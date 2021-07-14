#include "modules/io/spiral_file.h"
#include "base/command_line.h"
#include "modules/io/json_transfer.h"
#include "modules/io/make_unique.h"
#include "modules/io/mmap_buffer.h"
#include "modules/io/uuid.h"
#include "tools/build_stamp.h"

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <fstream>

std::string spiral_file_file_info::command_line_str() {
  std::string ret = "";
  for (auto c : command_line) {
    if (c.find(" ") != std::string::npos) {
      // "" parameters with embedded spaces
      ret += "\"" + c + "\" ";
    } else {
      ret += c + " ";
    }
  }
  return ret;
}

const char spiral_file::k_part_info_pathname[] = "part_info.json";
const char spiral_file::k_file_info_pathname[] = "file_info.json";

spiral_file::~spiral_file() {}

void spiral_file_create_state::create_membuf(const std::string &partname,
                                             const membuf &contents) const {
  create_membuf(partname, contents, m_options);
}
void spiral_file_create_state::create_membuf(const std::string &partname, const membuf &contents,
                                             const spiral_file_options &options) const {
  m_top->create_path_contents(m_dir + partname, contents, options);
}

mutable_membuf spiral_file_create_state::create_membuf(const std::string &partname,
                                                       size_t size) const {
  return create_membuf(partname, size, m_options);
}

mutable_membuf spiral_file_create_state::create_membuf(const std::string &partname, size_t size,
                                                       const spiral_file_options &options) const {
  return m_top->create_path(m_dir + partname, size, options);
}

spiral_file_create_state spiral_file_create_state::create_subpart(
    const std::string &partname) const {
  return create_subpart(partname, m_options);
}

spiral_file_create_state spiral_file_create_state::create_subpart(
    const std::string &partname, const spiral_file_options &options) const {
  return spiral_file_create_state(m_top, m_dir + partname + "/", options);
}

void spiral_file_create_state::set_version(const std::string &part_type,
                                           const product_version &version) const {
  spiral_file_part_info part_info;
  part_info.part_type = part_type;
  part_info.version = version;
  m_version_set = true;
  create_json<spiral_file_part_info>(spiral_file::k_part_info_pathname, part_info);
}

void spiral_file_create_state::set_ephemeral_version(const std::string &part_type) const {
  spiral_file_part_info part_info;
  part_info.part_type = part_type;
  m_version_set = true;
  create_json<spiral_file_part_info>(spiral_file::k_part_info_pathname, part_info);
}

std::string spiral_file_create_state::uuid() const { return m_top->uuid(); }

spiral_file_create_state::~spiral_file_create_state() {
  if (std::uncaught_exception()) {
    SPLOG("Create of part %s failed due to exception", m_dir.c_str())
  } else {
    CHECK(m_version_set) << "No version specified when creating file part";
  }
}

spiral_file_create_state::spiral_file_create_state(spiral_file_create *top, const std::string &dir,
                                                   const spiral_file_options &options)
    : m_top(top), m_dir(dir), m_options(options) {}

spiral_file_create_state::spiral_file_create_state(spiral_file_create_state &&old) {
  m_top = old.m_top;
  m_version_set = old.m_version_set;
  m_dir = old.m_dir;

  // If we moved, don't complain about the version when the old one goes away.
  old.m_version_set = true;
}

void spiral_file_create::create_path_contents(const std::string &path, const membuf &contents,
                                              const spiral_file_options &options) {
  mutable_membuf new_contents = create_path(path, contents.size(), options);
  CHECK_EQ(new_contents.size(), contents.size());
  memcpy(new_contents.mutable_data(), contents.data(), contents.size());
}

bool spiral_file_open_state::membuf_present(const std::string &partname) const {
  return m_top->path_is_present(m_dir + partname);
}
bool spiral_file_open_state::subpart_present(const std::string &partname) const {
  return m_top->path_is_present(m_dir + partname + "/" + spiral_file::k_part_info_pathname);
}

bool spiral_file_open_state::is_mutable() const { return m_top->is_mutable(); }

spiral_file_open_state spiral_file_open_state::open_subpart(const std::string &partname) const {
  return open_subpart(partname, m_options);
}

spiral_file_open_state spiral_file_open_state::open_subpart(
    const std::string &partname, const spiral_file_options &options) const {
  return spiral_file_open_state(m_top, m_dir + partname + "/", options);
}

spiral_file_open_state::~spiral_file_open_state() {
  if (!std::uncaught_exception()) {
    CHECK(m_version_checked) << "Must check version number when opening a file.";
  }
}

spiral_file_open_state::spiral_file_open_state(spiral_file_open_state &&old) {
  m_top = old.m_top;
  m_version_checked = old.m_version_checked;
  m_dir = old.m_dir;

  // If we moved, don't complain about the version when the old one goes away.
  old.m_version_checked = true;
}

void spiral_file_open_state::enforce_max_version(const std::string &part_type,
                                                 const product_version &enforce_version) const {
  spiral_file_part_info pi = part_info();
  if (part_type != pi.part_type) {
    throw(io_exception("Expecting type " + part_type + "; got " + pi.part_type));
  }
  if (!enforce_version.can_read(pi.version)) {
    throw(io_exception("Version " + pi.version.make_string() + " newer than supported version " +
                       enforce_version.make_string()));
  }
  m_version_checked = true;
}

void spiral_file_open_state::enforce_ephemeral_version(const std::string &part_type) const {
  spiral_file_part_info pi = part_info();
  CHECK_EQ(part_type, pi.part_type) << "Wrong file path provided";
  m_version_checked = true;
}

spiral_file_part_info spiral_file_open_state::part_info() const {
  return open_json_internal<spiral_file_part_info>(spiral_file::k_part_info_pathname, m_options);
}

spiral_file_file_info spiral_file_open_state::file_info() const { return m_top->file_info(); }

membuf spiral_file_open_state::open_membuf(const std::string &partname) const {
  CHECK(m_version_checked);
  return open_membuf_internal(partname, m_options);
}

membuf spiral_file_open_state::open_membuf(const std::string &partname,
                                           const spiral_file_options &options) const {
  CHECK(m_version_checked);
  return open_membuf_internal(partname, options);
}

membuf spiral_file_open_state::open_membuf_internal(const std::string &partname,
                                                    const spiral_file_options &options) const {
  return m_top->get_path(m_dir + partname, options);
}

mutable_membuf spiral_file_open_state::open_mutable_membuf(const std::string &partname) const {
  return open_mutable_membuf(partname, m_options);
}

mutable_membuf spiral_file_open_state::open_mutable_membuf(
    const std::string &partname, const spiral_file_options &options) const {
  CHECK(m_version_checked);
  return m_top->get_mutable_path(m_dir + partname, options);
}

spiral_file_open_state::spiral_file_open_state(spiral_file_open *top, const std::string &dir,
                                               const spiral_file_options &options)
    : m_top(top), m_dir(dir), m_options(options) {}

void spiral_file_create::write_file_info() {
  spiral_file_file_info file_info;
  file_info.create_timestamp = time(0);
  char ctime_buf[100];
  file_info.create_timestamp_text = ctime_r(&file_info.create_timestamp, ctime_buf);
  file_info.create_timestamp_text.pop_back();  // Remove \n
  file_info.build_revision = get_build_scm_revision();
  file_info.build_is_clean = build_is_clean();
  file_info.build_host = get_build_host();
  file_info.build_user = get_build_user();
  file_info.build_timestamp = get_build_timestamp();
  file_info.build_timestamp_text = ctime_r(&file_info.build_timestamp, ctime_buf);
  file_info.build_timestamp_text.pop_back();  // Remove \n
  file_info.uuid = m_uuid;
  file_info.command_line = original_program_args();
  if (!file_info.build_timestamp && !getenv("GTEST_TMP_DIR")) {
    // Missing build info is probably ok for tests, but make sure we
    // warn so we don't accidentally get a production binary produced
    // without a build stamp without people noticing.
    std::cerr << "WARNING: Binary was not compiled with build stamps enabled.\n"
                 "Version info in generated output file file will be missing.\n";
    SPLOG(
        "WARNING: Binary was not compiled with build stamps enabled. "
        "Version info in generated output file file will be missing.");
  }
  create_path_contents(spiral_file::k_file_info_pathname,
                       owned_membuf::from_str(json_serialize(file_info), "spiral_file_json"),
                       m_options);
}

spiral_file_create::spiral_file_create(const spiral_file_options &options)
    : m_uuid(make_uuid()), m_options(options) {}

spiral_file_create_state spiral_file_create::create() {
  write_file_info();

  return spiral_file_create_state(this, "", m_options);
}

spiral_file_open::spiral_file_open(const spiral_file_options &options) : m_options(options) {}

spiral_file_open_state spiral_file_open::open(const std::string &part_path) {
  return spiral_file_open_state(this, part_path, m_options);
}

spiral_file_file_info spiral_file_open::file_info() const {
  membuf file_info_buf = get_path(spiral_file::k_file_info_pathname, m_options);
  spiral_file_file_info file_info;
  json_deserialize(file_info, file_info_buf.str());
  return file_info;
}

std::string spiral_file_open::uuid() const { return file_info().uuid; }

std::string spiral_file_open_state::uuid() const { return file_info().uuid; }

std::ostream &operator<<(std::ostream &os, const spiral_file_options &opts) {
  os << "delayed_write=" << opts.delayed_write;
  os << ", small_object_threshold=" << opts.small_object_threshold;
  os << ", read_into_ram=" << opts.read_into_ram;
  return os;
}

spiral_file_options spiral_file_options::with_delayed_write(bool new_delayed_write) const {
  spiral_file_options modified = *this;
  modified.delayed_write = new_delayed_write;
  return modified;
}

spiral_file_options spiral_file_options::with_read_into_ram(bool new_read_into_ram) const {
  spiral_file_options modified = *this;
  modified.read_into_ram = new_read_into_ram;
  return modified;
}
