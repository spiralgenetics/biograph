#pragma once

#include "base/base.h"
#include "modules/io/json_transfer.h"
#include "modules/io/membuf.h"
#include "modules/io/transfer_object.h"
#include "modules/io/version.h"

// Utils to operate on spiral binary file formats.
//
// Goals of the file format:
// * Easily mmapable, both on input and output.
// * Easy to deal with versioning
//
// We use "zip" file formats using libarchive, with each file in the
// zip being uncompressed so it can be mmaped.
//
// In each directory in the archive, there is a "version.json" textual
// JSON file that contains a "product_version" structure.  This
// contains the version information for that particular part.
//
// Versioned subparts are in their own subdirectories, with their own
// version.json files.
//
// Top level there's also a "build_info.json" textual JSON file
// containing build stamp information about the binary that produced
// the file.

// File format for "part_info.json".
struct spiral_file_part_info {
  TRANSFER_OBJECT {
    VERSION(0);
    FIELD(part_type);
    FIELD(version);
  };

  std::string part_type;
  product_version version;
};

// File format for "file_info.json".
struct spiral_file_file_info {
  TRANSFER_OBJECT {
    VERSION(0);
    FIELD(create_timestamp);
    FIELD(create_timestamp_text);
    FIELD(uuid);
    FIELD(build_revision);
    FIELD(build_is_clean);
    FIELD(build_host);
    FIELD(build_user);
    FIELD(build_timestamp);
    FIELD(build_timestamp_text);
    FIELD(command_line);
  };

  time_t create_timestamp = 0;
  std::string create_timestamp_text;

  // Git revision the generating binary was built at.
  std::string build_revision;

  // True if compiled from an unmodified source tree at build_revision.
  bool build_is_clean = false;

  // Unique identifier for this file or part.
  std::string uuid;

  // Host and user that built this binary.
  std::string build_host;
  std::string build_user;

  // Time this binary was built.
  time_t build_timestamp = 0;
  std::string build_timestamp_text;

  // Command line used to generate this file.  First element is
  // argv[0].
  std::vector<std::string> command_line;

  std::string command_line_str();
};

struct spiral_file {
  // The name of the file in the archive that contains the version
  // information for that part.  Each part has this file present.
  static const char k_part_info_pathname[];

  // The name of the file in the archive that contains metadata for
  // the whole file.  Only one of these exists total.
  static const char k_file_info_pathname[];

  virtual ~spiral_file();
};

struct spiral_file_options {
  // If delayed_write is true, delay writing of all items until
  // they're closed.  This lets us avoid thrashing dirty buffers if
  // we're filling in a buffer non-sequentially.
  bool delayed_write = true;
  spiral_file_options with_delayed_write(bool new_delayed_write) const;

  // If delayed_write is false, we still delay write of objects
  // smaller than this amount to avoid having a bunch of mmap calls for tiny
  // regions
  size_t small_object_threshold = 4096;

  // If true, read buffers into RAM instead of using mmap.  This can
  // help when mmap performance is poor (like with gpfs).
  bool read_into_ram = false;
  spiral_file_options with_read_into_ram(bool new_read_into_ram) const;
};

std::ostream &operator<<(std::ostream &os, const spiral_file_options &opts);

class spiral_file_create;
class spiral_file_create_state {
 public:
  void create_membuf(const std::string &partname, const membuf &contents) const;
  void create_membuf(const std::string &partname, const membuf &contents,
                     const spiral_file_options &options) const;

  mutable_membuf create_membuf(const std::string &partname, size_t size) const;
  mutable_membuf create_membuf(const std::string &partname, size_t size,
                               const spiral_file_options &options) const;

  spiral_file_create_state create_subpart(const std::string &partname) const;
  spiral_file_create_state create_subpart(const std::string &partname,
                                          const spiral_file_options &options) const;
  ~spiral_file_create_state();

  void set_version(const std::string &part_type, const product_version &version) const;

  // Specifies that this part should not be used by processes other than the one
  // that created it.  I.e., it's a temp file.
  void set_ephemeral_version(const std::string &part_type) const;

  // Writes serialized JSON data.  If subpart already exists, it must
  // contain the same data.
  template <typename T>
  void create_json(const std::string &partname, const T &new_value) const {
    create_membuf(partname, owned_membuf::from_str(json_serialize(new_value), "spiral_file_json"),
                  m_options);
  }
  template <typename T>
  void create_json(const std::string &partname, const T &new_value,
                   const spiral_file_options &options) const {
    create_membuf(partname, owned_membuf::from_str(json_serialize(new_value), "spiral_file_json"),
                  options);
  }

  spiral_file_create_state(spiral_file_create_state &&);

  // Returns the UUID for the file being created.  This is the same for
  // all parts in the file.
  std::string uuid() const;

  const spiral_file_options &options() const { return m_options; }

 private:
  friend class spiral_file_create;
  spiral_file_create_state(spiral_file_create *top, const std::string &dir,
                           const spiral_file_options &options);
  spiral_file_create_state(const spiral_file_create_state &) = delete;
  spiral_file_create_state &operator=(const spiral_file_create_state &) = delete;
  void write_file_info();

  spiral_file_create *m_top = nullptr;
  mutable bool m_version_set = false;
  std::string m_dir;
  spiral_file_options m_options;
};

class spiral_file_open;
class spiral_file_open_state {
 public:
  // Returns true if the named membuf exists as part of this part.
  bool membuf_present(const std::string &partname) const;

  // Returns true if the named subpart exists already.
  bool subpart_present(const std::string &partname) const;

  // True if file isn't read-only.
  bool is_mutable() const;

  spiral_file_open_state open_subpart(const std::string &partname) const;
  spiral_file_open_state open_subpart(const std::string &partname,
                                      const spiral_file_options &options) const;
  ~spiral_file_open_state();

  // Raises an exception if this part has a version more recent than
  // the supplied version.
  void enforce_max_version(const std::string &part_type, const product_version &version) const;

  // Raises an exception if this part was not generated by this process.
  void enforce_ephemeral_version(const std::string &part_type) const;

  // Returns the stored version and type for this part.
  spiral_file_part_info part_info() const;

  // Returns the global file information for this file.  This is
  // shared among all parts.
  spiral_file_file_info file_info() const;

  // Returns the UUID for the file being opened.  This is the same for
  // all parts in the file.
  std::string uuid() const;

  // Provides access to a raw membuf subpart.
  membuf open_membuf(const std::string &partname) const;
  membuf open_membuf(const std::string &partname, const spiral_file_options &options) const;
  mutable_membuf open_mutable_membuf(const std::string &partname) const;
  mutable_membuf open_mutable_membuf(const std::string &partname,
                                     const spiral_file_options &options) const;

  // Provides access to a raw membuf subpart as serialized data.
  template <typename T>
  T open_json(const std::string &partname) const {
    CHECK(m_version_checked);
    return open_json_internal<T>(partname, m_options);
  }
  template <typename T>
  T open_json(const std::string &partname, const spiral_file_options &options) const {
    CHECK(m_version_checked);
    return open_json_internal<T>(partname, options);
  }

  spiral_file_open_state(spiral_file_open_state &&);

  const spiral_file_options &options() const { return m_options; }

 private:
  // Version of open_json that doesn't make sure we've checked the
  // part version.  That way we can use this as part of the part
  // versioning code.
  template <typename T>
  T open_json_internal(const std::string &partname, const spiral_file_options &options) const {
    T result;
    json_deserialize(result, open_membuf_internal(partname, options).str());
    return result;
  }
  membuf open_membuf_internal(const std::string &partname,
                              const spiral_file_options &options) const;

  friend class spiral_file_open;
  spiral_file_open_state(spiral_file_open *top, const std::string &dir,
                         const spiral_file_options &options);
  spiral_file_open_state(const spiral_file_open_state &) = delete;
  spiral_file_open_state &operator=(const spiral_file_open_state &) = delete;

  mutable bool m_version_checked = false;

  spiral_file_open *m_top = nullptr;
  std::string m_dir;
  spiral_file_options m_options;
};

class spiral_file_create : public spiral_file {
 public:
  spiral_file_create_state create();

  std::string uuid() const { return m_uuid; }

 protected:
  spiral_file_create(const spiral_file_options &options);

  friend class spiral_file_create_state;

  void write_file_info();
  virtual mutable_membuf create_path(const std::string &path, size_t size,
                                     const spiral_file_options &options) = 0;
  virtual void create_path_contents(const std::string &path, const membuf &contents,
                                    const spiral_file_options &options);

  std::string m_uuid;
  spiral_file_options m_options;
};

class spiral_file_open : public spiral_file {
 public:
  // Opens a given part.  Normally, a path is not specified and this
  // allows access to the top-level part.
  spiral_file_open_state open(const std::string &part_path = "");

  // Returns this file's build and run information.
  spiral_file_file_info file_info() const;

  // Returns the UUID associated with this file.
  std::string uuid() const;

  // Lists the files available in this archive.  This should not
  // normally be used except during testing or debugging.
  virtual std::set<std::string> contents() const = 0;

  // Returns true if this archive is open for writing.
  virtual bool is_mutable() const = 0;

  virtual bool path_is_present(const std::string &path) const = 0;

 protected:
  spiral_file_open(const spiral_file_options &options);

  friend class spiral_file_open_state;

  virtual membuf get_path(const std::string &path, const spiral_file_options &options) const = 0;
  virtual mutable_membuf get_mutable_path(const std::string &path,
                                          const spiral_file_options &options) = 0;

  spiral_file_options m_options;
};
