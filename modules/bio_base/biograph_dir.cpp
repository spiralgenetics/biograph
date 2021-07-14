#include <boost/filesystem.hpp>
#include <fstream>

#include "modules/bio_base/biograph_dir.h"
#include "modules/bio_base/seqset.h"
#include "modules/io/file_io.h"

namespace js = json_spirit;
namespace fs = boost::filesystem;

const std::vector<std::string> biograph_dir::m_subdirs = {"metadata", "coverage", "qc"};

biograph_dir::biograph_dir(const std::string& bg_dir, open_mode mode) {
  if (mode == READ_BGDIR) {
    m_path = fs::path(bg_dir);
    m_valid = check_bgdir() && load_metadata();
    if (not m_valid) {
      throw std::runtime_error("Attempted to open " + bg_dir +
                             " but the BioGraph was not valid. Cannot continue.");
    }
    return;
  }

  m_path = fs::path(bg_dir);
  for (auto p : m_subdirs) {
    fs::create_directories(m_path / fs::path(p));
  }
  fs::create_directories(m_path / fs::path("analysis"));
  m_valid = check_bgdir();
  if (not m_valid) {
    throw std::runtime_error("Attempted to create " + bg_dir +
                             " but the resulting biograph was not valid. Cannot continue.");
  }

  m_metadata.version = biograph_current_version.make_string();
}

bool biograph_dir::check_bgdir() {
  if (not fs::is_directory(m_path)) {
    return false;
  }

  // We can't count on any file (even /seqset) to exist, but the directory
  // structure should be the same.
  for (auto p : m_subdirs) {
    if (not fs::exists(m_path / fs::path(p))) {
      return false;
    }
  }

  // TODO: if /seqset exists, check the header
  // TODO: if metadata/bg_info.json exists, validate it
  return true;
}

void biograph_dir::save_metadata() {
  if (not check_bgdir()) {
    throw io_exception("Can't save to an invalid biograph: " + m_path.string());
  }

  // If biograph_id is not set, get it from the seqset (if available)
  if (m_metadata.biograph_id.empty()) {
    fs::path seqset_path(m_path / fs::path("seqset"));
    if (fs::exists(seqset_path)) {
      seqset_file the_seqset_file(seqset_path.string());
      m_metadata.biograph_id = the_seqset_file.get_seqset().uuid();
    }
  }

  fs::path metadata_path(m_path / m_metadata_file);
  std::ofstream os(metadata_path.string());
  os << json_serialize(m_metadata);
  if (not os.good()) {
    throw io_exception("Could not write to " + metadata_path.string());
  }
  os.close();
}

bool biograph_dir::force_check(const std::string& bg_dir) {
  fs::path my_path(bg_dir);

  // Other checks can be put in if we decide to change
  // the behavior
  return fs::exists(my_path);
}

bool biograph_dir::load_metadata() {
  fs::path metadata_path(m_path / m_metadata_file);
  if (not fs::exists(metadata_path)) {
    return false;
  }

  try {
    json_deserialize(m_metadata, slurp_file(metadata_path.string()));
  } catch (...) {
    throw std::runtime_error("Could not parse biograph metadata: " + metadata_path.string());
  }

  return true;
}

std::string biograph_dir::find_readmap(const std::string& in_readmap) {
  if (in_readmap.empty()) {
    if (m_metadata.samples.size() == 0) {
      throw std::runtime_error(
          "No samples found in BioGraph metadata. You must specify the full readmap path.");
    }
    if (m_metadata.samples.size() > 1) {
      throw std::runtime_error(
          "Multiple samples are present. You must specify a readmap or accession ID.");
    }
    for (auto sample : m_metadata.samples) {
      return readmap(sample.second);
    }
  }

  // Just the UUID?
  if (fs::exists(readmap(in_readmap))) {
    return readmap(in_readmap);
  }

  // Accession ID?
  for (auto sample : m_metadata.samples) {
    if (in_readmap == sample.first) {
      return readmap(sample.second);
    }
  }

  // Path to an existing file? Just return it.
  if (fs::exists(in_readmap) && fs::is_regular_file(in_readmap)) {
    return in_readmap;
  }

  // I give up.
  throw std::runtime_error("Couldn't find " + in_readmap + ". Cannot continue.");
}

std::string biograph_dir::find_readmap_accession(const std::string& rm) {
  if (rm.empty()) {
    if (m_metadata.samples.size() == 0) {
      return "SAMPLE";
    }
    if (m_metadata.samples.size() > 1) {
      throw std::runtime_error(
          "Multiple samples are present. You must specify a readmap or accession ID.");
    }
    for (auto sample : m_metadata.samples) {
      return sample.first;
    }
  }
  // Accession ID?
  for (auto sample : m_metadata.samples) {
    if (rm == sample.first) {
      return rm;
    }
    if (rm == sample.second) {
      return sample.first;
    }
  }
  // Fail
  throw std::runtime_error("Couldn't find accession id associated with" + rm +
                           ". Cannot continue.");
}
