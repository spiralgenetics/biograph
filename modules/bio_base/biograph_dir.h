#pragma once

#include <boost/filesystem.hpp>

#include "modules/io/transfer_object.h"
#include "json_spirit.h"

// {accession_id: md5_id_for_file}
typedef std::map<std::string, std::string> samples_t;

enum open_mode { READ_BGDIR, CREATE_BGDIR};

namespace fs = boost::filesystem;

struct biograph_metadata
{
	TRANSFER_OBJECT {
		VERSION(0);
		FIELD(version);
		FIELD(biograph_id);
		FIELD(accession_id);
		FIELD(samples);
		FIELD(command_history);
	};
	std::string version;
	std::string biograph_id;
	std::string accession_id;
	samples_t samples;
  std::vector<std::string> command_history;
};

class biograph_dir
{
public:
	// Default constructor
	biograph_dir() {};

	// Open for reading, create, or force create
	biograph_dir(const std::string& bg_dir, open_mode mode);

	// validate
	bool is_valid() const { return m_valid; };

	// Return the path to this biograph dir
	std::string path() const { return m_path.string(); };

	// Nominal path to the seqset. Does not check for existence.
	std::string seqset() const { return fs::path(m_path / fs::path("seqset")).string(); }

	// Nominal path to the given readmap. Does not check for existence.
	std::string readmap(const std::string& rm) const {
		return fs::path(m_path / fs::path("coverage") / fs::path(rm + ".readmap")).string();
	}

	// Return the (string) path to a readmap given a uuid or accession id.
	// Throws if no readmap can be found.
	std::string find_readmap(const std::string& rm);
	std::string find_readmap_accession(const std::string& rm);

	// Nominal path to the given seqpath. Does not check for existence.
	std::string seqpath(const std::string& sp) const {
		return fs::path(m_path / fs::path("assembly") / fs::path(sp + ".seqpath")).string();
	}

	// metadata
	const biograph_metadata& get_metadata() const { return m_metadata; };

	void set_metadata(biograph_metadata m) { m_metadata = m; };

	const std::string& biograph_id() const { return m_metadata.biograph_id; };
	const std::string& accession_id() const { return m_metadata.accession_id; };
	const samples_t& samples() const { return m_metadata.samples; };

	void save_metadata();

	//Checks if we'd be over writing existing files in mode CREATE_BGDIR
	static bool force_check(const std::string& bg_dir);

private:
	bool check_bgdir();
	bool load_metadata();

	fs::path m_path;
	bool m_valid = false;

	static const std::vector<std::string> m_subdirs;

	fs::path m_metadata_file = "metadata/bg_info.json";
	biograph_metadata m_metadata;
};
