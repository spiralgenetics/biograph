#pragma once

#include "modules/pipeline/dataset_meta.h" 
#include "modules/pipeline/direntry.h"
#include "modules/mapred/path.h" 
#include "modules/io/config.h"

// Manages path translation

class dataset_path
{
public:
	// Generate a path from a URL-encoded URL
	dataset_path(std::string url, bool decode = true);

	dataset_path append(const std::string& name) const;  // Get sub-path (invalid if not a directory)	
	dataset_path root() const; // Return the 'root' of this path
	static dataset_path root(const std::string& user); 

	bool is_reference() const { return m_reference; }
	const std::string& url() const { return m_url; } // Get as a URL again		
	const std::string& user() const { return m_user; } // Returns the user associated with this path, null for reference
	const std::string& parent() const { return m_parent; } // Returns the URL of the parent of this path, empty for root paths
	const std::string& name() const { return m_name; } // Returns the name of this path, will be empty for 'root' paths (ie, api/users/melvin or /api/reference)
	std::string friendly() const;
	const path& meta() const { return m_meta; }
	const path& data() const { return m_data_dir; } // Returns the 'path' where data should be placed
	const path& base() const { return m_base; }

	path::exist_enum exists() const; // Get data from cache
	void load(dataset_meta& meta) const { m_meta.json_get(meta); }
	direntry stat() const; // Get the direntry data only
	void create(const dataset_meta& meta) const; // Do a JSON put and update cache
	void update(const dataset_meta& meta) const; // Overwrite an existing file
	void remove(bool recursive = false) const; // Do a JSON put and update cache
	void mkdir() const;
	direntry create_remote(const dataset_meta& meta) const;  // Used to do write to meta-data on a remote machine
	static void update_cache(const direntry& de);  // Used to update cache for changes made remotely
	std::vector<direntry> list_dir() const;

private:
	dataset_path() = default;

	bool m_reference;
	std::string m_url;
	std::string m_user;  // Empty if reference
	std::string m_parent; // Parent url, empty if root
	std::string m_name;  
	path m_base;
	path m_meta;
	path m_data_dir;
};

void copy_dataset(const dataset_path& out, const dataset_path& in);
void create_ancestors(const dataset_path& p);
void set_directory_direntry(const dataset_path& path, const time_t t, direntry& de);
void set_file_direntry(const dataset_path& path, const dataset_meta& meta, const time_t t, direntry& de);

void gen_cache(const char* user);
