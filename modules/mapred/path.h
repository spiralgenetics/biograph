#pragma once

#include "modules/io/json_transfer.h"
#include <boost/operators.hpp>
#include <ctime>

class waiter
{
public:
	virtual ~waiter() = default;

	// returns the MD5 hash encoded in base64
	virtual std::string wait() = 0;
};

// Options to use when writing a path
struct path_write_options {
public:
	// Tags to use when writing a path; currently only supported in s3.
	std::map<std::string, std::string> tags;

	// Default write options.  Returns a reference so this is easy
	// to use for default arguments.
	static path_write_options& defaults() {
		static path_write_options default_options;
		return default_options;
	}
};

struct path_impl;
class path : public boost::totally_ordered<path>
{
	friend struct path_impl;
public:
	path();
	~path();

	path(const path& rhs);
	path& operator=(const path&);

	path(const std::string& url);
	operator std::string() const { return url(); } // Cast operator for serialization
	std::string url() const;
	std::string bare_path() const;
	std::string filename() const;
	
	bool valid() const { return (bool)m_impl; }

	// Return true if this path is excluded. Exclusion list is defined per product in config.json file
	// This is weak security feature, primarily to avoid accidentally running a destructive operations on a path
	bool excluded() const;

	typedef enum { UNKNOWN = 0, FILE, S3 } path_type;
	path_type type() const;

	bool operator<(const path& rhs) const;
	bool operator==(const path& rhs) const;

	// A helper so we dont have to do a much of hand manipulation
	// Adds another path component (ie += '/' + more)
	path append(const std::string& more) const;

	// Append a new 'unique' name to the path
	path append_unique(const std::string& prefix) const;

	// Get readers/writers
	// whoever calls these owns the returned pointer
	std::unique_ptr<readable> read() const;
	std::unique_ptr<writable> write(const path_write_options& options = path_write_options::defaults()) const;

	// Support for 'inverted' mode to allow fast s3 stuffs
	std::unique_ptr<waiter> write_inverted(reset_readable& source, size_t size, const path_write_options& options = path_write_options::defaults()) const;

	// Move and copy, move 'across' machines/types is not supported
	static void move(const path& src, const path& dest);
	static void copy(const path& src, const path& dest, const path_write_options& options = path_write_options::defaults());

	enum exist_enum
	{
		e_no_exist,
		e_file,
		e_directory
	};

	enum walk_state
	{
		w_dir_enter,  // Called when entering a directory
		w_dir_leave,  // Called when leaving a directory
		w_file    // Called for each file
	};

	static std::string str(const exist_enum& e);
	exist_enum exists() const;
	std::time_t modify_time() const;
	size_t size() const;
	std::vector<std::string> list() const;

	// Build a directory here
	void mkdir() const;

	// True if removed, false if wasn't there, throw on error
	bool remove() const;

	// True if removed, false if wasn't there, throw on error (including not-empty)
	bool rmdir(bool recursive = false) const;

	// Define the type for a walk callback
	struct walk_params
	{
		walk_params(
			walk_state state,
			const path& node,
			std::time_t last_modified,
			size_t size)
			: state(state)
			, node(node)
			, last_modified(last_modified)
			, size(size)
		{}

		walk_state  state;
		const path& node;
		std::time_t last_modified;
		size_t      size;
	};
	typedef std::function<void (const walk_params& params)> walker_f;

	// Walk this path recursively, calling the 'callback' for each subpath.  For directories, w_dir_enter is called
	// first, then the entries are walked, then w_dir_leave is called.  Files are called via w_file.
	// In addition, the modification time and size is passed for each entry.  The meaning of size for directories
	// is undefined. Stops when all paths have been visited.
	void walk(walker_f fn) const;

	// Convenience methods to allow slurping up a whole strings as a file
	std::string get() const;
	void put(const std::string& value, const path_write_options& options = path_write_options::defaults()) const;

	// Read and write JSON blobs into files with ease
	template<class X>
	void json_get(X& out) const
	{
		std::string s = get();
		json_deserialize(out, s);
	}

	template<class X>
	void json_put(const X& in, const path_write_options& options = path_write_options::defaults()) const
	{
		std::string s = json_serialize(in);
		put(s);
	}

private:
	path(std::unique_ptr<path_impl> impl);
	std::unique_ptr<path_impl> m_impl;
};

// Make path turn into a normal string for serialization

BASE_TYPE(path, std::string);
