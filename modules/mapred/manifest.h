#pragma once

#include "base/base.h"
#include "modules/mapred/path.h"
#include "modules/mapred/input_stream.h"
#include "modules/mapred/multi_reader.h"
#include "modules/mapred/metadata.h" 
#include "modules/io/keyvalue.h"
#include "modules/io/transfer_object.h"

#include <boost/iterator/iterator_facade.hpp>

// Manifests are an abstract class that represents a group of paths which are logically
// the same data set.  If the data is sorted, the manifest will contain the sort function
// along with first/last key for each part (which facilitates distributed merge sorts)
// Manifests are futher divided into 'paritions' which each contain only records who's
// key % num_paritions == parition number.  In the case of unsorted data, there is only
// one parition.

class file_info
{
public:
	file_info() = default;

	file_info(
		const path& file, 
		size_t size, 
		size_t num_records, 
		const std::string& first_key = "",	
		const std::string& last_key = "")
		: file(file)
		, size(size)
		, num_records(num_records)
		, first_key(first_key)
		, last_key(last_key) 
	{}

	TRANSFER_OBJECT
	{
		VERSION(0);
		FIELD(file, TF_STRICT); 
		FIELD(size, TF_STRICT); 
		FIELD(num_records, TF_STRICT); 
		FIELD(first_key); 
		FIELD(last_key); 
	}

	operator path() const { return file; }

	path file;
	size_t size = 0; // total bytes if one were to read AND decode the entire file. The encoding used is saved in the manifest
	size_t num_records = 0; // Number of records in this part
	std::string first_key; // First key in file (empty if not sorted)
	std::string last_key; // Last key in file (empty if not sorted)
};

class partition_info
{
public:
	partition_info() = default;

	TRANSFER_OBJECT
	{ 
		VERSION(0);
		FIELD(files, TF_STRICT); 
		FIELD(num_records, TF_STRICT); 
		FIELD(size, TF_STRICT); 
	}

	void add(const file_info& fi);
	void add(const partition_info& other);
	
	size_t size = 0;  // Precomputed accumulation of files
	size_t num_records = 0; // Precomputed accumulation of files
	std::vector<file_info> files;
};

class splitter;
class manifest
{
	class manifest_iterator;
	friend class manifest_iterator;
public:
	typedef std::function<void ()> notify_f;
	using is_manifest = std::true_type;

	TRANSFER_OBJECT 
	{ 
		VERSION(0);
		FIELD(partitions, TF_STRICT); 
		FIELD(size, TF_STRICT); 
		FIELD(num_records, TF_STRICT); 
		FIELD(sort); 
		FIELD(meta);
		FIELD(tags); // it's pronounced 'tehgs' :)
		FIELD(all_metadata);
	}

	manifest(const std::string& sort_ = "", size_t num_partitions = 1) 
		: size(0)
		, num_records(0)
		, sort(sort_)
		, partitions(num_partitions) 
	{}

	void set_meta(const std::string& meta_) { meta = meta_; }
	void set_sort(const std::string& sort_) { sort = sort_; }
	void set_notify(notify_f fn) { m_notify = fn; }

	std::string get_encoding() const;
	void set_encoding(const std::string& enc);

	std::string get_sort() { return sort; }
	const std::string& get_meta() { return meta; }
	void add(const manifest& other, bool unsorted = false);
	void add(const file_info& fi, size_t partition);
	size_t get_size() const { return size; }
	size_t get_num_records() const { return num_records;}
	const std::string& get_sort() const { return sort; }
	size_t get_num_partitions() const { return partitions.size(); }
	size_t get_num_chunks() const;
	size_t count_file_infos() const;

	typedef manifest_iterator const_iterator;
	typedef path value_type;

	void split_by_partition(std::vector<input_stream_params>& out) const;
	void split_by_goal_size(std::vector<input_stream_params>& out, size_t goal_size) const;
	void split_mergepart(std::vector<input_stream_params>& out, size_t goal_size, size_t max_files) const;
	bool split_sort(manifest& done, std::vector<input_stream_params>& to_sort, size_t max_merge, bool clean_break = true);
	void split_sort_reduce(std::vector<input_stream_params>& out, size_t goal_size, bool clean_break = true);
	// Partition file_infos if they need splitting.  Into to_split if yes, otherwise into target_manifest.
	void split_by_splitter(manifest& target_manifest, std::vector<input_stream_params>& to_split, const std::string& the_splitter) const;
	void sort_file_infos();
	size_t max_files() const; // Returns the maximum number of files per partition
	const_iterator begin() const;
	const_iterator end() const;
	
	void merge_tags(const manifest& in);

	meta::data& metadata() { return all_metadata; }
	const meta::data& metadata() const { return all_metadata; }

	// Functions that take an arbitrary number of manifests and merge the metadata from the
	// second parameter and on into the manifest in the first parameter.
	// Call like this: update_metadate(dest_manifest, source_manifest1, source_manifest2, ..., source_manifest_n);
	template<typename T>
	void update_metadata(const T& in)
	{
		static_assert(T::is_manifest::value, "The update_metadata function takes only manifests.  The input is not a manifest.");
		merge_tags(in);
	}

	template<typename T, typename... TInManifests>
	void update_metadata(const T& in, const TInManifests&... inputs)
	{
		static_assert(T::is_manifest::value, "The update_metadata function takes only manifests.  One of the inputs is not a manifest.");
		merge_tags(in);
		update_metadata(inputs...);
	}

private:
	class manifest_iterator : public boost::iterator_facade<manifest_iterator, file_info const, boost::forward_traversal_tag>
	{
	public:
		manifest_iterator();
		manifest_iterator(const manifest& self);
	private:
		friend class boost::iterator_core_access;

		void skip_to_valid();
		void increment();
		bool equal(manifest_iterator const& other) const;
		const file_info& dereference() const;

		const manifest* m_self;
		size_t m_part;
		size_t m_file;
	};

	size_t size;        // Precomputed accumulation of partitions
	size_t num_records;  // Precomputed accumulation of partitions
	std::string sort;  // If non-empty, elements are sorted
	std::vector<partition_info> partitions;
	std::string meta;
	meta::data all_metadata;
	notify_f m_notify;

public:
	std::map<std::string, std::string> tags;
};

class manifest_reader : public multi_reader<manifest::const_iterator>
{
public:
	manifest_reader(const manifest& m);
};
