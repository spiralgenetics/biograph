#include "modules/mapred/manifest.h"
#include "modules/mapred/splitter.h"
#include "modules/io/log.h"

#include "modules/bio_base/struct_var.h"
#include <boost/integer/common_factor_rt.hpp>

void partition_info::add(const file_info& fi)
{
	files.push_back(fi);
	size += fi.size;
	num_records += fi.num_records;
}

void partition_info::add(const partition_info& pi)
{
	files.insert(files.end(), pi.files.begin(), pi.files.end());
	size += pi.size;
	num_records += pi.num_records;
}

class by_start_sorter
{
public:
	by_start_sorter(std::shared_ptr<sorter> sorter_)
		: m_sorter(sorter_)
	{}
	
	bool operator()(const file_info& f1, const file_info& f2)
	{
		int cmp_first = m_sorter->compare(f1.first_key, f2.first_key);
		if (cmp_first < 0) return true;
		if (cmp_first > 0) return false;
		int cmp_last = m_sorter->compare(f1.last_key, f2.last_key);
		return cmp_last < 0;
	}
private:
	std::shared_ptr<sorter> m_sorter;
};

void manifest::merge_tags(const manifest& in) 
{ 
	all_metadata.merge(in.all_metadata); 
}

// Basic deal: Creates two outputs, things that need sorting still, and new manifest of already sorted data
// Outputs true this is the 'final' sort
bool manifest::split_sort(manifest& done, std::vector<input_stream_params>& to_sort, size_t max_merge, bool clean_break)
{
	bool final = true;
	if (partitions.size() > 1)
		throw io_exception("Sorting isn't supported on manifests with > 1 partition");
	std::vector<file_info>& files = partitions[0].files;
	std::shared_ptr<sorter> sorter_ = sorter_registry::get(sort, "");
	by_start_sorter order(sorter_);
	std::sort(files.begin(), files.end(), order);
	input_stream_params cur_param;
	cur_param.sort = sort;
	cur_param.clean_break = clean_break;
	std::string encoding = get_encoding();
	cur_param.encoding = encoding;
	done.set_encoding(encoding);
	std::string highest_key = "";
	//SPLOG("Sort: Doing split_sort");
	for(size_t i = 0; i < files.size(); i++)
	{
		//SPLOG("Sort: On file %d, first_key = '%s', last_key = '%s', highest_key = '%s', begin_on = '%s'", 
		//		(int) i, files[i].first_key.c_str(), files[i].last_key.c_str(),
		//		highest_key.c_str(), cur_param.begin_on.c_str());
		if (cur_param.inputs.size() > 0 &&
			(clean_break ? 
				sorter_->compare(highest_key , files[i].first_key) == -2 :
				sorter_->compare(highest_key , files[i].first_key) <= 0
			))
		{
			// Found a split
			if (cur_param.inputs.size() == 1 && cur_param.begin_on == "")  
			{
				//SPLOG("Sort: Previous file is sorted, adding to 'done'");
				done.add(cur_param.inputs[0], 0);   // Single element, pass straight through
			}
			else
			{
				//SPLOG("Sort: Found break, group of %d files onto list", (int) cur_param.inputs.size()); 
				to_sort.push_back(cur_param);  // Hit the end of a sort group > 1, add to to_sort
			}

			// Clear group
			cur_param.inputs.clear();
			cur_param.num_records = 0;
			cur_param.begin_on= "";
			cur_param.end_before = "";
			highest_key = "";
		}
		if (cur_param.inputs.size() >= max_merge)  // Not at break, but need to dump sort group
		{
			// Count how many files cross the artificial boundary
			size_t count_good = 0;
			input_stream_params carry_on;
			carry_on.sort = sort;
			carry_on.encoding = encoding;
			carry_on.clean_break = clean_break;
			for(size_t j = 0; j < cur_param.inputs.size(); j++)
			{
				//SPLOG("Sort: Checking if last (%s) < new first (%s)", 
				//		cur_param.inputs[j].last_key.c_str(),
				//		files[i].first_key.c_str());
				if (clean_break ?
					sorter_->compare(cur_param.inputs[j].last_key, files[i].first_key) <= -2 :
					sorter_->compare(cur_param.inputs[j].last_key, files[i].first_key) <= -1
				)
				{
					count_good++;
					//SPLOG("Sort: Yes");
				}
				else 
				{
					carry_on.inputs.push_back(cur_param.inputs[j]);
					carry_on.num_records += cur_param.inputs[j].num_records;
					//SPLOG("Sort: No, pushing as redo");
				}
			}
			//SPLOG("Sort: Good count = %d", (int) count_good);
			// See if it's worth making a psudo-split (>= 70%)
			if (count_good >= (max_merge * 7 / 10))
			{
				cur_param.end_before = files[i].first_key;
				//SPLOG("Sort: Doing psudo-break: pushing group of %d files, retaining %d", 
				//		(int) cur_param.inputs.size(), (int) carry_on.inputs.size());
				to_sort.push_back(cur_param);
				cur_param = carry_on;
				cur_param.begin_on = files[i].first_key;
			}
			else
			{
				//SPLOG("Sort: Breaking on non-clean boundray, pushing group of %d file", (int) cur_param.inputs.size());
				to_sort.push_back(cur_param);
				cur_param.inputs.clear();
				cur_param.num_records = 0;
				cur_param.begin_on= "";
				highest_key = "";
				final = false;
			}
		}	
		cur_param.inputs.push_back(files[i]);
		cur_param.num_records += files[i].num_records;
		if (highest_key == "" || sorter_->compare(highest_key, files[i].last_key) < 0)
			highest_key = files[i].last_key;
	}
	if(!cur_param.inputs.empty())
	{// we need this check, otherwise split_sort will keep on calling to_sort.push_back with an empty cur_param
		// which will put a sort_task in an infinite loop.
		if (cur_param.inputs.size() == 1)
		{
			//SPLOG("Sort: Adding single trailing file to done");
			done.add(files[files.size() - 1], 0);
		}
		else
		{
			//SPLOG("Sort: Pusing trailing group of %d files", (int) cur_param.inputs.size());
			to_sort.push_back(cur_param);
		}
	}

	return final;
}

void manifest::split_sort_reduce(std::vector<input_stream_params>& out, size_t goal_size, bool clean_break)
{
	if (partitions.size() > 1)
		throw io_exception("Sorting isn't supported on manifests with > 1 parition");
	std::vector<file_info>& files = partitions[0].files;
	if (files.size() == 0) return;
	std::shared_ptr<sorter> sorter_ = sorter_registry::get(sort, "");
        by_start_sorter order(sorter_);
        std::sort(files.begin(), files.end(), order);

	input_stream_params cur_param;
	cur_param.sort = sort;
	cur_param.num_records = 0;
	cur_param.encoding = get_encoding();
	cur_param.clean_break = clean_break;
	size_t cur_total = 0;
	for(size_t i = 0; i < files.size(); i++)
	{
		// Check if I should try a break
		if (cur_total >= goal_size)
		{
			// See if I've split from the pack
			std::string back_key = sorter_->bump_back(files[i].first_key);
			size_t count_disjoint = 0;
			std::vector<file_info> keep_it;
			for(size_t j = 0; j < cur_param.inputs.size(); j++)
			{
				if (sorter_->compare(cur_param.inputs[j].last_key, back_key) < (clean_break ? -1 : 0))
					count_disjoint++;
				else
					keep_it.push_back(cur_param.inputs[j]);
			}
			if (count_disjoint > cur_param.inputs.size()/2)
			{
				// Good for a split
				cur_param.end_before = files[i].first_key;
				out.push_back(cur_param);
				cur_param.begin_on = files[i].first_key;
				cur_param.inputs = keep_it;
				cur_param.num_records = 0;
				cur_total = 0;
				for(size_t j = 0; j < keep_it.size(); j++)
				{
					cur_param.num_records += keep_it[j].num_records;
					cur_total += keep_it[j].size;
				}
			}
		}
		cur_param.inputs.push_back(files[i]);
		cur_param.num_records += files[i].num_records;
		cur_total += files[i].size;
	}
	cur_param.end_before = "";
	out.push_back(cur_param);
}

void manifest::add(const file_info& fi, size_t partition)
{
	partitions[partition].add(fi);
	this->size += fi.size;
	this->num_records += fi.num_records;
	if (m_notify) {
		m_notify();
	}
}

void manifest::add(const manifest& other, bool unsorted)
{
	std::string encoding = get_encoding();
	std::string other_encoding = other.get_encoding();
	if (other.num_records == 0 && other.size == 0) {
		return;
	}

	if (other.sort != sort && (!unsorted || sort != "")) {
		throw io_exception("Manifests failed to match sort order");
	}

	if (num_records == 0 && size == 0) {
		*this = other;
		return;
	}

	if (other_encoding != encoding) {
		std::string msg("manifest::add mismatched encodings: this manifest encoding: '"+encoding+"' other.encoding: '"+other_encoding+"'");
		throw io_exception(msg);
	}

	size_t new_num_parts = boost::integer::gcd(other.partitions.size(), partitions.size());
	if (new_num_parts != partitions.size()) {
		std::vector<partition_info> new_parts(new_num_parts);
		for(size_t i = 0; i < partitions.size(); i++)
			new_parts[i % new_num_parts].add(partitions[i]);
		std::swap(partitions, new_parts);
	}
	
	for (size_t i = 0; i < other.partitions.size(); i++) {
		partitions[i % new_num_parts].add(other.partitions[i]);
	}

	size += other.size;
	num_records += other.num_records;

	merge_tags(other);
}

void manifest::split_by_partition(std::vector<input_stream_params>& out) const
{
	out.resize(partitions.size());
	std::string encoding = get_encoding();
	for(size_t i = 0; i < partitions.size(); i++)
	{
		out[i].num_records = partitions[i].num_records;
		for(size_t j = 0; j < partitions[i].files.size(); j++)
			out[i].inputs.push_back(partitions[i].files[j]);
		out[i].sort = sort;
		out[i].encoding = encoding;
	}
}

void manifest::split_by_goal_size(std::vector<input_stream_params>& out, size_t goal_size) const
{
	size_t cur_size = 0;
	input_stream_params cur_param;
	cur_param.num_records = 0;
	cur_param.sort = sort;
	cur_param.encoding = get_encoding();
	
	for(size_t i = 0; i < partitions.size(); i++)
	{
		for(size_t j = 0; j < partitions[i].files.size(); j++)
		{
			const file_info& fi = partitions[i].files[j];
			cur_param.inputs.push_back(fi);
			cur_param.num_records += fi.num_records;
			cur_size += fi.size;
			if (cur_size > goal_size)
			{
				out.push_back(cur_param);
				cur_param.num_records = 0;
				cur_param.inputs.clear();
				cur_size = 0;
			}
		}
	}
	out.push_back(cur_param);
}

void manifest::split_by_splitter(manifest& target_manifest, std::vector<input_stream_params>& to_split, const std::string& the_splitter) const
{
	bool is_first_run = true;
	std::unique_ptr<splitter> splitter_ptr;
	std::string encoding = get_encoding();
	target_manifest.set_encoding(encoding);
	
	for(size_t i = 0; i < partitions.size(); i++)
	{
		for(size_t j = 0; j < partitions[i].files.size(); j++)
		{
			const file_info& the_file_info = partitions[i].files[j];
			if (is_first_run)
			{
				splitter_ptr = splitter_registry::get_safe(the_splitter, the_file_info.first_key);
				is_first_run = false;
			}
			else
			{
				splitter_ptr->set_initial_key(the_file_info.first_key);
			}
			
			if ((*splitter_ptr)(the_file_info.last_key))
			{
				input_stream_params split_me_param;
				split_me_param.inputs.push_back(the_file_info);
				split_me_param.num_records = 0;
				split_me_param.sort = sort;
				split_me_param.begin_on = the_file_info.first_key;
				split_me_param.encoding = encoding;
				to_split.push_back(split_me_param);
			}
			else
			{
				target_manifest.add(the_file_info, i);
			}
		}
	}
}

void manifest::split_mergepart(std::vector<input_stream_params>& out, size_t goal_size, size_t max_files) const
{
	std::string encoding = get_encoding();
	if (encoding.empty())
		SPLOG_P(LOG_DEBUG, "manifest::split_mergepart : encoding is empty");

	for(size_t i = 0; i < partitions.size(); i++)
	{
		size_t cur_size = 0;
		input_stream_params cur_param;
		cur_param.num_records = 0;
		cur_param.sort = sort;
		cur_param.encoding = encoding;
		for(size_t j = 0; j < partitions[i].files.size(); j++)
		{
			const file_info& fi = partitions[i].files[j];
			cur_param.inputs.push_back(fi);
			cur_param.num_records += fi.num_records;
			cur_size += fi.size;
			if (cur_size > goal_size || cur_param.inputs.size() > max_files)
			{
				out.push_back(cur_param);
				cur_param.num_records = 0;
				cur_param.inputs.clear();
				cur_size = 0;
			}
		}
		if (cur_param.num_records > 0)
			out.push_back(cur_param);
	}
}

void manifest::sort_file_infos()
{
	std::shared_ptr<sorter> the_sorter_ptr(sorter_registry::get_safe(sort, ""));
	
	for(unsigned int i = 0; i < partitions.size(); i++)
	{
		std::stable_sort(partitions[i].files.begin(),
			partitions[i].files.end(),
			[the_sorter_ptr](const file_info& lhs, const file_info& rhs) -> bool { return the_sorter_ptr->lt(lhs.first_key, rhs.first_key); }
		);
	}
}

size_t manifest::max_files() const
{
	size_t max_val = 0;
	for(size_t i = 0; i < partitions.size(); i++)
		max_val = std::max(max_val, partitions[i].files.size());

	return max_val;
}

manifest::manifest_iterator::manifest_iterator() : m_self(NULL), m_part(0), m_file(0) {} 
manifest::manifest_iterator::manifest_iterator(const manifest& self) : m_self(&self), m_part(0), m_file(0) 
{
	skip_to_valid();
} 

void manifest::manifest_iterator::skip_to_valid()
{
	while(m_part < m_self->partitions.size()) {
		if(m_file < m_self->partitions[m_part].files.size())
			return;
		m_part++;
		m_file = 0;
	}
	m_self = NULL;
	m_part = 0;
	m_file = 0;
}
	
void manifest::manifest_iterator::increment() 
{ 
	m_file++;
	skip_to_valid();
}

bool manifest::manifest_iterator::equal(manifest_iterator const& other) const
{
	return (m_self == other.m_self && m_part == other.m_part && m_file == other.m_file);
}

const file_info& manifest::manifest_iterator::dereference() const { return m_self->partitions[m_part].files[m_file]; }

manifest::const_iterator manifest::begin() const { return manifest_iterator(*this); }
manifest::const_iterator manifest::end() const { return manifest_iterator(); }

manifest_reader::manifest_reader(const manifest& m)
	: multi_reader<manifest::const_iterator>(m.begin(), m.end(), m.get_encoding())
{}

std::string manifest::get_encoding() const
{
	return metadata().get(meta::ns::internal, "encoding", std::string(""));
}

void manifest::set_encoding(const std::string& enc)
{
	metadata().set(meta::ns::internal, "encoding", enc);
}

size_t manifest::get_num_chunks() const
{
	size_t total = 0;
	for (const auto& partition : partitions) {
		total += partition.files.size();
	}
	return total;
}

size_t manifest::count_file_infos() const
{
	size_t file_info_count = 0;
	for(const auto partition : partitions) {
		file_info_count += partition.files.size();
	}
	
	return file_info_count;
}
