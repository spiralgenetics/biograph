#pragma once

#include "modules/pipeline/datatype.h"
#include "modules/mapred/manifest.h"

// Represents the detailed map/reduce information associated with some dataset
class dataset_meta
{
public:
	dataset_meta() = default;

	dataset_meta(const std::string& sort)
		: the_manifest(sort)
	{}

	static
	void merge_add_element(dataset_meta& total, const std::string& part_url);

	TRANSFER_OBJECT
	{
		VERSION(0);
		FIELD(type, TF_STRICT);
		FIELD(the_manifest, TF_STRICT);
		FIELD(ref_name);
		FIELD(sort_keys);
		FIELD(in_progress, false);
		FIELD(blob, false);
	}

	datatype_ref type;
	std::string ref_name;
	manifest the_manifest;
	std::vector<std::string> sort_keys; // sequence of keys the dataset is sorted by.
	bool in_progress = false;
	bool blob = false;
	void validate();

	size_t get_output_index() const
	{
		return m_output_index;
	}

	void set_output_index(size_t output_index)
	{
		m_output_index = output_index;
	}

private:
	size_t m_output_index;
};
