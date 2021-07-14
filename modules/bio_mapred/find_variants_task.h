#pragma once

#include "modules/io/transfer_object.h"
#include "modules/mapred/task.h"
#include "modules/mapred/manifest.h"

struct find_variants_params
{
	TRANSFER_OBJECT
	{
		VERSION(0);
		FIELD(reference);
		FIELD(min_width, TF_STRICT);
		FIELD(min_depth, TF_STRICT);
		FIELD(min_edit_distance, TF_STRICT);
		FIELD(read_support, false);
		FIELD(coverage, false);
		FIELD(min_mask_length, TF_STRICT);
		FIELD(score_gap, TF_STRICT);
		FIELD(max_simple_alignment, TF_STRICT);
		FIELD(word_size, TF_STRICT);
		FIELD(evalue, TF_STRICT);
		FIELD(pair_insert_size, TF_STRICT);
		FIELD(trace_only);
		FIELD(skip_anchor);
	}

	std::string reference;
	double evalue;
	double max_simple_alignment;
	size_t min_width;
	size_t min_depth;
	size_t min_edit_distance;
	size_t min_mask_length;
	size_t score_gap;
	size_t word_size;
	unsigned pair_insert_size;

	bool read_support = false;
	bool coverage = false;

	bool trace_only;
	bool skip_anchor;

	void validate();
};

class find_variants_task : public task_impl<find_variants_task>
{
public:
	static std::string s_type() { return "find_variants_task"; }

	TRANSFER_OBJECT
	{
		VERSION(0);
		FIELD(input, TF_STRICT);
		FIELD(params, TF_STRICT);
		FIELD(m_state, TF_STRICT);
		FIELD(m_subtask, TF_STRICT);
		FIELD(m_coverage, TF_STRICT);
		FIELD(m_outputs, TF_STRICT);
		FIELD(m_assemblies, TF_STRICT);
		FIELD(m_supports, TF_STRICT);
		FIELD(m_anchors, TF_STRICT);
	}

	void run();

	manifest input;
	find_variants_params params;

private:
	int m_state = 0; // 0 = nothing done, 1 = anchoring done, 2 = assembly done, 3 = call_variants done
	subtask_id m_subtask;
	manifest m_assemblies;
	manifest m_coverage;
	manifest m_variants;
	manifest m_supports;
	manifest m_anchors;
	std::vector<manifest> m_outputs;
};
