#include <memory>

#include "modules/io/make_unique.h"
#include "modules/bio_mapred/kmerize_ref.h"
#include "modules/bio_base/reference.h"
#include "modules/mapred/sort_task.h"
#include "modules/mapred/output_stream.h"

REGISTER_TASK(kmerize_ref_task);
REGISTER_TASK(kmerize_supercontig_task);

void kmerize_ref_task::run()
{
	if (m_state == 0)
	{
		// Split progress
		split_progress(0.01, 0.4);

		// Make many kmerize_supercontig_tasks
		reference ref(params.reference);
		const reference_assembly& ref_assembly = ref.get_assembly();
		for (const supercontig& sc : ref_assembly.supercontigs) {
			std::unique_ptr<kmerize_supercontig_task> kst = make_unique<kmerize_supercontig_task>();
			kst->params = params;
			kst->the_supercontig = sc.name;
			m_subtasks.push_back(add_subtask(std::move(kst)));
		}
		// Move state forward
		m_state = 1;
	}
	else if (m_state == 1)
	{
		// Split progress
		split_progress(0.01, 0.01);

		// Collect results
		manifest out("lexical");
		for(size_t i = 0; i < m_subtasks.size(); i++)
		{
			manifest subout;
			get_output(subout, m_subtasks[i]);
			out.add(subout);
		}

		// Clear old subtasks and add sort
		m_subtasks.clear();
		std::unique_ptr<sort_task> st = make_unique<sort_task>();
		st->input = out;
		st->reduce = "sum";
		st->is_summary = true;
		m_subtasks.push_back(add_subtask(std::move(st)));

		// Move state forward
		m_state = 2;
	}
	else if (m_state == 2)
	{
		// Get sort result and return
		manifest out;
		get_output(out, m_subtasks[0]);
		set_output(out);
	}
}

void kmerize_supercontig_task::run()
{
	reference ref(params.reference);

	manifest out;
	output_stream_params osp;
	osp.sort = "lexical";
	osp.reduce = "sum";
	std::unique_ptr<kv_sink> sink = osp.build(get_root(), "ref_kmers", out);

	size_t ks = params.kmer_size;
	const supercontig& sc = ref.get_assembly().get_supercontig(the_supercontig);
	size_t c = 0;
	if (sc.len >= ks) {
		dna_const_iterator it = ref.get_dna(sc.tot_offset);
		dna_const_iterator itEnd = ref.get_dna(sc.tot_offset + sc.len - (ks - 1));
		while(it != itEnd) {
			update_progress(double(c++) / double(sc.len));
			dna_sequence seq(it, it + ks);
			seq.canonicalize();
			sink->write_msgpack(seq, 1);
			it++;
		}
	}

	sink->close();
	set_output(out);
}

