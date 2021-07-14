#include "modules/bio_format/struct_var.h"
#include "modules/bio_base/struct_var.h"
#include "modules/io/log.h"
#include "modules/io/registry.h"

#include <string.h>
#include <map>

REGISTER_3(exporter, struct_var, writable&, bool, const std::string&);

int sv_compute_edit_distance(const struct_var& sv, const reference& ref) {
	const int k_max_dist = 50;
	const reference_assembly& ref_assembly = ref.get_assembly();
	auto left_flat = long(ref_assembly.flatten(sv.ref_start));
	auto right_flat = long(ref_assembly.flatten(sv.ref_end));
	int ref_diff = -1;
	cost_matrix costs;

	if (sv.is_structural &&
		std::abs(right_flat - left_flat) < k_max_dist)
	{
		return k_max_dist;
	}

	if (sv.is_structural &&
		safe_range(left_flat, ref.size()) &&
		safe_range(right_flat, ref.size()))
	{
		dna_const_iterator left = ref.get_dna(left_flat);
		dna_const_iterator right = ref.get_dna(right_flat);
		if (sv.rev_start) left = left.rev_comp();
		if (sv.rev_end) right = right.rev_comp();
		std::vector<align_state> out;
		dna_sequence left_seq(left - k_max_dist, left + k_max_dist);
		dna_sequence right_seq(right - k_max_dist, right + k_max_dist);
		ref_diff = static_cast<int>(align_astar_exact(out, left_seq, right_seq, costs, 2.0 * k_max_dist));
	}
	return ref_diff;
}

void struct_var_exporter::write(const std::string& key, const std::string& value)
{
	struct_var sv;
	msgpack_deserialize(sv, value);
	std::string start_scaffold = m_reference_assembly.scaffold_order[sv.ref_start.scaffold_id];
	std::string end_scaffold = m_reference_assembly.scaffold_order[sv.ref_end.scaffold_id];
	int ref_diff = sv_compute_edit_distance(sv, m_reference);

	const double entropy = 2.0; // This field is deprecated, compute_entropy is actually broken.

	m_sink.print("%u\t", sv.var_id);
	m_sink.print("%s\t%ld\t%c\t", start_scaffold.c_str(), sv.ref_start.position, sv.rev_start ? '-' : '+');
	m_sink.print("%s\t%ld\t%c\t", end_scaffold.c_str(), sv.ref_end.position, sv.rev_end ? '-' : '+');
	m_sink.print("%s\t", sv.assembled.subseq(sv.var_start, sv.var_end - sv.var_start).as_string().c_str());
	m_sink.print("%s\t", sv.ref_seq.as_string().c_str());
	m_sink.print("%d\t%d\t", sv.is_structural, sv.is_ambig);
	m_sink.print("%ld\t%d\t", sv.depth, (int) sv.avg_depth);
	m_sink.print("%d\t%d\t", sv.min_overlap, (int) sv.avg_overlap);
	m_sink.print("%1.3f\t%d\t", entropy, ref_diff);
	m_sink.print("%s\n", sv.assembled.as_string().c_str());
}
