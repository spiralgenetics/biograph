#include "modules/bio_format/make_vars.h"
#include "modules/bio_base/call_structural.h"
#include "modules/bio_base/struct_var.h"
#include "modules/bio_base/pileup.h"
#include <vector>

#include <boost/format.hpp>
#include <boost/random.hpp>

using boost::format;
using boost::str;

namespace {
const unsigned k_max_var_a_start_cost = 200;

std::atomic<size_t> g_no_actual_depth{0};
std::atomic<size_t> g_no_variations_must_be_exact_match{0};
std::atomic<size_t> g_zero_seq_begin{0};
}

bool struct_var_adapter(
	const reference& ref,
    const std::function<void(const struct_var& var)>& sink_f,
	const dna_sequence& var_seq,
	const dna_slice& lower,
	const dna_slice& upper,
	ipileup* pile,
	const struct_var& per_assembly,
	size_t min_depth,
	bool dup_structural)
{
	bool is_var = false;

	auto vars = call_structural(var_seq, lower, upper, k_max_var_a_start_cost);
	if (vars.empty()) {
      g_no_variations_must_be_exact_match.fetch_add(1);
		return false;
	}
	bool has_anything = false;
	std::vector<size_t> depths(var_seq.size());
	std::vector<size_t> afwd(var_seq.size());
	std::vector<size_t> atot_qual(var_seq.size());
	// See if any of the vars are SV
	bool has_sv = false;
	for (size_t i = 0; i < vars.size(); i++) {
		if (vars[i].is_structural == true) {
			has_sv = true;
		}
	}
	// Give up if holes and SVs
	if (has_sv && per_assembly.has_holes) {
		return false;
	}

	if (pile) {
		for (size_t i = 0; i < vars.size(); i++) {
          if (vars[i].seq_begin == 0) {
            // TODO(nils): Why does this case happen?  Is something else wrong?
            if (g_zero_seq_begin.fetch_add(1) < 1000) {
              SPLOG("Seq begin is zero for sequence %s", var_seq.as_string().c_str());
            }
            return false;
          }
			size_t depth = pile->depth_at(vars[i].seq_begin - 1);
			for (int j = vars[i].seq_begin; j <= vars[i].seq_end; j++) {
				depth = std::min(depth, pile->depth_at(j));
			}
			if (depth > 0) {
				has_anything = true;
				break;
			}
		}
		if (!has_anything) {
			//pile->print();
          g_no_actual_depth.fetch_add(1);
			return false;
		}
		for (size_t i = 0; i < var_seq.size(); i++) {
			depths[i] = pile->depth_at(i);
			afwd[i] = pile->fwd_at(i);
			atot_qual[i] = pile->tot_qual_at(i);
		}
	}
	else {
		for (size_t i = 0; i < var_seq.size(); i++) {
			depths[i] = 0;
		}
	}

	int sub_id = 2 * per_assembly.var_id;
	for (size_t i = 0; i < vars.size(); i++) {
		const auto& var = vars[i];

		if (var.anchor_drop) {
			continue;
		}

		struct_var out = per_assembly;
		out.is_structural  = var.is_structural;
		out.align_failed   = var.align_failed;
		out.ref_start      = ref.get_seq_position(var.left_ref);
		out.rev_start      = var.left_ref.is_rev_comp();
		out.ref_end        = ref.get_seq_position(var.right_ref);
		out.rev_end        = var.right_ref.is_rev_comp();
		out.assembled      = var_seq;
		out.assembly_depth = depths;
		out.assembly_fwd   = afwd;
		out.assembly_tot_qual = atot_qual;
		out.var_start      = var.seq_begin;
		out.var_end        = var.seq_end;
		out.flipped        = false;
		out.sub_id         = sub_id;
		out.canonicalize();

		/*
		if (out.rev_start == false &&
			out.rev_end == false &&
			out.ref_start.scaffold_id == out.ref_end.scaffold_id &&
			out.ref_start != out.ref_end && // Hack to hide call_structural bug that returns equal reference ends.  TODO: Fix this.
			abs(
				int(out.ref_end.position - out.ref_start.position) -
				int(out.var_end - out.var_start)
			) < 200)
		{
		*/
		if (!out.is_structural) {
			if (out.flipped) {
				out.ref_seq = dna_sequence(var.right_ref.rev_comp() + 1, var.left_ref.rev_comp());
			}
			else {
				out.ref_seq = dna_sequence(var.left_ref + 1, var.right_ref);
			}
		}
		if (out.is_structural) {
			sub_id++;
		}
		if (pile) {
			size_t depth = pile->depth_at(out.var_start - 1);
			size_t tot_depth = depth;
			size_t tot_count = 1;
			for (size_t j = out.var_start; j <= out.var_end; j++) {
				tot_count++;
				tot_depth += pile->depth_at(j);
				depth = std::min(depth, pile->depth_at(j));
			}
			out.depth = depth;
			out.avg_depth = double(tot_depth) / double(tot_count);
		} else {
			out.depth = 0;
			out.avg_depth = 0;
		}
		if (out.depth >= min_depth) {
			// SPLOG("Writing variation: %s -> %s, depth = %d",
			// 	out.ref_seq.as_string().c_str(),
			// 	out.assembled.subseq(out.var_start, out.var_end - out.var_start).as_string().c_str(),
			// 	(int) out.depth);

			is_var = true;
			sink_f(out);
			if (dup_structural && out.is_structural) {
				out.flip();
				sink_f(out);
			}
		}
	}

	return is_var;
}

void log_struct_var_adapter_stats() {
  SPLOG("make_vars no actual depth: %ld", g_no_actual_depth.exchange(0));
  SPLOG("make_vars no variations must be exact match: %ld",
        g_no_variations_must_be_exact_match.exchange(0));
  SPLOG("make_vars seq_begin unexpectedly zero: %ld",
        g_zero_seq_begin.exchange(0));
}

make_vars::make_vars(const std::string& ref_name, size_t leeway, size_t read_length, const std::string& out, const std::string& scaf)
	: m_ref(ref_name)
	, m_gen(3)
	, m_leeway(leeway)
	, m_read_length(read_length)
	, m_next_id(0)
	, m_writer(out)
	, m_exporter(m_writer, ref_name)
	, m_scaffold(scaf)
{
	const reference_assembly& ref_assembly = m_ref.get_assembly();
	if (scaf == "") {
		m_scaffold_start = 0;
		m_scaffold_len = m_ref.size();
		return;
	}

	SPLOG("Looking for scaffold %s", scaf.c_str());
        auto sup_begin_it = ref_assembly.supercontigs.lower_bound(supercontig(scaf, 0, 0));
        auto sup_end_it = ref_assembly.supercontigs.upper_bound(supercontig(scaf, std::numeric_limits<std::size_t>::max(), 0));
	if (sup_begin_it == sup_end_it) {
		throw io_exception("Unknown scaffold: " + scaf);
	}
	m_scaffold_start = sup_begin_it->tot_offset;
	m_scaffold_len = 0;
	for(auto it = sup_begin_it; it != sup_end_it; ++it) {
		m_scaffold_len += it->len;
	}
	SPLOG("From at %zu, size = %zu", m_scaffold_start, m_scaffold_len);
}

size_t make_vars::random_loc(size_t space)
{
	int count = 0;
	while (true) {
		size_t off = m_scaffold_start + boost::uniform_int<size_t>(m_leeway, m_scaffold_len - space - m_leeway)(m_gen);
		//printf("Picking %d\n", (int) off);
		const supercontig& sc = m_ref.get_assembly().get_supercontig(off);
		if (sc.tot_offset + m_leeway > off || off + space > sc.tot_offset + sc.len - m_leeway) {
			count++;
			continue;
		}
		bool overlap = false;
		for (auto kvp : m_vars) {
			if (off + space + m_leeway > kvp.second.ref_start &&
				kvp.second.ref_end + m_leeway > off) {
				//printf("(%d, %d) too close to (%d, %d)\n",
				//	(int) off, (int) (off + space),
				//	(int) kvp.second.ref_start, (int) kvp.second.ref_end);
				overlap = true;
				break;
			}
		}
		if (!overlap)
			return off;
		count++;
		if (count == 1000) {
			throw io_exception(str(format("Giving up, can't find %d open bases") % space));
		}
	}
}

void make_vars::call(const var_info& vi)
{
	const int overlap = std::max((int) (5*m_read_length/4), 5);

	auto its = m_ref.get_dna(vi.ref_start);
	auto ite = m_ref.get_dna(vi.ref_end);

	auto new_seq = dna_sequence(its - overlap, its) + vi.var_seq + dna_sequence(ite, ite + overlap);
	its -= overlap;
	ite += overlap;

	//SPLOG("new_seq = %s", new_seq.as_string().c_str());
	//SPLOG("ref_seq = %s", dna_sequence(its, ite).as_string().c_str());

	dna_slice sbound = m_ref.get_supercontig(its - m_ref.get_dna(0));
	dna_slice ebound = m_ref.get_supercontig(ite - m_ref.get_dna(0));
	sbound = dna_slice(its, sbound.end());
    ebound = dna_slice(ebound.begin(), ite);

	struct_var svo;
	svo.var_id = m_next_id++;
	svo.is_ambig = false;
	svo.min_overlap = 0;
	svo.avg_overlap = 0;
	svo.has_holes = false;
	struct_var_adapter(m_ref, [&](const struct_var& var) {
        m_exporter.write_msgpack(var.ref_start,var);
      }, new_seq, sbound, ebound, nullptr, svo);
	m_vars.insert(std::make_pair(vi.ref_start, vi));
}

void make_vars::snp(const std::string& name)
{
	size_t loc = random_loc();
	dna_base base;
	do {
		base = dna_base(boost::uniform_int<int>(0, 3)(m_gen));
	} while (base == *m_ref.get_dna(loc));
	var_info vi;
	vi.name = name;
	vi.ref_start = loc;
	vi.ref_end = loc + 1;
	vi.ref_seq.push_back(*m_ref.get_dna(loc));
	vi.var_seq.push_back(base);
	call(vi);
}

void make_vars::random_insert(const std::string& name, size_t size)
{
	var_info vi;
	vi.name = name;
	vi.ref_start = random_loc();
	vi.ref_end = vi.ref_start;
	dna_base orig_start = *m_ref.get_dna(vi.ref_start);
	dna_base orig_end = *m_ref.get_dna(vi.ref_start - 1);
	do {
		vi.var_seq = dna_sequence();
		for (size_t i = 0; i < size; i++) {
			vi.var_seq.push_back(dna_base(boost::uniform_int<int>(0, 3)(m_gen)));
		}
	} while (vi.var_seq[0] == orig_start || vi.var_seq[size - 1] == orig_end);
	call(vi);
}

void make_vars::repeat_insert(const std::string& name, size_t size)
{
	var_info vi;
	vi.name = name;
	size_t loc;
	do {
		loc = random_loc();
	} while (*m_ref.get_dna(loc) == *m_ref.get_dna(loc + size - 1) ||
		*m_ref.get_dna(loc) == *m_ref.get_dna(loc + size) ||
		*m_ref.get_dna(loc + size - 1) == *m_ref.get_dna(loc - 1));
	vi.ref_start = loc;
	vi.ref_end = loc;
	vi.var_seq = dna_sequence(m_ref.get_dna(loc), m_ref.get_dna(loc + size));
	call(vi);
}

void make_vars::random_delete(const std::string& name, size_t size)
{
	var_info vi;
	size_t loc;
	do {
		loc = random_loc(size);
	} while (*m_ref.get_dna(loc) == *m_ref.get_dna(loc + size) ||
		*m_ref.get_dna(loc-1) == *m_ref.get_dna(loc + size -1 ));
	vi.name = name;
	vi.ref_start = loc;
	vi.ref_end = loc + size;
	vi.ref_seq = dna_sequence(m_ref.get_dna(loc), m_ref.get_dna(loc + size));
	call(vi);
}

void make_vars::transpose(const std::string& name, size_t size)
{
	var_info vi;
	size_t loc;
	do {
		loc = random_loc(size);
	} while (*m_ref.get_dna(loc) == (*m_ref.get_dna(loc + size - 1)).complement());
	vi.name = name;
	vi.ref_start = loc;
	vi.ref_end = loc + size;
	vi.ref_seq = dna_sequence(m_ref.get_dna(loc), m_ref.get_dna(loc + size));
	vi.var_seq = vi.ref_seq.rev_comp();

	if (size < m_read_length)
		call(vi);
	else {
		const int overlap = std::max((int) (5*m_read_length/4), 5);

		auto its = m_ref.get_dna(vi.ref_start);
		auto ite = m_ref.get_dna(vi.ref_end);
		auto join1 = dna_sequence(its - overlap, its) + dna_sequence(ite - overlap, ite).rev_comp();
		auto join2 = dna_sequence(its, its + overlap).rev_comp() + dna_sequence(ite, ite + overlap);
		dna_slice scaffold_bound = m_ref.get_supercontig(its - m_ref.get_dna(0));

		dna_slice sbound1 = scaffold_bound;
        sbound1 = dna_slice(its - overlap, sbound1.end());
		dna_slice sbound2 = scaffold_bound.rev_comp();
        sbound2 = dna_slice((its + overlap - 1).rev_comp(), sbound2.end());
		dna_slice ebound1 = scaffold_bound.rev_comp();
        ebound1 = dna_slice(ebound1.begin(), (ite - overlap -1).rev_comp());
		dna_slice ebound2 = scaffold_bound;
        ebound2 = dna_slice(ebound2.begin(), ite + overlap);

		struct_var svo;
		svo.var_id = m_next_id++;
		svo.is_ambig = false;
		svo.min_overlap = 0;
		svo.avg_overlap = 0;
		svo.has_holes = false;
        auto export_f = [&](const struct_var& var) {
          m_exporter.write_msgpack(var.ref_start, var);
        };
		struct_var_adapter(m_ref, export_f, join1, sbound1, ebound1, nullptr, svo);
		struct_var_adapter(m_ref, export_f, join2, sbound2, ebound2, nullptr, svo);
	}

	m_vars.insert(std::make_pair(vi.ref_start, vi));
}

void make_vars::print_sequence(FILE* out)
{
	const supercontig* sc = NULL;
	size_t in_line = 0;
	auto it = m_vars.begin();
	size_t i = 0;
	size_t pos = 0;
	std::string prev_scaffold;
	while (i < m_ref.size()) {
		if (it != m_vars.end() && i > it->first) {
			throw io_exception(str(format("WTF: %s %d %d %d %d")
				% it->second.name.c_str()
				% it->first
				% it->second.ref_start
				% it->second.ref_end
				% i
			));
		}
		if (it != m_vars.end() && (
			it->first != it->second.ref_start ||
			it->first > it->second.ref_end)) {
			throw io_exception(str(format("WTF: %s %d %d %d %d")
				% it->second.name.c_str()
				% it->first
				% it->second.ref_start
				% it->second.ref_end
				% i
			));
		}
		if (sc == NULL || i == sc->tot_offset + sc->len) {
			sc = &m_ref.get_assembly().get_supercontig(i);

			if (prev_scaffold != sc->scaffold_name) {
				SPLOG("sc->scaffold_name = %s, prev_scaffold = %s, i = %lu, pos = %lu", sc->scaffold_name.c_str(), prev_scaffold.c_str(), i, pos);
				prev_scaffold = sc->scaffold_name;
				pos = 0;

				if (in_line != 0) fprintf(out, "\n");
				if (m_scaffold == "" || (sc->scaffold_name == m_scaffold))
					fprintf(out, ">%s\n", sc->scaffold_name.c_str());
				in_line = 0;
			}

			if (m_scaffold == "" || (sc->scaffold_name == m_scaffold)) {
				while (pos < sc->offset) {
					fputc('N', out);
					in_line++;
					if (in_line == 80) { fprintf(out, "\n"); in_line = 0; }
					pos++;
				}
			}
		}
		if (it != m_vars.end() && it->first == i) {
			const var_info& vi = it->second;
			for (size_t j = 0; j < vi.var_seq.size(); j++) {
				fprintf(out, "%c", (char) vi.var_seq[j]);
				in_line++;
				if (in_line == 80) { fprintf(out, "\n"); in_line = 0; }
			}
			pos += vi.ref_end - i;
			i = vi.ref_end;
			++it;
		}
		else {
			if (m_scaffold == "" || (sc->scaffold_name == m_scaffold)) {
				fprintf(out, "%c", (char) *m_ref.get_dna(i));
				in_line++;
				if (in_line == 80) { fprintf(out, "\n"); in_line = 0; }
			}
			i++;
			pos++;
		}
	}
}
