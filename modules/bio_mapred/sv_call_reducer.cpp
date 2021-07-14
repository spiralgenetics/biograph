
#include "base/base.h"
#include "modules/bio_mapred/sv_call_reducer.h"
#include "modules/io/json_transfer.h"
#include "modules/io/make_unique.h"
#include <limits>

REGISTER_1(reducer, sv_call, const std::string&);

sv_call_reducer::sv_call_reducer(const std::string& params)
	: m_params(inline_json_deserialize<sv_call_params>(params))
	, m_ref(m_params.reference)
{
	//SPLOG("Constructing SV_CALL");
}

void sv_call_reducer::typed_start(const seq_position& key)
{
	//SPLOG("Got a start");
	if (m_bc_uniq == nullptr) {
		resource_manager rm;
		rm.read_resource(m_buf, m_params.coverage, [&](double) {
			m_watchdog();
		});
		size_t bit_size = bitcount::compute_size(m_ref.size());
		m_bc_uniq = make_unique<bitcount>(m_buf.buffer(), m_ref.size());
		m_bc_guess = make_unique<bitcount>(m_buf.buffer() + bit_size, m_ref.size());
	}
}

void sv_call_reducer::typed_add_value(const seq_position& key, const struct_var& value)
{
	if (m_pos_end.scaffold_id != key.scaffold_id ||
		m_pos_end.position <= key.position) {
		dump_current();
	}

	if (value.is_structural) {
		//SPLOG("Forwarding structural");
		sv_call svc;
		svc.sources.push_back(value);
		svc.position = key;
		svc.position.position++;
		size_t start = m_ref.flatten(value.ref_start);
		size_t end = m_ref.flatten(value.ref_end);
		if (start > end) std::swap(start, end);
		if (start == end) { end++; }
		if (m_ref.get_supercontig(start).begin() == m_ref.get_supercontig(end).begin()) {
          size_t super_start = m_ref.get_supercontig(start).begin() - m_ref.get_dna(0);
			size_t back_up = (start >= m_params.read_size ? start - m_params.read_size : 0);
			size_t final_start = std::max(super_start, back_up);
			size_t tot_coverage = 0;
			for(size_t i = final_start; i <= end; i++) {
				if (m_bc_uniq->get(i) || m_bc_guess->get(i)) {
					size_t os = std::max(i, start);
					size_t oe = std::min(i + m_params.read_size, end);
					tot_coverage += (oe - os);
				}
			}
			svc.sv_ref_depth = double(tot_coverage) / double(end - start);
		} else {
			svc.sv_ref_depth = -1;
		}

		output(svc.position, svc);
		return;
	}

	if (m_pos_end.scaffold_id != key.scaffold_id ||
		m_pos_end.position <= key.position) {
		m_pos_start = key;
	}

	if (value.ref_end > m_pos_end)
		m_pos_end = value.ref_end;
	m_current.sources.push_back(value);
}

void sv_call_reducer::typed_end()
{
	dump_current();
}

void sv_call_reducer::dump_current()
{
	if (m_current.sources.size() == 0) return;
	m_current.sv_ref_depth = -1;

	// Get data regarding the region of interest
	size_t start = m_ref.flatten(m_pos_start);
	size_t end = m_ref.flatten(m_pos_end);
	auto range = m_ref.get_supercontig(start);

	// Compute coverage data for the region
	size_t super_start = range.begin() - m_ref.get_dna(0);
	size_t back_up = (start >= m_params.read_size ? start - m_params.read_size : 0);
	size_t final_start = std::max(super_start, back_up);
	std::map<size_t, size_t> quick_pileup;
	for(size_t i = final_start; i <= end; i++) {
		if (m_bc_uniq->get(i) || m_bc_guess->get(i)) {
			for(size_t j = i; j < i + m_params.read_size; j++)
				quick_pileup[j]++;
		}
	}

	m_current.position = m_pos_start;
	m_current.position.position++;

	// Generate the reference allele
	allele ref_allele;
	ref_allele.seq = dna_sequence(m_ref.get_dna(start + 1), m_ref.get_dna(end));
	ref_allele.depth.resize(end - start + 1);
	for(size_t i = start; i <= end; i++) {
		ref_allele.depth[i - start] = quick_pileup[i];
	}
	m_current.alleles.push_back(ref_allele);

	// Group struct_vars by var_id and sort by position
	std::map<int, std::map<seq_position, struct_var> > by_var_id;
	for(const auto& v : m_current.sources) {
		by_var_id[v.sub_id][v.ref_start] = v;
	}

	std::map<dna_sequence, allele> uniq_alleles;

	// Make an allele for each unique sequence
	for(const auto& kvp : by_var_id) {
		const auto& by_pos = kvp.second;
		const dna_sequence& assembled = by_pos.begin()->second.assembled;
		const std::vector<size_t>& pileup = by_pos.begin()->second.assembly_depth;
		const std::vector<size_t>& fwd = by_pos.begin()->second.assembly_fwd;
		const std::vector<size_t>& tot_qual = by_pos.begin()->second.assembly_tot_qual;
		size_t tot_ref_start = m_ref.flatten(by_pos.begin()->second.ref_start);
		size_t tot_ref_end = m_ref.flatten(by_pos.rbegin()->second.ref_end);
		size_t tot_assembly_start = by_pos.begin()->second.var_start;
		size_t tot_assembly_end = by_pos.rbegin()->second.var_end;
		CHECK_LE(start + 1, tot_ref_start + 1);
		if (tot_assembly_start > tot_assembly_end) {
			for(const auto& kvp2 : by_pos) {
				SPLOG("%s", json_serialize(kvp2.second).c_str());
			}
		}
		CHECK_LE(tot_assembly_start, tot_assembly_end);
		CHECK_LE(tot_ref_end, end);
		dna_sequence allele_seq =
			dna_sequence(m_ref.get_dna(start + 1), m_ref.get_dna(tot_ref_start + 1)) +
			dna_sequence(assembled.begin() + tot_assembly_start, assembled.begin() + tot_assembly_end) +
			dna_sequence(m_ref.get_dna(tot_ref_end), m_ref.get_dna(end));
		allele& the_allele = uniq_alleles[allele_seq];
		the_allele.seq = allele_seq;
		the_allele.sub_ids.push_back(kvp.first);
		the_allele.depth.resize(tot_assembly_end - tot_assembly_start + 2);
		the_allele.fwd.resize(tot_assembly_end - tot_assembly_start + 2);
		the_allele.tot_qual.resize(tot_assembly_end - tot_assembly_start + 2);
		for(size_t i = tot_assembly_start - 1; i <= tot_assembly_end; i++) {
			the_allele.depth[i - (tot_assembly_start - 1)] += pileup[i];
			the_allele.fwd[i - (tot_assembly_start - 1)] += fwd[i];
			the_allele.tot_qual[i - (tot_assembly_start - 1)] += tot_qual[i];
		}
	}

	// Output all the alleles to current
	for(const auto& kvp : uniq_alleles) {
		m_current.alleles.push_back(kvp.second);
	}

	// Output results and clear
	output(m_pos_start, m_current);
	for(size_t i = 0; i < m_current.alleles.size(); i++) {
		const allele& the_allele = m_current.alleles[i];
		size_t min_depth = 100000;
		double avg_depth = 0.0;
		for(size_t i = 0; i < the_allele.depth.size(); i++) {
			min_depth = std::min(min_depth, the_allele.depth[i]);
			avg_depth += the_allele.depth[i];
		}
		avg_depth /= the_allele.depth.size();

		//SPLOG("%d: %s (%d, %d)", (int) i, the_allele.seq.as_string().c_str(),
		//	(int) min_depth, (int) avg_depth);
	}
	m_current.alleles.resize(0);
	m_current.sources.resize(0);
}
