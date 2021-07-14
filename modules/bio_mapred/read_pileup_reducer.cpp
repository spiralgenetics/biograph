
#include "base/base.h"
#include "modules/bio_mapred/read_pileup_reducer.h"
#include "modules/io/make_unique.h"
#include "modules/bio_format/make_vars.h"
#include <limits>

REGISTER_1(reducer, read_pileup, const std::string&);

void read_pileup_params::validate(){
	SPLOG_P(LOG_DEBUG, "read_pileup_params::validate> reference: %s, min_depth: %d, var_infos: %lu",
		reference.c_str(), min_depth, var_infos.get_size());
}

read_pileup_reducer::read_pileup_reducer(const std::string& params)
	: m_query_started(false)
{
	json_deserialize(m_params, params);
	m_params.validate();
}

void read_pileup_reducer::typed_start(const struct_var_key& key)
{
	//SPLOG("Doing, start, var id = %d", key.variation_id);
	if (!m_query_started) {
		m_ref = make_unique<reference>(m_params.reference);
		struct_var_key start = key;
		start.read_id = 0;
		struct_var_key end(std::numeric_limits<uint32_t>::max(), std::numeric_limits<uint32_t>::max());
		m_query.find_msgpack(m_params.var_infos, start, end);
		m_query_started = true;
	}
	struct_var_key kout;
	std::pair<dna_sequence, var_info> val;
	m_query.read_msgpack(kout, val);
	m_sequence = val.first;
	m_var_info = val.second;
	if (kout.variation_id != key.variation_id) {
		throw io_exception(printstring("var_info not synced with read support: %u, %u", kout.variation_id, key.variation_id));
	}
	m_pileup = make_unique<pileup>(m_sequence, 0);
	m_var_id = key.variation_id;
}

void read_pileup_reducer::typed_add_value(const struct_var_key& key, const read_support& support)
{
	dna_sequence seq = support.corrected;
	std::string qual = support.quality;
        if (support.flipped) {
		seq = seq.rev_comp();
		std::reverse(qual.begin(), qual.end());
	}
	//SPLOG("Doing add_value, seq = %s, pos = %lu", seq.as_string().c_str(), support.pos);
	m_pileup->add_read(support.name, seq, qual, !support.flipped, support.pos);
}

void read_pileup_reducer::typed_end()
{
	//SPLOG("Doing end");
	CHECK(m_ref);
	CHECK(m_pileup);
	//m_pileup->print();
	bool has_holes = false;
	for (size_t i = 0; i < m_sequence.size(); i++) {
		if (m_pileup->depth_at(i) == 0) has_holes = true;
	}

	//DEBUG_SPLOG("Computing sbound and ebound");
	dna_const_iterator sit = m_ref->get_dna(m_var_info.s_ref);
	dna_slice sbound = m_ref->get_supercontig(m_var_info.s_ref);
	dna_const_iterator eit = m_ref->get_dna(m_var_info.e_ref);
	dna_slice ebound = m_ref->get_supercontig(m_var_info.e_ref);
	if (m_var_info.s_flip) {
		//DEBUG_SPLOG("Flipping sbound");
		sit = sit.rev_comp();
		sbound = sbound.rev_comp();
	}
	if (m_var_info.e_flip) {
		//DEBUG_SPLOG("Flipping ebound");
		eit = eit.rev_comp();
		ebound = ebound.rev_comp();
	}
    sbound = dna_slice(sit, sbound.end());
    ebound = dna_slice(ebound.begin(), eit + 1);

	if (has_holes) {
		if (sit.is_rev_comp() != eit.is_rev_comp()) {
			SPLOG("read_graph_2::output_structural> Dropping %u due to holes and mismatching anchor complements", m_var_id);
			return;
		}
		else if (std::abs(int64_t(m_sequence.size()) - std::abs(eit - sit)) > 20)
		{
			SPLOG("read_graph_2::output_structural> Dropping due to holes with anchors %ld bases apart while assembly has %lu bases",
				std::abs(eit - sit), m_sequence.size());
			return;
		}
		else; // It's probably a compound heterozygote so we keep it and fall through.
	}

	struct_var base;
	base.var_id = m_var_id;
	base.is_ambig = m_var_info.is_ambig;
	base.min_overlap = m_var_info.min_overlap;
	base.avg_overlap = m_var_info.avg_overlap;
	base.has_holes = has_holes;
	struct_var_adapter(*m_ref,
                       [&](const struct_var& var) {
                         m_out_context->write_msgpack(var.ref_start,var);
                       },
                       m_sequence, sbound, ebound, m_pileup.get(), base,
                       m_params.min_depth, true);
}


