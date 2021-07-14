#pragma once

#include "base/base.h"
#include "modules/bio_base/seqset.h"

struct anchor
{
	anchor(dna_const_iterator _ref_pos, uint64_t _entry, uint8_t _overlap)
		: ref_pos(_ref_pos)
		, entry(_entry)
		, overlap(_overlap)
	{}
	// The location of the last base pair of 'entry' in ref seq 
	dna_const_iterator ref_pos;
	// The entry, with the non-ref region on the front
	uint64_t entry;
	// The overlap between then end of entry the the ref
	uint8_t overlap;
};

template<class OutIt>
bool seqset_anchor(OutIt out, const seqset& _seqset, dna_slice ref, 
		uint8_t min_overlap, uint32_t max_anchors, 
		const seqset_bitmap_base& bitmap = seqset_bitmap_true())
{
	seqset_range c = _seqset.ctx_begin();
	size_t count = 0;
	for(auto it = ref.begin(); it != ref.end(); ++it) {
		//CHECK(c.size() != read_len);
		dna_base ref_comp = it->complement();
		for(int b = 1; b < 4; b++) {
			seqset_range n = c.push_front_drop(dna_base((int(ref_comp) + b) % 4), min_overlap); 
			if (!n.valid()) continue;
			overlaps_t r;
			bool good = n.find_overlap_reads(r, max_anchors - count, min_overlap + 1, bitmap, false /* don't rely only on bitmap */, 1);
			for(const auto& kvp : r) {
				dna_const_iterator ref_pos = (it - kvp.second);
				*out++ = anchor(ref_pos, kvp.first, kvp.second);
				count++;
			}
			if (!good) { return false; }
		}
		c = c.push_front_drop(ref_comp);
	}
	return true;
}

