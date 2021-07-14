#pragma once

#include "modules/bio_base/seqset.h"
#include "modules/bio_base/seqset_anchor.h"
#include <unordered_map>

struct assembly
{
	assembly flip() const {
		struct assembly out;
		out.left = right.rev_comp();
		out.right = left.rev_comp();
		out.assembly = assembly.rev_comp();
		out.depth = depth;
		std::reverse(out.depth.begin(), out.depth.end());
		out.min_overlap = min_overlap;
		out.id = id;
		return out;
	}
	// dna iterator to reference left == assembled[0]
	dna_const_iterator left;
	// dna iterator to reference right == assembled[assembled.size() - 1]
	dna_const_iterator right;
	// The assembled sequence
	dna_sequence assembly;
	// The depth of each assembly position
	std::vector<uint8_t> depth;
	// The minimum overlap across the assembly
	uint8_t min_overlap;
	// A unique ID used by python
	uint32_t id;
};

template<typename OutIt, typename InCol>
bool seqset_assemble(OutIt out, const seqset& _seqset, const InCol& left, const InCol& right, 
		uint8_t min_overlap, uint32_t max_ops, bool skip_ambig,
		const seqset_bitmap_base& bitmap = seqset_bitmap_true())
{
	uint8_t read_len = _seqset.read_len();
	struct entry {
		uint64_t prev;
		uint8_t  tot_overlap;
		uint8_t  cur_overlap;
		bool     is_anchor;
		//branch points from original anchor
		uint32_t  branches;
	};
	//
	std::multimap<uint8_t, uint64_t> todo;
	std::unordered_map<uint64_t, entry> found;
	std::unordered_map<uint64_t, uint32_t> terms;
	std::vector<dna_const_iterator> refs;
	for(const anchor& a : left) {
		todo.emplace(0, a.entry);
		
		entry& e = found[a.entry];
		e.prev = refs.size();
		refs.push_back(a.ref_pos);
		e.tot_overlap = a.overlap;
		e.cur_overlap = a.overlap;
		e.is_anchor = true;
		e.branches = 0;
	}
	for(const anchor& a : right) {
		dna_sequence seq = _seqset.ctx_entry(a.entry).sequence();
		seqset_range rev_ctx = _seqset.find(seq.rev_comp());
		terms[rev_ctx.begin()] = refs.size();
		refs.push_back(a.ref_pos);
	}
	
	uint32_t op_count = 0;
	while(todo.size() && op_count < max_ops) {
		uint64_t src = todo.begin()->second;
		const entry& e = found[src];
		seqset_range c = _seqset.ctx_entry(src);
		todo.erase(todo.begin());

		overlaps_t results;
		//max_ops smaller - number of reads to return - single read to N 10 branch - skip it
		bool r = c.find_overlap_reads(results, max_ops - op_count, min_overlap, bitmap); //false on space

		if (!r) return false; 
		
		op_count += results.size();
		//false
		if (results.size() > 1 && skip_ambig) {
			continue;
		}
		
		for(const auto& kvp : results) {
			//SPLOG("   : %s", n.sequence().as_string().c_str());
			if (found.count(kvp.first)) {
				continue;
			}
			entry& ne = found[kvp.first];
			ne.prev = src;
			ne.cur_overlap = kvp.second;
			ne.tot_overlap = std::min(e.tot_overlap, ne.cur_overlap);
			ne.is_anchor = false;
			ne.branches += results.size() > 1 ? 1 : 0;
			if (terms.count(kvp.first)) {
				continue;
			}
			todo.emplace(ne.branches, kvp.first);
		}
	}
	
	uint32_t id = 0;
	for(const auto& kvp : terms) 
	{
		auto it = found.find(kvp.first);
		if (it == found.end()) {
			continue;
		}
		assembly r;
		r.right = refs[kvp.second].rev_comp();
		r.min_overlap = it->second.tot_overlap;
		std::vector<uint8_t> read_start;
		while(!it->second.is_anchor) {
			int uniq_size = read_len - it->second.cur_overlap;
			r.assembly += _seqset.ctx_entry(it->first).sequence(uniq_size);
			it = found.find(it->second.prev);
			read_start.push_back(1);
			for(int i = 0; i < uniq_size - 1; i++) {
				read_start.push_back(0);
			}
		}
		r.assembly += _seqset.ctx_entry(it->first).sequence();
		r.left = refs[it->second.prev];
		read_start.push_back(1);
		read_start.resize(r.assembly.size());
		r.depth.resize(r.assembly.size());
		int cur_depth = 0;
		for(size_t i = 0; i < r.assembly.size(); i++) {
			if (read_start[i]) cur_depth++;
			if (i >= read_len && read_start[i - read_len]) cur_depth--;
			r.depth[i] = cur_depth;
		}
		r.assembly = r.assembly.rev_comp();
		std::reverse(r.depth.begin(), r.depth.end());
		r.id = id++;
		*out++ = r;
	}
	return todo.size() == 0;
}

