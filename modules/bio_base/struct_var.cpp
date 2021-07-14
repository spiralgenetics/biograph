
#include "modules/bio_base/struct_var.h"

void read_support::flip()
{
	original = original.rev_comp();
	//corrected = corrected.complement();
	std::reverse(quality.begin(), quality.end());
	flipped = !flipped;
}

void struct_var::flip() 
{
	std::swap(ref_start, ref_end);
	std::swap(rev_start, rev_end);
	rev_start = !rev_start;
	rev_end = !rev_end;
	assembled = assembled.rev_comp();
	ref_seq = ref_seq.rev_comp();
	size_t old_start = var_start;
	var_start = assembled.size() - var_end;
	var_end = assembled.size() - old_start;
	flipped = !flipped;	
}

void struct_var::canonicalize() 
{
	if (ref_start < ref_end) return;
	flip();
}

bool safe_range(size_t pos, size_t tot_size) 
{
	return pos > 100 && pos + 100 < tot_size;
}
