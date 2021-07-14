
#include "modules/bio_mapred/display_var.h"
#include "modules/bio_base/call_structural.h"
#include <map>

static bool check_depth(const std::vector<struct_var>& vars, unsigned int min_depth);

static bool svar_less(const struct_var& a, const struct_var& b) {
	return a.var_start < b.var_start;
}

const static size_t ref_cont = 50;
void display_var(color_text_buffer& out,
	const reference& ref,
	const std::vector<struct_var>& _svars,
	const std::vector<read_support>& _reads)
{
	uint32_t highlight_color = color::green;

	std::vector<struct_var> svars = _svars;
	// Make 'flipped' state all the same
	bool shared_flip = svars[0].flipped;
	dna_sequence seq = svars[0].assembled;
	for(size_t i = 0; i < svars.size(); i++) {
		if (svars[i].flipped != shared_flip)
			svars[i].flip();
	}
	// Sort by order of assembled seq start
	std::sort(svars.begin(), svars.end(), svar_less);

	// Convert struct_vars back into sv_outs
	std::vector<sv_out> vars;
	const reference_assembly& ref_assembly = ref.get_assembly();
	for(size_t i = 0; i < svars.size(); i++)
	{
		sv_out sv;
		sv.is_structural = svars[i].is_structural;
		sv.seq_begin = svars[i].var_start;
		sv.seq_end = svars[i].var_end;
		size_t left_flat = ref_assembly.flatten(svars[i].ref_start);
		size_t right_flat = ref_assembly.flatten(svars[i].ref_end);
		sv.left_ref = ref.get_dna(left_flat);
		sv.right_ref = ref.get_dna(right_flat);
		if (svars[i].rev_start) sv.left_ref = sv.left_ref.rev_comp();
		if (svars[i].rev_end) sv.right_ref = sv.right_ref.rev_comp();
		vars.push_back(sv);
	}
	// Add a sentinal
	sv_out sentinal;
	sentinal.is_structural = 0;
	sentinal.seq_begin = seq.size();
	sentinal.seq_end = seq.size();
	vars.push_back(sentinal);

	int cur_var = 0;
	int which_ref = -1;
	int seq_pos = 0;
	int char_pos = 0;
	dna_const_iterator seqit = seq.begin();
	dna_const_iterator refit = vars[0].left_ref + 1 - vars[0].seq_begin;
	// Print ref1 name stuff
	seq_position ref_start = ref.get_seq_position(refit);
	std::string name1 = printstring("%s:%ld(%s) -> ",
		ref_assembly.scaffold_order[ref_start.scaffold_id].c_str(), ref_start.position + 1,
		refit.is_rev_comp() ? "rev" : "fwd");
	out.set_color(color::white);
	out.set_position(char_pos -int(name1.size()), which_ref);
	out.print("%s", name1.c_str());
	out.set_position(char_pos - strlen("assembled -> "), 0);
	out.print("assembled -> ");
	std::map<int, int> pos_to_char;
	while(true) {
		const sv_out& cv = vars[cur_var];

		// Print matching section
		size_t match_len = cv.seq_begin - seq_pos;
		out.set_color(color::white);
		out.set_position(char_pos, 0);
		out.print("%s", dna_sequence(seqit, seqit + match_len).as_string().c_str());
		out.set_position(char_pos, which_ref);
		out.print("%s", dna_sequence(refit, refit + match_len).as_string().c_str());
		for(size_t i = 0; i < match_len; i++)
			pos_to_char[seq_pos + i] = char_pos + i;
		seqit += match_len; refit += match_len; char_pos += match_len; seq_pos += match_len;

		if (cv.seq_begin == (int) seq.size()) break;

		// Print variation section
		out.set_color(highlight_color);
		out.set_position(char_pos, 0);
		out.print("%s", dna_sequence(seqit, seqit + cv.seq_end - cv.seq_begin).as_string().c_str());
		int seq_len = int(cv.seq_end - cv.seq_begin);
		for(size_t i = 0; i < size_t(seq_len); i++)
			pos_to_char[seq_pos + i] = char_pos + i;
		seqit += seq_len;
		seq_pos += seq_len;
		if (!cv.is_structural) {
			// Print ref region
			out.set_position(char_pos, which_ref);
			out.print("%s", dna_sequence(refit, cv.right_ref).as_string().c_str());
			int ref_len = (int) (cv.right_ref - cv.left_ref - 1);
			int max_len = std::max(seq_len, ref_len);
			// Draw grey regions as needed
			out.set_color(color::grey);
			out.set_position(char_pos + seq_len, 0);
			out.print("%s", std::string(max_len - seq_len, '.').c_str());
			out.set_position(char_pos + ref_len, which_ref);
			out.print("%s", std::string(max_len - ref_len, '.').c_str());
			char_pos += max_len;
		} else {
			// Print broken ends and switch references
			out.set_position(char_pos, which_ref);
			out.set_color(highlight_color);
			// TODO: 20 base pairs may go off the end of the supercontig
			out.print("%s...", dna_sequence(refit, refit + ref_cont).as_string().c_str());
			char_pos += cv.seq_end - cv.seq_begin;
			which_ref = 1;
			// Also 20 base problem
			dna_const_iterator refit2 = cv.right_ref - ref_cont;
			out.set_position(char_pos - (ref_cont + 3), which_ref);
			out.print("...%s", dna_sequence(refit2, refit2 + ref_cont).as_string().c_str());
			// Get ref name
			seq_position ref_start_2 = ref.get_seq_position(refit2);
			std::string name2 = printstring("%s:%ld(%s) - ",
				ref_assembly.scaffold_order[ref_start_2.scaffold_id].c_str(), ref_start_2.position + 1,
				refit2.is_rev_comp() ? "rev" : "fwd");
			out.set_color(color::white);
			out.set_position(char_pos -int(name2.size()) - (ref_cont + 3), which_ref);
			out.print("%s", name2.c_str());
		}
		refit = cv.right_ref;
		cur_var++;
	}
	// Map all the reads
	//pileup p(seq, 5000);
	std::multimap<int, read_support> reads;
	for(read_support rs : _reads) {
		if (rs.flipped != shared_flip)
		{
			rs.flip();
			rs.name += " rev";
		}
		int pos = rs.pos;
		if (shared_flip) pos = seq.size() - pos - rs.original.size();
		//p.add_read(rs.name, rs.original, rs.quality, flipped == shared_flip, pos);
		reads.insert(std::make_pair(pos, rs));
	}

	// Write them out in order
	int read_num = 3;
	for(const auto& kvp : reads)
	{
		int pos = kvp.first;
		const read_support& r = kvp.second;
		out.set_position(pos_to_char[pos] -int(r.name.size()) - 2, read_num);
		out.set_color(color::white);
		out.print("%s: ", r.name.c_str());
		int exp_char_pos = pos_to_char[pos];
		for(size_t i = 0; i < r.original.size(); i++) {
			int char_pos = pos_to_char[pos + i];
			while(exp_char_pos < char_pos) {
				out.set_color(color::grey);
				out.set_position(exp_char_pos, read_num);
				out.print(".");
				exp_char_pos++;
			}
			if (r.original[i] == seq[pos + i])
				out.set_color(color::white);
			else
				out.set_color(highlight_color);
			out.set_position(pos_to_char[pos + i], read_num);
			out.print("%c", (char) r.original[i]);
			exp_char_pos++;
		}
		read_num++;
	}
}

void display_vars(const path& out_path, const reference& ref, kv_source& svs, kv_source& rds, unsigned int min_depth)
{
	struct_var_key k1;
	struct_var_key k2;
	read_support rs;
	struct_var v;
	std::vector<struct_var> vars;
	std::vector<read_support> reads;

	uint32_t var_id = 0;
	bool var_valid = svs.read_msgpack(k1, v);
	bool read_valid= rds.read_msgpack(k2, rs);
	SPLOG("Reading loop");
	while(read_valid || var_valid || vars.size() > 0) {
		// If lookahead var is valid, add it and go around again
		if (var_valid && k1.variation_id == var_id) {
			vars.push_back(v);
			var_valid = svs.read_msgpack(k1, v);
			continue;
		}
		// If lookahead read is valid, add it and go around again
		if (read_valid && k2.variation_id == var_id) {
			if (reads.size() % 1000 == 0)
				SPLOG("Read: %d", (int) reads.size());
			reads.push_back(rs);
			read_valid = rds.read_msgpack(k2, rs);
			continue;
		}
		SPLOG("Done with var_id %d", (int) var_id);
		// Ok, do visualization if not empty
		if (check_depth(vars, min_depth))
		{
			color_text_buffer ctb;
			display_var(ctb, ref, vars, reads);
			path fpath = out_path.append(printstring("%d.html", var_id));
			std::unique_ptr<writable> out = fpath.write();
			out->print("<html><body>");
			out->print("<hr/>");
			ctb.render_as_html(*out);
			out->print("<hr/>");
			out->print("</body></html>");
			out->close();
		}
		// Clear buffers and move to next var_id
		vars.clear();
		reads.clear();
		if (var_valid) var_id = k1.variation_id;
		if ((! var_valid) || (read_valid && k2.variation_id < var_id)) var_id = k2.variation_id;
	}
}

static inline bool check_depth(const std::vector<struct_var>& vars, unsigned int min_depth)
{
	return std::any_of(vars.cbegin(), vars.cend(),
		[min_depth](const struct_var& the_sv) { return the_sv.depth >= min_depth; }
	);
}
