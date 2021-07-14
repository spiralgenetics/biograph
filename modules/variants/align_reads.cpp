#include "modules/variants/align_reads.h"

namespace variants {

namespace {

constexpr bool k_dbg = false;

}  // namespace

align_reads::~align_reads() { flush(); }

void align_reads::on_assembly_edges(optional_aoffset cur_pos,
                                    const std::vector<assembly_ptr>& left_edges,
                                    const std::vector<assembly_ptr>& inserts,
                                    const std::vector<assembly_ptr>& right_edges) {
  if (k_dbg) {
    std::cerr << "Edges at " << cur_pos << "\n";
  }
  propagate_t left_in;
  if (!m_propagate.empty()) {
    auto next_prop_it = m_propagate.begin();
    CHECK_GE(next_prop_it->first, cur_pos);
    if (next_prop_it->first == cur_pos) {
      left_in = std::move(next_prop_it->second);
      m_propagate.erase(next_prop_it);
    }
  }

  for (const auto& a : left_edges) {
    find_starts(a, &left_in);
  }

  propagate_t inserts_out;
  for (const auto& a : inserts) {
    find_starts(a, &inserts_out);
    propagate(a, left_in, &inserts_out);
  }
  inserts_out.merge(std::move(left_in));

  for (const auto& a : right_edges) {
    auto& right_out = m_propagate[a->right_offset];
    propagate(a, inserts_out, &right_out);
  }
}

void align_reads::propagate(const assembly_ptr& a, const propagate_t& prop_in,
                            propagate_t* prop_out) {
  if (k_dbg) {
    std::cerr << "Propagate through " << *a << "\n";
  }
  CHECK(a->read_coverage);
  const auto& cov = *a->read_coverage;
  for (const auto& in_elem : prop_in) {
    read_trace_key key = in_elem.first;
    const read_trace& read = in_elem.second;

    int read_len = read.read_len_tot;
    int offset = key.read_len_left - read_len;

    key.read_ids &= cov.get_read_ids_at(offset, read_len);
    if (key.read_ids.empty()) {
      continue;
    }

    if (k_dbg) {
      std::cerr << "Propagating through " << in_elem.first.read_ids << " -> " << key.read_ids
                << ": " << read << " in " << *a << "\n";
    }

    if (a->matches_reference) {
      propagate_ref(a, 0, std::move(key), read, prop_out);
    } else {
      propagate_var(a, 0, std::move(key), read, prop_out);
    }
  }
}

void align_reads::find_starts(const assembly_ptr& a, propagate_t* prop_out) {
  CHECK(a->read_coverage);
  for (const auto& cov_entry : a->read_coverage->reads()) {
    if (cov_entry.offset < 0) {
      continue;
    }

    read_trace_key new_key;
    new_key.read_ids = cov_entry.read_ids;
    new_key.read_len_left = cov_entry.read_len;

    if (k_dbg) {
      std::cerr << "New start: " << new_key << "\n";
    }

    if (a->matches_reference) {
      read_trace read;
      read.read_len_tot = cov_entry.read_len;
      read.left_offset = a->left_offset + cov_entry.offset;

      propagate_ref(a, cov_entry.offset, std::move(new_key), std::move(read), prop_out);
    } else {
      read_trace read;
      read.read_len_tot = cov_entry.read_len;
      read.left_offset = a->left_offset;
      if (cov_entry.offset > 0) {
        if (m_refskip_anchor && read.left_offset) {
          --read.left_offset;
          add_cigar(read, cigar_op::REF_SKIP, 1);
        }
        add_cigar(read, cigar_op::PAD, cov_entry.offset);
      }

      propagate_var(a, cov_entry.offset, std::move(new_key), std::move(read), prop_out);
    }
  }
}

void align_reads::propagate_ref(const assembly_ptr& a, aoffset_t offset, read_trace_key key,
                                read_trace read, propagate_t* prop_out) {
  aoffset_t end_offset = offset + key.read_len_left;
  aoffset_t ref_bases;

  if (end_offset > aoffset_t(a->seq.size())) {
    ref_bases = aoffset_t(a->seq.size()) - offset;
  } else {
    ref_bases = end_offset - offset;
  }
  CHECK_GT(ref_bases, 0);
  add_cigar(read, cigar_op::MATCH, ref_bases);
  read.seq += dna_slice(a->seq).subseq(offset, ref_bases);
  key.read_len_left -= ref_bases;
  CHECK_GE(key.read_len_left, 0);
  if (key.read_len_left) {
    // Propagate to next
    if (k_dbg) {
      std::cerr << "Ref to next: " << key << ", " << read << "\n";
    }
    prop_out->emplace(std::move(key), std::move(read));
  } else {
    // Ends here.
    read.right_offset = a->left_offset + end_offset;
    output_aligned(std::move(key), std::move(read));
  }
}

void align_reads::propagate_var(const assembly_ptr& a, aoffset_t offset, read_trace_key key,
                                read_trace read, propagate_t* prop_out) {
  aoffset_t end_offset = offset + key.read_len_left;
  aoffset_t var_bases;
  if (end_offset > aoffset_t(a->seq.size())) {
    var_bases = aoffset_t(a->seq.size()) - offset;
  } else {
    var_bases = end_offset - offset;
  }
  CHECK_GE(var_bases, 0);
  if (a->seq.size() > 0) {
    CHECK_GT(var_bases, 0);
  }

  aoffset_t ref_bases = a->right_offset - a->left_offset;
  if (ref_bases == var_bases) {
    add_cigar(read, cigar_op::MATCH, ref_bases);
  } else {
    add_cigar(read, cigar_op::DELETE, ref_bases);
    add_cigar(read, cigar_op::INSERT, var_bases);
  }

  CHECK_GE(key.read_len_left, var_bases);
  key.read_len_left -= var_bases;
  read.seq += dna_slice(a->seq).subseq(offset, var_bases);
  if (key.read_len_left) {
    // propagate to next
    if (k_dbg) {
      std::cerr << "Var to next: " << key << ", " << read << "\n";
    }
    prop_out->emplace(std::move(key), std::move(read));
  } else {
    // Save on output.
    read.right_offset = a->right_offset;

    aoffset_t pad = aoffset_t(a->seq.size()) - end_offset;
    CHECK_GE(pad, 0);
    add_cigar(read, cigar_op::PAD, pad);
    output_aligned(std::move(key), std::move(read));
  }
}

void align_reads::add_cigar(read_trace& read, cigar_op op, int num_bases) {
  if (k_dbg) {
    std::cerr << "Adding cigar to " << read << ": " << char(op) << ", " << num_bases << "\n";
  }
  CHECK_GE(num_bases, 0);
  if (!num_bases) {
    return;
  }
  if (op != read.cur_cigar_op) {
    flush_cigar(read);
  }
  read.cur_cigar_op = op;
  read.cur_cigar_count += num_bases;
}

void align_reads::flush_cigar(read_trace& read) {
  if (read.cur_cigar_count) {
    read.cigar += std::to_string(read.cur_cigar_count);
    read.cigar.push_back(char(read.cur_cigar_op));
    read.cur_cigar_count = 0;
  }
}

void align_reads::output_aligned(read_trace_key key, read_trace read) {
  flush_cigar(read);
  CHECK_EQ(read.read_len_tot, read.seq.size());
  m_on_aligned(std::move(key.read_ids), std::move(read));
}

align_reads::align_reads(const on_aligned_func_t& on_aligned, bool refskip_anchor,
                         pipeline_step_t output)
    : apply_edges_step(std::move(output)),
      m_on_aligned(on_aligned),
      m_refskip_anchor(refskip_anchor) {}

std::ostream& operator<<(std::ostream& os, const aligned_read& read) {
  return os << "AlignedRead([" << read.left_offset << "," << read.right_offset
            << "): " << read.cigar << " seq=" << read.seq << ")";
}

}  // namespace variants
