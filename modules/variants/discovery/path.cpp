#include "modules/variants/discovery/path.h"

namespace variants {
namespace discovery {

path::path(const readmap* rm, dna_slice seq, seqset_range r, unsigned path_overlap,
           unsigned bases_since_read, int anchor_len)
    : m_rm(rm),
      m_rc_seq(seq.rev_comp()),
      m_r(r),
      m_path_overlap(path_overlap),
      m_last_overlap(path_overlap),
      m_bases_since_read(bases_since_read),
      m_anchor_len(anchor_len) {
  CHECK_GE(m_bases_since_read, 0) << *this;
  CHECK_GE(m_anchor_len, 0) << *this;
  populate_longest_read();
  CHECK_GE(seq.size(), r.size()) << *this;
  CHECK_GE(seq.size(), m_anchor_len) << *this;
}

void path::push_front_drop(dna_base b) {
  CHECK(m_r.valid());
  return push_front_drop(b, m_r.push_front_drop(b));
}

void path::push_front_drop(dna_slice seq) {
  for (dna_base b : seq.rev_comp()) {
    push_front_drop(b.complement());
  }
}

void path::push_front_drop(dna_base b, const seqset_range& pushed) {
  CHECK(pushed.valid());
  CHECK_EQ(pushed.front(), b);

  if (m_bases_until_loop_check == 0) {
    if (pushed.begin() == m_loop_check_seqset_id) {
      m_loop_detected = true;
    }

    m_loop_check_seqset_id = pushed.begin();
    ++m_loop_check_interval;
    m_bases_until_loop_check = m_loop_check_interval;
  } else {
    --m_bases_until_loop_check;
  }

  ++m_bases_since_read;
  m_r = pushed;
  m_rc_seq.push_back(b.complement());
  populate_longest_read();
  int read_len = m_longest_read_id ? m_rm->get_readlength(*m_longest_read_id) : m_r.size();

  int overlap = read_len - m_bases_since_read;
  if (overlap < 0) {
    overlap = 0;
  }

  if (overlap < m_path_overlap) {
    m_path_overlap = overlap;
  }

  if (m_longest_read_id) {
    m_bases_since_read = 0;
    m_last_overlap = overlap;
  }

  m_tot_overlap += overlap;
  m_tot_overlap_bases++;
}

void path::populate_longest_read() { m_longest_read_id = m_rm->get_longest_prefix_read_id(m_r); }

void path::display_fwd(std::ostream& os) const {
  os << "path(r=" << m_r.sequence() << ", seq=" << m_rc_seq.rev_comp() << ", pol=" << m_path_overlap;
  if (m_tot_overlap_bases) {
    os << ", aol=" << (m_tot_overlap * 1.0 / m_tot_overlap_bases);
  } else {
    os << ", aol=(none)";
  }
  os << ", bases since read=" << m_bases_since_read;
  if (m_longest_read_id) {
    os << ", longest read="
       << m_rm->get_read_by_id(*m_longest_read_id).get_seqset_entry().sequence();
  }
  os << ")";
}

void path::display_rev(std::ostream& os) const {
  os << "path(rc_seq=" << m_rc_seq << ",rc_r=" << m_r.sequence().rev_comp()
     << ", ol=" << m_path_overlap << ", bases since read=" << m_bases_since_read;
  if (m_longest_read_id) {
    os << ", longest read rc="
       << m_rm->get_read_by_id(*m_longest_read_id).get_seqset_entry().sequence().rev_comp();
  }
  os << ")";
}

void path::check_invariants() const {
  CHECK_GE(seq().size(), anchor_len());
  dna_sequence r_seq = range().sequence();
  CHECK_EQ(seq().subseq(0, r_seq.size()), r_seq) << *this;
}

}  // namespace discovery
}  // namespace variants
