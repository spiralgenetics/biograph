#include "modules/variants/assembly_dot.h"

namespace variants {

namespace {

constexpr char k_dot_header[] = R"HEADER(
digraph G {
  mode="hier";  ranksep=1;  newrank="true";
  node [shape=record];
)HEADER";

constexpr char k_dot_footer[] = R"FOOTER(
}
)FOOTER";

constexpr size_t k_max_seq_size = 20;

std::string trimseq(std::string seq) {
  if (seq.size() <= k_max_seq_size) {
    return seq;
  }
  return seq.substr(0, k_max_seq_size / 2) + "..." +
         seq.substr(seq.size() - k_max_seq_size / 2);
}

}  // namespace

assembly_dot::assembly_dot(const scaffold& s) : m_scaffold(s) {
  m_result << k_dot_header;
}

void assembly_dot::add_assembly(const assembly& a) {
  CHECK(!m_finalized);

  if (a.matches_reference) {
    return;
  }
  ++m_id;

  std::string id = "A" + std::to_string(m_id);

  dna_slice left_anchor = dna_slice(a.seq).subseq(0, a.left_anchor_len);
  dna_slice variant = dna_slice(a.seq).subseq(
      a.left_anchor_len,
      dna_slice(a.seq).size() - a.right_anchor_len - a.left_anchor_len);
  dna_slice right_anchor = dna_slice(a.seq).subseq(
      a.seq.size() - a.right_anchor_len, a.right_anchor_len);

  m_result << "A" << m_id << " [";
  if (a.matches_reference) {
    m_result << "color=\"blue\", ";
  } else {
    m_result << "color=\"red\", ";
  }
  m_result << "vsort=10000, ";
  m_result << "label=\"{";
  m_result << "left=" << a.left_offset;

  if (!a.matches_reference) {
    m_result << "+" << a.left_anchor_len << "="
           << a.left_offset + a.left_anchor_len;
    m_result << " |" << trimseq(left_anchor.as_string());
  }

  if (a.matches_reference) {
    m_ref_right_edges.insert(std::make_pair(a.left_offset, id + ":ref"));
    m_ref_left_edges.insert(std::make_pair(a.right_offset, id + ":ref"));
    m_result << " | <ref> " << trimseq(a.seq.as_string());
  } else if (a.aligned_variants.empty()) {
    m_ref_right_edges.insert(
        std::make_pair(a.left_offset + a.left_anchor_len, id + ":variant"));
    m_ref_left_edges.insert(
        std::make_pair(a.right_offset - a.right_anchor_len, id + ":variant"));
    m_result << " | <variant> " << trimseq(variant.as_string());
  } else {
    aoffset_t var_id = 0;
    aoffset_t ref_offset = a.left_offset + a.left_anchor_len;
    aoffset_t seq_offset = a.left_anchor_len;
    for (const auto& var : a.aligned_variants) {
      CHECK_GE(var.left_offset, ref_offset);
      if (var.left_offset > ref_offset) {
        aoffset_t advance_dist = var.left_offset - ref_offset;
        m_result << " | "
                 << trimseq(a.seq.subseq(seq_offset, advance_dist).as_string());
        ref_offset += advance_dist;
        seq_offset += advance_dist;
      }
      std::string var_id_str = "v" + std::to_string(var_id);
      m_result << " | <" << var_id_str << "> " << trimseq(var.seq.as_string());
      CHECK_EQ(var.seq, a.seq.subseq(seq_offset, var.seq.size()));
      m_ref_right_edges.insert(
          std::make_pair(ref_offset, id + ":" + var_id_str));
      ref_offset = var.right_offset;
      seq_offset += var.seq.size();
      m_ref_left_edges.insert(
          std::make_pair(ref_offset, id + ":" + var_id_str));
      ++var_id;
    }
    if (a.right_offset - a.right_anchor_len > ref_offset) {
      aoffset_t advance_dist = a.right_offset - a.right_anchor_len - ref_offset;
      m_result << " | "
               << trimseq(a.seq.subseq(seq_offset, advance_dist).as_string());
      ref_offset += advance_dist;
      seq_offset += advance_dist;
    }
  }

  m_result << " | id=" << a.assembly_id << " ol= " << a.min_overlap
           << " score=" << a.score;
  if (a.matches_reference) {
    m_result << " REF";
  }
  m_result << " | " << trimseq(right_anchor.as_string());
  m_result << " | right=" << a.right_offset - a.right_anchor_len;
  if (!a.matches_reference) {
    m_result << "+"
             << a.right_anchor_len << "=" << a.right_offset;
  }
  m_result << "}\"]\n";
}

void assembly_dot::finalize() {
  m_finalized = true;

  if (m_ref_right_edges.empty() || m_ref_left_edges.empty()) {
    m_result << k_dot_footer;
    return;
  }

  std::set<aoffset_t> breakpoints;
  for (const auto& elist : {m_ref_right_edges, m_ref_left_edges}) {
    for (const auto& e : elist) {
      breakpoints.insert(e.first);
    }
  }

  if (m_ref_right_edges.begin()->first <= m_ref_left_edges.begin()->first) {
    breakpoints.insert(m_ref_right_edges.begin()->first - 10);
  }
  if (m_ref_left_edges.rbegin()->first >= m_ref_right_edges.rbegin()->first) {
    breakpoints.insert(m_ref_left_edges.rbegin()->first + 10);
  }

  auto it = breakpoints.begin();
  auto next = it;
  ++next;

  CHECK(next != breakpoints.end());
  std::string last_id;

  int pos = 0;

  while (next != breakpoints.end()) {
    aoffset_t start = *it;
    aoffset_t limit = *next;

    std::string id = "R" + std::to_string(start) + "to" + std::to_string(limit);

    pos++;
    m_result << id << " [color=\"green\", vsort=" << pos
             << ", weight=\"100\", pos=\"0," << (pos * 20) << "!\", label=\"{";

    m_result << trimseq(m_scaffold.subscaffold_str(start, limit - start));
    m_result << "|" << start << " - " << limit << "}\"]\n";

    for (std::string dir : {"n", "s"}) {
      const auto& edges = (dir == "s") ? m_ref_right_edges : m_ref_left_edges;
      aoffset_t pos = (dir == "s") ? limit : start;
      auto r = edges.equal_range(pos);
      for (auto it = r.first; it != r.second; ++it) {
        std::string dest = it->second;
        std::string src = id + ":" + dir;

        if (dir == "n") {
          std::swap(src, dest);
        }
        m_result << src << " -> " << dest;
        if (dir == "s") {
          m_result << " [constraint=false]";
        }
        m_result << "\n";
      }
    }

    if (last_id != "") {
      m_result << last_id << ":s -> " << id
               << ":n [color=\"green\", weight=500]\n";
    }
    last_id = id;

    ++it;
    ++next;
  }
  m_result << k_dot_footer;
}

}  // namespace variants
