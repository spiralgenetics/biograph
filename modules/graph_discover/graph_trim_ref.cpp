#include "modules/graph_discover/graph_trim_ref.h"
#include "modules/bio_base/kmer.h"

namespace variants {

static constexpr bool k_dbg = false;

graph_trim_ref::graph_trim_ref(const assemble_options& options, pipeline_step_t output)
    : sorted_output_pipeline_step(std::move(output)),
      m_options(options),
      m_scaffold(options.scaffold) {
  CHECK(m_scaffold);
}

graph_trim_ref::~graph_trim_ref() { flush(); }

void graph_trim_ref::flush() {
  advance_to(std::numeric_limits<aoffset_t>::max());
  flush_sorted();
  CHECK(m_ref_stops.empty());
  CHECK(m_ref_asms.empty());
}

void graph_trim_ref::advance_to(aoffset_t offset) {
  while (!m_ref_asms.empty()) {
    auto& a = m_ref_asms.front();
    if (a->right_offset > offset) {
      return;
    }

    split_and_output_ref(std::move(a));
    m_ref_asms.pop_front();
  }
  flush_sorted_to(offset);
  while (!m_ref_stops.empty() && *m_ref_stops.begin() < sort_flush_point()) {
    m_ref_stops.erase(m_ref_stops.begin());
  }
}

void graph_trim_ref::split_and_output_ref(assembly_ptr a) {
  CHECK(a->matches_reference);
  untrack_left_offset(a->left_offset);

  std::vector<aoffset_t> ref_stops(m_ref_stops.begin(), m_ref_stops.end());

  for (aoffset_t ref_stop : ref_stops) {
    if (ref_stop <= a->left_offset) {
      if (k_dbg) {
        std::cerr << "Discarding unused refstop: " << ref_stop << "\n";
      }
      continue;
    }

    if (ref_stop >= a->right_offset) {
      break;
    }
    if (k_dbg) {
      std::cerr << "Splitting ref at " << ref_stop << ": " << *a << "\n";
    }
    assembly_ptr left, right;
    aoffset_t left_offset = a->left_offset;
    std::tie(left, right) =
        split_assembly_absoffset(std::move(a), ref_stop - left_offset, ref_stop);

    sort_and_output(std::move(left));
    a = std::move(right);
  }
  sort_and_output(std::move(a));
}

void graph_trim_ref::on_assembly(assembly_ptr a) {
  static std::mutex g_mu;
  std::lock_guard<std::mutex> l(g_mu);
  if (k_dbg) {
    std::cerr << "graph_trim_ref processing assembly " << *a << "\n";
  }
  advance_to(min(a->left_offset, a->right_offset) - m_max_backtrack);

  if (a->matches_reference) {
    track_left_offset(a->left_offset);
    m_ref_asms.push_back(std::move(a));
    return;
  }

  int shared_left = 0;
  if (a->left_offset) {
    dna_slice left_seq = m_scaffold->split_extent_at(a->left_offset).second;
    shared_left = left_seq.shared_prefix_length(dna_slice(a->seq));
  }
  int shared_right = 0;
  if (a->right_offset) {
    dna_slice right_seq = m_scaffold->split_extent_at(a->right_offset).first;
    shared_right = right_seq.rev_comp().shared_prefix_length(dna_slice(a->seq).rev_comp());
  }

  aoffset_t max_anchor_len = a->seq.size();
  if (k_dbg) {
    std::cerr << "graph_trim shared left=" << shared_left << " right=" << shared_right << "\n";
  }
  if (a->left_offset && a->right_offset) {
    max_anchor_len = std::min<aoffset_t>(max_anchor_len, a->right_offset - a->left_offset);
  }
  if (!a->left_offset) {
    max_anchor_len = std::min<aoffset_t>(max_anchor_len, m_max_backtrack);
    if (shared_right > m_max_backtrack) {
      std::cerr << "Max backtrack exceeded on assembly: " << *a << "\n";
    }
  }

  if (shared_right > max_anchor_len) {
    shared_right = max_anchor_len;
  }

  if (shared_right + shared_left > max_anchor_len) {
    shared_left = max_anchor_len - shared_right;
  }

  if (k_dbg) {
    std::cerr << "Adjusted Shared left: " << shared_left << " right: " << shared_right
              << " max anchor len: " << max_anchor_len << "\n";
  }

  assembly_ptr discard;
  optional_aoffset left_offset = a->left_offset;
  std::tie(discard, a) = split_assembly_absoffset(
      std::move(a), shared_left,
      left_offset ? optional_aoffset(left_offset + shared_left) : optional_aoffset::none);

  if (k_dbg) {
    std::cerr << "After left discard: " << *a << "\n";
  }
  optional_aoffset right_offset = a->right_offset;
  aoffset_t seq_size = a->seq.size();
  std::tie(a, discard) = split_assembly_absoffset(
      std::move(a), seq_size - shared_right,
      right_offset ? optional_aoffset(right_offset - shared_right) : optional_aoffset::none);

  if (k_dbg) {
    std::cerr << "After right discard: " << *a << "\n";
  }

  if (a->seq.size() == 0 &&
      (a->left_offset == a->right_offset || !a->left_offset || !a->right_offset)) {
    if (k_dbg) {
      std::cout << "graph_trim_ref dropping variant that entirely matches reference: " << *a
                << "\n";
    }
    return;
  }
  if (a->left_offset) {
    m_ref_stops.insert(a->left_offset);
  }
  if (a->right_offset) {
    m_ref_stops.insert(a->right_offset);
  }

  sort_and_output(std::move(a));
}

}  // namespace variants
