#include "modules/variants/align.h"
#include "modules/bio_base/kmer.h"

#include <boost/range/irange.hpp>

namespace variants {

static constexpr bool k_align_debug = false;
static constexpr bool k_align_split_debug = false;

void aligner::on_assembly(assembly_ptr a) {
  flush_sorted_to(a->left_offset);

  if (a->matches_reference) {
    // No aligning to do.
    sort_and_output(std::move(a));
    return;
  }

  CHECK_GE(a->left_anchor_len, 0);
  CHECK_GE(a->right_anchor_len, 0);

  CHECK(m_scaffold);
  scaffold s =
      m_scaffold->subscaffold(a->left_offset, a->right_offset - a->left_offset);

  aoffset_t max_len =
      std::min<aoffset_t>(a->seq.size(), a->right_offset - a->left_offset);

  int shared_left = s.shared_prefix_length(a->seq);
  CHECK_GE(shared_left, 0);
  CHECK_GE(shared_left, a->left_anchor_len);

  int shared_right =
      s.rev_comp().shared_prefix_length(dna_slice(a->seq).rev_comp());
  CHECK_GE(shared_right, 0);
  CHECK_GE(shared_right, a->right_anchor_len);

  if (shared_left + shared_right > max_len) {
    shared_left = max_len - shared_right;
  }

  a->left_anchor_len = shared_left;
  a->right_anchor_len = shared_right;

  work_item w;
  w.v.left_offset = a->left_offset + a->left_anchor_len;
  w.v.right_offset = a->right_offset - a->right_anchor_len;
  w.v.seq = a->seq.subseq(
      a->left_anchor_len,
      aoffset_t(a->seq.size()) - a->right_anchor_len - a->left_anchor_len);
  w.s = s.subscaffold(a->left_anchor_len,
                      a->right_offset - a->right_anchor_len -
                          a->left_anchor_len - a->left_offset);

  if (k_align_debug) {
    std::cout << "Starting alignment of " << *a << " against " << s << "\n";
  }

  m_work.emplace_back(std::move(w));
  while (!m_work.empty()) {
    w = std::move(m_work.back());
    m_work.pop_back();
    process(a.get(), std::move(w));
  }

  std::sort(a->aligned_variants.begin(), a->aligned_variants.end());
  if (k_align_debug) {
    std::cout << "Done aligning, produced " << *a << " with "
              << a->aligned_variants.size() << " vars:\n";
    for (const auto& v : a->aligned_variants) {
      std::cout << "  " << v << "\n";
    }
  }

  if (a->aligned_variants.empty()) {
    CHECK_EQ(a->seq.size(), a->right_offset - a->left_offset);

    if (!m_options.trace_reference_assemblies) {
      // Discard this reference-only assembly.
      return;
    }
    
    a->matches_reference = true;
    a->left_anchor_len = 0;
    a->right_anchor_len = 0;
  }

  if (m_options.report_aligned_assemblies_func) {
    m_options.report_aligned_assemblies_func(m_options, *a);
  }

  sort_and_output(std::move(a));
}

void aligner::process(assembly* a, work_item w) {
  if (k_align_debug) {
    std::cout.setf(std::ios::unitbuf);
    std::cout << "Attempting to align " << w.v << " against scaffold " << w.s
              << "\n";
  }

  aoffset_t match_len, seq_match_start, scaffold_match_start;
  if (find_biggest_match_with_ends(m_options, w.v.seq, w.s, &match_len,
                                   &seq_match_start, &scaffold_match_start)) {
    CHECK_GE(match_len, 0);
    work_item left;
    left.v.left_offset = w.v.left_offset;
    left.v.right_offset = w.v.left_offset + scaffold_match_start;
    left.v.seq = w.v.seq.subseq(0, seq_match_start);
    CHECK_GE(scaffold_match_start, 0);
    left.s = w.s.subscaffold(0, scaffold_match_start);

    work_item right;
    right.v.left_offset = w.v.left_offset + scaffold_match_start + match_len;
    right.v.right_offset = w.v.right_offset;
    right.v.seq =
        w.v.seq.subseq(seq_match_start + match_len,
                       w.v.seq.size() - (seq_match_start + match_len));
    CHECK_LE(scaffold_match_start + match_len, w.s.end_pos())
        << "start: " << scaffold_match_start << "len: " << match_len;
    right.s =
        w.s.subscaffold(scaffold_match_start + match_len,
                        w.s.end_pos() - (scaffold_match_start + match_len));

    if (k_align_debug) {
      aligned_var ref_v;
      ref_v.left_offset = left.v.right_offset;
      ref_v.right_offset = right.v.left_offset;
      ref_v.seq = w.v.seq.subseq(seq_match_start, match_len);
      std::cout << "Aligned to " << left.v << " and " << right.v << "\n";
      std::cout << "Middle is: " << ref_v << "\n";
      CHECK_EQ(m_scaffold->subscaffold_str(
                   ref_v.left_offset, ref_v.right_offset - ref_v.left_offset),
               ref_v.seq.as_string());
    }

    if (!left.v.empty()) {
      m_work.emplace_back(std::move(left));
    }

    if (!right.v.empty()) {
      m_work.emplace_back(std::move(right));
    }
  } else {
    if (k_align_debug) {
      std::cout << "Nothing in common found; outputting var: " << w.v << " (empty=" << w.v.empty()
                << ")\n";
    }
    if (!w.v.empty()) {
      a->aligned_variants.emplace_back(std::move(w.v));
    }
  }
}

bool aligner::find_biggest_match_with_ends(const assemble_options& options,
                                           dna_slice seq, const scaffold& s,
                                           aoffset_t* match_len,
                                           aoffset_t* seq_match_start,
                                           aoffset_t* scaffold_match_start) {
  if (k_align_debug) {
    std::cout << "Searching for biggest match with ends\n";
  }
  aoffset_t min_match_size;
  if (find_biggest_match(options, seq, s, match_len, seq_match_start,
                         scaffold_match_start, &min_match_size,
                         anchor_type_t::ANCHORED_TO_BOTH)) {
    if (k_align_debug) {
      std::cout << "Found middle match of size " << *match_len << "\n";
    }
    return true;
  }

  return find_end_matches(options, seq, s, match_len, seq_match_start,
                          scaffold_match_start, min_match_size - 1);
}

bool aligner::find_end_matches(const assemble_options& options, dna_slice seq,
                               const scaffold& s, aoffset_t* match_len,
                               aoffset_t* seq_match_start,
                               aoffset_t* scaffold_match_start,
                               aoffset_t max_match_size) {
  aoffset_t shortest_len = std::min<aoffset_t>(seq.size(), s.end_pos());
  max_match_size = std::min<aoffset_t>(shortest_len, max_match_size);
  for (*match_len = max_match_size; *match_len >= 1; --*match_len) {
    aoffset_t search_len;
    switch (*match_len) {
      case 1:
        search_len = 2;
        break;
      case 2:
        search_len = 3;
        break;
      case 3:
        search_len = 5;
        break;
      case 4:
        search_len = 7;
        break;
      default:
        search_len = *match_len * options.ref_align_factor;
    }
    search_len = std::min<aoffset_t>(shortest_len, search_len);
    if (k_align_debug) {
      std::cout << "Searching for end matches of length " << *match_len
                << " up to " << search_len << " distance from ends\n";
    }
    if (find_match(seq.subseq(0, search_len), s.subscaffold(0, search_len),
                   *match_len, seq_match_start, scaffold_match_start,
                   anchor_type_t::ANCHORED_TO_LEFT)) {
      if (k_align_debug) {
        std::cout << "Found left end match\n";
      }
      return true;
    }
    if (find_match(seq.subseq(seq.size() - search_len, *match_len),
                   s.subscaffold(s.end_pos() - search_len, *match_len),
                   *match_len, seq_match_start, scaffold_match_start,
                   anchor_type_t::ANCHORED_TO_RIGHT)) {
      *seq_match_start += seq.size() - search_len;
      *scaffold_match_start += s.end_pos() - search_len;
      if (k_align_debug) {
        std::cout << "Found right end match\n";
      }
      return true;
    }
  }

  return false;
}

bool aligner::find_biggest_match(const assemble_options& options, dna_slice seq,
                                 const scaffold& s, aoffset_t* match_len,
                                 aoffset_t* seq_match_start,
                                 aoffset_t* scaffold_match_start,
                                 aoffset_t* min_match_size,
                                 anchor_type_t anchor) {
  *min_match_size =
      std::max<int>(seq.size(), s.end_pos()) / options.ref_align_factor;

  if (*min_match_size < 1) {
    *min_match_size = 1;
  }
  if (*min_match_size > options.max_ref_align_bases) {
    *min_match_size = options.max_ref_align_bases;
  }
  if (anchor != anchor_type_t::ANCHORED_TO_BOTH) {
    if (*min_match_size < int(options.min_anchor_drop_overlap)) {
      *min_match_size = options.min_anchor_drop_overlap;
    }
  }

  int max_match_size = std::min<int>(seq.size(), s.end_pos());
  if (max_match_size < *min_match_size) {
    return false;
  }

  if (k_align_debug) {
    std::cout << "Max match size: " << max_match_size
              << " min: " << *min_match_size << "\n";
  }

  *match_len = 0;
  auto match_sizes = boost::irange(*min_match_size, max_match_size + 1);
  auto match_sizes_it = std::upper_bound(
      match_sizes.begin(), match_sizes.end(), false, [&](bool found, int size) {
        aoffset_t seq_start, scaffold_start;

        if (find_match(seq, s, size, &seq_start, &scaffold_start, anchor)) {
          if (size > *match_len) {
            *match_len = size;
            *seq_match_start = seq_start;
            *scaffold_match_start = scaffold_start;
          }
          return found;
        }
        return !found;
      });

  if (!*match_len) {
    CHECK(match_sizes_it == match_sizes.begin());
    return false;
  }

  int matched_ub = match_sizes_it - match_sizes.begin() + *min_match_size;
  CHECK_EQ(matched_ub, *match_len + 1);
  return true;
}

bool aligner::find_match(dna_slice seq, const scaffold& s, int match_len,
                         aoffset_t* seq_match_start,
                         aoffset_t* scaffold_match_start,
                         anchor_type_t anchor) {
  if (match_len > s.end_pos() || match_len > int(seq.size())) {
    return false;
  }
  aoffset_t seq_anchor = 0;
  aoffset_t scaffold_anchor = 0;
  switch (anchor) {
    case anchor_type_t::ANCHORED_TO_BOTH:
      seq_anchor = (seq.size() - match_len) / 2;
      scaffold_anchor = (s.end_pos() - match_len) / 2;
      break;
    case anchor_type_t::ANCHORED_TO_LEFT:
      seq_anchor = 0;
      scaffold_anchor = 0;
      break;
    case anchor_type_t::ANCHORED_TO_RIGHT:
      seq_anchor = seq.size();
      scaffold_anchor = s.end_pos();
      break;
  }

  int kmer_size = match_len;
  if (kmer_size > 30) {
    kmer_size = 30;
  }

  std::unordered_multimap<kmer_t, aoffset_t, unsalted_hash> seq_kmers;
  aoffset_t seq_offset = 0;
  for (const auto& k : kmer_view(seq, kmer_size)) {
    if (seq_offset + match_len > int(seq.size())) {
      break;
    }
    seq_kmers.insert(std::make_pair(k, seq_offset));
    ++seq_offset;
  }
  if (int(seq.size()) >= match_len) {
    CHECK_EQ(seq_offset + match_len, seq.size() + 1);
  } else {
    CHECK_EQ(seq_offset, 0);
  }

  int64_t best_distance = std::numeric_limits<int64_t>::max();
  for (const auto& e : s.extents()) {
    aoffset_t scaffold_offset = e.offset;
    aoffset_t extent_end = e.offset + e.sequence.size();
    for (const auto& k : kmer_view(e.sequence, kmer_size)) {
      if (scaffold_offset + match_len > extent_end) {
        break;
      }

      auto eq_range = seq_kmers.equal_range(k);
      if (eq_range.first == eq_range.second) {
        ++scaffold_offset;
        continue;
      }

      int64_t scaffold_distance = scaffold_offset - scaffold_anchor;

      for (auto it = eq_range.first; it != eq_range.second; ++it) {
        aoffset_t seq_offset = it->second;
        int64_t seq_distance = seq_offset - seq_anchor;

        int64_t distance = llabs(seq_distance - scaffold_distance);

        if (distance > best_distance) {
          continue;
        }

        if (match_len > kmer_size) {
          scaffold sub = s.subscaffold(scaffold_offset, match_len);
          if (!sub.is_simple() ||
              sub.get_simple() != seq.subseq(seq_offset, match_len)) {
            continue;
          }
        }

        best_distance = distance;
        *seq_match_start = seq_offset;
        *scaffold_match_start = scaffold_offset;
      }
      ++scaffold_offset;
    }
    if (aoffset_t(e.sequence.size()) >= match_len) {
      CHECK_EQ(scaffold_offset + match_len, extent_end + 1);
    } else {
      CHECK_EQ(scaffold_offset, e.offset);
    }
  }

  return best_distance != std::numeric_limits<int64_t>::max();
}

void align_splitter::set_matches_reference(assembly& a) {
  CHECK_EQ(a.seq.size(), a.right_offset - a.left_offset);
  a.matches_reference = true;
  a.left_anchor_len = 0;
  a.right_anchor_len = 0;
}

void align_splitter::on_assembly(assembly_ptr a) {
  flush_sorted_to(a->left_offset);

  if (k_align_split_debug) {
    std::cout.setf(std::ios::unitbuf);
    std::cout << "Splitting alignment of " << *a << " with "
              << a->aligned_variants.size() << " vars:\n";
    for (const auto& v : a->aligned_variants) {
      std::cout << "  " << v << "\n";
    }
  }

  aoffset_t ref_offset = a->left_offset;

  aoffset_t orig_left = a->left_offset;
  aoffset_t orig_right = a->right_offset;

  std::vector<aligned_var> vars = a->aligned_variants;

  for (const auto& v : vars) {
    if (k_align_split_debug) {
      std::cout << "Splitting var '" << v
                << "', starting at ref offset=" << ref_offset << " and seq "
                << a->seq << "\n";
    }
    CHECK_GE(v.left_offset, orig_left);
    CHECK_LE(v.right_offset, orig_right);
    CHECK_GE(v.left_offset, ref_offset);

    assembly_ptr left, right;
    if (v.left_offset != ref_offset) {
      CHECK_GT(v.left_offset, ref_offset);
      aoffset_t ref_chunk_size = v.left_offset - ref_offset;
      CHECK_GE(a->seq.size(), ref_chunk_size);

      std::tie(left, right) =
          split_assembly(std::move(a), ref_chunk_size, ref_chunk_size);
      set_matches_reference(*left);
      if (k_align_split_debug) {
        std::cout << "Emitting ref section: " << *left << "\n";
      }
      sort_and_output(std::move(left));

      a = std::move(right);

      CHECK_EQ(ref_offset + ref_chunk_size, v.left_offset);
      ref_offset = v.left_offset;

      if (k_align_split_debug) {
        std::cout << "Emitted ref section, now starting at ref offset="
                  << ref_offset << " and seq " << a->seq << "\n";
        std::cout << "Assembly = " << *a << "\n";
      }
    }
    std::tie(left, right) = split_assembly(std::move(a), v.seq.size(),
                                           v.right_offset - v.left_offset);
    CHECK_EQ(left->seq, v.seq);
    if (k_align_split_debug) {
      std::cout << "Emitting var section: " << *left << "\n";
    }
    sort_and_output(std::move(left));
    a = std::move(right);
    ref_offset = a->left_offset;
  }

  CHECK_EQ(ref_offset, a->left_offset);
  if (a->left_offset == a->right_offset && a->seq.size() == 0) {
    return;
  }
  set_matches_reference(*a);
  if (k_align_split_debug) {
    std::cout << "Emitting final ref section: " << *a << "\n";
  }
  sort_and_output(std::move(a));
}

std::ostream& operator<<(std::ostream& os,
                         const aligner::anchor_type_t& anchor) {
  switch (anchor) {
    case aligner::anchor_type_t::ANCHORED_TO_LEFT:
      return os << "ANCHORED_TO_LEFT";
    case aligner::anchor_type_t::ANCHORED_TO_RIGHT:
      return os << "ANCHORED_TO_RIGHT";
    case aligner::anchor_type_t::ANCHORED_TO_BOTH:
      return os << "ANCHORED_TO_BOTH";
  }
  return os << "BAD ANCHOR VALUE " << int(anchor);
}

}  // namespace variants
