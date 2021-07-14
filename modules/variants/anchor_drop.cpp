#include "modules/variants/anchor_drop.h"
#include "modules/bio_base/kmer.h"
#include "modules/variants/align.h"

namespace variants {

constexpr bool k_anchor_drop_debug = false;

unsigned anchor_dropper::g_max_kmer_size = 31;

anchor_dropper::anchor_dropper(const assemble_options& options, pipeline_step_t output)
    : sorted_output_pipeline_step(std::move(output)),
      m_options(options),
      m_scaffold(options.scaffold) {
  set_expected_order(assembly::left_offset_less_than);

  m_read_ahead.scaffold_it = m_scaffold->begin();
  m_trail_behind.scaffold_it = m_scaffold->begin();

  m_kmer_size = std::min<unsigned>(m_options.min_overlap, g_max_kmer_size);
  m_kmer_mask = ~((~(kmer_t)0) << (m_kmer_size * 2));

  m_kmer_skip_len = m_options.min_overlap - m_kmer_size;

  if (k_anchor_drop_debug) {
    std::cout << "Kmer size " << m_kmer_size << " mask: " << kmer_str(m_kmer_mask, 32) << "\n";
  }
}

void anchor_dropper::advance_kmer(
    kmer_it& it, aoffset_t offset,
    const std::function<void(kmer_t, aoffset_t)>& process_kmer_f) const {
  while (it.scaffold_it.offset() < offset) {
    if (it.scaffold_it == m_scaffold->end()) {
      return;
    }

    if (process_kmer_f) {
      if (it.scaffold_it.first_in_extent()) {
        it.kmer_needed = m_kmer_size - 1;
      }

      it.kmer <<= 2;
      it.kmer |= kmer_t(int(*it.scaffold_it));

      if (it.kmer_needed) {
        --it.kmer_needed;
      } else {
        process_kmer_f(it.kmer & m_kmer_mask, it.scaffold_it.offset() + 1 - m_kmer_size);
        it.kmer_needed = m_kmer_skip_len;
      }
    }
    ++it.scaffold_it;
  }
}

void anchor_dropper::skip_to(kmer_it& it, aoffset_t offset) const {
  if (it.scaffold_it.offset() < offset && it.scaffold_it != m_scaffold->end()) {
    it.scaffold_it.skip_to(offset, "anchor_dropper");
    it.kmer_needed = m_kmer_size - 1;
  }
}

void anchor_dropper::advance_read_ahead_to(aoffset_t target_offset) {
  advance_kmer(m_read_ahead, target_offset, [&](kmer_t kmer, aoffset_t offset) {
    if (k_anchor_drop_debug) {
      std::cout << "Adding kmer ahead at " << offset << ": " << kmer_str(kmer, m_kmer_size) << "\n";
      std::cout.flush();
      CHECK_EQ(kmer_str(kmer, m_kmer_size), m_scaffold->subscaffold_str(offset, m_kmer_size));
    }
    m_kmers.insert(std::make_pair(kmer, offset));
  });
}

void anchor_dropper::advance_range_to(aoffset_t trail_behind, aoffset_t read_ahead) {
  // Advance over unneeded regions without kmerizing.  TODO(nils): We
  // could be even more efficient by clearning the trail behind
  // completely instead of advancing through it.
  skip_to(m_read_ahead, trail_behind - m_kmer_size);

  advance_read_ahead_to(read_ahead);
  advance_trail_behind_to(trail_behind);
}

void anchor_dropper::advance_trail_behind_to(aoffset_t target_offset) {
  advance_kmer(m_read_ahead, target_offset, [&](kmer_t kmer, aoffset_t offset) {
    if (k_anchor_drop_debug) {
      std::cout << "Removing kmer behind at " << offset << ": " << kmer_str(kmer, m_kmer_size)
                << "\n";
    }
    auto r = m_kmers.equal_range(kmer);
    bool found = false;
    for (auto it = r.first; it != r.second; ++it) {
      if (it->second == offset) {
        m_kmers.erase(it);
        found = true;
        break;
      }
    }
    CHECK(found);
  });
}

bool anchor_dropper::try_long_rejoin(assembly_ptr& a) {
  if (a->seq.size() < m_options.min_overlap + a->left_anchor_len) {
    return false;
  }

  aoffset_t seq_offset = 0;
  aoffset_t best_match_len = 0;
  aoffset_t best_match_seq_offset = 0;
  aoffset_t best_match_ref_offset = 0;

  dna_slice seq_slice(a->seq.begin(), a->seq.end());

  for (kmer_t kmer : kmer_view(seq_slice, m_kmer_size)) {
    if (seq_offset < a->left_anchor_len) {
      ++seq_offset;
      continue;
    }

    if (k_anchor_drop_debug || assembly_needs_trace(*a)) {
      std::cout << "Long drop looking for kmer " << kmer_str(kmer, m_kmer_size) << "\n";
    }

    auto r = m_kmers.equal_range(kmer);
    for (auto it = r.first; it != r.second; ++it) {
      if (it->second < a->left_offset + a->left_anchor_len) {
        continue;
      }
      aoffset_t max_before_len = std::min<aoffset_t>(
          seq_offset - a->left_anchor_len, it->second - a->left_offset - a->left_anchor_len);

      dna_slice ref_subseq_before, ref_subseq_after;
      std::tie(ref_subseq_before, ref_subseq_after) = m_scaffold->split_extent_at(it->second);

      dna_slice seq_subseq_before =
          seq_slice.subseq(a->left_anchor_len, seq_offset - a->left_anchor_len);
      dna_slice seq_subseq_after = seq_slice.subseq(seq_offset, seq_slice.size() - seq_offset);
      aoffset_t before_shared =
          ref_subseq_before.rev_comp().shared_prefix_length(seq_subseq_before.rev_comp());
      if (before_shared > max_before_len) {
        before_shared = max_before_len;
      }

      aoffset_t after_shared = ref_subseq_after.shared_prefix_length(seq_subseq_after);
      aoffset_t shared = before_shared + after_shared;
      if (k_anchor_drop_debug || assembly_needs_trace(*a)) {
        std::cout << "Comparing ref before seq at '" << ref_subseq_before << "' to '"
                  << seq_subseq_before << "'\n";
        std::cout << "Comparing ref after seq '" << ref_subseq_after << "' to '" << seq_subseq_after
                  << "'\n";
        std::cout << "Shared: " << before_shared << " + " << after_shared << " = " << shared
                  << "\n";
      }

      if (shared <= best_match_len || shared < aoffset_t(m_options.min_overlap)) {
        continue;
      }
      if (shared == best_match_len && (it->second - before_shared) > best_match_ref_offset) {
        // Closer drops are better.
        continue;
      }

      best_match_len = shared;
      best_match_seq_offset = seq_offset - before_shared;
      best_match_ref_offset = it->second - before_shared;
    }

    ++seq_offset;
  }

  if (best_match_len == 0) {
    return false;
  }
  CHECK_GE(best_match_len, m_options.min_overlap);
  auto new_seq_len = best_match_seq_offset + best_match_len;
  CHECK_GT(best_match_ref_offset, a->left_offset);
  a->right_offset = best_match_ref_offset + best_match_len;
  aoffset_t ref_len = a->right_offset - a->left_offset;
  assembly_ptr anchored =
      split_assembly(std::move(a), new_seq_len, ref_len).first;
  anchored->right_anchor_len = best_match_len;
  if (anchored->right_anchor_len + anchored->left_anchor_len >
      anchored->right_offset - anchored->left_offset) {
    anchored->right_anchor_len =
        anchored->right_offset - anchored->left_offset - anchored->left_anchor_len;
  }
  sort_and_output(std::move(anchored));
  return true;
}

void anchor_dropper::on_assembly(assembly_ptr a) {
  flush_sorted_to(a->left_offset);

  if (a->right_anchor_len || a->matches_reference) {
    if (!a->matches_reference) {
      CHECK(a->left_anchor_len);
    }
    sort_and_output(std::move(a));
    return;
  }

  if (k_anchor_drop_debug || assembly_needs_trace(*a)) {
    std::cout << "Anchor drop processing assembly " << *a << "\n";
  }

  CHECK(a->left_anchor_len);

  advance_range_to(a->left_offset, a->left_offset + m_options.read_ahead_distance);

  if (a->left_anchor_len >= aoffset_t(a->seq.size())) {
    if (k_anchor_drop_debug || assembly_needs_trace(*a)) {
      std::cout << "Left anchor len bigger than sequence; dropping: " << *a << "\n";
    }
    return;
  }
  CHECK_LT(a->left_anchor_len, a->seq.size());

  if (try_long_rejoin(a)) {
    return;
  }

  boost::optional<aoffset_t> new_ref_len;
  boost::optional<aoffset_t> new_seq_len;
  boost::optional<aoffset_t> new_match_len;

  // Ugh, https://gcc.gnu.org/bugzilla/show_bug.cgi?id=47679
  for (boost::optional<aoffset_t>* opt : {&new_ref_len, &new_seq_len, &new_match_len}) {
    opt->emplace(0);
    opt->reset();
  }

  aoffset_t seq_variant_start = 0;
  aoffset_t ref_variant_start = a->left_offset;

  for (;;) {
    // Advance through reference on left.
    dna_slice variant_seq =
        dna_slice(a->seq).subseq(seq_variant_start, a->seq.size() - seq_variant_start);
    scaffold s = m_scaffold->subscaffold(ref_variant_start, variant_seq.size());
    if (k_anchor_drop_debug || assembly_needs_trace(*a)) {
      std::cout << "Right anchor aligning " << variant_seq << " against " << s << "\n";
    }
    unsigned shared_start = s.shared_prefix_length(variant_seq);
    seq_variant_start += shared_start;
    CHECK_GE(seq_variant_start, a->left_anchor_len) << *a;
    ref_variant_start += shared_start;

    variant_seq = dna_slice(a->seq).subseq(seq_variant_start, a->seq.size() - seq_variant_start);
    if (variant_seq.size() == 0) {
      if (k_anchor_drop_debug || assembly_needs_trace(*a)) {
        std::cout << "Variant seq is empty\n";
      }
      break;
    }
    scaffold variant_s = m_scaffold->subscaffold(ref_variant_start, variant_seq.size() * 5 / 4 + 3);

    if (k_anchor_drop_debug || assembly_needs_trace(*a)) {
      std::cout << "After trim of " << shared_start << ":\n"
                << "Aligning:\n" << variant_seq << "\nagainst:\n" << variant_s << "\n";
    }

    int match_len;
    aoffset_t seq_match_start, scaffold_match_start;
    int min_match_size;
    if (!aligner::find_biggest_match(m_options, variant_seq, variant_s, &match_len,
                                     &seq_match_start, &scaffold_match_start, &min_match_size,
                                     aligner::anchor_type_t::ANCHORED_TO_LEFT)) {
      if (k_anchor_drop_debug || assembly_needs_trace(*a)) {
        std::cout << "No match found\n";
      }
      break;
    }

    if (k_anchor_drop_debug || assembly_needs_trace(*a)) {
      std::cout << "Common of length " << match_len << ": "
                << variant_seq.subseq(seq_match_start, match_len) << "\n";
    }

    seq_variant_start += seq_match_start + match_len;
    ref_variant_start += scaffold_match_start + match_len;

    new_seq_len.emplace(seq_variant_start);
    new_ref_len.emplace(ref_variant_start - a->left_offset);
    new_match_len.emplace(match_len);
  }

  if (new_seq_len && new_ref_len && new_match_len) {
    // Discard everything after the match.
    a->right_offset = a->left_offset + *new_ref_len;
    assembly_ptr anchored, discard;
    std::tie(anchored, discard) = split_assembly(std::move(a), *new_seq_len, *new_ref_len);
    anchored->right_anchor_len = *new_match_len;
    if (k_anchor_drop_debug || assembly_needs_trace(*anchored)) {
      std::cout << "After right anchor drop: " << *anchored << "\n";
      std::cout << "Discarding portion: " << *discard << "\n";
    }
    sort_and_output(std::move(anchored));

    if (m_options.report_half_aligned_func && discard->seq.size()) {
      CHECK_EQ(discard->left_anchor_len, 0);
      half_aligned_assembly ha;
      ha.scaffold_name = m_options.scaffold_name;
      ha.offset = discard->left_offset + discard->left_anchor_len;
      ha.right_anchor = false;  // anchored on the left
      ha.seq = discard->seq;
      ha.assembly_id = discard->assembly_id;
      ha.rc_read_ids = discard->rc_read_ids;
      m_options.report_half_aligned_func(ha);
    }
  } else {
    if (k_anchor_drop_debug || assembly_needs_trace(*a)) {
      std::cout << "Could not align right anchor.\n";
    }
    if (m_options.report_half_aligned_func) {
      dna_slice seq = a->seq;
      CHECK_GT(seq.size(), a->left_anchor_len);
      seq = seq.subseq(a->left_anchor_len, seq.size() - a->left_anchor_len);

      half_aligned_assembly ha;
      ha.scaffold_name = m_options.scaffold_name;
      ha.offset = a->left_offset + a->left_anchor_len;
      ha.right_anchor = false;  // anchored on the left
      ha.seq = dna_sequence(seq);
      ha.assembly_id = a->assembly_id;
      ha.rc_read_ids = a->rc_read_ids;
      m_options.report_half_aligned_func(ha);
    }
  }
}

}  // namespace variants
