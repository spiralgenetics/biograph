#include "modules/variants/split_variants.h"

#include "modules/bio_base/readmap.h"
#include "modules/bio_base/seqset.h"

namespace variants {

constexpr bool k_split_var_debug = false;

split_variants::split_variants(const assemble_options& options, pipeline_step_t output)
    : sorted_output_pipeline_step(std::move(output)), m_options(options) {
  CHECK(m_options.seqset);
  CHECK(m_options.readmap);
  CHECK(m_options.scaffold);

  m_ref.cov.cur = m_options.seqset->ctx_begin();
  m_ref.scaffold_it = m_options.scaffold->begin();
}

split_variants::~split_variants() {
  if (!m_active.empty()) {
    chunk_and_output();
  }
  flush_sorted();

  CHECK(m_active.empty());
}

void split_variants::on_assembly(assembly_ptr a) {
  if (a->matches_reference) {
    return;
  }
  CHECK(!a->aligned_variants.empty()) << *a << " must be aligned.";

  if (a->left_offset >= m_rightmost_offset) {
    if (!m_active.empty()) {
      chunk_and_output();
    }
    flush_sorted_to(a->left_offset);
    m_leftmost_offset = a->left_offset;
  }

  m_rightmost_offset = std::max<aoffset_t>(a->right_offset, m_rightmost_offset);

  CHECK_NE(a->seq.size(), 0) << *a;

  CHECK(!a->aligned_variants.empty());
  for (const auto& v : a->aligned_variants) {
    if (v.left_offset == v.right_offset) {
      auto var_int = boost::icl::interval<aoffset_t>::right_open(v.left_offset - 1, v.right_offset);
      m_variant_regions += var_int;
      if (k_split_var_debug) {
        std::cout << "Adding variant insert at " << v.left_offset << "\n";
      }
    } else {
      auto var_int = boost::icl::interval<aoffset_t>::right_open(v.left_offset, v.right_offset);
      m_variant_regions += var_int;
      if (k_split_var_debug) {
        std::cout << "Adding variant region " << var_int << " for " << v << "\n";
      }
    }
    if (k_split_var_debug) {
      std::cout << "New variant regions: " << m_variant_regions << "\n";
    }
  }

  asm_info i;
  i.a = std::move(a);
  i.cov.cur = m_options.seqset->ctx_begin();
  m_active.emplace_back(std::move(i));
}

void split_variants::chunk_and_output() {
  advance_ref_coverage_range(m_leftmost_offset, m_rightmost_offset + m_options.seqset->max_read_len());

  if (k_split_var_debug) {
    std::cout << "Variant regions: " << m_variant_regions << "\n";
  }

  for (const auto& interval : m_variant_regions) {
    aoffset_t start_offset = interval.lower();
    aoffset_t limit_offset = interval.upper();
    CHECK_GT(limit_offset, start_offset) << "0-ref-sized regions should have been expanded\n";
    assembly_ptr ref_asm = make_ref_coverage_assembly(start_offset, limit_offset);
    CHECK(ref_asm) << start_offset << " to " << limit_offset << " on " << m_options.scaffold_name;
    sort_and_output(std::move(ref_asm));
  }

  for (asm_info& i : m_active) {
    if (k_split_var_debug) {
      std::cout << "Processing assembly: " << *i.a << "\n";
    }
    // Pad it so it has everything we need for this variant region.
    pad_assembly(i.a.get(), m_leftmost_offset, m_rightmost_offset, m_options);

    if (k_split_var_debug) {
      std::cout << "Padded assembly: " << *i.a << "\n";
    }
    // Calculate coverage
    i.cov.cur = m_options.seqset->ctx_begin();
    i.cov.offset = 0;
    for (dna_base b : i.a->seq) {
      advance_coverage(&i.cov, b);
    }
    CHECK_EQ(i.cov.offset, i.a->seq.size());

    output_variant_regions(i);
  }

  m_variant_regions.clear();
  m_active.clear();
}

void split_variants::output_variant_regions(const asm_info& i) {
  aoffset_t ref_offset = i.a->left_offset;
  aoffset_t seq_offset = 0;
  auto vit = i.a->aligned_variants.begin();
  auto covit = i.cov.coverage.begin();

  for (const auto& interval : m_variant_regions) {
    if (k_split_var_debug) {
      std::cout << "Interval " << interval << "; ref_offset = " << ref_offset
                << " seq_offset = " << seq_offset << "\n";
    }

    aoffset_t start_offset = interval.lower();
    aoffset_t limit_offset = interval.upper();
    CHECK_GT(limit_offset, start_offset);
    dna_sequence var_seq;

    int min_depth = std::numeric_limits<int>::max();
    auto update_depth_here = [&]() {
      while (covit != i.cov.coverage.end() && covit->offset < seq_offset) {
        ++covit;
      }

      if (covit != i.cov.coverage.end()) {
        min_depth = std::min(min_depth, covit->depth);
      }
    };
    auto advance_seq = [&](aoffset_t adv) {
      aoffset_t adv_remaining = adv;
      while (adv_remaining) {
        update_depth_here();

        ++seq_offset;
        --adv_remaining;
      }
      update_depth_here();
      dna_slice to_add = dna_slice(i.a->seq).subseq(seq_offset - adv, adv);
      var_seq += to_add;
      if (k_split_var_debug) {
        std::cout << "Advanced " << adv << " bases: +" << to_add << " = " << var_seq << "\n";
      }
    };

    if (vit != i.a->aligned_variants.end()) {
      CHECK_LE(start_offset, vit->left_offset);
    }
    CHECK_GE(start_offset, ref_offset);

    if (ref_offset < start_offset) {
      aoffset_t adv = (start_offset - ref_offset);
      ref_offset += adv;
      if (k_split_var_debug) {
        std::cout << "Advancing start +" << adv << " to " << ref_offset << "\n";
      }
      advance_seq(adv);
    }

    // Starting variant region.
    min_depth = std::numeric_limits<int>::max();
    var_seq = dna_sequence();

    aoffset_t seq_adv = 0;
    while (vit != i.a->aligned_variants.end() && vit->right_offset <= limit_offset) {
      if (k_split_var_debug) {
        std::cout << "Advancing up to " << *vit << " from ref_offset=" << ref_offset
                  << " seq_offset=" << seq_offset << ", seq_adv= " << seq_adv << "\n";
      }
      aoffset_t adv = vit->left_offset - ref_offset;
      CHECK_GE(adv, 0);
      ref_offset += adv;
      seq_adv += adv;

      CHECK_EQ(ref_offset, vit->left_offset);
      CHECK_LE(ref_offset, limit_offset);

      if (k_split_var_debug) {
        std::cout << "Advancing through " << *vit << " from ref_offset=" << ref_offset
                  << " seq_offset=" << seq_offset << ", seq_adv= " << seq_adv << "\n";
      }
      ref_offset += vit->right_offset - vit->left_offset;
      seq_adv += vit->seq.size();
      ++vit;

      CHECK_LE(ref_offset, limit_offset) << dump_assembly_and_vars(*i.a);
    }

    aoffset_t final_adv = limit_offset - ref_offset;
    if (k_split_var_debug) {
      std::cout << "Advancing +" << final_adv << " to " << limit_offset
                << " after variants from ref_offset=" << ref_offset << " seq_offset=" << seq_offset
                << ", seq_adv= " << seq_adv << "\n";
    }
    CHECK_GE(final_adv, 0);
    ref_offset += final_adv;
    seq_adv += final_adv;

    advance_seq(seq_adv);
    CHECK_EQ(ref_offset, limit_offset);

    assembly_ptr var_asm = make_unique<assembly>();
    var_asm->assembly_id = i.a->assembly_id;
    var_asm->min_overlap = i.a->min_overlap;
    var_asm->seq = var_seq;
    var_asm->left_offset = start_offset;
    var_asm->right_offset = limit_offset;

    if (k_split_var_debug) {
      std::cout << "Resultant variant assembly: " << *var_asm << "\n";
    }
    sort_and_output(std::move(var_asm));
  }
}

void split_variants::advance_ref_coverage_range(aoffset_t flush_to, aoffset_t target) {
  if (k_split_var_debug) {
    std::cout << "Advancing ref coverage to [" << flush_to << ", " << target << ") from "
              << m_ref.cov.offset << "\n";
  }

  if (m_ref.cov.offset < flush_to - int(m_options.seqset->max_read_len())) {
    // Skip ahead and don't bother calculating all the coverage in between
    m_ref.cov.coverage.clear();
    m_ref.cov.cur = m_options.seqset->ctx_begin();
    m_ref.scaffold_it.skip_to(flush_to - int(m_options.seqset->max_read_len()), "split_variants");
    m_ref.cov.offset = m_ref.scaffold_it.offset();
    if (k_split_var_debug) {
      std::cout << "Skipping ahead to " << m_ref.cov.offset << "\n";
    }
  }

  unsigned nflush = 0;
  while (!m_ref.cov.coverage.empty() && m_ref.cov.coverage.front().offset < flush_to) {
    m_ref.cov.coverage.pop_front();
    ++nflush;
  }
  if (k_split_var_debug) {
    std::cout << "Flushed " << nflush << " old ref cov entries\n";
  }

  while (m_ref.scaffold_it != m_options.scaffold->end() && m_ref.cov.offset < target) {
    if (m_ref.scaffold_it.first_in_extent()) {
      m_ref.cov.cur = m_options.seqset->ctx_begin();
      m_ref.cov.offset = m_ref.scaffold_it.offset();
      if (k_split_var_debug) {
        std::cout << "First in extent: " << m_ref.cov.offset << "\n";
      }
    } else {
      CHECK_EQ(m_ref.cov.offset, m_ref.scaffold_it.offset());
    }

    advance_coverage(&m_ref.cov, *m_ref.scaffold_it);
    ++m_ref.scaffold_it;
  }
  if (k_split_var_debug) {
    std::cout << "Ref coverage now up to " << m_ref.cov.offset << "\n";
  }
}

void split_variants::advance_coverage(coverage_state* cov, dna_base base) {
  coverage_entry new_cov;
  new_cov.offset = cov->offset;
  cov->coverage.push_back(new_cov);
  cov->cur = cov->cur.push_front_drop(base.complement());
  if (cov->cur.begin() + 1 == cov->cur.end()) {
    auto reads = m_options.readmap->entry_to_index(cov->cur.begin());
    for (uint32_t read_id = reads.first; read_id != reads.second; ++read_id) {
      int read_len = m_options.readmap->get_readlength(read_id);
      if (read_len > int(cov->cur.size())) {
        continue;
      }
      aoffset_t coffset = cov->offset;
      for (auto cit = cov->coverage.rbegin();
           cit != cov->coverage.rend() && cit->offset == coffset && read_len > 1;
           --coffset, ++cit, --read_len) {
        ++cit->depth;
      }
    }
  }
  ++cov->offset;
}

assembly_ptr split_variants::make_ref_coverage_assembly(aoffset_t left_offset,
                                                        aoffset_t right_offset) {
  if (k_split_var_debug) {
    std::cout << "Making ref assembly for [" << left_offset << ", " << right_offset << ")\n";
  }
  std::deque<coverage_entry>::iterator start_it =
      std::lower_bound(m_ref.cov.coverage.begin(), m_ref.cov.coverage.end(), left_offset);
  std::deque<coverage_entry>::iterator end_it =
      std::upper_bound(m_ref.cov.coverage.begin(), m_ref.cov.coverage.end(), right_offset);

  int min_depth = std::numeric_limits<int>::max();
  aoffset_t offset = left_offset;
  if (start_it == end_it) {
    min_depth = 0;
  } else {
    for (auto it = start_it; it != end_it; ++it) {
      if (it->offset != offset) {
        min_depth = 0;
        break;
      }
      min_depth = std::min(min_depth, it->depth);
      ++it;
      ++offset;
    }
  }

  CHECK_NE(min_depth, std::numeric_limits<int>::max());

  assembly_ptr a = make_unique<assembly>();
  a->assembly_id = 0;
  a->left_offset = left_offset;
  a->right_offset = right_offset;
  a->other_depth = min_depth;
  a->matches_reference = true;

  scaffold s = m_options.scaffold->subscaffold(left_offset, right_offset - left_offset);
  if (!s.is_simple()) {
    // TODO(nils): Process missing sections of reference properly.  We
    // can't just set matches_reference, since that prohibits missing
    // sections in the middle.
    return nullptr;
  }
  dna_slice seq = s.get_simple();
  a->seq = dna_sequence(seq.begin(), seq.end());

  return a;
}

}  // namespace variants
