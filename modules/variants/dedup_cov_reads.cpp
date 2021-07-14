#include "modules/variants/dedup_cov_reads.h"

#include "modules/bio_base/readmap.h"
#include "modules/bio_base/seqset.h"

namespace variants {

static constexpr bool k_debug_dedup = false;

dedup_cov_reads::dedup_cov_reads(const assemble_options& opts, pipeline_step_t output)
    : sorted_output_pipeline_step(std::move(output)), m_options(opts) {}

dedup_cov_reads::~dedup_cov_reads() { flush(); }

void dedup_cov_reads::advance_to(aoffset_t offset) {
  while (m_cur_offset < offset) {
    if (m_active.empty()) {
      m_cur_offset = offset;
      return;
    }
    m_cur_offset = std::min<aoffset_t>(m_active.begin()->first, offset);
    while (!m_active.empty()) {
      if (m_active.begin()->first > m_cur_offset) {
        break;
      }
      CHECK_EQ(m_active.begin()->first, m_cur_offset);
      auto& a = m_active.begin()->second;

      track_reads(a.get(), false /* remove tracking */);
      aoffset_t left_offset = a->left_offset;
      sort_and_output(std::move(a));
      m_active.erase(m_active.begin());
      untrack_left_offset(left_offset);
    }
  }
}

void dedup_cov_reads::track_reads(assembly* a, bool track) {
  if (!a->edge_coverage) {
    throw(
        io_exception("dedup_cov_reads receieved an assembly without edge coverage; edge coverage "
                     "must be present for dedup_cov_reads"));
  }

  edge_coverage_t& ec = *a->edge_coverage;

  if (a->matches_reference) {
    track_ref(a, &ec.interior, track);
  } else {
    track_ref(a, &ec.reference_start, track);
    track_ref(a, &ec.reference_end, track);

    track_var(a, &ec.variant_end, track);
    track_var(a, &ec.variant_start, track);
    track_var(a, &ec.interior, track);
  }
}

void dedup_cov_reads::remove_read_if_seen_in_var(assembly* because_of, uint32_t read_id) {
  const bool debug_dedup = assembly_needs_trace(*because_of) || k_debug_dedup;
  auto eq_range = m_seen_var_reads.equal_range(read_id);
  auto var_next = eq_range.first;
  for (auto var_it = eq_range.first; var_it != eq_range.second; var_it = var_next) {
    ++var_next;
    assembly* var_a = var_it->second.a;
    if (var_a == because_of) {
      continue;
    }
    if (debug_dedup) {
      std::cout << "dedup_cov_reads: Deduplicating read " << read_id << " on both " << *because_of
                << " and " << *var_it->second.a << "\n";
      std::cout.flush();
    }
    var_it->second.collection->erase(read_id);
    m_seen_var_reads.erase(var_it);
  }
}

void dedup_cov_reads::track_ref(assembly* a, read_id_set* reads, bool track) {
  auto next = reads->begin();
  for (auto it = reads->begin(); it != reads->end(); it = next) {
    ++next;
    uint32_t read_id = *it;

    if (track) {
      auto ins_result = m_seen_ref_reads.emplace(read_id, 0);
      unsigned& count = ins_result.first->second;
      if (!count) {
        // Only check to see if we should remove it the first time through.
        remove_read_if_seen_in_var(a, read_id);
      }
      ++count;
    } else {
      auto it = m_seen_ref_reads.find(read_id);
      CHECK(it != m_seen_ref_reads.end());
      CHECK_GT(it->second, 0);
      --it->second;
      if (it->second == 0) {
        m_seen_ref_reads.erase(it);
      }
    }
  }
}

void dedup_cov_reads::track_var(assembly* a, read_id_set* reads, bool track) {
  const bool debug_dedup = assembly_needs_trace(*a) || k_debug_dedup;
  auto next = reads->begin();
  read_id_set to_erase_reads;
  for (auto it = reads->begin(); it != reads->end(); it = next) {
    ++next;
    uint32_t read_id = *it;

    if (track) {
      remove_read_if_seen_in_var(a, read_id);

      auto ref_it = m_seen_ref_reads.find(read_id);
      if (ref_it != m_seen_ref_reads.end()) {
        if (debug_dedup) {
          std::cout << "dedup_cov_reads: Deduplicating read " << read_id
                    << " both in reference and in " << *a << "\n";
          std::cout.flush();
        }

        to_erase_reads.insert(read_id);
        continue;
      }

      var_seen_t new_seen;
      new_seen.collection = reads;
      new_seen.a = a;
      m_seen_var_reads.emplace(read_id, new_seen);
    } else {
      auto eq_range = m_seen_var_reads.equal_range(read_id);
      bool found_eq = false;
      for (auto eq_it = eq_range.first; eq_it != eq_range.second; ++eq_it) {
        if (eq_it->second.a != a) {
          continue;
        }
        if (eq_it->second.collection != reads) {
          continue;
        }
        found_eq = true;

        m_seen_var_reads.erase(eq_it);
        break;
      }
      CHECK(found_eq) << "Missing tracking for " << read_id << "?";
    }
  }
  for (uint32_t read_id : to_erase_reads) {
    reads->erase(read_id);
  }
}

void dedup_cov_reads::on_assembly(assembly_ptr a) {
  track_left_offset(a->left_offset);
  advance_to(a->left_offset);
  track_reads(a.get(), true /* add tracking */);
  aoffset_t right_offset = a->right_offset;
  m_active.emplace(right_offset, std::move(a));
}

void dedup_cov_reads::flush() {
  advance_to(std::numeric_limits<aoffset_t>::max());
  CHECK(m_active.empty());
  CHECK(m_seen_ref_reads.empty());
  CHECK(m_seen_var_reads.empty());
  flush_sorted();
}

}  // namespace variants
