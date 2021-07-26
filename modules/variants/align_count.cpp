#include "modules/variants/align_count.h"

namespace variants {
align_count::align_count(const assemble_options& opts, pipeline_step_t output)
    : apply_edges_step(std::move(output)), m_opts(opts) {}

align_count::~align_count() { CHECK(m_active_counts.empty()); }

void align_count::on_assembly_edges(optional_aoffset reference_pos,
                                    const std::vector<assembly_ptr>& left_edges,
                                    const std::vector<assembly_ptr>& inserts,
                                    const std::vector<assembly_ptr>& right_edges) {
  for (const auto& a : left_edges) {
    end_assembly(a.get());
  }
  for (const auto& a : inserts) {
    start_assembly(a.get());
  }
  for (const auto& a : inserts) {
    end_assembly(a.get());
  }
  for (const auto& a : right_edges) {
    start_assembly(a.get());
  }
}

void align_count::start_assembly(assembly* a) {
  auto res = m_active.try_emplace(a);
  CHECK(res.second) << "Duplicate assembly?";
  auto it = res.first;

  auto& act = it->second;

  add_coverage(a, &act, true /* add */);
}

void align_count::add_coverage(assembly* a, active_assembly* act, bool add) {
  aoffset_t seqlen = a->seq.size();

  for (const auto& cov_entry : a->read_coverage->reads()) {
    aoffset_t aligned_bases;

    if (cov_entry.offset < 0) {
      if (cov_entry.offset + cov_entry.read_len > seqlen) {
        aligned_bases = seqlen;
      } else {
        aligned_bases = cov_entry.offset + cov_entry.read_len;
      }
    } else {
      if (cov_entry.offset + cov_entry.read_len > seqlen) {
        aligned_bases = seqlen - cov_entry.offset;
      } else {
        aligned_bases = cov_entry.read_len;
      }
    }

    auto count_it = m_active_counts.begin();

    for (uint32_t read_id : cov_entry.read_ids) {
      if (add) {
        if (!act->all_read_ids.contains(read_id)) {
          act->all_read_ids.insert(read_id);
          act->counts.local_read_lens += cov_entry.read_len;
          act->counts.local_aligned_bases += aligned_bases;
        }

        count_it = m_active_counts.try_emplace(count_it, read_id);
        // Store the aligned count so far
        act->counts.tot_aligned_bases += count_it->second;
        count_it->second += aligned_bases;

        // Save this alignment in any overlapping assemblies we've already seen.
        for (auto& other_act_ent : m_active) {
          active_assembly* other_act = &other_act_ent.second;
          if (other_act == act) {
            continue;
          }
          if (other_act->all_read_ids.contains(read_id)) {
            other_act->counts.tot_aligned_bases += aligned_bases;
          }
        }

        act->counts.tot_aligned_bases += aligned_bases;
      } else {
        auto count_it = m_active_counts.find(read_id);
        CHECK(count_it != m_active_counts.end());
        CHECK_GE(count_it->second, aligned_bases);
        count_it->second -= aligned_bases;
        if (count_it->second == 0) {
          m_active_counts.erase(count_it);
        }
      }
    }
  }
}

void align_count::end_assembly(assembly* a) {
  auto act_it = m_active.find(a);
  CHECK(act_it != m_active.end());
  auto* act = &act_it->second;

  a->align_count.emplace(act->counts);

  add_coverage(a, act, false /* remove */);

  m_active.erase(act_it);
}

}  // namespace variants
