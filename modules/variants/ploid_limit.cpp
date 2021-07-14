#include "modules/variants/ploid_limit.h"

#include "modules/variants/scaffold.h"

#include <boost/icl/interval_set.hpp>
#include <boost/icl/split_interval_map.hpp>

namespace variants {

constexpr int k_ploid_debug = 0;

ploid_limiter::~ploid_limiter() { ploid_flush(); }

void ploid_limiter::on_assembly(assembly_ptr a) {
  if (k_ploid_debug) {
    std::cout << "Got assembly: " << *a << "\n";
  }
  CHECK_GE(a->left_offset, m_cur_offset);
  m_cur_offset = a->left_offset;

  output_active();
  if (m_var_active == 0 && !m_deploid_scores.empty()) {
    do_deploid();
  }
  flush_queued();

  if (!a->matches_reference) {
    ++m_var_active;
  }
  track_left_offset(a->left_offset);
  aoffset_t right_offset = a->right_offset;
  m_active.emplace(std::make_pair(right_offset, std::move(a)));

  if (k_ploid_debug) {
    std::cout << "New active: " << m_active.size() << " vars: " << m_var_active
              << " deploids: " << m_deploid_scores.size() << "\n";
  }
}

void ploid_limiter::ploid_flush() {
  if (k_ploid_debug) {
    std::cout << "Ploid flush\n";
  }
  m_cur_offset = std::numeric_limits<aoffset_t>::max();
  output_active();
  CHECK_EQ(0, m_var_active);
  if (!m_deploid_scores.empty()) {
    do_deploid();
  }
  flush_queued();

  CHECK(m_active.empty());
  CHECK(m_deploid_scores.empty());
}

void ploid_limiter::output_active() {
  while (!m_active.empty() && m_active.begin()->first <= m_cur_offset) {
    assembly_ptr a = std::move(m_active.begin()->second);
    auto score = a->score;
    m_active.erase(m_active.begin());

    if (!a->matches_reference) {
      if (k_ploid_debug) {
        std::cout << "Adding var to deploid: " << *a << "\n";
      }
      CHECK_GT(m_var_active, 0);
      --m_var_active;
      m_deploid_scores.emplace(score, std::move(a));
      continue;
    }

    if (m_var_active || !m_deploid_scores.empty()) {
      // If we have any non-reference variants, ploid-limit them
      // against this reference assembly before continuing
      if (k_ploid_debug) {
        std::cout << "Adding ref to deploid: " << *a << "\n";
      }
      m_deploid_scores.emplace(score, std::move(a));
      continue;
    }

    untrack_left_offset(a->left_offset);
    if (k_ploid_debug) {
      std::cout << "Outputting: " << *a << "\n";
    }
    sort_and_output(std::move(a));
  }
}

void ploid_limiter::flush_queued() {
  if (m_var_active || !m_deploid_scores.empty()) {
    // Cant flush until we're finish deploiding.
    return;
  }

  flush_sorted_to(m_cur_offset);
}

void ploid_limiter::do_deploid() {
  // Include all currently active reference in the deplodification process.
  for (auto it = m_active.begin(); it != m_active.end(); ++it) {
    assembly_ptr a = std::move(it->second);
    auto score = a->score;

    CHECK(a->matches_reference);
    m_deploid_scores.emplace(std::make_pair(score, std::move(a)));
  }
  m_active.clear();
  CHECK_EQ(m_var_active, 0);

  std::vector<assembly_ptr> allele_assemblies;
  allele_assemblies.reserve(m_deploid_scores.size());

  if (k_ploid_debug) {
    std::cout << "Starting deploid at " << m_cur_offset << "\n";
  }
  // Go through assemblies, best score to worst score.
  for (auto it = m_deploid_scores.rbegin(); it != m_deploid_scores.rend(); ++it) {
    assembly_ptr a = std::move(it->second);
    if (k_ploid_debug) {
      std::cout << "Considering for deploid: " << dump_assembly_and_vars(*a) << "\n";
    }
    std::vector<size_t> merged_in;
    boost::icl::interval_map<aoffset_t, size_t> conflicts;
    unsigned max_conflict_count = 0;
    for (size_t merge_idx = 0; merge_idx != allele_assemblies.size(); ++merge_idx) {
      if (!allele_assemblies[merge_idx]) {
        continue;
      }
      auto& merge_with = allele_assemblies[merge_idx];
      if (merge_with->left_offset >= a->right_offset ||
          merge_with->right_offset <= a->left_offset) {
        continue;
      }
      CHECK(a);
      if (k_ploid_debug > 1) {
        std::cout << "Attempting to merge " << *a << " with " << *merge_with << "\n";
      }

      assembly_ptr merged = merge_assemblies(*a, *merge_with);
      if (k_ploid_debug > 1) {
        if (merged) {
          std::cout << "Merge successful: " << *merged << "\n";
        } else {
          std::cout << "Merge unsuccessful\n";
        }
      }

      if (merged) {
        merged_in.push_back(merge_idx);
        untrack_left_offset(a->left_offset);
        a = std::move(merged);
        track_left_offset(a->left_offset);
      } else {
        if (k_ploid_debug) {
          std::cout << "Conflicts with: " << dump_assembly_and_vars(*merge_with) << "\n";
        }
        auto conflict_i = boost::icl::interval<aoffset_t>::open(merge_with->left_offset,
                                                                merge_with->right_offset);
        conflicts += std::make_pair(conflict_i, size_t(1));

        for (const auto& c : (conflicts & conflict_i)) {
          if (c.second > max_conflict_count) {
            max_conflict_count = c.second;
          }
          if (max_conflict_count >= m_max_ploids) {
            break;
          }
        }
      }
      CHECK(a);
      if (max_conflict_count >= m_max_ploids) {
        break;
      }
    }
    CHECK(a);

    if (max_conflict_count >= m_max_ploids) {
      if (k_ploid_debug) {
        std::cout << "Too many conflicts; discarding: " << dump_assembly_and_vars(*a) << "\n";
      }
      untrack_left_offset(a->left_offset);
      continue;
    }

    size_t new_idx = allele_assemblies.size();
    if (k_ploid_debug) {
      std::cout << "Ploid saving " << *a << " at idx " << new_idx << "\n";
    }
    allele_assemblies.push_back(std::move(a));
    for (const size_t merge : merged_in) {
      CHECK(allele_assemblies[merge]);
      auto& old_a = allele_assemblies[merge];
      untrack_left_offset(old_a->left_offset);
      allele_assemblies[merge].release_and_discard();
    }
  }

  if (k_ploid_debug) {
    std::cout << "Done deploid\n";
  }

  for (auto& a : allele_assemblies) {
    if (!a) {
      continue;
    }
    auto right_offset = a->right_offset;
    if (right_offset <= m_cur_offset || !a->matches_reference) {
      // If it extends past cur but doesn't entirely match reference,
      // it got merged with some reference assemblies.  So the portion
      // that's still active is reference only.
      //
      // TODO(nils): It isn't quite right to remove from active and
      // output; it should really stay around so that it conflicts
      // with any compound hetrozygous variants.
      if (k_ploid_debug) {
        std::cout << "deploid output: " << *a << "\n";
      }
      untrack_left_offset(a->left_offset);
      sort_and_output(std::move(a));
    } else {
      if (k_ploid_debug) {
        std::cout << "deploid return to active: " << *a << "\n";
      }
      CHECK(a->matches_reference) << "cur offset=" << m_cur_offset << " assembly: " << *a;
      m_active.emplace(right_offset, std::move(a));
    }
  }
  m_deploid_scores.clear();
  CHECK_EQ(0, m_var_active);

  if (k_ploid_debug) {
    std::cout << "After deploiding, active:\n";
    for (const auto& i : m_active) {
      std::cout << "  " << *i.second << "\n";
      CHECK(i.second->matches_reference) << *i.second;
    }
  }
}

}  // namespace variants
