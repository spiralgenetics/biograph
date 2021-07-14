#include "modules/graph_discover/branch.h"

#include <boost/range/adaptor/reversed.hpp>
#include <chrono>
#include <ctime>
#include <random>

#include "modules/bio_base/readmap.h"
#include "modules/io/parallel.h"

namespace variants {

constexpr bool k_dbg = false;

branch_discover::branch_discover(const assemble_options& options, const std::string& tag,
                                 pipeline_step_t output)
    : graph_discover(options, std::move(output)), m_tag(tag) {}
branch_discover::~branch_discover() {}

void branch_discover::on_trace(const graph_discover::active_assembly* act) {
  if (!act->a->left_offset) {
    // We must have at least one anchor to trace from here.
    return;
  }

  aoffset_t min_overlap = opts().min_overlap;
  bool dbg = k_dbg;

  if (dbg) {
    std::cerr << "Branch tracing assembly " << *act->a << ", min overlap = " << min_overlap << "\n";
  }

  seqset_range_set cur = act->a->rc_seqset_entries.ends();

  dna_slice seq = dna_slice(act->a->seq);
  aoffset_t offset = 0;

  for (dna_base ref_b : seq) {
    dna_base rc_ref_b = ref_b.complement();

    for (dna_base b : dna_bases()) {
      dna_base rc_b = b.complement();
      if (rc_b == rc_ref_b) {
        continue;
      }

      bool extend_here = false;
      for (const auto& r : cur) {
        seqset_range next_r = r.push_front_drop(rc_b, min_overlap);
        if (!next_r.valid()) {
          continue;
        }

        extend_here = true;

        if (dbg) {
          std::cerr << "Extending at offset " << offset << "/" << seq.size() << " from "
                    << r.sequence() << " to " << next_r.sequence() << "\n";
        }
      }
      if (extend_here) {
        dna_sequence extension;
        extension.push_back(b);
        sort_and_output(discover_extend_right(act, offset, extension, m_tag));
        break;
      }
    }

    seqset_range_set new_cur;
    for (const auto& r : cur) {
      new_cur.insert(r.push_front_drop(rc_ref_b));
    }
    cur = std::move(new_cur);
    seqset_set_dedup_prefixes(cur);
    ++offset;
  }

  CHECK(cur == act->a->rc_seqset_entries.starts())
      << " cur=" << cur << " act: " << act->a->rc_seqset_entries.starts();
}

}  // namespace variants
