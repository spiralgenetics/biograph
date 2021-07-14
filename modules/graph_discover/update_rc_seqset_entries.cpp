#include "modules/graph_discover/update_rc_seqset_entries.h"

static constexpr bool k_dbg = false;

namespace variants {

update_rc_seqset_entries::update_rc_seqset_entries(const assemble_options &options,
                                                   pipeline_step_t output)
    : apply_edges_step(std::move(output)), m_options(options) {
  CHECK(m_options.seqset);
}

update_rc_seqset_entries::~update_rc_seqset_entries() {
  flush();
  CHECK(m_self_test_succeeded) << "Self test result should have been checked.";
}

void update_rc_seqset_entries::on_assembly_edges(optional_aoffset reference_pos,
                                                 const std::vector<assembly_ptr> &left_edges,
                                                 const std::vector<assembly_ptr> &inserts,
                                                 const std::vector<assembly_ptr> &right_edges) {
  seqset_range_set incoming;

  if (k_dbg) {
    std::cerr << this << " update rc seqset entries, reference_pos=" << reference_pos << "\n";
    std::cerr << "Left edges:\n";
    for (const auto &a : left_edges) {
      std::cerr << " " << *a << "\n";
    }
    std::cerr << "Inserts:\n";
    for (const auto &a : inserts) {
      std::cerr << " " << *a << "\n";
    }
    std::cerr << "Right edges:\n";
    for (const auto &a : right_edges) {
      std::cerr << " " << *a << "\n";
    }
  }

  for (const auto &a : left_edges) {
    CHECK(!a->rc_seqset_entries.empty())
        << "Should have already generated rc_seqset_entries for " << *a;

    const auto &entries = a->rc_seqset_entries.starts();
    incoming.insert(entries.begin(), entries.end());
  }

  if (incoming.empty()) {
    incoming.insert(m_options.seqset->ctx_begin());
  }

  if (!inserts.empty()) {
    propagate(incoming, inserts);

    for (const auto &a : inserts) {
      CHECK(!a->rc_seqset_entries.empty())
          << "Should have already generated rc_seqset_entries for " << *a;

      const auto &entries = a->rc_seqset_entries.starts();
      incoming.insert(entries.begin(), entries.end());
    }
  }

  propagate(incoming, right_edges);
}

void update_rc_seqset_entries::propagate(const seqset_range_set &incoming,
                                         const std::vector<assembly_ptr> &targets) {
  for (const auto &target : targets) {
    target->rc_seqset_entries.propagate_from_end(incoming, dna_slice(target->seq).rev_comp(),
                                                 m_options);
  }
}

}  // namespace variants
