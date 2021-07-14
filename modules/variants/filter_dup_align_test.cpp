#include "modules/variants/filter_dup_align.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "modules/bio_base/dna_testutil.h"
#include "modules/variants/add_ref.h"
#include "modules/variants/assemble_testutil.h"
#include "modules/variants/read_cov.h"

namespace variants {

using namespace testing;
using namespace dna_testutil;

class filter_dup_align_test : public assemble_test {
 public:
  struct pending_asm_t {
    assembly_ptr a;
    read_coverage_set reads;
  };
  void add_assembly(size_t priority, aoffset_t left_offset, aoffset_t right_offset) {
    assembly_ptr a = make_unique<assembly>();
    a->assembly_id = priority;
    a->tags.insert("filter_dup_align_test");
    a->left_offset = left_offset;
    a->right_offset = right_offset;
    a->seq = tseq(std::to_string(a->assembly_id));
    a->read_coverage.emplace();

    CHECK(!m_inputs.count(priority));
    m_inputs[priority].a = std::move(a);
  }

  void add_read(uint32_t read_id, std::vector<size_t> asm_ids) {
    int read_len = 0;
    std::vector<pending_asm_t*> pending;
    for (const auto& id : asm_ids) {
      CHECK(m_inputs.count(id));
      pending.push_back(&m_inputs[id]);
    }
    for (const auto* p : pending) {
      read_len += p->a->seq.size();
    }

    aoffset_t read_pos = 0;
    for (auto* p : pending) {
      p->reads.insert(read_pos, read_id, read_len);
      read_pos -= aoffset_t(p->a->seq.size());
    }
    CHECK_EQ(read_len, -read_pos);
  }

  std::map<size_t /* priority id */, pending_asm_t> m_inputs;

  static std::vector<assembly_ptr> sort_by_priority(std::vector<assembly_ptr> asms) {
    std::sort(asms.begin(), asms.end(), [](const assembly_ptr& a, const assembly_ptr& b) {
      return a->assembly_id < b->assembly_id;
    });
    return asms;
  }

  void run() {
    std::vector<assembly_ptr> asms;
    for (auto& input : m_inputs) {
      input.second.a->read_coverage.emplace(
          input.second.reads.build_and_clear(input.second.a->seq.size()));
      asms.push_back(input.second.a);
    }
    std::sort(asms.begin(), asms.end(), canon_assembly_order());

    std::unique_ptr<filter_dup_align> pipeline =
        make_unique<filter_dup_align>(sort_by_priority, test_output());
    for (auto& a : asms) {
      pipeline->add(std::move(a));
    }
    pipeline.reset();

    save_found_reads();
  }

  void save_found_reads() {
    for (const auto& a : m_assemblies) {
      for (const auto& cov : a.read_coverage->reads()) {
        for (uint32_t read_id : cov.read_ids) {
          m_filtered_reads[read_id].insert(a.assembly_id);
        }
      }
    }
  }

  std::map<uint32_t /* read id */, std::set<size_t> /* assembly ids */> m_filtered_reads;
};

TEST_F(filter_dup_align_test, singles) {
  add_assembly(1, 10, 20);
  add_assembly(2, 20, 30);
  add_assembly(3, 30, 40);
  add_read(1, {1});
  add_read(1, {2});
  add_read(1, {3});
  run();

  EXPECT_THAT(m_filtered_reads, ElementsAre(Pair(1, ElementsAre(1, 2, 3))));
}

TEST_F(filter_dup_align_test, conflict) {
  // Block 1
  add_assembly(1, 5, 10);
  // Block 2
  add_assembly(3, 10, 20);
  add_assembly(2, 20, 30);
  add_assembly(4, 20, 30);
  add_assembly(5, 30, 40);
  // Block 3
  add_assembly(6, 40, 45);

  // Read 1 has a conflict in block 2
  add_read(1, {3, 2, 5});  // Higher priority path
  add_read(1, {3, 4, 5});  // Lower priority path
  // But no conflicts in block 1 or 3
  add_read(1, {1});  // Should not be part of same block
  add_read(1, {6});  // Should not be part of same block.

  // Read 2 goes through all the blocks, and should still get placed.
  add_read(2, {1, 3, 4, 5, 6});
  run();

  EXPECT_THAT(           //
      m_filtered_reads,  //
      ElementsAre(       //
          Pair(1, ElementsAre(
                      // Placed in first block:
                      1,
                      // Placed in second block:
                      2, 3, 5,
                      // Placed in third block:
                      6)),
          Pair(2, ElementsAre(
                      // One placement spanning all three blocks.
                      // First block:
                      1,
                      // Second block:
                      3, 4, 5,
                      // Third block:
                      6))));
}

}  // namespace variants
