#include "modules/bio_base/dna_testutil.h"
#include "modules/bio_base/seqset_testutil.h"
#include "modules/bio_base/seqset_export.h"
#include "modules/io/config.h"

#include <gtest/gtest.h>
#include <gmock/gmock.h>

namespace {

using namespace dna_testutil;
using namespace testing;

}  // namespace

class seqset_export_test : public Test {
 public:
  void do_export(const std::vector<std::pair<dna_sequence, dna_sequence>>& paired_reads,
                 const std::vector<dna_sequence>& unpaired_reads);

  std::mutex m_mu;
  std::shared_ptr<seqset_file> m_seqset;
  std::unique_ptr<readmap> m_readmap;
  boost::optional<seqset_export> m_export;

  std::vector<std::pair<dna_sequence, dna_sequence>> m_exported_paired;
  std::vector<dna_sequence> m_exported_unpaired;
};

class test_export_worker : public seqset_export_worker {
 public:
  test_export_worker(seqset_export_test* test) : m_test(test) {}

  void output_paired(uint32_t read_id, dna_slice r1, dna_slice r2) override {
    std::lock_guard<std::mutex> l(m_test->m_mu);
    m_test->m_exported_paired.emplace_back(dna_sequence(r1.begin(), r1.end()),
                                           dna_sequence(r2.begin(), r2.end()));
  }
  void output_unpaired(uint32_t read_id, dna_slice r1) override {
    std::lock_guard<std::mutex> l(m_test->m_mu);
    m_test->m_exported_unpaired.emplace_back(dna_sequence(r1.begin(), r1.end()));
  }

 private:
  seqset_export_test* m_test = nullptr;
};

void seqset_export_test::do_export(
    const std::vector<std::pair<dna_sequence, dna_sequence>>& paired_reads,
    const std::vector<dna_sequence>& unpaired_reads) {
  std::vector<dna_sequence> all_reads = unpaired_reads;
  for (const auto& r : paired_reads) {
    all_reads.push_back(r.first);
    all_reads.push_back(r.second);
  }
  m_seqset = seqset_for_reads(all_reads);
  m_readmap = readmap_for_reads(m_seqset, paired_reads, unpaired_reads);

  m_export.emplace(m_seqset.get(), m_readmap.get(), CONF_S(storage_root));
  m_export->prepare();

  m_export->write_paired([&]() { return make_unique<test_export_worker>(this); });
  m_export->write_unpaired([&]() { return make_unique<test_export_worker>(this); });
}

MATCHER(IsPair, "") {
  auto expected = get<0>(arg);
  auto actual = get<1>(arg);

  if (expected == actual) {
    return true;
  }

  if (expected == std::make_pair(actual.second, actual.first)) {
    return true;
  }

  return false;
}

TEST_F(seqset_export_test, simple_unpaired) {
  std::vector<dna_sequence> unpaired = {tseq("a")};
  do_export({}, unpaired);

  EXPECT_THAT(m_exported_paired, IsEmpty());
  EXPECT_THAT(m_exported_unpaired, UnorderedElementsAreArray(unpaired));
}

TEST_F(seqset_export_test, simple_paired) {
  std::vector<std::pair<dna_sequence, dna_sequence>> paired = {
      {tseq("a"), tseq("b")}};
  do_export(paired, {});

  EXPECT_THAT(m_exported_unpaired, IsEmpty());
  EXPECT_THAT(m_exported_paired, UnorderedPointwise(IsPair(), paired));
}

TEST_F(seqset_export_test, multiple) {
  dna_sequence palindrome1 = tseq("a") + tseq_rc("a");
  dna_sequence palindrome2 = tseq("b") + tseq_rc("b");
  dna_sequence palindrome3 = tseq("c") + tseq_rc("c");

  std::vector<std::pair<dna_sequence, dna_sequence>> paired = {
      {palindrome1, tseq("d")},
      {tseq_rc("d"), palindrome2},
      {tseq("a"), tseq("xyz")}};

  std::vector<dna_sequence> unpaired = {palindrome1.subseq(0, 10), palindrome3};
  do_export(paired, unpaired);

  EXPECT_THAT(m_exported_unpaired, UnorderedElementsAreArray(unpaired));
  EXPECT_THAT(m_exported_paired, UnorderedPointwise(IsPair(), paired));
}
