#include "modules/bio_mapred/kmer_set.h"
#include "modules/bio_mapred/align_kmer.h"
#include "modules/bio_mapred/kmer_filter_mapper.h"
#include "modules/bio_mapred/kmerize_reads_mapper.h"
#include "modules/io/config.h"
#include "modules/mapred/map_reduce_task.h"
#include "modules/mapred/task_mgr.h"
#include "modules/test/fastq_test_utils.h"
#include "modules/test/test_utils.h"

#include <gtest/gtest.h>
#include <random>

using namespace testing;

class cmp_set_filler {
 public:
  cmp_set_filler(std::set<kmer_t>& set) : m_set(set) {}
  void operator()(size_t index, const kmer_t& k, size_t ks, const std::string& v) {
    if (index < 10) {
      printf("%ld: %s\n", index, dna_sequence(k, ks).as_string().c_str());
    }
    m_set.insert(k);
  }

 private:
  std::set<kmer_t>& m_set;
};

class kmer_set_test : public TestWithParam<unsigned /* kmer size */> {
 public:
  void SetUp() override {
    m_kmer_size = GetParam();
    populate_kmers();
  }

  void populate_kmers() {
    task_mgr_local tm;
    path out_path(make_path("align_kmer_test"));

    // Get ecoli data
    manifest e_coli_reads;
    make_fastq_kv("golden/e_coli_10000snp.fq", make_path("e_coli_10000.kvp"));
    e_coli_reads.add(file_info(path(make_path("e_coli_10000.kvp")), 1017780, 10000), 0);

    // Kmerize reads
    kmerize_reads_params kp;
    json_deserialize(kp, R"|(
			{
				"kmer_size" : )|" +
                             std::to_string(m_kmer_size) + R"|(,
				"trim" : 0,
				"use_score" : false
			}
		)|");

    std::unique_ptr<map_reduce_task> t = make_unique<map_reduce_task>();
    t->input = e_coli_reads;
    t->map = "kmerize_reads";
    t->map_param = json_serialize(kp);
    t->sort = "lexical";
    t->reduce = "kcount";
    t->is_summary = true;
    t->use_sort = true;
    manifest kmers;
    tm.run_task(kmers, out_path.append("kmers"), std::move(t));

    manifest_reader mr(kmers);
    printf("Reading\n");
    m_ks.emplace(mr, kmers.get_num_records(), kp.kmer_size, cmp_set_filler(m_cmp_set));
  }

 protected:
  unsigned m_kmer_size;
  boost::optional<kmer_set> m_ks;
  std::set<kmer_t> m_cmp_set;
};

TEST_P(kmer_set_test, basic) {
  kmer_set::const_iterator it1 = m_ks->begin();
  std::set<kmer_t>::const_iterator it2 = m_cmp_set.begin();
  size_t count = 0;
  size_t kmer_size = m_ks->kmer_size();
  ASSERT_GT(m_ks->size(), size_t(1000));
  ASSERT_EQ(m_ks->size(), m_cmp_set.size());
  while (it1 != m_ks->end() && it2 != m_cmp_set.end()) {
    if (*it1 != *it2) {
      throw io_exception(printstring("Mismatch: %s vs %s",
                                     dna_sequence(*it1, kmer_size).as_string().c_str(),
                                     dna_sequence(*it2, kmer_size).as_string().c_str()));
    }
    kmer_set::const_iterator it1x = m_ks->find(*it2);
    std::set<kmer_t>::const_iterator it2x = m_cmp_set.find(*it1);
    ASSERT_EQ(it1, it1x);
    ASSERT_EQ(it2, it2x);
    it1++;
    it2++;
    count++;
  }
  printf("Count = %ld\n", count);
}

// Tests populating a kmer set from an estimated size using a lambda.
TEST_P(kmer_set_test, from_kmer_source) {
  std::vector<kmer_t> unordered;
  std::random_device rand_dev;
  std::mt19937 rand_source(rand_dev());
  std::copy(m_cmp_set.begin(), m_cmp_set.end(), std::back_inserter(unordered));
  std::shuffle(unordered.begin(), unordered.end(), rand_source);

  for (int size_estimate_off : {0, 1, 10000}) {
    SCOPED_TRACE(size_estimate_off);
    size_t count = 0;
    kmer_set ks2(unordered.size() + size_estimate_off, m_kmer_size, get_maximum_mem_bytes(),
                 [&](const std::function<void(kmer_t, unsigned /* flags */)>& output_kmer,
                     progress_handler_t progress) {
                   for (auto it = unordered.begin(); it != unordered.end(); ++it) {
                     output_kmer(*it, 0 /* flags */);
                   }
                 });

    kmer_set::const_iterator it1 = ks2.begin();
    std::set<kmer_t>::const_iterator it2 = m_cmp_set.begin();
    ASSERT_EQ(ks2.size(), m_cmp_set.size());
    while (it1 != ks2.end() && it2 != m_cmp_set.end()) {
      dna_sequence s2(*it2, m_kmer_size);
      dna_sequence s1(*it1, m_kmer_size);
      EXPECT_EQ(s1.as_string(), s2.as_string()) << count;
      kmer_set::const_iterator it1x = ks2.find(*it2);
      std::set<kmer_t>::const_iterator it2x = m_cmp_set.find(*it1);
      EXPECT_EQ(*it1, *it1x) << count;
      EXPECT_EQ(*it2, *it2x) << count;
      it1++;
      it2++;
      count++;
    }
    EXPECT_TRUE(it1 == ks2.end());
    EXPECT_TRUE(it2 == m_cmp_set.end());
  }
}

TEST_P(kmer_set_test, flags) {
  std::vector<std::pair<kmer_t, unsigned>> unordered;
  std::random_device rand_dev;
  std::mt19937 rand_source(rand_dev());
  std::vector<unsigned> ordered_flags;
  for (kmer_t k : m_cmp_set) {
    unsigned flag_val = std::uniform_int_distribution<unsigned>(0, 3)(rand_source);
    unordered.emplace_back(k, flag_val);
    ordered_flags.push_back(flag_val);
  }
  std::shuffle(unordered.begin(), unordered.end(), rand_source);
  CHECK_EQ(unordered.size(), ordered_flags.size());

  size_t count = 0;
  kmer_set ks2(unordered.size(), m_kmer_size, get_maximum_mem_bytes(),
               [&](const std::function<void(kmer_t, unsigned /* flags */)>& output_kmer,
                   progress_handler_t progress) {
                 for (auto it = unordered.begin(); it != unordered.end(); ++it) {
                   output_kmer(it->first, it->second);
                 }
               });

  kmer_set::const_iterator it1 = ks2.begin();
  std::set<kmer_t>::const_iterator it2 = m_cmp_set.begin();
  ASSERT_EQ(ks2.size(), m_cmp_set.size());
  while (it1 != ks2.end() && it2 != m_cmp_set.end()) {
    dna_sequence s1(*it1, m_kmer_size);
    dna_sequence s2(*it2, m_kmer_size);
    EXPECT_EQ(s1.as_string(), s2.as_string()) << count;
    EXPECT_EQ(it1.get_flags(), ordered_flags[count]) << count;
    kmer_set::const_iterator it1x = ks2.find(*it2);
    std::set<kmer_t>::const_iterator it2x = m_cmp_set.find(*it1);
    EXPECT_EQ(*it1, *it1x) << count;
    EXPECT_EQ(*it2, *it2x) << count;
    it1++;
    it2++;
    count++;
  }
  EXPECT_TRUE(it1 == ks2.end());
  EXPECT_TRUE(it2 == m_cmp_set.end());
}

INSTANTIATE_TEST_CASE_P(kmer_set_test, kmer_set_test,
                        // Kmer sizes
                        ::testing::Values(20, 21, 22, 23, 30, 31, 32));
