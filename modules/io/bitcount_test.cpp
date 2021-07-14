#include "modules/io/bitcount.h"
#include "modules/io/log.h"
#include "modules/io/spiral_file_mem.h"
#include "modules/test/test_coverage.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <memory>
#include <random>
#include <vector>

using namespace testing;

TEST(bitcount, size) {
  uint64_t x = 0xfffe0000ffff0000;
  SPLOG("Sizeof unsigned int = %d", (int)sizeof(unsigned int));
  // __builtin_popcount should only count the lower 32 bits.
  EXPECT_EQ(16, __builtin_popcount(x));
  // __builtin_popcount should count all the bits.
  EXPECT_EQ(31, __builtin_popcountl(x));
}

class fake_bitcount {
 public:
  fake_bitcount(size_t size) : m_vec(size), m_tots(size + 1) {}
  void set(size_t i, bool v) { m_vec[i] = v; }
  void finalize() {
    size_t tot = 0;
    for (size_t i = 0; i < m_vec.size(); i++) {
      m_tots[i] = tot;
      tot += m_vec[i];
    }
    m_tots[m_vec.size()] = tot;
  }
  bool get(size_t i) { return m_vec[i]; }
  size_t count(size_t i) { return m_tots[i]; }

 private:
  std::vector<char> m_vec;
  std::vector<size_t> m_tots;
};

enum bitcount_constructor_type { OLD_STYLE_BUFFER, SPIRAL_FILE };

class bitcount_test : public testing::TestWithParam<bitcount_constructor_type> {
 public:
  // Since there is a lot of weird packing logic, do a lot of special size tests
  void test_size(size_t bc_size) {
    // SPLOG("Size test: %d", (int) bc_size);
    create_bc(bc_size);
    for (size_t i = 0; i < bc_size; i++) m_bc->set(i, 1);
    finalize_bc();
    for (size_t i = 0; i <= bc_size; i++) {
      // SPLOG("Input %d", (int) i);
      size_t out = m_bc_ro->count(i);
      EXPECT_EQ(out, i);
    }
    EXPECT_EQ(bc_size, m_bc_ro->size());
    EXPECT_EQ(bc_size, m_bc_ro->total_bits());
    EXPECT_EQ(bc_size, m_bc_ro->find_count(m_bc_ro->total_bits()));
    close_bc();
  }

  void create_bc(size_t nbits) {
    CHECK(!m_bc);
    CHECK(!m_bc_ro);
    m_nbits = nbits;
    switch (GetParam()) {
      case OLD_STYLE_BUFFER: {
        size_t bc_mem = bitcount::compute_size(nbits);
        m_buf.reset(new char[bc_mem]);
        m_bc.reset(new bitcount(m_buf.get(), nbits));
        m_bc->init();
        break;
      }
      case SPIRAL_FILE: {
        m_creator.reset(new spiral_file_create_mem());
        m_bc.reset(new bitcount(m_creator->create(), nbits));
        break;
      }
    }
  }

  void finalize_bc() {
    CHECK(m_bc);
    CHECK(!m_bc_ro);
    switch (GetParam()) {
      case OLD_STYLE_BUFFER: {
        m_bc->finalize();
        m_bc.reset();
        m_bc_ro.reset(new bitcount((const char *)m_buf.get(), m_nbits));
        break;
      }
      case SPIRAL_FILE: {
        CHECK(m_creator);
        CHECK(!m_opener);
        m_bc->finalize();
        m_bc.reset();
        spiral_file_mem_storage encoded = m_creator->close();
        m_creator.reset();

        m_opener.reset(new spiral_file_open_mem(encoded));
        m_bc_ro.reset(new bitcount(m_opener->open()));
      }
    }
  }

  void close_bc() {
    CHECK(!m_bc);
    CHECK(m_bc_ro);
    switch (GetParam()) {
      case OLD_STYLE_BUFFER: {
        m_bc_ro.reset();
        m_buf.reset();
        break;
      }
      case SPIRAL_FILE: {
        CHECK(!m_creator);
        CHECK(m_opener);
        m_bc_ro.reset();
        m_opener.reset();
        break;
      }
    }
  }

  size_t m_nbits = 0;

  std::unique_ptr<char[]> m_buf;
  std::unique_ptr<bitcount> m_bc;
  std::unique_ptr<bitcount> m_bc_ro;

  std::unique_ptr<spiral_file_create_mem> m_creator;
  std::unique_ptr<spiral_file_open_mem> m_opener;
};

TEST_P(bitcount_test, test) {
  size_t bc_size = 1000001;
  create_bc(bc_size);
  fake_bitcount bc2(bc_size);

  for (size_t i = 0; i < bc_size; i++) {
    bool x = random() % 2;
    m_bc->set(i, x);
    bc2.set(i, x);
  }
  finalize_bc();
  bc2.finalize();
  for (size_t i = 0; i < bc_size; i++) {
    size_t t1 = m_bc_ro->count(i);
    size_t t2 = bc2.count(i);
    // SPLOG("i = %d, t1 = %d, t2 = %d", (int) i, (int) t1, (int) t2);
    EXPECT_EQ(t1, t2);
    EXPECT_EQ(m_bc_ro->get(i), bc2.get(i));
  }

  size_t tot1 = m_bc_ro->count(bc_size);
  size_t tot2 = bc2.count(bc_size);
  // SPLOG("i = %d, t1 = %d, t2 = %d", (int) i, (int) t1, (int) t2);
  EXPECT_EQ(tot1, tot2);
}

TEST_P(bitcount_test, sizes) {
  // Yes, we even test size 0
  for (size_t i = 0; i <= 1024; i++) {
    test_size(i);
  }
}

TEST_P(bitcount_test, find_count) {
  size_t bitcount_size = 1024;
  create_bc(bitcount_size);
  for (size_t i = 0; i < bitcount_size; i++) {
    m_bc->set(i, i % 2);
  }
  finalize_bc();

  for (size_t i = 0; i < bitcount_size; i++) {
    EXPECT_EQ(m_bc_ro->count(i), i / 2);
    if (i < bitcount_size / 2) {
      EXPECT_EQ(m_bc_ro->find_count(i), 2 * i + 1) << i;
    }
  }

  EXPECT_EQ(m_bc_ro->size(), bitcount_size);
  EXPECT_EQ(m_bc_ro->total_bits(), bitcount_size / 2);
}

TEST_P(bitcount_test, find_count_with_index) {
  size_t bitcount_size = 1024;
  create_bc(bitcount_size);
  for (size_t i = 0; i < bitcount_size; i++) {
    m_bc->set(i, i % 2);
  }
  finalize_bc();
  m_bc_ro->make_find_count_index();

  for (size_t i = 0; i < bitcount_size; i++) {
    EXPECT_EQ(m_bc_ro->count(i), i / 2);
    if (i < bitcount_size / 2) {
      EXPECT_EQ(m_bc_ro->find_count(i), 2 * i + 1) << i;
    }
  }

  EXPECT_EQ(m_bc_ro->size(), bitcount_size);
  EXPECT_EQ(m_bc_ro->total_bits(), bitcount_size / 2);
}

INSTANTIATE_TEST_CASE_P(old_style_buffer_tests, bitcount_test, ::testing::Values(OLD_STYLE_BUFFER));
INSTANTIATE_TEST_CASE_P(spiral_file_tests, bitcount_test, ::testing::Values(SPIRAL_FILE));

TEST(bitcount_coverage_test, coverage) {
  scoped_test_coverage cov;

  size_t seed = size_t(time(0)) ^ (size_t(getpid()) << 32);
  std::cerr << "Generating random bitcount with seed " << seed << "\n";
  std::mt19937_64 random_source(seed);
  std::uniform_int_distribution<size_t> size_picker(1, 1 << 18);
  std::uniform_int_distribution<int> bool_picker(0, 1);

  std::string last_missing_coverage;

  while (!cov.missing("bitcount").empty()) {
    std::string missing_coverage = ::testing::PrintToString(cov.missing("bitcount"));
    if (missing_coverage != last_missing_coverage) {
      std::cerr << "Missing coverage: " << missing_coverage << "\n";
      last_missing_coverage = missing_coverage;
    }
    bool invert = bool_picker(random_source);
    size_t size = size_picker(random_source);
    std::uniform_int_distribution<size_t> stride_picker(0, 256);

    bitcount bc(size);

    size_t stride_left = stride_picker(random_source) * stride_picker(random_source);
    size_t pos = 0;
    size_t total_set_bits = 0;
    size_t last_set_bit = 0;
    while (pos < size) {
      bool val = invert;
      if (stride_left-- == 0) {
        val = !val;
        stride_left = stride_picker(random_source) * stride_picker(random_source);
      }
      bc.set(pos, val);
      if (val) {
        total_set_bits++;
        last_set_bit = pos;
      }
      pos++;
    }
    bc.finalize();
    if (bool_picker(random_source)) {
      bc.make_find_count_index();
    }

    for (pos = 0; pos < size; ++pos) {
      if (bc.get(pos)) {
        size_t c = bc.count(pos);
        EXPECT_EQ(pos, bc.find_count(c));
      }
    }

    EXPECT_EQ(total_set_bits, bc.total_bits());
    if (total_set_bits > 0) {
      EXPECT_EQ(total_set_bits - 1, bc.count(last_set_bit));
    } else {
      EXPECT_EQ(0, bc.count(size - 1));
    }
    EXPECT_EQ(bc.size(), bc.find_count(bc.total_bits()));
  }
}
