#include "benchmark/benchmark.h"
#include "modules/variants/read_cov.h"

#include <random>

constexpr size_t k_data_size = 1024 * 1024 * 64 /* 64 MB */;
constexpr size_t k_num_data = k_data_size / sizeof(uint32_t);

std::random_device random_dev;

namespace variants {

class read_cov_benchmark {
 public:
  read_cov_benchmark(size_t table_size) : m_rand_source(random_dev()) {
    init_lookup(table_size);
    init_data(table_size);
  }

  void init_lookup(size_t table_size) {
    // Generate a complete cycle of lookups so we don't get stuck in little loops.
    std::vector<uint32_t> nums;
    for (size_t i = 0; i != table_size; ++i) {
      nums.push_back(i);
    }
    std::shuffle(nums.begin(), nums.end(), m_rand_source);
    m_lookup.resize(table_size);
    for (size_t i = 0; i != table_size; ++i) {
      if (i == 0) {
        m_lookup[i] = nums[table_size - 1];
      } else {
        m_lookup[i] = nums[i - 1];
      }
    }
  }

  void init_data(size_t table_size) {
    std::uniform_int_distribution<uint32_t> rand_val(0, table_size - 1);

    m_data.reserve(k_num_data);

    for (size_t i = 0; i != k_num_data; ++i) {
      m_data.push_back(rand_val(m_rand_source));
    }
  }

  void do_translate() {
    read_cov::translate_uint32s(m_data.data(), m_data.data() + m_data.size(), m_lookup);
  }

 private:
  std::mt19937 m_rand_source;

  std::vector<uint32_t> m_lookup;
  std::vector<uint32_t> m_data;
};

}  // namespace variants

using namespace variants;

static void BM_translate(benchmark::State& state) {
  read_cov_benchmark b(state.range(0));

  while (state.KeepRunning()) {
    b.do_translate();
  }
  state.SetBytesProcessed(int64_t(state.iterations()) * k_num_data * sizeof(uint32_t));
}

BENCHMARK(BM_translate)->RangeMultiplier(2)->Range(1, 1024 * 1024);

BENCHMARK_MAIN();
