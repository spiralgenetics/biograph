#include "modules/io/packed_varbit_vector.h"
#include "modules/io/parallel.h"

#include "benchmark/benchmark.h"

#include <iostream>
#include <random>

namespace {

const size_t k_table_bytes = 4ULL * 1024 * 1024;
size_t seed = time(0);
std::mt19937_64 random_source(seed);

}  // namespace

static void AllBits(benchmark::internal::Benchmark* b) {
  for (int i = 0; i <= 64; ++i) {
    b->Arg(i);
  }
}

static void BM_get(benchmark::State& state) {
  unsigned bits_per_value = state.range(0);

  size_t element_count = k_table_bytes * 8 / std::max<size_t>(bits_per_value, 1);
  size_t max_value = (~uint64_t(0)) >> (64 - bits_per_value);

  mutable_packed_varbit_vector v(element_count, max_value, "packed_varbit_vector_benchmark");
  mutable_membuf mb = v.get_internal_elements();
  uint64_t* ptr = reinterpret_cast<uint64_t*>(mb.mutable_data());
  // Make sure all the pages get populated
  for (uint64_t i = 0; i < mb.size() / sizeof(uint64_t); i += 1000) {
    ptr[i] = random_source();
  }

  std::uniform_int_distribution<size_t> random_pos(0, element_count - 1);
  std::unique_ptr<int_map_interface> intf = v.get_int_map_interface();
  while (state.KeepRunning()) {
    benchmark::DoNotOptimize(intf->get(random_pos(random_source)));
  }
}

BENCHMARK(BM_get)->Apply(AllBits);

static void BM_set(benchmark::State& state) {
  unsigned bits_per_value = state.range(0);

  size_t element_count = k_table_bytes * 8 / std::max<size_t>(bits_per_value, 1);
  size_t max_value = (~uint64_t(0)) >> (64 - bits_per_value);

  mutable_packed_varbit_vector v(element_count, max_value, "packed_vabrit_vector_benchmark");

  std::uniform_int_distribution<size_t> random_pos(0, element_count - 1);
  std::uniform_int_distribution<size_t> random_val(0, max_value - 1);
  size_t val = 0;
  while (state.KeepRunning()) {
    v.set(random_pos(random_source), val);
    ++val;
    if (val >= max_value) {
      val = 0;
    }
  }
}

BENCHMARK(BM_set)->Apply(AllBits);

BENCHMARK_MAIN();
