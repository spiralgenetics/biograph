#include "modules/io/bitcount.h"

#include "benchmark/benchmark.h"

#include <iomanip>
#include <iostream>
#include <random>

namespace {

const size_t nbits = 16ULL * 1024 * 1024 * 1024;
size_t bufsize = 0;
char* bcbuf = nullptr;

std::unique_ptr<bitcount> bc;

size_t seed = time(0);
std::mt19937_64 random_source(seed);

void init_bc() {
  if (bcbuf) {
    return;
  }

  bufsize = bitcount::compute_size(nbits);
  bcbuf = new char[bufsize];

  uint64_t* bcbuf_64 = reinterpret_cast<uint64_t*>(bcbuf);

  std::cerr << "Populating bitcount with random data, seed " << seed << "\n";
  CHECK_EQ(0, random_source.min());
  CHECK_EQ(std::numeric_limits<uint64_t>::max(), random_source.max());

  size_t tot_uint64s = bufsize / sizeof(uint64_t);
  for (size_t i = 0; i < tot_uint64s; ++i) {
    if (!(i & ((1 << 20) - 1))) {
      std::cerr << " " << i;
    }
    bcbuf_64[i] = random_source();

    // Have some sections be sparser than others.
    if (i < (tot_uint64s/2)) {
      // Densify first half.
      size_t shifted = i >> 10;
      while (shifted) {
        bcbuf_64[i] |= random_source();
        shifted >>= 10;
      }
    } else {
      // Sparsify second half.
      size_t shifted = (tot_uint64s - i) >> 10;
      while (shifted) {
        bcbuf_64[i] &= random_source();
        shifted >>= 10;
      }
    }
  }

  bc.reset(new bitcount(bcbuf, nbits));
  std::cerr << "\nFinalizing...\n";
  bc->finalize();
  std::cerr << bc->total_bits() << " total bits present out of " << bc->size()
            << ": " << std::fixed << std::setw(6) << std::setprecision(2)
            << bc->total_bits() * 100. / bc->size() << "%.\n";
  std::cerr << "Done\n";
}

}  // namespace

static void BM_bitcount_count(benchmark::State& state) {
  init_bc();

  bitcount bc(bcbuf, nbits);

  std::uniform_int_distribution<size_t> random_pos(0, nbits - 1);
  size_t pos = random_pos(random_source);
  const size_t stride = random_pos(random_source);

  while (state.KeepRunning()) {
    size_t count = bc.count(pos);
    pos += stride ^ count;
    pos %= nbits;
  }
}

BENCHMARK(BM_bitcount_count);

static void BM_bitcount_find_count(benchmark::State& state) {
  init_bc();

  size_t total_bits_set = bc->total_bits();

  std::uniform_int_distribution<size_t> random_pos(0, total_bits_set - 1);
  size_t pos = random_pos(random_source);
  const size_t stride = random_pos(random_source);

  while (state.KeepRunning()) {
    size_t count = bc->find_count(pos);
    pos += stride ^ count;
    pos %= total_bits_set;
  }
}

BENCHMARK(BM_bitcount_find_count);

static void BM_bitcount_find_count_with_index(benchmark::State& state) {
  init_bc();

  bitcount bc(bcbuf, nbits);
  bc.make_find_count_index();

  size_t total_bits_set = bc.total_bits();

  std::uniform_int_distribution<size_t> random_pos(0, total_bits_set - 1);
  size_t pos = random_pos(random_source);
  const size_t stride = random_pos(random_source);

  while (state.KeepRunning()) {
    size_t count = bc.find_count(pos);
    pos += stride ^ count;
    pos %= total_bits_set;
  }
}

BENCHMARK(BM_bitcount_find_count_with_index);

static void BM_bitcount_make_find_count_index(benchmark::State& state) {
  init_bc();

  bitcount bc(bcbuf, nbits);
  while (state.KeepRunning()) {
    bc.make_find_count_index();
  }
}
BENCHMARK(BM_bitcount_make_find_count_index);

static void BM_bitcount_finalize(benchmark::State& state) {
  init_bc();

  bitcount bc(bcbuf, nbits);
  while (state.KeepRunning()) {
    bc.finalize();
  }
}
BENCHMARK(BM_bitcount_finalize);

BENCHMARK_MAIN();
