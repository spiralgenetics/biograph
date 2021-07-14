#include "benchmark/benchmark.h"
#include "modules/bio_base/shannon_entropy.h"
#include "modules/bio_base/dna_sequence.h"
#include "modules/bio_base/dna_testutil.h"

namespace {

constexpr unsigned k_seq_len = 1024*1024;

}  // namespace

static void BM_shannon_entropy(benchmark::State& state) {
  std::random_device true_rand;
  std::mt19937 rand_source(true_rand());
  
  dna_sequence seq = rand_dna_sequence(rand_source, k_seq_len);
  dna_const_iterator begin = seq.begin();
  dna_const_iterator it = seq.begin();
  dna_const_iterator end = seq.end();

  shannon_entropy e(state.range(0));
  
  while (state.KeepRunning()) {
    e.push_front(*it);
    ++it;
    if (it == end) {
      it = begin;
    }
    benchmark::DoNotOptimize(e.length_needed());
  }
}

BENCHMARK(BM_shannon_entropy)
    ->Arg(20)
    ->Arg(50)
    ->Arg(80)
    ->Arg(110)
    ->Arg(140)
    ->Arg(170)
    ->Arg(200)
    ->Arg(230);

BENCHMARK_MAIN();
