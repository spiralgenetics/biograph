#include "benchmark/benchmark.h"
#include "modules/bio_base/dna_sequence.h"
#include "modules/bio_base/dna_testutil.h"

#include <random>

namespace {

constexpr unsigned k_num_random_seqs = 10000;
constexpr unsigned k_random_seq_len = 100;
constexpr unsigned k_shared_prefix_len = 60;

std::vector<dna_sequence> g_random_seqs;
std::vector<dna_slice> g_random_slices;
std::vector<dna_slice> g_random_rc_slices;
std::vector<dna_slice> g_random_mixed_slices;

std::mt19937 g_rand_source;

void make_random_seqs() {
  if (g_random_seqs.empty()) {
    dna_sequence offset_seq("AAA");
    dna_sequence shared_prefix_seq = rand_dna_sequence(g_rand_source, k_shared_prefix_len);

    for (unsigned i = 0; i != k_num_random_seqs; ++i) {
      std::uniform_int_distribution<unsigned> offset_gen(0, 3);
      unsigned offset = offset_gen(g_rand_source);

      g_random_seqs.push_back(offset_seq.subseq(0, offset) + shared_prefix_seq +
                              rand_dna_sequence(g_rand_source, k_random_seq_len) +
                              shared_prefix_seq.rev_comp() + offset_seq.subseq(0, offset));

      const auto& seq = g_random_seqs.back();

      g_random_slices.emplace_back(seq.begin() + offset, seq.size() - offset);
      g_random_rc_slices.emplace_back(seq.rcbegin() + offset, seq.size() - offset);
      if (g_random_mixed_slices.size() & 1) {
        g_random_mixed_slices.emplace_back(seq.rcbegin() + offset, seq.size() - offset);

      } else {
        g_random_mixed_slices.emplace_back(seq.begin() + offset, seq.size() - offset);
      }
    }
    std::sort(g_random_slices.begin(), g_random_slices.end());
    std::sort(g_random_rc_slices.begin(), g_random_rc_slices.end());
    std::sort(g_random_mixed_slices.begin(), g_random_mixed_slices.end());
  }
}

}  // namespace

static void BM_compare_dna_slice(benchmark::State& state) {
  make_random_seqs();

  auto it = g_random_slices.begin();
  auto next = it + 1;

  while (state.KeepRunning()) {
    benchmark::DoNotOptimize(*it < *next);
    it = next;
    ++next;
    if (next == g_random_slices.end()) {
      next = g_random_slices.begin();
    }
  }
}

BENCHMARK(BM_compare_dna_slice);

static void BM_compare_dna_rc_slice(benchmark::State& state) {
  make_random_seqs();

  auto it = g_random_rc_slices.begin();
  auto next = it + 1;

  while (state.KeepRunning()) {
    benchmark::DoNotOptimize(*it < *next);
    it = next;
    ++next;
    if (next == g_random_rc_slices.end()) {
      next = g_random_rc_slices.begin();
    }
  }
}

BENCHMARK(BM_compare_dna_rc_slice);

static void BM_compare_dna_mixed_slice(benchmark::State& state) {
  make_random_seqs();

  auto it = g_random_mixed_slices.begin();
  auto next = it + 1;

  while (state.KeepRunning()) {
    benchmark::DoNotOptimize(*it < *next);
    it = next;
    ++next;
    if (next == g_random_mixed_slices.end()) {
      next = g_random_mixed_slices.begin();
    }
  }
}

BENCHMARK(BM_compare_dna_mixed_slice);

/*
  Run on (16 X 2199.99 MHz CPU s)
  2019-09-01 19:08:48
  -----------------------------------------------------------------------
  Benchmark                                Time           CPU Iterations
  -----------------------------------------------------------------------
  BM_copy_dna_sequence/1                  58 ns         58 ns   12119679
  BM_copy_dna_sequence/10093             195 ns        195 ns    3579309
  BM_assign_dna_sequence/1                64 ns         64 ns   10990189
  BM_assign_dna_sequence/10093           202 ns        202 ns    3474915
  BM_copy_dna_slice/1                     71 ns         71 ns    9885445
  BM_copy_dna_slice/10093              84392 ns      84390 ns       8354
  BM_assign_dna_slice/1                   69 ns         69 ns   10095238
  BM_assign_dna_slice/10093            80505 ns      80506 ns       8884
  BM_copy_dna_rc_slice/1                  96 ns         96 ns    6913265
  BM_copy_dna_rc_slice/10093           80205 ns      80205 ns       8556
  BM_assign_dna_rc_slice/1                69 ns         69 ns   10288428
  BM_assign_dna_rc_slice/10093         80466 ns      80466 ns       8548
  BM_copy_dna_subslice/1                  10 ns         10 ns   72920651
  BM_copy_dna_subslice/10093           86336 ns      86336 ns       8166
*/

namespace {

template <bool do_slice, bool do_rc, bool do_subslice, bool do_assign>
void run_copy_benchmark(benchmark::State& state) {
  dna_sequence seq = rand_dna_sequence(g_rand_source, state.range(0));
  dna_slice slice = seq;
  if (do_rc) {
    slice = slice.rev_comp();
  }
  if (do_subslice) {
    slice = slice.subseq(1, slice.size() - 1);
  }
  while (state.KeepRunning()) {
    if (do_slice) {
      if (do_assign) {
        dna_sequence assign_result;
        assign_result = slice;
        benchmark::DoNotOptimize(&assign_result);
      } else {
        benchmark::DoNotOptimize(dna_sequence(slice));
      }
    } else {
      if (do_assign) {
        dna_sequence assign_result;
        assign_result = seq;
        benchmark::DoNotOptimize(&assign_result);
      } else {
        benchmark::DoNotOptimize(dna_sequence(seq));
      }
    }
  }
}

}  // namespace

#define DNA_COPY_BENCHMARK(NAME, SLICE, RC, SUBSLICE, ASSIGN) \
  static void NAME(benchmark::State& state) {                 \
    run_copy_benchmark<SLICE, RC, SUBSLICE, ASSIGN>(state);   \
  }                                                           \
                                                              \
  BENCHMARK(NAME)->Arg(1)->Arg(32768)

DNA_COPY_BENCHMARK(BM_copy_dna_sequence, false, false, false, false);
DNA_COPY_BENCHMARK(BM_assign_dna_sequence, false, false, false, true);

DNA_COPY_BENCHMARK(BM_copy_dna_slice, true, false, false, false);
DNA_COPY_BENCHMARK(BM_assign_dna_slice, true, false, false, true);

DNA_COPY_BENCHMARK(BM_copy_dna_rc_slice, true, true, false, false);
DNA_COPY_BENCHMARK(BM_assign_dna_rc_slice, true, true, false, true);

DNA_COPY_BENCHMARK(BM_copy_dna_subslice, true, false, true, false);
DNA_COPY_BENCHMARK(BM_assign_dna_subslice, true, false, true, true);

DNA_COPY_BENCHMARK(BM_copy_dna_rc_subslice, true, true, true, false);
DNA_COPY_BENCHMARK(BM_assign_dna_rc_subslice, true, true, true, true);

#undef DNA_COPY_BENCHMARK

BENCHMARK_MAIN();
