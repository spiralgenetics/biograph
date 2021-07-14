#include "modules/bio_base/seqset.h"
#include "benchmark/benchmark.h"
#include "modules/io/utils.h"

#include <math.h>
#include <memory>
#include <random>
#include <vector>

namespace {

void update_progress(const float& new_progress) {
  static float prev_progress = 0;
  if (fabs(new_progress - prev_progress) > 0.0001) {
    prev_progress = new_progress;
    print_progress(new_progress);
  }
}

std::unique_ptr<seqset_file> g_seqset_file;
const seqset* g_seqset = nullptr;
std::mt19937 g_rand_source;
boost::optional<std::uniform_int_distribution<uint64_t>> g_get_seqset_id;

seqset_range rand_seqset_entry() { return g_seqset->ctx_entry((*g_get_seqset_id)(g_rand_source)); }

void init_seqset(bool pop_front_cache) {
  if (!g_seqset_file) {
    g_seqset_file.reset(new seqset_file("/scratch/HG001.hs37d5.50x.11197.bg/seqset"));
    g_seqset = &g_seqset_file->get_seqset();
    g_seqset->membufs().cache_in_memory(update_progress);
    g_get_seqset_id.emplace(0, g_seqset->size() - 1);
  }

  if (pop_front_cache && !g_seqset->is_pop_front_cached()) {
    g_seqset->populate_pop_front_cache(update_progress);
  } else if (!pop_front_cache && g_seqset->is_pop_front_cached()) {
    g_seqset->clear_pop_front_cache();
  }
}

constexpr size_t k_random_seq_chunk_size = 10000;

void fill_random_seqs(std::vector<dna_sequence>& seqs) {
  seqs.clear();
  seqs.reserve(k_random_seq_chunk_size);
  for (size_t i = 0; i != k_random_seq_chunk_size; ++i) {
    seqset_range e = rand_seqset_entry();
    seqs.emplace_back(e.sequence());
  }
}

}  // namespace

static void BM_seqset_sequence(benchmark::State& state) {
  init_seqset(false /* no pop front cache */);

  while (state.KeepRunning()) {
    benchmark::DoNotOptimize(rand_seqset_entry().sequence());
  }
}

BENCHMARK(BM_seqset_sequence);

static void BM_seqset_sequence_with_cache(benchmark::State& state) {
  init_seqset(true /* pop front cache */);

  while (state.KeepRunning()) {
    benchmark::DoNotOptimize(rand_seqset_entry().sequence());
  }
  g_seqset->clear_pop_front_cache();
}

BENCHMARK(BM_seqset_sequence_with_cache);

static void BM_seqset_is_maximal(benchmark::State& state) {
  init_seqset(false /* no pop front cache */);

  while (state.KeepRunning()) {
    benchmark::DoNotOptimize(rand_seqset_entry().is_maximal());
  }
}

BENCHMARK(BM_seqset_is_maximal);

static void BM_seqset_find_overlap_reads(benchmark::State& state) {
  init_seqset(false /* no pop front cache */);

  while (state.KeepRunning()) {
    overlaps_t overlaps;
    rand_seqset_entry().find_overlap_reads(overlaps, state.range(0), 60);
  }
}

BENCHMARK(BM_seqset_find_overlap_reads)->RangeMultiplier(2)->Range(1, 256);

static void BM_seqset_find_overlap_reads_fair(benchmark::State& state) {
  init_seqset(false /* no pop front cache */);

  while (state.KeepRunning()) {
    benchmark::DoNotOptimize(rand_seqset_entry().find_overlap_reads_fair(state.range(0), 60));
  }
}

BENCHMARK(BM_seqset_find_overlap_reads_fair)->RangeMultiplier(2)->Range(1, 256);

static void BM_seqset_find(benchmark::State& state) {
  init_seqset(false);

  std::vector<dna_sequence> seqs_to_find;
  auto seqs_to_find_it = seqs_to_find.end();

  while (state.KeepRunning()) {
    if (seqs_to_find_it == seqs_to_find.end()) {
      state.PauseTiming();
      fill_random_seqs(seqs_to_find);
      seqs_to_find_it = seqs_to_find.begin();
      state.ResumeTiming();
    }
    seqset_range r = g_seqset->find(*seqs_to_find_it);
    DCHECK_EQ(r.begin(), g_seqset->find_existing_unique(*seqs_to_find_it, 1));
    benchmark::DoNotOptimize(r);
    ++seqs_to_find_it;
  }
}

BENCHMARK(BM_seqset_find);

static void BM_seqset_find_unique(benchmark::State& state) {
  size_t expected_unique_len = state.range(0);

  init_seqset(false);

  std::vector<dna_sequence> seqs_to_find;
  auto seqs_to_find_it = seqs_to_find.end();

  while (state.KeepRunning()) {
    if (seqs_to_find_it == seqs_to_find.end()) {
      state.PauseTiming();
      fill_random_seqs(seqs_to_find);
      seqs_to_find_it = seqs_to_find.begin();
      state.ResumeTiming();
    }
    benchmark::DoNotOptimize(g_seqset->find_existing_unique(*seqs_to_find_it, expected_unique_len));
    ++seqs_to_find_it;
  }
}

BENCHMARK(BM_seqset_find_unique)
    ->Arg(1)
    ->Arg(5)
    ->Arg(10)
    ->Arg(20)
    ->Arg(23)
    ->Arg(26)
    ->Arg(29)
    ->Arg(32)
    ->Arg(35)
    ->Arg(38)
    ->Arg(41)
    ->Arg(50)
    ->Arg(100);

int main(int argc, char** argv) {
  ::benchmark::Initialize(&argc, argv);
  if (::benchmark::ReportUnrecognizedArguments(argc, argv)) {
    return 1;
  }
  ::benchmark::RunSpecifiedBenchmarks();
  g_seqset_file.reset();
}
