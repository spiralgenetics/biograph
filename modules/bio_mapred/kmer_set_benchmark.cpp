#include "benchmark/benchmark.h"
#include "modules/bio_mapred/kmer_set.h"
#include "modules/io/config.h"
#include "modules/io/log.h"
#include "modules/io/parallel.h"

namespace {

membuf raw_kmer_membuf;
const kmer_t* raw_kmers = nullptr;
size_t n_raw_kmers;
std::unique_ptr<kmer_set> loaded_ks;

void init_raw_kmers() {
  if (raw_kmer_membuf.size()) {
    return;
  }

  membuf mb = membuf(
      new mmap_buffer("/scratch/kmer_set_benchmark_kmers.dat", mmap_buffer::mode::read_only));

  raw_kmer_membuf = membuf(new owned_membuf(mb.data(), mb.size(), "raw_kmers"));
  raw_kmers = reinterpret_cast<const kmer_t*>(raw_kmer_membuf.data());
  n_raw_kmers = raw_kmer_membuf.size() / sizeof(kmer_t);
  CHECK_LT(n_raw_kmers, std::numeric_limits<uint32_t>::max());
  CHECK_GT(n_raw_kmers, 0);
}

void get_raw_kmers(const kmer_set::kmer_output_f& output_f, progress_handler_t progress) {
  CHECK(raw_kmer_membuf.size());

  parallel_for(
      0, n_raw_kmers,
      [&](size_t start, size_t limit) {
        for (size_t i = start; i != limit; ++i) {
          unsigned flags = 0;
          if (i & 1) {
            flags |= kmer_set::k_rev_starts_read;
          }
          if (i & 2) {
            flags |= kmer_set::k_fwd_starts_read;
          }
          output_f(raw_kmers[i], flags);
        }
      },
      progress);
}

std::unique_ptr<kmer_set> make_ks() {
  std::unique_ptr<kmer_set> ks =
      make_unique<kmer_set>(n_raw_kmers, 30 /* kmer size */, 100ULL*1024*1024, get_raw_kmers);
  CHECK_EQ(ks->size(), n_raw_kmers);
  return ks;
}

}  // namespace

static void BM_construct(benchmark::State& state) {
  init_raw_kmers();

  while (state.KeepRunning()) {
    std::unique_ptr<kmer_set> ks = make_ks();
    state.PauseTiming();
    loaded_ks = std::move(ks);
    state.ResumeTiming();
  }
}

BENCHMARK(BM_construct)->Unit(benchmark::kMillisecond);

static void BM_lookup(benchmark::State& state) {
  if (!loaded_ks) {
    init_raw_kmers();
    loaded_ks = make_ks();
  }

  const kmer_t* it = raw_kmers;
  const kmer_t* end = raw_kmers + n_raw_kmers;
  while (state.KeepRunning()) {
    benchmark::DoNotOptimize(loaded_ks->count(*it));
    ++it;
    if (it == end) {
      it = raw_kmers;
    }
  }
}

BENCHMARK(BM_lookup);

int main(int argc, char** argv) {
  log_init("kmer_set_benchmark", 2, 0);
  ::benchmark::Initialize(&argc, argv);
  if (::benchmark::ReportUnrecognizedArguments(argc, argv)) {
    return 1;
  }
  Config::load("etc/products/unittest.json");
  Config::set("resources_root", "/tmp/kmer_set_benchmark_storage");
  ::benchmark::RunSpecifiedBenchmarks();
  boost::filesystem::remove_all("/tmp/kmer_set_benchmark_storage");
  loaded_ks.reset();
  raw_kmer_membuf = membuf();
}
