#include "modules/build_seqset/repo_seq.h"
#include "benchmark/benchmark.h"
#include "modules/bio_base/dna_sequence.h"
#include "modules/bio_base/reference.h"
#include "modules/io/config.h"
#include "modules/io/progress.h"

#include <random>

constexpr size_t k_min_read_size = 8;
constexpr size_t k_max_read_size = 250;
constexpr size_t k_max_read_count = 300 * 1000 * 1000;
constexpr size_t k_suffix_count = 4;

using namespace build_seqset;

static void update_progress(double new_progress) {
  static float prev_progress = 0;
  if (fabs(new_progress - prev_progress) > 0.0001) {
    prev_progress = new_progress;
    print_progress(new_progress);
  }
}

class rand_repo {
 public:
  rand_repo() : m_rand_source(std::random_device()()) {
    static size_t n = 0;
    m_ref_path = "/scratch/ref";
    m_ref_path += std::to_string(n);
    unlink(m_ref_path.c_str());

    m_repo_path = "/scratch/repo";
    m_repo_path += std::to_string(n);
    unlink(m_repo_path.c_str());

    m_ref_builder.emplace(m_ref_path);
    m_repo_builder.emplace(m_repo_path);

    m_ref.emplace("", "/reference/hs37d5");
  }

  ~rand_repo() {
    unlink(m_ref_path.c_str());
    unlink(m_repo_path.c_str());
  }

  void add_rand_reads() {
    std::mutex mu;
    size_t ref_size = m_ref->size();
    std::cerr << "Adding ref repo...\n";
    size_t ref_start = m_repo_builder->write_seq(dna_slice(m_ref->get_dna(0), ref_size));
    std::cerr << "Done adding ref repo.  Adding reads...\n";
    CHECK_GT(ref_size, k_max_read_size);
    std::uniform_int_distribution<size_t> rand_ref_pos(ref_start,
                                                       ref_start + ref_size - k_max_read_size);
    std::uniform_int_distribution<size_t> rand_read_size_and_rc(k_min_read_size << 1,
                                                                k_max_read_size << 1);
    parallel_for(  //
        0, k_max_read_count,
        [&](size_t start, size_t limit) {
          std::unique_lock<std::mutex> l(mu);
          std::mt19937_64 rand_source(m_rand_source());
          l.unlock();

          std::vector<seq_repository::entry_data> chunk_reads;
          chunk_reads.reserve(limit - start);
          for (size_t i = start; i != limit; ++i) {
            size_t read_pos = rand_ref_pos(rand_source);
            size_t read_size_and_rc = rand_read_size_and_rc(rand_source);
            bool is_rc = read_size_and_rc & 1;
            size_t read_size = read_size_and_rc >> 1;
            size_t inline_size = std::min<size_t>(seq_repository::k_inline_bases, read_size);
            dna_slice seq(m_ref->get_dna(read_pos), read_size);

            if (is_rc) {
              seq = seq.rev_comp();
              read_pos += read_size;
            }

            dna_slice inline_part = seq.subseq(0, inline_size);

            if (read_size <= seq_repository::k_inline_bases) {
              read_pos = seq_repository::k_max_offset;
            }

            chunk_reads.emplace_back(read_size, inline_part, read_pos, is_rc);
          }

          l.lock();
          m_ref_builder->write_entries_and_clear(chunk_reads, true /* block on lock */);
          CHECK(chunk_reads.empty());
          m_read_count += limit - start;
        },
        update_progress);
    std::cerr << "\nDone adding reads\n";
  }

  void finalize() {
    m_ref_builder.reset();
    m_repo_builder.reset();
    m_entries.emplace(m_ref_path, m_repo_path);
    CHECK(m_entries->data_begin() + m_read_count == m_entries->data_end());
    size_t repo_size = boost::filesystem::file_size(m_repo_path);
    std::cerr << m_read_count << " reads added with a repo size totalling "
              << (repo_size / 1024 / 1024) << " MB\n";
  }

  void init_pass(size_t read_count) {
    m_data.reset();
    CHECK_LE(read_count, m_read_count);
    m_data.emplace(read_count, track_alloc("repo_seq_benchmark:entry_data"));
    std::copy(m_entries->data_begin(), m_entries->data_begin() + read_count, m_data->begin());
  }

  void do_sort() { m_entries->sort_entry_data(m_entries->data_begin(), m_entries->data_end()); }

 private:
  std::string m_ref_path;
  std::string m_repo_path;

  size_t m_read_count = 0;

  boost::optional<seq_repository::ref_builder> m_ref_builder;
  boost::optional<seq_repository::repo_builder> m_repo_builder;
  boost::optional<reference> m_ref;
  boost::optional<seq_repository> m_entries;
  std::mt19937_64 m_rand_source;
  boost::optional<tracked_vector<seq_repository::entry_data>> m_data;
};

boost::optional<rand_repo> g_rand_repo;

static void BM_sort(benchmark::State& state) {
  if (!g_rand_repo) {
    g_rand_repo.emplace();
    g_rand_repo->add_rand_reads();
    g_rand_repo->finalize();
  }
  CHECK_LE(state.range(0), k_max_read_count);
  size_t read_count = state.range(0);
  size_t bytes_processed = 0;
  while (state.KeepRunning()) {
    state.PauseTiming();
    g_rand_repo->init_pass(state.range(0));
    bytes_processed += read_count * sizeof(seq_repository::entry_data);
    state.ResumeTiming();
    g_rand_repo->do_sort();
  }
  state.SetBytesProcessed(bytes_processed);
}

BENCHMARK(BM_sort)
    ->Unit(benchmark::kMillisecond)
    ->RangeMultiplier(2)
    ->Range(1000, k_max_read_count);

int main(int argc, char** argv) {
  ::benchmark::Initialize(&argc, argv);
  if (::benchmark::ReportUnrecognizedArguments(argc, argv)) {
    return 1;
  }
  Config::load("etc/products/unittest.json");
  ::benchmark::RunSpecifiedBenchmarks();
  g_rand_repo.reset();
}
