#include "modules/bio_mapred/kmerize_bf.h"
#include "datavis/kmer_quality_report/kmer_quality_report.h"
#include "modules/bio_base/overrep.h"
#include "modules/bio_base/reference.h"
#include "modules/bio_base/unaligned_read.h"
#include "modules/bio_format/exporter.h"
#include "modules/bio_mapred/kmer_set.h"
#include "modules/build_seqset/kmer_counter.h"
#include "modules/io/config.h"
#include "modules/io/make_unique.h"
#include "modules/io/parallel.h"
#include "modules/io/progress.h"
#include "modules/io/semaphore.h"
#include "modules/io/stats.h"
#include "modules/io/stopwatch.h"
#include "modules/mapred/output_stream.h"
#include "modules/mapred/sort_task.h"

#include <boost/algorithm/string/predicate.hpp>

#include <fstream>
#include <future>
#include <list>
#include <thread>

REGISTER_TASK(kmerize_bf_task);
REGISTER_TASK(kmerize_bf_subtask);

void kmerize_bf_params::validate() {
  SPLOG_P(LOG_DEBUG,
          "kmerize_bf_params::validate> kmer_size: %lu ref: '%s' ref_size: %lu num_threads: %lu",
          kmer_size, reference.c_str(), ref_size, num_threads);

  if ((kmer_size < 16) or (kmer_size > 32)) {
    SPLOG_P(LOG_DEBUG, "Invalid kmer_size");
    throw io_exception("Invalid kmer_size");
  }

  if (read_length > 0 && read_length < kmer_size) {
    SPLOG_P(LOG_DEBUG, "read_length must be greater than kmer_size");
    throw io_exception("read_length must be greater than kmer_size");
  }

  if (reference.empty() && ref_size == 0) {
    SPLOG_P(LOG_DEBUG, "reference or ref_size must be specified.");
    throw io_exception("reference or ref_size must be specified.");
  }

  if (!num_threads) {
    num_threads = 4;
    SPLOG_P(LOG_DEBUG, "kmerize_bf_params::validate> threads unspecified, setting to %lu",
            num_threads);
  }
}

namespace {

struct reader_chain {
  explicit reader_chain(const std::string& encoding, const file_info& chunk)
      : raw(chunk.file.read()),
        decoded(make_decoder(encoding, *raw)),
        reader(make_unique<kv_reader>(*decoded)) {}

  std::unique_ptr<readable> raw;
  std::unique_ptr<readable> decoded;
  std::unique_ptr<kv_reader> reader;
};

struct output_chain {
  void build(const path& root, const std::string& prefix) {
    sink = std::unique_ptr<kv_sink>(osp.build(root, prefix, the_manifest));
  }

  manifest the_manifest;
  output_stream_params osp;
  std::unique_ptr<kv_sink> sink;
};

std::map<std::string, size_t> meminfo() {
  std::map<std::string, size_t> info;
  std::ifstream fin("/proc/meminfo");
  for (std::string line; std::getline(fin, line);) {
    auto pos = line.find(':');
    auto name = line.substr(0, pos);
    auto value = line.substr(pos + 1);
    info.insert(std::make_pair(name, std::stol(value)));
  }
  return info;
}

class kmerizer {
 public:
  kmerizer(const path& root, const manifest& input, build_seqset::kmer_counter* counter,
           const kmerize_bf_params& params, manifest& count_out, manifest& overrep_out,
           kv_sink& histogram_sink, const progress_handler_t& on_progress)
      : m_root(root),
        m_input(input),
        m_counter(counter),
        m_params(params),
        m_count_out(count_out),
        m_histogram_sink(histogram_sink),
        m_on_progress(on_progress) {}

  void prepare();
  void run();

  std::unique_ptr<kmer_set> release_kmer_set() { return std::move(m_ks); }

 private:
  template <typename PassProcessor>
  void run_pass(progress_handler_t progress);

  path m_root;

  // input
  const manifest& m_input;
  using kmer_counter = build_seqset::kmer_counter;
  kmer_counter* m_counter = nullptr;
  kmerize_bf_params m_params;
  // output
  manifest& m_count_out;
  kv_sink& m_histogram_sink;

  progress_handler_t m_on_progress;

  size_t m_overrep_threshold;

  build_seqset::count_kmer_options m_options;
  std::unique_ptr<kmer_counter> m_owned_counter;

  std::mutex m_mu;
  size_t m_overrep_filter_count = 0;
  std::map<size_t, size_t> m_histogram;
  std::unique_ptr<kmer_set> m_ks;
};

void kmerizer::prepare() {
  SPLOG("kmerizer::prepare> using kmer_size: %zu", m_params.kmer_size);
  SPLOG("kmerizer::prepare> using error_rate: %0.2f%%", m_params.error_rate * 100.0);
  SPLOG("kmerizer::prepare> using num_threads: %lu", m_params.num_threads);
  // read one sample to determine:
  // - read length
  // - number of reads in record (is it paired?)
  auto it_chunks = m_input.begin();
  if (it_chunks == m_input.end()) {
    throw io_exception("kmerizer::prepare> input manifest has no chunks. Check input files.");
  }
  reader_chain chain(m_input.get_encoding(), *it_chunks);

  read_id key;
  unaligned_reads value;
  if (!chain.reader->read_msgpack(key, value)) {
    throw io_exception("kmerizer::prepare> malformed input");
  }

  auto num_parts = value.size();
  if (num_parts < 1) {
    throw io_exception("kmerizer::prepare> input dataset missing reads");
  }

  if (m_params.read_parts == 0) {
    m_params.read_parts = num_parts;
    SPLOG("kmerizer::prepare> detected read_parts: %zu", m_params.read_parts);
  } else {
    SPLOG("kmerizer::prepare> user specified read_parts: %zu", m_params.read_parts);
  }

  if (m_params.read_length == 0) {
    m_params.read_length = value[0].sequence.size();
    SPLOG("kmerizer::prepare> detected read_length: %zu", m_params.read_length);
  } else {
    SPLOG("kmerizer::prepare> user specified read_length: %zu", m_params.read_length);
  }

  if (m_params.memory_bound == 0) {
    const auto ONE_KB = 1024ULL;
    const auto ONE_GB = ONE_KB * 1024 * 1024;

    auto info = meminfo();

    // meminfo() returns kb, so convert to bytes
    auto free_bytes = (info["MemFree"] + info["Buffers"] + info["Cached"]) * ONE_KB;
    m_params.memory_bound = (free_bytes - (free_bytes * 0.10));

    // display is more natural in gigabytes
    SPLOG("kmerizer::prepare> detected available system memory: %llu GB",
          (size_t(free_bytes) / ONE_GB));
    SPLOG("kmerizer::prepare> using memory_bound of: %llu GB",
          (size_t(m_params.memory_bound) / ONE_GB));
  } else {
    SPLOG("kmerizer::prepare> user specified memory_bound: %zu GB", m_params.memory_bound);

    // Convert gigabytes to bytes
    m_params.memory_bound *= 1024 * 1024 * 1024ll;
  }

  if (m_params.ref_size == 0) {
    SPLOG("kmerizer::prepare> loading ref");
    reference ref(m_params.reference);
    SPLOG("kmerizer::prepare> ref loaded.");
    m_params.ref_size = ref.size();
    SPLOG("kmerizer::prepare> detected ref_size: %zu from reference: %s", m_params.ref_size,
          m_params.reference.c_str());
  } else {
    SPLOG("kmerizer::prepare> user specified ref_size: %zu", m_params.ref_size);
  }

  m_options.kmer_size = m_params.kmer_size;
  m_options.max_memory_bytes = m_params.memory_bound;
  // Don't waste memory and storage making a probibilistic table
  // that's going to be less than 1% full.
  m_options.max_prob_table_entries = m_params.ref_size * 100;
  m_overrep_threshold = m_params.overrep;
  if (!m_counter) {
    m_owned_counter.reset(new build_seqset::kmer_counter(m_options));
  }
}

template <typename PassProcessor>
void kmerizer::run_pass(progress_handler_t progress) {
  std::vector<file_info> file_infos;
  for (auto the_file_info : m_input) {
    file_infos.push_back(the_file_info);
  }

  parallel_for(  //
      0, file_infos.size(),
      [&](size_t i) {
        reader_chain chain(m_input.get_encoding(), file_infos[i]);
        read_id key;
        unaligned_reads value;

        PassProcessor processor(*m_counter);
        while (chain.reader->read_msgpack(key, value)) {
          for (const auto& read : value) {
            processor.add(read.sequence);
          }
        }
      },
      progress);
}

enum class filter_result {
  PASSED,           // Passed kmer filtering
  BELOW_MIN_COUNT,  // Kmer is below minimum kmer count
  SKEWED,           // Too skewed between fwd and revcomp
  NEAR_OVERREP,     // Too little coverage next to an overrepresented kmer
  OVERREP,          // Not filtered, but an overrepresented kmer.
};

const char* filter_result_str(filter_result fr) {
  switch (fr) {
    case filter_result::PASSED:
      return "PASSED";
    case filter_result::BELOW_MIN_COUNT:
      return "BELOW_MIN_COUNT";
    case filter_result::SKEWED:
      return "SKEWED";
    case filter_result::NEAR_OVERREP:
      return "NEAR_OVERREP";
    case filter_result::OVERREP:
      return "OVERREP";
  }
  return nullptr;
}

void kmerizer::run() {
  if (!m_counter) {
    // No existing counter supplied; we have to create our own and run
    // the first pass.
    m_counter = m_owned_counter.get();
    m_counter->start_prob_pass();
    run_pass<kmer_counter::prob_pass_processor>(subprogress(m_on_progress, 0, 0.2));
    m_counter->close_prob_pass();
  }

  subprogress exact_progress(m_on_progress, 0.2, 0.8);

  for (unsigned pass_num = 0; pass_num < m_counter->exact_passes(); ++pass_num) {
    subprogress pass_progress(exact_progress, pass_num * 1. / m_counter->exact_passes(),
                              (pass_num + 1.) / m_counter->exact_passes());
    m_counter->start_exact_pass(pass_num);
    run_pass<kmer_counter::exact_pass_processor>(pass_progress);
  }

  m_counter->close_exact_passes();

  overrep_map overrep(m_params.kmer_size);

  auto kmer_passes = [&](const kmer_counter::element& elem) -> filter_result {
    size_t tot_count = elem.fwd_count + elem.rev_count;
    if (tot_count < m_params.min_count) {
      return filter_result::BELOW_MIN_COUNT;
    }
    int32_t min_count = std::min(elem.fwd_count, elem.rev_count);
    float low =
        float(min_count + m_params.prior_count) / float(tot_count + 2 * m_params.prior_count);
    if (low < m_params.skew_cutoff) {
      // Too out of wack, forget it
      return filter_result::SKEWED;
    }

    if (m_overrep_threshold) {
      overrep_t o;
      if (overrep.find_near(elem.kmer, o)) {
        uint32_t min_c = std::min(elem.fwd_count, elem.rev_count);
        uint32_t max_c = std::max(elem.fwd_count, elem.rev_count);
        if (min_c < o.second * m_params.rnd_err_thresh &&
            max_c < o.second * m_params.sys_err_thresh) {
          std::lock_guard<std::mutex> l(m_mu);
          ++m_overrep_filter_count;
          return filter_result::NEAR_OVERREP;
        }
      }
    }

    return filter_result::PASSED;
  };

  size_t approx_kmer_table_size = 0;
  std::map<filter_result, std::map<size_t, size_t>> per_fr_histo;
  m_counter->extract_exact_counts(  //
      [&](kmer_counter::extract_iterator start, kmer_counter::extract_iterator limit) {
        size_t count = 0;
        if (!m_overrep_threshold) {
          std::lock_guard<std::mutex> l(m_mu);
          approx_kmer_table_size += (limit - start);
          return;
        }
        for (auto it = start; it != limit; ++it) {
          ++count;

          auto elem = *it;
          size_t tot_count = elem.fwd_count + elem.rev_count;
          if (tot_count < m_overrep_threshold) {
            continue;
          }

          std::lock_guard<std::mutex> l(m_mu);
          overrep.add_overrep(std::make_pair(elem.kmer, tot_count));
          per_fr_histo[filter_result::OVERREP][tot_count]++;
        }
        std::lock_guard<std::mutex> l(m_mu);
        approx_kmer_table_size += count;
      });

  std::map<filter_result, size_t> filter_result_counts;
  bool need_collect_stats = true;
  boost::optional<file_writer> kmers_out;
  if (!m_params.dump_kmers_file.empty()) {
    kmers_out.emplace(m_params.dump_kmers_file);
  }
  m_ks = make_unique<kmer_set>(
      approx_kmer_table_size, m_params.kmer_size, m_params.memory_bound,
      [&](const kmer_set::kmer_output_f& output_f, progress_handler_t kmer_progress) {
        bool collect_stats_this_pass = false;
        {
          std::lock_guard<std::mutex> l(m_mu);
          if (need_collect_stats) {
            collect_stats_this_pass = true;
            need_collect_stats = false;
          }
        }
        m_counter->extract_exact_counts(
            [&](kmer_counter::extract_iterator start, kmer_counter::extract_iterator limit) {
              std::map<size_t, size_t> local_histo;
              std::map<filter_result, size_t> local_result_counts;
              std::map<filter_result, std::unordered_map<size_t, size_t>> local_per_fr_histo;
              std::stringstream local_kmers_out;

              for (auto it = start; it != limit; ++it) {
                auto elem = *it;
                filter_result result = kmer_passes(elem);
                if (collect_stats_this_pass) {
                  local_result_counts[result]++;
                }
                if (result != filter_result::PASSED) {
                  if (collect_stats_this_pass) {
                    local_per_fr_histo[result][elem.fwd_count + elem.rev_count]++;
                  }
                  continue;
                }

                if (collect_stats_this_pass) {
                  local_histo[elem.fwd_count + elem.rev_count]++;
                }

                unsigned flags = 0;
                if (it->fwd_starts_read) {
                  flags |= kmer_set::k_fwd_starts_read;
                }
                if (it->rev_starts_read) {
                  flags |= kmer_set::k_rev_starts_read;
                }

                output_f(it->kmer, flags);
                if (!m_params.dump_kmers_file.empty()) {
                  kmer_t k = it->kmer;
                  local_kmers_out.write(reinterpret_cast<const char*>(&k), sizeof(kmer_t));
                }
              }

              if (kmers_out) {
                std::string str = local_kmers_out.str();
                std::lock_guard<std::mutex> l(m_mu);
                kmers_out->write(str.data(), str.size());
              }

              if (collect_stats_this_pass) {
                std::lock_guard<std::mutex> l(m_mu);
                for (const auto& histo_entry : local_histo) {
                  m_histogram[histo_entry.first] += histo_entry.second;
                }
                for (const auto& result : local_result_counts) {
                  filter_result_counts[result.first] += result.second;
                }
                for (const auto& fr : local_per_fr_histo) {
                  auto& target_map = per_fr_histo[fr.first];
                  for (const auto& item : fr.second) {
                    target_map[item.first] += item.second;
                  }
                }
              }
            });
        if (kmers_out) {
          kmers_out.reset();
        }
      },
      subprogress(m_on_progress, 0.8, 0.9));
  m_counter->close();

  if (overrep.size()) {
    SPLOG("Found %ld overrepresented kmers (%.2f%%)", overrep.size(),
          overrep.size() * 100. / approx_kmer_table_size);
  }

  SPLOG("%ld total kmers before filtering, resulting in:", approx_kmer_table_size);

  for (const auto& result : filter_result_counts) {
    SPLOG("  %-15s %10ld (%6.2f%%)", filter_result_str(result.first), result.second,
          result.second * 100. / approx_kmer_table_size);
  }

  std::string tmp_dir = CONF_S(storage_root);
  if (boost::starts_with(tmp_dir, "file://")) {
    tmp_dir = tmp_dir.substr(7);
  }
  std::ofstream fout(tmp_dir + "/kmer_quality_report.html");
  fout << kmer_quality_report_header;

  for (const auto& item : m_histogram) {
    m_histogram_sink.write_msgpack(item.first, item.second);
    fout << "{'x':" << item.first << ",'y':" << item.second << "},";
  }

  fout << kmer_quality_report_footer;
  fout.close();

  for (const auto& per_fr : per_fr_histo) {
    std::ofstream fout(tmp_dir + "/kmer_quality_report-" + filter_result_str(per_fr.first) +
                       ".html");
    fout << kmer_quality_report_header;

    for (const auto& item : per_fr.second) {
      fout << "{'x':" << item.first << ",'y':" << item.second << "},";
    }
    fout << kmer_quality_report_footer;
    fout.close();
  }

  m_count_out.metadata().set(meta::ns::readonly, "kmer_size", m_params.kmer_size);
  SPLOG("kmerization complete");
}

}  // namespace

void kmerize_bf_subtask::run() {
  SPLOG_P(LOG_DEBUG, "kmerize_bf_subtask::run> Entry");

  manifest count;

  output_chain histogram;
  histogram.osp.sort = "uint64";
  histogram.osp.presorted = true;
  histogram.build(get_root(), "kmerize_bf_histogram");
  manifest overrep_manifest;

  auto duration = stopwatch([&] {
    kmerizer work(get_root(), input, nullptr, params, count, overrep_manifest, *histogram.sink,
                  [&](double progress) { update_progress(progress); });
    work.prepare();
    work.run();
  });

  histogram.sink->close();

  SPLOG_P(LOG_DEBUG, "kmerize_bf_subtask::run> Took %ld ms.", duration.count());

  SPLOG("kmerize_bf_subtask::run> Writing final data manifest")

  std::vector<manifest> outputs = {count, histogram.the_manifest, overrep_manifest};
  set_output(outputs);
}

std::pair<std::unique_ptr<kmer_set>, std::vector<manifest>> run_kmerize_subtask(
    const kmerize_bf_params& params, const manifest& input, build_seqset::kmer_counter* counter,
    progress_handler_t update_progress) {
  SPLOG_P(LOG_DEBUG, "run_kmerize_subtask> Entry");

  manifest count;

  output_chain histogram;
  histogram.osp.sort = "uint64";
  histogram.osp.presorted = true;
  histogram.build(CONF_S(path_bulkdata), "kmerize_bf_histogram");
  manifest overrep_manifest;

  std::unique_ptr<kmer_set> ks;
  auto duration = stopwatch([&] {
    kmerizer work(CONF_S(path_bulkdata), input, counter, params, count, overrep_manifest,
                  *histogram.sink, [&](double progress) { update_progress(progress); });
    work.prepare();
    work.run();

    ks = work.release_kmer_set();
  });

  CHECK(ks);

  histogram.sink->close();

  SPLOG_P(LOG_DEBUG, "run_kmerize_subtask> Took %ld ms.", duration.count());

  SPLOG("run_kmerize_subtask> Writing final data manifest")

  std::vector<manifest> outputs = {count, histogram.the_manifest, overrep_manifest};
  return std::make_pair(std::move(ks), outputs);
}

void kmerize_bf_task::run() {
  if (m_state == state::kmerize) {
    SPLOG("kmerize_bf_task::run> kmerize");

    split_progress(0.01, 0.50);

    auto kmerize = make_unique<kmerize_bf_subtask>();
    kmerize->input = input;
    kmerize->params = params;
    m_subtask = add_subtask(std::move(kmerize));

    m_state = state::sort_from_kmerize;
  } else if (m_state == state::sort) {
    // Sort from already-provided kmer counts
    SPLOG("kmerize_bf_task::run> sort");

    split_progress(0.01, 0.01);

    auto sort_kmers = make_unique<sort_task>();
    sort_kmers->input = m_kmer_counts;
    m_subtask = add_subtask(std::move(sort_kmers));

    m_state = state::final;
  } else if (m_state == state::sort_from_kmerize) {
    SPLOG("kmerize_bf_task::run> sort_from_kmerize");

    split_progress(0.01, 0.01);

    std::vector<manifest> outputs;
    get_output(outputs, m_subtask);
    m_histogram = outputs[1];
    m_overrep = outputs[2];
    auto sort_kmers = make_unique<sort_task>();
    sort_kmers->input = outputs[0];
    m_subtask = add_subtask(std::move(sort_kmers));

    m_state = state::final;
  } else if (m_state == state::final) {
    SPLOG("kmerize_bf_task::run> final");

    manifest sorted;
    get_output(sorted, m_subtask);
    sorted.metadata().set(meta::ns::readonly, "kmer_size", params.kmer_size);

    std::vector<manifest> outputs = {sorted, m_histogram, m_overrep};
    set_output(outputs);
  } else if (m_state == state::do_nothing) {
    SPLOG("kmerize_bf_task::run> do_nothing");
    std::vector<manifest> outputs = {m_kmer_counts, m_histogram, m_overrep};
    set_output(outputs);
  }
}

std::vector<std::string> get_kmer_filter_result_types() {
  std::vector<std::string> result;
  for (int i = 0;; ++i) {
    const char* s = filter_result_str(filter_result(i));
    if (!s) {
      return result;
    }
    result.push_back(s);
  }
};
