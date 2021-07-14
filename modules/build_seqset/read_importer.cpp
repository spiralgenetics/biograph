#include <htslib/bgzf.h>
#include <htslib/hfile.h>
#include <htslib/sam.h>
#include <sys/syscall.h>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/filesystem.hpp>
#include <boost/optional.hpp>
#include <boost/regex.hpp>

#include "modules/bio_base/dna_base_set.h"
#include "modules/bio_base/reference.h"
#include "modules/bio_format/fastq.h"
#include "modules/build_seqset/read_importer.h"
#include "modules/io/config.h"
#include "modules/io/defaults.h"
#include "modules/io/parallel.h"
#include "modules/io/zip.h"

namespace fs = boost::filesystem;
using boost::algorithm::ends_with;

static constexpr bool k_keep_quality_scores = false;

namespace build_seqset {

constexpr size_t read_importer_base::k_default_num_recs;
constexpr size_t read_importer_base::k_bam_line_batch_size;
constexpr unsigned read_importer_base::k_hts_threads;

namespace {

class fastq_batcher {
 public:
  fastq_batcher(fastq_reader& fq) : m_fq(fq) {
    m_thread = std::async(std::launch::async, [this]() {
      try {
        run_bg_thread();
      } catch (const io_exception& e) {
        std::cerr << e.what() << "\n";
        exit(1);
      }
    });
    fill_buffer();
  }
  ~fastq_batcher() {
    m_aborted = true;
    m_buffer_empty.notify_all();
    m_thread.get();
  }

  void run_bg_thread() {
    std::vector<std::pair<read_id, unaligned_read>> new_reads;
    new_reads.reserve(k_batch_size);
    bool got_eof = false;
    for (;;) {
      CHECK(new_reads.empty());
      while (new_reads.size() < k_batch_size) {
        new_reads.emplace_back();
        auto& r = new_reads.back();
        if (!m_fq.read(r.first, r.second)) {
          new_reads.pop_back();
          got_eof = true;
          break;
        }
      }

      std::unique_lock<std::mutex> l(m_mu);
      while (!m_aborted && !m_reads.empty()) {
        m_buffer_empty.wait(l);
      }
      std::swap(m_reads, new_reads);
      m_buffer_full.notify_one();
      if (got_eof) {
        m_eof = true;
        return;
      }
      if (m_aborted) {
        return;
      }
    }
  }

  bool read(read_id& id, unaligned_read& read) {
    if (m_read_pos == m_reads.end()) {
      return false;
    }
    id = std::move(m_read_pos->first);
    read = std::move(m_read_pos->second);
    ++m_read_pos;
    if (m_read_pos == m_reads.end()) {
      fill_buffer();
    }
    return true;
  }

  void fill_buffer() {
    std::unique_lock<std::mutex> l(m_mu);
    m_reads.clear();
    m_buffer_empty.notify_one();
    while (!m_eof && !m_aborted && m_reads.empty()) {
      m_buffer_full.wait(l);
    }
    m_read_pos = m_reads.begin();
  }

 private:
  static constexpr unsigned k_batch_size = 1024;
  fastq_reader& m_fq;

  std::future<void> m_thread;
  std::mutex m_mu;
  std::condition_variable m_buffer_empty;
  std::condition_variable m_buffer_full;
  bool m_aborted = false;
  bool m_eof = false;
  std::vector<std::pair<read_id, unaligned_read>> m_reads;
  std::vector<std::pair<read_id, unaligned_read>>::iterator m_read_pos;
};

}  // namespace

constexpr unsigned read_importer_base::read_batch::k_read_batch_size;

read_importer_base::read_importer_base(progress_handler_t progress) : m_progress(progress) {}

void read_importer_base::add_queued_import(double progress_part,
                                           const std::function<void(void)>& import_f) {
  thread_pool::work_t work{[import_f](parallel_state&) { import_f(); }};
  work.progress_part = progress_part;
  m_queued_imports.emplace_back(std::move(work));
}

size_t read_importer_base::import() {
  CHECK(!m_queued_imports.empty());

  std::sort(m_queued_imports.begin(), m_queued_imports.end(),
            [](const thread_pool::work_t& rhs, const thread_pool::work_t& lhs) {
              // Start bigger chunks first since they might take longer.
              return lhs.progress_part > rhs.progress_part;
            });

  parallel_pool().execute_worklist(m_queued_imports, subprogress(m_progress, 0, 0.9));
  m_queued_imports.clear();

  SPLOG("Imports done.  Cleaning up...");

  parallel_for(
      0, m_queued_cleanups.size(), [this](size_t idx) { m_queued_cleanups[idx](); },
      subprogress(m_progress, 0.9, 1));
  m_queued_cleanups.clear();

  m_progress(1);

  return m_total_reads;
}

void read_importer_base::read_batch::cut_reads(unsigned start, unsigned end) {
  CHECK_GT(end, start);
  for (auto& p : m_reads) {
    for (auto& r : p.second) {
      unsigned this_end = end;
      if (this_end > r.sequence.size()) {
        this_end = r.sequence.size();
      }
      CHECK_GT(this_end, start);
      r.sequence = r.sequence.substr(start, this_end - start);
    }
  }
}

void read_importer_base::flush_read_batch(read_batch& batch) {
  if (m_cut_reads_end) {
    // Trim down the reads before flushing the batch
    batch.cut_reads(m_cut_reads_start, m_cut_reads_end);
  }
  parallel_state* st = parallel_pool().get_state();
  CHECK(st);
  process_read_batch(*st, batch);
  batch.clear();
}

size_t read_importer_base::bam_process_line(bam_pair_cache_t& pair_cache, bam1_t& line,
                                            read_batch& batch, bool& got_paired) {
  if (line.core.flag & BAM_FSECONDARY || line.core.flag & BAM_FSUPPLEMENTARY) {
    return 0;
  }

  if (line.core.flag & BAM_FPAIRED) {
    got_paired = true;
    std::string qname(bam_get_qname(&line));
    auto it = pair_cache.find(qname);
    if (it != pair_cache.end()) {
      batch.add_paired_read(std::move(qname), bam1_to_unaligned_read(line), std::move(it->second));
      pair_cache.erase(it);
    } else {
      pair_cache.emplace(std::move(qname), bam1_to_unaligned_read(line));
    }
  } else {
    batch.add_unpaired_read(bam_get_qname(&line), bam1_to_unaligned_read(line));
  }

  return 1;
}

std::unique_ptr<read_importer_base::bam_pair_cache_t> read_importer_base::merge_bam_pair_cache(
    std::unique_lock<std::mutex>& l, std::unique_ptr<bam_pair_cache_t> local_pair_cache,
    std::unique_ptr<bam_pair_cache_t> total_pair_cache, read_batch& batch) {
  CHECK(!l);
  if (local_pair_cache->size() > total_pair_cache->size()) {
    std::swap(local_pair_cache, total_pair_cache);
  }
  if (!local_pair_cache->empty()) {
    auto local_it = local_pair_cache->begin();
    auto local_next = local_it;
    for (; local_it != local_pair_cache->end(); local_it = local_next) {
      ++local_next;
      const auto& qname = local_it->first;
      auto& ur = local_it->second;

      auto tot_it = total_pair_cache->find(qname);
      if (tot_it != total_pair_cache->end()) {
        batch.add_paired_read(qname, std::move(ur), std::move(tot_it->second));
        total_pair_cache->erase(tot_it);
      } else {
        total_pair_cache->emplace(qname, std::move(ur));
      }
      local_pair_cache->erase(local_it);
    }
  }
  l.lock();
  CHECK(local_pair_cache->empty());
  m_free_pair_caches.emplace_back(std::move(local_pair_cache));
  return total_pair_cache;
}

unaligned_read read_importer_base::bam1_to_unaligned_read(const bam1_t& line) {
  unaligned_read ret;
  if (line.core.flag & BAM_FREAD1) {
    ret.pair_number = 0;
  } else if (line.core.flag & BAM_FREAD2) {
    ret.pair_number = 1;
  }
  ret.sequence.resize(line.core.l_qseq);
  const uint8_t* seq = bam_get_seq(&line);
  for (int i = 0; i < line.core.l_qseq; i++) {
    char base = seq_nt16_str[bam_seqi(seq, i)];
    ret.sequence[i] = base;
  }

  if (k_keep_quality_scores) {
    const uint8_t* qual = bam_get_qual(&line);
    ;
    ret.quality.reserve(line.core.l_qseq);
    for (int i = 0; i < line.core.l_qseq; i++) {
      ret.quality.push_back(33 + qual[i]);
    }
  }

  if (line.core.flag & BAM_FREVERSE) {
    reverse_complement_iupac_string(ret.sequence);
    if (k_keep_quality_scores) {
      std::reverse(ret.quality.begin(), ret.quality.end());
    }
  }
  return ret;
}

void read_importer_base::submit_bam_line_batch(std::unique_ptr<bam_line_batch> line_batch,
                                               size_t num_lines, double progress_part,
                                               size_t batch_num, bam_file_state* st) {
  bam_line_batch* ptr = line_batch.release();

  thread_pool::work_t submit_work{{[this, ptr, num_lines, batch_num, st](parallel_state&) {
    std::unique_ptr<bam_line_batch> line_batch;
    line_batch.reset(ptr);

    read_batch batch(this);
    bool got_paired = false;

    std::unique_ptr<bam_pair_cache_t> local_pair_cache;
    {
      std::lock_guard<std::mutex> l(m_mu);
      if (!m_free_pair_caches.empty()) {
        local_pair_cache = std::move(m_free_pair_caches.back());
        m_free_pair_caches.pop_back();
      }
    }
    if (!local_pair_cache) {
      local_pair_cache = make_unique<bam_pair_cache_t>(track_alloc("bam_pair_cache"));
    }
    CHECK(local_pair_cache->empty());

    size_t read_count = 0;
    CHECK_LE(num_lines, line_batch->size());
    auto end = line_batch->begin() + num_lines;
    for (auto it = line_batch->begin(); it != end; ++it) {
      auto& line = *it;
      read_count += bam_process_line(*local_pair_cache, line, batch, got_paired);
    }

    std::unique_lock<std::mutex> l(m_mu);
    bool did_insert = st->pending.emplace(batch_num, std::move(local_pair_cache)).second;
    CHECK(did_insert);
    if (got_paired) {
      m_got_paired = true;
    }
    m_total_reads += read_count;
    m_free_bam_line_batches.push_back(std::move(line_batch));
    consolidate_pair_cache(l, st, batch);
  }}};
  submit_work.progress_part = progress_part;
  parallel_pool().add_work_async(std::move(submit_work));
}

void read_importer_base::consolidate_pair_cache(std::unique_lock<std::mutex>& l, bam_file_state* st,
                                                read_batch& batch) {
  CHECK(l);

  while (consolidate_pair_cache_once(l, st, batch)) {
  }
}

bool read_importer_base::consolidate_pair_cache_once(std::unique_lock<std::mutex>& l,
                                                     bam_file_state* st, read_batch& batch) {
  CHECK(l);

  if (st->pending.empty()) {
    return false;
  }

  auto first_it = st->pending.begin();
  if (first_it == st->pending.end() || first_it->first != st->merge_position) {
    return false;
  }

  auto second_it = first_it;
  ++second_it;
  if (second_it == st->pending.end() || second_it->first != (st->merge_position + 1)) {
    return false;
  }

  std::unique_ptr<bam_pair_cache_t> pc1 = std::move(first_it->second);
  st->pending.erase(first_it);
  std::unique_ptr<bam_pair_cache_t> pc2 = std::move(second_it->second);
  st->pending.erase(second_it);
  l.unlock();

  std::unique_ptr<bam_pair_cache_t> merged =
      merge_bam_pair_cache(l, std::move(pc1), std::move(pc2), batch);
  CHECK(l);

  ++st->merge_position;
  bool did_insert = st->pending.emplace(st->merge_position, std::move(merged)).second;
  CHECK(did_insert);
  return true;
}

std::unique_ptr<read_importer_base::bam_line_batch> read_importer_base::get_bam_line_batch() {
  std::unique_ptr<bam_line_batch> result;
  {
    std::lock_guard<std::mutex> l(m_mu);
    if (!m_free_bam_line_batches.empty()) {
      result = std::move(m_free_bam_line_batches.back());
      m_free_bam_line_batches.pop_back();
    }
  }
  if (!result) {
    result = make_unique<bam_line_batch>();
    result->resize(k_bam_line_batch_size);
    memset(result->data(), 0, sizeof(*result->data()) * result->size());
  }
  return result;
}

void read_importer_base::free_bam_line_batches() {
  std::lock_guard<std::mutex> l(m_mu);
  for (auto& line_batch : m_free_bam_line_batches) {
    for (auto& line : *line_batch) {
      if (line.data) {
        free(line.data);
      }
    }
    line_batch.reset();
  }
  m_free_bam_line_batches.clear();
}

void read_importer_base::bam_output_unpaired() {
  std::unique_lock<std::mutex> l(m_mu);
  std::vector<std::unique_ptr<bam_pair_cache_t>> tot_caches;
  for (const auto& st : m_bam_files) {
    CHECK_LT(st->num_batches, std::numeric_limits<size_t>::max());
    if (st->pending.empty()) {
      CHECK_EQ(0, st->num_batches);
      continue;
    }
    CHECK_EQ(st->merge_position + 1, st->num_batches);
    CHECK_EQ(1, st->pending.size());

    std::unique_ptr<bam_pair_cache_t> pc = std::move(st->pending.begin()->second);
    if (!pc->empty()) {
      SPLOG("File '%s' contains %lu pairs lacking mates", st->filename.c_str(), pc->size());
      tot_caches.emplace_back(std::move(pc));
    }
    st->pending.clear();
  }
  m_bam_files.clear();
  l.unlock();

  if (tot_caches.empty()) {
    return;
  }

  while (tot_caches.size() > 1) {
    size_t tot_entries = 0;
    for (const auto& pc : tot_caches) {
      tot_entries += pc->size();
    }
    SPLOG("Merging down final %lu pair caches with %lu entries", tot_caches.size(), tot_entries);

    std::vector<std::unique_ptr<bam_pair_cache_t>> new_tot_caches;
    size_t new_size = (tot_caches.size() + 1) / 2;
    new_tot_caches.resize(new_size);

    parallel_for(0, new_size, [&](size_t new_idx) {
      read_batch batch(this);
      std::unique_lock<std::mutex> local_l(m_mu, std::defer_lock);

      size_t old_idx1 = new_idx * 2;
      size_t old_idx2 = new_idx * 2 + 1;

      CHECK_LT(old_idx1, tot_caches.size());
      std::unique_ptr<bam_pair_cache_t> pc1 = std::move(tot_caches[old_idx1]);
      std::unique_ptr<bam_pair_cache_t> merged;
      if (old_idx2 == tot_caches.size()) {
        merged = std::move(pc1);
      } else {
        CHECK_LT(old_idx2, tot_caches.size());
        std::unique_ptr<bam_pair_cache_t> pc2 = std::move(tot_caches[old_idx2]);
        merged = merge_bam_pair_cache(local_l, std::move(pc1), std::move(pc2), batch);
        CHECK(local_l);
        local_l.unlock();
      }
      CHECK_LT(new_idx, new_size);
      new_tot_caches[new_idx] = std::move(merged);
    });
    std::swap(tot_caches, new_tot_caches);
  }

  CHECK_EQ(tot_caches.size(), 1);
  std::unique_ptr<bam_pair_cache_t> pair_cache = std::move(tot_caches.front());
  CHECK(pair_cache);
  if (pair_cache->empty()) {
    SPLOG("Final pair cache empty; all paired reads successfully matched with their mates.");
    return;
  }

  SPLOG(
      "WARNING: %lu entries remaining in pair cache from paired reads; treating as unpaired reads",
      pair_cache->size());
  size_t total_count = 0;
  parallel_for(0, pair_cache->bucket_count(),
               [this, &total_count, &pair_cache](size_t start, size_t limit) {
                 read_batch batch(this);
                 size_t local_count = 0;
                 for (size_t bucket = start; bucket != limit; ++bucket) {
                   auto it = pair_cache->begin(bucket);
                   auto end = pair_cache->end(bucket);

                   for (; it != end; ++it) {
                     ++local_count;
                     batch.add_unpaired_read(it->first, std::move(it->second));
                   }
                 }
                 std::lock_guard<std::mutex> l(m_mu);
                 total_count += local_count;
               });
  if (total_count) {
    SPLOG("Completed saving %lu unexpected unpaired reads", total_count);
  }
}

// Process a bam file, writing results to sink. Return the record count.
void read_importer_base::queue_bam(const std::string& in_file, const std::string& ref_dir) {
  size_t file_size = 0;
  if (fs::is_regular_file(in_file)) {
    file_size = fs::file_size(in_file);
  }
  // Open input bam
  hFILE* hfile_in = hopen(in_file.c_str(), "r");
  if (!hfile_in) {
    throw io_exception("Unable to open file " + in_file);
  }

  htsFile* bam_in = hts_hopen(hfile_in, in_file.c_str(), "r");
  if (bam_in == NULL) {
    throw io_exception("Unable to read file " + in_file);
  }

  // CRAM support
  const htsFormat* fmt = hts_get_format(bam_in);
  std::string cram_reference;
  if (fmt->format == cram) {
    cram_reference = std::string(ref_dir) + "/" + defaults.original_fasta;
  }

  if (!cram_reference.empty()) {
    hts_set_opt(bam_in, CRAM_OPT_REFERENCE, cram_reference.c_str());
  }

  // Read bam header
  bam_hdr_t* header = sam_hdr_read(bam_in);
  if (header == NULL) {
    throw std::runtime_error(in_file + " is not a valid BAM file.");
  }

  m_queued_cleanups.push_back([this]() { bam_output_unpaired(); });
  m_queued_cleanups.push_back([this]() { free_bam_line_batches(); });

  m_bam_files.emplace_back(make_unique<bam_file_state>());
  bam_file_state* file_state = m_bam_files.back().get();
  file_state->filename = in_file;

  add_queued_import(
      k_default_num_recs, [this, header, hfile_in, bam_in, in_file, file_size, file_state]() {
        unsigned hts_threads = get_thread_count() / 4;
        if (hts_threads > k_hts_threads) {
          hts_threads = k_hts_threads;
        }
        if (hts_threads > 0) {
          hts_set_threads(bam_in, hts_threads);
        }

        std::unique_ptr<bam_line_batch> line_batch = get_bam_line_batch();
        CHECK(!line_batch->empty());
        size_t record_count = 0;
        auto line_batch_it = line_batch->begin();
        DCHECK(line_batch_it != line_batch->end());

        size_t last_off = 0;
        size_t batch_num = 0;
        for (;;) {
          int read_result = sam_read1(bam_in, header, &*line_batch_it);

          if (read_result == -1) {
            // EOF
            break;
          }
          if (read_result < 0) {
            throw(io_exception(printstring("sam_read1 returned %d when reading %s", read_result,
                                           in_file.c_str())));
          }
          ++record_count;
          ++line_batch_it;
          if (line_batch_it == line_batch->end()) {
            double progress_part = double(k_bam_line_batch_size) / k_default_num_recs;
            size_t off = htell(hfile_in);
            if (off > last_off && off < file_size) {
              progress_part = double(off - last_off) / file_size;
            }
            last_off = off;
            size_t num_lines = line_batch_it - line_batch->begin();
            submit_bam_line_batch(std::move(line_batch), num_lines, progress_part, batch_num++,
                                  file_state);
            line_batch = get_bam_line_batch();
            line_batch_it = line_batch->begin();
            DCHECK(line_batch_it != line_batch->end());
          }
        }

        if (!line_batch->empty()) {
          size_t num_lines = line_batch_it - line_batch->begin();
          submit_bam_line_batch(std::move(line_batch), num_lines, 1. /* whole rest of the input */,
                                batch_num++, file_state);
        }

        bam_hdr_destroy(header);
        hts_close(bam_in);

        if (record_count) {
          SPLOG("%s: completed reading %lu records", in_file.c_str(), record_count);
        } else {
          SPLOG("WARNING: %s: no records present", in_file.c_str());
        }

        std::unique_lock<std::mutex> l(m_mu);
        file_state->num_batches = batch_num;
        read_batch batch(this);
        consolidate_pair_cache(l, file_state, batch);
      });
}

// Process FASTQ files, writing results to an import out handle. Return the
// record count.
void read_importer_base::queue_fastq(const std::string& in_file, const std::string& in_file2,
                                     const bool& interleaved) {
  add_queued_import(k_default_num_recs, [this, in_file, in_file2, interleaved]() {
    size_t read_count = read_fastq(in_file, in_file2, interleaved);
    std::lock_guard<std::mutex> l(m_mu);
    m_total_reads += read_count;
  });
}

size_t read_importer_base::read_fastq(const std::string& in_file, const std::string& in_file2,
                                      const bool& interleaved) {
  m_got_paired = interleaved or (not in_file2.empty());

  size_t count = 0;

  read_id id, id2;
  unaligned_reads ur, ur2;

  std::unique_ptr<file_reader> fin, fin2;
  std::unique_ptr<fastq_reader> fq, fq2;
  std::unique_ptr<zip_reader> unzip, unzip2;
  std::unique_ptr<fastq_batcher> fqb, fqb2;

  size_t file_size;

  if (fs::is_regular_file(in_file)) {
    try {
      file_size = fs::file_size(in_file);
    } catch (fs::filesystem_error& e) {
      file_size = std::numeric_limits<uint64_t>::max();
    }
  } else {
    file_size = std::numeric_limits<uint64_t>::max();
  }

  fin = make_unique<file_reader>(in_file);
  if (ends_with(in_file, ".gz")) {
    unzip = make_unique<zip_reader>(*fin);
    fq = make_unique<fastq_reader>(*unzip, false /* don't keep quality scores */);
  } else {
    fq = make_unique<fastq_reader>(*fin, false /* don't keep quality scores */);
  }
  fqb = make_unique<fastq_batcher>(*fq);

  if (not in_file2.empty()) {
    if (interleaved) {
      throw std::runtime_error("Interleaved reads must all be stored in one FASTQ file.");
    }
    fin2 = make_unique<file_reader>(in_file2);
    if (ends_with(in_file2, ".gz")) {
      unzip2 = make_unique<zip_reader>(*fin2);
      fq2 = make_unique<fastq_reader>(*unzip2, false /* don't keep quality scores */);
    } else {
      fq2 = make_unique<fastq_reader>(*fin2, false /* don't keep quality scores */);
    }
    fqb2 = make_unique<fastq_batcher>(*fq2);
  }

  std::unique_ptr<read_batch> batch;
  double last_progress = 0;
  auto flush_batch_if_needed = [&]() {
    if (batch && batch->full()) {
      read_batch* batch_ptr = batch.release();

      thread_pool::work_t work{[batch_ptr](parallel_state&) {
        // Batch is processed on destruct.
        delete batch_ptr;
      }};
      double cur_progress = fin->pos() / double(file_size);
      if (cur_progress > 1) {
        cur_progress = 1;
      }
      if (cur_progress > last_progress) {
        work.progress_part = (cur_progress - last_progress) / (1 - last_progress);
        last_progress = cur_progress;
      } else {
        work.progress_part = 0;
      }
      parallel_pool().add_work_async(work);
    }

    if (!batch) {
      batch = make_unique<read_batch>(this);
    }
  };

  for (;;) {
    flush_batch_if_needed();
    auto& id = batch->add_id();
    auto& ur = batch->add_read();
    if (!fqb->read(id, ur)) {
      batch->unadd_id();
      break;
    }
    count++;

    if (interleaved) {
      auto& ur2 = batch->add_read();
      if (!fqb->read(id2, ur2)) {
        SPLOG("Warning: interleaved fastq specified, but read an odd number of reads.");
        batch->unadd_id();
        break;
      }
      // We could check that id == id2, but there is no guarantee that the id
      // follows any convention whatsoever, even up to the first space.
      // Could be 0 1, 1 2, /1 /2, -1 -2, _0 _1, or no difference at all.
      count++;
    }

    else if (fqb2) {
      auto& ur2 = batch->add_read();
      if (fqb2->read(id2, ur2)) {
        count++;
      } else {
        batch->unadd_read();
        fqb2.reset();
      }
    }
  }

  // in_file2 may contain more unpaired reads
  if (fqb2) {
    for (;;) {
      flush_batch_if_needed();

      auto& id = batch->add_id();
      auto& ur = batch->add_read();
      if (!fqb2->read(id, ur)) {
        batch->unadd_id();
        break;
      }
      count++;
    }
  }

  return count;
}

void read_importer_base::map_bam_contigs_to_ref(const bam_hdr_t* bam_header,
                                                const reference& the_ref,
                                                const std::string& in_file,
                                                const std::string& ref_dir) {
  SPLOG("read_importer_base::map_bam_contigs_to_ref> Reference fasta path: %s",
        the_ref.fasta_path().c_str());
  for (int i = 0; i < bam_header->n_targets; i++) {
    scaffold key(bam_header->target_name[i], 0, 0);
    auto it = the_ref.get_assembly().scaffolds.find(key);
    if (it == the_ref.get_assembly().scaffolds.end()) {
      SPLOG("Did not find contig \"%s\" with index %d in reference.", bam_header->target_name[i],
            i);
      throw io_exception(boost::format("Contig \"%1%\" was found in the input, \"%2%\" "
                                       "but not in the reference, \"%3%\"") %
                         bam_header->target_name[i] % in_file % ref_dir);
      m_remap_contigs[i] = -1;
    } else {
      m_remap_contigs[i] = it->index;
      SPLOG(
          "Mapped input contig \"%s\" with index %d to reference contig index "
          "%d",
          bam_header->target_name[i], i, it->index);
    }
  }
}

void read_importer_base::read_batch::add_paired_read(std::string qname, unaligned_read rd1,
                                                     unaligned_read rd2) {
  if (full()) {
    m_importer->flush_read_batch(*this);
  }
  m_reads.emplace_back();
  auto& ur = m_reads.back();
  ur.first.pair_name = std::move(qname);
  ur.second.emplace_back(std::move(rd1));
  ur.second.emplace_back(std::move(rd2));
}

void read_importer_base::read_batch::add_unpaired_read(std::string qname, unaligned_read rd) {
  if (full()) {
    m_importer->flush_read_batch(*this);
  }
  m_reads.emplace_back();
  auto& ur = m_reads.back();
  ur.first.pair_name = std::move(qname);
  ur.second.emplace_back(std::move(rd));
}

}  // namespace build_seqset
