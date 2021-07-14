#pragma once

#include "modules/bio_base/unaligned_read.h"
#include "modules/bio_format/fastq.h"
#include "modules/io/keyvalue.h"
#include "modules/io/parallel.h"
#include "modules/io/progress.h"
#include "modules/io/track_mem.h"

#include <htslib/sam.h>

#include <map>

class reference;

namespace build_seqset {

class read_importer_base {
 public:
  // Number of records in a file to default to if we can't tell, for progress purposes.
  static constexpr size_t k_default_num_recs = 10000000;

  // Number of records to process in a bam line batch.  This should be
  // large enough that we get reasonable hits in the pair cache.
  static constexpr size_t k_bam_line_batch_size = 32768;

  // Maximum number of threads for htslib to use for decompression.  8
  // seems to be about the right amount to keep up with the single
  // thread bottleneck of calling sam_read1.
  static constexpr unsigned k_hts_threads = 8;

  read_importer_base(progress_handler_t progress = null_progress_handler);
  read_importer_base(const read_importer_base&) = delete;

  void set_cut_region(unsigned start, unsigned end) {
    CHECK_GT(end, start);
    m_cut_reads_start = start;
    m_cut_reads_end = end;
  }

  bool got_paired() const { return m_got_paired; }

  void queue_bam(const std::string& in_file, const std::string& ref_dir);
  void queue_fastq(const std::string& in_file, const std::string& in_file2,
                   const bool& interleaved = false);

  // Executes queued imports and returns the total number of reads imported.
  size_t import();

 protected:
  using bam_pair_cache_t = tracked_unordered_map<std::string, unaligned_read>;
  struct bam_file_state {
    std::string filename;

    // Pending pair caches
    std::map<size_t /* batch number */, std::unique_ptr<bam_pair_cache_t>> pending;

    // Total number of batches in file, if known.
    size_t num_batches = std::numeric_limits<size_t>::max();

    // Earliest batch in pending that hasn't been merged.
    size_t merge_position = 0;
  };

  class read_batch {
   public:
    static constexpr unsigned k_read_batch_size = 1024;

    read_batch() = delete;
    read_batch(const read_batch&) = delete;
    read_batch(read_importer_base* importer) : m_importer(importer) {
      m_reads.reserve(k_read_batch_size);
    }

    void flush_and_reserve(read_batch& rhs) {
      std::swap(m_reads, rhs.m_reads);
      CHECK_EQ(m_importer, rhs.m_importer);
      m_reads.reserve(rhs.m_reads.capacity());
    }

    ~read_batch() {
      if (!empty()) {
        CHECK(m_importer);
        m_importer->flush_read_batch(*this);
        CHECK(empty());
      }
    }

    read_id& add_id() {
      CHECK(m_importer);
      if (full()) {
        m_importer->flush_read_batch(*this);
        CHECK(empty());
      }
      m_reads.emplace_back();
      // Make sure we don't have to reallocate if we have pairing information.
      m_reads.back().second.reserve(2);
      return m_reads.back().first;
    }

    void unadd_id() {
      CHECK(!m_reads.empty());
      m_reads.pop_back();
    }

    unaligned_read& add_read() {
      CHECK(!m_reads.empty());
      CHECK_LT(m_reads.back().second.size(), 2);
      m_reads.back().second.emplace_back();
      return m_reads.back().second.back();
    }

    void unadd_read() {
      CHECK(!m_reads.empty());
      CHECK(!m_reads.back().second.empty());
      m_reads.back().second.pop_back();
    }

    void add_paired_read(std::string qname, unaligned_read rd1, unaligned_read rd2);
    void add_unpaired_read(std::string qname, unaligned_read rd);

    const std::vector<std::pair<read_id, unaligned_reads>>& reads() const { return m_reads; }

    void cut_reads(unsigned start, unsigned end);

    void clear() { m_reads.clear(); }

    bool empty() const {
      CHECK(m_importer);
      return m_reads.empty();
    }

    bool full() {
      CHECK_LE(m_reads.size(), k_read_batch_size);
      return m_reads.size() == k_read_batch_size;
    }
   private:

    std::vector<std::pair<read_id, unaligned_reads>> m_reads;
    read_importer_base* m_importer = nullptr;
  };

  void flush_read_batch(read_batch& batch);
  unaligned_read bam1_to_unaligned_read(const bam1_t& line);
  void map_bam_contigs_to_ref(const bam_hdr_t* bam_header, const reference& the_ref,
                              const std::string& in_file, const std::string& ref_dir);
  size_t bam_process_line(bam_pair_cache_t& pair_cache, bam1_t& line, read_batch& batch,
                          bool& got_paired);
  std::unique_ptr<bam_pair_cache_t> merge_bam_pair_cache(
      std::unique_lock<std::mutex>& l, std::unique_ptr<bam_pair_cache_t> local_pair_cache,
      std::unique_ptr<bam_pair_cache_t> total_pair_cache, read_batch& batch);
  void bam_output_unpaired();
  void consolidate_pair_cache(std::unique_lock<std::mutex>& l, bam_file_state* st,
                              read_batch& batch);
  // Returns true if any progress was made.
  bool consolidate_pair_cache_once(std::unique_lock<std::mutex>& l, bam_file_state* st,
                                   read_batch& batch);

  virtual void process_read_batch(parallel_state& st, read_batch&) = 0;

  size_t read_bam(const std::string& in_file, const std::string& ref_dir);
  size_t read_fastq(const std::string& in_file, const std::string& in_file2,
                    const bool& interleaved = false);

  void add_queued_import(double progress_part, const std::function<void(void)>& import_f);

  using bam_line_batch = std::vector<bam1_t>;
  std::unique_ptr<bam_line_batch> get_bam_line_batch();
  void submit_bam_line_batch(std::unique_ptr<bam_line_batch> line_batch, size_t num_lines,
                             double progress_part, size_t batch_num, bam_file_state* st);
  void free_bam_line_batches();

  progress_handler_t m_progress = null_progress_handler;

  // Main import tasks queued
  std::vector<thread_pool::work_t> m_queued_imports;

  // Cleanup tasks queued after all imports have been completed.
  std::vector<std::function<void(void)>> m_queued_cleanups;

  std::mutex m_mu;
  size_t m_total_reads = 0;

  // Current import states for all bam files.
  std::vector<std::unique_ptr<bam_file_state>> m_bam_files;

  bool m_got_paired = false;

  // Bam import batching
  std::vector<std::unique_ptr<bam_line_batch>> m_free_bam_line_batches;
  std::vector<std::unique_ptr<bam_pair_cache_t>> m_free_pair_caches;

  // If non-zero, the range of bases in each read that should be saved
  // to generate a trimmed dataset.
  unsigned m_cut_reads_start = 0;
  unsigned m_cut_reads_end = 0;

  // Contig remapping is not currently supported for seqsets, but this lets us
  // reuse code from anchored.
  using contig_remap_t = std::map<int, int>;
  contig_remap_t m_remap_contigs;
};

// Imports reads, using the given type T as a per-thread state object.
//
// T should be a subclass of parallel_local.
//
// T should have a typedef "init_type", which is the type of data
// supplied to initialize a new T.
//
// It should support the following methods:
//
// Process a batch of reads:
// void process(const std::vector<std::pair<read_id, unaligned_reads>>& reads):
//
// The destructor should finish any processing that needs to be done
// for previously provided reads.
//

template <typename T>
class read_importer : public read_importer_base {
 private:
  using init_type = typename T::init_type;

 public:
  read_importer(typename T::init_type init_data,
                progress_handler_t progress = null_progress_handler)
      : read_importer_base(progress), m_init_data(init_data) {}

  void process_read_batch(parallel_state& st, read_batch& batch) override {
    T* batch_state = st.get_local<T>(m_init_data);
    batch_state->process(batch.reads());
    batch.clear();
  }

 private:
  init_type m_init_data;
};

}  // namespace build_seqset
