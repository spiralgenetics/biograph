#pragma once

#include "modules/bio_base/dna_sequence.h"
#include "modules/bio_base/kmer.h"
#include "modules/build_seqset/repo_seq.h"

namespace build_seqset {

class part_repo {
 public:
  part_repo(unsigned partition_depth, const std::string& ref_path_prefix,
            const std::string& repo_path);

  struct partition_ref {
    kmer_t part_id;

    // Prefix of entries in this partition, e.g. "ACTG"
    dna_sequence prefix;

    // The repository containing all entries in this partition,
    // e.g. for "ACTG".
    std::shared_ptr<const seq_repository> main;

    // Repositories containing all entries for front pushes of entries
    // in this partition.  E.g. "AACT", "CACT", "GACT", "TACT".
    dna_base_array<
        std::pair<seq_repository::iterator, seq_repository::iterator>>
        pushed;

    dna_base_array<std::shared_ptr<const seq_repository>> pushed_repositories;

    // First sequence of next partition.
    dna_sequence next_entry;

    // Frees up references to this partition.
    void reset() {
      main.reset();
      for (dna_base b : dna_bases()) {
        pushed_repositories[b].reset();
      }
    }
  };

  class write_entry_buffer;
  void open_write_pass(const std::string& pass_name);
  std::unique_ptr<seq_repository::ref_builder> open_ref_builder(
      kmer_t part_id, const std::string& pass_name);

  void add_initial_repo(const dna_slice& reference_data);
  void write(const dna_slice& seq, unsigned fwd_suffixes, unsigned rc_suffixes);
  void write_using_repo(const dna_slice& seq, unsigned fwd_suffixes, unsigned rc_suffixes,
                        size_t repo_offset);
  void write(const seq_repository::entry_base& e);

  size_t write_with_expansions(const seq_repository::entry_base& e,
                               unsigned stride, unsigned count);
  void flush();

  std::unique_ptr<part_counts> release_part_counts(const std::string& pass_name);
  void reset_part_counts(const std::string& pass_name, std::unique_ptr<part_counts> counts);

  void for_each_partition(const std::string& pass_name,
                          const std::function<void(const partition_ref&)>& f,
                          progress_handler_t progress = null_progress_handler);
  // If "do_expand" is true, also populates pushed_repositories.
  std::vector<partition_ref> partitions(const std::string& pass_name,
                                        bool do_pushed = true,
                                        bool delete_on_close = false) const;

  const dna_slice& repo_slice() const { return m_repo_slice; }

  size_t partition_count() const { return 1 << (2 * m_depth); }
  size_t partition_depth() const { return m_depth; }

 private:
  // Keep statistics on 4^k_part_counts_depth chunks for each partition
  static constexpr unsigned k_part_counts_depth = 3;

  void dump_part_counts_if_needed() const;
  dna_sequence prefix_for_partition(kmer_t part_num) const;
  kmer_t partition_for_sequence(const dna_slice& seq) const;
  std::string ref_filename(kmer_t part_num, const std::string& pass_name) const;
  seq_repository::repo_builder* get_repo_builder();
  std::shared_ptr<seq_repository> open_part_repo(
      kmer_t part_num, const std::string& pass_name) const;
  void write_raw(const seq_repository::entry_data& e);
  std::tuple<std::shared_ptr<const seq_repository>, seq_repository::iterator,
             seq_repository::iterator>
  range_including_prefix(const std::vector<partition_ref>& parts,
                         const dna_slice& seq) const;

  unsigned m_depth;

  std::string m_ref_prefix;
  std::string m_repo_path;

  membuf m_repo;
  dna_slice m_repo_slice;

  std::mutex m_mu;
  std::vector<std::unique_ptr<seq_repository::ref_builder>> m_ref_builders;
  std::unique_ptr<seq_repository::repo_builder> m_repo_builder;
  std::unique_ptr<part_counts> m_part_counts;
  std::string m_part_counts_pass_name;
};

}  // namespace build_seqset
