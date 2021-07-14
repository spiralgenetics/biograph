#pragma once

#include "modules/bio_base/reference.h"
#include "modules/bio_base/seq_position.h"
#include "modules/bio_base/seqset.h"
#include "modules/io/packed_vector.h"
#include "modules/io/progress.h"
#include "modules/io/spiral_file.h"

namespace variants {

struct ref_anchor {
  // Location of anchor
  seq_position pos;
  // If true, anchor faces left in the reference.
  bool rev_comp;

  friend std::ostream& operator<<(std::ostream& os, const ref_anchor& a) {
    return os << "[@ " << a.pos.scaffold_id << ":" << a.pos.position << (a.rev_comp ? " RC" : "")
              << "]";
  }
  std::string as_string() const {
    std::stringstream os;
    os << (*this);
    return os.str();
  }

  ref_anchor& operator+=(int64_t offset) {
    if (rev_comp) {
      offset = -offset;
    }
    if (offset < 0) {
      CHECK_GE(pos.position, -offset);
    }
    pos.position = int64_t(pos.position) + offset;
    return *this;
  }
  ref_anchor& operator++() { return (*this) += 1; }
  ref_anchor operator+(int64_t offset) const {
    ref_anchor result = *this;
    result += offset;
    return result;
  }
  ref_anchor& operator-=(int64_t offset) { return (*this) += (-offset); }
};

// ref_map tracks which entries in a seqset match reference, and how
// ambiguous those matches are.
class ref_map {
 private:
  static constexpr unsigned k_fwd_flag = 1 << 7;
  static constexpr unsigned k_rev_flag = 1 << 6;
  static constexpr unsigned k_count_mask = (1 << 6) - 1;
  // Minimum chunk size of reference to process at once.
  friend class ref_map_test;
  static constexpr size_t k_min_chunk_size = 25600;

 public:
  class entry {
   public:
    explicit entry(uint8_t val) : m_val(val) {
      if (fwd_match() && rev_match()) {
        DCHECK_GT(match_count(), 1);
      } else if (!fwd_match() && !rev_match()) {
        DCHECK_EQ(match_count(), 0);
      }
    }

    // Returns true if this entry matches reference at all.
    bool is_match() const { return fwd_match() || rev_match(); }
    // Returns true if this entry matches reference in the forward direction
    bool fwd_match() const { return m_val & k_fwd_flag; }
    // Returns true if this entry matches reference in the reverse complement
    // direction
    bool rev_match() const { return m_val & k_rev_flag; }
    // Number of times this entry matches reference (bounded by a maximum of
    // k_count_mask)
    unsigned match_count() const { return m_val & k_count_mask; }

    bool is_unique() const { return match_count() == 1; }

    friend std::ostream& operator<<(std::ostream& os, const entry& e) {
      os << "[";
      if (e.fwd_match()) {
        os << "fwd ";
      }
      if (e.rev_match()) {
        os << "rev ";
      }
      os << "count=" << e.match_count() << "]";
      return os;
    }

   private:
    uint8_t m_val;
  };

  ref_map(const seqset* the_seqset, const reference* ref);
  ref_map(const seqset* the_seqset, const reference* ref, const spiral_file_create_state& state);
  ref_map(const seqset* the_seqset, const reference* ref, const spiral_file_open_state& state);

  void build(progress_handler_t progress = null_progress_handler);

  entry get(uint64_t seqset_id) const;

  boost::optional<ref_anchor> get_unique_ref_anchor(uint64_t seqset_id) const;

  dna_slice get_ref_slice(const ref_anchor& anchor) const;

 private:
  // Number of buckets to use to flush data to the ref map bitmap.
  static constexpr size_t k_num_flush_buckets = 32;
  static constexpr size_t k_flush_bucket_size = 64 * 1024;

  void flush_updates(std::vector<uint64_t>& seqset_ids, bool is_rev_comp);

  const seqset* m_seqset = nullptr;
  const reference* m_ref = nullptr;

  std::shared_ptr<mutable_packed_vector<unsigned, 8>> m_mutable_ref_map;
  std::shared_ptr<packed_vector<unsigned, 8>> m_ref_map;
  std::array<std::mutex, k_num_flush_buckets> m_flush_bucket_mu;
  size_t m_seqset_entries_per_flush_bucket = 0;
};

}  // namespace variants
