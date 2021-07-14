#pragma once

#include "modules/bio_base/reference.h"
#include "modules/bio_base/seqset.h"
#include "modules/io/packed_vector.h"
#include "modules/io/progress.h"
#include "modules/io/spiral_file.h"

// A refmatch allows looking up seqset entries and seeing how many
// different sections of reference they match, and in what directions.
//
// This allows one to easily see whether a seqset entry matches
// reference, and if so, whether it's ambiguous, and which
// orientations of it need to be looked up to find where it matches.
class refmatch {
 public:
  static constexpr unsigned k_fwd_flag = 1 << 7;
  static constexpr unsigned k_rev_flag = 1 << 6;
  static constexpr unsigned k_count_mask = (1 << 6) - 1;
  static const product_version k_refmatch_version;

  class entry {
   public:
    entry() = delete;
    entry(bool has_fwd, bool has_rev, unsigned tot_count)
        : m_has_fwd(has_fwd), m_has_rev(has_rev), m_count(tot_count) {
      if (has_fwd && has_rev) {
        CHECK_GT(tot_count, 1);
      } else if (has_fwd || has_rev) {
        CHECK_GE(tot_count, 1);
      } else {
        CHECK_EQ(tot_count, 0);
      }
    }
    entry(const entry&) = default;

    // True if this entry matches any section of reference in a forward
    // direction.
    bool has_fwd() const { return m_has_fwd; }
    // True if this entry matches any section of reference in a reverse
    // direction.
    bool has_rev() const { return m_has_rev; }
    // Total number of sections of reference this entry matches.
    unsigned matches() const { return m_count; }

    friend std::ostream& operator<<(std::ostream& os, const entry& e) {
      os << "[";
      if (e.has_fwd()) {
        os << "fwd ";
      }
      if (e.has_rev()) {
        os << "rev ";
      }
      os << "count=" << e.matches() << "]";
      return os;
    }

   private:
    bool m_has_fwd;
    bool m_has_rev;
    unsigned m_count;
  };

  refmatch(const seqset* the_seqset, const reference* ref,
           const spiral_file_open_state& state);

  entry get(uint64_t seqset_id) const;

 protected:
  refmatch(const seqset* the_seqset, const reference* ref);

  const seqset* m_seqset = nullptr;
  const reference* m_ref = nullptr;

  boost::optional<packed_vector<unsigned, 8>> m_per_entry;

  std::unordered_map<uint64_t /* seqset id */, unsigned /* count */> m_overflow;
};

class refmatch_builder : private refmatch {
 public:
  // Changable for testing
  static size_t g_min_chunk_size;

  refmatch_builder(const seqset* the_seqset, const reference* ref);

  void build(const spiral_file_create_state& state,
             progress_handler_t progress = null_progress_handler);

 private:
  void walk_reference(progress_handler_t progress);
  void save_overflow(const spiral_file_create_state& overflow_ids);

  boost::optional<mutable_packed_vector<unsigned, 8>> m_mutable_per_entry;
};
