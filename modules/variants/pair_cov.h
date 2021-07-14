#include "modules/variants/assemble.h"

#include <boost/container/flat_map.hpp>
#include <boost/container/flat_set.hpp>
#include <boost/container/small_vector.hpp>
#include <boost/icl/interval.hpp>
#include <boost/icl/interval_set.hpp>
#include "absl/container/btree_map.h"
#include "absl/container/btree_set.h"
#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"
#include "modules/bio_base/readmap.h"
#include "modules/io/ref_count.h"

namespace variants {

class pair_cov : public sorted_output_pipeline_step {
 public:
  pair_cov(const assemble_options& opts, pipeline_step_t output);
  ~pair_cov();

  void on_assembly(assembly_ptr a) override;

 private:
  // Simple interval set that just tracks the hull around all the
  // intervals that have been added.  This is less precise than using a
  // full boost::icl::interval_set, but is a whole lot faster.
  struct simple_interval_set {
   public:
    constexpr simple_interval_set() = default;
    constexpr simple_interval_set(aoffset_t lower, aoffset_t upper)
        : m_lower(lower), m_upper(upper) {}

    aoffset_t lower() const { return m_lower; }
    aoffset_t upper() const { return m_upper; }

    simple_interval_set& operator=(const simple_interval_set& rhs) {
      m_lower = rhs.m_lower;
      m_upper = rhs.m_upper;
      return *this;
    }

    simple_interval_set& operator+=(const simple_interval_set& rhs) {
      if (is_empty(*this)) {
        (*this) = rhs;
      } else {
        m_lower = std::min<aoffset_t>(m_lower, rhs.m_lower);
        m_upper = std::max<aoffset_t>(m_upper, rhs.m_upper);
      }
      return *this;
    }

    friend std::ostream& operator<<(std::ostream& os, const simple_interval_set& s) {
      return os << "[" << s.m_lower << "," << s.m_upper << ")";
    }

    bool operator==(const simple_interval_set& rhs) const {
      return m_lower == rhs.m_lower && m_upper == rhs.m_upper;
    }

    // ADL resolutions so this looks like a boost icl interval set to callers.
    friend bool contains(const simple_interval_set& s, aoffset_t pos) {
      return pos >= s.m_lower && pos < s.m_upper;
    }
    // Returns true if rhs is entirely contained in lhs.
    friend bool contains(const simple_interval_set& lhs, const simple_interval_set& rhs) {
      if (is_empty(lhs)) {
        return false;
      }
      if (is_empty(rhs)) {
        return false;
      }

      if (rhs.lower() < lhs.lower()) {
        return false;
      }
      if (rhs.upper() > lhs.upper()) {
        return false;
      }
      return true;
    }
    friend bool is_empty(const simple_interval_set& s) { return s.m_upper <= s.m_lower; }

    // Returns true if any part of lhs overlaps rhs.
    friend bool overlaps(const simple_interval_set& lhs, const simple_interval_set& rhs) {
      CHECK(!is_empty(lhs));
      CHECK(!is_empty(rhs));

      if (rhs.lower() >= lhs.upper()) {
        return false;
      }

      if (lhs.lower() >= rhs.upper()) {
        return false;
      }

      return true;
    }

   private:
    aoffset_t m_lower = 0;
    aoffset_t m_upper = 0;
  };

  using interval_t = simple_interval_set;
  using interval_set_t = simple_interval_set;

  struct result;
  using result_ptr =
      explicit_shared_ptr<result, false /* not atomic */, false /* no implicit copy */,
                          false /* no implicit destruct */>;

  struct result_offset {
    result* r = nullptr;

    // Read start, relative to the beginning of the assembly.
    aoffset_t read_start;
  } __attribute__((packed));

  // Only has to be unique among seqset ranges.  But the read id range
  // must be small enough that we can represent its bitmask in a
  // read_id_mask_t.
  struct multi_mid {
    uint32_t multi_id;
    int size;
    uint32_t read_id_chunk;

    bool operator==(const multi_mid& rhs) const noexcept {
      return multi_id == rhs.multi_id && size == rhs.size && read_id_chunk == rhs.read_id_chunk;
    }
  };

  struct multi_mid_hash {
    size_t operator()(const multi_mid& mm) const {
      size_t hash_val = boost::hash_value(mm.multi_id);
      boost::hash_combine(hash_val, mm.size);
      boost::hash_combine(hash_val, mm.read_id_chunk);
      return hash_val;
    }
  };

  struct read_id_masks {
    constexpr read_id_masks(read_id_mask_t total_arg, read_id_mask_t pending_arg)
        : total(total_arg), pending(pending_arg) {}

    read_id_mask_t total;
    read_id_mask_t pending;

    bool operator>(const read_id_masks& rhs) const {
      if (total != rhs.total) {
        return total > rhs.total;
      }
      return pending > rhs.pending;
    }
  };

  using roffs_t = boost::container::small_vector<result_offset, 3>;
  using pte_results_t = boost::container::flat_map<
      read_id_masks, roffs_t, std::greater<read_id_masks>,
      boost::container::small_vector<std::pair<read_id_masks, roffs_t>, 3>>;

  struct pair_table_entry {
    // Allow moving:
    pair_table_entry() = default;
    pair_table_entry(pair_table_entry&&) = default;
    pair_table_entry& operator=(pair_table_entry&&) = default;

    // But not copying:
    pair_table_entry(const pair_table_entry&) = delete;
    pair_table_entry& operator=(const pair_table_entry&) = delete;

    // Cache of union of all adjust(roff.read_starts, r->sv_adjust)
    interval_set_t tot_read_starts;

    // Union of all pending_masks in results.  This may contain
    // bits that are not in any result_offset::pending_mask entry if
    // entries get expired while still having pending mates.
    read_id_mask_t pending_mask;

    // Any results that contain this read.
    pte_results_t results;
  };

  struct pair_table : boost::noncopyable {
    absl::flat_hash_map<multi_mid, pair_table_entry, multi_mid_hash> entries;
  };

  struct pair_cov_pg;

  static simple_interval_set adjust_interval_set(const simple_interval_set& orig, aoffset_t adjust);
  static simple_interval_set adjust_interval_set(const simple_interval_set& orig,
                                                 const simple_interval_set& adjust);
  static simple_interval_set unadjust_interval_set(const simple_interval_set& orig,
                                                   const simple_interval_set& adjust);
  static interval_t closed_interval(aoffset_t start, aoffset_t end);
  static interval_t closed_interval(aoffset_t single_offset);
  static interval_set_t closed_interval_set(aoffset_t single_offset);
  static aoffset_t interval_set_upper_bound(const interval_set_t& s);

  // Add an assembly, save its reads, and advance.
  void add_assembly(pair_cov_pg* pg, assembly_ptr a);

  void advance_ref_to(aoffset_t pos);
  void advance_ref_towards(aoffset_t pos);
  void flush_active_to_here();
  void flush_pair_supported_offsets(result* r);
  void reap_result(std::unique_ptr<result> result);
  void save_reads(pair_cov_pg* pg);
  void match_mates(pair_cov_pg* pg);

  void join(std::vector<std::unique_ptr<pair_cov_pg>> pgs);
  void join(std::unique_ptr<pair_cov_pg> pg);
  void report_progress();
  void propagate_and_fill(pair_table& old_table, pair_table* new_table);
  std::string dump(pair_cov_pg* pg) const;
  std::string dump_read_id_mask(read_id_mask_t mask) const;
  std::string dump_pair_table(const pair_table& pt) const;

  void save_read(pair_cov_pg* pg, uint32_t read_id, const readmap::read& rd, const multi_mid& mm,
                 aoffset_t read_start);
  void save_pending_reads(pair_cov_pg* pg);
  bool find_mate(pair_cov_pg* pg, uint32_t mate_id, int mate_len, pair_table_entry& pte,
                 aoffset_t end_of_read_offset);
  absl::flat_hash_set<result*> find_expired_results() const;
  void flush_old();
  void update_adjusts();

  // Distance after which we can expire things.
  aoffset_t max_care_distance() const;

  std::string dump_multi_mid(const multi_mid& mm) const;

  assemble_options m_opts;

  std::vector<assembly_ptr> m_cur_inserts;
  std::vector<assembly_ptr> m_cur_non_inserts;

  aoffset_t m_cur_ref_offset = 0;
  aoffset_t m_next_flush_old = 0;

  // Path groups to join, and when they should be joined in to m_ref_pg.
  std::multimap<aoffset_t /* right offset */,
                std::unique_ptr<pair_cov_pg> /* path group to join at the given position */>
      m_active;

  pair_table m_main_pair_table;

  // Results which have rejoined but still need their pairing data kept around.
  absl::btree_set<result_ptr, result_ptr::less_than> m_pending_results;

  // Future adjustments
  absl::btree_map<
      aoffset_t /* reference offset to do adjustment? */,
      absl::flat_hash_map<result*, std::pair<result_ptr, interval_set_t /* adjustment */>>>
      m_future_adjusts;

  time_t m_last_report = 0;
  aoffset_t m_last_report_offset = 0;
  size_t m_last_report_asms = 0;
  size_t m_cur_asm_count = 0;
};

}  // namespace variants
