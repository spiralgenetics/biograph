#pragma once

// A class to provide vargraph-like coverage on assemblies.
//
// Inputs are:
//   A list of assemblies
//

#include "absl/container/btree_map.h"
#include "modules/variants/apply_edges.h"
#include "modules/variants/assemble.h"

namespace variants {

enum class cigar_op : char {
  MATCH = 'M',     // ref bases correspond to var bases
  DELETE = 'D',    // ref bases without var bases
  INSERT = 'I',    // var bases without ref bases
  REF_SKIP = 'N',  // Used for anchoring inserts.  Effectively the same as DELETE.
  PAD = 'P',       // Used for display purposes to align insert coverage
};

struct aligned_read {
  aoffset_t left_offset = std::numeric_limits<aoffset_t>::max();
  aoffset_t right_offset = std::numeric_limits<aoffset_t>::min();

  std::string cigar;
  dna_sequence seq;

  friend std::ostream& operator<<(std::ostream& os, const aligned_read& read);
};

class align_reads : public apply_edges_step {
 public:
  using on_aligned_func_t = std::function<void(const read_id_set& read_ids, aligned_read)>;
  align_reads(const on_aligned_func_t& on_aligned, bool refskip_anchor, pipeline_step_t output);
  ~align_reads();
  void on_assembly_edges(optional_aoffset reference_pos, const std::vector<assembly_ptr>& left_edges,
                         const std::vector<assembly_ptr>& inserts,
                         const std::vector<assembly_ptr>& right_edges) override;

 private:
  struct read_trace_key {
    read_id_set read_ids;
    int read_len_left = 0;

    bool operator<(const read_trace_key& rhs) const {
      int len_diff = rhs.read_len_left - read_len_left;
      if (len_diff) {
        return len_diff < 0;
      }

      return read_ids.total_order_lt(rhs.read_ids);
    }
    friend std::ostream& operator<<(std::ostream& os, const read_trace_key& key) {
      return os << "ReadTraceKey(left=" << key.read_len_left << ", ids=" << key.read_ids << ")";
    }
  };
  struct read_trace : public aligned_read {
    int read_len_tot = 0;
    cigar_op cur_cigar_op = cigar_op::MATCH;
    int cur_cigar_count = 0;

    friend std::ostream& operator<<(std::ostream& os, const read_trace& read) {
      os << static_cast<const aligned_read&>(read);
      if (read.cur_cigar_count) {
        os << ", cur op='" << char(read.cur_cigar_op) << "', count=" << read.cur_cigar_count << ".";
      }
      return os;
    }
  };
  using propagate_t = absl::btree_map<read_trace_key, read_trace>;

  void find_starts(const assembly_ptr& a, propagate_t* prop_out);
  void propagate(const assembly_ptr& a, const propagate_t& prop_in, propagate_t* prop_out);

  static void add_cigar(read_trace& read, cigar_op op, int num_bases);
  static void flush_cigar(read_trace& read);
  void output_aligned(read_trace_key key, read_trace read);
  void propagate_ref(const assembly_ptr& a, aoffset_t offset, read_trace_key key, read_trace read,
                     propagate_t* prop_out);
  void propagate_var(const assembly_ptr& a, aoffset_t offset, read_trace_key key, read_trace read,
                     propagate_t* prop_out);

  absl::btree_map<aoffset_t /* right offset */, propagate_t> m_propagate;
  on_aligned_func_t m_on_aligned;
  bool m_refskip_anchor;
};

}  // namespace variants
