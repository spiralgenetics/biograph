#pragma once

#include <boost/optional.hpp>
#include "absl/container/btree_map.h"
#include "modules/bio_base/reference.h"
#include "modules/bio_base/seq_position.h"
#include "modules/variants/assemble.h"
#include "modules/variants/ref_map.h"
#include "modules/variants/scaffold.h"

namespace variants {

using seqset_range_set = absl::btree_set<seqset_range>;
std::ostream& operator<<(std::ostream& os, const seqset_range_set& rs);

class graph_discover : public sorted_output_pipeline_step {
 public:
  graph_discover(const assemble_options& options, pipeline_step_t output);
  ~graph_discover() override;

  void on_assembly(assembly_ptr a) override;

  void flush() override;

  // Assembly that is within the window of discovery
  struct active_assembly {
    assembly_ptr a;

    virtual ~active_assembly() = default;
    virtual std::string to_string() const;
  };
  friend std::ostream& operator<<(std::ostream& os, const active_assembly& act);
  using active_assembly_ptr = std::unique_ptr<active_assembly>;
  virtual std::unique_ptr<active_assembly> make_active_assembly() {
    return make_unique<active_assembly>();
  }

  struct potential_anchor {
    // Assembly that this anchors to
    const active_assembly* act = nullptr;

    // Offset within the assembly.
    aoffset_t offset = 0;

    friend std::ostream& operator<<(std::ostream& os, const potential_anchor& anchor) {
      os << "anchor, offset=" << anchor.offset << " in ";
      if (anchor.act) {
        os << *anchor.act;
      } else {
        os << "null active ";
      }
      return os << " offset=" << anchor.offset;
    }
  };

  // Called to set up this active assembly after walking it.
  virtual void on_walk(active_assembly* a) {}

  // Called when we encounter the left side of the given assembly
  // during readahead, after adding it to our anchor lookups.
  virtual void on_readahead(const active_assembly* a) {}

  // Called when we encounter the left side of the given assembly during tracing,
  // after removing it from our anchor lookups.
  virtual void on_readahead_done(const active_assembly* a) {}

  // Called when we encounter the right side of the given assembly during tracing
  virtual void on_trace(const active_assembly* a) = 0;

  // Called when we advance to a certain point:
  virtual void on_advance_trace(aoffset_t offset) {}

  const assemble_options& opts() const { return m_options; }

 protected:
  assembly_ptr discover_extend_right(const active_assembly* act, aoffset_t offset, dna_slice seq,
                                     const std::string& tag,
                                     seqset_path new_rc_path = {}) WARN_UNUSED;
  assembly_ptr discover_anchor(const active_assembly* act, aoffset_t offset, dna_slice seq,
                               const potential_anchor& anchor, const std::string& tag,
                               seqset_path new_rc_path = {}) WARN_UNUSED;

 private:
  void advance_trace_to(aoffset_t offset);
  void advance_trace_towards(aoffset_t offset);

  void walk_readahead(active_assembly_ptr act);
  void process_readahead(const active_assembly* act);
  void process_readahead_done(const active_assembly* act);
  void process_trace(const active_assembly* act);

  void extend_assembly(const assembly* orig, assembly* a, dna_slice seq) const;

  absl::btree_multimap<aoffset_t /* min(right_offset, left_offset) */, active_assembly_ptr>
      m_readahead_done;

  absl::btree_multimap<aoffset_t /* max(right_offset, left_offset) */, active_assembly_ptr>
      m_trace_pending;

  assemble_options m_options;

  // Offset where we're tracing in m_trace_ready.
  aoffset_t m_trace_offset = std::numeric_limits<aoffset_t>::min();

  // Assembies that start at m_readahead_offset:
  std::vector<assembly_ptr> m_readahead_cur;
};

}  // namespace variants
