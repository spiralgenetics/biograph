#pragma once

#include "absl/container/btree_map.h"
#include "absl/container/btree_set.h"
#include "modules/io/ref_count.h"
#include "modules/variants/assemble.h"
#include "modules/variants/scaffold.h"

namespace variants {

// Joins together assemblies that share phase ids.
//
// Any assembly that has an empty phase_ids list or has matches_reference=true
// is outputted unchanged.
//
// Other assemblies are packaged up with assemblies with the same phase id and
// made as part of an encompassing assembly.
//
// max_phase_len specifies the maximum size of reference a phase
// may encompass in the middle, so we don't have to remember assemblies forever.
//
// max_phase_asm_len specifies the maximum length of the reference or
// variant sequence of an assembly before it will be considered
// separately from other things in its phase.

class join_phases : public sorted_output_pipeline_step {
 public:
  static const char k_join_phases_name[];

  join_phases(size_t max_phase_len, size_t max_phase_asm_len, pipeline_step_t output);
  ~join_phases() override;

  void on_assembly(assembly_ptr a) override;

  // For testing purposes:
  static bool g_check_invariants;

 private:
  struct active_t {
    // Constructed joined assembly that's the path traversed by this phased allele.
    assembly_ptr joined_a;

    // All reference assemblies seen since the last phased variant.
    std::vector<std::shared_ptr<assembly_ptr>> reference_after;

    aoffset_t right_offset = std::numeric_limits<aoffset_t>::max();
    aoffset_t var_right_offset = std::numeric_limits<aoffset_t>::max();
  };

  using active_ptr =
      explicit_shared_ptr<active_t, false /* not atomic */, false /* no implicit copy */,
                          false /* no implicit destruct */>;

  void check_invariants();
  void add_ref_asm(std::shared_ptr<assembly_ptr> shared_a);
  void add_to_active(active_t* act, std::shared_ptr<assembly_ptr> shared_a);
  void save_ref_asms(active_t* act);
  void split_active(active_t* act, const phase_set& keep_phases);
  active_ptr new_active(aoffset_t left_offset, const phase_set& new_phases);

  void add_var_asm(std::shared_ptr<assembly_ptr> shared_a);
  void flush();
  void advance_to(aoffset_t offset);
  void advance_towards(aoffset_t offset);
  void output_if_last_ref(std::shared_ptr<assembly_ptr> shared_a);
  void output_active(active_ptr act);

  void abort_phases(const phase_set& phases);
  void show_stats();

  static aoffset_t active_right_offset(const active_t& act);

  // In progress by phase id.
  absl::btree_map<std::string /* phase id */, active_ptr> m_active;

  // Phase ids to abort in the future
  absl::btree_map<aoffset_t, phase_set> m_abort_at;

  std::vector<assembly_ptr> m_cur_ref;
  aoffset_t m_cur_offset = 0;

  aoffset_t m_max_phase_len = 0;
  aoffset_t m_max_phase_asm_len = 0;
  size_t m_tot_seen = 0;
  size_t m_tot_seen_phases = 0;
  absl::btree_set<std::string> m_seen_phases;
};

void propagate_subassembly_coverage(assembly* a);

class split_phases : public sorted_output_pipeline_step {
 public:
  split_phases(pipeline_step_t output);
  ~split_phases() override;

  void on_assembly(assembly_ptr a) override;

 private:
  void advance_to(aoffset_t offset);

  std::map<aoffset_t, std::set<std::shared_ptr<assembly_ptr>>> m_active;
};

class resolve_phase_conflicts : public sorted_output_pipeline_step {
 public:
  using resolve_conflict_func_t =
      std::function<void(const assembly_ptr&, const assembly_ptr&, const phase_set& phase_ids)>;
  resolve_phase_conflicts(const resolve_conflict_func_t& resolve_conflict, pipeline_step_t output);
  ~resolve_phase_conflicts() override;

  void on_assembly(assembly_ptr a) override;

 private:
  void advance_to(aoffset_t offset);

  bool check_and_resolve_conflicts(const assembly_ptr& a, const assembly_ptr& b);

  resolve_conflict_func_t m_resolve_conflict;
  std::multimap<aoffset_t /* right offset */, assembly_ptr> m_active;
};

}  // namespace variants
