#pragma once

#include "modules/variants/assemble.h"

namespace variants {

class pair_counter : public assemble_pipeline_interface {
 public:
  pair_counter(const assemble_options& options, pipeline_step_t output);
  pair_counter() = delete;
  ~pair_counter() { flush(); }

  void on_assembly(assembly_ptr a) override;

 private:
  void flush();
  void advance_active();
  void advance_left();
  void advance_right();

  void calc_score(assembly& a) const;

  // Returns true if some advancement was done.
  bool advance_some(aoffset_t target_offset);

  std::vector<uint32_t> find_pair_matches(
      const assembly& a, const std::unordered_multiset<uint32_t, unsalted_hash>& read_ids,
      bool forward) const;

  assemble_options m_options;

  aoffset_t m_cur_offset = 0;
  std::multimap<aoffset_t /* offset of right anchor */, read_id_set /* read ids */>
      m_left;
  std::unordered_multiset<uint32_t, unsalted_hash> m_left_read_ids;

  std::multimap<aoffset_t /* offset of right anchor */, assembly_ptr> m_active;

  std::multimap<aoffset_t /* offset of left anchor */, assembly_ptr> m_right;
  std::unordered_multiset<uint32_t, unsalted_hash> m_right_read_ids;
  pipeline_step_t m_output;
};

}  // namespace variants
