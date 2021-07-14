#pragma once

#include "modules/variants/assemble.h"

#include <boost/icl/interval_map.hpp>

namespace variants {

// Filters assemblies to the number of allowed alleles based on a user supplied sort function.
class limit_alleles : public sorted_output_pipeline_step {
 public:
  // Sort function; this should sort the list of supplied assemblies, and return the highest
  // priority assemblies first.
  using sort_func_t = std::function<std::vector<assembly_ptr>(std::vector<assembly_ptr>)>;

  // Any assemblies that exceed the allele limit will be passed to this function to mark them as
  // over the limit.
  using on_limited_func_t = std::function<void(const assembly_ptr&)>;

  limit_alleles(size_t max_alleles, const sort_func_t& sort_func,
                const on_limited_func_t& on_limit_func, pipeline_step_t output)
      : sorted_output_pipeline_step(std::move(output)),
        m_max_alleles(max_alleles),
        m_sort_func(sort_func),
        m_on_limited_func(on_limit_func) {
    set_expected_order(assembly::left_offset_less_than);
  }
  limit_alleles() = delete;
  ~limit_alleles();

  void on_assembly(assembly_ptr a) override;

 private:
  using interval_t = typename boost::icl::right_open_interval<aoffset_t>;
  using depths_t =
      boost::icl::interval_map<aoffset_t, size_t, boost::icl::partial_absorber, std::less,
                               boost::icl::inplace_plus, boost::icl::inter_section, interval_t>;
  void advance_to(aoffset_t target);
  void advance_towards(aoffset_t target);
  void flush_block_contents();
  void sort_and_limit_block_contents();
  bool is_exceeded(const depths_t& depths) const;

  static interval_t interval_for_assembly(const assembly_ptr& a);

  std::multimap<aoffset_t /* right offset */, assembly_ptr> m_active;
  std::vector<assembly_ptr> m_block_contents;

  aoffset_t m_cur_offset = 0;

  unsigned m_max_alleles;
  sort_func_t m_sort_func;
  on_limited_func_t m_on_limited_func;
};

}  // namespace variants
