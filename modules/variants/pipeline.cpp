#include "modules/variants/pipeline.h"

#include "modules/bio_base/readmap.h"
#include "modules/variants/align.h"
#include "modules/variants/calc_coverage.h"
#include "modules/variants/dedup.h"
#include "modules/variants/normalize.h"
#include "modules/variants/pair_counter.h"
#include "modules/variants/ploid_limit.h"
#include "modules/variants/rvg_exclude.h"
#include "modules/variants/scaffold.h"
#include "modules/variants/simple_genotype_filter.h"
#include "modules/variants/sort.h"
#include "modules/variants/trim_ref.h"

#include <boost/range/adaptor/reversed.hpp>

namespace variants {

namespace {

class report_discovered_assemblies : public assemble_pipeline_interface {
 public:
  report_discovered_assemblies(const assemble_options& opts, pipeline_step_t output)
      : m_output(std::move(output)), m_options(opts) {}
  report_discovered_assemblies() = delete;

  void on_assembly(assembly_ptr a) override {
    m_options.report_discovered_assemblies_func(m_options, *a);
    m_output->add(std::move(a));
  }
  std::string description() const override { return "report_discovered_assemblies"; }

 private:
  pipeline_step_t m_output;
  assemble_options m_options;
};

}  // namespace

assemble_pipeline::~assemble_pipeline() {}

assemble_pipeline::assemble_pipeline(const assemble_options& options, pipeline_step_t output)
    : m_options(options) {
  if (m_options.use_pop_tracer || m_options.pop_trace_anchor_drop) {
    CHECK(options.readmap);
    options.readmap->calc_read_len_limits_if_needed();
  }

  // This is like left_offset_less_than, but canonicalizes the order as much as possible.
  add_step<sorter>(canon_assembly_order());
  if (m_options.report_discovered_assemblies_func) {
    add_step<report_discovered_assemblies>(m_options);
  }
  add_step<ref_trimmer>(m_options);
  add_step<deduper>();

  m_output = std::move(output);
}

pipeline_step_t assemble_pipeline::make_parallel_input() {
  {
    std::lock_guard<std::mutex> l(m_mu);

    if (!m_first_ser_step) {
      CHECK(m_output);
      m_first_ser_step = make_pipeline(m_ser, std::move(m_output));
    } else {
      CHECK(!m_output);
    }
  }

  return make_pipeline(m_par, make_unique<assemble_lambda_output>(
                                  [this](assembly_ptr a) {
                                    std::lock_guard<std::mutex> l(m_mu);
                                    m_first_ser_step->add(std::move(a));
                                  },
                                  "pipeline_parallel_guard"));
}

pipeline_step_t assemble_pipeline::make_pipeline(
    const std::vector<std::function<pipeline_step_t(pipeline_step_t)>>& steps,
    pipeline_step_t output) {
  pipeline_step_t cur = std::move(output);
  for (auto f : boost::adaptors::reverse(steps)) {
    cur = f(std::move(cur));
  }
  return cur;
}

void assemble_pipeline::add_standard_variants_pipeline() {
  CHECK(m_options.scaffold);
  add_step<aligner>(m_options);
  add_step<align_splitter>();
  add_step<normalizer>(m_options);
  add_step<exact_deduper>();
  add_step<vcf_padder>(m_options);
  if (!m_options.use_pop_tracer) {
    add_step<calc_coverage>(m_options);
    add_step<simple_genotype_filter>(m_options);
  }
  if (m_options.rvg_exclude) {
    add_step<rvg_exclude>(m_options);
  }
}

}  // namespace variants
