#include "python/biograph/variants/trim_ref.h"

#include "modules/variants/scaffold.h"
#include "modules/variants/trace_ref.h"
#include "modules/variants/trim_ref.h"
#include "python/biograph/variants/pipeline.h"
#include "python/common.h"

using namespace pybind11;
using namespace variants;

class __attribute__((visibility("hidden"))) trim_ref_generator : public asm_pipeline_wrapper {
 public:
  trim_ref_generator(const reference_wrapper& ref, const std::string& scaffold_name, object input,
                     bool rev_comp);

  void init();
  void on_input(assembly_ptr a) override;
  void on_input_done() override;

 private:
  assemble_options m_options;
  ::variants::scaffold m_scaffold;

  std::unique_ptr<ref_trimmer> m_trim_ref;
};

trim_ref_generator::trim_ref_generator(const reference_wrapper& ref,
                                       const std::string& scaffold_name, object input,
                                       bool rev_comp)
    : asm_pipeline_wrapper(input) {
  m_options.scaffold_name = scaffold_name;
  m_scaffold = trace_ref::ref_to_scaffold(ref.get_reference().get(), scaffold_name);
  if (rev_comp) {
    m_scaffold = m_scaffold.rev_comp();
  }
  m_options.scaffold = &m_scaffold;
  m_discard_reference_only = true;
}

void trim_ref_generator::init() {
  m_trim_ref.reset(new ref_trimmer(m_options, make_pipeline_output()));
}

void trim_ref_generator::on_input(assembly_ptr a) { m_trim_ref->add(std::move(a)); }

void trim_ref_generator::on_input_done() { m_trim_ref.reset(); }

void bind_trim_ref(module& m) {
  define_pipeline_generator<trim_ref_generator, const reference_wrapper&,
                            const std::string& /* scaffold_name */, object /* input */,
                            bool /* rev_comp */>(
      m, "trim_ref", "TrimRefGenreator", arg("ref"), arg("scaffold_name"), arg("input"),
      arg("rev_comp") = false, R"DOC(how does this work?)DOC");
}
