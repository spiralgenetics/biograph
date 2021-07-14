#include "python/biograph/variants/add_ref.h"

#include "modules/variants/add_ref.h"
#include "modules/variants/pair_cov.h"
#include "modules/variants/read_cov.h"
#include "modules/variants/scaffold.h"
#include "modules/variants/trace_ref.h"
#include "python/biograph/variants/assembly.h"
#include "python/biograph/variants/module.h"
#include "python/biograph/variants/pipeline.h"
#include "python/common.h"

using namespace pybind11;
using namespace variants;

class __attribute__((visibility("hidden"))) add_ref_generator : public asm_pipeline_wrapper {
 public:
  add_ref_generator(const reference_wrapper& ref, const std::string& scaffold_name, object input,
                    unsigned pad_bases, bool whole_ref, int max_len, bool rev_comp);
  void init();

  void on_input(assembly_ptr a) override;
  void on_input_done() override;

 private:
  assemble_options m_options;
  ::variants::scaffold m_scaffold;

  unsigned m_pad_bases;
  bool m_whole_ref;
  int m_max_len;
  bool m_rev_comp;

  std::unique_ptr<add_ref> m_add_ref;
};

add_ref_generator::add_ref_generator(const reference_wrapper& ref, const std::string& scaffold_name,
                                     object input, unsigned pad_bases, bool whole_ref, int max_len,
                                     bool rev_comp)
    : asm_pipeline_wrapper(input),
      m_pad_bases(pad_bases),
      m_whole_ref(whole_ref),
      m_max_len(max_len),
      m_rev_comp(rev_comp) {
  m_options.scaffold_name = scaffold_name;
  m_scaffold = trace_ref::ref_to_scaffold(ref.get_reference().get(), scaffold_name);
  if (m_rev_comp) {
    m_scaffold = m_scaffold.rev_comp();
  }
  m_options.scaffold = &m_scaffold;
}

void add_ref_generator::init() {
  m_add_ref.reset(
      new add_ref(m_options, m_pad_bases, m_whole_ref, m_max_len, make_pipeline_output()));
}

void add_ref_generator::on_input(assembly_ptr a) { m_add_ref->add(std::move(a)); }

void add_ref_generator::on_input_done() { m_add_ref.reset(); }

void bind_add_ref(module& m) {
  define_pipeline_generator<add_ref_generator, const reference_wrapper& /* ref */,
                            const std::string& /* scaffold_name */, object /* input */,
                            unsigned /* pad_bases*/, bool /* whole ref */, int /* max len */,
                            bool /* rev comp */>(
      m, "add_ref_assemblies", "AddRefGenerator", arg("ref"), arg("scaffold_name"), arg("input"),
      arg("pad_bases") = 0, arg("whole_ref") = false, arg("max_len") = 0, arg("rev_comp") = false,
      R"DOC(how does this work?)DOC");
}
