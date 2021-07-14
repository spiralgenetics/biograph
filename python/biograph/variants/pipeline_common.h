#pragma once

#include <pybind11/pybind11.h>
#include <boost/bind/apply.hpp>

pybind11::object just_return_self(pybind11::object self);

template <typename GeneratorType, typename... Args>
pybind11::object pipeline_generator_maker(Args... args) {
  std::shared_ptr<GeneratorType> gen = std::make_shared<GeneratorType>(args...);
  gen->init();
  return pybind11::cast(gen);
}

template <typename GeneratorType, typename... Args, typename... ArgsAndDocs>
void define_pipeline_generator(pybind11::module& m, const char* make_generator_name,
                               const char* generator_type_name, ArgsAndDocs&&... args_and_docs) {
  m.def(make_generator_name, pipeline_generator_maker<GeneratorType, Args...>,
        std::forward<ArgsAndDocs>(args_and_docs)...);

  pybind11::class_<GeneratorType, std::shared_ptr<GeneratorType>>(m, generator_type_name)
      .def("__iter__", just_return_self)
      // Python3 needs:
      .def("__next__", &GeneratorType::next)
      // Python2 needs:
      .def("next", &GeneratorType::next);
}

// Makes a wrapper class for a pipeline step which takes input and
// produces output.  This class can be used with
// define_pipeline_generator.
template <typename BaseType /* asm_pipeline_wrapper or par_asm_pipeline_wrapper */,
          typename... Args>
class __attribute__((visibility("hidden"))) pipeline_wrapper_with_args : public BaseType {
 public:
  pipeline_wrapper_with_args(pybind11::object input, const Args&... args)
      : BaseType(input), m_args{args...} {}

  void init() { m_step = make_pipeline_step(); }

  virtual variants::pipeline_step_t make_pipeline_step() = 0;

  void on_input(variants::assembly_ptr a) override { m_step->add(std::move(a)); }
  void on_input_done() override {
    m_step->flush();
    m_step.reset();
  }

 protected:
  std::tuple<Args...> m_args;
  variants::pipeline_step_t m_step;
};
