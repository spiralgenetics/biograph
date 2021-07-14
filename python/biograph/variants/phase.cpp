#include "python/biograph/variants/phase.h"

#include "modules/variants/phase.h"
#include "modules/variants/scaffold.h"
#include "python/biograph/variants/assembly.h"
#include "python/biograph/variants/module.h"
#include "python/biograph/variants/pipeline.h"
#include "python/common.h"

#include <pybind11/functional.h>

using namespace pybind11;
using namespace variants;

class __attribute__((visibility("hidden"))) join_phases_generator : public asm_pipeline_wrapper {
 public:
  join_phases_generator(object input, size_t max_phase_len, size_t max_phase_asm_len);
  void init();

  void on_input(assembly_ptr a) override;
  void on_input_done() override;

 private:
  size_t m_max_phase_len = 0;
  size_t m_max_phase_asm_len = 0;

  std::unique_ptr<join_phases> m_join_phases;
};

join_phases_generator::join_phases_generator(object input, size_t max_phase_len,
                                             size_t max_phase_asm_len)
    : asm_pipeline_wrapper(input),
      m_max_phase_len(max_phase_len),
      m_max_phase_asm_len(max_phase_asm_len) {}

void join_phases_generator::init() {
  m_join_phases.reset(
      new join_phases(m_max_phase_len, m_max_phase_asm_len, make_pipeline_output()));
}

void join_phases_generator::on_input(assembly_ptr a) { m_join_phases->add(std::move(a)); }

void join_phases_generator::on_input_done() { m_join_phases.reset(); }

class __attribute__((visibility("hidden"))) split_phases_generator : public asm_pipeline_wrapper {
 public:
  split_phases_generator(object input);
  void init();

  void on_input(assembly_ptr a) override;
  void on_input_done() override;

 private:
  std::unique_ptr<split_phases> m_split_phases;
};

split_phases_generator::split_phases_generator(object input) : asm_pipeline_wrapper(input) {}

void split_phases_generator::init() {
  m_split_phases.reset(new split_phases(make_pipeline_output()));
}

void split_phases_generator::on_input(assembly_ptr a) { m_split_phases->add(std::move(a)); }

void split_phases_generator::on_input_done() { m_split_phases.reset(); }

class __attribute__((visibility("hidden"))) resolve_phase_conflicts_generator
    : public asm_pipeline_wrapper {
 public:
  using resolve_conflict_func_t = resolve_phase_conflicts::resolve_conflict_func_t;
  resolve_phase_conflicts_generator(const resolve_conflict_func_t& resolve_func, object input);
  void init();

  void on_input(assembly_ptr a) override;
  void on_input_done() override;

 private:
  std::unique_ptr<resolve_phase_conflicts> m_resolve_phase_conflicts;
  resolve_conflict_func_t m_resolve_func;
};

resolve_phase_conflicts_generator::resolve_phase_conflicts_generator(
    const resolve_conflict_func_t& resolve_func, object input)
    : asm_pipeline_wrapper(input), m_resolve_func(resolve_func) {}

void resolve_phase_conflicts_generator::init() {
  m_resolve_phase_conflicts.reset(
      new resolve_phase_conflicts(m_resolve_func, make_pipeline_output()));
}

void resolve_phase_conflicts_generator::on_input(assembly_ptr a) {
  m_resolve_phase_conflicts->add(std::move(a));
}

void resolve_phase_conflicts_generator::on_input_done() { m_resolve_phase_conflicts.reset(); }

phase_set phase_set_from_format_fields(int format_index, iterable formats) {
  size_t sample_num = 0;
  phase_set ids;
  for (const auto& format_elem : formats) {
    PyObject* py_format_elem = format_elem.ptr();
    const char* utf8 = PyUnicode_AsUTF8(py_format_elem);
    if (!utf8) {
      throw error_already_set();
    }

    size_t index_left = format_index;
    const char* field_start = utf8;
    while (index_left) {
      field_start = strchr(field_start, ':');
      --index_left;
      if (!field_start) {
        throw(std::runtime_error(
            printstring("Unable to find format field #%d of '%s'", format_index, utf8)));
      }
      ++field_start;
    }

    const char* field_end = field_start;
    while (*field_end && *field_end != ':') {
      ++field_end;
    }

    if (field_start == field_end || (*field_start == '.' && field_start + 1 == field_end)) {
      // Empty format field
    } else {
      std::string id = std::to_string(sample_num);
      id += ":";
      id.insert(id.end(), field_start, field_end);
      ids.emplace(std::move(id));
    }

    ++sample_num;
  }
  return ids;
}

void bind_phases(module& m) {
  define_pipeline_generator<join_phases_generator, object /* input */, size_t /* max_phase_len */,
                            size_t /* max_phase_asm_len */>(
      m, "join_phases", "JoinPhasesGenerator", arg("input"), arg("max_phase_len") = 1000,
      arg("max_phase_asm_len") = 1000, R"DOC(how does this work?)DOC");
  m.def("propagate_subassembly_coverage", [](assembly_ptr a) {
    propagate_subassembly_coverage(a.get());
    return a;
  });
  define_pipeline_generator<split_phases_generator, object /* input */>(
      m, "split_phases", "SplitPhasesGenerator", arg("input"), R"DOC(how does this work?)DOC");
  define_pipeline_generator<resolve_phase_conflicts_generator,
                            resolve_phase_conflicts_generator::resolve_conflict_func_t,
                            object /* input */>(
      m, "resolve_phase_conflicts", "ResolvePhaseConflictsGenerator", arg("on_conflict"),
      arg("input"), R"DOC(how does this work?)DOC");

  class_<phase_set>(m, "PhaseSet", "Contains a set of phase ids")
      .def(init<>())
      .def(init([](iterable elems) -> phase_set {
        phase_set result;
        for (const auto& elem : elems) {
          result.insert(cast<std::string>(elem));
        }
        return result;
      }))
      .def(
          "__iter__",
          [](const phase_set& ids) {
            return make_iterator<return_value_policy::copy>(ids.begin(), ids.end());
          },
          keep_alive<0, 1>())
      .def("add", [](phase_set& ids, std::string id) { ids.insert(std::move(id)); })
      .def("__contains__",
           [](const phase_set& ids, const std::string& id) -> bool { return ids.count(id); })
      .def("__len__", &phase_set::size)
      .def("clear", &phase_set::clear)
      .def("difference", &phase_set::operator-)
      .def("union", &phase_set::operator+)
      .def("intersection", &phase_set::operator&)
      .def("__repr__",
           [](const phase_set& ids) {
             std::stringstream out;
             out << "PhaseSet([";
             bool first = true;
             for (const auto& id : ids) {
               if (first) {
                 first = false;
               } else {
                 out << ",";
               }
               out << "\"" << id << "\"";
             }
             out << "])";
             return out.str();
           })
      .def("__str__", str_from_ostream<phase_set>)
      .def_static("from_format_fields", &phase_set_from_format_fields, arg("field_index"),
                  arg("fields"),
                  "Builds a PhaseSet from VCF format fields.  field_index is the field number, and "
                  "fields a list of the format fields");
}
