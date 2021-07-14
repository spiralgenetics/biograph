#include "python/biograph/variants/pipeline.h"

#include "python/biograph/variants/assembly.h"
#include "python/common.h"

using namespace variants;
using namespace pybind11;

pipeline_step_t asm_pipeline_wrapper::make_pipeline_output() {
  // Make sure this object doesn't get deleted until we're done processing pipeline entries.
  std::shared_ptr<asm_pipeline_wrapper> shared_ref = shared_from_this();
  return make_unique<assemble_lambda_output>(
      [this, shared_ref](assembly_ptr a) {
        if (m_discard_reference_only && a->matches_reference) {
          return;
        }
        m_output_queue.emplace_back(std::move(a));
      },
      "python_pipeline_output");
}

assembly_ptr asm_pipeline_wrapper::next() {
  while (!m_input_done && m_output_queue.empty()) {
    read_more_input();
  }
  if (m_input_done && m_output_queue.empty()) {
    PyErr_SetObject(PyExc_StopIteration, Py_None);
    throw error_already_set();
  }
  CHECK(!m_output_queue.empty());

  assembly_ptr a = std::move(m_output_queue.front());
  m_output_queue.pop_front();

  return a;
}

void asm_pipeline_wrapper::read_more_input() {
  CHECK(!m_input_done);

  PyObject* py_next = PyIter_Next(m_input_iter_obj.ptr());
  if (PyErr_Occurred()) throw pybind11::error_already_set();
  if (py_next) {
    object next_obj = reinterpret_steal<object>(py_next);
    {
      assembly_ptr a = cast<assembly_ptr>(next_obj);
      if (!a) {
        throw(io_exception("Assemblies must not be blank"));
      }

      gil_scoped_release allow_threads;
      check_assembly_from_user(*a);
      aoffset_t min_left = min(a->left_offset, a->right_offset);
      if (min_left < m_last_left_offset) {
        throw(io_exception(printstring(  //
            "Assemblies must be sorted in order; "
            "got assembly with offset %d after assembly with offset %d: %s",
            aoffset_t(a->left_offset), m_last_left_offset, dump_assembly_and_vars(*a).c_str())));
      }
      m_last_left_offset = min_left;
      on_input(std::move(a));
    }
  } else {
    m_input_done = true;
    {
      gil_scoped_release allow_threads;
      on_input_done();
    }
  }
}

object just_return_self(object self) { return self; }
