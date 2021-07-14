#include "python/biograph/variants/par_pipeline.h"

#include "python/biograph/variants/assembly.h"
#include "python/common.h"

using namespace variants;
using namespace pybind11;

constexpr size_t par_asm_pipeline_wrapper::k_max_queue_size;

pipeline_step_t par_asm_pipeline_wrapper::make_pipeline_output() {
  // Make sure this object doesn't get deleted until we're done processing pipeline entries.
  std::shared_ptr<par_asm_pipeline_wrapper> shared_ref = shared_from_this();
  return make_unique<assemble_lambda_output>(
      [this, shared_ref](assembly_ptr a) {
        if (m_discard_reference_only && a->matches_reference) {
          return;
        }
        std::lock_guard<std::mutex> l(m_mu);
        m_output_queue.emplace_back(std::move(a));
        m_output_cv.notify_one();
      },
      "python_pipeline_output");
}

assembly_ptr par_asm_pipeline_wrapper::next() {
  assembly_ptr a;

  {
    gil_scoped_release allow_threads;

    std::unique_lock<std::mutex> l(m_mu);
    if (!m_thread_started) {
      m_thread_started = true;
      m_pipeline_step_thread = std::async(std::launch::async, [this]() { pipeline_step_thread(); });
    }
    while (!m_input_done && m_output_queue.empty()) {
      if (m_input_queue.size() >= k_max_queue_size) {
        m_output_cv.wait(l, [this]() { return m_input_queue.size() < k_max_queue_size; });
      }
      read_more_input(l);
      CHECK(l.owns_lock());
    }
    if (m_input_done && m_output_queue.empty()) {
      m_output_cv.wait(l, [this]() { return m_output_done || !m_output_queue.empty(); });
    }

    if (m_output_done && m_output_queue.empty()) {
      if (m_pipeline_step_thread.valid()) {
        // Sometimes next can be called multiple times on an iterator
        // even after it raises StopIteration.  We should only reap the thread once.
        m_pipeline_step_thread.get();
      }
    } else {
      CHECK(!m_output_queue.empty());

      a = std::move(m_output_queue.front());
      m_output_queue.pop_front();
    }
  }

  if (!a) {
    PyErr_SetObject(PyExc_StopIteration, Py_None);
    throw pybind11::error_already_set();
  }
  return a;
}

void par_asm_pipeline_wrapper::read_more_input(std::unique_lock<std::mutex>& l) {
  CHECK(!m_input_done);
  CHECK(l.owns_lock());
  l.unlock();

  assembly_ptr a;
  {
    gil_scoped_acquire gil;

    PyObject* py_next = PyIter_Next(m_input_iter_obj.ptr());
    if (PyErr_Occurred()) {
      throw pybind11::error_already_set();
    }
    if (py_next) {
      object next_obj = reinterpret_steal<object>(py_next);
      a = cast<assembly_ptr>(next_obj);
      if (!a) {
        throw(io_exception("Assemblies must not be blank"));
      }

      check_assembly_from_user(*a);
      aoffset_t min_left = min(a->left_offset, a->right_offset);
      if (min_left < m_last_left_offset) {
        throw(io_exception(printstring(  //
            "Assemblies must be sorted in order; "
            "got assembly with offset %d after assembly with offset %d: %s",
            aoffset_t(a->left_offset), m_last_left_offset, dump_assembly_and_vars(*a).c_str())));
      }
      m_last_left_offset = min_left;
    }
  }

  l.lock();
  if (a) {
    CHECK(!m_input_done);
    m_input_queue.emplace_back(std::move(a));
    m_input_cv.notify_one();
  } else {
    CHECK(!m_input_done);
    m_input_done = true;
    m_input_cv.notify_one();
  }
}

namespace {

class output_done_setter {
 public:
  output_done_setter(std::unique_lock<std::mutex>& l, bool& output_done,
                     std::condition_variable& output_cv)
      : m_l(l), m_output_done(output_done), m_output_cv(output_cv) {}
  // Make sure m_output_done is set even in the case of exceptions.
  ~output_done_setter() {
    if (!m_l) {
      m_l.lock();
    }
    m_output_done = true;
    m_output_cv.notify_one();
  }
  std::unique_lock<std::mutex>& m_l;
  bool& m_output_done;
  std::condition_variable& m_output_cv;
};

}  // namespace

void par_asm_pipeline_wrapper::pipeline_step_thread() {
  std::unique_lock<std::mutex> l(m_mu);
  // Make sure m_output_done is set even in the case of exceptions.
  output_done_setter set_output_done(l, m_output_done, m_output_cv);
  for (;;) {
    if (m_input_queue.empty()) {
      m_input_cv.wait(l, [this]() { return m_input_done || !m_input_queue.empty(); });
    }

    if (m_input_queue.empty()) {
      CHECK(m_input_done);
      break;
    }

    assembly_ptr a = std::move(m_input_queue.front());
    m_input_queue.pop_front();
    if (m_input_queue.empty()) {
      // Make sure we don't deadlock if we didn't output anything since the input queue was full.
      m_output_cv.notify_one();
    }
    l.unlock();
    on_input(std::move(a));
    l.lock();
  }

  CHECK(m_input_done);
  CHECK(m_input_queue.empty());
  CHECK(!m_output_done);
  l.unlock();
  on_input_done();
  l.lock();
  CHECK(!m_output_done);
}
