#include "python/biograph/variants/discover.h"

#include "modules/variants/dedup.h"
#include "modules/variants/scaffold.h"
#include "modules/variants/sort.h"
#include "modules/variants/trace_ref.h"
#include "modules/variants/trim_ref.h"
#include "python/biograph/variants/assembly.h"
#include "python/biograph/variants/module.h"
#include "python/biograph/variants/pipeline.h"
#include "python/common.h"

#include <pybind11/functional.h>

using namespace pybind11;
using namespace variants;

std::shared_ptr<ref_map> ref_map_generate_from(const std::shared_ptr<seqset>& ss,
                                               const reference_wrapper& ref, object py_progress) {
  std::shared_ptr<ref_map> result = std::make_shared<ref_map>(ss.get(), ref.get_reference().get());
  execute_with_py_progress(py_progress,
                           [&](progress_handler_t progress) { result->build(progress); });
  return result;
}

void ref_map_generate_and_save(const std::shared_ptr<seqset>& ss, const reference_wrapper& ref,
                               const std::string& pathname, object py_progress) {
  std::string new_refmap = pathname + ".new";
  unlink(new_refmap.c_str());
  {
    spiral_file_create_mmap c(new_refmap);
    ref_map build_rmap(ss.get(), ref.get_reference().get(), c.create());
    execute_with_py_progress(py_progress,
                             [&](progress_handler_t progress) { build_rmap.build(progress); });
  }
  rename(new_refmap.c_str(), pathname.c_str());
}
std::shared_ptr<ref_map> ref_map_load(const std::shared_ptr<seqset>& ss,
                                      const reference_wrapper& ref, const std::string& pathname) {
  spiral_file_open_mmap o(pathname);
  std::shared_ptr<ref_map> result =
      std::make_shared<ref_map>(ss.get(), ref.get_reference().get(), o.open());
  return result;
}

class parallel_discover;
class discover_pipeline_interface : public pipeline_interface {
 public:
  discover_pipeline_interface(parallel_discover* d, const std::string& scaffold_name,
                              const assemble_options& options);
  ~discover_pipeline_interface();

  pipeline_step_t make_parallel_input() override {
    return make_unique<assemble_lambda_output>(
        [this](assembly_ptr a) {
          std::lock_guard<std::mutex> l(m_mu);
          m_step->add(std::move(a));
        },
        "parallel_discover:make_parallel_input");
  }

 private:
  parallel_discover* m_d = nullptr;
  std::string m_scaffold_name;
  assemble_options m_options;
  pipeline_step_t m_step;

  std::mutex m_mu;
};

class __attribute__((visibility("hidden"))) parallel_discover : public scaffold_pipeline_interface {
 public:
  parallel_discover(const std::shared_ptr<readmap>& rm, const reference_wrapper& ref,
                    const std::shared_ptr<ref_map>& rmap)
      : m_readmap(rm), m_ref(ref.get_reference()), m_rmap(rmap) {
    m_options.seqset = rm->get_seqset().get();
    m_options.readmap = rm.get();
    m_options.ref = m_ref.get();
    m_options.rmap = m_rmap.get();
    m_trace_ref = make_unique<trace_ref>(m_options, this);
  }

  // TODO(nils): Why do we need this?
  parallel_discover(const parallel_discover& rhs) {
    throw(io_exception("Parallel discoverer may not be copied"));
  }

  ~parallel_discover() { abort_if_running(); }

  void add_entire_reference() {
    check_not_done();
    m_trace_ref->add_entire_reference();
  }
  void add_scaffold(const std::string& scaffold_name) {
    check_not_done();
    m_trace_ref->add_scaffold(scaffold_name);
  }
  void add_scaffold_range(const std::string& scaffold_name, aoffset_t start, aoffset_t limit) {
    check_not_done();
    m_trace_ref->add_scaffold_range(scaffold_name, start, limit);
  }

  std::unique_ptr<pipeline_interface> pipeline_for_scaffold(
      const assemble_options& options, const std::string& scaffold_name) override {
    return make_unique<discover_pipeline_interface>(this, scaffold_name, options);
  }

  void check_new_progress(std::unique_lock<std::mutex>& l) {
    if (m_new_progress) {
      boost::optional<double> new_progress = m_new_progress;
      m_new_progress.reset();
      if (new_progress && m_progress) {
        l.unlock();
        {
          gil_scoped_acquire acquire_gil;
          m_progress(*new_progress);
        }
        l.lock();
      }
    }
  }

  void process_all_output() {
    std::unique_lock<std::mutex> l(m_mu);
    CHECK(m_remaining);
    while (m_remaining) {
      m_cond.wait(l, [this]() -> bool {
        return m_remaining == 0 || !m_queued.empty() || m_aborted || m_new_progress;
      });

      if (m_aborted) {
        return;
      }

      check_new_progress(l);

      while (!m_queued.empty()) {
        std::string scaffold_name = m_queued.back().first;
        assembly_ptr a = std::move(m_queued.back().second);
        m_queued.pop_back();
        l.unlock();
        {
          gil_scoped_acquire acquire_gil;
          m_process_output(scaffold_name, std::move(a));
        }
        l.lock();

        check_new_progress(l);
      }
      if (m_queued.empty()) {
        m_cond.notify_all();
      }
    }
    CHECK(m_queued.empty());
  }

  void assemble(std::function<void(const std::string&, assembly_ptr)> process_output,
                object progress_handler) {
    check_not_done();
    m_process_output = process_output;
    m_progress = progress_handler;

    if (m_trace_ref->empty()) {
      throw(io_exception("ParallelDiscover: must specify regions to discover variants in."));
    }

    {
      gil_scoped_release release_gil;
      start_assemble();
      try {
        process_all_output();
      } catch (const error_already_set&) {
        abort_if_running();
        throw;
      }
      wait_for_finish();
    }
    if (m_progress) {
      m_progress(1);
    }
  }

  void abort_if_running() {
    {
      std::lock_guard<std::mutex> l(m_mu);
      if (!m_trace_ref) {
        return;
      }
      m_trace_ref->abort_trace();
      m_aborted = true;
      m_cond.notify_all();
    }

    wait_for_finish();
  }

  void add(const std::string& scaffold_name, assembly_ptr a) {
    // Don't let discovery run away wildly if python is being slow
    // processing assemblies.
    static constexpr size_t k_queue_hiwat = 1024;

    std::unique_lock<std::mutex> l(m_mu);
    m_cond.wait(l, [this]() -> bool { return m_queued.size() < k_queue_hiwat || m_aborted; });
    if (m_aborted) {
      return;
    }
    if (m_queued.empty()) {
      m_cond.notify_all();
    }
    m_queued.emplace_back(scaffold_name, std::move(a));
  }

  void incref() {
    std::lock_guard<std::mutex> l(m_mu);
    ++m_remaining;
  }

  void decref() {
    std::lock_guard<std::mutex> l(m_mu);
    CHECK_GT(m_remaining, 0);
    --m_remaining;
    if (!m_remaining) {
      m_cond.notify_all();
    }
  }

 private:
  void check_not_done() {
    if (m_assemble_started) {
      throw(io_exception("ParallelDiscover may not be reused"));
    }
  }

  void start_assemble() {
    std::function<void(double)> progress_handler = null_progress_handler;
    if (m_progress) {
      progress_handler = [this](double new_progress) {
        std::lock_guard<std::mutex> l(m_mu);
        if (new_progress < (m_last_progress + 0.0001)) {
          // Don't bother reporting progress changes under 0.01%.
          return;
        }
        m_last_progress = new_progress;
        if (!m_new_progress) {
          m_cond.notify_all();
        }
        m_new_progress = new_progress;
      };
    }
    m_assemble_started = true;
    m_run_assemble = std::async(std::launch::async, [this, progress_handler]() {
      m_trace_ref->assemble(progress_handler);
      m_trace_ref.reset();
    });
  }

  void wait_for_finish() {
    if (m_assemble_started) {
      m_run_assemble.get();
      CHECK(!m_trace_ref);
    } else {
      m_trace_ref.reset();
    }
    if (!m_aborted) {
      CHECK_EQ(0, m_queued.size());
    }
  }

  std::shared_ptr<readmap> m_readmap;
  std::shared_ptr<reference> m_ref;
  std::shared_ptr<ref_map> m_rmap;
  assemble_options m_options;

  std::function<void(const std::string&, assembly_ptr)> m_process_output;
  object m_progress;

  std::unique_ptr<trace_ref> m_trace_ref;
  std::future<void> m_run_assemble;
  bool m_assemble_started = false;

  std::mutex m_mu;
  std::condition_variable m_cond;
  size_t m_remaining = 0;
  std::vector<std::pair<std::string, assembly_ptr>> m_queued;
  double m_last_progress = 0;
  boost::optional<double> m_new_progress;
  bool m_aborted = false;
};

discover_pipeline_interface::discover_pipeline_interface(parallel_discover* d,
                                                         const std::string& scaffold_name,
                                                         const assemble_options& options)
    : m_d(d), m_scaffold_name(scaffold_name), m_options(options) {
  m_d->incref();
  auto saver = make_unique<assemble_lambda_output>(
      [this](assembly_ptr a) { m_d->add(m_scaffold_name, std::move(a)); }, "discover_generator");
  auto dedup = make_unique<exact_deduper>(std::move(saver));
  auto trim = make_unique<ref_trimmer>(m_options, std::move(dedup));
  m_step = make_unique<sorter>(canon_assembly_order(), std::move(trim));
}

discover_pipeline_interface::~discover_pipeline_interface() {
  // Flush all assemblies through our trim and dedup.
  m_step.reset();
  // Stop waiting for new assemblies from us.
  m_d->decref();
}

void bind_discover(module& m) {
  class_<ref_map, std::shared_ptr<ref_map>>(
      m, "RefMap", "Represents a mapping of which seqset entries have reference coverage.")
      .def_static("generate_from", ref_map_generate_from, arg("seqset"), arg("reference"),
                  arg("progress"), "Generates a reference map from the given seqset and reference")
      .def_static("generate_and_save", ref_map_generate_and_save, arg("seqset"), arg("reference"),
                  arg("path"), arg("progress"), "Generates a reference map and saves it")
      .def_static("load", ref_map_load, arg("seqset"), arg("reference"), arg("path"),
                  "Loads a refmap from the given file");

  class_<parallel_discover>(m, "ParallelDiscover", "Discovers variants in parallel")
      .def(init<const std::shared_ptr<readmap>&, const reference_wrapper&,
                const std::shared_ptr<ref_map>&>())
      .def("add_entire_reference", &parallel_discover::add_entire_reference,
           "Queue entire reference for variant discovery")
      .def("add_scaffold", &parallel_discover::add_scaffold,
           "Queue a whole scaffold (chromsome) for variant discovery")
      .def("add_scaffold_range", &parallel_discover::add_scaffold_range,
           "Queue a part of a scaffold for variant discovery")
      .def("assemble", &parallel_discover::assemble, arg("process_output"), arg("progress"),
           "Discover variants in queued regions in parallel. process_output is a function that is "
           "called with (scaffold name, assembly) for each assembly that is discovered.  "
           "Assemblies are not guaranteed to be produced in any particular order.  Progress is "
           "called periodically with a floating point value between 0 and 1 that estimates how "
           "much of the discovery process is complete.  If do_filtering is true, attempt to do "
           "some simple filtering on assemblies before outputting.");
}
