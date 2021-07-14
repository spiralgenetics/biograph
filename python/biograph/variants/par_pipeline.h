#pragma once

#include <pybind11/pybind11.h>
#include <memory>

#include "modules/variants/assemble.h"
#include "python/biograph/biograph.h"
#include "python/biograph/reference.h"
#include "python/biograph/variants/assembly.h"
#include "python/biograph/variants/pipeline_common.h"

// Version of asm_pipeline_wrapper that runs potentially slow pipeline
// steps in a separate thread.
class __attribute__((visibility("hidden"))) par_asm_pipeline_wrapper
    : public std::enable_shared_from_this<par_asm_pipeline_wrapper> {
 public:
  // Maximum number of entries in the input queue to a pipeline step
  // before blocking.
  static constexpr size_t k_max_queue_size = 4096;

  // Returns the next object that's output from this pipeline step.
  variants::assembly_ptr next();

  // Called when there's a new assembly available for this step to process
  virtual void on_input(variants::assembly_ptr a) = 0;

  // Called when there are no more assemblies for this step to process.
  virtual void on_input_done() = 0;

  virtual ~par_asm_pipeline_wrapper() = default;

 protected:
  // input is a generator that provides assemblies for us.
  par_asm_pipeline_wrapper(pybind11::object input) { m_input_iter_obj = input.attr("__iter__")(); }

  // Returns a pipeline step this step should output assemblies to.
  variants::pipeline_step_t make_pipeline_output();

 protected:
  // True if we discard reference-only assemblies when sending to python.
  bool m_discard_reference_only = false;

 private:
  // Gathers some input from m_input_iter_obj and sends it to the pipeline step.
  void read_more_input(std::unique_lock<std::mutex>& l);

  // Runs the asynchronous thread that processes this pipeline step.
  void pipeline_step_thread();

  std::mutex m_mu;
  std::condition_variable m_input_cv;
  std::condition_variable m_input_queue_drained_cv;
  std::condition_variable m_output_cv;
  bool m_thread_started = false;
  std::future<void> m_pipeline_step_thread;

  // Python generator object that provides input to this wrapper
  pybind11::object m_input_iter_obj;

  // Last left offset seen, to ensure we don't go backwards.
  variants::aoffset_t m_last_left_offset = std::numeric_limits<variants::aoffset_t>::min();

  // Assemblies that have been read from the input to this step but
  // have not yet been processed by the pipeline step.
  std::deque<variants::assembly_ptr> m_input_queue;

  // Assemblies that have been outputted from this pipeline step but
  // has not yet been consumed by callers to this generator.
  std::deque<variants::assembly_ptr> m_output_queue;

  // True if the input to this pipeline step has read all the input elements.
  bool m_input_done = false;

  // True if this pipeline step has finished emitting all its outputs.
  bool m_output_done = false;
};

void bind_pair_edge_cov();
