#include "python/common.h"

void execute_with_py_progress(pybind11::object py_progress_f,
                              const std::function<void(progress_handler_t progress)> f) {
  std::mutex mu;
  double last_progress = 0;
  boost::optional<double> pending_progress;
  bool done = false;
  std::condition_variable cond;

  progress_handler_t progress_handler = [&](double new_progress) {
    std::lock_guard<std::mutex> l(mu);

    if (new_progress < last_progress + 0.0001) {
      return;
    }
    last_progress = new_progress;
    pending_progress = new_progress;
    cond.notify_all();
  };

  std::future<void> run_thread = std::async(std::launch::async, [&]() {
    f(progress_handler);
    std::lock_guard<std::mutex> l(mu);
    done = true;
    cond.notify_all();
  });

  {
    std::unique_lock<std::mutex> l(mu);

    while (!done) {
      cond.wait(l, [&]() { return done || pending_progress; });

      if (done) {
        break;
      }

      CHECK(pending_progress);
      double new_progress = *pending_progress;
      pending_progress.reset();

      l.unlock();
      pybind11::gil_scoped_acquire acquire_gil;
      py_progress_f(new_progress);
      l.lock();
    }
  }

  run_thread.get();
}
