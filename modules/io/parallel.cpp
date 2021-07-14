#include "modules/io/parallel.h"
#include "base/base.h"
#include "modules/io/make_unique.h"

#include <boost/format.hpp>
#include <iostream>
#include <thread>

#ifdef _OPENMP
#include <omp.h>
#endif

size_t g_parallel_splits = 100000;
size_t g_subwork_parallel_splits = 100;

constexpr double thread_pool::k_subwork_progress_portion;

namespace {

void enforce_progress_bounds(double& p) {
  // Correct for potential floating point errors
  if (p < 0) {
    DCHECK_GT(p, -0.000001);
    p = 0;
  }
  if (p > 1) {
    DCHECK_LT(p, 1.000001);
    p = 1;
  }
}

}  // namespace

std::function<void(size_t, size_t, parallel_state& state)> parallel_for_func::get() {
  if (!m_chunk_f) {
    CHECK(m_individual_f);
    m_chunk_f = [this](size_t start, size_t limit, parallel_state& state) {
      while (start != limit) {
        CHECK_LT(start, limit);
        if (state.pool()->idle_threads() > 0 && start + 1 < limit) {
          // Some threads are idle; recursively break this up into more chunks.
          parallel_for(start, limit, m_chunk_f);
          return;
        }
        m_individual_f(start, state);
        ++start;
      }
    };
  }
  return m_chunk_f;
}

namespace {

struct chunk {
  size_t start;
  size_t limit;
};

std::vector<chunk> generate_chunks(size_t start, size_t limit, size_t nsplits) {
  CHECK_LE(start, limit);
  CHECK_GT(nsplits, 0);

  std::vector<chunk> chunks;

  size_t whole_range_size = limit - start;

  for (size_t i = 0; i < nsplits; ++i) {
    size_t chunk_start_offset = whole_range_size * i / nsplits;
    size_t chunk_limit_offset = whole_range_size * (i + 1) / nsplits;

    if (chunk_start_offset == chunk_limit_offset) {
      // Empty range; skip
      continue;
    }

    size_t chunk_start = start + chunk_start_offset;
    size_t chunk_limit = start + chunk_limit_offset;
    chunks.emplace_back(chunk{.start = chunk_start, limit = chunk_limit});
  }
  return chunks;
}

thread_local parallel_state* tl_state = nullptr;
thread_local std::unique_ptr<parallel_local> tl_local;

// two threads minimum
size_t g_num_threads = std::max(2U, std::thread::hardware_concurrency());

}  // namespace

std::unique_ptr<parallel_local>& thread_pool::get_untyped_local() { return tl_local; }

thread_pool& parallel_pool() {
  static std::once_flag once;
  // Don't use unique_ptr; we don't want to try to destruct it when we exit.
  static thread_pool* g_thread_pool = nullptr;

  std::call_once(once, []() {
    CHECK(!g_thread_pool);
    g_thread_pool = new thread_pool();
  });
  CHECK(g_thread_pool);

  return *g_thread_pool;
}

void parallel_for(size_t start, size_t limit, parallel_for_func process) {
  if (parallel_pool().get_state()) {
    // if we're in a work item, guess that we're doing 1/5th of it with this additional parallel
    // work?
    parallel_for_subprogress(start, limit, process, thread_pool::k_subwork_progress_portion);
  } else {
    parallel_for(start, limit, process, null_progress_handler);
  }
}

namespace {

std::vector<thread_pool::work_t> make_parallel_for_worklist_internal(
    size_t start, size_t limit, parallel_for_func& process,
    size_t max_num_chunks = std::numeric_limits<size_t>::max()) {
  auto process_range = process.get();
  std::vector<thread_pool::work_t> worklist;
  size_t nsplits = parallel_pool().get_state() ? g_subwork_parallel_splits : g_parallel_splits;
  if (nsplits > max_num_chunks) {
    CHECK_GT(max_num_chunks, 0);
    nsplits = max_num_chunks;
  }
  for (chunk ch : generate_chunks(start, limit, nsplits)) {
    worklist.emplace_back(
        [ch, process_range](parallel_state& st) { process_range(ch.start, ch.limit, st); });
  }
  return worklist;
}

}  // namespace

void parallel_for_subprogress(size_t start, size_t limit, parallel_for_func process,
                              double progress_subpart) {
  parallel_state* st = parallel_pool().get_state();
  CHECK(st) << "Cannot supply subprogress except inside of a parallel work item";
  auto process_range = process.get();
  auto worklist = make_parallel_for_worklist_internal(start, limit, process);

  parallel_pool().execute_worklist(worklist, parallel_pool().current_priority() + 1,
                                   progress_subpart);
}

void parallel_for(size_t start, size_t limit, parallel_for_func process,
                  progress_handler_t progress) {
  auto process_range = process.get();
  auto worklist = make_parallel_for_worklist_internal(start, limit, process);
  parallel_pool().execute_worklist(worklist, progress);
}

std::vector<thread_pool::work_t> make_parallel_for_worklist(size_t start, size_t limit,
                                                            parallel_for_func process,
                                                            size_t max_num_chunks) {
  return make_parallel_for_worklist_internal(start, limit, process, max_num_chunks);
}

void thread_pool::start_progress(progress_handler_t progress) {
  std::lock_guard<std::mutex> l(m_mu);
  CHECK(!m_progress_f) << "Only one progress allowed at once";
  m_progress_f.emplace(progress);
  m_tot_progress = 0;
}

void thread_pool::finish_progress() {
  std::lock_guard<std::mutex> l(m_mu);
  CHECK(m_progress_f) << "Missing progress handler";
  m_progress_f.reset();
}

void thread_pool::execute_worklist_internal(std::vector<work_t> new_worklist, int priority,
                                            double progress_subpart) {
  if (new_worklist.empty()) {
    return;
  }

  parallel_state* st = get_state();
  if (st) {
    progress_subpart *= st->m_progress_part;
    st->m_progress_part -= progress_subpart;
  }

  balance_worklist_progress(progress_subpart, new_worklist);

  std::unique_lock<std::mutex> l(m_mu);
  auto& worklist = m_work[priority];
  if (worklist.empty()) {
    m_state_changed.notify_all();
  }

  // If we get an exception we may abort processing this work and
  // return from this function before all the generated works finish.
  // So we want to supply them a pointer that will stick around.
  std::shared_ptr<size_t> work_left = std::make_shared<size_t>(new_worklist.size());

  for (const work_t& work : new_worklist) {
    work_t new_work = work;
    new_work.f = [this, work, work_left](parallel_state& state) {
      try {
        work.f(state);
      } catch (...) {
        std::lock_guard<std::mutex> l(m_mu);
        CHECK_GT(*work_left, 0);
        --*work_left;
        throw;
      }
      std::lock_guard<std::mutex> l(m_mu);
      CHECK_GT(*work_left, 0);
      if (0 == --*work_left) {
        m_state_changed.notify_all();
      }
    };

    worklist.emplace_back(new_work);
    ++m_queued_work;
  }

  m_more_work.notify_all();

  if (st) {
    // Worker thread or only thread; keep using this thread to work.
    CHECK(l.owns_lock());
    while (*work_left && !m_exception) {
      boost::optional<std::pair<int, work_t>> maybe_work = get_work(l, priority);
      if (maybe_work) {
        run_one_work_and_catch(l, maybe_work->first /* priority */, maybe_work->second);
        if (!l) {
          l.lock();
        }
      } else {
        CHECK(m_cur_threads);
        m_state_changed.wait(l);
      }
    }
  } else {
    // Main thread; wait all threads to finish and clean up
    size_t num_prio = m_work.size();
    if (num_prio) {
      double progress_per_prio = 1.0 / num_prio;
      for (auto& work_prio : m_work) {
        balance_worklist_progress(progress_per_prio, work_prio.second);
      }
    }
    start_threads(l);
    finish_threads(l);
  }
  if (!l) {
    l.lock();
  }
  if (!m_exception) {
    CHECK_EQ(*work_left, 0);
  }
}

parallel_state* thread_pool::get_state() const { return tl_state; }

size_t get_thread_count() { return g_num_threads; }

size_t set_thread_count(size_t num_threads) {
  return set_thread_count(boost::str(boost::format("%d") % num_threads));
}

size_t set_thread_count(const std::string& requested) {
  if (requested == "auto") {
    g_num_threads = std::max(2U, std::thread::hardware_concurrency());
  } else {
    try {
      g_num_threads = std::stoi(requested);
    } catch (const std::exception& ex) {
      throw std::runtime_error("--threads must specify an integer >= 1");
    }
  }

#ifdef _OPENMP
  if (g_num_threads < 1) {
    throw std::runtime_error("--threads must specify an integer >= 1");
  }

  if (g_num_threads >= 1) {
    omp_set_num_threads(g_num_threads);
  }
#endif

  return g_num_threads;
}

void thread_pool::start_threads(std::unique_lock<std::mutex>& l) {
  CHECK(l.owns_lock());
  CHECK_GT(g_num_threads, 0) << "Parallelism requested with no threads";
  CHECK_EQ(0, m_cur_threads) << "Threads already started?";
  CHECK_EQ(0, m_wanted_threads) << "Threads already started?";
  CHECK(m_threads.empty());

  m_wanted_threads = g_num_threads;
  while (m_cur_threads < m_wanted_threads) {
    ++m_cur_threads;
    std::thread new_thread([this]() {
      std::unique_lock<std::mutex> l(m_mu);
      run_worker(l);
      if (!l) {
        l.lock();
      }
      --m_cur_threads;
      m_state_changed.notify_all();
    });
    m_threads.emplace_back(std::move(new_thread));
  }
}

void thread_pool::run_one_work_and_catch(std::unique_lock<std::mutex>& l, int priority,
                                         const work_t& work) {
  CHECK(l.owns_lock());
  try {
    run_one_work(l, priority, work);
  } catch (...) {
    m_exception.emplace(std::current_exception());
    m_state_changed.notify_all();
    m_more_work.notify_all();
  }
}

boost::optional<std::pair<int /* priority */, thread_pool::work_t>> thread_pool::get_work(
    std::unique_lock<std::mutex>& l, int min_priority) {
  CHECK(l.owns_lock());
  boost::optional<std::pair<int, work_t>> work;
  if (m_work.empty()) {
    return boost::none;
  }

  auto it = m_work.begin();
  CHECK(it != m_work.end());
  int priority = it->first;
  if (priority < min_priority) {
    return boost::none;
  }
  auto& worklist = it->second;
  CHECK(!worklist.empty());

  auto& maybe_work = worklist.front();
  CHECK_LE(m_memory_reserved, m_memory_limit);
  if (maybe_work.reserve_memory + m_memory_reserved > m_memory_limit) {
    return boost::none;
  }

  work.emplace(priority, std::move(maybe_work));
  worklist.pop_front();
  if (worklist.empty()) {
    m_work.erase(it);
  }

  CHECK_GT(m_queued_work, 0);
  --m_queued_work;

  return work;
}

void thread_pool::run_one_work(std::unique_lock<std::mutex>& l, int priority, const work_t& work) {
  CHECK(l.owns_lock());

  ++m_active_work;

  parallel_state* orig_st = tl_state;
  parallel_state st;
  st.m_pool = this;
  st.m_cur_priority = priority;
  tl_state = &st;

  m_memory_reserved += work.reserve_memory;
  st.m_memory_reserved = work.reserve_memory;
  st.m_progress_part = work.progress_part;

  l.unlock();
  try {
    work.f(st);
  } catch (...) {
    note_work_finished(l, st, work);
    tl_state = orig_st;
    throw;
  }
  note_work_finished(l, st, work);
  tl_state = orig_st;
}

void thread_pool::note_work_finished(std::unique_lock<std::mutex>& l, parallel_state& st,
                                     const work_t& work) {
  if (!l) {
    l.lock();
  }
  if (st.m_memory_reserved) {
    CHECK_GE(m_memory_reserved, st.m_memory_reserved);
    m_memory_reserved -= st.m_memory_reserved;
    st.m_memory_reserved = 0;

    m_more_work.notify_all();
  }

  CHECK_GT(m_active_work, 0);
  --m_active_work;

  enforce_progress_bounds(st.m_progress_part);
  if (st.m_progress_part > 0) {
    m_tot_progress += st.m_progress_part;
    st.m_progress_part = 0;
    enforce_progress_bounds(m_tot_progress);
    m_new_progress = true;
  }
  CHECK(l.owns_lock());
  m_state_changed.notify_all();
}

void thread_pool::run_worker(std::unique_lock<std::mutex>& l) {
  CHECK(!tl_state);

  if (!l) {
    l.lock();
  }
  while (m_cur_threads <= m_wanted_threads) {
    CHECK(l.owns_lock());
    try {
      boost::optional<std::pair<int, work_t>> work;
      if (!m_exception) {
        work = get_work(l);
      }
      if (!work) {
        ++m_idle_threads;
        m_state_changed.notify_all();
        m_more_work.wait(l);
        CHECK_GT(m_idle_threads, 0);
        --m_idle_threads;
        continue;
      }

      run_one_work_and_catch(l, work->first /* priority */, work->second);
    } catch (const std::exception& e) {
      LOG(FATAL) << "Unhandled exception occured during parallel processing: " << e.what();
    }

    if (!l) {
      l.lock();
    }
  }

  auto& loc = get_untyped_local();
  if (loc) {
    l.unlock();
    try {
      loc->flush();
    } catch (...) {
      l.lock();
      m_exception.emplace(std::current_exception());
      m_state_changed.notify_all();
      m_more_work.notify_all();
    }
  }
  loc.reset();

  CHECK(!loc) << "State local should have been cleared already";
  CHECK(!tl_state) << "Parallel state should have been cleared already";
}

void thread_pool::throw_if_exception() {
  if (m_exception) {
    std::exception_ptr ex = *m_exception;
    m_exception.reset();
    std::rethrow_exception(ex);
  }
}

void thread_pool::check_progress_update(std::unique_lock<std::mutex>& l) {
  CHECK(l.owns_lock());
  CHECK(!get_state()) << "Progress updates should only happen in main thread";
  while (m_new_progress) {
    m_new_progress = false;

    double notify_progress = m_tot_progress;

    l.unlock();

    CHECK(m_progress_f) << "Missing progress handler";
    (*m_progress_f)(notify_progress);

    l.lock();
  }
}

void thread_pool::finish_threads(std::unique_lock<std::mutex>& l) {
  CHECK(!get_state()) << "finish_threads should only be executed in main thread";
  CHECK(l.owns_lock());

  CHECK(m_cur_threads) << "No threads present to process work?";

  m_state_changed.wait(l, [this, &l]() {
    check_progress_update(l);
    CHECK(l.owns_lock());

    if ((m_exception || m_work.empty()) && m_idle_threads == m_cur_threads) {
      return true;
    }

    CHECK(m_cur_threads) << "No threads present to process work?";
    return false;
  });

  check_progress_update(l);
  CHECK(l.owns_lock());

  // Now that all work is complete, take down our threads.  We can't
  // keep them along since we might fork or something...
  m_wanted_threads = 0;
  m_more_work.notify_all();
  m_state_changed.notify_all();

  std::vector<std::thread> reap = std::move(m_threads);
  m_threads.clear();
  l.unlock();

  while (!reap.empty()) {
    reap.back().join();
    reap.pop_back();
  }
  l.lock();

  CHECK(m_threads.empty());
  CHECK_EQ(0, m_cur_threads);
  CHECK_EQ(0, m_memory_reserved);
  if (m_exception) {
    m_work.clear();
    m_queued_work = 0;
  } else {
    CHECK_EQ(0, m_queued_work);
  }

  throw_if_exception();
}

void thread_pool::add_work_async(work_t work, int priority) {
  parallel_state* st = get_state();

  if (st) {
    enforce_progress_bounds(work.progress_part);

    work.progress_part *= st->m_progress_part;
    st->m_progress_part -= work.progress_part;

    enforce_progress_bounds(st->m_progress_part);
  }

  std::unique_lock<std::mutex> l(m_mu);
  auto& worklist = m_work[priority];
  m_more_work.notify_one();
  if (worklist.empty()) {
    m_state_changed.notify_all();
  }
  worklist.emplace_back(std::move(work));
  ++m_queued_work;

  if (!st) {
    // Haven't started processing yet; don't do anything.
    return;
  }

  // We're in a thread; execute one and make sure we don't build up a
  // huge backlog of async tasks.
  if (priority > st->priority() && worklist.size() > 2 * m_cur_threads) {
    boost::optional<std::pair<int, work_t>> maybe_work = get_work(l, priority);
    if (maybe_work) {
      run_one_work_and_catch(l, maybe_work->first /* priority */, maybe_work->second);
    }
  }
}

int thread_pool::current_priority() const {
  parallel_state* st = get_state();
  if (st) {
    return st->priority();
  } else {
    return 0;
  }
}

unsigned thread_pool::idle_threads() const { return m_idle_threads; }

void thread_pool::set_memory_limit(size_t new_limit) { m_memory_limit = new_limit; }

void parallel_state::unreserve_memory(size_t size) {
  if (!size) {
    return;
  }
  CHECK_GE(m_memory_reserved, size);
  m_memory_reserved -= size;
  pool()->unreserve_memory(size);
}

void thread_pool::unreserve_memory(size_t size) {
  if (!size) {
    return;
  }
  std::lock_guard<std::mutex> l(m_mu);

  CHECK_GE(m_memory_reserved, size);
  m_memory_reserved -= size;
  m_more_work.notify_all();
  m_state_changed.notify_all();
}
