#pragma once

#include "modules/io/progress.h"
#include "modules/io/utils.h"

#include <deque>
#include <map>
#include <memory>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#include <parallel/algorithm>
#else
#include <algorithm>
#endif

class parallel_state;
class parallel_local;
class thread_pool {
 public:
  // If unspecified, additional parallel operations within a work item
  // will be counted as this portion of the work item's progress.
  static constexpr double k_subwork_progress_portion = 0.2;

  struct work_t {
    // Bytes of memory that are reserved by this work.
    size_t reserve_memory = 0;

    // Weighted value for how much this piece of work contributes to
    // progress.  Default value is used for progress when processing
    // async work additions where we don't know how many of them will
    // occur.  When submitting a worklist, these values are normalized
    // so they add up to 1 for the wohle worklist.
    double progress_part = 0.01;

    std::function<void(parallel_state& state)> f;

    work_t(std::function<void(parallel_state& state)> f_) : f(std::move(f_)) {}
  };

  thread_pool() = default;

  // Enqueue additional work that doesn't need to be done right away;
  // just by the end of the stage.
  void add_work_async(thread_pool::work_t work) {
    add_work_async(std::move(work), current_priority() + 1);
  }
  void add_work_async(thread_pool::work_t work, int priority);

  // Returns the priority of our currently executing task, or 0 if
  // we're not inside a parallel state.
  int current_priority() const;

  template <typename Container>
  void execute_worklist(const Container& worklist) {
    if (get_state()) {
      execute_worklist(worklist, current_priority() + 1, k_subwork_progress_portion);
    } else {
      execute_worklist(worklist, null_progress_handler);
    }
  }

  template <typename Container>
  void execute_worklist(const Container& worklist, progress_handler_t progress) {
    CHECK(!get_state()) << "Cannot supply a new progress handler inside a work item";

    std::vector<work_t> converted(worklist.begin(), worklist.end());
    start_progress(progress);
    try {
      execute_worklist_internal(std::move(worklist), 0, 1.0);
      finish_progress();
    } catch (...) {
      std::lock_guard<std::mutex> l(m_mu);
      m_progress_f.reset();
      throw;
    }
  }

  template <typename Container>
  void execute_worklist(const Container& worklist, int priority, double progress_subpart) {
    CHECK(get_state()) << "Cannot supply a progress portion of a whole job";
    std::vector<work_t> converted(worklist.begin(), worklist.end());
    execute_worklist_internal(std::move(converted), priority, progress_subpart);
  }

  void set_progress(progress_handler_t progress) {
    std::lock_guard<std::mutex> l(m_mu);
    m_progress_f.emplace(progress);
    m_tot_progress = 0;
  }

  parallel_state* get_state() const;

  unsigned idle_threads() const;

  void set_memory_limit(size_t new_limit);

 private:
  friend class parallel_state;
  struct thread_done_exception : public std::exception {};

  void run_worker(std::unique_lock<std::mutex>& l);
  boost::optional<std::pair<int /* priority */, work_t>> get_work(
      std::unique_lock<std::mutex>& l, int min_priority = std::numeric_limits<int>::min());
  void run_one_work(std::unique_lock<std::mutex>& l, int priority, const work_t& work);
  void run_one_work_and_catch(std::unique_lock<std::mutex>& l, int priority, const work_t& work);
  void execute_worklist_internal(std::vector<work_t> worklist, int priority,
                                 double progress_subpart);
  void start_progress(progress_handler_t progress);
  void finish_progress();
  void throw_if_exception();
  void note_work_finished(std::unique_lock<std::mutex>& l, parallel_state& st, const work_t& work);
  void start_threads(std::unique_lock<std::mutex>& l);
  void finish_threads(std::unique_lock<std::mutex>& l);
  void check_progress_update(std::unique_lock<std::mutex>& l);
  void unreserve_memory(size_t size);
  template <typename Container>
  static void balance_worklist_progress(double tot_progress, Container& worklist);
  std::unique_ptr<parallel_local>& get_untyped_local();

  std::mutex m_mu;

  std::map<int /* priority */, std::deque<work_t>,
           std::greater<int> /* highest priority goes first */>
      m_work;

  // Count of number of work items in progress.
  unsigned m_active_work = 0;
  // Total work items queued
  unsigned m_queued_work = 0;
  // Number of threads currently waiting on work.
  std::atomic<unsigned> m_idle_threads{0};

  boost::optional<std::exception_ptr> m_exception;

  // Notified once whenever we add more work *or* there's a larger
  // change.  If we're just adding a single work, we can use
  // notify_one on here to not have to wake up all the threads.
  std::condition_variable m_more_work;

  // Notified whenever the current state changes in a way that's
  // different than a single work being added to an existing priority.
  std::condition_variable m_state_changed;

  std::vector<std::thread> m_threads;

  boost::optional<progress_handler_t> m_progress_f;
  double m_tot_progress = 0;
  // True if m_tot_progress has changed since m_progress_f was called.
  bool m_new_progress = false;

  size_t m_cur_threads = 0;
  size_t m_wanted_threads = 0;

  size_t m_memory_reserved = 0;
  size_t m_memory_limit = std::numeric_limits<size_t>::max();
};

class parallel_local {
 public:
  // Subclasses should override 'flush' to cleanup as opposed to doing
  // any work in the destructor.
  virtual void flush() {}

  virtual ~parallel_local() = default;

 protected:
  parallel_local() = default;
};

class parallel_state {
 public:
  // True if an exception has been thrown by some worker.
  bool exception_thrown() const { return m_pool->m_exception ? true : false; }

  // State shared between different chunks and the same worker.
  template <typename T, typename... Args>
  T* get_local(Args&&... args) {
    T* local = nullptr;
    std::unique_ptr<parallel_local>& untyped_local = m_pool->get_untyped_local();
    if (untyped_local) {
      local = dynamic_cast<T*>(untyped_local.get());
    }
    if (!local) {
      local = new T(std::forward<Args>(args)...);
      untyped_local.reset(local);
    }
    return local;
  }

  // Returns the current priority of the work that's being executed.
  int priority() const { return m_cur_priority; }

  thread_pool* pool() const { return m_pool; }

  void unreserve_memory(size_t size);
  size_t memory_reserved() const { return m_memory_reserved; }

 private:
  friend class thread_pool;

  void reset() { m_local.reset(); }

  int m_cur_priority;

  // Memory reserved by current task.
  size_t m_memory_reserved = 0;

  // Amount of progress remaining on the currently executing task.
 public:
  double m_progress_part = 0;

  thread_pool* m_pool = nullptr;

  std::unique_ptr<parallel_local> m_local;
};

template <typename Container>
void thread_pool::balance_worklist_progress(double tot_progress, Container& worklist) {
  CHECK(!worklist.empty());
  double sum_progress = 0;
  for (const auto& work : worklist) {
    sum_progress += work.progress_part;
  }
  if (sum_progress > 0) {
    for (auto& work : worklist) {
      work.progress_part *= tot_progress / sum_progress;
    }
  }
}

thread_pool& parallel_pool();

// Set the max number of threads
size_t set_thread_count(size_t num_threads);
size_t set_thread_count(const std::string& requested = "auto");

// Return the max number of threads
size_t get_thread_count();

class parallel_for_func {
 public:
  parallel_for_func() = delete;
  parallel_for_func(const parallel_for_func&) = default;

  // Process items individually, recursively chunking.
  template <typename T,
            typename std::enable_if<
                std::is_convertible<T, std::function<void(size_t, parallel_state&)>>::value,
                int>::type = 0>
  parallel_for_func(const T& f) {
    m_individual_f = f;
  }

  // Process items individually, recursively chunking; worker does not want state.
  template <typename T,
            typename std::enable_if<std::is_convertible<T, std::function<void(size_t)>>::value,
                                    int>::type = 0>
  parallel_for_func(const T& f) {
    m_individual_f = [f](size_t i, parallel_state&) { f(i); };
  }

  // Process items in chunks.
  template <typename T,
            typename std::enable_if<
                std::is_convertible<T, std::function<void(size_t, size_t, parallel_state&)>>::value,
                int>::type = 0>
  parallel_for_func(const T& f) {
    m_chunk_f = f;
  }

  // Process items in chunks, worker does not want state.
  template <typename T,
            typename std::enable_if<
                std::is_convertible<T, std::function<void(size_t, size_t)>>::value, int>::type = 0>
  parallel_for_func(const T& f) {
    m_chunk_f = [f](size_t start, size_t limit, parallel_state&) { f(start, limit); };
  }

  // Returns a canonical version of the function.  This
  // parallel_for_func object must outlast the return value.
  std::function<void(size_t, size_t, parallel_state& state)> get();

 private:
  std::function<void(size_t, size_t, parallel_state& state)> m_chunk_f;
  std::function<void(size_t, parallel_state& state)> m_individual_f;
};

void parallel_for(size_t start, size_t limit, parallel_for_func process);
void parallel_for(size_t start, size_t limit, parallel_for_func process,
                  progress_handler_t progress);

std::vector<thread_pool::work_t> make_parallel_for_worklist(
    size_t start, size_t limit, parallel_for_func process,
    size_t max_num_chunks = std::numeric_limits<size_t>::max());

// Same as parallel_for, but for work items that need to do additional
// parallel processing.  "progress_part" specifies how much of this
// work item's progress should be marked done by the work in the
// parallel for.
void parallel_for_subprogress(size_t start, size_t limit, parallel_for_func process,
                              double progress_subpart);

// Default number of parts to split into.  (Note: The number of parts
// that are processed at the same time is specified separately; see
// the parallel_execute implementation.)  If this is too high, we will
// spend all our time managing the paralellism.  If this is too low,
// progress will be rarely updated and we may have some CPUs that
// aren't very busy at the end for longer when most of the chunks are
// done.
extern size_t g_parallel_splits;
// Same, but for recursive calls to parallel_for within a parallel job.
extern size_t g_subwork_parallel_splits;

// Execute a sort in parallel, if enabled.  Otherwise, executes using
// normal std::sort.  Either way, it sorts it in place instead of
// allocating a significant amount of additional memory, like
// __gnu_parallel::sort can do in some run modes (e.g. multiway
// mergesort).
//
// However, there is a thread safety bug in std::partition
// (https://gcc.gnu.org/bugzilla/show_bug.cgi?id=69285), which is used
// by quicksort.  If the supplied compare function could potentially
// crash (as opposed to just returning a meaningless result) if the
// things it's comparing change while it's comparing them, it's best
// to use parallel_sort_thread_safe instead.
template <typename iterator_type,
          typename compare_type = std::less<typename iterator_type::value_type>>
inline void parallel_sort_in_place(
    iterator_type begin, iterator_type end,
    compare_type compare_f = std::less<typename iterator_type::value_type>()) {
#if defined(_OPENMP) && !defined(ADDRESS_SANITIZER) && !defined(THREAD_SANITIZER)
  // Quicksort is recursive; make sure it can use more than 2 threads
  // at once.
  int old_nested = omp_get_nested();
  omp_set_nested(1);
  __gnu_parallel::sort(begin, end, compare_f, __gnu_parallel::quicksort_tag());
  omp_set_nested(old_nested);
#else
  std::sort(begin, end, compare_f);
#endif
}

// parallel_sort_thread_safe uses multiway mergesort, which takes
// double the RAM as quicksort above.  However, multiway mergesort
// does not use parallel std::partition, so avoids the thread safety
// issues in https://gcc.gnu.org/bugzilla/show_bug.cgi?id=69285.
template <typename iterator_type,
          typename compare_type = std::less<typename iterator_type::value_type>>
inline void parallel_sort_thread_safe(
    iterator_type begin, iterator_type end,
    compare_type compare_f = std::less<typename iterator_type::value_type>()) {
#if defined(_OPENMP) && !defined(ADDRESS_SANITIZER) && !defined(THREAD_SANITIZER)
  __gnu_parallel::sort(begin, end, compare_f, __gnu_parallel::multiway_mergesort_tag());
#else
  std::sort(begin, end, compare_f);
#endif
}
