#pragma once

#include "base/base.h"

#include <atomic>
#include <cstddef>
#include <iostream>
#include <map>
#include <memory>
#include <mutex>
#include <set>
#include <unordered_map>
#include <vector>

#include <boost/program_options.hpp>

class track_mem {
 public:
  class mmap_tracker;
  class malloc_tracker;
  struct type_entry {
    std::string type;
    size_t current_usage = 0;
    size_t max_usage = 0;
  };
  class track_alloc;
  template <typename T>
  class allocator;

  static void initialize_if_needed();
  static std::vector<type_entry> get_stats();
  static void reset_stats();
  static void log_usage();

  // When used with e.g. tcmalloc's AddNewHook, this allows us to warn
  // and get a traceback if we allocate large chunks of untracked
  // memory.
  using malloc_hook_t = void (*)(const void* /* ptr */, size_t /* size */);
  static malloc_hook_t get_malloc_new_hook();
  static bool has_malloc_hook() { return g_has_malloc_hook; }
  static void set_reset_stats_hook(const std::function<void()>& new_hook) {
    g_reset_stats_hook = new_hook;
  }

  static mmap_tracker g_mmap_tracker;
  static malloc_tracker g_malloc_tracker;

 private:
  friend class track_mem_test;

  class tracked_type {
   public:
    tracked_type() = delete;

    const std::string& type() const { return m_type; }
    size_t current_usage() const {
      std::lock_guard<std::mutex> l(m_mu);
      return m_tot_alloc;
    }
    size_t max_usage() const {
      std::lock_guard<std::mutex> l(m_mu);
      return m_hiwat;
    }

    // Clear high water materk.
    void reset();

    // Write detail usage for this type to spiral log.
    void log_detail_usage(bool force_logging = false) const;
    void log_overview() const;
    virtual std::unordered_map<std::string, size_t> get_detail_usage() const = 0;

    void report_if_pending();

   protected:
    tracked_type(const std::string& type) : m_type(type) {}

    void note_allocation(size_t size, const std::string& description);
    void note_deallocation(size_t size, const std::string& description);

   private:
    virtual void on_new_hiwat(size_t new_hiwat) {}
    bool hiwat_is_interesting() const;
    void report_hiwat();

    // Type of item being tracked, e.g. "ALLOC" or "MMAP"
    std::string m_type;

    mutable std::mutex m_mu;
    size_t m_tot_alloc = 0;
    size_t m_hiwat = 0;

    // If k_min_hiwat_interval seconds from this time expires, report
    // even if we don't free anything.
    time_t m_hiwat_change_pending = 0;

    time_t m_last_hiwat_report = 0;
  };

  track_mem() = delete;
  track_mem(const track_mem&) = delete;

  static void malloc_new_hook(const void* result, size_t size);
  void reset_stats_internal();

  static std::string size_str(size_t sz);
  static void run_report_thread();

  // Don't itemize any allocations under this size.
  static constexpr size_t k_min_interesting_size = 1ULL * 1024 * 1024;
  // Only report new highwater marks this many bytes over previous
  static constexpr size_t k_min_interesting_hiwat_bytes = 1ULL * 1024 * 1024 * 1024;
  // Only report new highwater marks this many percent over previous
  static constexpr size_t k_min_interesting_hiwat_pct = 15;
  // Only report this many top memory users
  static constexpr size_t k_top_interesting = 5;
  // Report hiwater changes within this long, or after free.
  static constexpr time_t k_hiwat_report_interval = 30;
  // Check for hi water reporting every this many seconds
  static constexpr time_t k_hiwat_check_interval = 5;

  // Anything larger than this must use track_alloc instead of straight new/malloc.
  // TODO(nils): Crank this down to get better memory tracking.
  static constexpr size_t k_max_untracked_bytes = 50 * 1024 * 1024;

  // If true, kill the program when k_max_untracked_bytes is exceeded.
  // If false, only kill the problem when k_max_untracked_bytes is
  // exceeded if we're running in a test.
  static constexpr bool k_fatal_untracked = false;

  static constexpr size_t k_num_types = 2;
  static const std::array<tracked_type*, k_num_types> k_types;

  // Whether a malloc hook has been registered
  static bool g_has_malloc_hook;

  // Hook to be called on stats reset.
  static std::function<void()> g_reset_stats_hook;
};

// Returns program options suitable for configuring available memory.
boost::program_options::options_description track_mem_program_options();

// Returns the maximum number of bytes the user would like us to use.
size_t get_maximum_mem_bytes();

// Sets the maximum number of bytes the user would like to use.  Prefer
// using track_mem_program_options to calling this.
void set_maximum_mem_bytes(size_t new_max);

class track_mem::mmap_tracker : public track_mem::tracked_type {
 public:
  mmap_tracker() : tracked_type("MMAP") {}
  struct entry {
    std::string description;
    size_t size;
  };

  void* tracked_mmap(void* addr, size_t length, int prot, int flags, int fd, off_t offset,
                     const std::string& description);

  int tracked_munmap(void* addr, size_t length, const std::string& description);

  std::unordered_map<std::string, size_t> get_detail_usage() const override;

 private:
  mutable std::mutex m_mu;
  std::unordered_map<const void* /* pointer */, entry> m_mmaps;
};

class track_mem::malloc_tracker : public tracked_type {
 public:
  malloc_tracker() : tracked_type("ALLOC") {}

  struct entry {
    std::string description;
    std::atomic<size_t> ref_count{0};
    std::atomic<size_t> tot_used{0};
  };

  std::unordered_map<std::string, size_t> get_detail_usage() const override;
  void on_new_hiwat(size_t new_hiwat) override;

  // Increments the ref count before returning to caller.
  entry* get_entry(const std::string& description);

  void expect_big_malloc(size_t size) {
    if (track_mem::has_malloc_hook()) {
      std::lock_guard<std::mutex> l(g_expect_malloc_mu);
      g_expect_malloc.insert(size);
    }
  }
  void malloc_new_hook(size_t size);

  void note_malloc_allocation(entry* entry, const void* ptr, size_t size) {
    DCHECK_GE(entry->ref_count, 0U);
    entry->tot_used += size;
    note_allocation(size, entry->description);
  }

  void note_malloc_deallocation(entry* entry, const void* ptr, size_t size) {
    DCHECK_GE(entry->ref_count, 0U);
    note_deallocation(size, entry->description);
    // Can't just use -- since we want to be able to do error checking.
    for (;;) {
      size_t old_used = entry->tot_used;
      CHECK_GE(old_used, size);
      size_t new_used = old_used - size;
      if (entry->tot_used.compare_exchange_weak(old_used, new_used)) {
        break;
      }
    }
  }

  void inc_ref(entry* entry) { ++entry->ref_count; }
  void dec_ref(entry* entry) {
    // Can't just use -- since we want to be able to do error checking.
    std::unique_lock<std::mutex> l(m_mu, std::defer_lock);
    size_t new_count;
    for (;;) {
      size_t old_count = entry->ref_count.load();
      CHECK_GT(old_count, 0U) << entry->description;
      new_count = old_count - 1;
      if (new_count == 0) {
        l.lock();
      }
      if (entry->ref_count.compare_exchange_weak(old_count, new_count)) {
        break;
      }
      if (new_count == 0) {
        l.unlock();
      }
    }

    if (new_count == 0) {
      CHECK(l.owns_lock());
      CHECK_EQ(entry->tot_used, 0U) << entry->description;
      CHECK_EQ(entry->ref_count, 0U) << entry->description;
      auto it = m_allocators.find(entry->description);
      CHECK_EQ(it->second.get(), entry) << entry->description;
      m_allocators.erase(it);
    }
  }

 private:
  std::mutex g_expect_malloc_mu;
  std::multiset<size_t> g_expect_malloc;

  mutable std::mutex m_mu;
  std::unordered_map<std::string /* description */, std::unique_ptr<entry>> m_allocators;
};

using track_alloc = track_mem::track_alloc;
void* tracked_mmap(void* addr, size_t length, int prot, int flags, int fd, off_t offset,
                   const std::string& description);
int tracked_munmap(void* addr, size_t length, const std::string& description);

// Easy way to specify a description for tracking allocation
class track_mem::track_alloc {
 public:
  track_alloc() = delete;
  track_alloc(const std::string& description) : m_entry(g_malloc_tracker.get_entry(description)) {}

  template <typename T>
  operator allocator<T>() const;

  track_alloc(track_mem::malloc_tracker::entry* entry) : m_entry(entry) {
    g_malloc_tracker.inc_ref(m_entry);
  }
  ~track_alloc() { g_malloc_tracker.dec_ref(m_entry); }
  track_alloc(const track_alloc& alloc) : m_entry(alloc.m_entry) {
    g_malloc_tracker.inc_ref(m_entry);
  }
  track_alloc& operator=(const track_alloc& rhs) = delete;
  bool operator==(const track_alloc& rhs) const { return m_entry == rhs.m_entry; }

  track_mem::malloc_tracker::entry* get() const { return m_entry; }

 private:
  track_mem::malloc_tracker::entry* m_entry = nullptr;
};

// An allocator for STL containers that tracks memory usage.  This
// should only be used for long lived containers, since there is a
// performance penalty for tracking usage.
template <typename T>
class track_mem::allocator {
 public:
  using value_type = T;

  allocator() = delete;
  allocator(const allocator&) = default;
  allocator(const track_alloc& alloc) : m_alloc(alloc) {}
  template <typename U>
  allocator(const track_mem::allocator<U>& u) : m_alloc(u.get()) {}

  void note_external_allocation(T* ptr, std::size_t n) {
    size_t alloc_size = n * sizeof(T);
    if (alloc_size > k_max_untracked_bytes) {
      g_malloc_tracker.note_malloc_allocation(m_alloc.get(), ptr, alloc_size);
    }
  }

  void note_external_deallocation(T* ptr, std::size_t n) {
    size_t alloc_size = n * sizeof(T);
    if (alloc_size > k_max_untracked_bytes) {
      g_malloc_tracker.note_malloc_deallocation(m_alloc.get(), ptr, alloc_size);
    }
  }

  T* allocate(std::size_t n) {
    size_t alloc_size = n * sizeof(T);
    if (alloc_size > k_max_untracked_bytes) {
      g_malloc_tracker.expect_big_malloc(alloc_size);
    }
    T* ptr = m_std_alloc.allocate(n);
    if (alloc_size > k_max_untracked_bytes) {
      g_malloc_tracker.note_malloc_allocation(m_alloc.get(), ptr, alloc_size);
    }
    return ptr;
  }
  void deallocate(T* ptr, std::size_t n) {
    size_t alloc_size = n * sizeof(T);
    if (alloc_size > k_max_untracked_bytes) {
      g_malloc_tracker.note_malloc_deallocation(m_alloc.get(), ptr, alloc_size);
    }
    m_std_alloc.deallocate(ptr, n);
  }

  template <typename U>
  bool operator==(const track_mem::allocator<U>& a) const {
    return m_alloc == a.m_alloc;
  }

  const track_alloc& get() const { return m_alloc; }

#if (defined(__GNUC__) && __GNUC__ < 5)
  // Yuck, have to be compatible with old stl.
  template <typename U>
  struct rebind {
    using other = track_mem::allocator<U>;
  };
  using pointer = T*;
  using const_pointer = const T*;
  using reference = T&;
  using const_reference = const T&;
  using size_type = std::size_t;
  using difference_type = std::ptrdiff_t;
  using is_always_equal = std::false_type;

  size_type max_size() const { return m_std_alloc.max_size(); }
  template <class U>
  void destroy(U* p) {
    return m_std_alloc.destroy(p);
  }
  template <class U, class... Args>
  void construct(U* p, Args&&... args) {
    m_std_alloc.construct(p, std::forward<Args>(args)...);
  }
  pointer address(reference x) const { return m_std_alloc.reference(x); }
  const_pointer address(const_reference x) const { return m_std_alloc.const_reference(x); }
#endif

 private:
  track_alloc m_alloc;
  std::allocator<T> m_std_alloc;
};

template <typename T>
track_mem::track_alloc::operator allocator<T>() const {
  return allocator<T>(*this);
};

// Wrapper around std::vector that has easy initializers for using track_mem::allocator.
template <typename T>
class tracked_vector : public std::vector<T, track_mem::allocator<T>> {
 public:
  tracked_vector() = delete;
  tracked_vector(const tracked_vector&) = default;
  tracked_vector(const track_alloc& alloc)
      : std::vector<T, track_mem::allocator<T>>(track_mem::allocator<T>(alloc)) {}
  tracked_vector(size_t n, const track_alloc& alloc)
      : std::vector<T, track_mem::allocator<T>>(track_mem::allocator<T>(alloc)) {
    this->resize(n);
  }
};

// Wrapper around std::unordered_map that has easy initializers for using track_mem::allocator.
template <typename K, typename V, typename H = std::hash<K>, typename EQ = std::equal_to<K>,
          typename VALTYPE = typename std::unordered_multimap<K, V, H, EQ>::value_type>
class tracked_unordered_map
    : public std::unordered_map<K, V, H, EQ, track_mem::allocator<VALTYPE>> {
 public:
  tracked_unordered_map() = delete;
  tracked_unordered_map(const tracked_unordered_map&) = default;
#if (defined(__GNUC__) && __GNUC__ < 5)
  // Old std::unordered_map doesn't have the simple constructor.
  tracked_unordered_map(const track_alloc& alloc)
      : std::unordered_map<K, V, H, EQ, track_mem::allocator<VALTYPE>>(
            1000, std::hash<K>(), std::equal_to<K>(), track_mem::allocator<VALTYPE>(alloc)) {}
#else
  tracked_unordered_map(const track_alloc& alloc)
      : std::unordered_map<K, V, H, EQ, track_mem::allocator<VALTYPE>>(
            track_mem::allocator<VALTYPE>(alloc)) {}
#endif
};

// Wrapper around std::unordered_multimap that has easy initializers for using track_mem::allocator.
template <typename K, typename V, typename H = std::hash<K>, typename EQ = std::equal_to<K>,
          typename VALTYPE = typename std::unordered_multimap<K, V, H, EQ>::value_type>
class tracked_unordered_multimap
    : public std::unordered_multimap<K, V, H, EQ, track_mem::allocator<VALTYPE>> {
 public:
  tracked_unordered_multimap(const tracked_unordered_multimap&) = default;
#if (defined(__GNUC__) && __GNUC__ < 5)
  // Old std::unordered_multimap doesn't have the simple constructor.
  tracked_unordered_multimap(const track_alloc& alloc)
      : std::unordered_multimap<K, V, H, EQ, track_mem::allocator<VALTYPE>>(
            1000, std::hash<K>(), std::equal_to<K>(), track_mem::allocator<VALTYPE>(alloc)) {}
#else
  tracked_unordered_multimap(const track_alloc& alloc)
      : std::unordered_multimap<K, V, H, EQ, track_mem::allocator<VALTYPE>>(
            track_mem::allocator<VALTYPE>(alloc)) {}
#endif
};
