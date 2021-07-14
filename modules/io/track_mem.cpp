// Apparently utilities.h is needed first...
#include "vendor/glog/src/utilities.h"

#include "modules/io/io.h"
#include "modules/io/stopwatch.h"
#include "modules/io/track_mem.h"

#include "base/base.h"
#include "modules/io/log.h"
#include "modules/io/make_unique.h"
#include "modules/io/utils.h"

#include <sys/mman.h>
#include <mutex>
#include <set>
#include <thread>

constexpr size_t track_mem::k_min_interesting_size;
constexpr size_t track_mem::k_min_interesting_hiwat_bytes;
constexpr size_t track_mem::k_min_interesting_hiwat_pct;
constexpr size_t track_mem::k_top_interesting;
constexpr time_t track_mem::k_hiwat_report_interval;
constexpr time_t track_mem::k_hiwat_check_interval;
constexpr size_t track_mem::k_max_untracked_bytes;
constexpr bool track_mem::k_fatal_untracked;

track_mem::malloc_tracker track_mem::g_malloc_tracker;
track_mem::mmap_tracker track_mem::g_mmap_tracker;

bool track_mem::g_has_malloc_hook = false;

std::function<void()> track_mem::g_reset_stats_hook;

const std::array<track_mem::tracked_type*, track_mem::k_num_types> track_mem::k_types{
    {&track_mem::g_mmap_tracker, &track_mem::g_malloc_tracker}};

std::vector<track_mem::type_entry> track_mem::get_stats() {
  std::vector<track_mem::type_entry> result;
  for (const auto& type : k_types) {
    track_mem::type_entry e;
    e.type = type->type();
    e.current_usage = type->current_usage();
    e.max_usage = type->max_usage();
    result.push_back(e);
  }
  return result;
}

void* tracked_mmap(void* addr, size_t length, int prot, int flags, int fd, off_t offset,
                   const std::string& description) {
  return track_mem::g_mmap_tracker.tracked_mmap(addr, length, prot, flags, fd, offset, description);
}

int tracked_munmap(void* addr, size_t length, const std::string& description) {
  return track_mem::g_mmap_tracker.tracked_munmap(addr, length, description);
}

namespace {

std::once_flag g_track_mem_init_once;
boost::optional<std::thread> g_report_thread;

}  // namespace

void track_mem::initialize_if_needed() {
  std::call_once(g_track_mem_init_once, []() {
    g_report_thread.emplace(run_report_thread);
    g_report_thread->detach();
  });
}

void track_mem::run_report_thread() {
  for (;;) {
    std::this_thread::sleep_for(std::chrono::seconds(k_hiwat_check_interval));
    for (const auto& type : k_types) {
      type->report_if_pending();
    }
  }
}

void track_mem::tracked_type::report_if_pending() {
  try {
    time_t now = time(0);
    std::lock_guard<std::mutex> l(m_mu);
    if (m_hiwat_change_pending && m_hiwat_change_pending + k_hiwat_report_interval < now) {
      report_hiwat();
    }
  } catch (const std::exception& ex) {
    std::cerr << "Memory tracking error: " << ex.what() << std::endl << std::endl;
  }
}

bool track_mem::tracked_type::hiwat_is_interesting() const {
  if (m_tot_alloc < m_hiwat + k_min_interesting_hiwat_bytes) {
    return false;
  }
  if (m_tot_alloc < m_hiwat * (100 + k_min_interesting_hiwat_pct) / 100) {
    return false;
  }
  return true;
}

void track_mem::tracked_type::report_hiwat() {
  m_hiwat_change_pending = 0;
  if (!hiwat_is_interesting()) {
    return;
  }

  m_last_hiwat_report = time(0);
  SPLOG(
      "%s: BioGraph using significantly more memory than its previous maximum. "
      "Raised to %s from %s",
      m_type.c_str(), size_str(m_tot_alloc).c_str(), size_str(m_hiwat).c_str());
  m_hiwat = m_tot_alloc;
  log_detail_usage();
  on_new_hiwat(m_hiwat);
}

void track_mem::tracked_type::note_allocation(size_t nbytes, const std::string& description) {
  std::lock_guard<std::mutex> l(m_mu);
  m_tot_alloc += nbytes;
  if (!m_hiwat_change_pending && hiwat_is_interesting()) {
    m_hiwat_change_pending = time(0);
  }
}

void track_mem::tracked_type::note_deallocation(size_t nbytes, const std::string& description) {
  std::lock_guard<std::mutex> l(m_mu);

  if (hiwat_is_interesting()) {
    time_t now = time(0);
    if (m_last_hiwat_report + k_hiwat_report_interval < now) {
      report_hiwat();
    }
  }

  CHECK_GE(m_tot_alloc, nbytes) << "Attempted to deallocate more " << m_type
                                << " than was allocated, allocated: " << m_tot_alloc
                                << " Attempted to deallocate: " << nbytes;
  m_tot_alloc -= nbytes;
}

void track_mem::tracked_type::reset() {
  std::lock_guard<std::mutex> l(m_mu);
  if (m_hiwat > m_tot_alloc) {
    SPLOG_P(LOG_WARNING, "%s: Resetting high water mark from %s to %s", m_type.c_str(),
            size_str(m_hiwat).c_str(), size_str(m_tot_alloc).c_str());
    m_hiwat = m_tot_alloc;
    m_hiwat_change_pending = m_last_hiwat_report = 0;
  }
}

void track_mem::reset_stats() {
  for (auto& type : k_types) {
    type->reset();
  }

  if (g_reset_stats_hook) {
    g_reset_stats_hook();
  }
}

void track_mem::tracked_type::log_overview() const {
  std::lock_guard<std::mutex> l(m_mu);
  SPLOG("Type %s: %s used, high water %s", m_type.c_str(), size_str(m_tot_alloc).c_str(),
        size_str(m_hiwat).c_str());
}

void track_mem::log_usage() {
  SPLOG("Allocation stats:");
  for (auto& type : k_types) {
    type->log_detail_usage();
    type->log_overview();
  }
}

track_mem::malloc_hook_t track_mem::get_malloc_new_hook() {
  CHECK(!g_has_malloc_hook) << "Malloc should only be hooked once";
  g_has_malloc_hook = true;
  return malloc_new_hook;
}

namespace {

const std::string k_part_indicator = "-part-";

std::string description_for_path(const std::string& path) {
  auto pos = path.find(k_part_indicator);
  if (pos == std::string::npos) {
    return path;
  }

  // Path is of the form X-part-Y; change to X-part-*
  return path.substr(0, pos + k_part_indicator.size()) + "*";
}

}  // namespace

std::unordered_map<std::string, size_t> track_mem::mmap_tracker::get_detail_usage() const {
  std::unordered_map<std::string, size_t> dedup;
  {
    std::lock_guard<std::mutex> l(m_mu);

    // Deduplicate mmaps; there may be multiple maps of the same file.
    // Only use the largest of them.
    for (const auto& m : m_mmaps) {
      auto& d = dedup[m.second.description];
      d = std::max(d, m.second.size);
    }
  }

  // Group and sum by path patterns.
  std::unordered_map<std::string, size_t> result;
  for (const auto& d : dedup) {
    result[description_for_path(d.first)] += d.second;
  }
  return result;
}

void track_mem::tracked_type::log_detail_usage(bool force_logging) const {
  std::unordered_map<std::string, size_t> detail = get_detail_usage();
  std::multimap<size_t, std::string, std::greater<size_t>> top;
  for (const auto& d : detail) {
    top.insert(std::make_pair(d.second, d.first));
  }

  size_t count = 0;
  size_t other = 0;
  for (const auto& out : top) {
    if (count < k_top_interesting && out.first > k_min_interesting_size) {
      SPLOG_P(force_logging ? LOG_WARNING : LOG_DEBUG, "%s: %15s %s", m_type.c_str(),
              size_str(out.first).c_str(), out.second.c_str());
      ++count;
    } else {
      other += out.first;
    }
  }

  if (other > k_min_interesting_size) {
    SPLOG_P(force_logging ? LOG_WARNING : LOG_DEBUG, "%s: %15s %s", m_type.c_str(),
            size_str(other).c_str(), "Other");
  }
}

std::unordered_map<std::string, size_t> track_mem::malloc_tracker::get_detail_usage() const {
  std::unordered_map<std::string, size_t> result;

  std::lock_guard<std::mutex> l(m_mu);
  for (const auto& a : m_allocators) {
    result[a.first] = a.second->tot_used;
  }
  return result;
}

void track_mem::malloc_tracker::on_new_hiwat(size_t new_hiwat) {
  size_t configured_max = get_maximum_mem_bytes();
  if (new_hiwat > configured_max) {
#if NDEBUG
    // In release mode, only warn when we exceed the maximum.
    constexpr bool k_hard_limit = false;
#else
    constexpr bool k_hard_limit = true;
#endif
    std::string error_msg =
        printstring("New highwater %s %ld exceeds configured maximum %s %ld", size_str(new_hiwat).c_str(), new_hiwat,
                    size_str(configured_max).c_str(), configured_max);
    if (k_hard_limit && new_hiwat * 100 > configured_max * 101) {
      // Only error if it exceeds by more than 1%.
      SPLOG("ERROR: %s", error_msg.c_str());
      log_detail_usage(true /* force logging */);
      LOG(FATAL) << error_msg;
    } else {
      SPLOG("WARNING: %s", error_msg.c_str());
    }
  }
}

void track_mem::malloc_new_hook(const void* result, size_t size) {
  if (size <= k_max_untracked_bytes) {
    return;
  }

  g_malloc_tracker.malloc_new_hook(size);
}

void track_mem::malloc_tracker::malloc_new_hook(size_t size) {
  {
    std::lock_guard<std::mutex> l(g_expect_malloc_mu);
    auto it = g_expect_malloc.find(size);
    if (it != g_expect_malloc.end()) {
      g_expect_malloc.erase(it);
      return;
    }
  }

  static std::mutex g_untracked_report_mu;
  std::unique_lock<std::mutex> untracked_l(g_untracked_report_mu, std::try_to_lock);
  if (!untracked_l.owns_lock()) {
    // Only do one report at once, especially if generating a report
    // tries to call malloc again...
    return;
  }

  static size_t g_untracked_count = 0;
  size_t orig_count = g_untracked_count;
  ++g_untracked_count;
  // Report violations with exponential backoff.
  if (orig_count & g_untracked_count) {
    return;
  }

  static constexpr char k_untracked_msg[] =
      "Allocation exceeded size limit; allocation should be tracked using 'tracked_vector' or "
      "equivalent.";
  if (k_fatal_untracked || getenv("GTEST_TMP_DIR")) {
    CHECK_LE(size, k_max_untracked_bytes) << k_untracked_msg;
  }

  std::string stacktrace;
  DumpStackTraceToString(&stacktrace);

  std::cerr << "WARNING: " << k_untracked_msg << "\n";
  std::cerr << "Incident #" << g_untracked_count << " allocates " << size_str(size) << " (limit "
            << size_str(k_max_untracked_bytes) << ") at:\n"
            << stacktrace << "\n";
  std::cerr.flush();
}

std::string track_mem::size_str(size_t sz) {
  if (sz >= 100ULL * 1024 * 1024 * 1024) {
    // >100 gigs, report gigs
    return printstring("%ld G  ", sz / (1024 * 1024 * 1024));
  } else if (sz >= 10ULL * 1024 * 1024) {
    return printstring("%ld M ", sz / (1024 * 1024));
  } else if (sz >= 10ULL * 1024) {
    return printstring("%ld K", sz / 1024);
  } else {
    return printstring("%ld b", sz);
  }
}

void* track_mem::mmap_tracker::tracked_mmap(void* addr, size_t length, int prot, int flags, int fd,
                                            off_t offset, const std::string& description) {
  track_mem::initialize_if_needed();
  void* result = mmap(addr, length, prot, flags, fd, offset);
  if (result != MAP_FAILED) {
    {
      std::lock_guard<std::mutex> l(m_mu);
      entry new_entry;
      new_entry.description = description;
      new_entry.size = length;
      CHECK(m_mmaps.insert(std::make_pair(result, new_entry)).second) << description;
    }
    note_allocation(length, description);
  }

  return result;
}

int track_mem::mmap_tracker::tracked_munmap(void* addr, size_t length,
                                            const std::string& description) {
  note_deallocation(length, description);
  {
    std::lock_guard<std::mutex> l(m_mu);
    auto it = m_mmaps.find(addr);
    CHECK(it != m_mmaps.end()) << "Trying to munmap a region that was never mapped";
    CHECK_EQ(it->second.description, description);
    m_mmaps.erase(it);
  }
  return munmap(addr, length);
}

track_mem::malloc_tracker::entry* track_mem::malloc_tracker::get_entry(
    const std::string& description) {
  initialize_if_needed();
  std::lock_guard<std::mutex> l(m_mu);
  auto& e = m_allocators[description];
  if (e) {
    DCHECK_EQ(description, e->description);
  } else {
    e = make_unique<entry>();
    e->description = description;
  }
  inc_ref(e.get());
  return e.get();
}

namespace po = boost::program_options;

namespace {
// TODO(nils): There ought to be an easier way to do validation of
// program options.
struct max_mem_validator {
} g_max_mem_validator;

size_t g_max_mem_bytes = 0;

void validate(boost::any& v, const std::vector<std::string>& xs, max_mem_validator*, int) {
  po::validators::check_first_occurrence(v);
  std::string s = po::validators::get_single_string(xs);
  try {
    size_t max_mem_gb = boost::lexical_cast<size_t>(s);
    if (max_mem_gb < 1) {
      throw(io_exception("--max-mem must specify at least 1GB of RAM"));
    }
    size_t sys_gb = get_system_mem() / 1024 / 1024 / 1024;
    if (max_mem_gb > sys_gb) {
      throw(io_exception(printstring(
          "--max-mem must specify less than the total system memory of %ld GiB", sys_gb)));
    }
    g_max_mem_bytes = max_mem_gb * 1024 * 1024 * 1024;
  } catch (const boost::bad_lexical_cast&) {
    boost::throw_exception(po::invalid_option_value(s));
  }
}

}  // namespace

po::options_description track_mem_program_options() {
  po::options_description opts;
  opts.add_options()  //
      ("max-mem", po::value(&g_max_mem_validator),
       "Maximum memory to use, in GiB (48)")  //
      ;
  return opts;
}

size_t get_maximum_mem_bytes() {
  if (!g_max_mem_bytes) {
    g_max_mem_bytes = std::min<size_t>(get_system_mem(), 48ULL * 1024 * 1024 * 1024);
  }
  return g_max_mem_bytes;
}

void set_maximum_mem_bytes(size_t new_max) { g_max_mem_bytes = new_max; }
