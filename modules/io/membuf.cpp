#include "modules/io/membuf.h"
#include "base/base.h"
#include "modules/io/io.h"
#include "modules/io/parallel.h"
#include "modules/io/track_mem.h"

#include <string.h>
#include <sys/mman.h>
#include <unistd.h>

namespace {

static constexpr bool k_huge_pages_enabled = false;

}  // namespace

membuf membuf::subbuf(size_t offset, size_t new_size) const {
  CHECK_LE(offset + new_size, m_size) << "Offset: " << offset << " New size: " << new_size;
  return membuf(m_impl, m_data + offset, new_size);
}

membuf::membuf(const std::shared_ptr<membuf_impl>& impl, const char* data, size_t size)
    : m_impl(impl), m_data(data), m_size(size) {}

membuf::membuf(const std::shared_ptr<membuf_impl>& impl)
    : m_impl(impl), m_data(impl->data()), m_size(impl->size()) {}

mutable_membuf mutable_membuf::subbuf(size_t offset, size_t new_size) const {
  CHECK_LE(offset + new_size, size()) << "Offset: " << offset << " New size: " << new_size;
  return mutable_membuf(m_mutable_impl, m_mutable_data + offset, new_size);
}

mutable_membuf::mutable_membuf(const std::shared_ptr<mutable_membuf_impl>& impl, char* data,
                               size_t size)
    : membuf(impl, data, size), m_mutable_impl(impl), m_mutable_data(data) {}

mutable_membuf::mutable_membuf(const std::shared_ptr<mutable_membuf_impl>& impl)
    : mutable_membuf(impl, impl->mutable_data(), impl->size()) {}

void mutable_membuf::discard_region(char* start, size_t dsize) {
  CHECK(data());
  CHECK_GE(start, data());
  CHECK_LE(start, data() + size());
  CHECK_LE(start + dsize, data() + size());
  return m_mutable_impl->discard_region(start, dsize);
}

constexpr size_t mutable_membuf::k_populate_stride_size;

void mutable_membuf::populate_pages_for_write() {
  size_t pos = 0;
  while (pos < size()) {
    mutable_data()[pos]++;
    __sync_synchronize();
    mutable_data()[pos]--;
    pos += k_populate_stride_size;
  }
}

membuf_impl::membuf_impl() {}

mutable_membuf_impl::mutable_membuf_impl() {}

const char* mutable_membuf_impl::data() { return mutable_data(); }

borrowed_membuf::borrowed_membuf(const char* data, size_t size) : m_data(data), m_size(size) {}
const char* borrowed_membuf::data() { return m_data; }
size_t borrowed_membuf::size() { return m_size; }

borrowed_mutable_membuf::borrowed_mutable_membuf(char* data, size_t size)
    : borrowed_membuf(data, size), m_mutable_data(data) {}
char* borrowed_mutable_membuf::mutable_data() { return m_mutable_data; }
size_t borrowed_mutable_membuf::size() { return borrowed_membuf::size(); }

owned_membuf::owned_membuf(size_t size, const std::string& description)
    : m_alloc((size < k_mmap_threshold) ? description : (description + "(mmap)")),
      m_size(size),
      m_adjusted_size(size) {
  if (size < k_mmap_threshold) {
    m_data = m_alloc.allocate(m_adjusted_size);
    if (!m_data) {
      track_mem::log_usage();
      throw(io_exception("Unable to allocate " + std::to_string(size) + " bytes for " +
                         description + ": " + strerror(errno)));
    }
    memset(m_data, 0, m_adjusted_size);
  } else {
    m_mmapped = true;

    constexpr size_t k_gb = 1024 * 1024 * 1024;

    if (size >= k_gb * 2 && k_huge_pages_enabled) {
      m_adjusted_size = size;
      m_adjusted_size += k_gb - 1;
      m_adjusted_size -= m_adjusted_size % k_gb;
// Use hugepages
#ifdef MAP_HUGE_SHIFT
      m_data = reinterpret_cast<char*>(
          mmap(nullptr /* don't specify base address */, m_adjusted_size, PROT_READ | PROT_WRITE,
               MAP_PRIVATE | MAP_ANONYMOUS | MAP_HUGETLB | (30 << MAP_HUGE_SHIFT), -1 /* fd */,
               0 /* offset */));
#else
      // MAP_HUGE_SHIFT is missing on Linux pre-3.8, so use the default instead.
      m_data = reinterpret_cast<char*>(
          mmap(nullptr /* don't specify base address */, m_adjusted_size, PROT_READ | PROT_WRITE,
               MAP_PRIVATE | MAP_ANONYMOUS | MAP_HUGETLB, -1 /* fd */, 0 /* offset */));
#endif
      static bool did_complain = false;
      if (m_data != MAP_FAILED) {
        did_complain = false;
        return;
      }

      if (!did_complain) {
        SPLOG("Unable to allocate huge pages for size %ld (to support %ld): %s", m_adjusted_size,
              size, strerror(errno));
        did_complain = true;
        track_mem::log_usage();
      }
      m_adjusted_size = size;
    }

    m_data = reinterpret_cast<char*>(mmap(nullptr /* don't specify base address */, size,
                                          PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS,
                                          -1 /* fd */, 0 /* offset */));
    if (m_data == MAP_FAILED) {
      track_mem::log_usage();
      throw(io_exception("Unable to allocate (via mmap) " + std::to_string(size) +
                         " bytes: " + strerror(errno)));
    }

    m_alloc.note_external_allocation(m_data, m_adjusted_size);
  }
}

mutable_membuf owned_membuf::from_str(const std::string& str, const std::string& description) {
  return mutable_membuf(std::make_shared<owned_membuf>(str.data(), str.size(), description));
}

owned_membuf::owned_membuf(const char* data, size_t size, const std::string& description)
    : owned_membuf(size, description) {
  memcpy(m_data, data, size);
}

owned_membuf::~owned_membuf() {
  CHECK(m_data);

  if (m_mmapped) {
    m_alloc.note_external_deallocation(m_data, m_adjusted_size);
    CHECK_EQ(0, munmap(m_data, m_adjusted_size))
        << ((void*)m_data) << ": " << m_adjusted_size << ": " << strerror(errno);
  } else {
    m_alloc.deallocate(m_data, m_adjusted_size);
  }
}

char* owned_membuf::mutable_data() {
  CHECK(m_data);
  return m_data;
}

size_t owned_membuf::size() { return m_size; }

void owned_membuf::discard_region(char* start, size_t size) {
  if (!size || !m_mmapped) {
    return;
  }

  CHECK_GT(start, m_data);
  CHECK_LT(start, m_data + m_size);
  CHECK_LE(start + size, m_data + m_size);

  int pagesize = sysconf(_SC_PAGESIZE);

  CHECK_GT(pagesize, 0) << "Couldn't determine page size: " << strerror(errno);

  size_t start_offset = reinterpret_cast<intptr_t>(start) % pagesize;
  if (start_offset) {
    size_t start_advance = pagesize - start_offset;

    if (size < start_advance) {
      return;
    }

    start += start_advance;
    size -= start_advance;
  }

  if (size < size_t(pagesize)) {
    return;
  }

  size_t end_trim = size % pagesize;
  size -= end_trim;

  if (madvise(start, size, MADV_DONTNEED) < 0) {
    SPLOG("Discard of region starting at %p for %ld bytes failed: %s", start, size,
          strerror(errno));
  }
}

constexpr size_t membuf_cachelist::k_cache_stride_size;
constexpr double membuf_cachelist::k_gigabytes_per_second;
constexpr size_t membuf_cachelist::k_cache_chunk_size;

membuf_cachelist::membuf_cachelist(membuf b) {
  for (size_t i = 0; i < b.size(); i += k_cache_chunk_size) {
    m_membufs.push_back(b.subbuf(i, std::min(b.size() - i, k_cache_chunk_size)));
  }
}

void membuf_cachelist::cache_in_memory(progress_handler_t progress) const {
  size_t sum = 0;

  parallel_for(  //
      0, m_membufs.size(),
      [this, &sum](size_t i, parallel_state& st) {
        size_t local_sum = 0;
        const auto* data = m_membufs[i].data();
        size_t size = m_membufs[i].size();
        for (size_t pos = 0; pos < size; pos += k_cache_stride_size) {
          local_sum += data[pos];
        }
        // Doesn't need to be thread safe because we don't catually care what
        // the sum is.
        sum += local_sum;
      },
      progress);

  // Make sure we don't optimize out the sum variable we uselessly
  // calculated.
  asm volatile("" ::"r"(sum));
}

bool membuf_cachelist::is_cached_in_memory() const {
  size_t tot_size = 0;
  for (const auto& mb : m_membufs) {
    tot_size += mb.size();
  }

  if (tot_size < k_cache_stride_size) {
    return true;
  }

  unsigned sum = 0;

  auto wait_time = std::chrono::nanoseconds(size_t(tot_size / k_gigabytes_per_second)) +
                   std::chrono::milliseconds(5);

  std::atomic<unsigned> taking_too_long{0};
  std::mutex mu;
  std::condition_variable cv;

  std::future<void> watchdog = std::async(std::launch::async, [&]() {
    std::unique_lock<std::mutex> l(mu);
    if (cv.wait_for(l, wait_time) == std::cv_status::timeout) {
      taking_too_long.store(1);
    }
  });

  size_t pos = 0;
  auto start_time = std::chrono::steady_clock::now();
  for (const auto& mb : m_membufs) {
    for (size_t pos = 0; pos < mb.size() && !taking_too_long.load(); pos += k_cache_stride_size) {
      sum += mb.data()[pos];
      pos += k_cache_stride_size;
    }
  }
  auto end_time = std::chrono::steady_clock::now();

  bool did_take_too_long = taking_too_long.load();
  cv.notify_all();

  if (did_take_too_long) {
    std::chrono::duration<double> duration = end_time - start_time;
    SPLOG(
        "Slow membuf read: was only able to read %ld bytes in %.8f ms "
        "(%.3f ms/gigabyte)",
        pos, duration.count() * 1000., duration.count() * 1000. * 1024 * 1024 * 1024 / pos);
  }

  // Make sure we don't optimize out the sum variable we uselessly
  // calculated.
  asm volatile("" ::"r"(sum));

  return !did_take_too_long;
}

membuf_cachelist& membuf_cachelist::operator+=(const membuf_cachelist& rhs) {
  m_membufs.insert(m_membufs.end(), rhs.m_membufs.begin(), rhs.m_membufs.end());
  return *this;
}
