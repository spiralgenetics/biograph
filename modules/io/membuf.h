#pragma once

#include <sys/types.h>
#include <map>
#include <memory>
#include <mutex>

#include "modules/io/progress.h"
#include "modules/io/track_mem.h"

class membuf_impl;
class mutable_membuf_impl;

// A read only reference to memory.  Copying just copies a reference.
class membuf {
 public:
  membuf() {}
  membuf(const std::shared_ptr<membuf_impl>& impl);
  membuf(const membuf& old) = default;
  membuf& operator=(const membuf& old) = default;

  // Takes ownership of a membuf impl:
  membuf(membuf_impl* impl) : membuf(std::shared_ptr<membuf_impl>(impl)) {}

  const char* data() const { return m_data; }
  size_t size() const { return m_size; }

  membuf subbuf(size_t offset, size_t length) const;

  // Returns contents as a std::string.  Do not use this on large
  // objects or in performance critical code as it makes a complete
  // copy of the contents of this membuf.
  std::string str() const { return std::string(m_data, m_size); }

 protected:
  membuf(const std::shared_ptr<membuf_impl>& impl, const char* data, size_t size);

 private:
  std::shared_ptr<membuf_impl> m_impl;
  const char* m_data = 0;
  size_t m_size = 0;
};

// A mutable reference to a  memory region.  Copying just copies a reference.
class mutable_membuf : public membuf {
 public:
  mutable_membuf() {}
  mutable_membuf(const std::shared_ptr<mutable_membuf_impl>& impl);
  mutable_membuf(const mutable_membuf& old) = default;
  mutable_membuf& operator=(const mutable_membuf& old) = default;

  // Takes ownership of a mutable_membuf impl:
  mutable_membuf(mutable_membuf_impl* impl)
      : mutable_membuf(std::shared_ptr<mutable_membuf_impl>(impl)) {}

  // Releases, if possible, memory required to hold the data in the
  // range [start, start + size).  The data may be zeroed out, but is
  // not guaranteed to be.
  void discard_region(char* start, size_t size);

  // Write to all the RAM in this membuf to ensure it's all writable,
  // and all mapped.  This can decrease memory fragmentation when
  // populating sparse arrays via random accesses.
  void populate_pages_for_write();

  char* mutable_data() const { return m_mutable_data; }

  mutable_membuf subbuf(size_t offset, size_t length) const;

 private:
  // Stride size for populate_pages_for_write. This number should not
  // exceed the page size.
  static constexpr size_t k_populate_stride_size = 4096;

  mutable_membuf(const std::shared_ptr<mutable_membuf_impl>& impl, char* data, size_t size);

  std::shared_ptr<mutable_membuf_impl> m_mutable_impl;
  char* m_mutable_data = 0;
};

// A provider of a read-only memory buffer.
class membuf_impl {
 public:
  virtual ~membuf_impl() = default;
  virtual const char* data() = 0;
  virtual size_t size() = 0;

 protected:
  membuf_impl();
};

// A provider of a writable memory buffer.
class mutable_membuf_impl : public membuf_impl {
 public:
  virtual ~mutable_membuf_impl() = default;
  // By default, discard does nothing.
  virtual void discard_region(char* start, size_t size) {}
  virtual char* mutable_data() = 0;
  const char* data() override;

 protected:
  mutable_membuf_impl();
};

// Refers to existing memory that's managed externally to this buffer.
// Do not use this in new code.
class borrowed_membuf : public membuf_impl {
 public:
  borrowed_membuf(const char* data, size_t size);
  const char* data() override;
  size_t size() override;

 private:
  const char* m_data = nullptr;
  size_t m_size = 0;
};

class borrowed_mutable_membuf : public borrowed_membuf, public mutable_membuf_impl {
 public:
  borrowed_mutable_membuf(char* data, size_t size);
  char* mutable_data() override;
  size_t size() override;

 private:
  char* m_mutable_data = nullptr;
};

// A membuf that owns its own data, and frees when done.
class owned_membuf : public mutable_membuf_impl {
 public:
  owned_membuf(size_t size, const std::string& description);
  ~owned_membuf() override;

  // Creates a membuf with the contents of the given string.  This
  // copies all the data, so should not be used in performance
  // critical code.
  static mutable_membuf from_str(const std::string& str, const std::string& description);

  // Creates an owned membuf with a copy of the given data.  Do not
  // use this in performance critical code, as it makes a copy of all
  // the data.
  owned_membuf(const char* data, size_t size, const std::string& description);

  void discard_region(char* start, size_t size) override;
  char* mutable_data() override;
  size_t size() override;

  // Buffers larger than this size will be allocated by calling mmap
  // instead of calloc.  Setting this too small may cause the system limit on
  // number of maps (vm.max_map_count) to be exceeded.
  static constexpr size_t k_mmap_threshold = 64ULL * 1024 * 1024;  // 64 MB

 private:
  track_mem::allocator<char> m_alloc;

  char* m_data = nullptr;
  size_t m_size = 0;
  size_t m_adjusted_size = 0;

  // True if this buffer was allocated via mmap, false if allocated via calloc.
  bool m_mmapped = false;
};

class membuf_cachelist {
 public:
  membuf_cachelist() = default;
  membuf_cachelist(const membuf_cachelist&) = default;
  membuf_cachelist(std::initializer_list<membuf_cachelist> membufs) {
    for (const auto& mb : membufs) {
      (*this) += mb;
    }
  }
  // Implicitly convertable from membuf:
  membuf_cachelist(membuf b);

  // Requires that these membufs (which might be mmaped) be cached in
  // RAM so it performs well even when hosted on devices that have
  // slow seeks.
  void cache_in_memory(progress_handler_t progress = null_progress_handler) const;

  // Returns true if this membuf appears to all be present in RAM.
  bool is_cached_in_memory() const;

  membuf_cachelist& operator+=(const membuf_cachelist& rhs);

 private:
  // When checking for whether a mapped file is cached in RAM, only
  // look at 1 byte every k_cache_stride_size bytes.  This number
  // should not exceed the page size.
  static constexpr size_t k_cache_stride_size = 4096;

  // Minimum speed we should be able to access every
  // k_cache_stride_size byte from a buffer, per gigabyte, if it's
  // already cached in RAM.
  static constexpr double k_gigabytes_per_second = 10;

  // Largest chunk to attempt to cache at once.
  static constexpr size_t k_cache_chunk_size = 1024 * 1024 * 128; // 128 MB

  // Membufs to cache.
  std::vector<membuf> m_membufs;
};

// TODO(nils): mmap_buffer should be an instance of membuf.
//
// TODO(nils): the following should use a membuf inside:
//  * packed_vector
//  * bitcount
