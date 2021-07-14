#include "modules/io/mmap_buffer.h"

#include "modules/io/io.h"
#include "modules/io/track_mem.h"
#include "modules/io/utils.h"

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <unistd.h>
#include <boost/format.hpp>
#include <cstdlib>
#include <cstring>

mmap_buffer::mmap_buffer(const std::string& path, size_t size) { open(path, size); }

mmap_buffer::mmap_buffer(const std::string& path, mode the_mode) { open(path, the_mode); }

mmap_buffer::~mmap_buffer() { close(); }

bool mmap_buffer::is_open() const { return m_file.is_open(); }

void mmap_buffer::open(const std::string& path, size_t size) {
  if (is_open()) {
    throw io_exception(boost::format("Trying to open already open mmap_buffer %1%") % path);
  }

  m_file.open(path, O_CREAT | O_RDWR | O_EXCL, 0644);

  if (ftruncate(m_file.get_fd(), size) < 0) {
    m_file.close();
    throw io_exception("Unable to truncate mmap file");
  }

  m_buffer = (char*)tracked_mmap(nullptr, size, PROT_READ | PROT_WRITE, MAP_SHARED, m_file.get_fd(),
                                 0, path);
  if (m_buffer == MAP_FAILED) {
    m_file.close();
    m_buffer = nullptr;
    throw io_exception(boost::format("Unable to mmap new file %1%. Make sure there is sufficient "
                                     "free memory and try again.") %
                       path);
  }

  m_size = size;
  m_mode = mode::read_write;
}

void mmap_buffer::truncate(size_t new_size) {
  if (ftruncate(m_file.get_fd(), new_size) < 0) {
    throw io_exception("Unable to truncate mmap file");
  }
  m_size = new_size;
}

void mmap_buffer::open(const std::string& path, mode the_mode) {
  if (is_open()) {
    throw io_exception(boost::format("Trying to open already open mmap_buffer %1%") % path);
  }

  int file_mode = O_RDONLY;
  int prot_mode = PROT_READ;
  int flags = MAP_SHARED;

  switch (the_mode) {
    case mode::read_only:
      break;
    case mode::read_write:
      file_mode = O_RDWR;
      prot_mode |= PROT_WRITE;
      break;
    case mode::copy_on_write:
      prot_mode |= PROT_WRITE;
      flags = MAP_PRIVATE;
      break;
    case mode::read_populate:
      flags = MAP_SHARED | MAP_POPULATE;
      break;
    default:
      throw io_exception(boost::format("Invalid mode in mmap_buffer::open %1%") %
                         static_cast<int>(the_mode));
  }
  m_mode = the_mode;
  m_file.open(path, file_mode);

  off_t seek_off = lseek(m_file.get_fd(), 0, SEEK_END);
  if (seek_off == -1) {
    m_file.close();
    throw io_exception("Unable to seek in mmap");
  }

  m_buffer = (char*)tracked_mmap(nullptr, seek_off, prot_mode, flags, m_file.get_fd(), 0, m_file.path());
  if (m_buffer == MAP_FAILED) {
    m_file.close();
    m_buffer = nullptr;
    throw io_exception("Unable to mmap. Make sure there is sufficient free memory and try again.");
  }

  m_size = (size_t)seek_off;
}

void mmap_buffer::sync() {
  if (!m_file.is_open()) {
    throw io_exception("Trying to sync unopened mmap file");
  }

  msync(m_buffer, m_size, MS_ASYNC);
}

void mmap_buffer::close() {
  if (!m_file.is_open()) {
    return;
  }

  tracked_munmap(m_buffer, m_size, m_file.path());
  m_file.close();
  m_buffer = nullptr;
  m_size = 0;
}
