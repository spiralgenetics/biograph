#include "modules/io/spiral_file_mmap.h"

#include "modules/io/mmap_buffer.h"
#include "vendor/minizip/ioapi_mem.h"
#include "vendor/minizip/unzip.h"
#include "vendor/minizip/zip.h"

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>

namespace {

void throw_if_unz_error(const std::string& what, int err) {
  if (err != UNZ_OK) {
    throw(io_exception(what + ": Got UNZIP error " + std::to_string(err)));
  }
}

void throw_if_zip_error(const std::string& what, int err) {
  if (err != ZIP_OK) {
    throw(io_exception(what + ": Got ZIP error " + std::to_string(err)));
  }
}

using autoclosing_fd = spiral_file_mmap_autoclosing_fd;

}  // namespace

class spiral_file_mmap_autoclosing_fd {
 public:
  spiral_file_mmap_autoclosing_fd() = delete;
  spiral_file_mmap_autoclosing_fd(int fd) : m_fd(fd) {}
  spiral_file_mmap_autoclosing_fd(const spiral_file_mmap_autoclosing_fd&) =
      delete;
  ~spiral_file_mmap_autoclosing_fd() {
    if (m_fd >= 0) ::close(m_fd);
  }

  int fd() { return m_fd; }

 private:
  int m_fd;
};

spiral_file_open_mmap::spiral_file_open_mmap(const std::string& filename,
                                             mmap_buffer::mode mode,
                                             const spiral_file_options& options)
    : spiral_file_open(options) {

  int fd;
  switch (mode) {
    case mmap_buffer::mode::read_write:
      m_mutable = true;
      fd = ::open(filename.c_str(), O_RDONLY);
      break;
    case mmap_buffer::mode::copy_on_write:
      m_mutable = true;
      fd = ::open(filename.c_str(), O_RDWR);
      break;
    default:
      m_mutable = false;
      fd = ::open(filename.c_str(), O_RDONLY);
      break;
  }
  m_fd = std::make_shared<autoclosing_fd>(fd);

  if (m_mutable) {
    m_mutable_mmap_buffer.emplace(new mmap_buffer(filename, mode));
    m_mmap_buffer = m_mutable_mmap_buffer;
  } else {
    m_mmap_buffer.emplace(new mmap_buffer(filename, mode));
  }
  import_zip_contents(filename);
}

bool spiral_file_open_mmap::is_mutable() const { return m_mutable; }

spiral_file_open_mmap::~spiral_file_open_mmap() {}

void spiral_file_open_mmap::import_zip_contents(const std::string& filename) {
  unzFile uf = unzOpen64(filename.c_str());
  if (!uf) {
    throw(io_exception("Could not open membuf zip file " + filename +
                       " for reading"));
  }
  int err = unzGoToFirstFile(uf);
  throw_if_unz_error("unzGoToFirstFile", err);

  do {
    // TODO(nils): Should we dynamically allocate filename_inzip?
    char filename_inzip[8192] = {0};
    unz_file_info64 file_info = {0};
    err = unzGetCurrentFileInfo64(uf, &file_info, filename_inzip,
                                  sizeof(filename_inzip), NULL, 0, NULL, 0);
    throw_if_unz_error("unzGetCurrentFileInfo64", err);

    // Make sure it's not compressed.
    if (file_info.compression_method != 0) {
      throw(io_exception("Compression not supported in spiral files(1)"));
    }

    if (file_info.compressed_size != file_info.uncompressed_size) {
      throw(io_exception("Compression not supported in spiral files(2)"));
    }

    int method;
    err = unzOpenCurrentFile2(uf, &method, 0, 1);
    throw_if_unz_error("unzOpenCurrentFile2", err);
    size_t pos = unzGetCurrentFileZStreamPos64(uf);
    CHECK_GT(pos, 0);

    CHECK(!m_paths.count(filename_inzip)) << filename_inzip;
    path_info& info = m_paths[filename_inzip];
    info.offset = pos;
    info.size = file_info.uncompressed_size;

    err = unzCloseCurrentFile(uf);
    throw_if_unz_error("unzCloseCurrentFile", err);
    err = unzGoToNextFile(uf);
  } while (err == UNZ_OK);
  if (err != UNZ_END_OF_LIST_OF_FILE) {
    throw_if_unz_error("unzGoToNetFile", err);
  }
  err = unzClose(uf);
  throw_if_unz_error("unzClose", err);
}

membuf spiral_file_open_mmap::get_path(
    const std::string& path, const spiral_file_options& options) const {
  auto it = m_paths.find(path);
  CHECK(it != m_paths.end()) << path;
  const path_info& i = it->second;

  std::lock_guard<std::mutex> l(i.mu);

  if (options.read_into_ram) {
    if (!i.in_mem_buf) {
      mutable_membuf in_ram =
          new owned_membuf(i.size, "spiral_file_mmap: " + path);
      char* bufptr = in_ram.mutable_data();
      size_t size_left = i.size;
      size_t offset = i.offset;
      while (size_left) {
        ssize_t nread =
            pread(m_fd->fd(), bufptr,
                  std::min<size_t>(size_left, 512ULL * 1024 * 1024), offset);
        if (nread <= 0) {
          throw(io_exception("Incomplete read into memory of " +
                             std::to_string(nread) + " bytes of " + path));
        }
        CHECK_LE(nread, size_left);
        size_left -= nread;
        offset += nread;
        bufptr += nread;
      }
      i.in_mem_buf.emplace(in_ram);
    }
    return *i.in_mem_buf;
  } else {
    return m_mmap_buffer->subbuf(i.offset, i.size);
  }
}

mutable_membuf spiral_file_open_mmap::get_mutable_path(
    const std::string& path, const spiral_file_options& options) {
  CHECK(path_is_present(path));
  path_info& i = m_paths[path];

  std::lock_guard<std::mutex> l(i.mu);

  // TODO(nils): Should we wire up a write-on-close here os we can
  // support options.read_into_ram?
  return m_mutable_mmap_buffer->subbuf(i.offset, i.size);
}

bool spiral_file_open_mmap::path_is_present(const std::string& path) const {
  return m_paths.count(path);
}

std::set<std::string> spiral_file_open_mmap::contents() const {
  std::set<std::string> result;
  for (const auto& path : m_paths) {
    result.insert(path.first);
  }
  return result;
}

struct spiral_file_mmap_internal {
  zipFile zf = nullptr;
  std::shared_ptr<autoclosing_fd> fd;
};

namespace {

void* fopen_func(void* opaque, const void* /* filename */, int /* mode */) {
  return opaque;
}

void* fopendisk_func(void* /* opaque */, void* /* stream */,
                     uint32_t /* number_disk */, int /* mode */) {
  return 0;
};

uint32_t fread_func(void* opaque, void* stream, void* buf, uint32_t size) {
  spiral_file_mmap_internal* internal =
      reinterpret_cast<spiral_file_mmap_internal*>(opaque);
  size_t tot_nread = 0;
  while (tot_nread != size) {
    int nread =
        read(internal->fd->fd(), ((char*)buf) + tot_nread, size - tot_nread);
    if (nread < 0) {
      return -1;
    }
    if (nread == 0) {
      return 0;
    }
    tot_nread += size;
    CHECK_LE(tot_nread, size);
  }

  return tot_nread;
}

uint32_t fwrite_func(void* opaque, void* stream, const void* buf,
                     uint32_t size) {
  spiral_file_mmap_internal* internal =
      reinterpret_cast<spiral_file_mmap_internal*>(opaque);
  size_t tot_nwrote = 0;
  while (tot_nwrote != size) {
    int nwrote =
        write(internal->fd->fd(), ((char*)buf) + tot_nwrote, size - tot_nwrote);
    if (nwrote <= 0) {
      return -1;
    }
    tot_nwrote += nwrote;
    CHECK_LE(tot_nwrote, size);
  }

  return tot_nwrote;
}

uint64_t ftell_func(void* opaque, void* stream) {
  spiral_file_mmap_internal* internal =
      reinterpret_cast<spiral_file_mmap_internal*>(opaque);
  return lseek(internal->fd->fd(), 0, SEEK_CUR);
}

long fseek_func(void* opaque, void* stream, uint64_t offset, int origin) {
  spiral_file_mmap_internal* internal =
      reinterpret_cast<spiral_file_mmap_internal*>(opaque);
  int lseek_origin = 0;
  switch (origin) {
    case ZLIB_FILEFUNC_SEEK_CUR:
      lseek_origin = SEEK_CUR;
      break;
    case ZLIB_FILEFUNC_SEEK_END:
      lseek_origin = SEEK_END;
      break;
    case ZLIB_FILEFUNC_SEEK_SET:
      lseek_origin = SEEK_SET;
      break;
    default:
      return -1;
  }

  if (lseek(internal->fd->fd(), offset, lseek_origin) < 0) {
    return -1;
  }
  return 0;
}

int fclose_func(void* opaque, void* stream) { return 0; }
int ferror_func(void* opaque, void* stream) { return 0; }

void fill_filefunc64(zlib_filefunc64_def* filefunc64) {
  filefunc64->zopen64_file = fopen_func;
  filefunc64->zopendisk64_file = fopendisk_func;
  filefunc64->zread_file = fread_func;
  filefunc64->zwrite_file = fwrite_func;
  filefunc64->ztell64_file = ftell_func;
  filefunc64->zseek64_file = fseek_func;
  filefunc64->zclose_file = fclose_func;
  filefunc64->zerror_file = ferror_func;
}

}  // namespace

// A mutable membuf that automatically writes the contents to a file
// after construction is done.  If we mmap the file directly while
// constructing a large structure, we run into a couple of problems:
//
// 1) the OS will often try to flush out dirty buffers to disk before
// we're done writing them.  If we have a 30 GB structure that we
// construct with random access and the OS is only willing to keep 10
// GB of dirty pages, it will spend all its time rewriting the
// structure to disk over and over again.
//
// 2) Under at least Ubuntu's Trusty, we can't use huge TLB pages when
// mapping files.  If we're constructing something with random access,
// we'll spend a lot of time pulling TLB entries from main memory for
// each of the 4k pages.  TODO(nils): It seems that we automatically
// get huge TLB pages... sometimes?  when using an anonymous MMAP.  It
// would be better to understand under what circumstances this
// happens.
//
// We also use mmap to generate an anonymous mapping if allocating
// large regions.  This is much faster than calloc since we're
// guaranteed that the pages are already zero and don't have to spend
// time filling them with zeros ourself.
class file_writing_membuf : public mutable_membuf_impl {
 public:
  // Maximum size to use calloc up to instead of mmapping our own
  // anonymous pages.
  const size_t k_max_malloc_size = 1024 * 1024 * 4;  // 4 MB

  file_writing_membuf(std::shared_ptr<autoclosing_fd> fd, size_t offset,
                      size_t size, const std::string& description)
      : m_fd(fd), m_offset(offset) {
    m_owned = mutable_membuf(new owned_membuf(size, description));
  }

  ~file_writing_membuf() {
    size_t to_write = m_owned.size();
    size_t offset = m_offset;
    const char* buf = m_owned.data();

    while (to_write) {
      // Write a maximum of 64 megs per write call.
      size_t to_write_this_chunk = std::min<size_t>(to_write, 64 * 1024 * 1024);
      ssize_t nwrote = pwrite(m_fd->fd(), buf, to_write_this_chunk, offset);
      CHECK_GT(nwrote, 0) << strerror(errno);
      CHECK_LE(nwrote, to_write_this_chunk);
      buf += nwrote;
      offset += nwrote;
      to_write -= nwrote;
    }
  }

  char* mutable_data() override { return m_owned.mutable_data(); }
  size_t size() override { return m_owned.size(); }

 private:
  std::shared_ptr<autoclosing_fd> m_fd;
  mutable_membuf m_owned;
  size_t m_offset = 0;
};

spiral_file_create_mmap::spiral_file_create_mmap(
    const std::string& filename, const spiral_file_options& options)
    : spiral_file_create(options), m_filename(filename) {
  m_internal.reset(new spiral_file_mmap_internal);

  zlib_filefunc64_def filefunc64 = {0};
  fill_filefunc64(&filefunc64);
  filefunc64.opaque = m_internal.get();

  int fd = open(filename.c_str(), O_CREAT | O_RDWR | O_TRUNC, 0666);
  if (fd < 0) {
    throw(io_exception("Could not open zip for writing: " + filename + ": " +
                       strerror(errno)));
  }
  m_internal->fd.reset(new autoclosing_fd(fd));

  m_internal->zf = zipOpen2_64("", 0, 0, &filefunc64);
  if (!m_internal->zf) {
    throw(io_exception("zipOpen2_64"));
  }
}

void spiral_file_create_mmap::create_path_contents(
    const std::string& path, const membuf& contents,
    const spiral_file_options& options) {
  CHECK(m_internal);
  zip_fileinfo file_info = {0};
  int err = zipOpenNewFileInZip2_64(
      m_internal->zf, path.c_str(), &file_info, NULL, 0, NULL,
      0, /* extrafields */
      0 /* comment */, 0 /* method: no compression */,
      Z_NO_COMPRESSION /* level */, 0 /* raw */, 1 /* zip64 */);
  throw_if_zip_error("zipOpenNewFileInZip64", err);

  membuf contents_left = contents;

  while (contents_left.size() > 0) {
    // zip library takes 32 bit size, so only write in 512 MB chunks.
    size_t to_write = std::min<size_t>(contents_left.size(), 512 * 1024 * 1024);
    CHECK_LT(to_write, 0x8FFFFFFF);
    err = zipWriteInFileInZip(m_internal->zf, contents_left.data(), to_write);
    throw_if_zip_error("zipWriteFInileInZip", err);
    contents_left =
        contents_left.subbuf(to_write, contents_left.size() - to_write);
  }
}

mutable_membuf spiral_file_create_mmap::create_path(
    const std::string& path, size_t part_size,
    const spiral_file_options& options) {
  zip_fileinfo file_info = {0};
  int err = zipOpenNewFileInZip2_64(
      m_internal->zf, path.c_str(), &file_info, NULL, 0, NULL,
      0, /* extrafields */
      0 /* comment */, 0 /* method: no compression */,
      Z_NO_COMPRESSION /* level */, 1 /* raw */, 1 /* zip64 */);
  throw_if_zip_error("zipOpenNewFileInZip64", err);

  err = zipFlush(m_internal->zf);
  throw_if_zip_error("zipFlush", err);

  off_t part_offset = lseek(m_internal->fd->fd(), 0, SEEK_END);
  if (part_offset <= 0) {
    throw(io_exception("lseek to end of zip: " + std::string(strerror(errno))));
  }
  CHECK_GT(part_offset, 0);
  size_t new_file_size = part_offset + part_size;
  if (ftruncate(m_internal->fd->fd(), new_file_size) < 0) {
    throw(io_exception("ftruncate to extend zip: " +
                       std::string(strerror(errno))));
  }

  size_t end_offset = lseek(m_internal->fd->fd(), 0, SEEK_END);
  CHECK_EQ(new_file_size, end_offset);

  err = zipPretendWriteInZip64(m_internal->zf, part_size);
  throw_if_zip_error("zipPretendWriteInZip64", err);

  err = zipCloseFileInZipRaw64(m_internal->zf, part_size,
                               0 /* don't bother to fill crc */);
  throw_if_zip_error("zipCloseFileInZipRaw64", err);

  err = zipFlush(m_internal->zf);
  throw_if_zip_error("zipFlush(2)", err);

  if (options.delayed_write || part_size < options.small_object_threshold) {
    return mutable_membuf(new file_writing_membuf(
        m_internal->fd, part_offset, part_size, "spiral_file: " + path));
  } else {
    mutable_membuf whole_file(
        new mmap_buffer(m_filename, mmap_buffer::mode::read_write));
    return whole_file.subbuf(part_offset, part_size);
  }
}

size_t spiral_file_create_mmap::close() {
  CHECK(m_internal);

  int err = zipClose_64(m_internal->zf, 0);
  throw_if_zip_error("zipClose_64", err);

  size_t end_offset = lseek(m_internal->fd->fd(), 0, SEEK_END);
  CHECK_GT(end_offset, 0);
  m_internal.reset();
  return end_offset;
}

spiral_file_create_mmap::~spiral_file_create_mmap() {
  if (m_internal) {
    close();
  }
}
