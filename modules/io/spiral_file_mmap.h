#pragma once

#include "modules/io/mmap_buffer.h"
#include "modules/io/spiral_file.h"

#include <mutex>

class spiral_file_mmap_autoclosing_fd;

class spiral_file_open_mmap : public spiral_file_open {
 public:
  spiral_file_open_mmap(const std::string& filename, mmap_buffer::mode mode,
                        const spiral_file_options& options);
  spiral_file_open_mmap(const std::string& filename,
                        mmap_buffer::mode mode = mmap_buffer::mode::read_only)
      : spiral_file_open_mmap(filename, mode, m_spiral_file_opts) {}
  spiral_file_open_mmap(const std::string& filename, const spiral_file_options& options)
      : spiral_file_open_mmap(filename, mmap_buffer::mode::read_only, options) {}

  ~spiral_file_open_mmap() override;

  bool is_mutable() const override;

 private:
  struct path_info {
    size_t offset;
    size_t size;

    mutable std::mutex mu;
    mutable boost::optional<mutable_membuf> in_mem_buf;
  };

  spiral_file_options m_spiral_file_opts;
  membuf get_path(const std::string& path, const spiral_file_options& options) const override;
  mutable_membuf get_mutable_path(const std::string& path,
                                  const spiral_file_options& options) override;
  bool path_is_present(const std::string& path) const override;
  std::set<std::string> contents() const override;

  // Reads the zip format from m_buf into m_paths.
  void import_zip_contents(const std::string& filename);

  std::shared_ptr<spiral_file_mmap_autoclosing_fd> m_fd;
  boost::optional<mutable_membuf> m_mutable_mmap_buffer;
  boost::optional<membuf> m_mmap_buffer;
  std::map<std::string, path_info> m_paths;
  bool m_mutable = false;
};

struct spiral_file_mmap_internal;
class spiral_file_create_mmap : public spiral_file_create {
 public:
  spiral_file_create_mmap(const std::string& filename, const spiral_file_options& options);
  spiral_file_create_mmap(const std::string& filename)
      : spiral_file_create_mmap(filename, m_spiral_file_opts) {}
  ~spiral_file_create_mmap() override;

  // Returns the size of the resultant file.
  size_t close();

 private:
  friend class file_writing_membuf;
  friend struct spiral_file_mmap_internal;

  mutable_membuf create_path(const std::string& path, size_t size,
                             const spiral_file_options& options) override;
  void create_path_contents(const std::string& path, const membuf& contents,
                            const spiral_file_options& options) override;

  std::string m_filename;
  std::unique_ptr<spiral_file_mmap_internal> m_internal;
  spiral_file_options m_spiral_file_opts;
};
