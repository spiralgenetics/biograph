#include "modules/io/spiral_file_mem.h"

spiral_file_open_mem::spiral_file_open_mem(const spiral_file_mem_storage& storage)
    : spiral_file_open(m_spiral_file_options), m_storage(storage) {}

spiral_file_open_mem::~spiral_file_open_mem() {}

spiral_file_create_mem::spiral_file_create_mem() : spiral_file_create(m_spiral_file_options) {}
spiral_file_create_mem::~spiral_file_create_mem() {}

mutable_membuf spiral_file_create_mem::create_path(const std::string& path, size_t size,
                                                   const spiral_file_options& options) {
  CHECK(!m_storage.paths.count(path)) << path;
  mutable_membuf result = new owned_membuf(size, path);
  m_storage.paths[path] = result;
  return result;
}

spiral_file_mem_storage spiral_file_create_mem::close() { return std::move(m_storage); }

membuf spiral_file_open_mem::get_path(const std::string& path,
                                      const spiral_file_options& options) const {
  CHECK(path_is_present(path)) << path;
  auto it = m_storage.paths.find(path);
  CHECK(it != m_storage.paths.end());
  return it->second;
}
mutable_membuf spiral_file_open_mem::get_mutable_path(const std::string& path,
                                                      const spiral_file_options& options) {
  CHECK(path_is_present(path)) << path;
  return m_storage.paths[path];
}

std::set<std::string> spiral_file_open_mem::contents() const {
  std::set<std::string> result;
  for (const auto& p : m_storage.paths) {
    result.insert(p.first);
  }
  return result;
}

bool spiral_file_open_mem::path_is_present(const std::string& path) const {
  return m_storage.paths.count(path);
}
