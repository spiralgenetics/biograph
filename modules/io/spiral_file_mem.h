#pragma once

#include "modules/io/spiral_file.h"

struct spiral_file_mem_storage {
  std::map<std::string, mutable_membuf> paths;
};

class spiral_file_open_mem : public spiral_file_open {
 public:
  spiral_file_open_mem(const spiral_file_mem_storage &storage);
  ~spiral_file_open_mem();

  bool is_mutable() const override { return true; }

 private:
  membuf get_path(const std::string &path, const spiral_file_options &options) const override;
  mutable_membuf get_mutable_path(const std::string &path,
                                  const spiral_file_options &options) override;
  std::set<std::string> contents() const override;
  bool path_is_present(const std::string &path) const override;

  spiral_file_mem_storage m_storage;
  spiral_file_options m_spiral_file_options;
};

class spiral_file_create_mem : public spiral_file_create {
 public:
  spiral_file_create_mem();
  ~spiral_file_create_mem();

  spiral_file_mem_storage close();

 private:
  mutable_membuf create_path(const std::string &path, size_t size,
                             const spiral_file_options &options) override;

  spiral_file_mem_storage m_storage;
  spiral_file_options m_spiral_file_options;
};
