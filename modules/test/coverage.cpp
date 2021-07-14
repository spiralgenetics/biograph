#include "modules/test/coverage.h"
#include "base/base.h"

#include <iostream>
#include <map>

#ifndef NDEBUG

namespace coverage_internal {

bool test_coverage_enabled = false;

namespace {

std::map<coverage_entry, coverage_info>* g_coverage = nullptr;

std::mutex g_mutex;

}  // namespace

test_coverage_registerer::test_coverage_registerer(const char* module_name,
                                                   const char* file_name,
                                                   int line) {
  std::lock_guard<std::mutex> l(g_mutex);
  if (!g_coverage) {
    g_coverage = new std::map<coverage_entry, coverage_info>();
  }
  CHECK(
      g_coverage
          ->insert(std::make_pair(coverage_entry{module_name, file_name, line},
                                  coverage_info{nullptr, false}))
          .second);
}

void test_coverage_registerer::register_marked(const char* module_name,
                                               const char* file_name, int line,
                                               bool* marked) {
  std::lock_guard<std::mutex> l(g_mutex);

  auto i = g_coverage->find(coverage_entry{module_name, file_name, line});
  CHECK(i != g_coverage->end());
  CHECK(!i->second.marked);
  i->second = coverage_info{marked, false};
}

void reset_test_coverage() {
  std::lock_guard<std::mutex> l(g_mutex);

  for (auto& marked : *g_coverage) {
    if (marked.second.marked) {
      if (*marked.second.marked) {
        marked.second.ever_marked = true;
      }
      *marked.second.marked = false;
    }
  }
}

std::map<coverage_entry, coverage_info> get_coverage_map() {
  std::lock_guard<std::mutex> l(g_mutex);
  return *g_coverage;
}

}  // namespace coverage_internal

#endif
