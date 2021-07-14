#pragma once

#include "modules/test/coverage.h"

// Support for making sure required code paths are tested.  See
// coverage.h for details.

#include <iostream>
#include <gtest/gtest.h>

#ifdef NDEBUG

class scoped_test_coverage {
 public:
  scoped_test_coverage() {}
  std::set<std::string> marked() const {
    std::cerr << "WARNING: NDEBUG build; test coverage markers not enabled.\n";
    return std::set<std::string>();
  }

  std::set<std::string> missing(const char* module_name) const {
    std::cerr << "WARNING: NDEBUG build; coverage testing not available.";
    return std::set<std::string>();
  }
};

#else

class scoped_test_coverage {
 public:
  scoped_test_coverage() {
    CHECK(!coverage_internal::test_coverage_enabled);

    coverage_internal::test_coverage_enabled = true;
  }

  std::set<std::string> missing(const char* module_name) const {
    std::set<std::string> result;
    auto m = coverage_internal::get_coverage_map();

    unsigned coverage_in_module = 0;

    for (const auto& entry : m) {
      if (entry.first.module_name == module_name) {
        coverage_in_module++;
        if (!entry.second.ever_marked &&
            !(entry.second.marked && *entry.second.marked)) {
          result.insert(entry.first.file_name + ":" +
                        std::to_string(entry.first.line));
        }
      }
    }

    // Make sure we have some coverage markers that might be missing
    EXPECT_GT(coverage_in_module, 0)
        << module_name
        << " test coverage requested but no test coverage markers present.";
    return result;
  }

  ~scoped_test_coverage() {
    CHECK(coverage_internal::test_coverage_enabled);
    coverage_internal::test_coverage_enabled = false;
    coverage_internal::reset_test_coverage();
  }
};

#endif
