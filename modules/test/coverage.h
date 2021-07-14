#pragma once

#include <map>
#include <mutex>
#include <set>
#include <string>

// To assist with making sure edge cases have good test coverage,
// coverage.h provides macros to mark certain parts of the code as
// required to be exercised by tests.
//
// If "NDEBUG" is defined, these macros are compiled away completely.
// Otherwise, if coverage testing isn't enabled, only a single
// variable read is done at runtime so it has minimal overhead.
//
// To enable coverage testing, put DECLARE_TEST_COVERAGE(my_module) at
// the beginning of the module source file.  Then, for each section
// you want to ensure is covered by the test, call the
// "NOTE_TEST_COVERAGE" or "NOTE_TEST_COVERAGE_IF" macros to note
// blocks which must be exercised by test.
//
// To make sure all the noted sections are reached, you can do
// something like the following in the test:
//
// #include "modules/test/test_coverage.h"
//
// scoped_test_coverage cov;
// while (!cov.missing("my_module").empty()) {
//    std::cerr << "Missing coverage: "
//              << ::testing::PrintToString(cov.missing("my_module")) << "\n";
//    // Run my_module with random parameters to try to exercise
//    // all the marked blocks.
// }
//
// The strings returned by "missing" will contain the file name
// and the line number of the markers that were not exercised.

#ifndef NDEBUG

// At the beginning of a file that needs to note whether sections are
// covered by tests, call this to enable checking.
#define DECLARE_TEST_COVERAGE(MODULE)                              \
  struct test_coverage_##MODULE {                                  \
    static constexpr const char* module_name() { return #MODULE; } \
    static constexpr const char* file_name() { return __FILE__; }  \
  };

// Notes that this block should be exercised by a test when CONDITION
// is true.
#define NOTE_TEST_COVERAGE_IF(MODULE, CONDITION)                           \
  do {                                                                     \
    if (coverage_internal::test_coverage_enabled) {                        \
      static coverage_internal::test_coverage_tmpl<test_coverage_##MODULE, \
                                                   __LINE__>               \
          cov;                                                             \
      if (!cov.marked() && (CONDITION)) {                                  \
        cov.mark();                                                        \
      }                                                                    \
    }                                                                      \
  } while (0)

// Notes that this block should be exercised by a test.
#define NOTE_TEST_COVERAGE(MODULE) NOTE_TEST_COVERAGE_IF(MODULE, true)

namespace coverage_internal {

// The following should not be used directly.
extern bool test_coverage_enabled;

class test_coverage_registerer {
 public:
  test_coverage_registerer(const char* module_name, const char* file_name,
                           int line);
  static void register_marked(const char* module_name, const char* file_name,
                              int line, bool* marked);
};

template <typename T, int line>
class test_coverage_tmpl {
 public:
  test_coverage_tmpl() {
    registerer.register_marked(T::module_name(), T::file_name(), line,
                               &m_marked);
  }

  bool marked() const { return m_marked; }

  void mark() { m_marked = true; }

 private:
  bool m_marked = false;

  static test_coverage_registerer registerer;
};

template <typename T, int line>
test_coverage_registerer test_coverage_tmpl<T, line>::registerer(
    T::module_name(), T::file_name(), line);

struct coverage_entry {
  std::string module_name;
  std::string file_name;
  int line;

  bool operator<(const coverage_entry& entry) const {
    if (module_name != entry.module_name) {
      return module_name < entry.module_name;
    }
    if (file_name != entry.file_name) {
      return file_name < entry.file_name;
    }
    return line < entry.line;
  }
};

struct coverage_info {
  bool* marked;
  bool ever_marked;
};

std::map<coverage_entry, coverage_info> get_coverage_map();

void reset_test_coverage();

}  // namespace coverage_internal

#else

#define DECLARE_TEST_COVERAGE(MODULE) \
  static_assert(true,"")
#define NOTE_TEST_COVERAGE(MODULE) \
  static_assert(true,"")
#define NOTE_TEST_COVERAGE_IF(MODULE, CONDITION) \
  static_assert(true,"")

#endif  // NDEBUG
