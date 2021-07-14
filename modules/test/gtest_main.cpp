#include <gtest/gtest.h>
#include <string>
#include "base/base.h"

#include <iostream>

int main(int argc, char** argv) {
  spiral_init(&argc, &argv);
  const char* tmpdir = getenv("TEST_TMPDIR");
  if (tmpdir) {
    setenv("TMPDIR", tmpdir, 1 /* overwrite */);
  }

  // Make sure standard output is line buffered when we're running
  // tests, even if we're capturing the output.
  std::cout << std::unitbuf;

  testing::InitGoogleTest(&argc, argv);

  return RUN_ALL_TESTS();
}
