#include "base/command_line.h"

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <sys/types.h>
#include <unistd.h>
#include <fstream>

using namespace testing;

TEST(command_line_test, setproctitle) {
  setproctitle("Here is my new process title!");

  int mypid = getpid();
  
  std::string cmd = "ps ww " + std::to_string(mypid);
  FILE *f = popen(cmd.c_str(), "r");
  char buf[1024];
  std::string ps_output;
  while (fgets(buf, 1024, f)) {
    ps_output += buf;
  }
  fclose(f);

  EXPECT_THAT(ps_output, HasSubstr("Here is my new process title!"));

  // Make sure our original argv and environment are still entact.
  EXPECT_THAT(getenv("PATH"), HasSubstr(":"));
  EXPECT_THAT(original_program_args()[0], HasSubstr("command_line_test"));
}

TEST(command_line_test, setproctitle_long) {
  std::string long_title = "Here is my new long process title!";
  
  // 1 MB
  while (long_title.size() < 1024*1024*1024) {
    long_title = long_title + long_title;
  }

  setproctitle(long_title);

  int mypid = getpid();
  
  std::string cmd = "ps ww " + std::to_string(mypid);
  FILE *f = popen(cmd.c_str(), "r");
  char buf[1024];
  std::string ps_output;
  while (fgets(buf, 1024, f)) {
    ps_output += buf;
  }
  fclose(f);

  EXPECT_THAT(ps_output, HasSubstr("Here is my new long process title!"));

  // Make sure our original argv and environment are still entact.
  EXPECT_THAT(getenv("PATH"), HasSubstr(":"));
  EXPECT_THAT(original_program_args()[0], HasSubstr("command_line_test"));
}
