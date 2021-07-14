#include "base/command_line.h"

#include <iostream>
#include <string.h>
#include <sys/prctl.h>
#include <unistd.h>

namespace {

// A copy of the original argv.
std::vector<std::string> g_original_args;

// Pointer to the start of the buffer we can use to overwrite our
// process title.
char *g_proctitle_buffer = nullptr;
size_t g_proctitle_buffer_len = 0;

// Pointer to the original argv;
char **g_orig_argv = nullptr;

}  // namespace

void save_command_line(int argc, char **argv) {
  // Duplicate argv with one allocated from the heap.
  for (int i = 0; i < argc; i++) {
    g_original_args.push_back(argv[i]);
  }

  if (getenv("HEAPCHECK")) {
    return;
  }
  // Duplicate the original environment pointer with one allocated from the
  // HEAP.
  char **orig_environ = environ;
  size_t environ_count = 0;
  while (orig_environ[environ_count])
    environ_count++;

  char **new_environ = (char **)malloc(sizeof(char *) * (environ_count + 1));
  if (!new_environ) {
    return;
  }

  for (size_t i = 0; i < environ_count; ++i) {
    new_environ[i] = strdup(orig_environ[i]);
  }
  new_environ[environ_count] = nullptr;

  environ = new_environ;

  g_orig_argv = argv;
  g_proctitle_buffer = argv[0];

  char* proctitle_buffer_end;
  proctitle_buffer_end = argv[argc - 1];

  proctitle_buffer_end += strlen(proctitle_buffer_end);
  g_proctitle_buffer_len = proctitle_buffer_end - g_proctitle_buffer;
}

void setproctitle(const std::string &new_proctitle) {
  if (getenv("HEAPCHECK")) {
    return;
  }
  memset(g_proctitle_buffer, 0, g_proctitle_buffer_len);
  memcpy(g_proctitle_buffer, new_proctitle.data(),
         std::min<size_t>(new_proctitle.size(), g_proctitle_buffer_len - 1));

  g_orig_argv[0] = g_proctitle_buffer;
  g_orig_argv[1] = nullptr;
  prctl(PR_SET_NAME, g_proctitle_buffer); // this is required for 'top' to work
}

const std::vector<std::string>& original_program_args() {
  return g_original_args;
}
