#include "base/base.h"
#include "base/command_line.h"

namespace {

bool g_spiral_initted = false;

}  // namespace

void spiral_init(int* argc, char*** argv) {
  CHECK(!g_spiral_initted) << "Must not call spiral_init more than once.";
  save_command_line(*argc, *argv);
  google::InitGoogleLogging((*argv)[0]);
  g_spiral_initted = true;
}

bool spiral_initted() {
  return g_spiral_initted;
}
