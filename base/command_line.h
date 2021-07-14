#pragma once

#include <string>
#include <vector>

// Saves the current command line, initializes proctitle, and replaces
// "environ" with a copy allocated from the heap. This should
// generally only be called from spiral_init.
void save_command_line(int argc, char** argv);

// Sets the current process title.
void setproctitle(const std::string& new_proctitle);

// Returns the full command line and args originally used to invoke
// this program, before any processing.
const std::vector<std::string>& original_program_args();
