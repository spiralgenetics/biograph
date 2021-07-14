// Base.h includes things that are needed in most spiral code in once
// place, so they don't have to be repeatedly specified.
#pragma once

// google log utilities: CHECK macro, etc.
#include "glog/logging.h"

// Runs one-time initialization common to spiral programs.
// This does the following:
//
// * Initializes GLOG so that CHECK failures get a stack trace
//
// * Saves the full command line so it can be used for e.g. output
//   file metadata
//
// * Parses any gflags present, and removes them from argc/argv.
//
// * Saves process title so it can be changed later with setproctitle.
void spiral_init(int* argc, char*** argv);

// Returns true if spiral_init has been called.
bool spiral_initted();

// Provide for functions that do magic that isn't fully understood by
// address/thread sanitizer.
#if defined(__clang__) || (defined (__GNUC__) && __GNUC__ > 4)
# define ATTRIBUTE_NO_SANITIZE_ADDRESS __attribute__((no_sanitize_address))
# define ATTRIBUTE_NO_SANITIZE_THREAD __attribute__((no_sanitize_thread))
#else
# define ATTRIBUTE_NO_SANITIZE_ADDRESS
# define ATTRIBUTE_NO_SANITIZE_THREAD
#endif

#if defined(__clang__) || (defined (__GNUC__) && __GNUC__ > 4)
#define WARN_UNUSED __attribute__((warn_unused_result))
#else
#define WARN_UNUSED
#endif
