#include "tools/malloc_select.h"

#if TCMALLOC
// Workaround for https://sourceware.org/bugzilla/show_bug.cgi?id=20432
// Workaround suggested in https://github.com/jemalloc/jemalloc/issues/442
// It should be safe to remove these after glibc 2.25
extern "C" {
void __malloc_fork_lock_parent() {}
void __malloc_fork_unlock_parent() {}
void __malloc_fork_unlock_child() {}
};
#endif
