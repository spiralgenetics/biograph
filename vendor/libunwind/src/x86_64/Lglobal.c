#define UNW_LOCAL_ONLY
#include "libunwind-config.h"
#include <libunwind.h>
#if defined(UNW_LOCAL_ONLY) && !defined(UNW_REMOTE_ONLY)
#include "Gglobal.c"
#endif
