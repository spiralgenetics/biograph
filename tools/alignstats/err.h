#ifndef _ERR_H
#define _ERR_H

#include <stdio.h>
#include <stdlib.h>

#define die_on_alloc_fail(obj)                                                 \
    if ((obj) == NULL) {                                                       \
        fprintf(stderr, __FILE__ ":%d memory allocation failed!\n", __LINE__); \
        exit(EXIT_FAILURE);                                                    \
    }

#endif /* _ERR_H */
