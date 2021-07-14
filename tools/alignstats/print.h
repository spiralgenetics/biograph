#ifndef _PRINT_H
#define _PRINT_H

#include <stdint.h>
#include <stdio.h>

int print_pct(char *buffer, size_t num, uint64_t numerator, uint64_t denominator);
char *copy_to_buffer(char *dest, const char *src, size_t num);

#endif /* _PRINT_H */
