#ifndef _LOGGING_H
#define _LOGGING_H

#include <stdio.h>
#include <stdlib.h>

/**
 * Simple macros for logging. Argument format is the same as that of fprintf().
 */

/**
 * Log level in increasing order of verbosity
 */
#define LL_NONE  0
#define LL_ERROR 1
#define LL_WARN  2
#define LL_INFO  3
#define LL_DEBUG 4

#define LOGGING_NAME  "AlignStats"
#define LOGGING_DEST  stderr
#define LOGGING_LEVEL LL_INFO
/*
#define LOGGING_LEVEL LL_DEBUG
*/
#define LOGGING_NEWLINE fputc('\n', LOGGING_DEST)

#define log_error(...)                                                         \
    do {                                                                       \
        if (LOGGING_LEVEL >= LL_ERROR) {                                       \
            fprintf(LOGGING_DEST, "[" LOGGING_NAME " ERROR] " __VA_ARGS__);    \
            LOGGING_NEWLINE;                                                   \
        }                                                                      \
    } while (0)

#define log_error_r(...)                                                       \
    do {                                                                       \
        if (LOGGING_LEVEL >= LL_ERROR) {                                       \
            fprintf(LOGGING_DEST, __VA_ARGS__);                                \
        }                                                                      \
    } while (0)

#define log_warning(...)                                                       \
    do {                                                                       \
        if (LOGGING_LEVEL >= LL_WARN) {                                        \
            fprintf(LOGGING_DEST, "[" LOGGING_NAME " WARNING] " __VA_ARGS__);  \
            LOGGING_NEWLINE;                                                   \
        }                                                                      \
    } while (0)

#define log_warning_r(...)                                                     \
    do {                                                                       \
        if (LOGGING_LEVEL >= LL_WARN) {                                        \
            fprintf(LOGGING_DEST, __VA_ARGS__);                                \
        }                                                                      \
    } while (0)

#define log_info(...)                                                          \
    do {                                                                       \
        if (LOGGING_LEVEL >= LL_INFO) {                                        \
            fprintf(LOGGING_DEST, "[" LOGGING_NAME " INFO] " __VA_ARGS__);     \
            LOGGING_NEWLINE;                                                   \
        }                                                                      \
    } while (0)

#define log_info_r(...)                                                        \
    do {                                                                       \
        if (LOGGING_LEVEL >= LL_INFO) {                                        \
            fprintf(LOGGING_DEST, __VA_ARGS__);                                \
        }                                                                      \
    } while (0)

#define log_debug(...)                                                         \
    do {                                                                       \
        if (LOGGING_LEVEL >= LL_DEBUG) {                                       \
            fprintf(LOGGING_DEST, "[" LOGGING_NAME " DEBUG] " __VA_ARGS__);    \
            LOGGING_NEWLINE;                                                   \
        }                                                                      \
    } while (0)

#define log_debug_r(...)                                                       \
    do {                                                                       \
        if (LOGGING_LEVEL >= LL_DEBUG) {                                       \
            fprintf(LOGGING_DEST, __VA_ARGS__);                                \
        }                                                                      \
    } while (0)

#define log(...) fprintf(LOGGING_DEST, __VA_ARGS__)

#endif /* _LOGGING_H */
