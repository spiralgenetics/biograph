#ifndef _REPORT_H
#define _REPORT_H

#include <stdio.h>
#include <stdlib.h>

#define REPORT_BUFFER_SIZE 4096
#define REPORT_FORMAT "%s: %s\n"

/* Report element structure */
struct report_element {
    char *key;                   /* Elem's key (not NULL) */
    char *value;                 /* Elem's value (not NULL) */
    struct report_element *next; /* Elem's next sibling (linked-list) */
};
typedef struct report_element report_element_t;

report_element_t *report_element_init(char *key, char *value);
void report_element_destroy(report_element_t *element);

/* Report document structure */
struct report {
    report_element_t *head_elem; /* Document's head element in LL */
    report_element_t *tail_elem; /* Document's tail element in LL */
};
typedef struct report report_t;

report_t *report_init();
void report_destroy(report_t *report);

/* Report functions */
void report_add_element(report_t *report, report_element_t *element);
void report_add_key_value(report_t *report, char *key, char *value);
void report_print(FILE *stream, report_t *report);

#endif /* _REPORT_H */
