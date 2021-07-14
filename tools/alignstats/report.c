#include "report.h"
#include "err.h"
#include <stdbool.h>
#include <string.h>

/* Report element structure */

/**
 * Create and return new report_element_t* with key:value.
 */
report_element_t *report_element_init(char *key, char *value)
{
    /* key and value cannot be NULL or empty */
    if (key == NULL || *key == '\0' || value == NULL || *value == '\0') {
        return NULL;
    }

    char *new_key, *new_value;
    size_t len;

    /* Key */
    len = strlen(key);
    new_key = malloc((size_t)((len + 1) * sizeof(char)));
    die_on_alloc_fail(new_key);
    strncpy(new_key, key, len);
    new_key[len] = '\0';

    /* Value */
    len = strlen(value);
    new_value = malloc((size_t)((len + 1) * sizeof(char)));
    die_on_alloc_fail(new_value);
    strncpy(new_value, value, len);
    new_value[len] = '\0';

    /* Report element */
    report_element_t *element = calloc(1, sizeof(report_element_t));
    die_on_alloc_fail(element);

    element->key = new_key;
    element->value = new_value;
    element->next = NULL;

    return element;
}

/**
 * Free *report_element_t from memory.
 */
void report_element_destroy(report_element_t *element)
{
    if (element != NULL) {
        free(element->key);
        free(element->value);
        free(element);
    }
}

/* Report document structure */

/**
 * Create and return new *report_t.
 */
report_t *report_init()
{
    report_t *report = calloc(1, sizeof(report_t));
    die_on_alloc_fail(report);

    report->head_elem = report->tail_elem = NULL;

    return report;
}

/**
 * Free *report_t from memory.
 */
void report_destroy(report_t *report)
{
    report_element_t *curr_elem, *next_elem;

    if (report != NULL) {
        curr_elem = report->head_elem;

        /* free each element */
        while (curr_elem != NULL) {
            next_elem = curr_elem->next;
            report_element_destroy(curr_elem);
            curr_elem = next_elem;
        }

        free(report);
    }
}

/* Report functions */

/**
 * Add element to report.
 */
void report_add_element(report_t *report, report_element_t *element)
{

    if (report->tail_elem == NULL) { /* Empty report */
        report->tail_elem = report->head_elem = element;
    } else {                         /* Non-empty report */
        report->tail_elem = report->tail_elem->next = element;
    }
}

/**
 * Add a new element with key:value to report.
 */
void report_add_key_value(report_t *report, char *key, char *value)
{
    report_add_element(report, report_element_init(key, value));
}

/**
 * Print report elements to stream.
 */
void report_print(FILE *stream, report_t *report)
{
    report_element_t *curr_elem = report->head_elem;

    /* for each element print key, value */
    while (curr_elem != NULL) {
        fprintf(stream, REPORT_FORMAT, curr_elem->key, curr_elem->value);
        curr_elem = curr_elem->next;
    }
}
