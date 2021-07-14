#include "filter.h"
#include "err.h"
#include "logging.h"
#include "print.h"
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

/* Counter for preliminarily filtered reads */

/**
 * Create new *filter_counter_t and initialize values
 */
filter_counter_t *filter_counter_init(uint8_t min_qual, uint16_t filter_incl, uint16_t filter_excl)
{
    filter_counter_t *fc = calloc(1, sizeof(filter_counter_t));
    die_on_alloc_fail(fc);

    /* global record filter */
    fc->min_qual = min_qual;
    fc->filter_incl = filter_incl;
    fc->filter_excl = filter_excl;

    fc->r_filtered = 0;
    fc->r_unfiltered = 0;

    return fc;
}

/**
 * Free *fc from memory.
 */
void filter_counter_destroy(filter_counter_t *fc)
{
    free(fc);
}

/**
 * Process a bam1_t record for filter counter.
 */
void filter_counter_process_record(bam1_t *rec, filter_counter_t *fc)
{
    if (filter_test_qual(rec->core.qual, fc->min_qual) ||
        filter_test_flag(rec->core.flag, fc->filter_incl, fc->filter_excl))
    {
        ++fc->r_filtered;
    } else {
        ++fc->r_unfiltered;
    }
}

/**
 * Write metrics in fc to report.
 */
void filter_counter_report(report_t *report, filter_counter_t *fc)
{
    char *key_buffer = malloc(REPORT_BUFFER_SIZE * sizeof(char));
    die_on_alloc_fail(key_buffer);
    char *value_buffer = malloc(REPORT_BUFFER_SIZE * sizeof(char));
    die_on_alloc_fail(value_buffer);

    uint64_t r_total = fc->r_filtered + fc->r_unfiltered;

    copy_to_buffer(key_buffer, "Total_Reads", REPORT_BUFFER_SIZE);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", r_total);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_buffer, "Unfiltered_Reads", REPORT_BUFFER_SIZE);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", fc->r_unfiltered);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_buffer, "Unfiltered_Reads_Pct", REPORT_BUFFER_SIZE);
    print_pct(value_buffer, REPORT_BUFFER_SIZE, fc->r_unfiltered, r_total);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_buffer, "Filtered_Reads", REPORT_BUFFER_SIZE);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", fc->r_filtered);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_buffer, "Filtered_Reads_Pct", REPORT_BUFFER_SIZE);
    print_pct(value_buffer, REPORT_BUFFER_SIZE, fc->r_filtered, r_total);
    report_add_key_value(report, key_buffer, value_buffer);

    free(key_buffer);
    free(value_buffer);
}
