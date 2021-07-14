#ifndef _FILTER_H
#define _FILTER_H

#include "htslib/sam.h"
#include "report.h"
#include <stdint.h>

/* Filter counter structure */
struct filter_counter {
    uint8_t min_qual;
    uint16_t filter_incl;
    uint16_t filter_excl;

    uint64_t r_filtered;
    uint64_t r_unfiltered;
};
typedef struct filter_counter filter_counter_t;

/* Is record filtered by qual? */
#define filter_test_qual(qual, min_qual)    \
    ((qual) < (min_qual))

/* Is record filtered by flag? */
#define filter_test_flag(flag, filter_incl, filter_excl)    \
    (((flag) & (filter_incl)) != (filter_incl) ||           \
     ((flag) & (filter_excl)) != 0)

filter_counter_t *filter_counter_init(uint8_t min_qual, uint16_t filter_incl, uint16_t filter_excl);
void filter_counter_destroy(filter_counter_t *fc);

/* Filter counter calculation */
void filter_counter_process_record(bam1_t *rec, filter_counter_t *fc);
void filter_counter_report(report_t *report, filter_counter_t *fc);

#endif /* _FILTER_H */
