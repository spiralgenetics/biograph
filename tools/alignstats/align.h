#ifndef _ALIGN_H
#define _ALIGN_H

#include "htslib/sam.h"
#include "readtype.h"
#include "report.h"
#include <stdbool.h>
#include <stdint.h>

/* Alignment metrics structure */
struct align_metrics {
    /* Read metrics */
    uint64_t r_total;
    uint64_t r_aligned;
    uint64_t r_dup;
    uint64_t r_mapped;
    uint64_t r_unmapped;
    uint64_t r_soft_clipped;
    uint64_t r_exact_match;
    uint64_t r_mapq20;

    /* Base metrics */
    uint64_t b_total;
    uint64_t b_dup;
    uint64_t b_mapped;
    uint64_t b_unmapped;
    uint64_t b_aligned;
    uint64_t b_matched;
    uint64_t b_mismatched;
    uint64_t b_inserted;
    uint64_t b_deleted;
    uint64_t b_soft_clipped;
    uint64_t b_exact_match;
    uint64_t b_q20;

    /* Filter reads */
    uint16_t filter;
};
typedef struct align_metrics align_metrics_t;

align_metrics_t *align_metrics_init();
void align_metrics_destroy(align_metrics_t *am);

/* Alignment metrics calculation */
void align_process_record(bam1_t *rec, align_metrics_t *am, bool process_cigar);
void align_report(report_t *report, align_metrics_t *am, read_type_t rt);

#endif /* _ALIGN_H */
