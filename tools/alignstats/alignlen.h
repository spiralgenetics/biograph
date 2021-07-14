#ifndef _ALIGNLEN_H
#define _ALIGNLEN_H

#include "htslib/sam.h"
#include "readtype.h"
#include "report.h"
#include "treemap.h"
#include <stdint.h>

/* Alignment length metrics structure */
struct align_len_metrics {
    double mean;            /* Mean alignment length */
    uint64_t median;        /* Median alignment length */
    uint64_t mode;          /* Mode alignment length */
    tree_map_t *length_map; /* Treemap of alignment lengths */
};
typedef struct align_len_metrics align_len_metrics_t;

align_len_metrics_t *align_len_metrics_init();
void align_len_metrics_destroy(align_len_metrics_t *alm);

/* Alignment length metrics calculation */
void align_len_process_record(bam1_t *rec, align_len_metrics_t *alm);
void align_len_finalize(align_len_metrics_t *alm);
void align_len_report(report_t *report, align_len_metrics_t *alm, read_type_t rt);

#endif /* _ALIGNLEN_H */
