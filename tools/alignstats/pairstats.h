#ifndef _PAIRSTATS_H
#define _PAIRSTATS_H

#include "htslib/sam.h"
#include "report.h"

/* Pair statistics metrics */
struct pair_stats_metrics {
    /* Pair mapping metrics */
    uint64_t pairs_total;
    uint64_t pairs_mapped;
    uint64_t pairs_mapped_same_chr;
    uint64_t read1_mapped;
    uint64_t read2_mapped;
    uint64_t pairs_unmapped;

    /* Chimeric rate */
    uint64_t cr_mapped;        /* chimeric rate mapped */
    uint64_t cr_improper_pair; /* chimeric rate mapped + improper pair */
    uint16_t cr_filter_mapped;
    uint16_t cr_filter_improper_pair;
};
typedef struct pair_stats_metrics pair_stats_metrics_t;

pair_stats_metrics_t *pair_stats_metrics_init();
void pair_stats_metrics_destroy(pair_stats_metrics_t *psm);

/* Pair statistics metrics calculation */
void pair_stats_process_record(bam1_t *rec, pair_stats_metrics_t *psm);
void pair_stats_report(report_t *report, pair_stats_metrics_t *psm);

#endif /* _PAIRSTATS_H */
