#include "pairstats.h"
#include "err.h"
#include "print.h"
#include <stdlib.h>
#include <string.h>

/**
 * Create new *pair_stats_metrics_t and initialize values.
 */
pair_stats_metrics_t *pair_stats_metrics_init()
{
    pair_stats_metrics_t *psm = calloc(1, sizeof(pair_stats_metrics_t));
    die_on_alloc_fail(psm);

    psm->cr_filter_mapped = BAM_FUNMAP |
                            BAM_FSECONDARY |
                            BAM_FQCFAIL |
                            BAM_FDUP;
    psm->cr_filter_improper_pair = psm->cr_filter_mapped |
                                   BAM_FPROPER_PAIR |
                                   BAM_FMUNMAP;

    return psm;
}

/**
 * Free *psm from memory.
 */
void pair_stats_metrics_destroy(pair_stats_metrics_t *psm)
{
    free(psm);
}

/**
 * Process a bam1_t record for pair statistics metrics.
 */
void pair_stats_process_record(bam1_t *rec, pair_stats_metrics_t *psm)
{
    /* Chimeric rate */
    if (!(rec->core.flag & psm->cr_filter_mapped)) {
        ++psm->cr_mapped;
        if (!(rec->core.flag & psm->cr_filter_improper_pair)) {
            ++psm->cr_improper_pair;
        }
    }

    /* Ignore unpaired and 1st reads for pair stats */
    if (!(rec->core.flag & BAM_FREAD1) && rec->core.flag & BAM_FPAIRED) {
        /* Count 2nd reads for total pairs */
        if (rec->core.flag & BAM_FREAD2) {
            ++psm->pairs_total;
        }

        if (rec->core.flag & BAM_FUNMAP) {
            if (rec->core.flag & BAM_FMUNMAP) {
                ++psm->pairs_unmapped; /* If neither read in pair is mapped */
            } else {
                ++psm->read1_mapped; /* Read 1 mapped, read 2 unmapped */
            }
        } else {
            if (rec->core.flag & BAM_FMUNMAP) {
                ++psm->read2_mapped; /* Read 2 mapped, read 1 unmapped */
            } else {
                ++psm->pairs_mapped; /* If both reads in pair are mapped */

                /* Mapped to same chromosome */
                if (rec->core.tid == rec->core.mtid) {
                    ++psm->pairs_mapped_same_chr;
                }
            }
        }
    }
}

/**
 * Write metrics in psm to report.
 */
void pair_stats_report(report_t *report, pair_stats_metrics_t *psm)
{
    char *key_buffer = malloc(REPORT_BUFFER_SIZE * sizeof(char));
    die_on_alloc_fail(key_buffer);
    char *value_buffer = malloc(REPORT_BUFFER_SIZE * sizeof(char));
    die_on_alloc_fail(value_buffer);

    copy_to_buffer(key_buffer, "Total_Pairs", REPORT_BUFFER_SIZE);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", psm->pairs_mapped);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_buffer, "Total_Same_Chr_Pairs", REPORT_BUFFER_SIZE);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", psm->pairs_mapped_same_chr);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_buffer, "Total_Same_Chr_Pairs_Pct", REPORT_BUFFER_SIZE);
    print_pct(value_buffer, REPORT_BUFFER_SIZE, psm->pairs_mapped_same_chr, psm->pairs_total);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_buffer, "Unpaired_Reads", REPORT_BUFFER_SIZE);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", psm->pairs_unmapped);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_buffer, "Unpaired_Reads_Pct", REPORT_BUFFER_SIZE);
    print_pct(value_buffer, REPORT_BUFFER_SIZE, psm->pairs_unmapped, psm->pairs_total);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_buffer, "R1_Unpaired_Reads", REPORT_BUFFER_SIZE);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", psm->read1_mapped);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_buffer, "R1_Unpaired_Reads_Pct", REPORT_BUFFER_SIZE);
    print_pct(value_buffer, REPORT_BUFFER_SIZE, psm->read1_mapped, psm->pairs_total);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_buffer, "R2_Unpaired_Reads", REPORT_BUFFER_SIZE);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", psm->read2_mapped);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_buffer, "R2_Unpaired_Reads_Pct", REPORT_BUFFER_SIZE);
    print_pct(value_buffer, REPORT_BUFFER_SIZE, psm->read2_mapped, psm->pairs_total);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_buffer, "Chimeric_Rate", REPORT_BUFFER_SIZE);
    print_pct(value_buffer, REPORT_BUFFER_SIZE, psm->cr_improper_pair, psm->cr_mapped);
    report_add_key_value(report, key_buffer, value_buffer);

    free(key_buffer);
    free(value_buffer);
}
