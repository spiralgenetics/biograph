#include "align.h"
#include "err.h"
#include "print.h"
#include "logging.h"
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Alignment metrics structure */

/**
 * Create and return new *align_metrics_t.
 */
align_metrics_t *align_metrics_init()
{
    align_metrics_t *am = calloc(1, sizeof(align_metrics_t));
    die_on_alloc_fail(am);

    /* Secondary read or QC failure */
    am->filter = BAM_FSECONDARY | BAM_FQCFAIL;

    return am;
}

/**
 * Free *align_metrics_t from memory.
 */
void align_metrics_destroy(align_metrics_t *am)
{
    free(am);
}

/* Alignment metrics calculation */

/**
 * Process a bam1_t record for alignment metrics.
 */
void align_process_record(bam1_t *rec, align_metrics_t *am, bool process_cigar)
{
    char op;
    uint8_t /* *nm_tag,*/ *md_tag, *qual;
    int32_t pos, end_pos;
    uint32_t oplen, *cigar;
    uint64_t num_m, num_i, num_d, num_s, num_eq, num_x;
    uint64_t num_matches, num_mismatches, num_aligned;

    /* Total reads */
    ++am->r_total;
    am->b_total += rec->core.l_qseq;

    /* Unmapped reads */
    if (rec->core.flag & BAM_FUNMAP) {
        ++am->r_unmapped;
        am->b_unmapped += rec->core.l_qseq;
    } else {
        /* Mapped reads */
        ++am->r_mapped;
        am->b_mapped += rec->core.l_qseq;
        
        /* MAPQ20 reads */
        if (rec->core.qual >= 20) {
            ++am->r_mapq20;
        }
        
        /* Aligned reads: filter out reads with filtered flags */
        if (!(rec->core.flag & am->filter)) {
            ++am->r_aligned;
        }

        /* Duplicate reads */
        if (rec->core.flag & BAM_FDUP) {
            ++am->r_dup;
            am->b_dup += rec->core.l_qseq;
        }

        if (process_cigar) {
            qual = bam_get_qual(rec);
            cigar = bam_get_cigar(rec);
            pos = 0;
            num_m = num_i = num_d = num_s = num_eq = num_x = 0;

            /* Count M, I, D, S in CIGAR string (processCIGAR) */

            /* for each CIGAR operator */
            for (uint16_t i = 0; i < rec->core.n_cigar; ++i) {
                op = bam_cigar_op(cigar[i]);
                oplen = bam_cigar_oplen(cigar[i]);

                /* Sum the length of these operators */
                switch (op) {
                case BAM_CMATCH:
                    num_m += oplen;
                    break;
                case BAM_CEQUAL:
                    num_eq += oplen;
                    break;
                case BAM_CDIFF:
                    num_x += oplen;
                    break;
                case BAM_CINS:
                    num_i += oplen;
                    break;
                case BAM_CSOFT_CLIP:
                    num_s += oplen;
                    pos += (int32_t)oplen;
                    break;
                case BAM_CDEL:
                    num_d += oplen;
                    break;
                default:
                    break;
                }

                /* Q20 aligned bases */
                switch (op) {
                case BAM_CMATCH:
                case BAM_CEQUAL:
                case BAM_CDIFF:
                case BAM_CINS:
                    end_pos = pos + (int32_t)oplen;

                    /* for each aligned base check qual and increment pos */
                    while (pos < end_pos) {
                        if (qual[pos++] >= 20) {
                            ++am->b_q20;
                        }
                    }
                    break;
                default:
                    break;
                }
            }

            /*
             * Count mismatches from MD tag if present.
             * if num_mismatches was calculated from an MD tag, subtract dels
             * otherwise num_mismatches becomes num_x
             */
            /*
            if ((nm_tag = bam_aux_get(rec, "NM")) != NULL) {
                num_mismatches = bam_aux2i(nm_tag);
                num_mismatches -= num_d;
                num_matches = num_m - num_mismatches;
            } else
            */
            if ((md_tag = bam_aux_get(rec, "MD")) != NULL) {
                num_mismatches = 0;
                for (char *md_str = bam_aux2Z(md_tag); *md_str != '\0'; ++md_str) {
                    if (isalpha(*md_str)) {
                        ++num_mismatches;
                    }
                }
                num_mismatches -= num_d;
                num_matches = num_m - num_mismatches;
            } else {
                num_mismatches = num_x;
                num_matches = (num_m > num_eq) ? num_m - num_mismatches : num_eq;
            }

            /*
             * Set alignment metrics
             * b_matched + b_mismatched == num_m
             * num_m + b_inserted == b_aligned
             * b_aligned + b_soft_clipped == b_mapped == l_qseq
             */
            num_aligned = rec->core.l_qseq;
            if (num_s > 0) {
                ++am->r_soft_clipped;
                am->b_soft_clipped += num_s;
                num_aligned -= num_s;
            }
            am->b_aligned += num_aligned;
            am->b_matched += num_matches;
            am->b_mismatched += num_mismatches;
            am->b_inserted += num_i;
            am->b_deleted += num_d;

            /* Exact matches: no mismatches, insertions, or deletions */
            if (num_mismatches + num_d + num_i == 0) {
                ++am->r_exact_match;
                am->b_exact_match += num_aligned;
            }
        }
    }
}

/**
 * Add alignment metrics to report.
 */
void align_report(report_t *report, align_metrics_t *am, read_type_t rt)
{
    char *key_buffer = malloc(REPORT_BUFFER_SIZE * sizeof(char));
    die_on_alloc_fail(key_buffer);
    char *value_buffer = malloc(REPORT_BUFFER_SIZE * sizeof(char));
    die_on_alloc_fail(value_buffer);

    /*
     * Add key prefix for read1 and read2
     * All reads: prefix = ""
     */
    const char *prefix = (rt == RT_READ1) ? "R1_" : rt == RT_READ2 ? "R2_" : "";
    size_t prefix_len = strlen(prefix);
    size_t copy_size = REPORT_BUFFER_SIZE - prefix_len;
    char *key_start = key_buffer + prefix_len;

    copy_to_buffer(key_buffer, prefix, REPORT_BUFFER_SIZE);

    copy_to_buffer(key_start, "Yield_Reads", copy_size);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", am->r_total);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "Yield_Bases", copy_size);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", am->b_total);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "Unmapped_Reads", copy_size);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", am->r_unmapped);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "Unmapped_Reads_Pct", copy_size);
    print_pct(value_buffer, REPORT_BUFFER_SIZE, am->r_unmapped, am->r_total);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "Unmapped_Bases", copy_size);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", am->b_unmapped);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "Unmapped_Bases_Pct", copy_size);
    print_pct(value_buffer, REPORT_BUFFER_SIZE, am->b_unmapped, am->b_total);
    report_add_key_value(report, key_buffer, value_buffer);

    if (rt == RT_ALL) {
        copy_to_buffer(key_start, "Duplicate_Reads", copy_size);
        snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", am->r_dup);
        report_add_key_value(report, key_buffer, value_buffer);

        copy_to_buffer(key_start, "Duplicate_Reads_Pct", copy_size);
        print_pct(value_buffer, REPORT_BUFFER_SIZE, am->r_dup, am->r_aligned);
        report_add_key_value(report, key_buffer, value_buffer);

        copy_to_buffer(key_start, "Duplicate_Bases", copy_size);
        snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", am->b_dup);
        report_add_key_value(report, key_buffer, value_buffer);

        copy_to_buffer(key_start, "Duplicate_Bases_Pct", copy_size);
        print_pct(value_buffer, REPORT_BUFFER_SIZE, am->b_dup, am->b_aligned);
        report_add_key_value(report, key_buffer, value_buffer);
    }

    copy_to_buffer(key_start, "Mapped_Reads", copy_size);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", am->r_mapped);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "Mapped_Reads_Pct", copy_size);
    print_pct(value_buffer, REPORT_BUFFER_SIZE, am->r_mapped, am->r_total);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "Mapped_Bases", copy_size);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", am->b_mapped);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "Mapped_Bases_Pct", copy_size);
    print_pct(value_buffer, REPORT_BUFFER_SIZE, am->b_mapped, am->b_total);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "Aligned_Bases", copy_size);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", am->b_aligned);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "Aligned_Bases_Pct", copy_size);
    print_pct(value_buffer, REPORT_BUFFER_SIZE, am->b_aligned, am->b_total);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "Matched_Bases", copy_size);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", am->b_matched);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "Matched_Bases_Pct", copy_size);
    print_pct(value_buffer, REPORT_BUFFER_SIZE, am->b_matched, am->b_aligned);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "Mismatched_Bases", copy_size);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", am->b_mismatched);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "Mismatched_Bases_Pct", copy_size);
    print_pct(value_buffer, REPORT_BUFFER_SIZE, am->b_mismatched, am->b_aligned);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "Inserted_Bases", copy_size);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", am->b_inserted);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "Inserted_Bases_Pct", copy_size);
    print_pct(value_buffer, REPORT_BUFFER_SIZE, am->b_inserted, am->b_aligned);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "Deleted_Bases", copy_size);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", am->b_deleted);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "Deleted_Bases_Pct", copy_size);
    print_pct(value_buffer, REPORT_BUFFER_SIZE, am->b_deleted, am->b_aligned);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "SoftClipped_Reads", copy_size);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", am->r_soft_clipped);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "SoftClipped_Reads_Pct", copy_size);
    print_pct(value_buffer, REPORT_BUFFER_SIZE, am->r_soft_clipped, am->r_mapped);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "SoftClipped_Bases", copy_size);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", am->b_soft_clipped);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "SoftClipped_Bases_Pct", copy_size);
    print_pct(value_buffer, REPORT_BUFFER_SIZE, am->b_soft_clipped, am->b_mapped);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "Perfect_Reads", copy_size);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", am->r_exact_match);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "Perfect_Reads_Pct", copy_size);
    print_pct(value_buffer, REPORT_BUFFER_SIZE, am->r_exact_match, am->r_mapped);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "Perfect_Bases", copy_size);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", am->b_exact_match);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "Perfect_Bases_Pct", copy_size);
    print_pct(value_buffer, REPORT_BUFFER_SIZE, am->b_exact_match, am->b_mapped);
    report_add_key_value(report, key_buffer, value_buffer);
    
    copy_to_buffer(key_start, "MAPQ20_Reads", copy_size);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", am->r_mapq20);
    report_add_key_value(report, key_buffer, value_buffer);
    
    copy_to_buffer(key_start, "MAPQ20_Reads_Pct", copy_size);
    print_pct(value_buffer, REPORT_BUFFER_SIZE, am->r_mapq20, am->r_mapped);
    report_add_key_value(report, key_buffer, value_buffer);
    
    copy_to_buffer(key_start, "Q20_Bases", copy_size);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", am->b_q20);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "Q20_Bases_Pct", copy_size);
    print_pct(value_buffer, REPORT_BUFFER_SIZE, am->b_q20, am->b_aligned);
    report_add_key_value(report, key_buffer, value_buffer);

    free(key_buffer);
    free(value_buffer);
}
