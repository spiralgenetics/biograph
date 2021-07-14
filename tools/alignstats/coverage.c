#include "coverage.h"
#include "err.h"
#include "logging.h"
#include "print.h"
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

/**
 * Capture and whole genome statistics
 */

/* Target_Coverage helpers */
/*
char* base_cov_to_str(uint32_t *coverage) {
    //Currently an unused function.
    size_t length = sizeof(coverage)/sizeof(coverage[0]);
    size_t cur_int_buffer = 20;
    size_t cur_buffer_pos = 0;
    size_t buffer_length = length;
    char *ret = calloc(length, sizeof(char));

    for (size_t i = 0 ; i < length; i++)
    {
        if ((cur_buffer_pos + cur_int_buffer) > buffer_length) {
            char* new_str = (char*) malloc(strlen(ret) * sizeof(char) * 2);
            strcpy(new_str, ret);
            ret = new_str;
            buffer_length = sizeof(ret)/sizeof(ret[0]);
        }
        snprintf(ret + cur_buffer_pos, cur_int_buffer+1, "%d,", coverage[i]);
        cur_buffer_pos += cur_int_buffer + 1;
    }
    return ret;
    }*/

/**
 * Create and return new *target_coverage objects
 */
target_coverage_block_t *target_coverage_block_init(size_t target_count) {
    target_coverage_block_t *my_block = calloc(1, sizeof(target_coverage_block_t));
    die_on_alloc_fail(my_block);
    my_block->start_pos = NULL;
    my_block->start_pos = calloc(target_count, sizeof(uint32_t));
    my_block->end_pos = NULL;
    my_block->end_pos = calloc(target_count, sizeof(uint32_t));
    my_block->mean = NULL;
    my_block->mean = calloc(target_count, sizeof(float));
    my_block->min = NULL;
    my_block->min = calloc(target_count, sizeof(uint32_t));
    my_block->cov_lt5 = NULL;
    my_block->cov_lt5 = calloc(target_count, sizeof(uint32_t));
    my_block->cov_lt10 = NULL;
    my_block->cov_lt10 = calloc(target_count, sizeof(uint32_t));
    my_block->cov_lt20 = NULL;
    my_block->cov_lt20 = calloc(target_count, sizeof(uint32_t));
    my_block->base_coverage = NULL; 
    my_block->base_coverage = (uint32_t**) calloc(target_count, sizeof(uint32_t*));
    die_on_alloc_fail(my_block->base_coverage);

    return my_block;
}

target_coverage_t *target_coverage_init(uint32_t chrom_count) { 
    target_coverage_t *my_target = calloc(1, sizeof(target_coverage_t));
    die_on_alloc_fail(my_target);
    my_target->chroms = NULL;
    my_target->chroms = calloc(chrom_count, sizeof(target_coverage_block_t));
    die_on_alloc_fail(my_target->chroms);

    return my_target;
}
 

/* Coverage info structure */

/**
 * Create and return new *coverage_info_t.
 */
coverage_info_t *coverage_info_init()
{
    coverage_info_t *ci = calloc(1, sizeof(coverage_info_t));
    die_on_alloc_fail(ci);

    ci->cov_histo_len = 0x40000; /* == 2^18 == 262144 */
    ci->cov_histo = calloc(ci->cov_histo_len, sizeof(uint64_t));
    die_on_alloc_fail(ci->cov_histo);

    return ci;
}

/**
 * Free *coverage_info_t from memory.
 */
void coverage_info_destroy(coverage_info_t *ci)
{
    if (ci != NULL) {
        free(ci->cov_histo);
        free(ci);
    }
}

/* Capture metrics structure */

/**
 * Create and return new *capture_metrics_t.
 */
capture_metrics_t *capture_metrics_init(bed_t *target_design)
{
    capture_metrics_t *cm = calloc(1, sizeof(capture_metrics_t));
    cm->t_total = target_design->num_targets;
    cm->t_target_cov = target_coverage_init(target_design->num_chroms);
    cm->t_target_cov->chrom_names = target_design->chrom_names;
    
    // For every chrom in target_design and every target per chrom, initialize the object
    for (int cur_chrom = 0; cur_chrom < target_design->num_chroms; cur_chrom++) {
        bed_chrom_t *cur_chrom_bed = target_design->chroms[cur_chrom];
        cm->t_target_cov->chroms[cur_chrom] = target_coverage_block_init(cur_chrom_bed->num_targets);
        for (int chrom_target = 0; chrom_target < cur_chrom_bed->num_targets; chrom_target++) {
            //off by one?
            cm->t_target_cov->chroms[cur_chrom]->start_pos[chrom_target] = cur_chrom_bed->start_pos[chrom_target];
            cm->t_target_cov->chroms[cur_chrom]->end_pos[chrom_target] = cur_chrom_bed->end_pos[chrom_target];
            size_t span = cur_chrom_bed->end_pos[chrom_target] - cur_chrom_bed->start_pos[chrom_target];
            cm->t_target_cov->chroms[cur_chrom]->base_coverage[chrom_target] = calloc(span, sizeof(uint32_t));
        }
    }
    
    die_on_alloc_fail(cm);

    return cm;
}

/**
 * Finalize capture metrics once all records are processed.
 * Calculates median coverage.
 */
void capture_metrics_finalize(capture_metrics_t *cm, coverage_info_t *ci, bed_t *ti)
{
    uint64_t sum = 0, mid = (ti == NULL ? cm->b_total : cm->b_targeted) / 2;

    /* Set median coverage */
    for (size_t i = 0; i < ci->cov_histo_len; ++i) {
        if ((sum += ci->cov_histo[i]) >= mid) {
            cm->c_median = i;
            break;
        }
    }

    /* Masked bases */
    cm->b_total -= cm->b_masked;
}

/**
 * Free *capture_metrics_t from memory.
 */
void capture_metrics_destroy(capture_metrics_t *cm)
{
    free(cm);
}

/* Capture target and coverage statistics calculation */

/**
 * Increment ci->cov_histo[cov], realloc if needed
 */
void incr_cov_histo(coverage_info_t *ci, uint32_t cov)
{
    uint64_t *tmp_cov_histo;

    if (cov + 1 > ci->cov_histo_len) {
        /* Buffer a small amount (256) to reduce # of reallocs for slowly increasing coverage */
        ci->cov_histo_len = (size_t)(cov + 1 + 256);
        tmp_cov_histo = realloc(ci->cov_histo, ci->cov_histo_len * sizeof(uint64_t));

        if (tmp_cov_histo != NULL) {
            ci->cov_histo = tmp_cov_histo;
        } else {
            free(ci->cov_histo);
            die_on_alloc_fail(tmp_cov_histo);
        }
    }

    ++ci->cov_histo[cov];
}

/**
 * Record whole genome coverage metrics for chromosome cname.
 */
void handle_wgs_coverage(const uint32_t *coverage, capture_metrics_t *cm,
                         coverage_info_t *ci, int32_t chrom_len)
{
    uint32_t cov;

    /* for each base position in chromosome */
    for (int32_t i = 0; i < chrom_len; ++i) {
        cov = coverage[i];

        /* Bases with coverage of at least 1, 10, 20, etc. */
        if (cov >= 1) {
            ++cm->b_1_plus_hits;
            if (cov >= 10) {
                ++cm->b_10_plus_hits;
                if (cov >= 20) {
                    ++cm->b_20_plus_hits;
                    if (cov >= 30) {
                        ++cm->b_30_plus_hits;
                        if (cov >= 40) {
                            ++cm->b_40_plus_hits;
                            if (cov >= 50) {
                                ++cm->b_50_plus_hits;
                                if (cov >= 100) {
                                    ++cm->b_100_plus_hits;
                                    if (cov >= 500) {
                                        ++cm->b_500_plus_hits;
                                        if (cov >= 1000) {
                                            ++cm->b_1000_plus_hits;
        }   }   }   }   }   }   }   }   }

        incr_cov_histo(ci, cov);
        cm->c_total += cov;
    }
}

/**
 * (getTargetsAndWriteCoverage)
 * Record capture coverage metrics for chromosome cname.
 */
void handle_target_coverage(const uint32_t *coverage, capture_metrics_t *cm,
                            coverage_info_t *ci, bed_t *ti, int32_t chrom_idx,
                            const char *chrom, int32_t chrom_len)
{
    /* Target coverage statistics */
    bool target_hit, buffer_hit; /*, space_fasta = false;*/
    uint32_t cov;
    int32_t start, end, j, buffer_end;
    bed_chrom_t *tic = ti->chroms[chrom_idx];
    
    target_coverage_block_t *cur_chrom_cov = cm->t_target_cov->chroms[chrom_idx];
    /* for each target */
    for (size_t tgt_idx = 0; tgt_idx < tic->num_targets; ++tgt_idx) {
        start = tic->start_pos[tgt_idx];
        end = tic->end_pos[tgt_idx];
        target_hit = false;
        uint32_t *cur_target_cov = cur_chrom_cov->base_coverage[tgt_idx];
        float sum_cov = 0;
        uint32_t min_cov = UINT_MAX;
        /* for each base position */
        for (int32_t j = start; j <= end; ++j) {
            cov = coverage[j];
            /* Copy coverage over */
            cur_target_cov[j - start] = cov;
            sum_cov += cov;
            if (cov < min_cov) {
                min_cov = cov;
            }
            
            if (cov < 20) {
                ++cur_chrom_cov->cov_lt20[tgt_idx];
                if (cov < 10) {
                    ++cur_chrom_cov->cov_lt10[tgt_idx];
                    if (cov < 5) {
                        ++cur_chrom_cov->cov_lt5[tgt_idx];
                    
            } } } 
            
            /* Bases with coverage of at least 1, 10, 20, etc. */
            if (cov >= 1) {
                ++cm->b_1_plus_hits;
                target_hit = true;
                if (cov >= 10) {
                    ++cm->b_10_plus_hits;
                    if (cov >= 20) {
                        ++cm->b_20_plus_hits;
                        if (cov >= 30) {
                            ++cm->b_30_plus_hits;
                            if (cov >= 40) {
                                ++cm->b_40_plus_hits;
                                if (cov >= 50) {
                                    ++cm->b_50_plus_hits;
                                    if (cov >= 100) {
                                        ++cm->b_100_plus_hits;
                                        if (cov >= 500) {
                                            ++cm->b_500_plus_hits;
                                            if (cov >= 1000) {
                                                ++cm->b_1000_plus_hits;
            }   }   }   }   }   }   }   }   }

            incr_cov_histo(ci, cov);
            cm->c_total += cov;
        }
        //And record average
        cur_chrom_cov->mean[tgt_idx] = sum_cov / (end - start);
        cur_chrom_cov->min[tgt_idx] = min_cov;
        
        if (target_hit) {
            ++cm->t_hit;
        } else {
            /* Check if buffers were hit */
            buffer_hit = false;

            /* Left buffer */
            j = (start > BUFFER) ? start - BUFFER : 0;
            buffer_end = (start < chrom_len) ? start : chrom_len - 1;

            while (j < buffer_end) {
                if (coverage[j++] > 0) {
                    buffer_hit = true;
                    break;
                }
            }

            if (buffer_hit) {
                ++cm->t_buffers_hit;
            } else {
                /* Right buffer */
                ++end;
                j = (end > 0) ? end : 0;
                buffer_end =
                    (end + BUFFER < chrom_len) ? end + BUFFER : chrom_len - 1;

                while (j < buffer_end) {
                    if (coverage[j++] > 0) {
                        ++cm->t_buffers_hit;
                        break;
                    }
                }
            }
        }
    }
}

/**
 * Set values from start to end in coverage to 0.
 */
void clear_coverage(uint32_t *coverage, int32_t start, int32_t end, int32_t chrom_len)
{
    if (start < 0) {
        start = 0;
    }
    if (end >= chrom_len) {
        end = chrom_len - 1;
    }

    /* Clear target and buffer bases */
    memset(coverage + start, 0, (end - start + 1) * sizeof(uint32_t));
}

/**
 * (findWhereReadsHit)
 * Record contiguous regions of >= 20X coverage outside of target and buffer
 * regions. This function is destructive to coverage, only perform once all
 * records for a chromosome have been processed.
 */
void handle_miss_reads(uint32_t *coverage, capture_metrics_t *cm, bed_t *ti,
                       int32_t chrom_idx, int32_t chrom_len)
{
    if (ti->num_targets > 0) {
        bed_chrom_t *tic = ti->chroms[chrom_idx];

        /* for each target set coverage[start-buffer:end+buffer] values to 0 */
        for (size_t i = 0; i < tic->num_targets; ++i) {
            clear_coverage(coverage, tic->start_pos[i] - MISS_BUFFER,
                           tic->end_pos[i] + MISS_BUFFER, chrom_len);
        }

        /* for each base in coverage with targets cleared */
        for (int32_t j = 0; j < chrom_len; ++j) {
            if (coverage[j] >= 20) {
                ++cm->t_non_target_good_hits;

                /* Scoot j past this contiguous region of coverage */
                while (j < chrom_len && coverage[j] > 0) {
                    ++j;
                }
            }
        }
    }
}

/**
 * Set coverage values for regions of known N bases in reference to 0.
 * Regions are defined as targets in cov_mask_ti.
 */
void handle_coverage_mask(uint32_t *coverage, bed_t *cov_mask_ti,
                          int32_t chrom_idx, int32_t chrom_len)
{
    if (cov_mask_ti->num_targets > 0) {
        bed_chrom_t *tic = cov_mask_ti->chroms[chrom_idx];

        /* for each target set coverage[start:end] values to 0 */
        for (size_t i = 0; i < tic->num_targets; ++i) {
            clear_coverage(coverage, tic->start_pos[i], tic->end_pos[i], chrom_len);
        }
    }
}

/**
 * Erase target regions overlapping masked regions
 * Note - Cap_Q20_Bases cannot be masked, so if there's overlap between the mask
 * and target, Cap_Q20_Bases will be inflated
 */
void handle_coverage_mask_target(uint8_t *target_cov, capture_metrics_t *cm,
                                 bed_t *cov_mask_ti, int32_t chrom_idx,
                                 int32_t chrom_len)
{
    uint8_t *curr_pos, *end_pos;
    int32_t start, end;
    bed_chrom_t *tic;
    
    if (cov_mask_ti->num_targets > 0) {
        tic = cov_mask_ti->chroms[chrom_idx];

        /* for each target */
        for (size_t i = 0; i < tic->num_targets; ++i) {
            start = tic->start_pos[i];
            end = tic->end_pos[i];

            for (curr_pos = target_cov + start, end_pos = target_cov + end;
                 curr_pos <= end_pos;
                 ++curr_pos)
            {
                switch (*curr_pos) {
                case TARGET_IN:
                    --cm->b_targeted;
                    break;
                case TARGET_BUFFER:
                    --cm->b_buffer;
                    break;
                default:
                    break;
                }
            }

            /* Set masked region to TARGET_OUT */
            memset(target_cov + start, TARGET_OUT, (end - start + 1) * sizeof(uint8_t));
        }
    }
}

/**
 * Set target positions in target_cov.
 */
void set_target_cov(uint8_t *target_cov, capture_metrics_t *cm, bed_t *ti,
                    int32_t chrom_idx, int32_t chrom_len)
{
    uint8_t *curr_pos, *end_pos;
    int32_t start, end, start_, end_, target_start;
    bed_chrom_t *tic;

    if (ti->num_targets > 0) {
        tic = ti->chroms[chrom_idx];

        /* Clear target_cov */
        memset(target_cov, TARGET_OUT, chrom_len * sizeof(uint8_t));

        /* for each target */
        for (size_t j = 0; j < tic->num_targets; ++j) {
            start = tic->start_pos[j];
            end = tic->end_pos[j];

            /* Left buffer */
            start_ = (start > BUFFER) ? start - BUFFER : 0;
            end_ = (start < chrom_len) ? start : chrom_len - 1;

            for (curr_pos = target_cov + start_, end_pos = target_cov + end_;
                 curr_pos < end_pos;
                 ++curr_pos)
            {
                if (*curr_pos == TARGET_OUT) {
                    *curr_pos = TARGET_BUFFER;
                }
            }

            /* Target */
            target_start = (start > 0) ? start : 0;
            memset(target_cov + target_start, TARGET_IN,
                   ((end < chrom_len) ? end : chrom_len - 1) - target_start + 1);

            /* Right buffer */
            ++end;
            start_ = (end > 0) ? end : 0;
            end_ = (end + BUFFER < chrom_len) ? end + BUFFER : chrom_len - 1;

            for (curr_pos = target_cov + start_, end_pos = target_cov + end_;
                 curr_pos < end_pos;
                 ++curr_pos)
            {
                if (*curr_pos == TARGET_OUT) {
                    *curr_pos = TARGET_BUFFER;
                }
            }
        }

        /* Record number of buffer and targeted bases in target_cov */
        for (curr_pos = target_cov, end_pos = target_cov + chrom_len;
             curr_pos < end_pos;
             ++curr_pos)
        {
            switch (*curr_pos) {
            case TARGET_IN:
                ++cm->b_targeted;
                break;
            case TARGET_BUFFER:
                ++cm->b_buffer;
                break;
            default:
                break;
            }
        }
    }
}

void _capture_process_record1(bam1_t *rec, capture_metrics_t *cm)
{
    ++cm->r_aligned;

    /* if paired read */
    if (rec->core.flag & BAM_FPAIRED) {
        ++cm->r_paired;

        /* if mate is mapped */
        if (!(rec->core.flag & BAM_FMUNMAP)) {
            ++cm->r_paired_w_mate;
        }
    }

    /* if duplicate read */
    if (rec->core.flag & BAM_FDUP) {
        ++cm->r_dup;
    }
}

/**
 * Process coverage metrics specific to whole genome or capture.
 */
void _capture_process_record2(bam1_t *rec, capture_metrics_t *cm,
                              target_state_t target_status)
{
    cm->b_aligned += rec->core.l_qseq;

    /* Record read as out of target, in target buffer, or in target */
    switch (target_status) {
    case TARGET_IN:
        ++cm->r_in_target;
        cm->b_on_target += rec->core.l_qseq;
        if (rec->core.qual >= 20) {
            ++cm->r_in_target_mapq20;
            cm->b_in_target_mapq20 += rec->core.l_qseq;
        }
        break;
    case TARGET_BUFFER:
        ++cm->r_in_buffer;
        break;
    default:
        ++cm->r_out_target;
        break;
    }
}

/**
 * Process record for target and coverage info.
 * cm_wgs: Capture metrics for whole genome statistics
 * cm_cap: Capture metrics for capture statistics
 */
void capture_process_record(bam1_t *rec, uint32_t *coverage,
                            const uint8_t *target_cov,
                            capture_metrics_t *cm_wgs,
                            capture_metrics_t *cm_cap, int32_t chrom_len,
                            bool remove_dups)
{
    bool in_target, in_buffer;
    uint8_t *qual, *bqual;
    int32_t pos, start, start_pos, end_pos, ref_pos;
    uint32_t oplen, *cigar;
    const uint16_t FILTER = BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL;
    const uint32_t COV_MAX = UINT32_MAX - 1;


    if (cm_wgs != NULL) {
        ++cm_wgs->r_total;
    }
    if (cm_cap != NULL) {
        ++cm_cap->r_total;
    }

    /* Filter out reads with filtered flags */
    if (rec->core.flag & FILTER) {
        return;
    }

    if (cm_wgs != NULL) {
        _capture_process_record1(rec, cm_wgs);
    }
    if (cm_cap != NULL) {
        _capture_process_record1(rec, cm_cap);
    }

    /* remove duplicate reads */
    if (remove_dups && (rec->core.flag & BAM_FDUP)) {
        return;
    }

    /*
     * Record coverage and whether aligned bases hit target or buffer
     * Confirmed this now matches up with the correct ref pos and # CIGAR Ms
     */
    in_target = false;
    in_buffer = false;
    ref_pos = 0;
    start = rec->core.pos;
    cigar = bam_get_cigar(rec);
    bqual = bam_get_qual(rec);

    /* for each CIGAR op */
    for (uint32_t i = 0; i < rec->core.n_cigar; ++i) {
        oplen = bam_cigar_oplen(cigar[i]);

        switch (bam_cigar_op(cigar[i])) {
        case BAM_CHARD_CLIP:
            break;
        /*
         * M, =, X: record coverage for CIGAR M bases (matches + mismatches)
         * or = (matches) and X (mismatches)
         */
        case BAM_CMATCH:
        case BAM_CEQUAL:
        case BAM_CDIFF:
            start_pos = pos = start + ref_pos;
            end_pos = start_pos + oplen - 1;

            if (start_pos < 0) {
                start_pos = 0;
                pos = 0;
            }
            if (end_pos >= chrom_len) {
                end_pos = chrom_len - 1;
            }

            /* Record coverage */
            qual = bqual + ref_pos;
            while (pos <= end_pos) {
                if (coverage[pos] < COV_MAX) {
                    ++coverage[pos];
                } else {
                    /* Y'know just in case */
                    log_warning("Coverage of greater than %u detected. "
                                "Coverage statistics may not be accurate.", COV_MAX);
                }
                ++pos;
                ++qual;
            }
            /* pos == end_pos */

            if (target_cov != NULL) {
                while (--pos >= start_pos) {
                    if (target_cov[pos] == TARGET_BUFFER) {
                        in_buffer = true;
                    } else if (target_cov[pos] == TARGET_IN) {
                        in_target = true;
                        break;
                    }
                }
            }
            ref_pos += oplen;
            break;
        /* D: advance ref_pos past deletion */
        case BAM_CDEL:
            ref_pos += oplen;
            break;
        default:
            break;
        }
    }

    if (cm_wgs != NULL) {
        _capture_process_record2(rec, cm_wgs, TARGET_OUT);
    }
    if (cm_cap != NULL) {
        _capture_process_record2(rec, cm_cap,
            in_target ? TARGET_IN : in_buffer ? TARGET_BUFFER : TARGET_OUT);
    }
}

/**
 * Write capture metrics to report.
 */
void capture_report(report_t *report, capture_metrics_t *cm, bed_t *ti)
{
    char *key_buffer = malloc(REPORT_BUFFER_SIZE * sizeof(char));
    die_on_alloc_fail(key_buffer);
    char *value_buffer = malloc(REPORT_BUFFER_SIZE * sizeof(char));
    die_on_alloc_fail(value_buffer);

    const char *prefix = (ti == NULL) ? "Wgs_" : "Cap_";
    size_t prefix_len = strlen(prefix);
    size_t copy_size = REPORT_BUFFER_SIZE - prefix_len;
    char *key_start = key_buffer + prefix_len;

    /* Capture (targets) or whole genome (no targets) percentage denominator */
    uint64_t denominator = (ti == NULL) ? cm->b_total : cm->b_targeted;

    copy_to_buffer(key_buffer, prefix, REPORT_BUFFER_SIZE);

    copy_to_buffer(key_start, "Total_Reads", copy_size);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", cm->r_total);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "Cov_Duplicate_Reads", copy_size);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", cm->r_dup);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "Cov_Duplicate_Reads_Pct", copy_size);
    print_pct(value_buffer, REPORT_BUFFER_SIZE, cm->r_dup, cm->r_total);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "Aligned_Reads", copy_size);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", cm->r_aligned);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "Aligned_Reads_Pct", copy_size);
    print_pct(value_buffer, REPORT_BUFFER_SIZE, cm->r_aligned, cm->r_total);
    report_add_key_value(report, key_buffer, value_buffer);
    

    
    copy_to_buffer(key_start, "Reads_Paired", copy_size);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", cm->r_paired);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "Reads_Paired_With_Mates", copy_size);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", cm->r_paired_w_mate);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "Average_Coverage", copy_size);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%.2f",
             (denominator != 0) ? (double)cm->c_total / (double)denominator
                                : 0.0);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "Median_Coverage", copy_size);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", cm->c_median);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "Expected_Aligned_Reads", copy_size);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", cm->r_aligned);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "Calculated_Aligned_Reads", copy_size);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", cm->r_in_target + cm->r_in_buffer + cm->r_out_target);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "Coverage_Bases_1", copy_size);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", cm->b_1_plus_hits);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "Coverage_Bases_1_Pct", copy_size);
    print_pct(value_buffer, REPORT_BUFFER_SIZE, cm->b_1_plus_hits, denominator);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "Coverage_Bases_10", copy_size);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", cm->b_10_plus_hits);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "Coverage_Bases_10_Pct", copy_size);
    print_pct(value_buffer, REPORT_BUFFER_SIZE, cm->b_10_plus_hits, denominator);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "Coverage_Bases_20", copy_size);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", cm->b_20_plus_hits);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "Coverage_Bases_20_Pct", copy_size);
    print_pct(value_buffer, REPORT_BUFFER_SIZE, cm->b_20_plus_hits, denominator);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "Coverage_Bases_30", copy_size);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", cm->b_30_plus_hits);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "Coverage_Bases_30_Pct", copy_size);
    print_pct(value_buffer, REPORT_BUFFER_SIZE, cm->b_30_plus_hits, denominator);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "Coverage_Bases_40", copy_size);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", cm->b_40_plus_hits);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "Coverage_Bases_40_Pct", copy_size);
    print_pct(value_buffer, REPORT_BUFFER_SIZE, cm->b_40_plus_hits, denominator);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "Coverage_Bases_50", copy_size);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", cm->b_50_plus_hits);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "Coverage_Bases_50_Pct", copy_size);
    print_pct(value_buffer, REPORT_BUFFER_SIZE, cm->b_50_plus_hits, denominator);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "Coverage_Bases_100", copy_size);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", cm->b_100_plus_hits);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "Coverage_Bases_100_Pct", copy_size);
    print_pct(value_buffer, REPORT_BUFFER_SIZE, cm->b_100_plus_hits, denominator);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "Coverage_Bases_1000", copy_size);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", cm->b_1000_plus_hits);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "Coverage_Bases_1000_Pct", copy_size);
    print_pct(value_buffer, REPORT_BUFFER_SIZE, cm->b_1000_plus_hits, denominator);
    report_add_key_value(report, key_buffer, value_buffer);

    /* ti != NULL? capture stats: wgs */
    if (ti != NULL) {
        copy_to_buffer(key_start, "Buffer_Aligned_Reads", copy_size);
        snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", cm->r_in_buffer);
        report_add_key_value(report, key_buffer, value_buffer);

        copy_to_buffer(key_start, "Buffer_Aligned_Reads_Pct", copy_size);
        print_pct(value_buffer, REPORT_BUFFER_SIZE, cm->r_in_buffer, cm->r_aligned);
        report_add_key_value(report, key_buffer, value_buffer);

        copy_to_buffer(key_start, "Target_Aligned_Reads", copy_size);
        snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", cm->r_in_target);
        report_add_key_value(report, key_buffer, value_buffer);

        copy_to_buffer(key_start, "Target_Aligned_Reads_Pct", copy_size);
        print_pct(value_buffer, REPORT_BUFFER_SIZE, cm->r_in_target, cm->r_aligned);
        report_add_key_value(report, key_buffer, value_buffer);
        
        copy_to_buffer(key_start, "Target_Aligned_Bases", copy_size);
        snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", cm->b_on_target);
        report_add_key_value(report, key_buffer, value_buffer);

        copy_to_buffer(key_start, "Target_Aligned_Bases_Pct", copy_size);
        print_pct(value_buffer, REPORT_BUFFER_SIZE, cm->b_on_target, cm->b_aligned);
        report_add_key_value(report, key_buffer, value_buffer);

        copy_to_buffer(key_start, "Target_MAPQ20_Reads", copy_size);
        snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", cm->r_in_target_mapq20);
        report_add_key_value(report, key_buffer, value_buffer);

        copy_to_buffer(key_start, "Target_MAPQ20_Reads_Pct", copy_size);
        print_pct(value_buffer, REPORT_BUFFER_SIZE, cm->r_in_target_mapq20, cm->r_aligned);
        report_add_key_value(report, key_buffer, value_buffer);
        
        copy_to_buffer(key_start, "Target_MAPQ20_Bases", copy_size);
        snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", cm->b_in_target_mapq20);
        report_add_key_value(report, key_buffer, value_buffer);
        
        copy_to_buffer(key_start, "Target_MAPQ20_Bases_Pct", copy_size);
        print_pct(value_buffer, REPORT_BUFFER_SIZE, cm->b_in_target_mapq20, cm->b_on_target);
        report_add_key_value(report, key_buffer, value_buffer);       
        
        copy_to_buffer(key_start, "Targets_Hit", copy_size);
        snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", cm->t_hit);
        report_add_key_value(report, key_buffer, value_buffer);

        copy_to_buffer(key_start, "Targets_Hit_Pct", copy_size);
        print_pct(value_buffer, REPORT_BUFFER_SIZE, cm->t_hit, cm->t_total);
        report_add_key_value(report, key_buffer, value_buffer);

        copy_to_buffer(key_start, "Target_Buffers_Hit", copy_size);
        snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", cm->t_buffers_hit);
        report_add_key_value(report, key_buffer, value_buffer);

        copy_to_buffer(key_start, "Target_Buffers_Hit_Pct", copy_size);
        print_pct(value_buffer, REPORT_BUFFER_SIZE, cm->t_buffers_hit, cm->t_total);
        report_add_key_value(report, key_buffer, value_buffer);

        copy_to_buffer(key_start, "Total_Targets", copy_size);
        snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", cm->t_total);
        report_add_key_value(report, key_buffer, value_buffer);

        copy_to_buffer(key_start, "High_Coverage_Non_Target_Hits", copy_size);
        snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", cm->t_non_target_good_hits);
        report_add_key_value(report, key_buffer, value_buffer);

        copy_to_buffer(key_start, "Bases_On_Target", copy_size);
        snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", cm->b_targeted);
        report_add_key_value(report, key_buffer, value_buffer);
        

        
        copy_to_buffer(key_start, "Bases_On_Buffer", copy_size);
        snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", cm->b_buffer);
        report_add_key_value(report, key_buffer, value_buffer);

        copy_to_buffer(key_start, "Reads_On_Target_Or_Buffer", copy_size);
        snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", cm->r_in_target + cm->r_in_buffer);
        report_add_key_value(report, key_buffer, value_buffer);

        copy_to_buffer(key_start, "Reads_On_Target_Or_Buffer_Pct", copy_size);
        print_pct(value_buffer, REPORT_BUFFER_SIZE, cm->r_in_target + cm->r_in_buffer, cm->r_aligned);
        report_add_key_value(report, key_buffer, value_buffer);
        
        //For each target - Report coverage
        copy_to_buffer(key_start, "Target_Coverage_Header", copy_size);
        snprintf(value_buffer, REPORT_BUFFER_SIZE, "chrom start end mean_cov min_cov cov_lt5 cov_lt10 cov_lt20");
        report_add_key_value(report, key_buffer, value_buffer);
        for(int i = 0; i < ti->num_chroms; i++) {
            for (int j = 0; j < ti->chroms[i]->num_targets; j++) {
                copy_to_buffer(key_start, "Target_Coverage", copy_size);
                snprintf(value_buffer, REPORT_BUFFER_SIZE, "%s %u %u %.2f %u %u %u %u", 
                         cm->t_target_cov->chrom_names[i],
                         cm->t_target_cov->chroms[i]->start_pos[j],
                         cm->t_target_cov->chroms[i]->end_pos[j],
                         cm->t_target_cov->chroms[i]->mean[j],
                         cm->t_target_cov->chroms[i]->min[j], 
                         cm->t_target_cov->chroms[i]->cov_lt5[j], 
                         cm->t_target_cov->chroms[i]->cov_lt10[j], 
                         cm->t_target_cov->chroms[i]->cov_lt20[j]);
                report_add_key_value(report, key_buffer, value_buffer);
                
                //need a way to output the base_coverage metrics...
                // base_cov_to_str(cm->t_target_cov->chroms[i]->base_coverage[j]));
            }
        }
    }

    free(key_buffer);
    free(value_buffer);
}
