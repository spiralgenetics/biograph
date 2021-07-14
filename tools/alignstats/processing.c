#include "processing.h"
#include "align.h"
#include "alignlen.h"
#include "alignstats.h"
#include "coverage.h"
#include "filter.h"
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "insertsize.h"
#include "logging.h"
#include "pairstats.h"
#include <stdbool.h>
#include <stdint.h>

#ifdef USE_PTHREAD
#include <pthread.h>
#endif

/**
 * Functions for reading input alignment files and processing alignment
 * records.
 */

/**
 * Set args->iter to new region after calling either move_to_first_region()
 * or move_to_first_region().
 */
void set_iter(args_t *args)
{
    bed_chrom_t *chrom = args->regions->chroms[args->regions_curr_chrom_idx];

    /* destroy previous iter if set */
    if (args->iter != NULL) {
        sam_itr_destroy(args->iter);
    }

    /* set new iter */
    args->iter = sam_itr_queryi(args->index, args->regions_curr_chrom_idx,
                                chrom->start_pos[args->regions_curr_target_idx],
                                chrom->end_pos[args->regions_curr_target_idx]);

    if (args->iter == NULL) {
        log_warning("Unable to set iterator in region %zu", args->regions_curr_chrom_idx);
    }
}

/**
 * Move to first region in args->region->chroms.
 * Returns whether moved to next region.
 */
bool move_to_first_region(args_t *args)
{
    bool rv = true;
    bed_chrom_t **chroms = args->regions->chroms;

    if (args->regions_curr_chrom_idx >= args->regions->num_chroms) {
        rv = false;
    } else {
        while (rv && chroms[args->regions_curr_chrom_idx]->num_targets == 0) {
            if (args->regions_curr_chrom_idx < args->regions->num_chroms) {
                ++args->regions_curr_chrom_idx;
                args->regions_curr_target_idx = 0;
            } else {
                rv = false;
            }
        }
        set_iter(args);
    }

    return rv;
}

/**
 * Move to first region in args->region->chroms. Only use this after moving
 * to first regions using move_to_first_region().
 * Returns whether moved to next region.
 */
bool move_to_next_region(args_t *args)
{
    bool rv = true;
    bed_chrom_t **chroms = args->regions->chroms;

    /* no more regions */
    if (args->regions_curr_chrom_idx >= args->regions->num_chroms) {
        rv = false;
    } else {
        ++args->regions_curr_target_idx;

        /* find next chrom */
        if (args->regions_curr_target_idx >= chroms[args->regions_curr_chrom_idx]->num_targets) {
            do {
                ++args->regions_curr_chrom_idx;
                if (args->regions_curr_chrom_idx < args->regions->num_chroms) {
                    args->regions_curr_target_idx = 0;
                } else {
                    rv = false;
                }
            } while (rv && chroms[args->regions_curr_chrom_idx]->num_targets == 0);
        }

        if (rv) {
            set_iter(args);
        }
    }

    return rv;
}

/**
 * Base function for reading bam1_t records into args->read_buff
 * Reads records using sam_itr_next
 * Returns number of records read.
 */
uint32_t read_bam1(args_t *args)
{
    /* for each of up to the next RECORD_BUFFER_SIZE reads */
    uint32_t record_idx = 0;

    while (record_idx < args->reads_per_buffer &&
           sam_read1(args->input_sf, args->hdr, args->read_buff[record_idx]) >= 0)
    {
        ++record_idx;
    }

    return record_idx;
}

/**
 * Base function for reading bam1_t records into args->read_buff.
 * Reads records using sam_itr_next
 * Returns number of records read.
 */
uint32_t read_bam_itr(args_t *args)
{
    uint32_t record_idx = 0;

    while (record_idx < args->reads_per_buffer) {
        if (sam_itr_next(args->input_sf, args->iter, args->read_buff[record_idx]) < 0) {
            if (!move_to_next_region(args)) {
                /* process '*' chrom reads */
                if (args->process_unmapped && !args->process_unmapped_done) {
                    if (args->iter != NULL) {
                        sam_itr_destroy(args->iter);
                    }
                    args->iter = sam_itr_querys(args->index, args->hdr, "*");
                    args->process_unmapped_done = true;
                } else {
                    break;
                }
            }
        } else {
            ++record_idx;
        }
    }

    return record_idx;
}

/**
 * Evalute whether current record tid is in coordinate-sorted order with
 * respect to tid for previous mapped chromosome (i.e. not "*").
 */
#define check_order_tid(prev_mapped_tid, curr_tid) \
    (((curr_tid) == -1) || (prev_mapped_tid) <= (curr_tid))

/**
 * Evalute whether current record pos is in coordinate-sorted order with
 * respect to previous record pos.
 */
#define check_order_pos(prev_pos, curr_pos) \
    ((prev_pos) <= (curr_pos))

/**
 * Base function for processing bam1_t records in args->curr_buff.
 * Returs number of records processed
 */
uint32_t process_records(args_t *args)
{
    bool is_read_filtered = false;
    int32_t prev_rec_pos = -2;
    uint32_t record_idx = 0;
    bam1_t *rec;

    /* for each read in read buffer reads */
    while (record_idx < args->read_buff_size) {
        rec = args->curr_buff[record_idx++];

        /* if new chromosome */
        if (rec->core.tid != args->curr_chrom_idx) {
            args->new_chrom = true;
            args->prev_chrom_idx = args->curr_chrom_idx;
            if (args->curr_chrom_idx >= 0) {
                args->prev_mapped_chrom_idx = args->curr_chrom_idx;
            }
            args->prev_chrom_name = args->curr_chrom_name;
            args->prev_chrom_len = args->curr_chrom_len;

            if ((args->curr_chrom_idx = rec->core.tid) >= 0) {
                args->curr_chrom_name = args->hdr->target_name[args->curr_chrom_idx];
                args->curr_chrom_len = args->hdr->target_len[args->curr_chrom_idx];
            } else {
                args->curr_chrom_name = NULL;
                args->curr_chrom_len = 0;
            }

            if (args->order_warn && !check_order_tid(args->prev_mapped_chrom_idx, args->curr_chrom_idx)) {
                log_warning("Record not in sorted order. Results may not be accurate: "
                            "tid %d > %d", args->prev_mapped_chrom_idx,
                                           args->curr_chrom_idx);
                args->order_warn = false;
            }
            prev_rec_pos = -2;

            /* Don't process CIGAR if RNAME is "*" */
            args->process_cigar = (args->curr_chrom_name == NULL ||
                                   *args->curr_chrom_name != '*');
        }

        if (args->order_warn && !check_order_pos(prev_rec_pos, rec->core.pos)) {
            log_warning("Record not in sorted order. Results may not be accurate: "
                        "tid %d, pos %d > %d", args->curr_chrom_idx,
                                               prev_rec_pos + 1,
                                               rec->core.pos + 1);
            args->order_warn = false;
        }

        /* Preliminary filtering */
        filter_counter_process_record(rec, args->fc);
        is_read_filtered = (filter_test_qual(rec->core.qual, args->fc->min_qual) ||
                            filter_test_flag(rec->core.flag, args->fc->filter_incl, args->fc->filter_excl));

        if (!is_read_filtered)
        {
            prev_rec_pos = rec->core.pos;

            /* Alignment stats */
            if (args->do_alignment) {
                align_process_record(rec, args->am_all, args->process_cigar);
                align_len_process_record(rec, args->alm_all);

                if (rec->core.flag & BAM_FREAD1) {
                    align_process_record(rec, args->am_read1, args->process_cigar);
                    align_len_process_record(rec, args->alm_read1);
                } else if (rec->core.flag & BAM_FREAD2) {
                    align_process_record(rec, args->am_read2, args->process_cigar);
                    align_len_process_record(rec, args->alm_read2);
                }
                /*
                align_process_record(rec, args->am_fragment, args->process_cigar);
                align_len_process_record(rec, args->alm_fragment);
                */

                pair_stats_process_record(rec, args->psm);
                insert_size_process_record(rec, args->ism);
            }
        }

        /* Whole genome or capture stats */
        if (args->do_wgs || args->do_capture) {
            /* if new chromosome */
            if (args->new_chrom) {
                /* Non-target good hits for previous chromosome */
                if (args->prev_chrom_idx >= 0) {
                    if (args->do_cov_mask) {
                        handle_coverage_mask(args->coverage, args->cov_mask_ti,
                                             args->prev_chrom_idx,
                                             args->prev_chrom_len);
                    }

                    if (args->do_wgs) {
                        handle_wgs_coverage(args->coverage, args->cm_wgs,
                                            args->ci_wgs, args->prev_chrom_len);
                    }

                    if (args->do_capture) {
                        handle_target_coverage(
                            args->coverage, args->cm, args->ci, args->ti,
                            args->prev_chrom_idx, args->prev_chrom_name,
                            args->prev_chrom_len);
                        /*
                         * Note: target and buffer regions in coverage are
                         * wiped out by handle_miss_reads(). Only call
                         * handle_miss_reads() once coverage no longer
                         * needs to be processed.
                         */
                        handle_miss_reads(args->coverage, args->cm, args->ti,
                                          args->prev_chrom_idx,
                                          args->prev_chrom_len);
                    }
                }

                /* Next chromosome id, if available */
                if (args->curr_chrom_idx >= 0) {
                    /* Reset entire coverage array */
                    clear_coverage(args->coverage, 0, args->curr_chrom_len,
                                   args->curr_chrom_len);

                    if (args->do_wgs) {
                        args->cm_wgs->b_total += args->curr_chrom_len;
                    }
                    if (args->do_capture) {
                        args->cm->b_total += args->curr_chrom_len;
                        set_target_cov(args->target_cov, args->cm, args->ti,
                                       args->curr_chrom_idx,
                                       args->curr_chrom_len);

                        if (args->do_cov_mask) {
                            handle_coverage_mask_target(
                                args->target_cov, args->cm, args->cov_mask_ti,
                                args->curr_chrom_idx, args->curr_chrom_len);
                        }
                    }
                }
                args->new_chrom = false;
            }

            if (!is_read_filtered)
            {
                /* Process read for target and coverage info */
                capture_process_record(rec, args->coverage, args->target_cov,
                                       args->cm_wgs, args->cm, args->curr_chrom_len,
                                       args->remove_dups);
            }
        }

        if (args->verbose && ++args->num_records_processed % args->interval == 0) {
            log_info("%lu records processed...", args->num_records_processed);
        }
    }

    return record_idx;
}

/**
 * Handle coverage operations for last processed chromosome and write final
 * report to file.
 */
void finalize_results(args_t *args)
{
    /* Finalize results */
    if (args->do_alignment) {
        align_len_finalize(args->alm_all);
        align_len_finalize(args->alm_read1);
        align_len_finalize(args->alm_read2);
        /*align_len_finalize(args->alm_fragment);*/
        insert_size_finalize(args->ism);
    }

    /* Target coverage and miss reads for last chromosome */
    if ((args->do_wgs || args->do_capture) && args->curr_chrom_idx >= 0) {
        if (args->do_cov_mask) {
            handle_coverage_mask(args->coverage, args->cov_mask_ti,
                                 args->curr_chrom_idx, args->curr_chrom_len);
        }

        if (args->do_wgs) {
            handle_wgs_coverage(args->coverage, args->cm_wgs, args->ci_wgs,
                                args->curr_chrom_len);
        }

        if (args->do_capture) {
            handle_target_coverage(args->coverage, args->cm, args->ci, args->ti,
                                   args->curr_chrom_idx, args->curr_chrom_name,
                                   args->curr_chrom_len);
            handle_miss_reads(args->coverage, args->cm, args->ti,
                              args->curr_chrom_idx, args->curr_chrom_len);
        }
    }

    if (args->do_wgs) {
        capture_metrics_finalize(args->cm_wgs, args->ci_wgs, NULL);
    }
    if (args->do_capture) {
        capture_metrics_finalize(args->cm, args->ci, args->ti);
    }

    if (args->verbose) {
        log_info("Finished processing records.");
    }

    /* Write report */

    if (args->verbose) {
        log_info("Writing report.");
    }

    /* Report filtered/unfiltered read counts */
    fprintf(args->output_fp, "## AlignStats Report ##\n");
    filter_counter_report(args->fc_report, args->fc);
    report_print(args->output_fp, args->fc_report);

    /* Write alignment stats report */
    if (args->do_alignment) {
        align_report(args->a_report, args->am_all, RT_ALL);
        align_report(args->a_report, args->am_read1, RT_READ1);
        align_report(args->a_report, args->am_read2, RT_READ2);
        align_len_report(args->a_report, args->alm_all, RT_ALL);
        align_len_report(args->a_report, args->alm_read1, RT_READ1);
        align_len_report(args->a_report, args->alm_read2, RT_READ2);
        pair_stats_report(args->a_report, args->psm);
        insert_size_report(args->a_report, args->ism);

        fprintf(args->output_fp, "## Alignment Statistics ##\n");
        report_print(args->output_fp, args->a_report);
    }

    /* Write whole genome and/or coverage stats reports */
    if (args->do_capture || args->do_wgs) {
        fprintf(args->output_fp, "## Coverage Statistics ##\n");
    }

    /* Write whole genome stats report */
    if (args->do_wgs) {
        capture_report(args->wgs_report, args->cm_wgs, NULL);
        fprintf(args->output_fp, "### Whole Genome Metrics ###\n");
        report_print(args->output_fp, args->wgs_report);
    }

    /* Write capture stats report */
    if (args->do_capture) {
        capture_report(args->cap_report, args->cm, args->ti);
        fprintf(args->output_fp, "### Capture Metrics ###\n");
        report_print(args->output_fp, args->cap_report);
    }
}

#ifdef USE_PTHREAD
/**
 * Pthread barrier-synced function to read bam1_t records into buffers.
 */
void *pt_read_bam(void *arg)
{
    args_t *args = (args_t *)arg;
    uint32_t record_idx = 0;

    do {
        record_idx = args->read_bam_func(args);

        /* sync point 1 */
        pthread_barrier_wait(&args->barrier1);

        /* set read buff size for processing */
        args->read_buff_size = record_idx;

        /* swap buffer in use */
        args->read_buff = (args->read_buff == args->rec_buff_arr[0])
                              ? args->rec_buff_arr[1]
                              : args->rec_buff_arr[0];

        /* sync point 2 */
        pthread_barrier_wait(&args->barrier2);
    } while (record_idx > 0);

    pthread_exit(NULL);
}

/**
 * Pthread barrier-synced function to process bam1_t records read into buffers.
 */
void *pt_process_records(void *arg)
{
    bool more_to_process = true;
    args_t *args = (args_t *)arg;

    while (more_to_process) {
        process_records(args);

        /* sync point 1 */
        pthread_barrier_wait(&args->barrier1);

        /* swap buffer in use */
        args->curr_buff = (args->curr_buff == args->rec_buff_arr[0])
                              ? args->rec_buff_arr[1]
                              : args->rec_buff_arr[0];

        /* sync point 2 */
        pthread_barrier_wait(&args->barrier2);

        more_to_process = (args->read_buff_size > 0);
    }

    finalize_results(args);

    pthread_exit(NULL);
}
#endif /* USE_PTHREAD */

/**
 * Read bam1_t records into buffers and process.
 */
void read_and_process(args_t *args)
{
    args->read_buff = args->curr_buff = args->rec_buff_arr[0];

    while ((args->read_buff_size = args->read_bam_func(args)) > 0) {
        process_records(args);
        args->curr_buff = args->read_buff =
            (args->read_buff == args->rec_buff_arr[0])
                ? args->rec_buff_arr[1]
                : args->rec_buff_arr[0];
    }

    finalize_results(args);
}
