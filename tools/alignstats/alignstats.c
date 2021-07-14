#include "alignstats.h"
#include "err.h"
#include "logging.h"
#include "processing.h"
#include "report.h"
#include "treemap.h"
#include <ctype.h>
#include <getopt.h>
#include <stdlib.h>
#include <string.h>

/**
 * ------===< AlignStats >===------
 * Produce a report of alignment and coverage metrics for an input SAM, BAM, or CRAM file.
 *
 * Copyright 2015-2017 Baylor College of Medicine Human Genome Sequencing Center
 *
 * Author: Jesse Farek (farek at bcm dot edu)
 * License: BSD 3-clause (see LICENSE or https://opensource.org/licenses/BSD-3-Clause)
 */

/**
 * Todo list/Wish list
 *
 * long options, should be easy
 * bed support using htslib regidx
 * option for multithreading by chunking from read buffer
 * target_cov[] can just be made the highest 2 or so bits of coverage[]
 */

void usage()
{
    fprintf(stderr, "Usage: alignstats [-i INPUT] [-j FORMAT] [-o OUTPUT]\n");
    fprintf(stderr, "                  [-h] [-v] [-n NUMREADS] [-p] [-z]\n");
    fprintf(stderr, "                  [-r REGIONS] [-t TARGET] [-m COVMASK] [-T REFFASTA]\n");
    fprintf(stderr, "                  [-q INT] [-f INT] [-F INT]\n");
    fprintf(stderr, "                  [-D] [-U] [-A] [-C] [-W]\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Runtime options:\n");
    fprintf(stderr, "    -h          Print usage information.\n");
    fprintf(stderr, "    -v          Print verbose runtime information output to stderr.\n");
    fprintf(stderr, "    -n INT      Maximum number of records to keep in memory.\n");
    fprintf(stderr, "    -p          Use separate threads for reading and processing records\n");
    fprintf(stderr, "                (requires builtin pthread support).\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "File options:\n");
    fprintf(stderr, "    -i INPUT    Read INPUT as the input SAM, BAM, or CRAM file (stdin). Input\n");
    fprintf(stderr, "                must be coordinate-sorted for accurate results.\n");
    fprintf(stderr, "    -j FORMAT   Specify file format of input alignment file (\"sam\", \"bam\", or\n");
    fprintf(stderr, "                \"cram\" available, default guessed from filename or \"sam\").\n");
    fprintf(stderr, "    -o OUTPUT   Write report to OUTPUT (stdout).\n");
    fprintf(stderr, "    -r REGIONS  File in BED format listing which regions to process. By\n");
    fprintf(stderr, "                default, all available records are processed. This option\n");
    fprintf(stderr, "                requires the alignment file to be indexed.\n");
    fprintf(stderr, "    -t TARGET   File in BED format listing capture coverage regions. Required\n");
    fprintf(stderr, "                if capture coverage statistics are enabled.\n");
    fprintf(stderr, "    -m COVMASK  File in BED format listing regions of N bases in reference.\n");
    fprintf(stderr, "                Coverage counts will be suppressed for these regions.\n");
    fprintf(stderr, "    -T REFFASTA Indexed FASTA reference file for CRAM input alignment.\n");
    fprintf(stderr, "    -z           Translate BED files from 1-based to 0-based format\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Processing options:\n");
    fprintf(stderr, "    -q INT      Only process records with minimum read quality of INT.\n");
    fprintf(stderr, "    -f INT      Only process records with all bits in INT set in FLAG.\n");
    fprintf(stderr, "    -F INT      Only process records with none of bits in INT set in FLAG.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Reporting options:\n");
    fprintf(stderr, "    -D          Disable excluding duplicate reads from coverage statistics.\n");
    fprintf(stderr, "    -U          Disable processing unplaced unmapped reads (CHROM \"*\") when\n");
    fprintf(stderr, "                using the -r option.\n");
    fprintf(stderr, "    -A          Disable reporting alignment statistics.\n");
    fprintf(stderr, "    -C          Disable reporting capture coverage statistics.\n");
    fprintf(stderr, "    -W          Disable reporting whole genome coverage statistics.\n");
}

void header()
{
    fprintf(stderr, "AlignStats v" ALIGNSTATS_VERSION " (using HTSlib %s)\n", hts_version());
}

enum input_format { INPUT_SAM = 0, INPUT_BAM, INPUT_CRAM };

#define MODE_LEN 3
#define FMT_LEN  5

int main(int argc, char **argv)
{
    bool input_from_stdin, rec_buff_allocated;
    char *input_fn, *output_fn, *target_fn, *regions_fn;
    char *cov_mask_fn, *reference_fn, *ref_buff, *end;
    char mode[MODE_LEN], fmt[FMT_LEN];
    int c, exit_val;
    uint8_t min_qual;
    uint16_t filter_incl, filter_excl;
    uint32_t max_chrom_len, min_buffer_reads, max_reads_tmp;
    size_t fn_len, rec_buff_size;
    void *status;
    FILE *target_fp, *cov_mask_fp, *regions_fp;
    enum input_format format;
    const uint16_t FLAG_MAX = 4095; /* == (1<<12)-1 */

    args_t *args;
#ifdef USE_PTHREAD
    pthread_t *read_bam_thread, *process_records_thread;
#endif

    exit_val = EXIT_SUCCESS;
    input_from_stdin = false;
    rec_buff_allocated = false;
    input_fn = NULL;
    output_fn = NULL;
    target_fn = NULL;
    cov_mask_fn = NULL;
    regions_fn = NULL;
    reference_fn = NULL;
    target_fp = NULL;
    cov_mask_fp = NULL;
    regions_fp = NULL;
    ref_buff = NULL;
    min_buffer_reads = 200;
    memset(mode, '\0', MODE_LEN);
    memset(fmt, '\0', FMT_LEN);
    format = INPUT_SAM;

    args = calloc(1, sizeof(args_t));
    die_on_alloc_fail(args);

    args->verbose = false;
    args->do_alignment = true;
    args->do_capture = true;
    args->do_wgs = true;
    args->do_cov_mask = false;
    args->do_pthread = false;
    args->remove_dups = true;
    args->process_unmapped = true;
    args->process_unmapped_done = false;
    args->zero_based = false;

    args->input_sf = NULL;
    args->output_fp = NULL;

    args->prev_chrom_name = NULL;
    args->curr_chrom_name = NULL;
    args->prev_chrom_idx = -999;
    args->prev_mapped_chrom_idx = -999;
    args->curr_chrom_idx = -999;
    args->prev_chrom_len = 0;
    args->curr_chrom_len = 0;
    args->num_records_processed = 0;
    args->interval = 10000000;
    args->new_chrom = true;
    args->process_cigar = false;
    args->order_warn = true;
    args->reads_per_buffer = RECORD_BUFFER_SIZE;

    args->iter = NULL;
    args->index = NULL;
    args->regions = NULL;
    args->regions_curr_chrom_idx = 0;
    args->regions_curr_target_idx = 0;
    args->read_bam_func = NULL;

    args->am_all = NULL;
    args->am_read1 = NULL;
    args->am_read2 = NULL;
    /*args->am_fragment = NULL;*/
    args->alm_all = NULL;
    args->alm_read1 = NULL;
    args->alm_read2 = NULL;
    /*args->alm_fragment = NULL;*/
    args->psm = NULL;
    args->ism = NULL;
    args->cm = NULL;
    args->cm_wgs = NULL;
    args->ci = NULL;
    args->ci_wgs = NULL;
    args->ti = NULL;
    args->cov_mask_ti = NULL;
    args->coverage = NULL;
    args->target_cov = NULL;
    args->a_report = NULL;
    args->wgs_report = NULL;
    args->cap_report = NULL;

    /* No parameters */
    if (argc == 1) {
        header();
        usage();
        exit_val = EXIT_FAILURE;
        goto end;
    }

    min_qual = 0;
    filter_incl = 0;
    filter_excl = 0;

    /* Read parameters */
    while ((c = getopt(argc, argv, "ACDF:T:UWf:hi:j:m:n:o:pq:r:t:vz")) != -1) {
        switch (c) {
        case 'A': /* Turn off alignment stats */
            args->do_alignment = false;
            break;
        case 'C': /* Turn off capture stats */
            args->do_capture = false;
            break;
        case 'D': /* Don't remove duplicate reads */
            args->remove_dups = false;
            break;
        case 'F': /* Filter exclude */
            filter_excl = (uint16_t)strtol(optarg, &end, 10);
            if (filter_excl > FLAG_MAX) {
                log_warning("Invalid flag received from -F option. Setting to 0.");
                filter_excl = 0;
            }
            break;
        case 'T': /* Reference filename for cram */
            reference_fn = optarg;
            break;
        case 'U': /* Don't process '*' chrom */
            args->process_unmapped = false;
            break;
        case 'W': /* Turn off whole genome stats */
            args->do_wgs = false;
            break;
        case 'f': /* Filter include */
            filter_incl = (uint16_t)strtol(optarg, &end, 10);
            if (filter_incl > FLAG_MAX) {
                log_warning("Invalid flag received from -f option. Setting to 0.");
                filter_incl = 0;
            }
            break;
        case 'h': /* Usage */
            header();
            usage();
            exit_val = EXIT_SUCCESS;
            goto end;
            break;
        case 'i': /* Input filename */
            input_fn = optarg;
            break;
        case 'j': /* Input file format */
            strncpy(fmt, optarg, FMT_LEN);
            break;
        case 'm': /* Coverage mask filename */
            cov_mask_fn = optarg;
            break;
        case 'n': /* Maximum number of reads in memory */
            max_reads_tmp = (uint32_t)strtol(optarg, &end, 10);
            if (max_reads_tmp < min_buffer_reads) {
                log_warning("Given value for -r is too small, using %u", min_buffer_reads);
                max_reads_tmp = min_buffer_reads;
            }
            args->reads_per_buffer = max_reads_tmp / 2;
            break;
        case 'o': /* Output filename */
            output_fn = optarg;
            break;
        case 'p': /* Read input in separate thread */
            args->do_pthread = true;
            break;
        case 'q': /* Minimum read quality */
            min_qual = (uint8_t)strtol(optarg, &end, 10);
            break;
        case 'r': /* Regions */
            regions_fn = optarg;
            break;
        case 't': /* Target filename */
            target_fn = optarg;
            break;
        case 'z': /* Translate bed to zero-based */
            args->zero_based = true;
            break;
        case 'v': /* Verbose runtime information */
            args->verbose = true;
            break;
        default:
            header();
            usage();
            exit_val = EXIT_FAILURE;
            goto end;
            break;
        }
    }

    if (args->verbose) {
        log_info("Running AlignStats v" ALIGNSTATS_VERSION);
        log_info("Command:");
        for (int i = 0; i < argc; ++i) {
            log_info_r((i == 0 ? "%s" : " %s"), argv[i]);
        }
        log_info_r("\n");
    }

    /* Check parameters */

    if (!args->do_alignment && !args->do_capture && !args->do_wgs) {
        log_error("At least one of -A, -C, or -W must be unset.");
        exit_val = EXIT_FAILURE;
        goto end;
    }

    if (input_fn == NULL) {
        input_fn = "-";
    }

    if (strcmp(input_fn, "-") == 0) {
        input_from_stdin = true;
    }

    /* Open files */

    /* Format specified with -j option */
    if (*fmt != '\0') {
        for (size_t i = 0; fmt[i] != '\0'; ++i) {
            fmt[i] = tolower(fmt[i]);
        }

        if (strcmp(fmt, "bam") == 0) {
            format = INPUT_BAM;
        } else if (strcmp(fmt, "cram") == 0) {
            format = INPUT_CRAM;
        } else if (strcmp(fmt, "sam") != 0) {
            log_warning("Unrecognized input format received from -j option.");
        }
    /* Guess format from filename */
    } else if (!input_from_stdin) {
        fn_len = strlen(input_fn);

        if (fn_len >= 3 && strcmp(input_fn + fn_len - 3, "bam") == 0) {
            format = INPUT_BAM;
        } else if (fn_len >= 4 && strcmp(input_fn + fn_len - 4, "cram") == 0) {
            format = INPUT_CRAM;
        } else if (!(fn_len >= 3 && strcmp(input_fn + fn_len - 3, "sam") == 0)) {
            log_warning("Input filename has an unrecognized file extension.");
        }
    }

    switch (format) {
    case INPUT_BAM:
        strncpy(mode, "rb", 3);
        if (args->verbose) {
            log_info("Opening input file as BAM file.");
        }
        break;
    case INPUT_CRAM:
        strncpy(mode, "rc", 3);
        if (args->verbose) {
            log_info("Opening input file as CRAM file.");
        }
        break;
    case INPUT_SAM:
    default:
        strncpy(mode, "r", 2);
        if (args->verbose) {
            log_info("Opening input file as SAM file.");
        }
        break;
    }

    if (format == INPUT_CRAM && reference_fn == NULL) {
        log_error("No reference specified for input CRAM alignment file.");
        exit_val = EXIT_FAILURE;
        goto end;
    }

    /* Open input file */
    if ((args->input_sf = sam_open(input_fn, mode)) == NULL) {
        log_error("Failed to open input file \"%s\".", input_fn);
        perror(NULL);
        exit_val = EXIT_FAILURE;
        goto end;
    }

    /* Set reference .fa.fai for cram input */
    if (format == INPUT_CRAM) {
        rec_buff_size = strlen(reference_fn) + 5;
        ref_buff = calloc(rec_buff_size, sizeof(char));
        die_on_alloc_fail(ref_buff);

        strncpy(ref_buff, reference_fn, rec_buff_size);
        strncat(ref_buff, ".fai", strlen(".fai"));

        if (hts_set_fai_filename(args->input_sf, ref_buff) != 0) {
            log_error("hts_set_fai_filename() failed for \"%s\".", ref_buff);
            perror(NULL);
            exit_val = EXIT_FAILURE;
            goto end;
        }

        free(ref_buff);
    }

    if ((args->hdr = sam_hdr_read(args->input_sf)) == NULL) {
        log_error("Failed to read header for input file \"%s\".", input_fn);
        perror(NULL);
        exit_val = EXIT_FAILURE;
        goto end;
    }

    /* Alignment index file if processing by region */
    if (regions_fn != NULL) {
        if (input_from_stdin) {
            log_error("Cannot process regions on input alignment from stdin.");
            exit_val = EXIT_FAILURE;
            goto end;
        } else if ((args->index = sam_index_load(args->input_sf, input_fn)) == NULL) {
            log_error("Failed to read index for input file \"%s\".", input_fn);
            perror(NULL);
            exit_val = EXIT_FAILURE;
            goto end;
        }
    }

    /* Target file for capture */
    if (args->do_capture) {
        if (args->verbose) {
            log_info("Opening target file.");
        }

        /* Open target file */
        if (target_fn != NULL) {
            if ((target_fp = fopen(target_fn, "r")) == NULL) {
                log_error("Failed to open target file \"%s\".", target_fn);
                perror(NULL);
                exit_val = EXIT_FAILURE;
                goto end;
            }
        } else {
            log_error("No target file specified for capture statistics.");
            exit_val = EXIT_FAILURE;
            goto end;
        }
    }

    /* Output file */
    if (output_fn == NULL) {
        log_warning("No output file specified, using stdout for report.");
        args->output_fp = stdout;
    } else {
        if (args->verbose) {
            log_info("Opening output file.");
        }
        if ((args->output_fp = fopen(output_fn, "w")) == NULL) {
            log_error("Failed to open output file \"%s\".", output_fn);
            perror(NULL);
            exit_val = EXIT_FAILURE;
            goto end;
        }
    }

    /* N bases coverage values mask */
    if (cov_mask_fn != NULL) {
        if (args->verbose) {
            log_info("%s", "Opening coverage mask file.");
        }
        if ((cov_mask_fp = fopen(cov_mask_fn, "r")) == NULL) {
            log_error("Failed to open coverage mask file \"%s\".", cov_mask_fn);
            perror(NULL);
            exit_val = EXIT_FAILURE;
            goto end;
        }

        args->do_cov_mask = true;
    }

    /* Region BED file */
    if (regions_fn != NULL) {
        if (args->verbose) {
            log_info("Opening regions file.");
        }
        if ((regions_fp = fopen(regions_fn, "r")) == NULL) {
            log_error("Failed to open regions file \"%s\".", regions_fn);
            perror(NULL);
            exit_val = EXIT_FAILURE;
            goto end;
        }
    }

    if (args->verbose) {
        log_info("Files opened successfully. Files in use:");
        log_info("Input file: %s", input_fn);
        if (target_fp != NULL) {
            log_info("Target file: %s", target_fn);
        }
        log_info("Output file: %s",
                 (args->output_fp == stdout) ? "(stdout)" : output_fn);
        if (regions_fn != NULL) {
            log_info("Regions file: %s", regions_fn);
        }
        if (cov_mask_fn != NULL) {
            log_info("Coverage mask file: %s", cov_mask_fn);
        }
    }

    /* Prepare all the data structures */

    args->fc = filter_counter_init(min_qual, filter_incl, filter_excl);

    /* Create metrics calculators */
    if (args->do_alignment) {
        args->am_all = align_metrics_init();
        args->am_read1 = align_metrics_init();
        args->am_read2 = align_metrics_init();
        /*args->am_fragment = align_metrics_init();*/
        args->alm_all = align_len_metrics_init();
        args->alm_read1 = align_len_metrics_init();
        args->alm_read2 = align_len_metrics_init();
        /*args->alm_fragment = align_len_metrics_init();*/
        args->psm = pair_stats_metrics_init();
        args->ism = insert_size_metrics_init();
    }

    if (regions_fn != NULL) {
        if (args->verbose) {
            log_info("Loading regions file...");
        }

        args->regions = bed_init();
        if (load_bed(regions_fp, args->regions, args->hdr, args->zero_based) == 0) {
            log_error("No usable regions in regions file.");
            exit_val = EXIT_FAILURE;
            goto end;
        }
        args->regions_curr_chrom_idx = 0;
        args->regions_curr_target_idx = 0;

        if (move_to_first_region(args)) {
            if (args->verbose) {
                log_info("Regions file loaded successfully.");
            }
        } else {
            log_error("No usable regions in regions file.");
            exit_val = EXIT_FAILURE;
            goto end;
        }

        if (fclose(regions_fp) != 0) {
            log_error("Error closing regions file \"%s\".", regions_fn);
            perror(NULL);
        }
    }

    /* Load targets for capture stats */
    if (args->do_capture) {
        /* Setup capture stats metrics */
        args->ti = bed_init();
        args->ci = coverage_info_init();

        if (args->verbose) {
            log_info("Loading targets...");
        }

        if (load_bed(target_fp, args->ti, args->hdr, args->zero_based) == 0) {
            log_error("No usable targets in target file.");
            exit_val = EXIT_FAILURE;
            goto end;
        } else {
            if (args->verbose) {
                log_info("Targets loaded successfully.");
            }
            /* Can't populate targets info until we've got it loaded */
            args->cm = capture_metrics_init(args->ti);
        }

        if (fclose(target_fp) != 0) {
            log_error("Error closing target file \"%s\".", target_fn);
            perror(NULL);
        }
    }

    /* Load coverage mask regions as targets */
    if (args->do_cov_mask) {
        args->cov_mask_ti = bed_init();

        if (args->verbose) {
            log_info("Loading coverage mask targets...");
        }

        if (load_bed(cov_mask_fp, args->cov_mask_ti, args->hdr, args->zero_based) == 0) {
            log_warning("No usable targets in coverage mask file.");
            exit_val = EXIT_FAILURE;
            goto end;
        } else if (args->verbose) {
            log_info("Coverage mask targets loaded successfully.");
        }

        if (fclose(cov_mask_fp) != 0) {
            log_error("Error closing coverage mask file \"%s\".", cov_mask_fn);
            perror(NULL);
        }
    }

    if (args->do_wgs) {
        //blank bed...?
        args->cm_wgs = capture_metrics_init(bed_init());
        args->cm_wgs->b_masked = bed_sum_bases(args->cov_mask_ti);
        args->ci_wgs = coverage_info_init();
    }

    args->read_bam_func = (regions_fn == NULL) ? read_bam1 : read_bam_itr;

    /* Reports */
    args->fc_report = report_init();
    if (args->do_alignment) {
        args->a_report = report_init();
    }
    if (args->do_wgs) {
        args->wgs_report = report_init();
    }
    if (args->do_capture) {
        args->cap_report = report_init();
    }

    /* Set sizes of coverage and target_cov to size of largest chromosome */
    max_chrom_len = 0;

    for (int32_t i = 0; i < args->hdr->n_targets; ++i) {
        if (args->hdr->target_len[i] > max_chrom_len) {
            max_chrom_len = args->hdr->target_len[i];
        }
    }

    if (max_chrom_len == 0) {
        log_error("No usable reference sequences in SAM header.");
        exit_val = EXIT_FAILURE;
        goto end;
    }

    /* coverage[] remains NULL otherwise */
    if (args->do_wgs || args->do_capture) {
        args->coverage = calloc(max_chrom_len, sizeof(uint32_t));
        die_on_alloc_fail(args->coverage);
    }

    /* target_cov[] remains NULL otherwise */
    if (args->do_capture) {
        args->target_cov = calloc(max_chrom_len, sizeof(uint8_t));
        die_on_alloc_fail(args->target_cov);
    }

    args->rec_buff_arr[0] = malloc(args->reads_per_buffer * sizeof(bam1_t *));
    die_on_alloc_fail(args->rec_buff_arr[0]);
    args->rec_buff_arr[1] = malloc(args->reads_per_buffer * sizeof(bam1_t *));
    die_on_alloc_fail(args->rec_buff_arr[1]);
    rec_buff_allocated = true;

    for (uint32_t i = 0; i < args->reads_per_buffer; ++i) {
        args->rec_buff_arr[0][i] = bam_init1();
        args->rec_buff_arr[1][i] = bam_init1();
    }

    args->read_buff = args->rec_buff_arr[0];
    args->curr_buff = args->rec_buff_arr[1];
    args->read_buff_size = 0;

    if (args->do_pthread) {
#ifdef USE_PTHREAD
        pthread_barrier_init(&args->barrier1, NULL, 2);
        pthread_barrier_init(&args->barrier2, NULL, 2);

        read_bam_thread = calloc(1, sizeof(pthread_t));
        die_on_alloc_fail(read_bam_thread);
        process_records_thread = calloc(1, sizeof(pthread_t));
        die_on_alloc_fail(process_records_thread);

        /* Create and start the reading and processing threads */
        if (pthread_create(read_bam_thread, NULL, pt_read_bam, args)) {
            log_error("Error creating read_bam_thread.");
            exit_val = EXIT_FAILURE;
            goto end;
        }

        if (pthread_create(process_records_thread, NULL, pt_process_records, args)) {
            log_error("Error creating process_records_thread.");
            exit_val = EXIT_FAILURE;
            goto end;
        }

        if (args->verbose) {
            log_info("Start processing records.");
        }

        /* Join main thread with the two threads */
        if (pthread_join(*read_bam_thread, &status) != 0) {
            log_error("Error joining read_bam_thread.");
            exit_val = EXIT_FAILURE;
            goto end;
        }

        if (pthread_join(*process_records_thread, &status) != 0) {
            log_error("Error joining process_records_thread.");
            exit_val = EXIT_FAILURE;
            goto end;
        }

        free(read_bam_thread);
        free(process_records_thread);

        pthread_barrier_destroy(&args->barrier1);
        pthread_barrier_destroy(&args->barrier2);
    } else {
#else
        log_warning("AlignStats not built with pthread, multithreading disabled.");
    }
#endif
        if (args->verbose) {
            log_info("Start processing records.");
        }

        read_and_process(args);
#ifdef USE_PTHREAD
    }
#endif

    /* Clean up */
end:
    bam_hdr_destroy(args->hdr);
    if (args->input_sf != NULL) {
        sam_close(args->input_sf);
    }

    if (rec_buff_allocated) {
        for (uint32_t i = 0; i < args->reads_per_buffer; ++i) {
            bam_destroy1(args->rec_buff_arr[0][i]);
            bam_destroy1(args->rec_buff_arr[1][i]);
        }
    }

    report_destroy(args->fc_report);
    report_destroy(args->a_report);
    report_destroy(args->wgs_report);
    report_destroy(args->cap_report);
    align_metrics_destroy(args->am_all);
    align_metrics_destroy(args->am_read1);
    align_metrics_destroy(args->am_read2);
    /*align_metrics_destroy(args->am_fragment);*/
    align_len_metrics_destroy(args->alm_all);
    align_len_metrics_destroy(args->alm_read1);
    align_len_metrics_destroy(args->alm_read2);
    /*align_len_metrics_destroy(args->alm_fragment);*/
    pair_stats_metrics_destroy(args->psm);
    insert_size_metrics_destroy(args->ism);
    capture_metrics_destroy(args->cm);
    capture_metrics_destroy(args->cm_wgs);
    coverage_info_destroy(args->ci);
    coverage_info_destroy(args->ci_wgs);
    filter_counter_destroy(args->fc);
    bed_destroy(args->ti);
    bed_destroy(args->cov_mask_ti);
    free(args->coverage);
    free(args->target_cov);
    free(args->rec_buff_arr[0]);
    free(args->rec_buff_arr[1]);

    if (regions_fn != NULL) {
        bed_destroy(args->regions);
        hts_itr_destroy(args->iter);
        hts_idx_destroy(args->index);
    }
    if (args->output_fp != NULL && fclose(args->output_fp) != 0) {
        log_error("Error closing output file \"%s\".", output_fn);
        perror(NULL);
    }
    if (args->verbose && exit_val == EXIT_SUCCESS) {
        log_info("Finished successfully.");
    }

    free(args);

    return exit_val;
}
