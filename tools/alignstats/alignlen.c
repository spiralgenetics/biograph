#include "alignlen.h"
#include "err.h"
#include "print.h"
#include <string.h>

/* Alignment length metrics structure */

/**
 * Create and return new *align_len_t.
 */
align_len_metrics_t *align_len_metrics_init()
{
    align_len_metrics_t *alm = calloc(1, sizeof(align_len_metrics_t));
    die_on_alloc_fail(alm);

    alm->length_map = tree_map_init();

    return alm;
}

/**
 * Free *align_len_t from memory.
 */
void align_len_metrics_destroy(align_len_metrics_t *alm)
{
    if (alm != NULL) {
        tree_map_destroy(alm->length_map);
        free(alm);
    }
}

/* Alignment length metrics calculation */

/**
 * Process a bam1_t record for alignment length metrics.
 */
void align_len_process_record(bam1_t *rec, align_len_metrics_t *alm)
{
    int32_t key, num_s;
    uint32_t *cigar;
    tree_node_t *node;

    if (!(rec->core.flag & BAM_FUNMAP)) {
        num_s = 0;
        cigar = bam_get_cigar(rec);

        /* Count number of CIGAR S bases in record */
        for (uint16_t i = 0; i < rec->core.n_cigar; ++i) {
            if (bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) {
                num_s += bam_cigar_oplen(cigar[i]);
            }
        }

        /* Increment alignment length count */
        key = rec->core.l_qseq - num_s;
        node = tree_map_get(alm->length_map, key);

        tree_map_set(alm->length_map, key, (node == NULL) ? 1 : node->value + 1);
    }
}

/**
 * Finalize alignment length metrics once all records are processed.
 * Calculates mean, median, and mode alignment lengths.
 */
void align_len_finalize(align_len_metrics_t *alm)
{
    size_t i;
    uint64_t num_lengths, sum_lengths, mode_length, curr_length, median_idx;
    tree_node_key_t *keyset;
    tree_node_t *node;

    if (tree_map_set_keyset(&keyset, alm->length_map)) {
        num_lengths = sum_lengths = mode_length = 0;

        /* Mean and mode alignment lengths */
        for (i = 0; i < alm->length_map->num_nodes; ++i) {

            if (!tree_map_set_node(&node, alm->length_map, keyset[i])) {
                goto fail;
            }

            if ((curr_length = (uint64_t)node->value) > mode_length) {
                mode_length = curr_length;
                alm->mode = (uint64_t)keyset[i];
            }

            num_lengths += curr_length;
            sum_lengths += (uint64_t)keyset[i] * curr_length;
        }

        if (num_lengths != 0) {
            alm->mean = (double)sum_lengths / (double)num_lengths;

            /* Median alignment length */
            median_idx = num_lengths / 2;
            num_lengths = 0;

            for (i = 0; i < alm->length_map->num_nodes; ++i) {
                if (!tree_map_set_node(&node, alm->length_map, keyset[i])) {
                    goto fail;
                }

                if ((num_lengths += (uint64_t)node->value) >= median_idx) {
                    alm->median = (uint64_t)keyset[i];
                    break;
                }
            }
        } else {
            alm->mean = 0.0;
            alm->median = 0;
        }

fail:
        free(keyset);
    }
}

/**
 * Add alignment length metrics to report.
 */
void align_len_report(report_t *report, align_len_metrics_t *alm, read_type_t rt)
{
    char *key_buffer = malloc(REPORT_BUFFER_SIZE * sizeof(char));
    die_on_alloc_fail(key_buffer);
    char *value_buffer = malloc(REPORT_BUFFER_SIZE * sizeof(char));
    die_on_alloc_fail(value_buffer);

    char *key_start;
    size_t prefix_len, copy_size;

    /*
     * Add key prefix for read1 and read2
     * All reads: prefix = ""
     */
    const char *prefix = (rt == RT_READ1) ? "R1_" : rt == RT_READ2 ? "R2_" : "";
    prefix_len = strlen(prefix);
    copy_to_buffer(key_buffer, prefix, REPORT_BUFFER_SIZE);
    key_start = key_buffer + strlen(prefix);
    copy_size = REPORT_BUFFER_SIZE - prefix_len;

    copy_to_buffer(key_start, "Mean_Aligned_Read_Length", copy_size);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%.2f", alm->mean);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "Median_Aligned_Read_Length", copy_size);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", alm->median);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "Mode_Aligned_Read_Length", copy_size);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", alm->mode);
    report_add_key_value(report, key_buffer, value_buffer);

    free(key_buffer);
    free(value_buffer);
}
