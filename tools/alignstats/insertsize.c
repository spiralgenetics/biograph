#include "insertsize.h"
#include "err.h"
#include "print.h"
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

/* Insert size metrics structure */

/**
 * Create new *insert_size_metrics_t and initialize values
 */
insert_size_metrics_t *insert_size_metrics_init()
{
    insert_size_metrics_t *ism = calloc(1, sizeof(insert_size_metrics_t));
    die_on_alloc_fail(ism);

    ism->insert_size_map = tree_map_init();
    ism->filter = BAM_FREAD1 |
                  BAM_FSECONDARY |
                  BAM_FUNMAP |
                  BAM_FMUNMAP |
                  BAM_FDUP;

    return ism;
}

/**
 * Free *ism from memory.
 */
void insert_size_metrics_destroy(insert_size_metrics_t *ism)
{
    if (ism != NULL) {
        tree_map_destroy(ism->insert_size_map);
        free(ism);
    }
}

/* Insert size metrics calculation */

/**
 * Process a bam1_t record for insert size metrics.
 */
void insert_size_process_record(bam1_t *rec, insert_size_metrics_t *ism)
{
    int32_t key;
    tree_node_t *node;

    if (rec->core.flag & BAM_FPAIRED &&    /* reads not in pair */
        !(rec->core.flag & ism->filter) && /* reads with filtered flags */
        rec->core.tid == rec->core.mtid)   /* read and mate not on same chr */
    {
        key = abs(rec->core.isize);
        node = tree_map_get(ism->insert_size_map, key);

        /* Increment insert size in map */
        tree_map_set(ism->insert_size_map, key,
                     (node == NULL) ? 1 : node->value + 1);
    }
}

/**
 * Finalize insert size metrics once all records are processed.
 * Calculates mean, median, and mode insert sizes.
 */
void insert_size_finalize(insert_size_metrics_t *ism)
{
    tree_node_key_t *keyset;
    tree_node_t *node;
    uint64_t num_sizes, sum_sizes, mode_size, curr_size, median_idx;

    if (tree_map_set_keyset(&keyset, ism->insert_size_map)) {
        num_sizes = sum_sizes = mode_size = 0;

        /* Mean and mode alignment sizes */
        for (size_t i = 0; i < ism->insert_size_map->num_nodes; ++i) {
            if (!tree_map_set_node(&node, ism->insert_size_map, keyset[i])) {
                goto fail;
            }

            curr_size = (uint64_t)node->value;
            if (curr_size > mode_size) {
                mode_size = curr_size;
                ism->mode = (uint64_t)keyset[i];
            }

            num_sizes += curr_size;
            sum_sizes += (uint64_t)keyset[i] * curr_size;
        }

        if (num_sizes != 0) {
            ism->mean = (double)sum_sizes / (double)num_sizes;

            /* Median alignment size */
            median_idx = num_sizes / 2;
            num_sizes = 0;

            for (size_t i = 0; i < ism->insert_size_map->num_nodes; ++i) {
                if (!tree_map_set_node(&node, ism->insert_size_map, keyset[i])) {
                    goto fail;
                }

                if ((num_sizes += (uint64_t)node->value) >= median_idx) {
                    ism->median = (uint64_t)keyset[i];
                    break;
                }
            }
        } else {
            ism->mean = 0.0;
            ism->median = 0;
        }

fail:
        free(keyset);
    }
}

/**
 * Write metrics in ism to report.
 */
void insert_size_report(report_t *report, insert_size_metrics_t *ism)
{
    char *key_buffer = malloc(REPORT_BUFFER_SIZE * sizeof(char));
    die_on_alloc_fail(key_buffer);
    char *value_buffer = malloc(REPORT_BUFFER_SIZE * sizeof(char));
    die_on_alloc_fail(value_buffer);

    copy_to_buffer(key_buffer, "Mean_Insert_Size", REPORT_BUFFER_SIZE);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%.2f", ism->mean);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_buffer, "Median_Insert_Size", REPORT_BUFFER_SIZE);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", ism->median);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_buffer, "Mode_Insert_Size", REPORT_BUFFER_SIZE);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", ism->mode);
    report_add_key_value(report, key_buffer, value_buffer);

    free(key_buffer);
    free(value_buffer);
}
