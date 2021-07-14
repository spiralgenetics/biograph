#ifndef _BED_H
#define _BED_H

#include "htslib/sam.h"
#include "report.h"
#include <stdbool.h>
#include <stdint.h>

#define CHROM_BUFFER_SIZE 4096
#define CHROM_BUFFER_SIZE_SCAN "%4095s"

/* Bed structure */
struct bed_chrom {
    size_t num_targets;  /* Number of targets read from file */
    int32_t *start_pos;  /* Start positions array */
    int32_t *end_pos;    /* End positions array */
};
typedef struct bed_chrom bed_chrom_t;

bed_chrom_t *bed_chrom_init(size_t count);
void bed_chrom_destroy(bed_chrom_t *tic);

struct bed {
    char **chrom_names;
    size_t num_targets;
    size_t num_chroms;
    bed_chrom_t **chroms;
};
typedef struct bed bed_t;

bed_t *bed_init();
uint32_t bed_sum_bases(bed_t *ti);
void bed_destroy(bed_t *ti);
size_t load_bed(FILE *fp, bed_t *ti, bam_hdr_t *hdr, bool zero_based);

#endif /* _BED_H */
