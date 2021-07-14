#ifndef _PROCESSING_H
#define _PROCESSING_H

#include "alignstats.h"

void set_iter(args_t *args);
bool move_to_first_region(args_t *args);
bool move_to_next_region(args_t *args);
uint32_t read_bam1(args_t *args);
uint32_t read_bam_itr(args_t *args);
uint32_t process_records(args_t *args);
void finalize_results(args_t *args);
#ifdef USE_PTHREAD
void *pt_read_bam(void *arg);
void *pt_process_records(void *arg);
#endif
void read_and_process(args_t *args);
void clean_up_files(args_t *args);

#endif /* _PROCESSING_H */
