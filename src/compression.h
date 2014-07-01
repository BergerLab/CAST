#ifndef __CABLAST_COMPRESSION_H__
#define __CABLAST_COMPRESSION_H__

/* Apparently this is required to make pthread_rwlock* stuff available. */
#define __USE_UNIX98

#include <pthread.h>
#include <stdint.h>

#include "ds.h"

#include "align.h"
#include "coarse.h"
#include "compressed.h"
#include "database.h"
#include "seeds.h"
#include "seq.h"

struct cb_compress_workers {
    pthread_t *threads;
    int32_t num_workers;
    struct DSQueue *jobs;
    void *args;
};

struct cb_compress_workers *cb_compress_start_workers(struct cb_database *db,
                                                      int32_t num_workers);

void cb_compress_join_workers(struct cb_compress_workers *workers);

void cb_compress_free_workers(struct cb_compress_workers *workers);

void cb_compress_send_job(struct cb_compress_workers *workers,
                          struct cb_seq *org_seq);

struct cb_compressed_seq *cb_compress(struct cb_coarse *coarse_db,
                                      struct cb_seq *org_seq,
                                      struct cb_align_nw_memory *mem,
                                      struct cb_seeds_add_memory *seeds_mem);

#endif
