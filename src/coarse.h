#ifndef __CABLAST_COARSE_H__
#define __CABLAST_COARSE_H__

/* Apparently this is required to make pthread_rwlock* stuff available. */
#define __USE_UNIX98

#include <pthread.h>
#include <stdint.h>
#include <stdio.h>

#include "ds.h"

#include "link_to_compressed.h"
#include "seeds.h"
#include "seq.h"
#include "stdbool.h"

struct cb_link_to_compressed *
cb_link_to_compressed_init(int32_t org_seq_id, int16_t coarse_start,
                            int16_t coarse_end, uint64_t original_start,
                            uint64_t original_end, bool dir);

void
cb_link_to_compressed_free(struct cb_link_to_compressed *link);

struct cb_coarse_seq {
    int32_t id;
    struct cb_seq *seq;
    struct cb_link_to_compressed *links;
    pthread_rwlock_t lock_links;
};

struct cb_coarse_seq *
cb_coarse_seq_init(int32_t id, char *residues, int32_t start, int32_t end);

void
cb_coarse_seq_free(struct cb_coarse_seq *seq);

void
cb_coarse_seq_addlink(struct cb_coarse_seq *seq,
                       struct cb_link_to_compressed *newlink);

struct cb_coarse {
    struct DSVector *seqs;
    struct cb_seeds *seeds;
    uint64_t dbsize;
    FILE *file_fasta;
    FILE *file_seeds;
    FILE *file_links;
    FILE *file_links_index;
    FILE *file_fasta_index;
    FILE *file_params;
    pthread_rwlock_t lock_seq;
};

struct cb_coarse *
cb_coarse_init(int32_t seed_size,
                FILE *file_fasta, FILE *file_seeds, FILE *file_links,
                FILE *file_links_index, FILE *file_fasta_index,
                FILE *file_params);

void
cb_coarse_free(struct cb_coarse *coarse_db);

struct cb_coarse_seq *
cb_coarse_add(struct cb_coarse *coarse_db,
               char *residues, int32_t start, int32_t end);

struct cb_coarse_seq *
cb_coarse_get(struct cb_coarse *coarse_db, int32_t i);

void
cb_coarse_save_binary(struct cb_coarse *coarse_db);

void
cb_coarse_save_plain(struct cb_coarse *coarse_db);

void
cb_coarse_save_seeds_binary(struct cb_coarse *coarse_db);

void
cb_coarse_save_seeds_plain(struct cb_coarse *coarse_db);


char *get_coarse_header(FILE *f);
struct cb_link_to_compressed *read_coarse_link(FILE *f);
struct DSVector *get_coarse_sequence_links(FILE *f);
struct DSVector *get_coarse_sequence_links_at(FILE *links, FILE *index,
                                                           int32_t id);
int64_t cb_coarse_find_offset(FILE *index_file, int id);
struct fasta_seq *cb_coarse_read_fasta_seq(struct cb_coarse *coarsedb,
                                            int id);

#endif
