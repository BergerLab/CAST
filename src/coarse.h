#ifndef __CABLAST_COARSE_H__
#define __CABLAST_COARSE_H__

//Apparently this is required to make pthread_rwlock* stuff available.
#define __USE_UNIX98

#include <pthread.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>

#include "clibs/include/ds.h"

#include "link_to_compressed.h"
#include "seeds.h"
#include "seq.h"

struct cb_coarse_seq {
    int32_t id;
    struct cb_seq *seq;
    struct cb_link_to_compressed *links;
    struct cb_link_to_compressed *last_link;
    pthread_rwlock_t lock_links;
};

struct cb_coarse_seq *cb_coarse_seq_init(int32_t id, char *residues,
                                         int32_t start, int32_t end);

void cb_coarse_seq_free(struct cb_coarse_seq *seq);

void cb_coarse_seq_addlink(struct cb_coarse_seq *seq,
                           struct cb_link_to_compressed *newlink);

struct cb_coarse {
    struct DSVector *seqs;
    struct cb_seeds *seeds;
    uint64_t dbsize;
    pthread_rwlock_t lock_seq;

    //The FASTA file of the residues of each sequence in the coarse database
    FILE *file_fasta;

    //The links to the compressed sequences in each coarse sequence
    FILE *file_links;

    //The byte index of the start of each coarse sequence in coarse.links
    FILE *file_links_index;

    /*The indices into the bases of the coarse FASTA file for the coarse start
      and end of each link*/
    FILE *file_links_base_index;

    //The number of links to compressed sequences in each coarse sequence
    FILE *file_links_count_index;

    /*The index of the start of the header of each FASTA sequence in the coarse
      FASTA file*/
    FILE *file_fasta_index;

    //The base index of the start of each sequence in the coarse FASTA file
    FILE *file_fasta_base_index;

    //Contains the size of the original FASTA file
    FILE *file_params;
};

struct cb_coarse *
cb_coarse_init(int32_t seed_size,
               FILE *file_fasta, FILE *file_links,
               FILE *file_links_index, FILE *file_links_base_index,
               FILE *file_links_count_index, FILE *file_fasta_index,
               FILE *file_fasta_base_index, FILE *file_params, bool read);

void cb_coarse_free(struct cb_coarse *coarse_db);

struct cb_coarse_seq *cb_coarse_add(struct cb_coarse *coarse_db, char *residues,
                                    int32_t start, int32_t end,
                                    struct cb_seeds_add_memory *seeds_mem);

struct cb_coarse_seq *cb_coarse_get(struct cb_coarse *coarse_db, int32_t i);

void cb_coarse_db_update_dbsize(struct cb_coarse *coarse_db, int32_t size);

void cb_coarse_save_binary(struct cb_coarse *coarse_db);
void cb_coarse_save_plain(struct cb_coarse *coarse_db);

struct cb_link_to_compressed_data *read_coarse_link_data(FILE *f);
struct DSVector *read_coarse_links(FILE *f, int64_t num_links);
int64_t cb_coarse_find_offset(FILE *index_file, int id);
struct fasta_seq *cb_coarse_read_fasta_seq(struct cb_coarse *coarsedb,
                                           int id);

struct cb_coarse_r {
    struct cb_coarse *db;
    struct DSVector *links;
    struct DSVector *link_inds_by_block;
    struct DSVector *link_inds_by_seq;
    int32_t link_block_size;
    int64_t num_coarse_seqs;
    int64_t *seq_link_counts;
    int64_t *seq_base_indices;
    char *all_residues;
};

struct cb_coarse_r *
cb_coarse_r_init(int32_t seed_size,
                 FILE *file_fasta, FILE *file_links,
                 FILE *file_links_index, FILE *file_links_base_index,
                 FILE *file_links_count_index, FILE *file_fasta_index,
                 FILE *file_fasta_base_index, FILE *file_params,
                 bool load_coarse_residues, bool load_coarse_links,
                 int32_t link_block_size);

void cb_coarse_r_free(struct cb_coarse_r *coarse_db);
void cb_coarse_r_init_blocks(struct cb_coarse_r *coarse_db);
struct DSVector *cb_coarse_r_get_block(struct cb_coarse_r *coarse_db,
                                       int32_t index);

void cb_coarse_r_read_all_residues(struct cb_coarse_r *coarse_db);
void cb_coarse_r_read_all_links(struct cb_coarse_r *coarse_db);

char *cb_coarse_r_get_seq_residues(struct cb_coarse_r *coarse_db,
                                   int64_t id);

/*Coarse database functions ending in _r are used on cb_coarse_r structs
 *and are used as wrapper functions for the regular coarse database functions
 *being called on the cb_coarse_r struct's coarsedb.
 */
struct cb_coarse_seq *cb_coarse_get_r(struct cb_coarse_r *coarse_db,
                                      int32_t i);

struct fasta_seq *cb_coarse_read_fasta_seq_r(struct cb_coarse_r *coarsedb,
                                             int64_t id);

#endif
