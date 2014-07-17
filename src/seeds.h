#ifndef __CABLAST_SEEDS_H__
#define __CABLAST_SEEDS_H__

//Apparently this is required to make pthread_rwlock* stuff available.
#define __USE_UNIX98

#include <pthread.h>
#include <stdint.h>

#define CABLAST_SEEDS_ALPHA_SIZE 4

const int8_t cb_seeds_alpha_size[26];

struct cb_seed_loc {
    uint32_t coarse_seq_id;
    uint16_t residue_index;
    struct cb_seed_loc *next;
};

struct cb_seed_loc *cb_seed_loc_init(uint32_t coarse_seq_id,
                                     uint16_t residue_index);

struct cb_seeds {
    int32_t seed_size;
    int32_t *loc_counts;
    struct cb_seed_loc ***locs;
    int32_t locs_length;
    int32_t *powers;
    int32_t powers_length;
    pthread_rwlock_t lock;
    pthread_rwlock_t *alloc_locks;
};

struct cb_seeds_add_memory {
    int32_t *hashes;
    struct cb_seed_loc **locs;
};

struct cb_seeds_add_memory *cb_seeds_add_memory_init();

struct cb_seeds *cb_seeds_init(int32_t seed_size);

void cb_seeds_free(struct cb_seeds *seeds);

struct cb_coarse_seq;

void cb_seeds_add(struct cb_seeds *seeds, struct cb_coarse_seq *seq,
                  struct cb_seeds_add_memory *mem);

int32_t cb_seeds_lookup(struct cb_seeds *seeds, char *kmer);

int32_t hash_kmer(struct cb_seeds *seeds, char *kmer);
int32_t update_kmer(struct cb_seeds *seeds, char *kmer, int32_t key);
char *unhash_kmer(struct cb_seeds *seeds, int hash);

#endif
