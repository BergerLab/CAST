#include <stdio.h>
#include <stdlib.h>

#include "coarse.h"
#include "flags.h"
#include "seeds.h"

#include "ds.h"

const int8_t cb_seeds_alpha_size[] = {
    0,   /* 'A' */
    -1,  /* 'B' */
    1,   /* 'C' */
    -1,  /* 'D' */
    -1,  /* 'E' */
    -1,  /* 'F' */
    2,   /* 'G' */
    -1,  /* 'H' */
    -1,  /* 'I' */
    -1,  /* 'J' */
    -1,  /* 'K' */
    -1,  /* 'L' */
    -1,  /* 'M' */
    -1,  /* 'N' */
    -1,  /* 'O' */
    -1,  /* 'P' */
    -1,  /* 'Q' */
    -1,  /* 'R' */
    -1,  /* 'S' */
    3,   /* 'T' */
    -1,  /* 'U' */
    -1,  /* 'V' */
    -1,  /* 'W' */
    -1,  /* 'X' */
    -1,  /* 'Y' */
    -1   /* 'Z' */
};

static int32_t residue_value(char residue);

struct cb_seeds *cb_seeds_init(int32_t seed_size){
    struct cb_seeds *seeds;
    int32_t errno, p;

    seeds = malloc(sizeof(*seeds));
    assert(seeds);

    if (0 != (errno = pthread_rwlock_init(&seeds->add_lock, NULL))) {
        fprintf(stderr, "Could not create rwlock. Errno: %d\n", errno);
        exit(1);
    }

    seeds->quad_locks = malloc(256*sizeof(*(seeds->quad_locks)));
    assert(seeds->quad_locks);

    for (int i = 0; i < 256; i++) {
        pthread_rwlock_t lock;
        if (0 != (errno = pthread_rwlock_init(&lock, NULL))) {
            fprintf(stderr, "Could not create rwlock. Errno: %d\n", errno);
            exit(1);
        }
        seeds->quad_locks[i] = lock;
    }

    seeds->seed_size = seed_size;
    seeds->powers_length = seed_size + 1;
    
    seeds->powers = malloc((seeds->powers_length)*sizeof(*seeds->powers));
    assert(seeds->powers);

    p = 1;
    for (int i = 0; i < seeds->powers_length; i++) {
        seeds->powers[i] = p;
        p *= CABLAST_SEEDS_ALPHA_SIZE;
    }

    seeds->locs_length = seeds->powers[seed_size];

    seeds->locs = malloc(seeds->locs_length*sizeof(*seeds->locs));
    assert(seeds->locs);

    for (int i = 0; i < seeds->locs_length; i++)
        seeds->locs[i] =
          malloc(compress_flags.max_kmer_freq*sizeof(seeds->locs));

    seeds->loc_counts = malloc(seeds->locs_length*sizeof(*seeds->loc_counts));
    assert(seeds->loc_counts);

    int locs_length = seeds->locs_length;
    int32_t *loc_counts = seeds->loc_counts;

    for (int i = 0; i < locs_length; i++)
        loc_counts[i] = 0;

    return seeds;
}

void cb_seeds_free(struct cb_seeds *seeds){
    int32_t errno;

    if (0 != (errno = pthread_rwlock_destroy(&seeds->add_lock))) {
        fprintf(stderr, "Could not destroy rwlock. Errno: %d\n", errno);
        exit(1);
    }
    for (int i = 0; i < 256; i++)
        if (0 != (errno = pthread_rwlock_destroy(&(seeds->quad_locks[i])))) {
            fprintf(stderr, "Could not destroy rwlock. Errno: %d\n", errno);
            exit(1);
        }

    for (int i = 0; i < seeds->locs_length; i++)
        free(seeds->locs[i]);

    free(seeds->locs);
    free(seeds->loc_counts);
    free(seeds->powers);
    free(seeds);
}

struct cb_seeds_add_memory *cb_seeds_add_memory_init(){
    struct cb_seeds_add_memory *mem = malloc(sizeof(*mem));
    assert(mem);

    mem->loc_counts = malloc(256*sizeof(*(mem->loc_counts)));
    assert(mem->loc_counts);

    mem->loc_caps = malloc(256*sizeof(*(mem->loc_caps)));
    assert(mem->loc_caps);

    mem->locs = malloc(256*sizeof(*(mem->locs)));
    assert(mem->locs);

    for (int i = 0; i < 256; i++) {
        mem->locs[i] = malloc(256*sizeof(*(mem->locs[i])));
        assert(mem->locs[i]);

        mem->loc_counts[i] = 0;
        mem->loc_caps[i]   = 256;
    }

    return mem;
}

void cb_seeds_add_memory_free(struct cb_seeds_add_memory *mem) {
    for (int i = 0; i < 256; i++)
        free(mem->locs[i]);
    free(mem->loc_counts);
    free(mem->loc_caps);
    free(mem);
}

void cb_seeds_add(struct cb_seeds *seeds, struct cb_coarse_seq *seq,
                  struct cb_seeds_add_memory *seeds_mem){
    struct cb_seed_loc *loc;
    int32_t hash, q, seed_size = seeds->seed_size+1,
            seq_length = seq->seq->length, *loc_counts = seeds->loc_counts,
            id = seq->id;
    struct cb_seed_loc ***locs = seeds->locs;
    char *kmer, *residues = seq->seq->residues;

    for (int i = 0; i < 256; i++)
        seeds_mem->loc_counts[i] = 0;

    hash = hash_kmer(seeds, residues);

    if (hash != -1) {
        q = hash % 256;
        loc = cb_seed_loc_init(id, 0);

        struct cb_seed_h_loc *h_loc = malloc(sizeof(*h_loc));
        assert(h_loc);

        h_loc->loc  = loc;
        h_loc->hash = hash;

        seeds_mem->locs[q][seeds_mem->loc_counts[q]] = h_loc;
        (seeds_mem->loc_counts[q])++;
    }

    for (int i = 1; i <= seq_length - seed_size; i++) {
        kmer = residues + i;
        hash = update_kmer(seeds, kmer, hash);

        if (hash != -1) {
            q = hash % 256;
            loc = cb_seed_loc_init(id, i);

            struct cb_seed_h_loc *h_loc = malloc(sizeof(*h_loc));
            assert(h_loc);

            h_loc->loc = loc;
            h_loc->hash = hash;

            seeds_mem->locs[q][seeds_mem->loc_counts[q]] = h_loc;
            (seeds_mem->loc_counts[q])++;

            if (seeds_mem->loc_counts[q] == seeds_mem->loc_caps[q]) {
                seeds_mem->locs[q] =
                  realloc(seeds_mem->locs[q],
                          2*seeds_mem->loc_caps[q]
                           *sizeof(seeds_mem->locs[q]));
                assert(seeds_mem->locs[q]);

                seeds_mem->loc_caps[q] *= 2;
            }
        }
    }

    pthread_rwlock_wrlock(&seeds->add_lock);

    for (int i = 0; i < 256; i++)
        if (seeds_mem->loc_counts[i] > 0) {
            pthread_rwlock_t l = seeds->quad_locks[i];
            pthread_rwlock_wrlock(&l);
            for (int j = 0; j < seeds_mem->loc_counts[i]; j++) {
                struct cb_seed_h_loc *h_loc = seeds_mem->locs[i][j];
                int h = h_loc->hash;
                struct cb_seed_loc *loc = h_loc->loc;
                locs[h][loc_counts[h]] = loc;
                (loc_counts[h])++;
            }
            pthread_rwlock_unlock(&l);
        }

    pthread_rwlock_unlock(&seeds->add_lock);
}

int32_t cb_seeds_lookup(struct cb_seeds *seeds, char *kmer){
    int32_t hash = hash_kmer(seeds, kmer);
    int32_t count = 0;

    if (hash < 0)
        return -1;

    int l = hash % 256;
    pthread_rwlock_t q = seeds->quad_locks[l];

    pthread_rwlock_rdlock(&q);
    count = seeds->loc_counts[hash];
    pthread_rwlock_unlock(&q);

    return count;
}

struct cb_seed_loc *cb_seed_loc_init(uint32_t coarse_seq_id,
                                     uint16_t residue_index){
    struct cb_seed_loc *seedLoc;

    seedLoc = malloc(sizeof(*seedLoc));
    assert(seedLoc);

    seedLoc->coarse_seq_id = coarse_seq_id;
    seedLoc->residue_index = residue_index;
    seedLoc->next          = NULL;

    return seedLoc;
}

/*void cb_seed_loc_free(struct cb_seed_loc *seedLoc){
    struct cb_seed_loc *seed1, *seed2;

    for (seed1 = seedLoc; seed1 != NULL;) {
        seed2 = seed1->next;
        free(seed1);
        seed1 = seed2;
    }
}*/

static int32_t residue_value(char residue){
    int32_t i = residue - 'A', val;

    if (i < 0 || i >= 26) {
        fprintf(stderr, "Invalid nucleotide residue: %c\n", residue);
        exit(1);
    }
    val = cb_seeds_alpha_size[i];

    return val;
}

/*Takes in as input a seeds table and a k-mer and returns the k-mer's index
  in the seeds table*/
int32_t hash_kmer(struct cb_seeds *seeds, char *kmer){
    int32_t i = 0, key = 0, val = 0, seed_size = seeds->seed_size,
            *powers = seeds->powers;

    for (i = 0; i < seed_size; i++) {
        val = residue_value(kmer[i]);
        if (val == -1)
            return -1;

        key += val * powers[i];
    }
    return key;
}

/*Takes in as input a seeds table and a k-mer and returns the k-mer's index
  in the seeds table*/
int32_t update_kmer(struct cb_seeds *seeds, char *kmer, int32_t key){
    int32_t seed_size = seeds->seed_size, *powers = seeds->powers,
            val = residue_value(kmer[seed_size-1]);

    if (val == -1)
        return -1;

    key /= CABLAST_SEEDS_ALPHA_SIZE;
    key += val*powers[seed_size-1];
    
    return key;
}


/*Convert an integer to the k-mer that it represents.  Currently only works for
  size k = 10*/
char *unhash_kmer(struct cb_seeds *seeds, int hash){
    int i;
    char nucleotides[4] = {'A','C','G','T'};

    char *kmer = malloc(11*sizeof(*kmer));
    assert(kmer);

    kmer[10] = '\0';
    for (i = 0; i < seeds->seed_size; i++) {
        kmer[i] = nucleotides[hash%4];
        hash /= 4;
    }
    return kmer;
}

/*Output the seeds table in plain text format for debugging*/
/*void print_seeds(struct cb_seeds *seeds){
    int32_t i, j;

    char *kmer = malloc(seeds->seed_size*sizeof(*kmer));
    assert(kmer);

    for (i = 0; i < seeds->locs_length; i++) {
        struct cb_seed_loc *s = seeds->locs[i];
        uint32_t new_kmer = (uint32_t)0;
        char *kmer = unhash_kmer(seeds, i);

        printf("%s\n", kmer);
        for (j = 0; j < seeds->seed_size; j++) {
            new_kmer <<= 2;
            new_kmer |= ((i >> (2*j)) & ((uint32_t)3));
        }
        printf("%s\n", kmer);
        free(kmer);
        while (s) {
            printf("(%d %d) > ", s->coarse_seq_id, s->residue_index);
            s = s->next;
        }
        printf("\n");
    }
}*/
