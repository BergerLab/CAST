#include <stdio.h>
#include <stdlib.h>

#include "coarse.h"
#include "flags.h"
#include "seeds.h"

const int8_t cb_seeds_alpha_size[] = {
    0,  /* 'A' */
    -1, /* 'B' */
    1,  /* 'C' */
    -1, /* 'D' */
    -1, /* 'E' */
    -1, /* 'F' */
    2,  /* 'G' */
    -1, /* 'H' */
    -1, /* 'I' */
    -1, /* 'J' */
    -1, /* 'K' */
    -1, /* 'L' */
    -1, /* 'M' */
    -1, /* 'N' */
    -1, /* 'O' */
    -1, /* 'P' */
    -1, /* 'Q' */
    -1, /* 'R' */
    -1, /* 'S' */
    3,  /* 'T' */
    -1, /* 'U' */
    -1, /* 'V' */
    -1, /* 'W' */
    -1, /* 'X' */
    -1, /* 'Y' */
    -1  /* 'Z' */
};

static int32_t residue_value(char residue);

struct cb_seeds *cb_seeds_init(int32_t seed_size){
    struct cb_seeds *seeds;
    int32_t errno, p;

    seeds = malloc(sizeof(*seeds));
    assert(seeds);

    if (0 != (errno = pthread_rwlock_init(&seeds->lock, NULL))) {
        fprintf(stderr, "Could not create rwlock. Errno: %d\n", errno);
        exit(1);
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
        seeds->locs[i] = NULL;

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

    if (0 != (errno = pthread_rwlock_destroy(&seeds->lock))) {
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

void cb_seeds_add(struct cb_seeds *seeds, struct cb_coarse_seq *seq){
    struct cb_seed_loc *loc;
    int32_t hash, seed_size = seeds->seed_size+1,
            seq_length = seq->seq->length, *loc_counts = seeds->loc_counts,
            id = seq->id;
    struct cb_seed_loc ***locs = seeds->locs;
    char *kmer, *residues = seq->seq->residues;

    pthread_rwlock_wrlock(&seeds->lock);

    hash = hash_kmer(seeds, residues);
    if (hash != -1 && loc_counts[hash] < compress_flags.max_kmer_freq) {
        loc = cb_seed_loc_init(id, 0);

        if (locs[hash] == NULL) {
            locs[hash]=malloc(compress_flags.max_kmer_freq*sizeof(locs[hash]));
            assert(locs[hash]);
        }

        locs[hash][loc_counts[hash]] = loc;
        (loc_counts[hash])++;
    }
   
    for (int i = 1; i <= seq_length - seed_size; i++) {
        kmer = residues + i;

        hash = update_kmer(seeds, kmer, hash);
        if (hash == -1 ||
            loc_counts[hash] >= compress_flags.max_kmer_freq)
            continue;

        loc = cb_seed_loc_init(id, i);

        if (locs[hash] == NULL) {
            locs[hash]=malloc(compress_flags.max_kmer_freq*sizeof(locs[hash]));
            assert(locs[hash]);
        }

        locs[hash][loc_counts[hash]] = loc;
        (loc_counts[hash])++;
    }

    pthread_rwlock_unlock(&seeds->lock);
}

int32_t cb_seeds_lookup(struct cb_seeds *seeds, char *kmer){
    int32_t hash = hash_kmer(seeds, kmer);
    int32_t count = 0;

    if (hash < 0)
        return -1;

    pthread_rwlock_rdlock(&seeds->lock);
    count = seeds->loc_counts[hash];
    pthread_rwlock_unlock(&seeds->lock);

    return count;
}

struct cb_seed_loc *cb_seed_loc_init(uint32_t coarse_seq_id,
                                     uint16_t residue_index){
    struct cb_seed_loc *seedLoc;

    seedLoc = malloc(sizeof(*seedLoc));
    assert(seedLoc);

    seedLoc->coarse_seq_id = coarse_seq_id;
    seedLoc->residue_index = residue_index;
    seedLoc->next = NULL;

    return seedLoc;
}

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


/*Convert an integer to the k-mer that it represents. Currently only works for
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
