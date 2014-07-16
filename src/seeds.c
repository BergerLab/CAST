#include <stdio.h>
#include <stdlib.h>

#include "coarse.h"
#include "flags.h"
#include "seeds.h"

/*The residue values for each letter of the alphabet, used for hashing k-mers
  -1 = invalid residue.*/
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

/*Takes in the length of the k-mers used and creates a new seeds table with one
  array of seed locations for each k-mer.*/
struct cb_seeds *cb_seeds_init(int32_t seed_size){
    struct cb_seeds *seeds;
    int32_t errno, p;

    seeds = malloc(sizeof(*seeds));
    assert(seeds);

    seeds->alloc_locks = malloc(256*sizeof(*(seeds->alloc_locks)));
    assert(seeds->alloc_locks);

    /*Create one alloc lock for k-mers ending in each of the 256 possible sets
     *of 4 bases.  These are used for allocating an array of seed locations for
     *a k-mer.
     */
    for (int i = 0; i < 256; i++)
        if (0 != (errno = pthread_rwlock_init(&(seeds->alloc_locks[i]),NULL))) {
            fprintf(stderr, "Could not create rwlock. Errno: %d\n", errno);
            exit(1);
        }

    //Create the lock for adding seeds to the seeds table.
    if (0 != (errno = pthread_rwlock_init(&seeds->lock, NULL))) {
        fprintf(stderr, "Could not create rwlock. Errno: %d\n", errno);
        exit(1);
    }

    seeds->seed_size     = seed_size;
    seeds->powers_length = seed_size + 1;
    
    seeds->powers = malloc((seeds->powers_length)*sizeof(*seeds->powers));
    assert(seeds->powers);

    p = 1;
    for (int i = 0; i < seeds->powers_length; i++) {
        seeds->powers[i] = p;
        p *= CABLAST_SEEDS_ALPHA_SIZE;
    }

    seeds->locs_length = seeds->powers[seed_size];

    seeds->locs        = malloc(seeds->locs_length*sizeof(*seeds->locs));
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

/*Frees the seeds table and all of the seed locations in the table.*/
void cb_seeds_free(struct cb_seeds *seeds){
    int32_t errno;

    //Destroy the locks and free the array of alloc locks.
    if (0 != (errno = pthread_rwlock_destroy(&seeds->lock))) {
        fprintf(stderr, "Could not destroy rwlock. Errno: %d\n", errno);
        exit(1);
    }

    for (int i = 0; i < 256; i++)
        if (0 != (errno = pthread_rwlock_destroy(&(seeds->alloc_locks[i])))) {
            fprintf(stderr, "Could not create rwlock. Errno: %d\n", errno);
            exit(1);
        }
    free(seeds->alloc_locks);

    for (int i = 0; i < seeds->locs_length; i++)
        free(seeds->locs[i]);

    free(seeds->locs);
    free(seeds->loc_counts);
    free(seeds->powers);
    free(seeds);
}

/*Takes in the seeds table and a coarse sequence and adds a new seed location
 *in the seeds table for each coarse sequence.  If a k-mer has more seed
 *locations than the number specified with the --max-kmer-freq flag
 *(default 500), then no new seed locations will be added for that k-mer.
 */
void cb_seeds_add(struct cb_seeds *seeds, struct cb_coarse_seq *seq,
                  struct cb_seeds_add_memory *mem){
    struct cb_seed_loc **locs_to_add = mem->locs;
    int32_t hash, seed_size = seeds->seed_size+1,
            seq_length = seq->seq->length, *loc_counts = seeds->loc_counts,
            id = seq->id, *hashes = mem->hashes;

    struct cb_seed_loc ***locs = seeds->locs;
    char *kmer, *residues = seq->seq->residues;

    hash = hash_kmer(seeds, residues);
    hashes[0] = hash;
    locs_to_add[0] = hash != -1 ? cb_seed_loc_init(id, 0) : NULL;

    //Fill the array of hashes and create the seed locations to add 
    for (int i = 1; i <= seq_length - seed_size; i++) {
        kmer = residues + i;
        hash = update_kmer(seeds, kmer, hash);
        hashes[i] = hash;
        locs_to_add[i] = hash != -1 ? cb_seed_loc_init(id, i) : NULL;
    }

    /*Allocate a seed locations array for any k-mers being added that currently
      don't have a locations array in the seeds table*/
    for (int i = 0; i <= seq_length - seed_size; i++) {
        hash = hashes[i];
        if (locs[hash] == NULL) {
            pthread_rwlock_wrlock(&(seeds->alloc_locks[hash%256]));
            if (locs[hash] == NULL)
                locs[hash]=malloc(compress_flags.max_kmer_freq*sizeof(locs[hash]));
            assert(locs[hash]);
            pthread_rwlock_unlock(&(seeds->alloc_locks[hash%256]));
        }
    }

    /*For each k-mer, add its seed locations to the seeds table and update that
      k-mer's entry in loc_counts.*/
    pthread_rwlock_wrlock(&seeds->lock);
    for (int i = 0; i <= seq_length - seed_size; i++) {
        hash = hashes[i];
        if (locs_to_add[i] != NULL &&
              loc_counts[hash] < compress_flags.max_kmer_freq) {
            locs[hash][loc_counts[hash]] = locs_to_add[i];
            (loc_counts[hash])++;
        }
    }
    pthread_rwlock_unlock(&seeds->lock);
}

/*Takes in the seeds table and a k-mer and returns the number of seed locations
  in the table for that k-mer.*/
int32_t cb_seeds_lookup(struct cb_seeds *seeds, char *kmer){
    int32_t hash = hash_kmer(seeds, kmer);
    return hash >= 0 ? seeds->loc_counts[hash] : -1;
}

/*Takes in a coarse sequence ID number and an index into the coarse sequence
  with that ID number and uses them to create a new seed location struct.*/
struct cb_seed_loc *cb_seed_loc_init(uint32_t coarse_seq_id,
                                     uint16_t residue_index){
    struct cb_seed_loc *seedLoc = malloc(sizeof(*seedLoc));
    assert(seedLoc);

    seedLoc->coarse_seq_id = coarse_seq_id;
    seedLoc->residue_index = residue_index;

    return seedLoc;
}

/*Takes in a residue and returns its value in the cb_seeds_alpha_size array,
  which is used for hashing k-mers.*/
static inline int32_t residue_value(char residue){
    if (residue < 'A' || residue > 'Z') {
        fprintf(stderr, "Invalid nucleotide residue: %c\n", residue);
        exit(1);
    }

    return cb_seeds_alpha_size[residue-'A'];
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

/*Allocates a new seeds_add_memory for a compression worker, which contains
  an array of k-mers as integers and an array of seed locations.*/
struct cb_seeds_add_memory *cb_seeds_add_memory_init(){
    struct cb_seeds_add_memory *mem = malloc(sizeof(*mem));
    assert(mem);

    mem->hashes =
      malloc(2*compress_flags.max_chunk_size*sizeof(*(mem->hashes)));
    assert(mem->hashes);

    mem->locs   = malloc(2*compress_flags.max_chunk_size*sizeof(*(mem->locs)));
    assert(mem->locs);

    return mem;
}

