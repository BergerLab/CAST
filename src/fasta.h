#ifndef __CABLAST_FASTA_H__
#define __CABLAST_FASTA_H__

#include <pthread.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>

#include "clibs/include/ds.h"

struct fasta_seq {
    char *name;
    char *seq;
};

struct fasta_seq *fasta_read_next(FILE *f, const char *exclude);

void fasta_free_seq(struct fasta_seq *seq);


struct fasta_seq_gen {
    pthread_t thread;
    FILE *fp;
    struct DSQueue *seqs;
    const char *exclude;
};

struct fasta_seq_gen *
fasta_generator_start(const char *file_name, const char *exclude,
                      int buffer_capacity);

void fasta_generator_free(struct fasta_seq_gen *fsg);

struct fasta_seq *fasta_generator_next(struct fasta_seq_gen *fsg);

#endif
