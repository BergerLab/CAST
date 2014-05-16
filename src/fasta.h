#ifndef __CABLAST_FASTA_H__
#define __CABLAST_FASTA_H__

#include <pthread.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>

#include "ds.h"

#define FASTA_INITIAL_SEQUENCE_LENGTH 1000
#define FASTA_MAX_LINE 1024
#define FASTA_EXCLUDE_NCBI_BLOSUM62 "JOU"

struct fasta_file {
    struct fasta_seq **seqs;
    int32_t length;
};

struct fasta_seq {
    char *name;
    char *seq;
};

struct fasta_seq_gen {
    pthread_t thread;
    FILE *fp;
    struct DSQueue *seqs;
    const char *exclude;
};

struct fasta_file *
fasta_read_all(const char *file_name, const char *exclude);

struct fasta_seq *
fasta_read_next(FILE *f, const char *exclude);

void
fasta_free_all(struct fasta_file *ff);

void
fasta_free_seq(struct fasta_seq *seq);

struct fasta_seq_gen *
fasta_generator_start(const char *file_name, const char *exclude,
                      int buffer_capacity);

void
fasta_generator_free(struct fasta_seq_gen *fsg);

struct fasta_seq *
fasta_generator_next(struct fasta_seq_gen *fsg);

#endif
