#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fasta.h"
#include "util.h"

static bool
is_new_sequence_start(FILE *f);

static void
exclude_residues(char *seq, const char *exclude);

static void *
fasta_generator(void *gen);

struct fasta_file *
fasta_read_all(const char *file_name, const char *exclude)
{
    FILE *f;
    struct fasta_file *ff;
    struct fasta_seq *new_seq;
    int allocated = 0;

    printf("Reading all sequences from: %s\n", file_name);

    f = fopen(file_name, "r");
    if (f == NULL) {
        perror("fasta_read_all");
        exit(1);
    }

    ff = malloc(sizeof(*ff));
    assert(ff);

    ff->length = 0;

    allocated = FASTA_INITIAL_SEQUENCE_LENGTH;

    ff->seqs = malloc(allocated * sizeof(*ff->seqs));
    assert(ff->seqs);

    while (NULL != (new_seq = fasta_read_next(f, exclude))) {
        if (ff->length == allocated) {
            allocated *= 1.5;

            ff->seqs = realloc(ff->seqs, allocated * sizeof(*ff->seqs));
            assert(ff->seqs);
        }
        ff->seqs[ff->length++] = new_seq;
    }

    fclose(f);

    return ff;
}

void
fasta_free_all(struct fasta_file *ff)
{
    int i;

    for (i = 0; i < ff->length; i++)
        fasta_free_seq(ff->seqs[i]);
    free(ff->seqs);
    free(ff);
}

struct fasta_seq *
fasta_read_next(FILE *f, const char *exclude)
{
    struct fasta_seq *fs;
    char *line = NULL;

    /* check to make sure the next line starts a new sequence record */
    if (!is_new_sequence_start(f)){
        return NULL;}

    /* read in the sequence id */
    if (0 == readline(f, &line)) {
        free(line);
        return NULL;
    }

    line = trim(line, "> \n\r\t");

    fs = malloc(sizeof(*fs));
    assert(fs);

    fs->name = malloc((strlen(line) + 1) * sizeof(*fs->name));
    assert(fs->name);

    strcpy(fs->name, line);
    free(line);

    /* Now read all of the sequence data for this record */
    fs->seq = malloc(sizeof(*fs->seq));
    assert(fs->seq);

    fs->seq[0] = '\0';
    while (!is_new_sequence_start(f)) {
        if (0 == readline(f, &line)) {
            free(line);
            return fs;
        }

        line = trim(line, "* \n\r\t");
        exclude_residues(line, exclude);
        exclude_residues(line, "*");

        fs->seq = realloc(
            fs->seq, sizeof(*fs->seq) * (1 + strlen(line) + strlen(fs->seq)));
        assert(fs->seq);

        strcat(fs->seq, line);
        free(line);
    }

    return fs;
}

void
fasta_free_seq(struct fasta_seq *seq)
{
    free(seq->name);
    free(seq->seq);
    free(seq);
}

struct fasta_seq_gen *
fasta_generator_start(const char *file_name, const char *exclude,
                      int buffer_capacity)
{
    FILE *fp;
    struct fasta_seq_gen *fsg;
    int errno;

    assert(buffer_capacity > 0);

    if (NULL == (fp = fopen(file_name, "r"))) {
        perror("fasta_start_generator");
        exit(1);
    }

    fsg = malloc(sizeof(*fsg));
    assert(fsg);

    fsg->seqs = ds_queue_create(buffer_capacity);
    fsg->fp = fp;
    fsg->exclude = exclude;

    errno = pthread_create(&fsg->thread, NULL, fasta_generator, (void*) fsg);
    if (0 != errno) {
        fprintf(stderr, "Could not create pthread. Errno: %d\n", errno);
        exit(1);
    }

    return fsg;
}

void
fasta_generator_free(struct fasta_seq_gen *fsg)
{
    int errno;

    if (0 != (errno = pthread_join(fsg->thread, NULL))) {
        fprintf(stderr, "Could not join thread. Errno: %d\n", errno);
        exit(1);
    }
    fclose(fsg->fp);
    ds_queue_free(fsg->seqs);
    free(fsg);
}

static void *
fasta_generator(void *gen)
{
    struct fasta_seq_gen *fsg;
    struct fasta_seq *seq;

    fsg = (struct fasta_seq_gen *) gen;

    while (NULL != (seq = fasta_read_next(fsg->fp, fsg->exclude))) {
        ds_queue_put(fsg->seqs, seq);
    }

    ds_queue_close(fsg->seqs);
    return NULL;
}

struct fasta_seq *
fasta_generator_next(struct fasta_seq_gen *fsg)
{
    struct fasta_seq *seq;

    seq = (struct fasta_seq *) ds_queue_get(fsg->seqs);
    return seq;
}

static bool
is_new_sequence_start(FILE *f)
{
    char next;
    bool is_new_seq;

    next = fgetc(f);
    is_new_seq = next == '>';
    if (EOF == ungetc(next, f)) {
        return false;
        perror("is_new_sequence_start");
        exit(1);
    }

    return is_new_seq;
}

static void
exclude_residues(char *seq, const char *exclude)
{
    int i, j;

    for (i = 0; seq[i] != '\0'; i++)
        for (j = 0; exclude[j] != '\0'; j++)
            if (seq[i] == exclude[j])
                seq[i] = 'X';
}
