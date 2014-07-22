#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fasta.h"
#include "util.h"

static bool is_new_sequence_start(FILE *f);

static void exclude_residues(char *seq, const char *exclude);

static void *fasta_generator(void *gen);

/*Takes in a FASTA file pointer and an array of residue characters to exclude
 *and reads in the next FASTA sequence to create a FASTA sequence struct
 *that excludes the characters passed into "exclude", returning it as a
 *fasta_seq struct.  If the file pointer is not at the start of a FASTA sequence
 *or the sequence id is not read in, a NULL pointer is returned. 
 */
struct fasta_seq *fasta_read_next(FILE *f, const char *exclude){
    struct fasta_seq *fs;
    char *line = NULL, *seq = NULL;

    //Check to make sure the next line starts a new sequence record.
    if (!is_new_sequence_start(f))
        return NULL;

    //Read in the sequence id.
    if (0 == readline(f, &line)) {
        free(line);
        return NULL;
    }

    line = trim(line, "> \n\r\t");

    fs = malloc(sizeof(*fs));
    assert(fs);

    fs->name = malloc((strlen(line)+1)*sizeof(*fs->name));
    assert(fs->name);

    strcpy(fs->name, line);
    free(line);

    //Now read all of the sequence data for this record.
    seq = malloc(sizeof(*seq));
    assert(seq);

    seq[0] = '\0';
    int seq_length = 0;
    while (!is_new_sequence_start(f)) {
        if (0 == readline(f, &line)) {
            free(line);
            seq[seq_length] = '\0';
            fs->seq = seq;
            return fs;
        }

        line = trim(line, "* \n\r\t");
        exclude_residues(line, exclude);
        exclude_residues(line, "*");
        int line_length = strlen(line);

        seq = realloc(seq, sizeof(*seq)*(1+line_length+seq_length));
        assert(seq);

        for (int i = 0; line_length - i != 0; i++)
            seq[seq_length+i] = line[i];
        seq_length += line_length;

        free(line);
    }
    seq[seq_length] = '\0';

    fs->seq = seq;
    return fs;
}

//Frees a FASTA sequence struct.
void fasta_free_seq(struct fasta_seq *seq){
    free(seq->name);
    free(seq->seq);
    free(seq);
}

/*Takes in the name of a FASTA file, an array of residue characters to exclude,
 *and a buffer capacity and creates a new FASTA generator struct with its queue
 *having the capacity passed into buffer_capacity, also starting a pthread for
 *running fasta_generator on the file.
 */
struct fasta_seq_gen *
fasta_generator_start(const char *file_name, const char *exclude,
                      int buffer_capacity){
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

    fsg->seqs    = ds_queue_create(buffer_capacity);
    fsg->fp      = fp;
    fsg->exclude = exclude;

    errno = pthread_create(&fsg->thread, NULL, fasta_generator, (void *)fsg);
    if (0 != errno) {
        fprintf(stderr, "Could not create pthread. Errno: %d\n", errno);
        exit(1);
    }

    return fsg;
}

//Frees the FASTA generator, joins its thread, and frees its thread-safe queue.
void fasta_generator_free(struct fasta_seq_gen *fsg){
    int errno;

    if (0 != (errno = pthread_join(fsg->thread, NULL))) {
        fprintf(stderr, "Could not join thread. Errno: %d\n", errno);
        exit(1);
    }
    fclose(fsg->fp);
    ds_queue_free(fsg->seqs);
    free(fsg);
}

/*The FASTA generator function, which is run in a pthread and reads FASTA
 *sequences from the FASTA file and puts them into the generator's thread-safe
 *queue.
 */
static void *fasta_generator(void *gen){
    struct fasta_seq_gen *fsg;
    struct fasta_seq *seq;

    fsg = (struct fasta_seq_gen *)gen;

    while (NULL != (seq = fasta_read_next(fsg->fp, fsg->exclude))) {
        ds_queue_put(fsg->seqs, seq);
    }

    ds_queue_close(fsg->seqs);
    return NULL;
}

/*Takes in a FASTA generator and returns the next FASTA sequence in the
  generator's thread-safe queue.*/
struct fasta_seq *fasta_generator_next(struct fasta_seq_gen *fsg){
    struct fasta_seq *seq;

    seq = (struct fasta_seq *)ds_queue_get(fsg->seqs);
    return seq;
}

/*Takes in a file pointer and checks if the pointer is at the start of a new
  FASTA sequence.*/
static bool is_new_sequence_start(FILE *f){
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

/*Takes in a sequence and an array of residues to exclude and replaces all
  excluded residues in the sequence with X's.*/
static void exclude_residues(char *seq, const char *exclude){
    for (int i = 0; seq[i] != '\0'; i++)
        for (int j = 0; exclude[j] != '\0'; j++)
            if (seq[i] == exclude[j])
                seq[i] = 'X';
}
