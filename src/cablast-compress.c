#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <unistd.h>

#include "clibs/include/opt.h"

#include "coarse.h"
#include "compression.h"
#include "database.h"
#include "fasta.h"
#include "flags.h"
#include "seq.h"
#include "util.h"

//Joins two strings to create a new file path string.
static char *path_join(char *a, char *b){
    char *joined;

    joined = malloc((1+strlen(a)+1+strlen(b))*sizeof(*joined));
    assert(joined);

    sprintf(joined, "%s/%s", a, b);
    return joined;
}

int main(int argc, char **argv){
    struct cb_database *db;
    struct cb_compress_workers *workers;
    struct fasta_seq_gen *fsg;
    struct fasta_seq *seq;
    struct cb_seq *org_seq;
    int org_seq_id, len_command;
    struct timeval start, current;
    long double elapsed;
    char *coarse_filename, *makeblastdb;
    struct opt_config *conf = load_compress_args();
    struct opt_args *args = opt_config_parse(conf, argc, argv);

    if (args->nargs < 2) {
        fprintf(stderr, 
          "Usage: %s [flags] database-dir fasta-file [ fasta-file ... ]\n",
          argv[0]);
        opt_config_print_usage(conf);
        exit(1);
    }

    assert(compress_flags.map_seed_size > 0 &&
           compress_flags.map_seed_size < 16);
    db = cb_database_init(args->args[0], compress_flags.map_seed_size);
    workers = cb_compress_start_workers(db, compress_flags.procs);

    org_seq_id = 0;
    gettimeofday(&start, NULL);
    for (int i = 1; i < args->nargs; i++) {
        fsg = fasta_generator_start(args->args[i], "", 100);

        while (NULL != (seq = fasta_generator_next(fsg))) {
            org_seq = cb_seq_init(org_seq_id, seq->name, seq->seq);
            cb_compress_send_job(workers, org_seq);

            fasta_free_seq(seq);

            org_seq_id++;
            if (compress_flags.notify_every > 0 &&
                org_seq_id % compress_flags.notify_every == 0) {
                gettimeofday(&current, NULL);
                elapsed = (long double)(current.tv_sec - start.tv_sec);
                fprintf(stderr, "%d sequences compressed (%0.4Lf seqs/sec)\n",
                                org_seq_id, ((long double)org_seq_id)/elapsed);
            }
        }

        fasta_generator_free(fsg);
    }

    cb_compress_join_workers(workers);
    cb_coarse_save_binary(db->coarse_db);

    coarse_filename = path_join(args->args[0], "coarse.fasta");
    len_command = strlen("makeblastdb -dbtype nucl -in  -out") + 2
                  * strlen(coarse_filename) + 1;

    makeblastdb = malloc(len_command*sizeof(makeblastdb));
    assert(makeblastdb);

    sprintf(makeblastdb, "makeblastdb -dbtype nucl -in %s -out %s",
                         coarse_filename, coarse_filename);

    system(makeblastdb);

    free(coarse_filename);
    free(makeblastdb);
    cb_database_free(db);
    opt_config_free(conf);
    opt_args_free(args);

    cb_compress_free_workers(workers);

    return 0;
}
