#include <assert.h>
#include <errno.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <unistd.h>

#include "clibs/include/opt.h"

#include "coarse.h"
#include "compressed.h"
#include "database.h"
#include "edit_scripts.h"
#include "flags.h"
#include "fasta.h"
#include "seq.h"

//Joins two strings to create a new file path string.
static char *path_join(char *a, char *b){
    char *joined;

    joined = malloc((1+strlen(a)+1+strlen(b))*sizeof(*joined));
    assert(joined);

    sprintf(joined, "%s/%s", a, b);
    return joined;
}

int main(int argc, char **argv){ 
    struct cb_database_r *db;
    struct fasta_seq_gen *fsg;
    struct fasta_seq *seq;
    struct timeval start;
    struct fasta_seq **coarse_sequences;
    struct cb_link_to_coarse *link;
    struct opt_config *conf = load_compress_args();
    struct cb_compressed_seq **compressed;
    uint64_t last_end, num_coarse_sequences = 0;
    int i, org_seq_id, overlap;
    char *fasta_filename, *compressed_filename;

    FILE *compressed_file;

    struct opt_args *args = opt_config_parse(conf, argc, argv);
    if (args->nargs < 2) {
        fprintf(stderr, "Usage: %s [flags] database-dir output-fasta-file\n",
                argv[0]);
        exit(1);
    }

    db = cb_database_r_init(args->args[0], false, false, false, 30000);

    org_seq_id = 0;
    gettimeofday(&start, NULL);
    fasta_filename      = path_join(args->args[0], CABLAST_COARSE_FASTA);
    compressed_filename = path_join(args->args[0], CABLAST_COMPRESSED);
    if (NULL == (compressed_file = fopen(compressed_filename, "r"))) {
        fprintf(stderr, "fopen: 'fopen %s' failed: %s\n",
                compressed_filename, strerror(errno));
        exit(1);
    }
    compressed = read_compressed(compressed_file);

    coarse_sequences = malloc(10000*sizeof(*coarse_sequences));
    assert(coarse_sequences);

    fsg = fasta_generator_start(fasta_filename, "", 100);
    while (NULL != (seq = fasta_generator_next(fsg)))
        coarse_sequences[num_coarse_sequences++] = seq;
    for (i = 0; compressed[i] != NULL; i++) {
        int current_chunk = 0;
        last_end = 0;
        printf(">%s\n", compressed[i]->name);
        for (link = (compressed[i])->links; link != NULL; link = link->next) {
            struct cb_seq *chunk =
              cb_seq_init_range(-1, "",
                                coarse_sequences[link->coarse_seq_id]->seq,
                                link->coarse_start, link->coarse_end + 1);
            int length;
            char *decompressed;

            /*overlap represents the length of the overlap of the parts of the
              decompressed sequence that has been printed and the parts of the
              decompressed sequence currently being decompressed.*/
            overlap = last_end - link->original_start;

            for (length = 0; chunk->residues[length] != '\0'; length++);

            decompressed = read_edit_script(link->diff,chunk->residues,length);

            /*Print all characters of the decompressed chunk past the index
             *"overlap" unless overlap is greater than the length of the
             *decompressed chunk.
             */
            decompressed += overlap;
            if (overlap < link->original_end - link->original_start ||
                 (overlap == link->original_end - link->original_start &&
                  !link->next))
                printf("%s", decompressed);
            decompressed -= overlap;

            free(decompressed);

            if (link->original_end > last_end)
                last_end = link->original_end + 1;

            cb_seq_free(chunk);

            current_chunk++;
        }
        putc('\n', stdout);
    }

    free(fasta_filename);
    free(compressed_filename);

    for (i = 0; compressed[i] != NULL; i++)
        cb_compressed_seq_free(compressed[i]);
    free(compressed);

    fasta_generator_free(fsg);
    for (i = 0; i < num_coarse_sequences; i++)
        free(coarse_sequences[i]);
    free(coarse_sequences);

    cb_database_r_free(db);
    opt_config_free(conf);
    opt_args_free(args);

    return 0;
}
