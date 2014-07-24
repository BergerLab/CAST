#include <assert.h>
#include <libxml/parser.h>
#include <libxml/tree.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <unistd.h>

#include "clibs/include/ds.h"
#include "clibs/include/opt.h"

#include "bitpack.h"
#include "coarse.h"
#include "compressed.h"
#include "compression.h"
#include "database.h"
#include "decompression.h"
#include "fasta.h"
#include "flags.h"
#include "progress-bar.h"
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
    FILE *query_file = NULL;
    struct cb_database_r *db = NULL;
    struct opt_config *conf;
    struct opt_args *args;
    struct DSVector *queries = NULL;
    struct fasta_seq *query = NULL;
    uint64_t dbsize = 0;
    int i = 0, j = 0;

    conf = load_search_args();
    args = opt_config_parse(conf, argc, argv);

    if (args->nargs < 3) {
        fprintf(stderr, 
                "Usage: %s [flags] database-dir fasta-file results-file\n",
                argv[0]);
        opt_config_print_usage(conf);
        exit(1);
    }

    /*if (!search_flags.hide_progress)
        fprintf(stderr, "Loading database data\n\n");*/


    db = cb_database_r_init(args->args[0],
                            (search_flags.load_coarse_db ||
                             search_flags.load_coarse_residues),
                            (search_flags.load_coarse_db ||
                             search_flags.load_coarse_links),
                            search_flags.load_compressed_db,
                            search_flags.link_block_size);
    dbsize = read_int_from_file(8, db->coarse_db->db->file_params);

    /*if (!search_flags.hide_progress)
        fprintf(stderr, "Running coarse BLAST\n\n");*/
    //blat_coarse(args, dbsize);

    query_file = fopen(args->args[1], "r");
    queries = ds_vector_create();
    query = fasta_read_next(query_file, "");
    while (query) {
        ds_vector_append(queries, (void *)query);
        query = fasta_read_next(query_file, "");
    }

    fclose(query_file);

    opt_args_free(args);
    opt_config_free(conf);

    return 0;
}
