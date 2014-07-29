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
#include "psl.h"
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

/*Takes in the arguments for the program and returns a string of all of the
 *arguments for fine BLAT (the arguments after --blat-args) concatenated and
 *separated by a space.
 */
char *get_blat_args(struct opt_args *args){
    int i = 0, index = -1, length = 1;
    char *blat_args = NULL;

    for (i = 0; i < args->nargs; i++)
        if (index >= 0)
            length += strlen(args->args[i]) + 1;
        else
            index = strcmp(args->args[i], "--blat_args") == 0 ? i : -1;

    blat_args = malloc(length*sizeof(*args));
    assert(blat_args);

    *blat_args = '\0';
    if (index != -1)
        for (i = index + 1; i < args->nargs; i++) {
            blat_args = strcat(blat_args, args->args[i]);
            if (i < args->nargs - 1)
                blat_args = strcat(blat_args, " ");
        }
    return blat_args;
}

/*Takes in the vector of hits from coarse BLAT and the database we are using
 *for CaBLAT and returns a vector of every original sequence section re-created
 *from the calls to cb_coarse_expand for the hits we are expanding.
 */
struct DSVector *expand_blat_hits(struct DSVector *hits,
                                  struct cb_database_r *db){
    struct DSVector *expanded_hits = ds_vector_create();

    for (int i = 0; i < hits->size; i++) {
        struct psl_entry *h = (struct psl_entry *)ds_vector_get(hits, i);

        int32_t coarse_start  = h->t_start, coarse_end = h->t_end,
                coarse_seq_id = atoi(h->t_name);

        struct DSVector *oseqs =
          cb_coarse_expand(db->coarse_db, db->com_db, coarse_seq_id,
                           coarse_start, coarse_end, 10);

        for (int j = 0; j < oseqs->size; j++)
            ds_vector_append(expanded_hits, ds_vector_get(oseqs, j));

        ds_vector_free_no_data(oseqs);
    }

    return expanded_hits;
}


/*Runs BLAT on the coarse FASTA file and stores the results in a temporary
  psl file.*/
void blat_coarse(struct opt_args *args){
    char *target_path = path_join(args->args[0], CABLAST_COARSE_FASTA),
         *coarse_blat_command =
            malloc(
              (strlen("$HOME/bin/$MACHTYPE/blat    -noHead -minIdentity=80")
               +strlen(target_path)+strlen(args->args[1])
               +strlen("coarse-blat.psl"))*sizeof(*coarse_blat_command));

    sprintf(coarse_blat_command,
            "$HOME/bin/$MACHTYPE/blat %s %s %s -noHead -minIdentity=80",
            target_path, args->args[1], "coarse-blat.psl");

    if (!cablat_flags.hide_progress)
        fprintf(stderr, "%s\n", coarse_blat_command);

    system(coarse_blat_command);
}

//Runs BLAT on the fine FASTA file
void blat_fine(struct opt_args *args, uint64_t dbsize){
    char *blat, *blat_args = get_blat_args(args);
    int command_length = 1024;

    blat = malloc(command_length*sizeof(*blat));
    assert(blat);

    if (blat_args[0] == '\0')
        sprintf(blat, "$HOME/bin/$MACHTYPE/blat CaBLAT_fine.fasta %s %s",
                args->args[1], args->args[2]);
    else
        sprintf(blat, "$HOME/bin/$MACHTYPE/blat %s CaBLAT_fine.fasta %s %s",
                blat_args, args->args[1], args->args[2]);

    if (!cablat_flags.hide_progress)
        fprintf(stderr, "\n%s\n", blat);

    system(blat); //Run fine BLAT

    free(blat_args);
    free(blat);
}


/*Takes in a vector of expansion structs for the expanded BLAT hits from a
  query and outputs them to the FASTA file CaBLAT_fine.fasta.*/
void write_fine_fasta(struct DSVector *oseqs){
    FILE *temp = fopen("CaBLAT_fine.fasta", "w");
    int i;

    if (!temp) {
        fprintf(stderr, "Could not open CaBLAT_fine.fasta for writing\n");
        return;
    }
    for (i = 0; i < oseqs->size; i++) {
        struct cb_seq *current_seq =
          ((struct cb_hit_expansion *)ds_vector_get(oseqs, i))->seq;
        fprintf(temp, "> %s\n%s\n", current_seq->name, current_seq->residues);
    }
    fclose(temp);

    if (cablat_flags.fine_blat_db)
        system("makeblastdb -dbtype nucl -in CaBLAT_fine.fasta -out "
               "CaBLAST_fine.fasta");
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

    conf = load_cablat_args();
    args = opt_config_parse(conf, argc, argv);

    if (args->nargs < 3) {
        fprintf(stderr, 
                "Usage: %s [flags] database-dir fasta-file results-file\n",
                argv[0]);
        opt_config_print_usage(conf);
        exit(1);
    }

    db = cb_database_r_init(args->args[0],
                            (cablat_flags.load_coarse_db ||
                             cablat_flags.load_coarse_residues),
                            (cablat_flags.load_coarse_db ||
                             cablat_flags.load_coarse_links),
                            cablat_flags.load_compressed_db,
                            cablat_flags.link_block_size);
    dbsize = read_int_from_file(8, db->coarse_db->db->file_params);

    blat_coarse(args);

    FILE *coarse_blat_output = fopen("coarse-blat.psl", "r");
    struct DSVector *coarse_hits = psl_read(coarse_blat_output);
    fclose(coarse_blat_output);

    struct DSVector *expanded_hits = expand_blat_hits(coarse_hits, db);
    write_fine_fasta(expanded_hits);

    blat_fine(args, dbsize);

    opt_args_free(args);
    opt_config_free(conf);

    return 0;
}
