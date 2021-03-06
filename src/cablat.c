#include <assert.h>
#include <errno.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "clibs/include/ds.h"
#include "clibs/include/opt.h"

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
 *for CaBLAT and returns a vector containing every original sequence section
 *re-created from the calls to cb_coarse_expand for the hits we are expanding.
 */
struct DSVector *expand_blat_hits(struct DSVector *hits,
                                  struct cb_database_r *db){
    struct DSVector *expanded_hits = ds_vector_create();

    for (int i = 0; i < hits->size; i++) {
        if (!cablat_flags.hide_progress) { //Display the progress bar
            int32_t digits_full = floor(log10((double)hits->size)),
                    digits_i    = floor(log10((double)i)),
                    spaces      = digits_full - digits_i;
            char *bar = progress_bar(i, hits->size);
            fprintf(stderr, "\r");
            fprintf(stderr, "Expanding coarse BLAT hit: %d/%d",
                            i+1, hits->size);
            for (int j = 0; j < spaces; j++)
                putc(' ', stderr);
            fprintf(stderr, " %s ", bar);
            free(bar);
        }

        struct psl_entry *h = (struct psl_entry *)ds_vector_get(hits, i);

        if (h->t_end - h->t_start > 1000)
            continue;

        int32_t coarse_start  = h->t_start, coarse_end = h->t_end,
                coarse_seq_id = atoi(h->t_name);

        struct DSVector *oseqs =
          cb_coarse_expand(db->coarse_db, db->com_db, coarse_seq_id,
                           coarse_start, coarse_end, 10);

        for (int j = 0; j < oseqs->size; j++)
            ds_vector_append(expanded_hits,
                             (struct cb_hit_expansion *)ds_vector_get(oseqs,j));

        ds_vector_free_no_data(oseqs);
    }
    if (!cablat_flags.hide_progress)
        fprintf(stderr, "\n");

    return expanded_hits;
}


/*Runs BLAT on the coarse FASTA file and stores the results in a temporary
  psl file.*/
void blat_coarse(char *db, char *queries){
    char *target_path = path_join(db, CABLAST_COARSE_FASTA),
         *coarse_blat_command =
            malloc(
              (strlen("blat    -noHead -minIdentity=80")
               +strlen(target_path)+strlen(queries)+strlen("coarse-blat.psl")+1)
              *sizeof(*coarse_blat_command));

    sprintf(coarse_blat_command,
            "blat %s %s %s -noHead -minIdentity=80",
            target_path, queries, "coarse-blat.psl");

    if (!cablat_flags.hide_progress)
        fprintf(stderr, "Running coarse BLAT\n%s\n\n", coarse_blat_command);

    system(coarse_blat_command);

    free(target_path);
    free(coarse_blat_command);
}

/*Takes in the name of the queries file, the destination filename for the fine
 *BLAT output, and a string of any BLAT arguments passed in and uses those
 *arguments to run BLAT on the fine FASTA file.
 */
void blat_fine(char *queries, char *out, char *blat_args){
    char *blat;
    int command_length = 1024;
    bool complete_psl = cablat_flags.complete_psl;

    blat = malloc(command_length*sizeof(*blat));
    assert(blat);

    if (blat_args[0] == '\0')
        sprintf(blat, "blat CaBLAT_fine.fasta %s %s %s",
                complete_psl ? "-noHead" : "", queries, out);
    else
        sprintf(blat, "blat %s %s CaBLAT_fine.fasta %s %s",
                blat_args, complete_psl ? "-out=psl -noHead" : "",
                queries, out);

    if (!cablat_flags.hide_progress)
        fprintf(stderr, "\nRunning fine BLAT\n%s\n", blat);

    system(blat); //Run fine BLAT

    free(blat_args);
    free(blat);
}

/*Takes in a vector of hit expansion structs for the expanded BLAT hits from
 *each coarse BLAT hit, a destination filename, and a bool telling whether or
 *not to include the offsets in the original sequence of each expanded hit in
 *the file and outputs the expanded hits to the specified file in FASTA format,
 *putting the offset at the end of the FASTA header if show_offsets is true.
 */
void write_fine_fasta(struct DSVector *oseqs, char *dest, bool show_offsets){
    FILE *temp;
    int sequences = 0;
    bool number_targets =
      cablat_flags.number_targets || cablat_flags.complete_psl;

    if (NULL == (temp = fopen(dest, "w"))) {
        fprintf(stderr,"fopen: 'fopen %s' failed: %s\n",dest,strerror(errno));
        exit(1);
    }

    for (int i = 0; i < oseqs->size; i++) {
        struct cb_hit_expansion *current_expansion =
            (struct cb_hit_expansion *)ds_vector_get(oseqs,i);

        struct cb_seq *current_seq = current_expansion->seq;
        int64_t offset             = current_expansion->offset;

        if (!number_targets) {
            if (show_offsets)
                fprintf(temp, "> %s (offset %ld)\n%s\n",
                               current_seq->name, offset,
                               current_seq->residues);
            else
                fprintf(temp, "> %s\n%s\n", current_seq->name,
                              current_seq->residues);
        }
        else {
            if (show_offsets)
                fprintf(temp, "> %d %s (offset %ld)\n%s\n", ++sequences,
                              current_seq->name, offset,
                              current_seq->residues);
            else
                fprintf(temp, "> %d %s\n%s\n", ++sequences,
                              current_seq->name, current_seq->residues);
        }
    }

    fclose(temp);
}

/*Takes in the name of a FASTA file and a destination filename and converts the
  FASTA file to a FASTA file at that destination with its headers numbered.*/
void write_numbered_fasta(char *fname, char *dest){
    assert(fname);
    assert(dest);

    FILE *f, *d;
    if (NULL == (f = fopen(fname, "r"))) {
        fprintf(stderr, "fopen: 'fopen %s' failed: %s\n",
                fname, strerror(errno));
        exit(1);
    }
    if (NULL == (d = fopen(dest, "w"))) {
        fprintf(stderr, "fopen: 'fopen %s' failed: %s\n",
                dest, strerror(errno));
        exit(1);
    }

    struct fasta_seq *s = NULL;
    int index = 1;

    s = fasta_read_next(f, "");
    while (s) {
        fprintf(d, ">%d %s\n%s\n", index, s->name, s->seq);
        fasta_free_seq(s);
        s = fasta_read_next(f, "");
        index++;
    }

    fclose(f);
    fclose(d);
}

/*Takes in the query FASTA file's name and returns a vector of fasta_seq
  structs for each query.*/
struct DSVector *read_queries(char *filename){
    FILE *query_file;
    if (NULL == (query_file = fopen(filename, "r"))) {
        fprintf(stderr, "fopen: 'fopen %s' failed: %s\n",
                filename, strerror(errno));
        exit(1);
    }

    struct DSVector *queries = ds_vector_create();
    struct fasta_seq *query = NULL;

    query = fasta_read_next(query_file, "");
    while (query) {
        ds_vector_append(queries, (void *)query);
        query = fasta_read_next(query_file, "");
    }

    fclose(query_file);

    return queries;
}

int main(int argc, char **argv){
    struct cb_database_r *db = NULL;
    struct opt_config *conf;
    struct opt_args *args;

    conf = load_cablat_args();
    args = opt_config_parse(conf, argc, argv);

    char *file_nq = "CaBLAT_numbered_queries.fasta",
         *file_fr = "CaBLAT_fine_results.psl";

    bool complete_psl = cablat_flags.complete_psl,
         number_queries = cablat_flags.number_queries || complete_psl;

    if (args->nargs < 3) {
        fprintf(stderr, 
                "Usage: %s [flags] database-dir fasta-file results-file "
                "[ --blat_args BLAT_ARGUMENTS ]\n", argv[0]);
        opt_config_print_usage(conf);
        exit(1);
    }

    //Load in the database
    db = cb_database_r_init(args->args[0],
                            (cablat_flags.load_coarse_db ||
                             cablat_flags.load_coarse_residues),
                            (cablat_flags.load_coarse_db ||
                             cablat_flags.load_coarse_links),
                            cablat_flags.load_compressed_db,
                            cablat_flags.link_block_size);

    //Number the queries if --complete-psl is passed in
    if (number_queries)
        write_numbered_fasta(args->args[1], file_nq);

    //Run coarse BLAT
    blat_coarse(args->args[0], number_queries ? file_nq : args->args[1]);

    struct DSVector *queries =
      complete_psl ? read_queries(args->args[1]) : NULL;

    FILE *coarse_blat_output;

    //Convert the coarse BLAT output to psl_entry structs
    if (NULL == (coarse_blat_output = fopen("coarse-blat.psl", "r"))) {
        fprintf(stderr, "fopen: 'fopen %s' failed: %s\n",
                        "coarse-blat.psl", strerror(errno));
        exit(1);
    }

    struct DSVector *coarse_hits = psl_read(coarse_blat_output);
    fclose(coarse_blat_output);

    /*Run hit expansion on the coarse BLAT hits to create the fine BLAT target
      file.*/
    struct DSVector *expanded_hits = expand_blat_hits(coarse_hits, db);
    write_fine_fasta(expanded_hits, "CaBLAT_fine.fasta", false);

    /*If an expanded FASTA file with offsets is requested, also output the
      expanded FASTA file*/
    if (strcmp(cablat_flags.output_expanded_fasta, "") != 0)
        write_fine_fasta(
          expanded_hits, cablat_flags.output_expanded_fasta, true);

    //Run fine BLAT
    blat_fine(number_queries ? file_nq : args->args[1],
              complete_psl ? file_fr : args->args[2], get_blat_args(args));

    //Convert the output to a complete .psl file if --complete-psl is passed
    if (cablat_flags.complete_psl) {
        FILE *fine_blat_output, *output_file;

        if (NULL == (fine_blat_output = fopen(file_fr, "r"))) {
            fprintf(stderr, "fopen: 'fopen %s' failed: %s\n",
                            file_fr, strerror(errno));
            exit(1);
        }

        if (NULL == (output_file = fopen(args->args[2], "w"))) {
            fprintf(stderr, "fopen: 'fopen %s' failed: %s\n",
                            args->args[2], strerror(errno));
            exit(1);
        }

        struct DSVector *fine_hits = psl_read(fine_blat_output);
        int64_t *seq_lengths = cb_compressed_get_lengths(db->com_db);
        for (int i = 0; i < fine_hits->size; i++){
            struct psl_entry *hit =
              (struct psl_entry *)ds_vector_get(fine_hits, i);
            int target_i = atoi(hit->t_name)-1,
                query_i  = atoi(hit->q_name)-1;
            struct cb_hit_expansion *target_expansion =
              (struct cb_hit_expansion *)ds_vector_get(expanded_hits, target_i);

            free(hit->t_name);
            hit->t_name =
              malloc((strlen(target_expansion->seq->name)+1)
                     *sizeof(*(hit->t_name)));
            assert(hit->t_name);
            strcpy(hit->t_name, target_expansion->seq->name);

            struct fasta_seq *query_seq =
              (struct fasta_seq *)ds_vector_get(queries, query_i);
            char *query_name = query_seq->name;
            free(hit->q_name);
            hit->q_name = malloc((strlen(query_name)+1)*sizeof(*(hit->q_name)));
            assert(hit->q_name);
            strcpy(hit->q_name, query_name);

            hit->t_size = (unsigned)(seq_lengths[target_expansion->seq->id]);
            hit->t_start += target_expansion->offset;
            hit->t_end += target_expansion->offset;

            for (int j = 0; j < hit->block_count; j++)
                hit->t_starts[j] += target_expansion->offset;

            psl_entry_print(hit, output_file);
            psl_entry_free(hit);
        }

        ds_vector_free_no_data(fine_hits);
        free(seq_lengths);
        fclose(fine_blat_output);
        fclose(output_file);
    }

    if (complete_psl) {
        for (int i = 0; i < queries->size; i++)
            fasta_free_seq((struct fasta_seq *)ds_vector_get(queries, i));
        ds_vector_free_no_data(queries);
    }

    /*Free the expanded hits and the database and delete the intermediate files
      unless --no-cleanup is passed.*/
    for (int i = 0; i < expanded_hits->size; i++)
        cb_hit_expansion_free(
          (struct cb_hit_expansion *)ds_vector_get(expanded_hits, i));
    ds_vector_free_no_data(expanded_hits);

    for (int i = 0; i < coarse_hits->size; i++)
        psl_entry_free((struct psl_entry *)ds_vector_get(coarse_hits, i));
    ds_vector_free_no_data(coarse_hits);

    cb_database_r_free(db);
    opt_args_free(args);
    opt_config_free(conf);

    if (!cablat_flags.no_cleanup) {
        system("rm coarse-blat.psl");
        system("rm CaBLAT_fine.fasta");

        if (!cablat_flags.no_cleanup) {
            if (number_queries)
                system("rm CaBLAT_numbered_queries.fasta");
            if (complete_psl)
                system("rm CaBLAT_fine_results.psl");
        }
    }

    return 0;
}
