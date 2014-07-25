#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "clibs/include/ds.h"

#include "psl.h"
#include "util.h"

/*Load a psl entry from a row parsed with split_space, based on psl_load in Jim
  Kent's BLAT.*/
struct psl_entry *psl_load(char **row){

    struct psl_entry *ret = malloc(sizeof(*ret));
    assert(ret);

    ret->matches       = atoi(row[0]);
    ret->mismatches    = atoi(row[1]);
    ret->rep_matches   = atoi(row[2]);
    ret->n_count       = atoi(row[3]);
    ret->q_num_insert  = atoi(row[4]);
    ret->q_base_insert = atoi(row[5]);
    ret->t_num_insert  = atoi(row[6]);
    ret->t_base_insert = atoi(row[7]);

    ret->strand        = malloc((strlen(row[8])+1)*sizeof(*(ret->strand)));
    assert(ret->strand);

    strcpy(ret->strand, row[8]);

    ret->q_name        = malloc((strlen(row[9])+1)*sizeof(*(ret->q_name)));
    assert(ret->q_name);

    strcpy(ret->q_name, row[9]);

    ret->q_size        = atoi(row[10]);
    ret->q_start       = atoi(row[11]);
    ret->q_end         = atoi(row[12]);

    ret->t_name        = row[13];
    ret->t_name        = malloc((strlen(row[13])+1)*sizeof(*(ret->q_name)));
    assert(ret->t_name);

    strcpy(ret->t_name, row[13]);

    ret->t_size        = atoi(row[14]);
    ret->t_start       = atoi(row[15]);
    ret->t_end         = atoi(row[16]);
    ret->block_count   = atoi(row[17]);
    ret->block_size    = atoi(row[18]);
    ret->q_starts      = atoi(row[19]);
    ret->t_starts      = atoi(row[20]);

    return NULL;
}

struct DSVector *psl_read(FILE *f){
    struct DSVector *entries = ds_vector_create();
    char *line = NULL;

    while (0 != (readline(f, &line))) {
        char *no_newline = trim_space(line);
        struct psl_entry *entry = psl_load(split_spaces(no_newline));
        ds_vector_append(entries, (void *)entry);
    }

    return entries;
}
