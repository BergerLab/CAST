#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "clibs/include/ds.h"

#include "psl.h"
#include "util.h"

struct psl_entry *
psl_entry_init(int matches, int mismatches, int rep_matches, int n_count,
               int q_num_insert, int q_base_insert, int t_num_insert,
               int t_base_insert, char *strand, char *q_name, unsigned q_size,
               int q_start, int q_end, char *t_name, unsigned t_size,
               int t_start, int t_end, unsigned block_count,
               unsigned block_size, unsigned q_starts, unsigned t_starts){
    struct psl_entry *entry = malloc(sizeof(*entry));
    assert(entry);

    entry->matches       = matches;
    entry->mismatches    = mismatches;
    entry->rep_matches   = rep_matches;
    entry->n_count       = n_count;
    entry->q_num_insert  = q_num_insert;
    entry->q_base_insert = q_base_insert;
    entry->t_num_insert  = t_num_insert;
    entry->t_base_insert = t_base_insert;

    entry->strand        = malloc((strlen(strand)+1)*sizeof(*(entry->strand)));
    assert(entry->strand);

    strcpy(entry->strand, strand);

    entry->q_name        = malloc((strlen(q_name)+1)*sizeof(*(entry->q_name)));
    assert(entry->q_name);

    strcpy(entry->q_name, q_name);

    entry->q_size        = q_size;
    entry->q_start       = q_start;
    entry->q_end         = q_end;

    entry->t_name        = t_name;
    entry->t_name        = malloc((strlen(t_name)+1)*sizeof(*(entry->t_name)));
    assert(entry->t_name);

    strcpy(entry->t_name, t_name);

    entry->t_size        = t_size;
    entry->t_start       = t_start;
    entry->t_end         = t_end;
    entry->block_count   = block_count;
    entry->block_size    = block_size;
    entry->q_starts      = q_starts;
    entry->t_starts      = t_starts;

    return entry;
}

//Frees a psl entry struct
void psl_entry_free(struct psl_entry *entry){
    free(entry->strand);
    free(entry->q_name);
    free(entry->t_name);
    free(entry);
}

/*Load a psl entry from a row parsed with split_space, based on psl_load in Jim
  Kent's BLAT.*/
struct psl_entry *psl_load(char **row){
    return psl_entry_init(
      atoi(row[0]), atoi(row[1]), atoi(row[2]), atoi(row[3]), atoi(row[4]),
      atoi(row[5]), atoi(row[6]), atoi(row[7]), row[8], row[9], atoi(row[10]),
      atoi(row[11]), atoi(row[12]), row[13], atoi(row[14]), atoi(row[15]),
      atoi(row[16]), atoi(row[17]), atoi(row[18]), atoi(row[19]),
      atoi(row[20]));
}

/*Takes in a file pointer for a .psl file and parses each of its lines with
  psl_load to create a psl_entry struct, which is stored in a vector.*/
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
