#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "clibs/include/ds.h"

#include "psl.h"
#include "util.h"

//Takes in the data for a .psl entry and uses it to create a psl_entry struct. 
struct psl_entry *
psl_entry_init(int matches, int mismatches, int rep_matches, int n_count,
               int q_num_insert, int q_base_insert, int t_num_insert,
               int t_base_insert, char *strand, char *q_name, unsigned q_size,
               int q_start, int q_end, char *t_name, unsigned t_size,
               int t_start, int t_end, unsigned block_count, char *block_sizes,
               char *q_starts, char *t_starts){
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

    entry->t_name        = malloc((strlen(t_name)+1)*sizeof(*(entry->t_name)));
    assert(entry->t_name);
    strcpy(entry->t_name, t_name);

    entry->t_size        = t_size;
    entry->t_start       = t_start;
    entry->t_end         = t_end;
    entry->block_count   = block_count;

    char **block_size_array = split_char(block_sizes, ','),
         **q_starts_array   = split_char(q_starts, ','),
         **t_starts_array   = split_char(t_starts, ',');

    unsigned *block_size_ints = malloc(block_count*sizeof(*block_size_ints));
    assert(block_size_ints);

    unsigned *q_starts_ints = malloc(block_count*sizeof(*q_starts_ints));
    assert(q_starts_ints);

    unsigned *t_starts_ints = malloc(block_count*sizeof(*t_starts_ints));
    assert(t_starts_ints);

    for (unsigned i = 0; i < block_count; i++) {
        block_size_ints[i] = (unsigned)atoi(block_size_array[i]);
        q_starts_ints[i]   = (unsigned)atoi(q_starts_array[i]);
        t_starts_ints[i]   = (unsigned)atoi(t_starts_array[i]);
    }

    entry->block_sizes   = block_size_ints;
    entry->q_starts      = q_starts_ints;
    entry->t_starts      = t_starts_ints;

    for (int i = 0; block_size_array[i] != NULL; i++)
        free(block_size_array[i]);
    free(block_size_array);

    for (int i = 0; q_starts_array[i] != NULL; i++)
        free(q_starts_array[i]);
    free(q_starts_array);

    for (int i = 0; t_starts_array[i] != NULL; i++)
        free(t_starts_array[i]);
    free(t_starts_array);

    return entry;
}

//Frees a psl entry struct
void psl_entry_free(struct psl_entry *entry){
    free(entry->strand);
    free(entry->q_name);
    free(entry->t_name);
    free(entry->block_sizes);
    free(entry->q_starts);
    free(entry->t_starts);
    free(entry);
}

/*Load a psl entry from a row parsed with split_space, based on psl_load in Jim
  Kent's BLAT.*/
struct psl_entry *psl_load(char **row){
    return psl_entry_init(
      atoi(row[0]), atoi(row[1]), atoi(row[2]), atoi(row[3]), atoi(row[4]),
      atoi(row[5]), atoi(row[6]), atoi(row[7]), row[8], row[9], atoi(row[10]),
      atoi(row[11]), atoi(row[12]), row[13], atoi(row[14]), atoi(row[15]),
      atoi(row[16]), atoi(row[17]), row[18], row[19], row[20]);
}

/*Takes in a file pointer for a .psl file and parses each of its lines with
  psl_load to create a psl_entry struct, which is stored in a vector.*/
struct DSVector *psl_read(FILE *f){
    struct DSVector *entries = ds_vector_create();
    char *line = NULL;

    while (0 != (readline(f, &line))) {
        char *no_newline = trim_space(line);
        char **psl_data = split_char(no_newline, '\t');
        struct psl_entry *entry = psl_load(psl_data);

        ds_vector_append(entries, (void *)entry);

        for (int i = 0; i < 21; i++)
            free(psl_data[i]);
        free(psl_data);
        free(no_newline);
    }

    free(line);

    return entries;
}
