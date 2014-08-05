#ifndef __CABLAT_PSL_H__
#define __CABLAT_PSL_H__

#include <stdio.h>
#include <stdlib.h>

#include "clibs/include/ds.h"

/*Summary info about a patSpace alignment, modified from psl.h in Jim Kent's
  BLAT.*/
struct psl_entry {
    int matches;
    int mismatches;
    int rep_matches;
    int n_count;
    int q_num_insert;
    int q_base_insert;
    int t_num_insert;
    int t_base_insert;
    char *strand;
    char *q_name;
    unsigned q_size;
    int q_start;
    int q_end;
    char *t_name;
    unsigned t_size;
    int t_start;
    int t_end;
    unsigned block_count;
    unsigned *block_sizes;
    unsigned *q_starts;
    unsigned *t_starts;
};

struct psl_entry *
psl_entry_init(int matches, int mismatches, int rep_matches, int n_count,
               int q_num_insert, int q_base_insert, int t_num_insert,
               int t_base_insert, char *strand, char *q_name, unsigned q_size,
               int q_start, int q_end, char *t_name, unsigned t_size,
               int t_start, int t_end, unsigned block_count, char *block_sizes,
               char *q_starts, char *t_starts);

void psl_entry_free(struct psl_entry *entry);

void psl_entry_print(struct psl_entry *p, FILE *f);

struct psl_entry *psl_load(char **row);
struct DSVector *psl_read(FILE *f);

#endif
