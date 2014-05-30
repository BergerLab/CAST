#ifndef __CABLAST_ALIGN_H__
#define __CABLAST_ALIGN_H__

/*
This is adapted from Dan Kortschak's Needleman-Wunsch algorithm from the biogo
package: code.google.com/p/biogo.

It's mostly copied from its original form, but it is optimized specifically
for cablast to limit allocations and to absolve the need for the biogo/seq.Seq
type.
*/

#include <stdint.h>

#define CABLAST_ALIGN_SEQ_SIZE 10000

struct ungapped_alignment{
    int32_t length;
    bool found_bad_window;
};

struct ungapped_alignment
cb_align_ungapped(char *rseq, int32_t rstart, int32_t rend, int32_t dir1,
                  int32_t i1, char *oseq, int32_t ostart, int32_t oend,
                  int32_t dir2, int32_t i2, bool *matches,
                  bool *matches_past_clump, int *matches_index);

struct cb_align_nw_memory {
    int32_t *table;
    int32_t *zeroes;
    char *ref;
    char *org;
};

struct cb_nw_tables {
    int **dp_score;
    int **dp_from;
};

int *best_edge(int **dp_score, int dp_len1, int dp_len2);
int *backtrack_to_clump(struct cb_nw_tables tables, int *pos);

struct cb_nw_tables
make_nw_tables(char *rseq, int dp_len1, int i1, int dir1,
               char *oseq, int dp_len2, int i2, int dir2);


struct cb_align_nw_memory *
cb_align_nw_memory_init();

void
cb_align_nw_memory_free(struct cb_align_nw_memory *mem);

struct cb_alignment {
    char *ref;
    char *org;
    int32_t length;
};

struct cb_alignment
cb_align_nw(struct cb_align_nw_memory *mem,
            char *rseq, int dp_len1, int i1, int dir1,
            char *oseq, int dp_len2, int i2, int dir2,
            bool *matches, int *matches_index);

int32_t
cb_align_length_nogaps(char *residues);

int32_t
attempt_ext(int32_t i1, const int32_t dir1, const char *s1, int32_t len1,
            int32_t start1, int32_t i2, const int32_t dir2, const char *s2,
            int32_t len2, int32_t start2);

int check_and_update(bool *matches, int *matches_index, int *num_matches,
                     bool *temp, int temp_index);

int max_dp_len(int i, int dir, int len);
#endif
