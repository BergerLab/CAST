#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "DNAutils.h"
#include "align.h"
#include "DNAalphabet.h"
#include "flags.h"

struct ungapped_alignment
cb_align_ungapped(char *rseq, int32_t rstart, int32_t rend, int32_t dir1,
                   int32_t i1, char *oseq, int32_t ostart, int32_t oend,
                   int32_t dir2, int32_t i2, bool *matches,
                   bool *matches_past_clump, int *matches_index)
{
    int32_t length, scanned, successive;
    int32_t rlen, olen;
    int32_t i;
    int32_t matches_since_last_consec;
    int consec_match_clump_size;
    int32_t dir_prod;
    int matches_count;
    int temp_index;
    struct ungapped_alignment ungapped;

    ungapped.length = -1;
    ungapped.found_bad_window = false;
    length = 0;
    scanned = 0;
    consec_match_clump_size = compress_flags.consec_match_clump_size;
    successive = consec_match_clump_size;
    rlen = rend - rstart;
    olen = oend - ostart;
    dir_prod = dir1 * dir2;
    temp_index = 0;
    matches_count = 0;
    matches_since_last_consec = 0;

    for (i = *matches_index - 100; i < *matches_index; i++)
        if (matches[i])
            matches_count++;
    while (i1 >= rstart && i1 < rend && i2 >= ostart && i2 < oend) {
        int cur_ismatch = bases_match(rseq[i1], oseq[i2], dir_prod);
        i1 += dir1;
        i2 += dir2;
        scanned++;
        if (cur_ismatch == 1) {
            matches_past_clump[temp_index] = true;
            temp_index++;
            successive++;
            if (successive >= consec_match_clump_size) {
                int update = check_and_update(matches, matches_index, 
                                              &matches_count,
                                              matches_past_clump, temp_index);
                    length += update;
                    if (update != temp_index) {
                        ungapped.length = length;
                        ungapped.found_bad_window = true;
                        return ungapped;
                    }
                temp_index = 0;
                matches_since_last_consec = 0;
            }
            else
                matches_since_last_consec++;
        }
        else {
            matches_past_clump[temp_index] = false;
            temp_index++;
            successive = 0;
            if (scanned - length >= compress_flags.btwn_match_min_dist_check) {
                if ((double)matches_since_last_consec < (scanned-length)*0.5) {
                    ungapped.length = length;
                    return ungapped;
                }
            }
        }
    }
    ungapped.length = scanned;
    return ungapped;
}

int32_t
cb_align_identity(char *rseq, int32_t rstart, int32_t rend,
                   char *oseq, int32_t ostart, int32_t oend)
{
    int32_t rlen, olen, same, i;

    rlen = rend - rstart;
    olen = oend - ostart;
    assert(rlen == olen);

    if (rlen == 0 && olen == 0)
        return 0;

    same = 0;
    for (i = 0; i < rlen; i++)
        if (rseq[rstart + i] == oseq[ostart + i])
            same++;

    return (same * 100) / rlen;
}

struct cb_align_nw_memory *
cb_align_nw_memory_init()
{
    struct cb_align_nw_memory *mem;
    int seq_size = CABLAST_ALIGN_SEQ_SIZE;

    mem = malloc(sizeof(*mem));
    assert(mem);

    mem->table = malloc(seq_size * seq_size * sizeof(*mem->table));
    assert(mem->table);

    mem->zeroes = malloc(seq_size * seq_size * sizeof(*mem->zeroes));
    assert(mem->zeroes);

    mem->ref = malloc(seq_size * sizeof(*mem->ref));
    assert(mem->ref);

    mem->org = malloc(seq_size * sizeof(*mem->org));
    assert(mem->org);

    return mem;
}

void
cb_align_nw_memory_free(struct cb_align_nw_memory *mem)
{
    free(mem->zeroes);
    free(mem->table);
    free(mem->ref);
    free(mem->org);
    free(mem);
}

/*Makes the tables used in Needleman-Wunsch alignment; takes in two sequences,
 *the lengths of the sections of the sequences that we are aligning, indices
 *into these sequences, and the directions of the sequences and returns a
 *struct containing a table of the scores in the alignment and a table of
 *directions for backtracking to the start of the alignment.
 */
struct cb_nw_tables
make_nw_tables(char *rseq, int dp_len1, int i1, int dir1,
               char *oseq, int dp_len2, int i2, int dir2)
{
    struct cb_nw_tables tables;
    int i, j1, j2;
    int dir_prod = dir1*dir2;
    int **dp_score = malloc((dp_len1+1)*sizeof(*dp_score));
    int **dp_from = malloc((dp_len1+1)*sizeof(*dp_from));
    for (i = 0; i < dp_len1+1; i++) {
        dp_score[i] = malloc((dp_len2+1)*sizeof(**dp_score));
        dp_from[i] = malloc((dp_len2+1)*sizeof(**dp_from));
    }
    for (i = 0; i <= dp_len2; i++) {
        dp_score[0][i] = -3*i;
        dp_from[0][i] = 2;
    }
    for (i = 1; i <= dp_len1; i++) {
        dp_score[i][0] = -3*i;
        dp_from[i][0] = 1;
    }
    for (j1 = 1; j1 <= dp_len1; j1++)
        for (j2 = 1; j2 <= dp_len2; j2++){
            int score0, score1, score2;
            score0 = dp_score[j1-1][j2-1] +
                     (bases_match(rseq[i1+dir1*(j1-1)], oseq[i2+dir2*(j2-1)],
                                                         dir_prod) ? 1 : -3);
            score1 = dp_score[j1-1][j2] - 3;
            score2 = dp_score[j1][j2-1] - 3;
            if (score0 >= score1 && score0 >= score2) {
                dp_score[j1][j2] = score0;
                dp_from[j1][j2] = 0;
            }
            else if (score2 >= score1) {
                dp_score[j1][j2] = score2;
                dp_from[j1][j2] = 2;
            }
            else {
                dp_score[j1][j2] = score1;
                dp_from[j1][j2] = 1;
            }
        }
    tables.dp_score = dp_score;
    tables.dp_from = dp_from;
    return tables;
}

/*Finds the space on the bottom and right edges of a Needleman-Wunsch score
 *table with the best score, returning it as an array of two ints.
 */
int *best_edge(int **dp_score, int dp_len1, int dp_len2){
    int j1, j2;
    int max_dp_score = -1000;
    int *best = malloc(2*sizeof(*best));
    best[0] = -1;
    best[1] = -1;
    for (j2 = 0; j2 <= dp_len2; j2++){
        if (dp_score[dp_len1][j2] >= max_dp_score) {
            max_dp_score = dp_score[dp_len1][j2];
            best[0] = dp_len1; best[1] = j2;
        }}
    for (j1 = 0; j1 <= dp_len1; j1++)
        if (dp_score[j1][dp_len2] >= max_dp_score) {
            max_dp_score = dp_score[j1][dp_len2];
            best[0] = j1; best[1] = dp_len2;
        }
   return best; 
}

int *backtrack_to_clump(struct cb_nw_tables tables, int *pos){
    int consec_matches = 0;
    int **dp_score = tables.dp_score;
    int **dp_from = tables.dp_from;
    int consec_match_clump_size = compress_flags.consec_match_clump_size;
    while (!(pos[0] == 0 && pos[1] == 0)) {
        int prev_j1, prev_j2;
        if (consec_matches == consec_match_clump_size) { /*found chunk; stop*/
            pos[0] += consec_match_clump_size;
            pos[1] += consec_match_clump_size;
            break;
        }

        switch (dp_from[pos[0]][pos[1]]) { /*backtrack to previous cell*/
            case 0: prev_j1 = pos[0]-1; prev_j2 = pos[1]-1; break;
            case 2: prev_j1 = pos[0]; prev_j2 = pos[1]-1;break;
            default: prev_j1 = pos[0]-1; prev_j2 = pos[1];
        }
        if (dp_from[pos[0]][pos[1]] == 0)
            if (dp_score[pos[0]][pos[1]] > dp_score[prev_j1][prev_j2]) /*match*/
                consec_matches++;
            else
                consec_matches = 0;
        else
            consec_matches = 0;
        pos[0] = prev_j1;
        pos[1] = prev_j2;
    }
    /*Couldn't find a 4-mer clump*/
    if (consec_matches < consec_match_clump_size) {
        pos[0] = 0;
        pos[1] = 0;
    }
    return pos;
}

struct cb_alignment
cb_align_nw(struct cb_align_nw_memory *mem,
             char *rseq, int dp_len1, int i1, int dir1,
             char *oseq, int dp_len2, int i2, int dir2,
             bool *matches, int *matches_index)
{
    struct cb_alignment align;
    int matches_count = 0, i = 0;
    struct cb_nw_tables tables = make_nw_tables(rseq, dp_len1, i1, dir1,
                                                 oseq, dp_len2, i2, dir2);
    int *best = best_edge(tables.dp_score, dp_len1, dp_len2);

    best = backtrack_to_clump(tables, best);

    if (best[0] <= 0) {
        align.ref = "\0";
        align.org = "\0";
        align.length = -1;
        free(best);
        for (i = 0; i <= dp_len1; i++) {
            free(tables.dp_score[i]);
            free(tables.dp_from[i]);
        }
        free(tables.dp_score);
        free(tables.dp_from);
       
        return align;
    }
    int cur_j1 = best[0];
    int cur_j2 = best[1];
    int dir_prod = dir1 * dir2;
    int **dp_score = tables.dp_score;
    int **dp_from = tables.dp_from;
    bool *matches_to_add = malloc((cur_j1 + cur_j2)*sizeof(*matches_to_add));
    char *subs1_dp = malloc((cur_j1 + cur_j2)*sizeof(*subs1_dp));
    char *subs2_dp = malloc((cur_j1 + cur_j2)*sizeof(*subs2_dp));
    int num_steps = 0;
    align.ref = "\0";
    align.org = "\0";
    align.length = -1;
    while (!(cur_j1 == 0 && cur_j2 == 0)) {
        int prev_j1, prev_j2;
        switch (dp_from[cur_j1][cur_j2]) {
            char c1, c2;
        case 0:
            prev_j1 = cur_j1-1; prev_j2 = cur_j2-1; /*match or substitution*/
            c1 = rseq[i1+dir1*prev_j1]; /*comp if antisense*/
            c2 = oseq[i2+dir2*prev_j2];
            if (dir_prod == -1) c2 = base_complement(c2);
            subs1_dp[num_steps] = c1;
            subs2_dp[num_steps] = c2;
            break;
        case 2: prev_j1 = cur_j1; prev_j2 = cur_j2-1; /*advance 2; gap in 1*/
            c2 = oseq[i2+dir2*prev_j2];
            if (dir_prod == -1) c2 = base_complement(c2); /*comp if antisense*/
            subs1_dp[num_steps] = '-';
            subs2_dp[num_steps] = c2;
            break;
        default: prev_j1 = cur_j1-1; prev_j2 = cur_j2; /*advance 1; gap in 2*/
            c1 = rseq[i1+dir1*prev_j1];
            subs1_dp[num_steps] = c1;
            subs2_dp[num_steps] = '-';
        }
        matches_to_add[num_steps] = dp_score[cur_j1][cur_j2] >
                                    dp_score[prev_j1][prev_j2];
        num_steps++;
        cur_j1 = prev_j1; cur_j2 = prev_j2;
    }
    for (i = 0; i < num_steps/2; i++) { /* flip order */
        bool temp = matches_to_add[num_steps-i-1];
        matches_to_add[num_steps-1-i] = matches_to_add[i];
        matches_to_add[i] = temp;
    }

    /*note: need to flip order*/
    if (dp_len1 < compress_flags.min_match_len &&
        dp_len2 < compress_flags.min_match_len)
        for (i = *matches_index - 100; i < *matches_index; i++)
            if (matches[i])
                matches_count++;

    /*Make sure we don't have a bad window unless we are running
      Needleman-Wunsch alignment on a match.  If we have a bad window, then
      throw out this alignment.  Otherwise, copy the alignment into align.org
      and align.ref.*/
    if (dp_len1 < compress_flags.min_match_len &&
        dp_len2 < compress_flags.min_match_len &&
        check_and_update(matches, matches_index, &matches_count,
                         matches_to_add, num_steps) != num_steps)
        align.length = -1;
    else {
        align.length = num_steps;
        align.org = malloc((align.length+1)*sizeof(*(align.org)));
        align.ref = malloc((align.length+1)*sizeof(*(align.ref)));
        for (i = 0; i < align.length; i++) {
            /*Don't update the matches array if we are running Needleman-Wunsch
              alignment on a match.*/
            if (dp_len1 < compress_flags.min_match_len &&
                dp_len2 < compress_flags.min_match_len)
                matches[(*matches_index)+i] = matches_to_add[i];

            align.ref[i] = subs1_dp[align.length-i-1];
            align.org[i] = subs2_dp[align.length-i-1];
        }
        align.org[align.length] = '\0';
        align.ref[align.length] = '\0';
    }
    free(best);
    for (i = 0; i <= dp_len1; i++) {
        free(tables.dp_score[i]);
        free(tables.dp_from[i]);
    }
    free(tables.dp_score);
    free(tables.dp_from);
    free(subs1_dp);
    free(subs2_dp);
    free(matches_to_add);
    return align;
}

/*Returns the number of non-gap characters in a string*/
int32_t
cb_align_length_nogaps(char *residues)
{
    int i, len, rlen;

    len = 0;
    rlen = strlen(residues);
    for (i = 0; i < rlen; i++)
        if (residues[i] != '-')
            len++;

    return len;
}

/* returns extension distance (but does not move pointers)*/
int32_t
attempt_ext(int32_t i1, const int32_t dir1, const char *s1, int32_t len1,
            int32_t start1, int32_t i2, const int32_t dir2, const char *s2,
            int32_t len2, int32_t start2)
{
    const int32_t dir_prod = dir1*dir2;
    int32_t progress = 0;
    int32_t consec_mismatch = 0;
    i1 += dir1;
    i2 += dir2;
/*printf("attempt_ext: i1: %d, i2: %d, progress: ", i1-start1, i2-start2);*/
    /*Replace this 3 with the flag for max_consec_mismatch*/
    while (consec_mismatch < 3 &&
           i1 >= start1 && i1 < start1+len1 &&
           i2 >= start2 && i2 < start2+len2) {
        if (!bases_match(s1[i1], s2[i2], dir_prod))
            consec_mismatch++;
        else
            consec_mismatch = 0;
        i1 += dir1; i2 += dir2;
        progress++;
    }
/*fprintf(stderr, "attempt_ext, progress = %d\n", progress);*/
    return progress;
}

/*Takes in as input an array of bools representing data on whether or not
 *previous pairs of bases matched, an index into that array, the current number
 *of matches, an array of bools to add to the array of matches, and the number
 *of bools to add.  check_and_update adds bools from the temp array to the
 *matches array until it either has added all of the matches or a bad window
 *was found.  If a bad window was found, check_and_update returns false.
 *Otherwise, check_and_update returns true.
 */
int check_and_update(bool *matches, int *matches_index, int *num_matches,
                                             bool *temp, int temp_index){
    int i;
    for (i = 0; i < temp_index; i++) {
        int hundred_bases_ago = *matches_index - 100;
        matches[(*matches_index)] = temp[i];
        if (temp[i])
            (*num_matches)++;
        if (matches[hundred_bases_ago])
            (*num_matches)--;
        (*matches_index)++;
        if (*num_matches < 85)
            return i;
    }
    return temp_index;
}

int min(int a, int b){return a<b?a:b;}

int max_dp_len(int i, int dir, int len){
    return dir == 1 ? min(25, len-i) : min(25, i+1);
}
