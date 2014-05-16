#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "align.h"

struct {
    char *rseq;
    char *oseq;
    int32_t answer;
} tests[] = {
    { "A", "A", 0 },
    {"AB", "AB", 0},
    {"ABC", "ABC", 3},
    {"ABCD", "ABCD", 3},
    {"ABCYEFG", "ABCZEFG", 3},
    {"ABCYEFGH", "ABCZEFGH", 8},
    {"ABCDEFGHIJKLMNOP", "ABCDEFGHIJKLMNOP", 15},
    {"ABCDEF", "ABC", 3},
    {"ABC", "ABCDEF", 3},
    {"ABCDEFGHIKLMNPQR", "ABCDEFGHIKLMNPQR", 15},
    { NULL, NULL, 0 }
};

int main(void)
{
    int32_t window_size, kmer_size, id_threshold;
    int32_t tval;
    int i;

    window_size = 10;
    kmer_size = 3;
    id_threshold = 50;

    for (i = 0; tests[i].rseq != NULL; i++) {
        tval = cbp_align_ungapped(
            window_size, kmer_size, id_threshold,
            tests[i].rseq, 0, strlen(tests[i].rseq),
            tests[i].oseq, 0, strlen(tests[i].oseq));

        if (tval != tests[i].answer) {
            printf("TEST %d FAILED\n", i);
            printf("Ungapped extension on '%s' and '%s' should yield a "
                   "length of %d, but cbp_align_ungapped returned %d.\n",
                   tests[i].rseq, tests[i].oseq, tests[i].answer, tval);
            exit(1);
        }
    }

    printf("ALL TESTS PASSED\n");

    return 0;
}

