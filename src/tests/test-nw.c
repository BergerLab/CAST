#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "align.h"

struct {
    char *ref;
    char *org;
    char *aref;
    char *aorg;
    int length;
} tests[] = {
    {
        "GHIKLMNPQR",
        "GAAAHIKLMN",
        "---GHIKLMNPQR",
        "GAAAHIKLMN---",
        13
    },
    {
        "ABCD",
        "ABCD",
        "ABCD",
        "ABCD",
        4
    },
    {
        "GHIKLMNPQRSTVW",
        "GAAAHIKLMNPQRSTVW",
        "---GHIKLMNPQRSTVW",
        "GAAAHIKLMNPQRSTVW",
        17
    },
    {
        "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
        "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
        "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
        "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
        62
    },
    {
        "GHIKLMNPQRSTVW",
        "GAAAHIKLMNPQRSTVW",
        "---GHIKLMNPQRSTVW",
        "GAAAHIKLMNPQRSTVW",
        17
    },
    /* {  */
        /* "ABCDEFGWXYZ",  */
        /* "ABCDEFMNPQRSTZABEGWXYZ",  */
        /* "ABCDEF-----------GWXYZ",  */
        /* "ABCDEFMNPQRSTZABEGWXYZ",  */
        /* 22 */
    /* },  */
    { NULL, NULL, NULL, NULL, 0 }
};

int main(void)
{
    struct cbp_align_nw_memory *mem;
    struct cbp_alignment alignment;
    int i;

    mem = cbp_align_nw_memory_init();

    for (i = 0; tests[i].ref != NULL; i++) {
        alignment = cbp_align_nw(mem,
            tests[i].ref, 0, strlen(tests[i].ref),
            tests[i].org, 0, strlen(tests[i].org));

        if (0 != strcmp(alignment.ref, tests[i].aref)
            || 0 != strcmp(alignment.org, tests[i].aorg)
            || alignment.length != tests[i].length) {
            printf("TEST %d FAILED\n", i);
            printf("The alignment for:\n");
            printf("---------------------------\n");
            printf("\t%s\n", tests[i].ref);
            printf("\t%s\n", tests[i].org);
            printf("---------------------------\n");
            printf("produced:\n");
            printf("---------------------------\n");
            printf("\tlength: %d\n", alignment.length);
            printf("\t%s\n", alignment.ref);
            printf("\t%s\n", alignment.org);
            printf("---------------------------\n");
            printf("but should have been:\n");
            printf("---------------------------\n");
            printf("\tlength: %d\n", tests[i].length);
            printf("\t%s\n", tests[i].aref);
            printf("\t%s\n", tests[i].aorg);
            printf("---------------------------\n");
            exit(1);
        }
    }

    cbp_align_nw_memory_free(mem);

    printf("ALL TESTS PASSED\n");

    return 0;
}

