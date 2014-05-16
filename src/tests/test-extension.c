#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "align.h"
#include "flags.h"

struct extend_match {
    int32_t rlen;
    int32_t olen;
};

extern struct extend_match
extend_match(struct cbp_align_nw_memory *mem,
             char *rseq, int32_t rstart, int32_t rend,
             char *oseq, int32_t ostart, int32_t oend);

struct {
    char *ref;
    char *org;
    int32_t rlen;
    int32_t olen;
} tests[] = {
    {
        "ABCDEFGHIKLMNPQR",
        "ABCDEFGHIKLMNPQR",
        16,
        16
    },
    {
        "ABCDEFGHIKLMNPQRSTVW",
        "ABCDEFGAAAHIKLMNPQRSTVW",
        20,
        23
    },
    {
        "ABCDEFGHIKLMNPQRSTVW",
        "ABCDEFGAAAHIKLMNPQRSTBBBBBBBBBBBBBBBBBBBVW",
        6,
        6
    },
    { NULL, NULL, 0, 0 }
};

int32_t main(int32_t argc, char **argv)
{
    struct extend_match ext;
    struct cbp_align_nw_memory *mem;
    struct opt_config *conf;
    struct opt_args *args;
    int i;

    mem = cbp_align_nw_memory_init();

    conf = load_compress_args();
    args = opt_config_parse(conf, argc, argv);

    compress_flags.ungapped_window_size = 10;
    compress_flags.match_kmer_size = 3;
    compress_flags.ext_seq_id_threshold = 50;
    compress_flags.gapped_window_size = 25;

    for (i = 0; tests[i].ref != NULL; i++) {
        ext = extend_match(
            mem,
            tests[i].ref, 0, strlen(tests[i].ref),
            tests[i].org, 0, strlen(tests[i].org));

        if (ext.rlen != tests[i].rlen || ext.olen != tests[i].olen) {
            printf("TEST %d FAILED\n", i);
            printf("The match extension for:\n");
            printf("---------------------------\n");
            printf("\t%s\n", tests[i].ref);
            printf("\t%s\n", tests[i].org);
            printf("---------------------------\n");
            printf("yield lengths: (%d, %d)\n", ext.rlen, ext.olen);
            printf("but expected:  (%d, %d)\n", tests[i].rlen, tests[i].olen);
            exit(1);
        }
    }

    cbp_align_nw_memory_free(mem);
    opt_config_free(conf);
    opt_args_free(args);

    printf("ALL TESTS PASSED\n");

    return 0;
}

