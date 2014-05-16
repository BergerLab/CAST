#ifndef __CABLAST_SEQ_H__
#define __CABLAST_SEQ_H__

#include <stdint.h>

struct cb_seq {
    int32_t id;
    char *name;
    char *residues;
    int32_t length;
};

struct cb_hit_expansion {
    struct cb_seq *seq;
    int64_t offset;
};

struct cb_seq *
cb_seq_init(int32_t id, char *name, char *residues);

struct cb_seq *
cb_seq_init_range(int32_t id, char *name, char *residues,
                   int32_t start, int32_t end);

void
cb_seq_free(struct cb_seq *seq);

void
cb_hit_expansion_free(struct cb_hit_expansion *expansion);

#endif
