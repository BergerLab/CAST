#ifndef __CABLAST_COMPRESSED_H__
#define __CABLAST_COMPRESSED_H__

#include <stdint.h>
#include <stdio.h>

#include "ds.h"

#include "align.h"
#include "bitpack.h"
#include "edit_scripts.h"
#include "link_to_coarse.h"
#include "seq.h"

#include "stdbool.h"

struct cb_link_to_coarse *
cb_link_to_coarse_init(int32_t coarse_seq_id,
                        uint64_t original_start, uint64_t original_end,
                        uint16_t coarse_start, uint16_t coarse_end,
                        struct cb_alignment alignment, bool dir);

struct cb_link_to_coarse *
cb_link_to_coarse_init_nodiff(int32_t coarse_seq_id,
                               uint64_t original_start, uint64_t original_end,
                               uint16_t coarse_start, uint16_t coarse_end,
                               bool dir);

void
cb_link_to_coarse_free(struct cb_link_to_coarse *link);

struct cb_compressed_seq {
    uint64_t id;
    char *name;
    struct cb_link_to_coarse *links;
};

struct cb_compressed_seq *
cb_compressed_seq_init(int32_t id, char *name);

void
cb_compressed_seq_free(struct cb_compressed_seq *seq);

void
cb_compressed_seq_addlink(struct cb_compressed_seq *seq,
                           struct cb_link_to_coarse *link);

struct cb_coarse;

struct cb_compressed {
    struct DSVector *seqs;
    FILE *file_compressed;
    FILE *file_index;
};

struct cb_compressed *
cb_compressed_init(FILE *file_compressed, FILE *file_index);

void
cb_compressed_free(struct cb_compressed *com_db);

int32_t
cb_compressed_size(struct cb_compressed *com_db);

void
cb_compressed_add(struct cb_compressed *com_db,
                   struct cb_compressed_seq *seq);

void
cb_compressed_save_binary(struct cb_compressed *com_db);

void
cb_compressed_save_plain(struct cb_compressed *com_db);

void
cb_compressed_write(struct cb_compressed *com_db,
                     struct cb_compressed_seq *seq);

void
cb_compressed_write_binary(struct cb_compressed *com_db,
                            struct cb_compressed_seq *seq);

struct cb_compressed_seq *
cb_compressed_seq_at(struct cb_compressed *com_db, int32_t i);

char *get_compressed_header(FILE *f);
struct cb_link_to_coarse *read_link(FILE *f);
struct cb_compressed_seq *get_compressed_seq(FILE *f, int id);
struct cb_compressed_seq **read_compressed(FILE *f);
int64_t cb_compressed_link_offset(struct cb_compressed *comdb, int id);
struct cb_compressed_seq *cb_compressed_read_seq_at(
                                                 struct cb_compressed *comdb,
                                                 int32_t id);
int64_t cb_compressed_get_seq_length(FILE *f);
int64_t *cb_compressed_get_lengths(struct cb_compressed *comdb);
#endif
