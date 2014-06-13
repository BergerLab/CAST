#ifndef __CABLAST_LINK_TO_COARSE_H__
#define __CABLAST_LINK_TO_COARSE_H__

#include <stdbool.h>
#include <stdint.h>

#include "align.h"

struct cb_link_to_coarse {
    char *diff;
    uint64_t coarse_seq_id;
    uint64_t original_start;
    uint64_t original_end;
    uint16_t coarse_start;
    uint16_t coarse_end;
    struct cb_link_to_coarse *next;
};

struct cb_link_to_coarse_data {
    uint64_t coarse_seq_id;
    uint64_t original_start;
    uint64_t original_end;
    uint16_t coarse_start;
    uint16_t coarse_end;
    uint16_t script_length;
};

struct cb_link_to_coarse *
cb_link_to_coarse_init(uint64_t coarse_seq_id,
                       uint64_t original_start, uint64_t original_end,
                       uint16_t coarse_start, uint16_t coarse_end,
                       struct cb_alignment alignment, bool dir);

struct cb_link_to_coarse *
cb_link_to_coarse_init_nodiff(uint64_t coarse_seq_id,
                              uint64_t original_start, uint64_t original_end,
                              uint16_t coarse_start, uint16_t coarse_end,
                              bool dir);

void cb_link_to_coarse_free(struct cb_link_to_coarse *link);

struct cb_link_to_coarse_data *
cb_link_to_coarse_get_data(struct cb_link_to_coarse *link);

struct cb_link_to_coarse *
cb_link_to_coarse_from_data(struct cb_link_to_coarse_data *data, char *diff){

#endif
