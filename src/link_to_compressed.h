#ifndef __CABLAST_LINK_TO_COMPRESSED_H__
#define __CABLAST_LINK_TO_COMPRESSED_H__

#include <stdbool.h>
#include <stdint.h>

struct cb_link_to_compressed_data {
    bool dir;
    int32_t org_seq_id;
    int16_t coarse_start;
    int16_t coarse_end;
    uint64_t original_start;
    uint64_t original_end;
};

struct cb_link_to_compressed {
    bool dir;
    int32_t org_seq_id;
    int16_t coarse_start;
    int16_t coarse_end;
    uint64_t original_start;
    uint64_t original_end;
    struct cb_link_to_compressed *next;
};

struct cb_link_to_compressed *
cb_link_to_compressed_init(int32_t org_seq_id, int16_t coarse_start,
                           int16_t coarse_end, uint64_t original_start,
                           uint64_t original_end, bool dir);

void cb_link_to_compressed_free(struct cb_link_to_compressed *link);

struct cb_link_to_compressed_data *
cb_link_to_compressed_get_data(struct cb_link_to_compressed *link);

#endif
