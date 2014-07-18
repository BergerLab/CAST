#ifndef __CABLAST_DECOMPRESSION_H__
#define __CABLAST_DECOMPRESSION_H__

#include <stdint.h>
#include <stdio.h>

#include "clibs/include/ds.h"
#include "stdbool.h"

#include "coarse.h"
#include "compressed.h"
#include "database.h"
#include "link_to_coarse.h"
#include "seq.h"

struct DSVector *cb_coarse_expand(struct cb_coarse_r *coarse_db,
                                  struct cb_compressed *comdb,
                                  int32_t id, int32_t hit_from, int32_t hit_to,
                                  int32_t hit_pad_length);

void decode_edit_script(char *orig, int dest_len, int original_start,
                        struct cb_coarse_r *coarsedb,
                        struct cb_link_to_coarse *link);
#endif
