#ifndef __CABLAST_DECOMPRESSION_H__
#define __CABLAST_DECOMPRESSION_H__

#include <stdint.h>
#include <stdio.h>

#include "ds.h"
#include "stdbool.h"

#include "coarse.h"
#include "compressed.h"
#include "seq.h"

struct cb_seq *cb_decompress_seq(struct cb_compressed_seq *cseq,
                                   struct cb_coarse *coarsedb);
struct DSVector *
cb_coarse_expand(struct cb_coarse *coarsedb, struct cb_compressed *comdb,
                  int32_t id, int32_t start, int32_t end,
                  int32_t hit_pad_length);

#endif
