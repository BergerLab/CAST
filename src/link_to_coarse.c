#include <assert.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

#include "align.h"
#include "edit_scripts.h"
#include "link_to_coarse.h"

struct cb_link_to_coarse *
cb_link_to_coarse_init(int32_t coarse_seq_id,
                        uint64_t original_start, uint64_t original_end,
                        uint16_t coarse_start, uint16_t coarse_end,
                        struct cb_alignment alignment, bool dir){
    struct cb_link_to_coarse *link;

    link = malloc(sizeof(*link));
    assert(link);

    link->coarse_seq_id  = coarse_seq_id;
    link->original_start = original_start;
    link->original_end   = original_end;
    link->coarse_start   = coarse_start;
    link->coarse_end     = coarse_end;
    link->next           = NULL;
    link->diff           = make_edit_script(alignment.org, alignment.ref, dir,
                                            alignment.length);
    assert(link->diff);

    return link;
}

struct cb_link_to_coarse *
cb_link_to_coarse_init_nodiff(int32_t coarse_seq_id,
                               uint64_t original_start, uint64_t original_end,
                               uint16_t coarse_start, uint16_t coarse_end,
                               bool dir){
    struct cb_link_to_coarse *link;

    link = malloc(sizeof(*link));
    assert(link);

    link->diff = malloc(2*sizeof(*(link->diff)));
    assert(link->diff);

    link->diff[0]        = dir ? '0' : '1';
    link->diff[1]        = '\0';
    link->coarse_seq_id  = coarse_seq_id;
    link->original_start = original_start;
    link->original_end   = original_end;
    link->coarse_start   = coarse_start;
    link->coarse_end     = coarse_end;
    link->next           = NULL;

    return link;
}

void
cb_link_to_coarse_free(struct cb_link_to_coarse *link)
{
    free(link->diff);
    free(link);
}

