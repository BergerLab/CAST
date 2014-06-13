#include <assert.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "align.h"
#include "edit_scripts.h"
#include "link_to_coarse.h"

/*Initializes a link to a sequence in the coarse database.  Creates the
  link's edit script using the alignment passed in.*/
struct cb_link_to_coarse *
cb_link_to_coarse_init(uint64_t coarse_seq_id,
                       uint64_t original_start, uint64_t original_end,
                       uint16_t coarse_start, uint16_t coarse_end,
                       struct cb_alignment alignment, bool dir){
    struct cb_link_to_coarse *link = malloc(sizeof(*link));
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

/*Initializes a link to a coarse sequence in which the coarse and original
  sequences are identical and therefore a blank edit script can be used.*/
struct cb_link_to_coarse *
cb_link_to_coarse_init_nodiff(uint64_t coarse_seq_id,
                              uint64_t original_start, uint64_t original_end,
                              uint16_t coarse_start, uint16_t coarse_end,
                              bool dir){
    struct cb_link_to_coarse *link = malloc(sizeof(*link));
    assert(link);

    link->coarse_seq_id  = coarse_seq_id;
    link->original_start = original_start;
    link->original_end   = original_end;
    link->coarse_start   = coarse_start;
    link->coarse_end     = coarse_end;

    /*Give the link a blank edit script*/
    link->diff = malloc(2*sizeof(*(link->diff)));
    assert(link->diff);

    link->diff[0]        = dir ? '0' : '1';
    link->diff[1]        = '\0';

    link->next           = NULL;

    return link;
}

/*Freeing function for a cb_link_to_coarse*/
void cb_link_to_coarse_free(struct cb_link_to_coarse *link){
    if (link != NULL) {
        cb_link_to_coarse_free(link->next);
        free(link->diff);
        free(link);
    }
}

/*Takes in a link_to_coarse and gets its coarse sequence ID, its indices, and
 *the length of its edit script, which are returned in a
 *cb_link_to_coarse_get_data struct.
 */
struct cb_link_to_coarse_data *
cb_link_to_coarse_get_data(struct cb_link_to_coarse *link){
    struct cb_link_to_coarse_data *data = malloc(sizeof(*data));
    assert(data);

    data->coarse_seq_id  = link->coarse_seq_id;
    data->original_start = link->original_start;
    data->original_end   = link->original_end;
    data->coarse_start   = link->coarse_start;
    data->coarse_end     = link->coarse_end;

    data->script_length = (uint16_t)0;
    while (link->diff[data->script_length] != '\0')
        data->script_length++;

    return data;
}

/*Takes in a cb_link_to_coarse_get_data struct and an edit script and uses them
  to create a new cb_link_to_coarse struct.*/
struct cb_link_to_coarse *
cb_link_to_coarse_from_data(struct cb_link_to_coarse_data *data, char *diff){
    struct cb_link_to_coarse *link = malloc(sizeof(*link));
    assert(link);

    link->coarse_seq_id  = data->coarse_seq_id;
    link->original_start = data->original_start;
    link->original_end   = data->original_end;
    link->coarse_start   = data->coarse_start;
    link->coarse_end     = data->coarse_end;
    link->next           = NULL;
    link->diff           = diff;

    return link;
}

