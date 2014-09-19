#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "seq.h"

/*Calls cb_seq_init_range to initialize a sequence with the data passed in
  that copies all of the residues passed in.*/
struct cb_seq *cb_seq_init(int32_t id, char *name, char *residues){
    return cb_seq_init_range(id, name, residues, 0, strlen(residues));
}

/*@param id, name, residues: The ID number, name, and residues of a sequence
 *@param start, end: The starting and ending indices of the part of the sequence
 *to copy
 *
 *@return: A new sequence data structure for the sequence copied
 */
struct cb_seq *cb_seq_init_range(int32_t id, char *name, char *residues,
                                 int32_t start, int32_t end){
    int len;

    struct cb_seq *seq = malloc(sizeof(*seq));
    assert(seq);

    seq->id       = id;
    seq->length   = end - start;
    seq->name     = NULL;
    seq->residues = NULL;

    assert(seq->length > 0);

    len = strlen(name);
    if (len > 0) {
        seq->name = malloc((1+len)*sizeof(*seq->name));
        assert(seq->name);

        strcpy(seq->name, name);
    }
    if (seq->length > 0) {
        assert(start >= 0 && start < end);

        seq->residues = malloc((1+seq->length)*sizeof(*(seq->residues)));
        assert(seq->residues);

        strncpy(seq->residues, residues + start, seq->length);
        seq->residues[seq->length] = '\0';
        assert(strlen(seq->residues) > 0);
    }

    assert(seq->residues != NULL);

    return seq;
}

//Frees a cb_seq struct
void cb_seq_free(struct cb_seq *seq){
    free(seq->name);
    free(seq->residues);
    free(seq);
}

struct cb_hit_expansion *
cb_hit_expansion_init(int64_t offset, bool strand,
                      int32_t id, char *name, char *residues){
    struct cb_hit_expansion *expansion = malloc(sizeof(*expansion));
    assert(expansion);

    expansion->offset = offset;
    expansion->strand = strand ? '+' : '-';
    expansion->seq    = cb_seq_init(id, name, residues);

    return expansion;
}


//Frees a hit expansion
void cb_hit_expansion_free(struct cb_hit_expansion *expansion){
    cb_seq_free(expansion->seq);
    free(expansion);
}
