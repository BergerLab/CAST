#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "ds.h"
#include "stdbool.h"
#include "stdlib.h"

#include "coarse.h"
#include "compressed.h"
#include "decompression.h"
#include "fasta.h"
#include "seq.h"
#include "util.h"

/*Takes in an entry in a compressed database and a coarse database and returns
  the decompressed sequence that the database entry came from as a pointer to
  a struct cb_seq.*/
struct cb_seq *cb_decompress_seq(struct cb_compressed_seq *cseq,
                                   struct cb_coarse *coarsedb){
    struct DSVector *decompressed_chunks = ds_vector_create();
    struct cb_link_to_coarse *link = NULL;
    struct cb_seq *seq = NULL;
    uint64_t last_end = 0;
    int32_t overlap, decompressed_length = 0, i = 0, j = 0, copied = 0;
    char *residues = NULL;

    for (link = cseq->links; link != NULL; link = link->next) {
        struct fasta_seq *chunk = cb_coarse_read_fasta_seq(coarsedb,
                                                           link->coarse_seq_id);
        int length, coarse_len = link->coarse_end - link->coarse_start;
        char *dec_chunk;

        char *coarse_sub = malloc((coarse_len+1)*sizeof(*coarse_sub));
        assert(coarse_sub);

        memcpy(coarse_sub, chunk->seq + link->coarse_start, coarse_len);
        coarse_sub[coarse_len] = '\0';

        for (length = 0; chunk->seq[length] != '\0'; length++);

        /*overlap represents the length of the overlap of the parts of the
          decompressed sequence that has been printed and the parts of the
          decompressed sequence currently being decompressed.*/
        overlap = last_end - link->original_start + 1;
        dec_chunk = read_edit_script(link->diff, coarse_sub, length);

        /*Print all characters of the decompressed chunk past the index
          "overlap" unless overlap is greater than the length of the
          decompressed chunk.*/
        dec_chunk += overlap;
        if ((unsigned int)overlap < link->original_end - link->original_start) {
            int chunk_length = 0;
            char *section = NULL;

            for (chunk_length=0; dec_chunk[chunk_length]!='\0'; chunk_length++);

            section = malloc((chunk_length + 1)*sizeof(*section));
            assert(section);

            for (i = 0; i <= chunk_length; i++)
                section[i] = dec_chunk[i];
            ds_vector_append(decompressed_chunks, (void *)section);
            decompressed_length += chunk_length;
            if (link->original_end > last_end)
                last_end = link->original_end;
        }

        dec_chunk -= overlap;
        free(dec_chunk);
        free(coarse_sub);
        fasta_free_seq(chunk);
    }

    residues = malloc((decompressed_length+1)*sizeof(*residues));
    assert(residues);

    for (i = 0; i < decompressed_chunks->size; i++) {
        char *current_chunk = (char *)ds_vector_get(decompressed_chunks, i);
        for (j = 0; current_chunk[j] != '\0'; j++)
            residues[copied++] = current_chunk[j];
    }
    residues[copied] = '\0';
    seq = cb_seq_init(cseq->id, cseq->name, residues);
    ds_vector_free(decompressed_chunks);
    free(residues);
    return seq;
}

int get_min(int a, int b){return a<b?a:b;}
int get_max(int a, int b){return a>b?a:b;}

/*Takes in a coarse database, a compressed database, the accession number,
 *hit_from, and hit_to values of a BLAST Hsp, and a hit pad length and returns
 *a vector containing a cb_hit_expansion struct with an expanded original
 *sequence section for each link_to_compressed from the coarse sequence that is
 *in the range between the indices hit_from and hit_to.
 */
struct DSVector *
cb_coarse_expand(struct cb_coarse *coarsedb, struct cb_compressed *comdb,
                  int32_t id, int32_t hit_from, int32_t hit_to,
                  int32_t hit_pad_length){
    FILE *links = coarsedb->file_links,
         *coarse_links_index = coarsedb->file_links_index,
         *fasta = coarsedb->file_fasta,
         *compressed = comdb->file_compressed;
    struct DSVector *oseqs = ds_vector_create(), *coarse_seq_links;
    struct fasta_seq *residues = cb_coarse_read_fasta_seq(coarsedb, id);
    int64_t *seq_lengths = cb_compressed_get_lengths(comdb);
    int32_t fasta_length = strlen(residues->seq), i = 0, j = 0;

    /*Get all links_to_compressed for the coarse sequence we are expanding.*/
    coarse_seq_links =
        get_coarse_sequence_links_at(links, coarse_links_index, id);

    /*Get the residues of the coarse sequence we are expanding.*/
    for (i = 0; i < coarse_seq_links->size; i++) {
        struct cb_link_to_compressed *link =
            (struct cb_link_to_compressed *)ds_vector_get(coarse_seq_links, i);

        /*Only expand the link if it overlaps the range for the BLAST Hsp we
          are expanding from.*/
        if (link->coarse_start <= hit_to && link->coarse_end >= hit_from) {
            struct cb_link_to_coarse *current = NULL;
            struct cb_compressed_seq *seq;
            struct cb_hit_expansion *expansion;
            uint64_t original_start, original_end, original_range;
            char *orig_str;
            bool dir = link->dir;

            /*Calculate the range in the original sequence for the section of
              the original sequence we want to re-create with this expansion.*/
            original_start =
                get_max(0, (dir ?
                            get_min(hit_from + (link->original_start -
                                                link->coarse_start),
                                    hit_from + (link->original_end -
                                                link->coarse_end)) :
                            get_min(link->original_start +
                                    link->coarse_end - hit_to,
                                    link->original_end -
                                    (hit_to-link->coarse_start)))
                             - hit_pad_length);
            original_end =
                get_min((dir ? get_max(hit_to + (link->original_start -
                                                 link->coarse_start),
                                       hit_to + (link->original_end -
                                                 link->coarse_end)) :
                               get_max(link->original_end -
                                       (hit_from-link->coarse_start),
                                       link->original_start +
                                       link->coarse_end-hit_from))
                         + hit_pad_length, seq_lengths[link->org_seq_id] - 1);
            original_range = original_end - original_start + 1;

            seq = cb_compressed_read_seq_at(comdb, link->org_seq_id);

            orig_str = malloc((original_end-original_start+2) *
                              sizeof(*orig_str));
            assert(orig_str);

            for (j = 0; j < original_end-original_start+1; orig_str[j++]='?');

            /*Run decode_edit_script for each link_to_coarse in the compressed
              sequence to re-create the section of the original string.*/ 
            current = seq->links;
            for (; current; current = current->next) {
                int coarse_range = current->coarse_end - current->coarse_start,
                    init_i0 = current->original_start-(int32_t)original_start,
                    last_i0 = init_i0 + coarse_range;

                if (0 < last_i0 && (int32_t)original_range > init_i0)
                    decode_edit_script(orig_str, original_end-original_start+1,
                                        original_start, coarsedb, current);
            }
            orig_str[original_end-original_start+1] = '\0';

/*printf("%s\n", orig_str);*/
            expansion = malloc(sizeof(*expansion));
            assert(expansion);

            expansion->offset = (int64_t)original_start;
            expansion->seq = cb_seq_init(link->org_seq_id,
                                         seq->name, orig_str);
            ds_vector_append(oseqs, (void *)expansion);
            free(orig_str);
            cb_compressed_seq_free(seq);
        }
    }

    ds_vector_free(coarse_seq_links);
    fasta_free_seq(residues);
    return oseqs;
}
