#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "ds.h"
#include "stdbool.h"
#include "stdlib.h"

#include "coarse.h"
#include "compressed.h"
#include "decompression.h"
#include "DNAutils.h"
#include "edit_scripts.h"
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
struct DSVector *cb_coarse_expand(struct cb_coarse_r *coarse_db,
                                  struct cb_compressed *comdb,
                                  int32_t id, int32_t hit_from, int32_t hit_to,
                                  int32_t hit_pad_length){
    struct cb_coarse *coarsedb = coarse_db->db;
    FILE *links              = coarsedb->file_links,
         *coarse_links_index = coarsedb->file_links_index,
         *fasta              = coarsedb->file_fasta,
         *compressed         = comdb->file_compressed;
    struct DSVector *oseqs = ds_vector_create(), *coarse_seq_links = NULL;
    int64_t *seq_lengths      = cb_compressed_get_lengths(comdb),
            *seq_base_indices = coarse_db->seq_base_indices,
            hit_from_ind      = hit_from + seq_base_indices[id],
            hit_to_ind        = hit_to + seq_base_indices[id];
    int32_t first_block = hit_from_ind / coarse_db->link_block_size,
            last_block  = hit_to_ind / coarse_db->link_block_size,
            i = 0, j = 0, k = 0;

    for (i = first_block; i <= last_block; i++) {
        /*Current block of links*/
        struct DSVector *link_block = cb_coarse_r_get_block(coarse_db, i);
        /*Indices in the current block*/
        struct DSVector *ind_block =
          (struct DSVector *)ds_vector_get(coarse_db->link_inds_by_block, i);

        for (j = 0; j < link_block->size; j++){
            struct cb_link_to_compressed_data *link =
              (struct cb_link_to_compressed_data *)ds_vector_get(link_block, j);
            int64_t link_ind     = *(int64_t *)ds_vector_get(ind_block, j),
                    coarse_start = link->coarse_start,
                    coarse_end   = link->coarse_end;

            /*Determine which sequence the link belongs to*/
            for (k = 0; k <= coarse_db->num_coarse_seqs; k++)
                if (coarse_db->seq_link_counts[k] > link_ind)
                    break;

            coarse_start += seq_base_indices[k];
            coarse_end += seq_base_indices[k];

            /*Only expand the link if it overlaps the range for the BLAST Hsp we
              are expanding from.*/
            if (coarse_start <= hit_to_ind && coarse_end >= hit_from_ind) {
                struct cb_link_to_coarse *current = NULL;
                struct cb_compressed_seq *seq;
                struct cb_hit_expansion *expansion;
                uint64_t original_start, original_end, original_range;
                char *orig_str;
                bool dir = link->dir;

                /*Calculate the range in the original sequence for the section
                 *of the original sequence we want to re-create with this
                 *expansion.
                 */
                original_start =
                  get_max(0, (dir ? get_min(hit_from + (link->original_start -
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

                seq = comdb->seqs->size == 0 ?
                        cb_compressed_read_seq_at(comdb, link->org_seq_id) :
                        cb_compressed_seq_at(comdb, link->org_seq_id);

                orig_str = malloc((original_range+1)*sizeof(*orig_str));
                assert(orig_str);

                for (k = 0; k < original_range; orig_str[k++]='?');

                /*Run decode_edit_script for each link_to_coarse in the compressed
                  sequence to re-create the section of the original string.*/
                current = seq->links;
                for (; current; current = current->next) {
                    int coarse_range = current->coarse_end - current->coarse_start,
                        init_i0 = current->original_start-(int32_t)original_start,
                        last_i0 = init_i0 + coarse_range;

                    if (0 < last_i0 && (int32_t)original_range > init_i0)
                        decode_edit_script(orig_str, original_range,
                                           original_start, coarse_db, current);
                }
                orig_str[original_range] = '\0';
printf("%s\n", orig_str);
                expansion = malloc(sizeof(*expansion));
                assert(expansion);

                expansion->offset = (int64_t)original_start;
                expansion->seq = cb_seq_init(link->org_seq_id, seq->name, orig_str);
                ds_vector_append(oseqs, (void *)expansion);

                free(orig_str);

                /*Free the compressed sequence if it was loaded using
                  cb_compressed_read_seq_at*/
                if (comdb->seqs->size == 0)
                    cb_compressed_seq_free(seq);
            }
        }

        ds_vector_free_no_data(link_block);
    }

    /*if (!coarse_db->links)
        ds_vector_free(coarse_seq_links);*/

    return oseqs;
}

/*An implementation of decode_edit_script that is used in Po-Ru's C++ version
 *used in search for expanding coarse BLAST hits.
 *
 *Takes in a string orig, which is the section of the original sequence we want
 *to re-create, the length of this section of the string, the index in the
 *original sequence of the start of the section, the coarse database, and a
 *link_to_coarse we want to apply and applies the link to re-create part of or
 *all of the section of the original sequence.
 */
void decode_edit_script(char *orig, int dest_len, int original_start,
                        struct cb_coarse_r *coarsedb,
                        struct cb_link_to_coarse *link){
    struct fasta_seq *fasta = NULL;
    struct edit_info *edit = NULL;
    int i = 0, i0 = 0, i1 = 0, coarse_pos, last_edit_str_len, script_pos;
    char *diff = link->diff, *residues;
    bool fwd = (diff[0] & ((char)0x7f)) == '0';

    if (coarsedb->all_residues == NULL) {
        fasta = cb_coarse_read_fasta_seq_r(coarsedb, link->coarse_seq_id); 
        residues = fasta->seq;
    }
    else
        residues = coarsedb->all_residues +
                   coarsedb->seq_base_indices[link->coarse_seq_id];

    /*The link represents an exact match so there are no edits to make*/
    if (diff[1] == '\0' && fwd) {
        int starting_i0 = -1, last_i0 = -1;

        i0 = link->original_start - original_start;
        for (i1 = link->coarse_start; i1 < link->coarse_end; i0++, i1++)
            if (0 <= i0 && i0 < dest_len) {
                starting_i0 = (starting_i0 == -1 ? i0 : starting_i0);
                last_i0 = i0;
                orig[i0] = residues[i1];
            }
        /*If the link is from a reverse-complement match, convert the original
          string to its reverse complement.*/
        if (fasta != NULL)
            fasta_free_seq(fasta);

        return;
    }

    coarse_pos = link->coarse_start;
    last_edit_str_len = 0;

    edit = malloc(sizeof(*edit));
    assert(edit);

    script_pos = 1;

    /*We are decompressing a link from a forward match*/
    if (fwd) {
        i0 = link->original_start - original_start;
        while (next_edit(diff, &script_pos, edit)) {
            int x = 0, xmin = -i0, xmax = dest_len - i0;

            for (x = maximum(0, xmin);
                 x < minimum(edit->last_dist-last_edit_str_len, xmax); x++)
                orig[i0+x] = residues[x+coarse_pos];

            i0 += edit->last_dist - last_edit_str_len;
            coarse_pos += edit->last_dist - last_edit_str_len;
            for (i = 0; i < edit->str_length; i++)
                if (edit->str[i] != '-') {
                    if (0 <= i0 && i0 < dest_len)
                        orig[i0] = edit->str[i];
                    i0++;
                }
            if (edit->is_subdel) coarse_pos += edit->str_length;

            last_edit_str_len = edit->str_length;

            if (i0 >= dest_len) {
                free(edit);
                if (fasta != NULL)
                    fasta_free_seq(fasta);

                return;
            }
        }
    }
    else {
        i0 = link->original_end - original_start;
        while (next_edit(diff, &script_pos, edit)) {
            int x = 0, xmin = i0 - dest_len + 1, xmax = i0 + 1;

            for (x = maximum(0, xmin);
                 x < minimum(edit->last_dist-last_edit_str_len, xmax); x++)
                orig[i0-x] = base_complement(residues[x+coarse_pos]);

            i0 -= edit->last_dist - last_edit_str_len;
            coarse_pos += edit->last_dist - last_edit_str_len;
            for (i = 0; i < edit->str_length; i++)
                if (edit->str[i] != '-') {
                    if (0 <= i0 && i0 < dest_len)
                        orig[i0] = base_complement(edit->str[i]);
                    i0--;
                }

            if (edit->is_subdel) coarse_pos += edit->str_length;

            last_edit_str_len = edit->str_length;

            if (i0 < 0) {
                free(edit->str);
                free(edit);
                if (fasta != NULL)
                    fasta_free_seq(fasta);

                return;
            }
        }
    }
    if ((fwd && i0 < dest_len) || (!fwd && i0 >= 0)) {
        int dir = fwd ? 1 : -1;
        for (i1 = coarse_pos; i1 <= link->coarse_end; i1++) {
            if (0 <= i0 && i0 < dest_len)
                orig[i0] = fwd ? residues[i1] : base_complement(residues[i1]);
            i0 += dir;
        }
    }
    free(edit);
    if (fasta != NULL)
        fasta_free_seq(fasta);
}

