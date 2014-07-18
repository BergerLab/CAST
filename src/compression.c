#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "clibs/include/ds.h"
#include "compression.h"
#include "DNAutils.h"
#include "edit_scripts.h"
#include "flags.h"
#include "seeds.h"

struct worker_args {
    struct cb_database *db;
    struct DSQueue *jobs;
    int id;
};

struct extend_match {
    char *rseq;
    char *oseq;
    int32_t rlen;
    int32_t olen;
};

struct extend_match
extend_match(struct cb_align_nw_memory *mem,
             char *rseq, int32_t rstart, int32_t rend, int32_t resind,
             int32_t dir1, char *oseq, int32_t ostart, int32_t oend,
             int32_t current, int32_t dir2);

static int32_t add_without_match(struct cb_coarse *coarse_db,
                                 struct cb_seq *org_seq,
                                 int32_t ostart, int32_t oend,
                                 struct cb_seeds_add_memory *seeds_mem);

static void *cb_compress_worker(void *data);

static int32_t min(int32_t a, int32_t b);

/*Takes in the database and the number of threads being used and creates a new
 *compress worker for each processor being used and starts a new pthread for
 *each one to start compression.
 */
struct cb_compress_workers *cb_compress_start_workers(struct cb_database *db,
                                                      int32_t num_workers){
    struct DSQueue *jobs;
    struct cb_compress_workers *workers;
    struct worker_args **wargs;
    int32_t i, errno;

    jobs = ds_queue_create(20);

    wargs = malloc(num_workers*sizeof(*wargs));
    assert(wargs);

    for (i = 0; i < num_workers; i++) {
        wargs[i] = malloc(sizeof(*(wargs[i])));
        assert(wargs[i]);

        wargs[i]->db   = db;
        wargs[i]->jobs = jobs;
        wargs[i]->id   = i;
    }

    workers = malloc(sizeof(*workers));
    assert(workers);

    workers->threads     = malloc(num_workers*sizeof(*workers->threads));
    assert(workers->threads);

    workers->num_workers = num_workers;
    workers->jobs        = jobs;
    workers->args        = (void *)wargs;

    for (i = 0; i < num_workers; i++) {
        errno = pthread_create(&workers->threads[i], NULL,
                               cb_compress_worker, (void *)wargs[i]);
        if (errno != 0) {
            fprintf(stderr,
              "cb_compress_start_workers: Could not start thread. Errno: %d",
              errno);
            exit(1);
        }
    }

    return workers;
}

//Closes the jobs queue and joins the compression workers.
void cb_compress_join_workers(struct cb_compress_workers *workers){
    int errno;

    ds_queue_close(workers->jobs);

    for (int i = 0; i < workers->num_workers; i++) {
        errno = pthread_join(workers->threads[i], NULL);
        if (errno != 0) {
            fprintf(stderr,
              "cb_compress_join_workers: Could not join thread. Errno: %d",
              errno);
            exit(1);
        }
    }
}

//Frees the compression workers.
void cb_compress_free_workers(struct cb_compress_workers *workers){
    ds_queue_free(workers->jobs);
    free(workers->args);
    free(workers->threads);
    free(workers);
}

/*Takes in the compression workers struct and an original sequence and puts the
  original sequence on the jobs queue.*/
void cb_compress_send_job(struct cb_compress_workers *workers,
                          struct cb_seq *org_seq){
    ds_queue_put(workers->jobs, (void *)org_seq);
}

/*The compress worker function run on a pthread for running compression.  While
 *the jobs queue is open and contains sequences, the compress worker gets
 *sequences from the queue, compresses them, and writes them to the compressed
 *database file in the order they are in in the original FASTA file.
 */
static void *cb_compress_worker(void *data){
    struct worker_args *args;
    struct cb_align_nw_memory *mem;
    struct cb_seeds_add_memory *seeds_mem;
    struct cb_seq *s;
    struct cb_compressed_seq *cseq = NULL;
    struct DSLinkedList *seqs_to_write = ds_list_create();

    args = (struct worker_args *)data;
    mem = cb_align_nw_memory_init();
    seeds_mem = cb_seeds_add_memory_init();

    while (NULL != (s = (struct cb_seq *)ds_queue_get(args->jobs)) ||
           NULL != seqs_to_write->first) {
        if (NULL != s) {
            cseq = cb_compress(args->db->coarse_db, s, mem, seeds_mem);
            cb_coarse_db_update_dbsize(args->db->coarse_db, s->length);
            ds_list_append(seqs_to_write, (struct cb_compressed_seq *)cseq);
            cb_seq_free(s);
        }

        struct DSListNode *current_node = seqs_to_write->first;
        struct cb_compressed_seq *current_seq = current_node ?
          (struct cb_compressed_seq *)seqs_to_write->first->data : NULL;
        if (current_seq != NULL &&
              current_seq->id == args->db->com_db->next_seq_to_write) {
            cb_compressed_write_binary(args->db->com_db, current_seq);
            cb_compressed_seq_free(current_seq);
            ds_list_remove(seqs_to_write, seqs_to_write->first);

            //Atmoically update which compressed sequence is the next to write.
            __sync_fetch_and_add(&(args->db->com_db->next_seq_to_write), 1);
        }
    }

    cb_align_nw_memory_free(mem);

    return NULL;
}

/*@param coarse_db: The coarse database
 *@param org_seq: The original sequence to compress
 *@param mem: The compress worker's alignment memory
 *@param seeds_mem: The compress worker's seeds_add memory
 *
 *@return: A compressed sequence struct with links to the the coarse database
 *  representing parts of the original sequence that aligned well with a
 *  sequence in the coarse database.
 *
 *Goes through a sequence from the FASTA file passed in and finds locations in
 *the sequence that align well with sequences in the coarse database, adding
 *links to them in the compressed sequence that is returned.  New coarse
 *sequences are added to the coarse database if more than max_chunk_size
 *(default 10000) have been traversed without finding a match or when there are
 *bases before a match that do not align well with a coarse sequence.  For more
 *detail read the comments in the function.
 */
struct cb_compressed_seq *
cb_compress(struct cb_coarse *coarse_db, struct cb_seq *org_seq,
            struct cb_align_nw_memory *mem,
            struct cb_seeds_add_memory *seeds_mem){
    struct cb_seeds *coarse_seeds = coarse_db->seeds;
    struct extend_match mseqs_fwd, mseqs_rev;
    struct cb_coarse_seq *coarse_seq;
    struct cb_compressed_seq *cseq =
             cb_compressed_seq_init(org_seq->id, org_seq->name);
    struct cb_seed_loc *seedLoc;
    struct cb_alignment alignment;
    int32_t seed_size = coarse_db->seeds->seed_size,
            resind = -1, current = 0, i = 0, new_coarse_seq_id = -1,
            fwd_rlen, rev_rlen, fwd_olen, rev_olen, index,
            chunks = 0,
            max_chunk_size = compress_flags.max_chunk_size,
            max_section_size = max_chunk_size * 2,
            overlap = compress_flags.overlap,
            start_of_section = 0,
            org_seq_len = org_seq->length,
            end_of_chunk = min(start_of_section + max_chunk_size,
                                   org_seq_len),
            end_of_section = min(start_of_section + max_section_size,
                                   org_seq_len),
            org_seq_id = org_seq->id, *loc_counts = coarse_seeds->loc_counts;
    char *kmer, revcomp[seed_size+1], *org_seq_residues = org_seq->residues;
    bool *matches, *matches_temp, found_match;

    //Initialize the matches and matches_temp arrays.
    matches = malloc(max_section_size*sizeof(*matches));
    assert(matches);

    matches_temp = malloc(max_section_size*sizeof(*matches_temp));
    assert(matches_temp);

    for (i = 0; i < max_section_size; i++) {
        matches[i] = true;
        matches_temp[i] = true;
    }

    for (current = 0; current <= org_seq_len-seed_size; current++) {
        found_match = false;

        /*If we are at the beginning of the first chunk of the first sequence,
         *add the first chunk without a match and skip ahead to the start of
         *the second chunk.
         */
        if (current == 0 && org_seq_id == 0) {
            new_coarse_seq_id =
              add_without_match(coarse_db, org_seq, 0,
                                minimum(org_seq_len, max_chunk_size),
                                seeds_mem);

            cb_compressed_seq_addlink(cseq, cb_link_to_coarse_init_nodiff(
                                               new_coarse_seq_id, 0,
                                               end_of_chunk - 1, 0,
                                               end_of_chunk - 1, true));

            if (end_of_chunk < org_seq_len - seed_size) {
                start_of_section += max_chunk_size - overlap;
                end_of_chunk = min(start_of_section + max_chunk_size,
                                   org_seq_len);
                end_of_section = min(start_of_section + max_section_size,
                                     org_seq_len);
                current = start_of_section-1;
            }
            chunks++;
            continue;
        }

        //Get the k-mer and allocate a copy of its reverse complement
        kmer = org_seq_residues+current;
        kmer_revcomp(revcomp, kmer, seed_size);

        struct cb_seed_loc ***locs = coarse_seeds->locs;
        int32_t hash         = hash_kmer(coarse_seeds, kmer),
                hash_revcomp = hash_kmer(coarse_seeds, revcomp);

        /*The locations of all seeds in the database that start with the
          current k-mer.*/
        int32_t seeds_count = loc_counts[hash];

        /*The locations of all seeds in the database that start with the
          current k-mer's reverse complement.*/
        int32_t seeds_r_count = loc_counts[hash_revcomp];

        for (i = 0; i < seeds_count; i++) {
            if (found_match)
                break;

            seedLoc = locs[hash][i];

            resind = seedLoc->residue_index;
            coarse_seq = cb_coarse_get(coarse_db, seedLoc->coarse_seq_id);

            if (resind + seed_size > coarse_seq->seq->length)
                continue;

            if ((attempt_ext(current, -1, org_seq_residues,
                             end_of_section - start_of_section,
                             start_of_section, resind, -1,
                             coarse_seq->seq->residues,
                             coarse_seq->seq->length, 0) +
                 attempt_ext(current+seed_size-1, 1, org_seq_residues,
                             end_of_section - start_of_section,
                             start_of_section, resind+seed_size-1, 1,
                             coarse_seq->seq->residues,
                             coarse_seq->seq->length, 0)) >
                 compress_flags.min_progress) {
                mseqs_rev = extend_match(mem,
                                         coarse_seq->seq->residues, 0,
                                         coarse_seq->seq->length, resind, -1,
                                         org_seq_residues, start_of_section,
                                         end_of_section, current, -1);

                mseqs_fwd = extend_match(mem,
                                         coarse_seq->seq->residues, 0,
                                         coarse_seq->seq->length,
                                         resind + seed_size - 1, 1,
                                         org_seq_residues, start_of_section,
                                         end_of_section,
                                         current + seed_size - 1, 1);

                fwd_rlen = cb_align_length_nogaps(mseqs_fwd.rseq);
                rev_rlen = cb_align_length_nogaps(mseqs_rev.rseq);
                fwd_olen = cb_align_length_nogaps(mseqs_fwd.oseq);
                rev_olen = cb_align_length_nogaps(mseqs_rev.oseq);

                //If the match was too short, try the next seed.
                if (rev_olen+seed_size+fwd_olen-1 < compress_flags.min_match_len) {
                    free(mseqs_fwd.rseq);
                    free(mseqs_fwd.oseq);
                    free(mseqs_rev.rseq);
                    free(mseqs_rev.oseq);

                    continue;
                }

                found_match = true;

                /*Concatenate the extensions and the k-mer for the coarse
                  sequence*/
                index = 0;

                alignment.ref = malloc((strlen(mseqs_rev.rseq)+seed_size+
                                        strlen(mseqs_fwd.rseq)+1)*
                                        sizeof(*(alignment.ref)));
                assert(alignment.ref);

                for (i = 0; mseqs_rev.rseq[i] != '\0'; i++)
                    alignment.ref[index++] = mseqs_rev.rseq[i];
                for (i = 0; i < seed_size; i++)
                    alignment.ref[index++] = kmer[i];
                for (i = 0; mseqs_fwd.rseq[i] != '\0'; i++)
                    alignment.ref[index++] = mseqs_fwd.rseq[i];
                alignment.ref[index] = '\0';

                /*Concatenate the extensions and the k-mer for the original
                  sequence*/
                index = 0;

                alignment.org = malloc((strlen(mseqs_rev.oseq)+seed_size+
                                        strlen(mseqs_fwd.oseq)+1)*
                                        sizeof(*(alignment.org)));
                assert(alignment.org);

                for (i = 0; mseqs_rev.oseq[i] != '\0'; i++)
                    alignment.org[index++] = mseqs_rev.oseq[i];
                for (i = 0; i < seed_size; i++)
                    alignment.org[index++] = kmer[i];
                for (i = 0; mseqs_fwd.oseq[i] != '\0'; i++)
                    alignment.org[index++] = mseqs_fwd.oseq[i];
                alignment.org[index] = '\0';

                alignment.length = mseqs_rev.olen + seed_size + mseqs_fwd.olen;

                //Make a new chunk for the parts of the chunk before the match
                if (current - rev_olen - start_of_section > 0) {
                    new_coarse_seq_id = add_without_match(coarse_db, org_seq,
                                                    start_of_section,
                                                    current - rev_olen +
                                                      compress_flags.overlap,
                                                    seeds_mem);
                    cb_compressed_seq_addlink(cseq,
                      cb_link_to_coarse_init_nodiff(
                        new_coarse_seq_id,
                        start_of_section,
                        current - rev_olen
                                + compress_flags.overlap - 1,
                        0, current - rev_olen
                                   + compress_flags.overlap
                                   - start_of_section - 1,
                        true));
                    chunks++;
                }

                //Add a link to the coarse sequence in the compressed sequence.
                cb_compressed_seq_addlink(cseq,
                                          cb_link_to_coarse_init(coarse_seq->id,
                                            current - rev_olen,
                                            current + seed_size + fwd_olen - 1,
                                            resind - rev_rlen,
                                            resind + seed_size + fwd_rlen - 1,
                                            alignment, true));

                //Add a link to the compressed sequence in the coarse sequence.
                cb_coarse_seq_addlink(coarse_seq,
                                      cb_link_to_compressed_init(
                                        org_seq->id,
                                        resind - rev_rlen,
                                        resind + seed_size + fwd_rlen - 1,
                                        current - rev_olen,
                                        current + seed_size + fwd_olen - 1,
                                        true));

                //Update the current position in the sequence.
                if (current + fwd_olen <= org_seq_len-seed_size-1)
                    start_of_section = current + fwd_olen -
                                       compress_flags.overlap + seed_size;
                else
                    start_of_section = current + fwd_olen + seed_size;

                current = start_of_section - 1;
                end_of_chunk = min(start_of_section + max_chunk_size,
                                     org_seq_len);
                end_of_section = min(start_of_section + max_section_size,
                                     org_seq_len);

                chunks++;

                free(mseqs_fwd.rseq);
                free(mseqs_fwd.oseq);
                free(mseqs_rev.rseq);
                free(mseqs_rev.oseq);

                free(alignment.ref);
                free(alignment.org);
            }
        }
        for (i = 0; i < seeds_r_count; i++) {
            /*If we found a match in the seed locations for the k-mer, then
             *there is no need to check the locations for the reverse
             *complement.
             */
            seedLoc = locs[hash_revcomp][i];
            if (found_match)
                break;

            resind = seedLoc->residue_index;
            coarse_seq = cb_coarse_get(coarse_db, seedLoc->coarse_seq_id);

            if (resind + seed_size > coarse_seq->seq->length)
                continue;

            if ((attempt_ext(current, -1, org_seq_residues,
                             end_of_section - start_of_section,
                             start_of_section, resind + seed_size - 1, 1,
                             coarse_seq->seq->residues,
                             coarse_seq->seq->length, 0) +
                 attempt_ext(current+seed_size-1, 1, org_seq_residues,
                             end_of_section - start_of_section,
                             start_of_section, resind, -1,
                             coarse_seq->seq->residues,
                             coarse_seq->seq->length, 0)) >
                 compress_flags.min_progress) {
                mseqs_rev = extend_match(mem,
                                         coarse_seq->seq->residues, 0,
                                         coarse_seq->seq->length, resind, -1,
                                         org_seq_residues, start_of_section,
                                         end_of_section,
                                         current + seed_size - 1, 1);

                mseqs_fwd = extend_match(mem,
                                         coarse_seq->seq->residues, 0,
                                         coarse_seq->seq->length,
                                         resind+seed_size-1, 1,
                                         org_seq_residues, start_of_section,
                                         end_of_section, current, -1);

                fwd_rlen = cb_align_length_nogaps(mseqs_fwd.rseq);
                rev_rlen = cb_align_length_nogaps(mseqs_rev.rseq);
                fwd_olen = cb_align_length_nogaps(mseqs_fwd.oseq);
                rev_olen = cb_align_length_nogaps(mseqs_rev.oseq);

                //If the match was too short, try the next seed.
                if (rev_olen+seed_size+fwd_olen-1 <
                      compress_flags.min_match_len) {
                    free(mseqs_fwd.rseq);
                    free(mseqs_fwd.oseq);
                    free(mseqs_rev.rseq);
                    free(mseqs_rev.oseq);
                    continue;
                }

                found_match = true;

                /*Concatenate the extensions and the k-mer's reverse complement
                  for the coarse sequence*/
                index = 0;

                alignment.ref = malloc((strlen(mseqs_rev.rseq)+seed_size+
                                        strlen(mseqs_fwd.rseq)+1)*
                                        sizeof(*(alignment.ref)));
                assert(alignment.ref);

                for (i = 0; mseqs_rev.rseq[i] != '\0'; i++)
                    alignment.ref[index++] = mseqs_rev.rseq[i];
                for (i = 0; i < seed_size; i++)
                    alignment.ref[index++] = revcomp[i];
                for (i = 0; mseqs_fwd.rseq[i] != '\0'; i++)
                    alignment.ref[index++] = mseqs_fwd.rseq[i];
                alignment.ref[index] = '\0';

                /*Concatenate the extensions and the k-mer's reverse complement
                  for the original sequence*/
                index = 0;

                alignment.org = malloc((strlen(mseqs_rev.oseq)+seed_size+
                                        strlen(mseqs_fwd.oseq)+1)*
                                        sizeof(*(alignment.org)));
                assert(alignment.org);

                for (i = 0; mseqs_rev.oseq[i] != '\0'; i++)
                    alignment.org[index++] = mseqs_rev.oseq[i];
                for (i = 0; i < seed_size; i++)
                    alignment.org[index++] = revcomp[i];
                for (i = 0; mseqs_fwd.oseq[i] != '\0'; i++)
                    alignment.org[index++] = mseqs_fwd.oseq[i];
                alignment.org[index] = '\0';

                alignment.length = mseqs_rev.olen + seed_size + mseqs_fwd.olen;


                //Make a new chunk for the parts of the chunk before the match.
                if (current - fwd_olen - start_of_section > 0) {
                    new_coarse_seq_id = add_without_match(coarse_db, org_seq,
                                                    start_of_section,
                                                    current - fwd_olen +
                                                      compress_flags.overlap,
                                                    seeds_mem);
                    cb_compressed_seq_addlink(cseq,
                      cb_link_to_coarse_init_nodiff(
                        new_coarse_seq_id,
                        start_of_section,
                        current - fwd_olen + compress_flags.overlap - 1,
                        0, current - fwd_olen + compress_flags.overlap
                                   - start_of_section - 1,
                        true));
                    chunks++;
                }

                //Add a link to the coarse sequence in the compressed sequence.
                cb_compressed_seq_addlink(cseq,
                                          cb_link_to_coarse_init(
                                            coarse_seq->id,
                                            current - fwd_olen,
                                            current + seed_size + rev_olen - 1,
                                            resind - rev_rlen,
                                            resind + seed_size + fwd_rlen - 1,
                                            alignment, false));

                //Add a link to the compressed sequence in the coarse sequence.
                cb_coarse_seq_addlink(coarse_seq,
                                      cb_link_to_compressed_init(
                                        org_seq->id,
                                        resind - rev_rlen,
                                        resind + seed_size + fwd_rlen - 1,
                                        current - fwd_olen,
                                        current + seed_size + rev_olen - 1,
                                        false));

                //Update the current position in the sequence.
                if (current + rev_olen <= org_seq_len-seed_size-1)
                    start_of_section = current + rev_olen -
                                       compress_flags.overlap + seed_size;
                else
                    start_of_section = current + rev_olen + seed_size;
                current = start_of_section - 1;
                end_of_chunk = min(start_of_section + max_chunk_size,
                                     org_seq_len);
                end_of_section = min(start_of_section + max_section_size,
                                     org_seq_len);

                chunks++;

                free(mseqs_fwd.rseq);
                free(mseqs_fwd.oseq);
                free(mseqs_rev.rseq);
                free(mseqs_rev.oseq);

                free(alignment.ref);
                free(alignment.org);
            }
        }

        /*If we have traversed an entire chunk of bases without finding a match,
         *then add the whole chunk as a sequence in the database and update
         *start_of_section, end_of_chunk, and end_of_section
         */
        if (current >= end_of_chunk - seed_size && !found_match) {
            new_coarse_seq_id = add_without_match(coarse_db, org_seq,
                                                  start_of_section,
                                                  end_of_chunk, seeds_mem);

            cb_compressed_seq_addlink(cseq, cb_link_to_coarse_init_nodiff(
                                              new_coarse_seq_id,
                                              start_of_section,
                                              end_of_chunk - 1, 0,
                                              end_of_chunk -
                                                start_of_section - 1,
                                              true));

            if (end_of_chunk < org_seq_len - seed_size - 1) {
                start_of_section = end_of_chunk - overlap;
                end_of_chunk = min(start_of_section + max_chunk_size,
                                   org_seq_len);
                end_of_section = min(start_of_section + max_section_size,
                                     org_seq_len);
                current = start_of_section - 1;
            }
            chunks++;
        }
    }

    free(matches);
    free(matches_temp);

    return cseq;
}

/*@param mem: A compress worker's alignment memory
 *@param rseq, oseq: The coarse and original sequences
 *@param rstart, rend, ostart, oend: The start and end of the sections of the
 *  sequences being aligned.
 *@param resind, current: The starting base indices in the coarse and original
 *  sequences.
 *@param dir1, dir2: The direction of the match extension for the coarse and
 *  original sequences (1 = forward, -1 = reverse)
 *
 *@return: An extend_match struct containing the aligned sequences and their
 *  lengths, including gaps from Needleman-Wunsch alignment.
 *
 *Extends a match found with attempt_ext by making alternating calls to
 *cb_align_ungapped and cb_align_nw until a bad window is found or the end of
 *one of the sequences is reached.  For more detail read the comments in the
 *function.
 */
struct extend_match
extend_match(struct cb_align_nw_memory *mem,
             char *rseq, int32_t rstart, int32_t rend, int32_t resind,
             int32_t dir1, char *oseq, int32_t ostart, int32_t oend,
             int32_t current, int32_t dir2){
    struct DSVector *rseq_segments = ds_vector_create(),
                    *oseq_segments = ds_vector_create();
    struct cb_alignment alignment;
    struct extend_match mseqs;
    struct ungapped_alignment ungapped;
    int32_t rlen, olen, i, j, m, matches_count, matches_index, max_section_size,
            rseq_len = 0, oseq_len = 0, dir_prod = dir1 * dir2;
    bool *matches, *matches_past_clump;
    bool found_bad_window;
    char *base_complement = "TNGNNNCNNNNNNNNNNNNANNNNNN";

    max_section_size = 2 * compress_flags.max_chunk_size;

    //Initialize the matches and matches_past_clump arrays.
    matches = mem->matches;
    matches_past_clump = mem->matches_past_clump;

    matches_index = compress_flags.gapped_window_size;
    for (i = 0; i < max_section_size; i++) {
        matches[i] = true;
        matches_past_clump[i] = true;
    }

    resind += dir1;
    current += dir2;

    rlen = rend - rstart;
    olen = oend - ostart;

    mseqs.rlen = 0;
    mseqs.olen = 0;
    mseqs.rseq = "";
    mseqs.oseq = "";

    while (true) {
        int dp_len1, dp_len2, i, r_align_len, o_align_len;
        char *r_segment, *o_segment;

        if (mseqs.rlen == rlen || mseqs.olen == olen)
            break;

        /*Get the maximum length for ungapped alignment and extend the match
          by that distance.*/
        ungapped = cb_align_ungapped(rseq, rstart, rend, dir1, resind,
                                     oseq, ostart, oend, dir2, current,
                                     matches,matches_past_clump,&matches_index);

        m = ungapped.length;
        found_bad_window = ungapped.found_bad_window;

        mseqs.rlen += m;
        mseqs.olen += m;

        if (m > 0) {
            r_segment = malloc((m+1)*sizeof(*r_segment));
            assert(r_segment);

            o_segment = malloc((m+1)*sizeof(*o_segment));
            assert(o_segment);

            for (i = 0; i < m; i++) {
                r_segment[i] = rseq[resind+dir1*i];
                o_segment[i] = dir_prod > 0 ? oseq[current+dir2*i] :
                                      base_complement[oseq[current+dir2*i]-'A'];
            }

            r_segment[m] = '\0';
            o_segment[m] = '\0';
            ds_vector_append(rseq_segments, (void *)r_segment);
            ds_vector_append(oseq_segments, (void *)o_segment);
        }

        resind += m * dir1;
        current += m * dir2;

        //End the extension if we found a bad window in ungapped alignment.
        if (found_bad_window)
            break;

        /*Carry out Needleman-Wunsch alignment and end the extension if we
          found a bad window or couldn't find a 4-mer match in the alignment.*/
        dp_len1 = max_dp_len(resind-rstart, dir1, rend-rstart);
        dp_len2 = max_dp_len(current-ostart, dir2, oend-ostart);

        alignment = cb_align_nw(mem, rseq, dp_len1, resind, dir1,
                                oseq, dp_len2, current, dir2,
                                matches, &matches_index);

        //End the extension if we found a bad window in ungapped alignment.
        if (alignment.length == -1)
            break;

        matches_count = 0;

        /*Check for a bad window manually and end the extension if a bad
          window is found.*/
        for (i = matches_index - 100; i < matches_index; i++)
            if (matches[i])
                matches_count++;
        if (matches_count < compress_flags.window_ident_thresh)
            break;

        ds_vector_append(rseq_segments, (void *)alignment.ref);
        ds_vector_append(oseq_segments, (void *)alignment.org);

        r_align_len = strlen(alignment.ref);
        o_align_len = strlen(alignment.org);

        /*Update the lengths of the alignments and the indices of the
          sequences.*/
        mseqs.rlen += r_align_len;
        mseqs.olen += o_align_len;
        resind += cb_align_length_nogaps(alignment.ref) * dir1;
        current += cb_align_length_nogaps(alignment.org) * dir2;
    }

    mseqs.rseq = malloc((mseqs.rlen+1)*sizeof(*(mseqs.rseq)));
    assert(mseqs.rseq);

    mseqs.oseq = malloc((mseqs.olen+1)*sizeof(*(mseqs.oseq)));
    assert(mseqs.oseq);


    //Copy each segment into the overall extension.
    for (i = 0; i < rseq_segments->size; i++) {
        char *r_segment = (char *)ds_vector_get(rseq_segments, i),
             *o_segment = (char *)ds_vector_get(oseq_segments, i);
        for (j = 0; r_segment[j] != '\0'; j++)
            mseqs.rseq[rseq_len++] = r_segment[j];
        for (j = 0; o_segment[j] != '\0'; j++)
            mseqs.oseq[oseq_len++] = o_segment[j];
    }

    mseqs.rseq[mseqs.rlen] = '\0';
    mseqs.oseq[mseqs.olen] = '\0';

    /*If the direction of the coarse sequence is reversed, reverse both
      extensions.*/
    if (dir1 < 0) {
        char *temp = malloc((mseqs.rlen+1)*sizeof(*temp));
        assert(temp);

        for (i = 0; i < mseqs.rlen; i++)
            temp[i] = mseqs.rseq[mseqs.rlen-i-1];
        temp[mseqs.rlen] = '\0';
        free(mseqs.rseq);
        mseqs.rseq = temp;

        temp = malloc((mseqs.olen+1)*sizeof(*temp));
        assert(temp);

        for (i = 0; i < mseqs.olen; i++)
            temp[i] = mseqs.oseq[mseqs.olen-i-1];
        temp[mseqs.olen] = '\0';
        free(mseqs.oseq);
        mseqs.oseq = temp;
    }

    ds_vector_free(rseq_segments);
    ds_vector_free(oseq_segments);
    return mseqs;
}

/*Creates a new coarse sequence for a section of the original DNA sequence that
*does not have a match, such as if the maximum chunk length was traversed in
*the sequence without finding a match or if there is any DNA in the sequence
*before the latest match that isn't part of the match.
*
*In addition to creating a new coarse sequence, add_without_match also adds a
*a link to the sequence being compressed to the new coarse sequence.
*/
static int32_t add_without_match(struct cb_coarse *coarse_db,
                                 struct cb_seq *org_seq,
                                 int32_t ostart, int32_t oend,
                                 struct cb_seeds_add_memory *seeds_mem){
    struct cb_link_to_compressed *link = NULL;
    struct cb_coarse_seq *coarse_seq =
      cb_coarse_add(coarse_db, org_seq->residues, ostart, oend, seeds_mem);
    cb_coarse_seq_addlink(coarse_seq,
      cb_link_to_compressed_init(org_seq->id, 0, oend - ostart - 1,
                                 ostart, oend - 1, true));
    link = coarse_seq->links;
    return coarse_seq->id;
}

static int32_t min(int32_t a, int32_t b){
    if (a < b)
        return a;
    return b;
}
