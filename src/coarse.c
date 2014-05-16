#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bitpack.h"
#include "coarse.h"
#include "DNAutils.h"
#include "seeds.h"
#include "seq.h"

/*Takes in the size of the k-mers that will be used in compression and file
  pointers for the database and returns a newly-created coarse database.*/
struct cb_coarse *
cb_coarse_init(int32_t seed_size,
                FILE *file_fasta, FILE *file_seeds, FILE *file_links,
                FILE *file_links_index, FILE *file_fasta_index,
                FILE *file_params)
{
    struct cb_coarse *coarse_db;
    int32_t errno;

    coarse_db = malloc(sizeof(*coarse_db));
    assert(coarse_db);

    coarse_db->seqs = ds_vector_create_capacity(10000000);
    coarse_db->seeds = cb_seeds_init(seed_size);
    coarse_db->dbsize = (uint64_t)0;

    /*Initialize the file pointers*/
    coarse_db->file_fasta = file_fasta;
    coarse_db->file_seeds = file_seeds;
    coarse_db->file_links = file_links;
    coarse_db->file_fasta_index = file_fasta_index;
    coarse_db->file_links_index = file_links_index;
    coarse_db->file_params = file_params;

    if (0 != (errno = pthread_rwlock_init(&coarse_db->lock_seq, NULL))) {
        fprintf(stderr, "Could not create rwlock. Errno: %d\n", errno);
        exit(1);
    }

    return coarse_db;
}

/*Takes in a coarse database, closes its file pointers, and frees the struct
  for the database.*/
void
cb_coarse_free(struct cb_coarse *coarse_db)
{
    int32_t errno;
    int32_t i;

    fclose(coarse_db->file_fasta);
    fclose(coarse_db->file_seeds);
    fclose(coarse_db->file_links);
    fclose(coarse_db->file_links_index);
    fclose(coarse_db->file_fasta_index);
    fclose(coarse_db->file_params);

    if (0 != (errno = pthread_rwlock_destroy(&coarse_db->lock_seq))) {
        fprintf(stderr, "Could not destroy rwlock. Errno: %d\n", errno);
        exit(1);
    }
    for (i = 0; i < coarse_db->seqs->size; i++)
        cb_coarse_seq_free(
            (struct cb_coarse_seq *) ds_vector_get(coarse_db->seqs, i));

    ds_vector_free_no_data(coarse_db->seqs);
    cb_seeds_free(coarse_db->seeds);
    free(coarse_db);
}

struct cb_coarse_seq *
cb_coarse_add(struct cb_coarse *coarse_db,
               char *residues, int32_t start, int32_t end)
{
    struct cb_coarse_seq *seq;
    int32_t id;

    pthread_rwlock_wrlock(&coarse_db->lock_seq);
    id = coarse_db->seqs->size;
    seq = cb_coarse_seq_init(id, residues, start, end);
    ds_vector_append(coarse_db->seqs, (void*) seq);
    pthread_rwlock_unlock(&coarse_db->lock_seq);

    cb_seeds_add(coarse_db->seeds, seq);

    return seq;
}

struct cb_coarse_seq *
cb_coarse_get(struct cb_coarse *coarse_db, int32_t i)
{
    struct cb_coarse_seq *seq;

    pthread_rwlock_rdlock(&coarse_db->lock_seq);
    seq = (struct cb_coarse_seq *) ds_vector_get(coarse_db->seqs, i);
    pthread_rwlock_unlock(&coarse_db->lock_seq);

    return seq;
}

/*Outputs the sequences in the coarse database to a FASTA file in plain text,
 *outputs the links to the compressed database in a binary format, and outputs
 *the size of the database to the params file.
 */
void
cb_coarse_save_binary(struct cb_coarse *coarse_db)
{
    struct cb_coarse_seq *seq;
    struct cb_link_to_compressed *link;
    int64_t i;
    int j;

    /*Keeps track of the index to be printed to coarse.links.index*/
    uint64_t index = (uint64_t)0;
    
    int16_t mask = (((int16_t)1)<<8)-1;

    for (i = 0; i < coarse_db->seqs->size; i++) {
        uint64_t coarse_fasta_loc;
        char *fasta_output = malloc(30000*sizeof(*fasta_output));
        char *link_header = malloc(30*sizeof(*fasta_output));
        seq = (struct cb_coarse_seq *) ds_vector_get(coarse_db->seqs, i);

        /*At the start of outputting each sequence, output the indices for the
          coarse links and FASTA files to their index files.*/
        output_int_to_file(index, 8, coarse_db->file_links_index);
        coarse_fasta_loc = ftell(coarse_db->file_fasta);
        output_int_to_file(coarse_fasta_loc, 8, coarse_db->file_fasta_index);

        /*Output the FASTA sequence to the coarse FASTA file*/
        sprintf(fasta_output, "> %ld\n%s\n", i, seq->seq->residues);
        for(j = 0; fasta_output[j] != '\0'; j++)
            putc(fasta_output[j], coarse_db->file_fasta);

        /*Output the link header for the sequence to the coarse links file*/
        sprintf(link_header, "> %ld\n", i);
        for(j = 0; link_header[j] != '\0'; j++)
            putc(link_header[j], coarse_db->file_links);

        /*For each character outputted to the coarse links file, increment
          "index"*/
        index += strlen(link_header);

        free(fasta_output);
        free(link_header);

        /*Output all links for the current sequence to the coarse links file*/
        for (link = seq->links; link != NULL; link = link->next) {
            int j;
            char *id_bytes = read_int_to_bytes(link->org_seq_id, 8);
            for (j = 0; j < 8; j++)
                putc(id_bytes[j], coarse_db->file_links);
            /*Convert the start and end indices for the link to two
              characters.*/
            int16_t coarse_start = (int16_t)link->coarse_start;
            int16_t coarse_end   = (int16_t)link->coarse_end;
            char coarse_start_left  = (coarse_start >> 8) & mask;
            char coarse_start_right = coarse_start & mask;
            char coarse_end_left    = (coarse_end >> 8) & mask;
            char coarse_end_right   = coarse_end & mask;
            /*Prints the binary representations of the indices and the
              direction of the link to the links file*/
            putc(coarse_start_left, coarse_db->file_links);
            putc(coarse_start_right, coarse_db->file_links);
            putc(coarse_end_left, coarse_db->file_links);
            putc(coarse_end_right, coarse_db->file_links);
            output_int_to_file(link->original_start, 8, coarse_db->file_links);
            output_int_to_file(link->original_end, 8, coarse_db->file_links);
            putc((link->dir?'0':'1'), coarse_db->file_links);

            index += 29;

            /*0 is used as a delimiter to signify that there are more links
              for this sequence*/
            if (link->next != NULL) {
                putc(0, coarse_db->file_links);
                index++;
            }
        }
        /*'#' is used as a delimiter to signify the last link of the sequence*/
        if (i+1 < coarse_db->seqs->size){
            putc('#', coarse_db->file_links);
            index++;
        }
    }
    output_int_to_file(coarse_db->dbsize, 8, coarse_db->file_params);
    putc('\n', coarse_db->file_params);
    putc('\n', coarse_db->file_links);
}

void
cb_coarse_save_plain(struct cb_coarse *coarse_db)
{
    struct cb_coarse_seq *seq;
    struct cb_link_to_compressed *link;
    int32_t i;

    for (i = 0; i < coarse_db->seqs->size; i++) {
        seq = (struct cb_coarse_seq *) ds_vector_get(coarse_db->seqs, i);
        fprintf(coarse_db->file_fasta, "> %d\n%s\n", i, seq->seq->residues);

        fprintf(coarse_db->file_links, "> %d\n", i);
        for (link = seq->links; link != NULL; link = link->next)
            fprintf(coarse_db->file_links_index,
                "original sequence id: %d, reference range: (%d, %d), "
                  "direction: %c\n",
                link->org_seq_id, link->coarse_start, link->coarse_end,
                (link->dir?'0':'1'));
    }
}

void
cb_coarse_save_seeds_binary(struct cb_coarse *coarse_db)
{
    struct cb_coarse_seq *seq;
    struct cb_link_to_compressed *link;
    int32_t i, j;
    char *kmer;
    uint32_t mask = (uint32_t)3;

    for (i = 0; i < coarse_db->seeds->locs_length; i++) {
        kmer = unhash_kmer(coarse_db->seeds, i);
        struct cb_seed_loc *loc = cb_seeds_lookup(coarse_db->seeds, kmer);
        if (loc) {
            output_int_to_file(i, 4, coarse_db->file_seeds);    
            struct cb_seed_loc *loc_first = loc;
            while (loc) {
                output_int_to_file(loc->coarse_seq_id,4,coarse_db->file_seeds);
                output_int_to_file(loc->residue_index,2,coarse_db->file_seeds);
                loc = loc->next;
                if (loc) putc((char)0, coarse_db->file_seeds);
            }
            putc((char)1, coarse_db->file_seeds);
            cb_seed_loc_free(loc_first);
        }
        free(kmer);
    }
    putc('\n', coarse_db->file_seeds);
}

void
cb_coarse_save_seeds_plain(struct cb_coarse *coarse_db)
{
    struct cb_coarse_seq *seq;
    struct cb_link_to_compressed *link;
    int32_t i, j;
    uint32_t i2;
    uint32_t mask = (uint32_t)3;
    char *kmer;

    for (i = 0; i < coarse_db->seeds->locs_length; i++) {
        struct cb_seed_loc *loc;
        i2 = (uint32_t)0;
        for (j = 0; j < coarse_db->seeds->seed_size; j++) {
            i2 <<= 2;
            i2 |= ((i >> (2*j)) & mask);
        }
        kmer = unhash_kmer(coarse_db->seeds, i2);
        fprintf(coarse_db->file_seeds, "%s\n", kmer);
        loc = cb_seeds_lookup(coarse_db->seeds, kmer);
        struct cb_seed_loc *loc_first = loc;
        while (loc) {
            if (loc->coarse_seq_id < 500)
                fprintf(coarse_db->file_seeds,"(%d, %d) > ",
                        loc->coarse_seq_id, loc->residue_index);
            loc = loc->next;
        }
        cb_seed_loc_free(loc_first);
        fprintf(coarse_db->file_seeds, "\n");
        free(kmer);
    }
}

struct cb_coarse_seq *
cb_coarse_seq_init(int32_t id, char *residues, int32_t start, int32_t end)
{
    struct cb_coarse_seq *seq;
    int32_t errno;

    seq = malloc(sizeof(*seq));
    assert(seq);

    seq->id = id;
    seq->seq = cb_seq_init_range(id, "", residues, start, end);
    seq->links = NULL;

    if (0 != (errno = pthread_rwlock_init(&seq->lock_links, NULL))) {
        fprintf(stderr, "Could not create rwlock. Errno: %d\n", errno);
        exit(1);
    }

    return seq;
}

void
cb_coarse_seq_free(struct cb_coarse_seq *seq)
{
    struct cb_link_to_compressed *link1, *link2;
    int32_t errno;

    if (0 != (errno = pthread_rwlock_destroy(&seq->lock_links))) {
        fprintf(stderr, "Could not destroy rwlock. Errno: %d\n", errno);
        exit(1);
    }
    for (link1 = seq->links; link1 != NULL; ) {
        link2 = link1->next;
        cb_link_to_compressed_free(link1);
        link1 = link2;
    }

    cb_seq_free(seq->seq);
    free(seq);
}

void
cb_coarse_seq_addlink(struct cb_coarse_seq *seq,
                       struct cb_link_to_compressed *newlink)
{
    struct cb_link_to_compressed *link;

    assert(newlink->next == NULL);
    pthread_rwlock_wrlock(&seq->lock_links);
    if (seq->links == NULL) {
        seq->links = newlink;
        pthread_rwlock_unlock(&seq->lock_links);
        return;
    }

    for (link = seq->links; link->next != NULL; link = link->next);
    link->next = newlink;
    pthread_rwlock_unlock(&seq->lock_links);
}

struct cb_link_to_compressed *
cb_link_to_compressed_init(int32_t org_seq_id, int16_t coarse_start,
                            int16_t coarse_end, uint64_t original_start,
                            uint64_t original_end, bool dir)
{
    struct cb_link_to_compressed *link;

    link = malloc(sizeof(*link));
    assert(link);

    link->org_seq_id = org_seq_id;
    link->coarse_start = coarse_start;
    link->coarse_end = coarse_end;
    link->original_start = original_start;
    link->original_end = original_end;
    link->dir = dir;
    link->next = NULL;

    return link;
}

void
cb_link_to_compressed_free(struct cb_link_to_compressed *link)
{
    free(link);
}

/*A function for getting the header for an entry in the coarse links file.
  Returns NULL if EOF is found before a newline.*/
char *get_coarse_header(FILE *f){
    int c = 0;
    char *header = malloc(30*sizeof(*header));
    int header_length = 30;
    int i = 0;
    while (c != EOF && c != '\n') {
        c = getc(f);
        if (c != EOF) {
            header[i] = (char)c;
            i++;
            if (i == header_length - 1) {
                header_length *= 2;
                header = realloc(header, header_length*sizeof(*header));
            }
        }
        else
            return NULL;
    }
    header[i] = '\0';
    header = realloc(header, (i+1)*sizeof(*header));
    return header;
}

/*Reads one link from a file with the links to the compressed database and
  converts its data to a struct cb_link_to_compressed*/
struct cb_link_to_compressed *read_coarse_link(FILE *f){
    struct cb_link_to_compressed *link = malloc(sizeof(*link));
    link->org_seq_id = (uint64_t)read_int_from_file(8, f);
    if (feof(f)) {
        free(link);
        return NULL;
    }

    link->coarse_start = (uint16_t)read_int_from_file(2, f);
    if (feof(f)) {
        free(link);
        return NULL;
    }

    link->coarse_end = (uint16_t)read_int_from_file(2, f);
    if (feof(f)) {
        free(link);
        return NULL;
    }

    link->original_start = (uint64_t)read_int_from_file(8, f);
    if (feof(f)) {
        free(link);
        return NULL;
    }

    link->original_end = (uint64_t)read_int_from_file(8, f);
    if (feof(f)) {
        free(link);
        return NULL;
    }

    link->dir = getc(f) == '0';
    link->next = NULL;
    return link;
}

/*Takes in a pointer to the coarse.links file generated by cablast-compress and
 *returns a vector containing all of the links in the coarse database entry for
 *the sequence the file pointer currently points to.  For this function to work
 *properly the file pointer must be pointing to the start of the header of a
 *sequence entry in the links file.
 */
struct DSVector *get_coarse_sequence_links(FILE *f){
    struct DSVector *links = ds_vector_create();
    char *h = get_coarse_header(f);
    if (h == NULL) {
        ds_vector_free(links);
        return NULL;
    }
    free(h);

    while (true) {
        char c = 1;
        struct cb_link_to_compressed *current_link = read_coarse_link(f);
        if (current_link == NULL)
            break;
        ds_vector_append(links, (void *)current_link);
        c = getc(f);
        if (c == '#')
            break;
    }
    return links;
}

/*A wrapper function for get_coarse_sequence_links that handles fseek calls;
 *seeks to the index in the coarse.links file for the sequence at index id
 *and then calls get_coarse_sequence_links.  If fseek is successful, then
 *return the coarse sequence links for the coarse sequence at index id.
 *Otherwise, return NULL.
 */
struct DSVector *get_coarse_sequence_links_at(FILE *links, FILE *index,
                                                           int32_t id){
    int64_t offset = cb_coarse_find_offset(index, id);
    if (offset < 0)
        return NULL;
    bool fseek_success = fseek(links, offset, SEEK_SET) == 0;
    if (!fseek_success) { 
        fprintf(stderr, "Error in seeking to offset %lu\n", offset);
        return NULL;
    }
    return get_coarse_sequence_links(links);
}


/*Takes in an index file from a coarse database and the ID number of the
  sequence in the corresponding database file that the user wants to seek to
  and returns the byte offset of the sequence in its database file.*/
int64_t cb_coarse_find_offset(FILE *index_file, int id){
    int i;
    int try_off = id * 8;
    bool fseek_success = fseek(index_file, try_off, SEEK_SET) == 0;
    int64_t mask = make_mask(8);
    int64_t offset = (int64_t)(-1);
    if (!fseek_success) {
        fprintf(stderr, "Error in seeking to offset %d\n", try_off);
        return (int64_t)(-1);
    }
    for (i = 0; i < 8; i++) {
        int64_t current_byte=((int64_t)(getc(index_file))&mask);
        offset <<= 8;
        offset |= current_byte;
    }
    return offset;
}

/*Takes in as arguments a coarse database and the ID number of the sequence in
 *the coarse FASTA file to read in and gets a struct fasta_seq for that
 *sequence.
 */
struct fasta_seq *cb_coarse_read_fasta_seq(struct cb_coarse *coarsedb,
                                                               int id){
    int64_t offset = cb_coarse_find_offset(coarsedb->file_fasta_index, id);
    if (offset < 0)
        return NULL;
    bool fseek_success = fseek(coarsedb->file_fasta, offset, SEEK_SET) == 0;
    if (!fseek_success) {
        fprintf(stderr, "Error in seeking to offset %d\n", offset);
        return NULL;
    }
    return fasta_read_next(coarsedb->file_fasta, "");
}
