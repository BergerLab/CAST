#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bitpack.h"
#include "coarse.h"
#include "DNAutils.h"
#include "fasta.h"
#include "seeds.h"
#include "seq.h"

/*Takes in the size of the k-mers that will be used in compression, the file
 *pointers for the database, and whether or not to load in the coarse residues
 *from the coarse database's coarse FASTA file and returns a newly-created
 *coarse database.
 */
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

    coarse_db->seqs         = ds_vector_create_capacity(10000000);
    coarse_db->seeds        = cb_seeds_init(seed_size);
    coarse_db->dbsize       = (uint64_t)0;

    /*Initialize the file pointers*/
    coarse_db->file_fasta       = file_fasta;
    coarse_db->file_seeds       = file_seeds;
    coarse_db->file_links       = file_links;
    coarse_db->file_fasta_index = file_fasta_index;
    coarse_db->file_links_index = file_links_index;
    coarse_db->file_params      = file_params;

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

    /*if (coarse_db->all_residues != NULL)
        free(coarse_db->all_residues);*/

    free(coarse_db);
}

/*Takes in the coarse database and the residues and original start and end
 *indices for a sequence to be added to a coarse database and adds a sequence
 *created from those residues while also adding the sequence to the seeds
 *table.
 */
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

/*Get the coarse sequence in the coarse database at index i.*/
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
    /*Keeps track of the index to be printed to coarse.links.index*/
    uint64_t index = (uint64_t)0;
    int16_t mask = (int16_t)0xff;
    int j;

    for (i = 0; i < coarse_db->seqs->size; i++) {
        uint64_t coarse_fasta_loc;
        char *fasta_output, *link_header;

        fasta_output = malloc(30000*sizeof(*fasta_output));
        assert(fasta_output);

        link_header = malloc(30*sizeof(*fasta_output));
        assert(link_header);

        seq = (struct cb_coarse_seq *) ds_vector_get(coarse_db->seqs, i);

        /*At the start of outputting each sequence, output the indices for the
          coarse links and FASTA files to their index files.*/
        output_int_to_file(index, 8, coarse_db->file_links_index);
        coarse_fasta_loc = ftell(coarse_db->file_fasta);
        output_int_to_file(coarse_fasta_loc, 8, coarse_db->file_fasta_index);

        /*Output the FASTA sequence to the coarse FASTA file*/
        sprintf(fasta_output, ">%ld\n%s\n", i, seq->seq->residues);
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
            int32_t j;
            /*Convert the start and end indices for the link to two
              characters.*/
            int16_t coarse_start = (int16_t)link->coarse_start,
                    coarse_end   = (int16_t)link->coarse_end;
            char *id_bytes          = read_int_to_bytes(link->org_seq_id, 8),
                 coarse_start_left  = (coarse_start >> 8) & mask,
                 coarse_start_right = coarse_start & mask,
                 coarse_end_left    = (coarse_end >> 8) & mask,
                 coarse_end_right   = coarse_end & mask;

            for (j = 0; j < 8; j++)
                putc(id_bytes[j], coarse_db->file_links);

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

/*Output the coarse database's sequences and links in a human-readable format.
  Used for debugging purposes.*/
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

/*Takes in the coarse database and saves its seeds table in a binary format.*/
void
cb_coarse_save_seeds_binary(struct cb_coarse *coarse_db)
{
    struct cb_coarse_seq *seq;
    struct cb_link_to_compressed *link;
    int32_t i, j;
    uint32_t mask = (uint32_t)3;
    char *kmer;

    for (i = 0; i < coarse_db->seeds->locs_length; i++) {
        struct cb_seed_loc *loc;
        kmer = unhash_kmer(coarse_db->seeds, i);

        loc = cb_seeds_lookup(coarse_db->seeds, kmer);

        if (loc) {
            struct cb_seed_loc *loc_first = loc;
            output_int_to_file(i, 4, coarse_db->file_seeds);    
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

/*Takes in the coarse database and saves its seeds table in a human-readable
  format.  Used for debugging.*/
void
cb_coarse_save_seeds_plain(struct cb_coarse *coarse_db)
{
    struct cb_coarse_seq *seq;
    struct cb_link_to_compressed *link;
    int32_t i, j;
    uint32_t i2 = (uint32_t)0, mask = (uint32_t)3;
    char *kmer;

    for (i = 0; i < coarse_db->seeds->locs_length; i++) {
        struct cb_seed_loc *loc, *loc_first = loc;

        i2 = (uint32_t)0;
        for (j = 0; j < coarse_db->seeds->seed_size; j++) {
            i2 <<= 2;
            i2 |= ((i >> (2*j)) & mask);
        }

        kmer = unhash_kmer(coarse_db->seeds, i2);
        fprintf(coarse_db->file_seeds, "%s\n", kmer);
        loc = cb_seeds_lookup(coarse_db->seeds, kmer);
        loc_first = loc;

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

/*Loads all of the residues in the coarse database's links file into the coarse
  database's links vector.*/
void cb_coarse_get_all_links(struct cb_coarse_db_read *coarse_db){
    struct DSVector *link_vectors = ds_vector_create();
    int32_t num_link_vectors = 0, i = 0, j = 0;

    bool fseek_success =
      fseek(coarse_db->coarsedb->file_links_index, 0, SEEK_END) == 0;
    if (!fseek_success)
        fprintf(stderr, "Error in seeking to end of FASTA index file\n");
    fseek(coarse_db->coarsedb->file_links_index, 0, SEEK_SET) == 0;

    coarse_db->links = ds_vector_create();

    num_link_vectors = ftell(coarse_db->coarsedb->file_fasta_index)/8;
    for (i = 0; i < num_link_vectors; i++) {
        struct DSVector *links =
          get_coarse_sequence_links_at(coarse_db->coarsedb->file_links,
                                       coarse_db->coarsedb->file_links_index,i);
        for (j = 0; j < links->size; j++)
            ds_vector_append(coarse_db->links, ds_vector_get(links, j));
        ds_vector_free_no_data(links);
    }
}

/*Loads all of the residues in the coarse database's FASTA file into the coarse
  database's all_residues string*/
void cb_coarse_get_all_residues(struct cb_coarse_db_read *coarse_db){
    struct DSVector *fasta_seqs = ds_vector_create();
    int32_t num_fasta_entries = 0, num_bases = 0,
            i = 0, j = 0, bases_copied = 0;

    bool fseek_success =
      fseek(coarse_db->coarsedb->file_fasta_index, 0, SEEK_END) == 0;
    if (!fseek_success)
        fprintf(stderr, "Error in seeking to end of FASTA index file\n");

    num_fasta_entries = ftell(coarse_db->coarsedb->file_fasta_index)/8;
    fseek(coarse_db->coarsedb->file_fasta_index, 0, SEEK_SET) == 0;

    for (i = 0; i < num_fasta_entries; i++) {
        struct fasta_seq *current_seq = cb_coarse_read_fasta_seq_r(coarse_db, i);
        if (!current_seq)
            fprintf(stderr, "Error getting FASTA sequence #%d in "
                            "cb_coarse_get_all_residues\n", i);
        ds_vector_append(fasta_seqs, (void *)current_seq);
        num_bases += strlen(current_seq->seq);
    }

    coarse_db->all_residues =
      malloc((num_bases+1)*sizeof(*(coarse_db->all_residues)));
    assert(coarse_db->all_residues);

    for (i = 0; i < num_fasta_entries; i++) {
        struct fasta_seq *sequence =
            (struct fasta_seq *)ds_vector_get(fasta_seqs, i);
        for (j = 0; j < strlen(sequence->seq); j++) {
            coarse_db->all_residues[bases_copied++] = sequence->seq[j];
        }
        fasta_free_seq(sequence);
    }
    coarse_db->all_residues[num_bases] = '\0';
    ds_vector_free_no_data(fasta_seqs);
}

/*Takes in an ID number and the residues and original start and end indices
 *for a sequence to be added to a coarse database and creates a coarse
 *sequence from that information.
 */
struct cb_coarse_seq *
cb_coarse_seq_init(int32_t id, char *residues, int32_t start, int32_t end)
{
    struct cb_coarse_seq *seq;
    int32_t errno;

    seq = malloc(sizeof(*seq));
    assert(seq);

    seq->id    = id;
    seq->seq   = cb_seq_init_range(id, "", residues, start, end);
    seq->links = NULL;

    if (0 != (errno = pthread_rwlock_init(&seq->lock_links, NULL))) {
        fprintf(stderr, "Could not create rwlock. Errno: %d\n", errno);
        exit(1);
    }

    return seq;
}

/*Frees a coarse sequence.*/
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

/*Takes in a coarse sequence and a link_to_compressed and adds the link to the
  sequence at the end of the linked list of links.*/
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

/*A function for getting the header for an entry in the coarse links file.
  Returns NULL if EOF is found before a newline.*/
char *get_coarse_header(FILE *f){
    int header_length = 30, c = 0, i = 0;
    char *header;

    header = malloc(30*sizeof(*header));
    assert(header);

    while (c != EOF && c != '\n') {
        c = getc(f);
        if (c != EOF) {
            header[i] = (char)c;
            i++;
            if (i == header_length - 1) {
                header_length *= 2;

                header = realloc(header, header_length*sizeof(*header));
                assert(header);
            }
        }
        else
            return NULL;
    }
    header[i] = '\0';

    header = realloc(header, (i+1)*sizeof(*header));
    assert(header);

    return header;
}

/*Reads one link from a file with the links to the compressed database and
  converts its data to a struct cb_link_to_compressed*/
struct cb_link_to_compressed *read_coarse_link(FILE *f){
    struct cb_link_to_compressed *link = malloc(sizeof(*link));
    assert(link);

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
        struct cb_link_to_compressed *current_link = read_coarse_link(f);
        char c = 1;
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
    bool fseek_success;

    if (offset < 0)
        return NULL;
    fseek_success = fseek(links, offset, SEEK_SET) == 0;
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
    int64_t mask = make_mask(8),
            offset = (int64_t)(-1);
    int i, try_off = id * 8;
    bool fseek_success = fseek(index_file, try_off, SEEK_SET) == 0;

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
    bool fseek_success;

    if (offset < 0)
        return NULL;

    fseek_success = fseek(coarsedb->file_fasta, offset, SEEK_SET) == 0;
    if (!fseek_success) {
        fprintf(stderr, "Error in seeking to offset %d\n", offset);
        return NULL;
    }
    return fasta_read_next(coarsedb->file_fasta, "");
}

/*Loads a cb_coarse_db_read struct, which includes a coarse database as well
  as data from the files being read*/
struct cb_coarse_db_read *
cb_coarse_read_init(int32_t seed_size,
                    FILE *file_fasta, FILE *file_seeds, FILE *file_links,
                    FILE *file_links_index, FILE *file_fasta_index,
                    FILE *file_params, bool load_coarse_residues,
                    bool load_coarse_links){
    struct cb_coarse_db_read *coarsedb = malloc(sizeof(*coarsedb));
    assert(coarsedb);

    coarsedb->coarsedb = cb_coarse_init(seed_size, file_fasta, file_seeds,
                                        file_links, file_links_index,
                                        file_fasta_index, file_params);
    coarsedb->link_inds_by_block = ds_vector_create();

    coarsedb->all_residues = NULL;
    coarsedb->links        = NULL;

    /*If the --load-coarse-residues search flag was passed in, load the coarse
      residues into coarse_db->all_residues.*/
    if (load_coarse_residues)
        cb_coarse_get_all_residues(coarsedb);

    if (load_coarse_links)
        cb_coarse_get_all_links(coarsedb);

    return coarsedb;
}

/*Frees a cb_coarse_db_read struct.*/
void
cb_coarse_db_read_free(struct cb_coarse_db_read *coarsedb){
    int i;

    cb_coarse_free(coarsedb->coarsedb);
    if (coarsedb->all_residues)
        free(coarsedb->all_residues);
    if (coarsedb->links) {
        for (i = 0; i < coarsedb->links->size; i++)
            cb_link_to_compressed_free(
              (struct cb_link_to_compressed *)ds_vector_get(coarsedb->links,i));
        ds_vector_free_no_data(coarsedb->links);
    }
    ds_vector_free(coarsedb->link_inds_by_block);
    free(coarsedb);
}

struct cb_coarse_seq * cb_coarse_get_r(struct cb_coarse_db_read *coarse_db,
                                       int32_t i){
    return cb_coarse_get(coarse_db->coarsedb, i);
}

struct fasta_seq *cb_coarse_read_fasta_seq_r(struct cb_coarse_db_read *coarsedb,
                                             int id){
    return cb_coarse_read_fasta_seq(coarsedb->coarsedb, id);
}

