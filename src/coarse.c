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
#include "util.h"

/*Takes in the size of the k-mers that will be used in compression, the file
 *pointers for the database, and whether or not to load in the coarse residues
 *from the coarse database's coarse FASTA file and returns a newly-created
 *coarse database.
 */
struct cb_coarse *
cb_coarse_init(int32_t seed_size,
               FILE *file_fasta, FILE *file_seeds, FILE *file_links,
               FILE *file_links_index, FILE *file_links_base_index,
               FILE *file_links_count_index, FILE *file_fasta_index,
               FILE *file_fasta_base_index, FILE *file_params, bool read){
    int32_t errno;

    struct cb_coarse *coarse_db = malloc(sizeof(*coarse_db));
    assert(coarse_db);

    coarse_db->seqs   = ds_vector_create_capacity(10000000);
    coarse_db->seeds  = read ? NULL : cb_seeds_init(seed_size);
    coarse_db->dbsize = (uint64_t)0;

    /*Initialize the file pointers*/
    coarse_db->file_fasta             = file_fasta;
    coarse_db->file_seeds             = file_seeds;
    coarse_db->file_links             = file_links;
    coarse_db->file_links_base_index  = file_links_base_index;
    coarse_db->file_links_count_index = file_links_count_index;
    coarse_db->file_fasta_index       = file_fasta_index;
    coarse_db->file_fasta_base_index  = file_fasta_base_index;
    coarse_db->file_links_index       = file_links_index;
    coarse_db->file_params            = file_params;

    if (0 != (errno = pthread_rwlock_init(&coarse_db->lock_seq, NULL))) {
        fprintf(stderr, "Could not create rwlock. Errno: %d\n", errno);
        exit(1);
    }

    return coarse_db;
}

/*Takes in a coarse database, closes its file pointers, and frees the struct
  for the database.*/
void cb_coarse_free(struct cb_coarse *coarse_db){
    int32_t errno;
    int32_t i;

    fclose(coarse_db->file_fasta);
    fclose(coarse_db->file_seeds);
    fclose(coarse_db->file_links);
    fclose(coarse_db->file_links_index);
    fclose(coarse_db->file_links_base_index);
    fclose(coarse_db->file_links_count_index);
    fclose(coarse_db->file_fasta_index);
    fclose(coarse_db->file_fasta_base_index);
    fclose(coarse_db->file_params);

    if (0 != (errno = pthread_rwlock_destroy(&coarse_db->lock_seq))) {
        fprintf(stderr, "Could not destroy rwlock. Errno: %d\n", errno);
        exit(1);
    }
    for (i = 0; i < coarse_db->seqs->size; i++)
        cb_coarse_seq_free(
            (struct cb_coarse_seq *)ds_vector_get(coarse_db->seqs, i));

    ds_vector_free_no_data(coarse_db->seqs);
    if (coarse_db->seeds != NULL)
        cb_seeds_free(coarse_db->seeds);

    free(coarse_db);
}

/*Takes in the coarse database and the residues and original start and end
 *indices for a sequence to be added to a coarse database and adds a sequence
 *created from those residues while also adding the sequence to the seeds
 *table.
 */
struct cb_coarse_seq *cb_coarse_add(struct cb_coarse *coarse_db, char *residues,
                                    int32_t start, int32_t end){
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
struct cb_coarse_seq *cb_coarse_get(struct cb_coarse *coarse_db, int32_t i){
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
void cb_coarse_save_binary(struct cb_coarse *coarse_db){
    struct cb_coarse_seq *seq;
    struct cb_link_to_compressed *link;
    int64_t i;

    /*Keeps track of the indices to be printed to coarse.links.index and
      coarse_fasta_base_index.*/
    uint64_t link_index = (uint64_t)0, base_index = (uint64_t)0;
    int16_t mask = (int16_t)0xff;
    int32_t j;

    for (i = 0; i < coarse_db->seqs->size; i++) {
        uint64_t coarse_fasta_index, link_count = 0;
        char *fasta_output;

        fasta_output = malloc(30000*sizeof(*fasta_output));
        assert(fasta_output);

        seq = (struct cb_coarse_seq *)ds_vector_get(coarse_db->seqs, i);

        coarse_fasta_index = ftell(coarse_db->file_fasta);

        link_index = ftell(coarse_db->file_links);
        /*At the start of outputting each sequence, output the indices for the
          coarse links and FASTA files to their index files.*/
        fwrite(&link_index, sizeof(link_index), 1, coarse_db->file_links_index);
        fwrite(&base_index, sizeof(base_index), 1,
               coarse_db->file_fasta_base_index);
        fwrite(&coarse_fasta_index, sizeof(coarse_fasta_index), 1,
               coarse_db->file_fasta_index);

        /*Output the FASTA sequence to the coarse FASTA file*/
        sprintf(fasta_output, ">%ld\n%s\n", i, seq->seq->residues);
        fputs(fasta_output, coarse_db->file_fasta);
        free(fasta_output);

        /*Output all links for the current sequence to the coarse links file*/
        for (link = seq->links; link != NULL; link = link->next) {
            uint64_t start_base        = base_index +
                                           (uint64_t)link->data->coarse_start,
                     end_base          = base_index +
                                           (uint64_t)link->data->coarse_end,
                     coarse_base_start = link->data->coarse_start + base_index,
                     coarse_base_end   = link->data->coarse_end + base_index;

            /*Convert the start and end indices for the link to two
              characters.*/
            int16_t coarse_start = (int16_t)link->data->coarse_start,
                    coarse_end   = (int16_t)link->data->coarse_end;

            fwrite(link->data, sizeof(link->data), 1, coarse_db->file_links);

            fwrite(&coarse_base_start, sizeof(coarse_base_start),
                   1, coarse_db->file_links_base_index);
            fwrite(&coarse_base_end, sizeof(coarse_base_end),
                   1, coarse_db->file_links_base_index);

            link_count++;
        }

        fwrite(&link_count, sizeof(link_count), 1,
               coarse_db->file_links_count_index);
        base_index += strlen(seq->seq->residues);
    }

    fwrite(&base_index, sizeof(base_index), 1,
           coarse_db->file_fasta_base_index);

    output_int_to_file(coarse_db->dbsize, 8, coarse_db->file_params);
    putc('\n', coarse_db->file_params);
}

/*Output the coarse database's sequences and links in a human-readable format.
  Used for debugging purposes.*/
void cb_coarse_save_plain(struct cb_coarse *coarse_db){
    struct cb_coarse_seq *seq;
    struct cb_link_to_compressed *link;
    int32_t i;

    for (i = 0; i < coarse_db->seqs->size; i++) {
        seq = (struct cb_coarse_seq *)ds_vector_get(coarse_db->seqs, i);
        fprintf(coarse_db->file_fasta, "> %d\n%s\n", i, seq->seq->residues);

        fprintf(coarse_db->file_links, "> %d\n", i);
        for (link = seq->links; link != NULL; link = link->next) {
            fprintf(coarse_db->file_links,
              "    original sequence id: %d, reference range: (%d, %d), "
              "direction: %c\n",
              link->data->org_seq_id, link->data->coarse_start, link->data->coarse_end,
              (link->data->dir?'0':'1'));
        }
    }
}

/*Takes in the coarse database and saves its seeds table in a binary format.*/
void cb_coarse_save_seeds_binary(struct cb_coarse *coarse_db){
    struct cb_coarse_seq *seq;
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
                if (loc)
                    putc((char)0, coarse_db->file_seeds);
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
void cb_coarse_save_seeds_plain(struct cb_coarse *coarse_db){
    struct cb_coarse_seq *seq;
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

        for (; loc; loc = loc->next)
            if (loc->coarse_seq_id < 500)
                fprintf(coarse_db->file_seeds,"(%d, %d) > ",
                        loc->coarse_seq_id, loc->residue_index);

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
      fseek(coarse_db->db->file_links_index, 0, SEEK_END) == 0;
    if (!fseek_success)
        fprintf(stderr, "Error in seeking to end of FASTA index file\n");
    num_link_vectors = ftell(coarse_db->db->file_links_index)/8;
    fseek(coarse_db->db->file_links_index, 0, SEEK_SET);

    coarse_db->links = ds_vector_create();

    for (i = 0; i < num_link_vectors; i++) {
        struct DSVector *links =
          get_coarse_sequence_links(coarse_db->db->file_links);
        for (j = 0; j < links->size; j++)
            ds_vector_append(coarse_db->links, ds_vector_get(links, j));
        ds_vector_free_no_data(links);
    }
}

/*Loads all of the residues in the coarse database's FASTA file into the coarse
  database's all_residues string*/
void cb_coarse_get_all_residues(struct cb_coarse_db_read *coarse_db){
    int64_t num_bases = 0, i = 0;
    int32_t num_fasta_entries = 0, bases_copied = 0, j = 0;
    char *line = NULL;

    /*Get the number of entries in the coarse FASTA file*/
    bool fseek_success = fseek(coarse_db->db->file_fasta_index,0,SEEK_END) == 0;
    if (!fseek_success)
        fprintf(stderr, "Error in seeking to end of FASTA index file\n");
    num_fasta_entries = ftell(coarse_db->db->file_fasta_index)/8;
    fseek(coarse_db->db->file_fasta_index, 0, SEEK_SET);

    /*Get the number of bases in the coarse FASTA file*/
    fseek_success =
      fseek(coarse_db->db->file_fasta_base_index, -8, SEEK_END) == 0;
    if (!fseek_success)
        fprintf(stderr, "Error in seeking to end of FASTA index file\n");
    fread(&num_bases, sizeof(num_bases), 1,
          coarse_db->db->file_fasta_base_index);
    fseek(coarse_db->db->file_fasta_base_index, 0, SEEK_SET);

    coarse_db->all_residues =
      malloc((num_bases+1)*sizeof(*(coarse_db->all_residues)));
    assert(coarse_db->all_residues);

    for (i = 0; i < num_fasta_entries; i++) {
        int32_t line_length = coarse_db->seq_base_indices[i+1] -
                              coarse_db->seq_base_indices[i];
        FILE *fasta = coarse_db->db->file_fasta;
        if (0 == readline(fasta, &line)) {
            free(line);
            return;
        }
        else
            free(line);
        line = malloc((line_length+1)*sizeof(*line));
        fread(line, sizeof(*line), line_length+1, coarse_db->db->file_fasta);
        /*strcat(coarse_db->all_residues, line);*/
        for (j = 0; line[j] != '\0' && line[j] != '\n'; j++)
            coarse_db->all_residues[bases_copied++] = line[j];
        free(line);
    }

    coarse_db->all_residues[num_bases] = '\0';
}

/*Takes in an ID number and the residues and original start and end indices
 *for a sequence to be added to a coarse database and creates a coarse
 *sequence from that information.
 */
struct cb_coarse_seq *cb_coarse_seq_init(int32_t id, char *residues,
                                         int32_t start, int32_t end){
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
void cb_coarse_seq_free(struct cb_coarse_seq *seq){
    int32_t errno;

    if (0 != (errno = pthread_rwlock_destroy(&seq->lock_links))) {
        fprintf(stderr, "Could not destroy rwlock. Errno: %d\n", errno);
        exit(1);
    }
    cb_link_to_compressed_free(seq->links);

    cb_seq_free(seq->seq);
    free(seq);
}

/*Takes in a coarse sequence and a link_to_compressed and adds the link to the
  sequence at the end of the linked list of links.*/
void cb_coarse_seq_addlink(struct cb_coarse_seq *seq,
                           struct cb_link_to_compressed *newlink){
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

    link->data = malloc(sizeof(*(link->data)));
    assert(link->data);

    fread(&(link->data->org_seq_id), sizeof(link->data->org_seq_id), 1, f);
    fread(&(link->data->coarse_start), sizeof(link->data->coarse_start), 1, f);
    fread(&(link->data->coarse_end), sizeof(link->data->coarse_end), 1, f);
    fread(&(link->data->original_start), sizeof(link->data->original_start), 1, f);
    fread(&(link->data->original_end), sizeof(link->data->original_end), 1, f);
    link->data->dir = getc(f) == '0';
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
        if (feof(f) || c == '#')
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
    int64_t mask = make_mask(8), offset = (int64_t)(-1), try_off = id * 8;
    int32_t i;
    bool fseek_success = fseek(index_file, try_off, SEEK_SET) == 0;

    if (!fseek_success) {
        fprintf(stderr, "Error in seeking to offset %ld\n", try_off);
        return (int64_t)(-1);
    }
    fread(&offset, sizeof(offset), 1, index_file);

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
        fprintf(stderr, "Error in seeking to offset %ld\n", offset);
        return NULL;
    }
    return fasta_read_next(coarsedb->file_fasta, "");
}

/*Loads a cb_coarse_db_read struct, which includes a coarse database as well
  as data from the files being read*/
struct cb_coarse_db_read *
cb_coarse_read_init(int32_t seed_size,
                    FILE *file_fasta, FILE *file_seeds, FILE *file_links,
                    FILE *file_links_index, FILE *file_links_base_index,
                    FILE *file_links_count_index, FILE *file_fasta_index,
                    FILE *file_fasta_base_index, FILE *file_params,
                    bool load_coarse_residues, bool load_coarse_links,
                    int32_t link_block_size){
    int64_t link_count = (int64_t)0, links_in_sequence;
    int32_t i;
    bool fseek_success;

    struct cb_coarse_db_read *coarsedb = malloc(sizeof(*coarsedb));
    assert(coarsedb);

    /*Initialize the main cb_coarse structure*/
    coarsedb->db = cb_coarse_init(seed_size, file_fasta, file_seeds,
                                  file_links, file_links_index,
                                  file_links_base_index, file_links_count_index,
                                  file_fasta_index, file_fasta_base_index,
                                  file_params, true);

    /*Get the number of sequences in the coarse database and get the base index
      for each sequence.*/
    fseek_success = (fseek(file_fasta_base_index, 0, SEEK_END)) == 0;
    if (!fseek_success)
        fprintf(stderr, "Error in seeking to end of FASTA base index file\n");
    coarsedb->num_coarse_seqs = ftell(file_fasta_base_index)/8 - 1;
    fseek(file_fasta_base_index, 0, SEEK_SET);

    /*Initialize the seq_link_counts and seq_link_indices arrays.*/
    coarsedb->seq_link_counts = malloc(coarsedb->num_coarse_seqs
                                       *sizeof(*(coarsedb->seq_link_counts)));
    assert(coarsedb->seq_link_counts);

    coarsedb->seq_base_indices = malloc((coarsedb->num_coarse_seqs+1)
                                        *sizeof(*(coarsedb->seq_base_indices)));
    assert(coarsedb->seq_base_indices);

    for (i = 0; i < coarsedb->num_coarse_seqs; i++) {
        fread(&links_in_sequence, sizeof(links_in_sequence),
              1, file_links_count_index);
        link_count += links_in_sequence;
        coarsedb->seq_link_counts[i] = link_count;
    }

    fread(coarsedb->seq_base_indices, sizeof(*(coarsedb->seq_base_indices)),
          (coarsedb->num_coarse_seqs+1), file_fasta_base_index);

    coarsedb->link_block_size = link_block_size;

    /*Get the blocks of link indices.*/
    cb_coarse_db_read_init_blocks(coarsedb);

    coarsedb->all_residues = NULL;
    coarsedb->links        = NULL;

    /*If the --load-coarse-residues or --load-coarse-db search flag was passed
      in, load the coarse residues into coarse_db->all_residues.*/
    if (load_coarse_residues)
        cb_coarse_get_all_residues(coarsedb);
    /*If the --load-coarse-links or --load-coarse-db search flag was passed in,
      load the coarse residues into coarse_db->all_residues.*/
    if (load_coarse_links)
        cb_coarse_get_all_links(coarsedb);

    return coarsedb;
}

/*Frees a cb_coarse_db_read struct.*/
void cb_coarse_db_read_free(struct cb_coarse_db_read *coarsedb){
    int32_t i, j;

    cb_coarse_free(coarsedb->db);

    for (i = 0; i < coarsedb->link_inds_by_block->size; i++)
        ds_vector_free(ds_vector_get(coarsedb->link_inds_by_block, i));
    ds_vector_free_no_data(coarsedb->link_inds_by_block);

    free(coarsedb->seq_base_indices);
    free(coarsedb->seq_link_counts);

    /*Only free the residues and links if they were allocated*/
    if (coarsedb->all_residues)
        free(coarsedb->all_residues);

    if (coarsedb->links)
        ds_vector_free(coarsedb->links);

    free(coarsedb);
}

/*Initializes a cb_coarse_db_read's link_inds_by_block to have an empty vector
  of link indices for each block of bases in the coarse FASTA file.*/
void cb_coarse_db_read_init_blocks(struct cb_coarse_db_read *coarse_db){
    FILE *file_fasta_base_index  = coarse_db->db->file_fasta_base_index,
         *file_links_base_index  = coarse_db->db->file_links_base_index,
         *file_links_count_index = coarse_db->db->file_links_count_index;
    int64_t num_link_blocks, current_link, *current_link_ptr = NULL;
    int32_t current_seq = 0, link_count = 0, i = 0,
            block_size = coarse_db->link_block_size;
    bool fseek_success;

    /*Initialize link_inds_by_block*/
    coarse_db->link_inds_by_block = ds_vector_create();

    /*Get the number of bases in the coarse FASTA file to determine the number
      of blocks of link indices to make.*/
    fseek_success = (fseek(file_fasta_base_index, -8, SEEK_END)) == 0;
    if (!fseek_success)
        fprintf(stderr, "Error in seeking to end of FASTA base index file\n");
    fread(&num_link_blocks, sizeof(num_link_blocks), 1, file_fasta_base_index);
    num_link_blocks = num_link_blocks / block_size + 1;
    fseek(file_fasta_base_index, 0, SEEK_SET);

    /*Initialize the blocks*/
    for (i = 0; i < num_link_blocks; i++)
        ds_vector_append(coarse_db->link_inds_by_block,
                         (void *)ds_vector_create());

    fseek(file_fasta_base_index, 0, SEEK_SET);

    current_link = 0;

    while (!feof(file_links_base_index)) {
        int64_t current_start, current_end;
        fread(&current_start, sizeof(current_start), 1, file_links_base_index);
        fread(&current_end, sizeof(current_end), 1, file_links_base_index);

        if (feof(file_links_base_index))
            break;

        current_start /= block_size;
        current_end /= block_size;

        for (i = current_start; i <= current_end; i++) {
            current_link_ptr = malloc(sizeof(current_link_ptr));
            assert(current_link_ptr);

            *current_link_ptr = current_link;

            ds_vector_append(ds_vector_get(coarse_db->link_inds_by_block, i),
                             current_link_ptr);
        }

        current_link++;
    }
}

struct DSVector *cb_coarse_get_block(struct cb_coarse_db_read *coarse_db,
                                     int32_t index){
    struct DSVector *links = ds_vector_create(),
                    *block = (struct DSVector *)ds_vector_get(
                               coarse_db->link_inds_by_block,
                               index);
    int32_t i;

    if (coarse_db->links != NULL)
        for (i = 0; i < block->size; i++) {
            int32_t link_index = *(int32_t *)ds_vector_get(block, i);
            struct cb_link_to_compressed *link =
              (struct cb_link_to_compressed *)ds_vector_get(coarse_db->links,
                                                            link_index);
            ds_vector_append(links, link);
        }
    return links;
}

char *cb_coarse_get_seq_residues(struct cb_coarse_db_read *coarse_db,
                                 int64_t id){
    int64_t start = coarse_db->seq_base_indices[id],
            end   = coarse_db->seq_base_indices[id+1], i;
    char *residues     = malloc((end-start+1)*sizeof(*residues)),
         *all_residues = coarse_db->all_residues;

    if (all_residues == NULL) {
        fprintf(stderr, "cb_coarse_get_seq_residues only works if "
                        "all_residues was initialized in "
                        "cb_coarse_read_init.\n");
        return NULL;
    }

    for (i = start; i < end; i++)
        residues[i-start] = all_residues[i];
    residues[end-start] = '\0';

    return residues;
}

/*Coarse database functions ending in _r are used on cb_coarse_db_read structs
 *and are used as wrapper functions for the regular coarse database functions
 *being called on the cb_coarse_db_read struct's coarsedb.
 */
struct cb_coarse_seq *cb_coarse_get_r(struct cb_coarse_db_read *coarse_db,
                                      int32_t i){
    return cb_coarse_get(coarse_db->db, i);
}
struct fasta_seq *cb_coarse_read_fasta_seq_r(struct cb_coarse_db_read *coarsedb,
                                             int64_t id){
    return cb_coarse_read_fasta_seq(coarsedb->db, id);
}
