#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "clibs/include/ds.h"

#include "bitpack.h"
#include "compressed.h"
#include "edit_scripts.h"

/*Takes in pointers to a compressed database file and its index file and a
 *bool telling whether or not to populate the database with the files passed in
 *and creates a new compressed database struct, reading in the compressed
 *sequences if true is passed into "populate".
 */
struct cb_compressed *cb_compressed_init(FILE *file_compressed,
                                         FILE *file_index, bool populate){
    struct cb_compressed *com_db = malloc(sizeof(*com_db));
    assert(com_db);

    com_db->file_compressed   = file_compressed;
    com_db->file_index        = file_index;
    com_db->seqs              = ds_vector_create();
    com_db->next_seq_to_write = 0;

    if (populate) {
        int64_t num_sequences;
        int32_t i;

        bool fseek_success = fseek(file_index, 0, SEEK_END) == 0;
        if (!fseek_success) {
            fprintf(stderr, "Error in seeking to end of compressed.cb.index\n");
            return NULL;
        }
        num_sequences = ftell(file_index)/8;
        fseek_success = fseek(file_index, 0, SEEK_SET) == 0;
        if (!fseek_success) {
            fprintf(stderr,
                    "Error in seeking to start of compressed.cb.index\n");
            return NULL;
        }

        for (i = 0; i < num_sequences; i++)
            ds_vector_append(com_db->seqs, cb_compressed_read_seq_at(com_db,i));
    }

    return com_db;
}

//Frees the compressed database and its sequences and closes its files.
void cb_compressed_free(struct cb_compressed *com_db){
    fclose(com_db->file_compressed);
    fclose(com_db->file_index);

    for (int i = 0; i < com_db->seqs->size; i++)
        cb_compressed_seq_free(cb_compressed_seq_at(com_db, i));

    ds_vector_free_no_data(com_db->seqs);

    free(com_db);
}

//Returns the number of sequences in the compressed database
int32_t cb_compressed_size(struct cb_compressed *com_db){
    return com_db->seqs->size;
}

//Add a sequence to the compressed database
void cb_compressed_add(struct cb_compressed *com_db,
                       struct cb_compressed_seq *seq){
    ds_vector_append(com_db->seqs, (void *)seq);
}

/*Outputs a compressed sequence in the compressed database to the database's
  compressed file in human-readable text.  Used for debugging.*/
void cb_compressed_write(struct cb_compressed *com_db,
                         struct cb_compressed_seq *seq){
    struct cb_link_to_coarse *link;

    fprintf(com_db->file_compressed, "> %ld; %s\n", seq->id, seq->name);
    for (link = seq->links; link != NULL; link = link->next)
        fprintf(com_db->file_compressed,
                "reference sequence id: %ld, reference range: (%d, %d)\n%s\n",
                link->coarse_seq_id, link->coarse_start, link->coarse_end,
                link->diff);
}

/*Outputs a compressed sequence in the compressed database to the database's
  compressed file in binary format.*/
uint64_t cb_compressed_write_binary(struct cb_compressed *com_db,
                                    struct cb_compressed_seq *seq){
    struct cb_link_to_coarse *link, *find_length;
    uint64_t index = ftell(com_db->file_compressed), original_length = 0;
    char *id_string;

    id_string = malloc((20+strlen(seq->name))*sizeof(*id_string));
    assert(id_string);

    fwrite(&index, sizeof(index), 1, com_db->file_index);

    //Output the header for the sequence
    sprintf(id_string, ">%s\n", seq->name);
    fputs(id_string, com_db->file_compressed);

    free(id_string);

    find_length = seq->links;
    for (; find_length; find_length = find_length->next)
        original_length = find_length->original_end + 1;
    fwrite(&original_length,sizeof(original_length),1,com_db->file_compressed);

    for (link = seq->links; link != NULL; link = link->next){
        struct cb_link_to_coarse_data *link_data =
          cb_link_to_coarse_get_data(link);
        int32_t odd;
        uint16_t script_length = (uint16_t)0;
        char *edit_script = link->diff,
             *script      = edit_script_to_half_bytes(edit_script);

        while (edit_script[script_length] != '\0')
            script_length++;
        odd = script_length % 2 == 1 ? 1 : 0;

        //Output link_data in binary format
        fwrite(link_data, sizeof(*link_data), 1, com_db->file_compressed);

        //Output all of the characters of the edit script as half-bytes
        fwrite(script, sizeof(*script),
               script_length/2+odd, com_db->file_compressed);

        /*If there are more links for this sequence, the character after
         *the edit script is a space.  Otherwise, the character after the
         *edit script is the > of the FASTA header.  If we are printing the
         *last link of the last sequence, print a newline.
         */
        if (link->next)
            putc(' ', com_db->file_compressed);

        free(link_data);
        free(script);
    }
    putc('\n', com_db->file_compressed);
    return seq->id;
}

/*Takes in the ID number and the name of a sequence and creates a new compressed
  sequence with its list of links initialized to NULL.*/
struct cb_compressed_seq *cb_compressed_seq_init(int32_t id, char *name){
    struct cb_compressed_seq *seq = malloc(sizeof(*seq));
    assert(seq);

    seq->id        = id;
    seq->links     = NULL;
    seq->last_link = NULL;

    seq->name      = malloc((1+strlen(name))*sizeof(*seq->name));
    assert(seq->name);

    strcpy(seq->name, name);

    return seq;
}

//Frees a compressed sequence
void cb_compressed_seq_free(struct cb_compressed_seq *seq){
    cb_link_to_coarse_free(seq->links);
    free(seq->name);
    free(seq);
}

/*Takes in a compressed sequence and a cb_link_to_coarse and adds the link to
  the tail of the sequence's list of links*/
void cb_compressed_seq_addlink(struct cb_compressed_seq *seq,
                               struct cb_link_to_coarse *newlink){
    assert(newlink->next == NULL);

    if (seq->links == NULL) {
        seq->links     = newlink;
        seq->last_link = newlink;
        return;
    }

    seq->last_link->next = newlink;
    seq->last_link       = newlink;
}

//Gets the i-th compressed sequence from the compressed database.
struct cb_compressed_seq *cb_compressed_seq_at(struct cb_compressed *com_db,
                                               int32_t i){
    return (struct cb_compressed_seq *)ds_vector_get(com_db->seqs, i);
}

/*A function for getting the header for an entry in the compressed links file.
  Returns NULL if EOF is found before a newline.*/
char *get_compressed_header(FILE *f){
    int header_length = 10000, chars_read = 0 , c = 0 , i = 0;
    char *header;

    header = malloc(10000*sizeof(*header));
    assert(header);

    while (c != EOF && (char)c != '\n') {
        chars_read++;
        c = getc(f);
        if (c != EOF && (char)c != '\n' && chars_read > 1) {
            header[i] = (char)c;
            i++;
            if (i == header_length-1) {
                header_length *= 2;
                header = realloc(header, header_length*sizeof(*header));
                assert(header);
            }
        }
        if (c == EOF) {
            free(header);
            return NULL;
        }
    }
    header[i] = '\0';

    header = realloc(header, (i+1)*sizeof(*header));
    assert(header);

    return header;
}

/*Reads one link from a file with the links to the coarse database and converts
  its data to a struct cb_link_to_coarse*/
struct cb_link_to_coarse *read_compressed_link(FILE *f){
    struct cb_link_to_coarse_data *link_data;
    int32_t chars_to_read;
    uint16_t script_length = (uint16_t)0;
    char *half_bytes, *diff;
    bool fread_success;

    //Read in the data of the current link
    link_data = malloc(sizeof(*link_data));
    assert(link_data);

    fread_success = fread(link_data, sizeof(*link_data), 1, f) == 1;
    assert(fread_success);

    //Read in the link's edit script and convert it to ASCII format
    script_length = link_data->script_length;

    chars_to_read = script_length / 2;
    if (script_length % 2 == 1)
        chars_to_read++;

    half_bytes = malloc(chars_to_read*sizeof(*half_bytes));
    assert(half_bytes);

    fread_success =
      fread(half_bytes, sizeof(*half_bytes), chars_to_read, f) == chars_to_read;
    assert(fread_success);

    diff = half_bytes_to_ASCII(half_bytes, script_length);

    free(half_bytes);

    //Create a new cb_link_to_coarse with the link's data and edit script.
    struct cb_link_to_coarse *link =
      cb_link_to_coarse_from_data(link_data, diff);
    free(link_data);

    return link;
}

/*Takes in a pointer to the coarse.links file generated by cablast-compress and
 *the ID number of the sequence and returns the compressed sequence that the
 *data being read represents.  For this function to work properly the file
 *pointer must be pointing to the start of the header of a sequence entry in the
 *links file.
 */
struct cb_compressed_seq *get_compressed_seq(FILE *f, int id){
    struct cb_link_to_coarse *first_link = NULL, *last_link = NULL;
    struct cb_compressed_seq *seq;
    int64_t l;
    char *h = get_compressed_header(f);
    bool fread_success;

    seq = cb_compressed_seq_init(id, h);

    if (h == NULL) {
        fprintf(stderr, "Could not get compressed sequence\n");
        return NULL;
    }

    free(h);

    fread_success = fread(&l, sizeof(l), 1, f) == 1;
    assert(fread_success);

    while (true) {
        char c = 1;
        struct cb_link_to_coarse *current_link = read_compressed_link(f);

        if (current_link == NULL)
            break;

        if (!first_link)
            first_link = current_link;
        else
            last_link->next = current_link;
        last_link = current_link;
        c = getc(f);
        if (c == '\n')
            break;
    }
    seq->links = first_link;

    return seq;
}

/*Takes in a pointer to the coarse.links file generated by cablast-compress and
 *the ID number of the sequence and returns the compressed sequence that the
 *data being read represents.  For this function to work properly the file
 *pointer must be pointing to the start of the header of a sequence entry in the
 *links file.
 */
struct cb_compressed_seq *cb_compressed_read_seq_at(struct cb_compressed *comdb,
                                                    int32_t id){
    FILE *links = comdb->file_compressed;
    int64_t offset = cb_compressed_link_offset(comdb, id);
    bool fseek_success;

    if (offset < 0)
        return NULL;

    fseek_success = fseek(links, offset, SEEK_SET) == 0;
    if (!fseek_success) { 
        fprintf(stderr, "Error in seeking to offset %lu\n", offset);
        return NULL;
    }

    return get_compressed_seq(links, id);
}

/*Takes in a file pointer to a compressed database's compressed.cb file and
 *returns the length of the original sequence that was compressed to produce
 *the sequence that the file pointer is pointing to.  For this function to work
 *properly, f must be pointing to the start of a compressed sequence in the
 *compressed.cb file.
 */
int64_t cb_compressed_get_seq_length(FILE *f){
    int64_t length;
    char *h = get_compressed_header(f);
    bool fread_success;

    if (h == NULL) {
        fprintf(stderr, "Could not get sequence length\n");
        return -1;
    }
    free(h);

    fread_success = fread(&length, sizeof(length), 1, f) == 1;
    assert(fread_success);

    return length;
}

//Gets the lengths in bases for all sequences in the database
int64_t *cb_compressed_get_lengths(struct cb_compressed *comdb){
    FILE *links = comdb->file_compressed, *index = comdb->file_index;
    int64_t *lengths = NULL, num_sequences;
    bool fseek_success;

    fseek_success = fseek(index, 0, SEEK_END) == 0;
    if (!fseek_success) {
        fprintf(stderr, "Error in seeking to end of compressed.cb.index\n");
        return NULL;
    }

    num_sequences = ftell(index) / 8;

    fseek_success = fseek(index, 0, SEEK_SET) == 0;
    if (!fseek_success) {
        fprintf(stderr, "Error in seeking to start of compressed.cb.index\n");
        return NULL;
    }

    if (num_sequences > 0) {
        int32_t i = 0;

        lengths = malloc(num_sequences*sizeof(*lengths));
        assert(lengths);

        for (i = 0; i < num_sequences; i++) {
            int64_t offset = cb_compressed_link_offset(comdb, i);

            fseek_success = fseek(links, offset, SEEK_SET) == 0;
            if (!fseek_success) {
                fprintf(stderr, "error in seeking to offset %ld\n", offset);
                free(lengths);
                return NULL;
            }

            lengths[i] = cb_compressed_get_seq_length(links);
        }
    }

    fseek_success = fseek(links, 0, SEEK_SET) == 0;
    if (!fseek_success) {
        fprintf(stderr, "error in seeking to start of compressed.cb.index");
        free(lengths);
        return NULL;
    }

    fseek_success = fseek(index, 0, SEEK_SET) == 0;
    if (!fseek_success) {
        fprintf(stderr, "error in seeking to start of compressed.cb");
        free(lengths);
        return NULL;
    }

    return lengths;
}

/*Takes the compressed.cb generated by cablast-compress and parses it to get
  an array of compressed sequences.*/
struct cb_compressed_seq **read_compressed(FILE *f){
    int length = 0, capacity = 1000;

    struct cb_compressed_seq **compressed_seqs =
      malloc(1000*sizeof(*compressed_seqs));
    assert(compressed_seqs);

    //Read each sequence
    while (true) {
        struct cb_link_to_coarse *links = NULL;
        char *header = get_compressed_header(f);

        if (header == NULL)
            break;

        //Skip over the length of the original sequence
        read_int_from_file(8, f);

        //Read each link in the sequence
        while (true) {
            struct cb_link_to_coarse *current_link = read_compressed_link(f);
            char c = 1;

            if (current_link == NULL)
                break;
            current_link->next = links;
            if (current_link->next == NULL)
                links = current_link;
            else {
                while (current_link->next->next != NULL)
                    current_link->next = current_link->next->next;
                current_link->next->next = current_link;
                current_link->next = NULL;
            }
            c = getc(f);

            /*If we find a newline at the end of the link we are at the end of
              the sequence.*/
            if (c == '\n')
                break;
        }

        compressed_seqs[length] = malloc(sizeof(*(compressed_seqs[length])));
        assert(compressed_seqs[length]);

        compressed_seqs[length]->name  = header;
        compressed_seqs[length]->links = links;
        length++;

        if (length == capacity) {
            capacity *= 2;
            compressed_seqs =
              realloc(compressed_seqs, capacity*sizeof(*compressed_seqs));
            assert(compressed_seqs);
        }
    }

    compressed_seqs =
      realloc(compressed_seqs, (length+1)*sizeof(*compressed_seqs));
    assert(compressed_seqs);

    compressed_seqs[length] = NULL;
    return compressed_seqs;
}

/*Gets the offset in the compressed database's compressed.cb file for the
  link whose index is passed into id.*/
int64_t cb_compressed_link_offset(struct cb_compressed *comdb, int id){
    int64_t offset;
    int32_t try_off = id * 8;
    bool fread_success,
         fseek_success = fseek(comdb->file_index, try_off, SEEK_SET) == 0;

    if (!fseek_success) {
        fprintf(stderr, "Error in seeking to offset %d", try_off);
        return (int64_t)(-1);
    }

    fread_success = fread(&offset, sizeof(offset), 1, comdb->file_index) == 1;
    assert(fread_success);

    return offset;
}
