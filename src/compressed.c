#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "ds.h"

#include "bitpack.h"
#include "compressed.h"
#include "edit_scripts.h"

struct cb_compressed *cb_compressed_init(FILE *file_compressed,
                                         FILE *file_index, bool populate){
    struct cb_compressed *com_db;

    com_db = malloc(sizeof(*com_db));
    assert(com_db);

    com_db->file_compressed = file_compressed;
    com_db->file_index      = file_index;
    com_db->seqs            = ds_vector_create();

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

        for (i = 0; i < num_sequences; i++) {
            ds_vector_append(com_db->seqs,
                             cb_compressed_read_seq_at(com_db, i));
        }
    }

    return com_db;
}

void cb_compressed_free(struct cb_compressed *com_db){
    int i;

    fclose(com_db->file_compressed);
    fclose(com_db->file_index);

    for (i = 0; i < com_db->seqs->size; i++)
        cb_compressed_seq_free(cb_compressed_seq_at(com_db, i));

    ds_vector_free_no_data(com_db->seqs);
    free(com_db);
}

int32_t cb_compressed_size(struct cb_compressed *com_db){
    return com_db->seqs->size;
}

void cb_compressed_add(struct cb_compressed *com_db,
                       struct cb_compressed_seq *seq){
    ds_vector_append(com_db->seqs, (void *)seq);
}

/*Takes in as input the compressed database and converts its sequences to a
 *binary output format, which is printed to the file pointed to by
 *com_db->file_compressed.
 */
void cb_compressed_save_binary(struct cb_compressed *com_db){
    struct cb_compressed_seq *seq;
    struct cb_link_to_coarse *link;

    uint64_t index = (uint64_t)0;
    int32_t i;
    int16_t mask = (((int16_t)1)<<8)-1;

    for (i = 0; i < com_db->seqs->size; i++) {
        struct cb_link_to_coarse *find_length;
        uint64_t original_length = 0;
        int32_t j;
        char *id_bytes;

        output_int_to_file(index, 8, com_db->file_index);
        seq = cb_compressed_seq_at(com_db, i);
        id_bytes = read_int_to_bytes(seq->id, 8);

        putc('>', com_db->file_compressed);
        putc(' ', com_db->file_compressed);
        for (j = 0; j < 8; j++)
            putc(id_bytes[j], com_db->file_compressed);
        putc(';', com_db->file_compressed);
        putc(' ', com_db->file_compressed);

        index += 12;

        for(j = 0; seq->name[j] != '\0'; j++)
            putc(seq->name[j], com_db->file_compressed);
        index += j;
        putc('\n', com_db->file_compressed);
        index++;

        find_length = seq->links;
        for (; find_length; find_length = find_length->next)
            original_length = find_length->original_end + 1;

        output_int_to_file(original_length, 8, com_db->file_compressed);

        for (link = seq->links; link != NULL; link = link->next){
            /*Convert the start and end indices for the link to two
              characters.*/
            uint64_t org_start = link->original_start,
                     org_end   = link->original_end,
                     coarse_seq_id;
            int16_t cor_start = (int16_t)link->coarse_start,
                    cor_end   = (int16_t)link->coarse_end;
            char cor_start_left  = (cor_start >> 8) & mask,
                 cor_start_right = cor_start & mask,
                 cor_end_left    = (cor_end >> 8) & mask,
                 cor_end_right   = cor_end & mask;

            char script_left, script_right,
                 *script, *edit_script = link->diff;
            int16_t script_length = (int16_t)0;
            int32_t j;

            /*Output the ID of the current chunk as 8 characters*/
            coarse_seq_id = link->coarse_seq_id;
            for (j = 7; j >= 0; j--) {
                char b = shift_right(coarse_seq_id, j*8) & mask;
                putc(b, com_db->file_compressed);
            }
            index += 8;

            /*Represent the length of the edit script as two characters and get
              the edit script as a sequence of half-bytes*/
            while (edit_script[script_length] != '\0')
                script_length++;

            script_left  = (script_length >> 8) & mask;
            script_right = script_length & mask;
            script = edit_script_to_half_bytes(edit_script);

            /*Output the indices of the start and end of the sequence being
             *linked to and the length of the edit script represented in 16
             *bits, and the length of the edit script.
             */
            output_int_to_file(org_start, 8, com_db->file_compressed);
            output_int_to_file(org_end, 8, com_db->file_compressed);
            putc(cor_start_left, com_db->file_compressed);
            putc(cor_start_right, com_db->file_compressed);
            putc(cor_end_left, com_db->file_compressed);
            putc(cor_end_right, com_db->file_compressed);
            putc(script_left, com_db->file_compressed);
            putc(script_right, com_db->file_compressed);
            index += 6;

            /*Output all of the characters of the edit script as half-bytes*/
            for (j = 0; j < script_length/2; j++)
                putc(script[j], com_db->file_compressed);
            index += j;

            /*If there are more links for this sequence, the character after
             *the edit script is a space.  Otherwise, the character after the
             *edit script is the > of the FASTA header.  If we are printing the
             *last link of the last sequence, print a newline.
             */
            if (link->next) {
                putc(' ', com_db->file_compressed);
                index++;
            }
            else
                if (i == com_db->seqs->size-1) {
                    putc('\n', com_db->file_compressed);
                    index++;
                }
        }
    }
}

void cb_compressed_save_plain(struct cb_compressed *com_db){
    struct cb_compressed_seq *seq;
    struct cb_link_to_coarse *link;
    int i;

    for (i = 0; i < com_db->seqs->size; i++) {
        seq = cb_compressed_seq_at(com_db, i);
        fprintf(com_db->file_compressed, "> %ld; %s\n", seq->id, seq->name);
        for (link = seq->links; link != NULL; link = link->next)
            fprintf(com_db->file_compressed,
              "reference sequence id: %ld, reference range: (%d, %d), "
              "original sequence range: (%ld %ld)\n%s\n",
              link->coarse_seq_id, link->coarse_start, link->coarse_end,
              link->original_start, link->original_end,
              link->diff);
    }
}

void cb_compressed_write(struct cb_compressed *com_db,
                         struct cb_compressed_seq *seq){
    struct cb_link_to_coarse *link;

    fprintf(com_db->file_compressed, "> %ld; %s\n", seq->id, seq->name);
    for (link = seq->links; link != NULL; link = link->next) {
        fprintf(com_db->file_compressed,
          "reference sequence id: %ld, reference range: (%d, %d)\n%s\n",
          link->coarse_seq_id, link->coarse_start, link->coarse_end,
          link->diff);
    }
}

/*Outputs a compressed sequence in the compressed database to the database's
  compressed file in binary format.*/
void cb_compressed_write_binary(struct cb_compressed *com_db,
                                struct cb_compressed_seq *seq){
    struct cb_link_to_coarse *link, *find_length;
    uint64_t index = ftell(com_db->file_compressed), original_length = 0;
    int32_t i;
    int16_t mask = (((int16_t)1)<<8)-1;
    char *id_string;

    id_string = malloc((20+strlen(seq->name))*sizeof(*id_string));
    assert(id_string);

    fwrite(&index, sizeof(index), 1, com_db->file_index);
    /*output_int_to_file(index, 8, com_db->file_index);*/

    sprintf(id_string, "> %ld; %s\n", seq->id, seq->name);

    /*Output the header for the sequence*/
    fputs(id_string, com_db->file_compressed);
    /*putc('>', com_db->file_compressed);
    putc(' ', com_db->file_compressed);
    for (i = 0; id_string[i] != '\0'; i++)
        putc(id_string[i], com_db->file_compressed);
    putc(';', com_db->file_compressed);
    putc(' ', com_db->file_compressed);
    for (i = 0; seq->name[i] != '\0'; i++)
        putc(seq->name[i], com_db->file_compressed);
    putc('\n', com_db->file_compressed);*/

    free(id_string);

    find_length = seq->links;
    for (; find_length; find_length = find_length->next)
        original_length = find_length->original_end + 1;

    /*output_int_to_file(original_length, 8, com_db->file_compressed);*/
    fwrite(&original_length,sizeof(original_length),1,com_db->file_compressed);

    for (link = seq->links; link != NULL; link = link->next){
        /*Convert the start and end indices for the link to two characters.*/
        uint64_t org_start = link->original_start,
                 org_end   = link->original_end,
                 coarse_seq_id;
        int32_t odd;
        uint16_t cor_start     = link->coarse_start,
                 cor_end       = link->coarse_end,
                 script_length = (uint16_t)0;
        char cor_start_left  = (cor_start >> 8) & mask,
             cor_start_right = cor_start & mask,
             cor_end_left    = (cor_end >> 8) & mask,
             cor_end_right   = cor_end & mask;

        char *edit_script = link->diff,
             *script = edit_script_to_half_bytes(edit_script),
             script_left, script_right;

        /*Output the ID of the current chunk as 8 characters*/
        coarse_seq_id = link->coarse_seq_id;
fprintf(stderr, "coarse_seq_id: %lu\n", coarse_seq_id);
        fwrite(&coarse_seq_id,sizeof(coarse_seq_id),1,com_db->file_compressed);
        /*for (i = 7; i >= 0; i--) {
            char b = (char)(shift_right(coarse_seq_id, i*8) & mask);
            putc(b, com_db->file_compressed);
        }*/

        /*Represent the length of the edit script as two characters and get
          the edit script as a sequence of half-bytes*/
        while (edit_script[script_length] != '\0')
            script_length++;
        odd = script_length % 2 == 1 ? 1 : 0;
        script_left  = (script_length >> 8) & mask;
        script_right = script_length & mask;

        /*Output the indices of the start and end of the sequence being
         *linked to and the length of the edit script represented in 16
         *bits, and the length of the edit script.
         */
        /*output_int_to_file(org_start, 8, com_db->file_compressed);*/
        fwrite(&org_start, sizeof(org_start), 1, com_db->file_compressed);
        /*output_int_to_file(org_end, 8, com_db->file_compressed);*/
        fwrite(&org_end, sizeof(org_end), 1, com_db->file_compressed);
        /*putc(cor_start_left, com_db->file_compressed);
        putc(cor_start_right, com_db->file_compressed);*/
        fwrite(&cor_start, sizeof(cor_start), 1, com_db->file_compressed);
        /*putc(cor_end_left, com_db->file_compressed);
        putc(cor_end_right, com_db->file_compressed);*/
        fwrite(&cor_end, sizeof(cor_end), 1, com_db->file_compressed);
        /*putc(script_left, com_db->file_compressed);
        putc(script_right, com_db->file_compressed);*/
        fwrite(&script_length,sizeof(script_length),1,com_db->file_compressed);

        /*Output all of the characters of the edit script as half-bytes*/
        fwrite(script, sizeof(*script),
               script_length/2+odd, com_db->file_compressed);
        /*for (i = 0; i < script_length/2+odd; i++)
            putc(script[i], com_db->file_compressed);*/

        /*If there are more links for this sequence, the character after
         *the edit script is a space.  Otherwise, the character after the
         *edit script is the > of the FASTA header.  If we are printing the
         *last link of the last sequence, print a newline.
         */
        if (link->next)
            putc(' ', com_db->file_compressed);

        free(script);
    }
    putc('\n', com_db->file_compressed);
}

struct cb_compressed_seq *cb_compressed_seq_init(int32_t id, char *name){
    struct cb_compressed_seq *seq = malloc(sizeof(*seq));
    assert(seq);

    seq->id = id;
    seq->links = NULL;

    seq->name = malloc((1+strlen(name))*sizeof(*seq->name));
    assert(seq->name);

    strcpy(seq->name, name);

    return seq;
}

void cb_compressed_seq_free(struct cb_compressed_seq *seq){
    cb_link_to_coarse_free(seq->links);
    free(seq->name);
    free(seq);
}

void cb_compressed_seq_addlink(struct cb_compressed_seq *seq,
                               struct cb_link_to_coarse *newlink){
    struct cb_link_to_coarse *link;

    assert(newlink->next == NULL);

    if (seq->links == NULL) {
        seq->links = newlink;
        return;
    }

    for (link = seq->links; link->next != NULL; link = link->next);
    link->next = newlink;
}

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
        if (c != EOF && (char)c != '\n' && chars_read > 2) {
            header[i] = (char)c;
            i++;
            if (i == header_length-1) {
                header_length *= 2;
                header = realloc(header, header_length*sizeof(*header));
                assert(header);
            }
        }
        if (c == EOF)
            return NULL;
    }
    header[i] = '\0';

    header = realloc(header, (i+1)*sizeof(*header));
    assert(header);

    return header;
}

/*Reads one link from a file with the links to the coarse database and converts
  its data to a struct cb_link_to_coarse*/
struct cb_link_to_coarse *read_compressed_link(FILE *f){
    struct cb_link_to_coarse *link;
    int32_t chars_to_read;
    uint16_t script_length = (uint16_t)0;
    char *link_bytes, *half_bytes;

    link = malloc(sizeof(*link));
    assert(link);

    link_bytes = malloc(30*sizeof(link_bytes));
    assert(link_bytes);

    /*fread(link_bytes, 1, 30, f);
    if (feof(f)) {
        free(link);
        return NULL;
    }*/

    /*link->coarse_seq_id  = (uint64_t)bytes_to_int(link_bytes, 0, 8);*/
    fread(&(link->coarse_seq_id), sizeof(link->coarse_seq_id), 1, f);
printf("coarse_seq_id %ld\n", link->coarse_seq_id);
    /*link->original_start = (uint64_t)bytes_to_int(link_bytes, 8, 8);*/
    fread(&(link->original_start), sizeof(link->original_start), 1, f);
printf("original_start %ld\n", link->original_start);
    /*link->original_end   = (uint64_t)bytes_to_int(link_bytes, 16, 8);*/
    fread(&(link->original_end), sizeof(link->original_end), 1, f);
printf("original_end %ld\n", link->original_end);
    /*link->coarse_start   = (int16_t)bytes_to_int(link_bytes, 24, 2);*/
    fread(&(link->coarse_start), sizeof(link->coarse_start), 1, f);
printf("coarse_start %hd\n", link->coarse_start);
    /*link->coarse_end     = (int16_t)bytes_to_int(link_bytes, 26, 2);*/
    fread(&(link->coarse_end), sizeof(link->coarse_end), 1, f);
printf("coarse_end %hd\n", link->coarse_end);
    /*script_length        = (int16_t)bytes_to_int(link_bytes, 28, 2);*/
    fread(&script_length, sizeof(script_length), 1, f);
printf("script_length %hd\n\n", script_length);

    chars_to_read = script_length / 2;
    if (script_length % 2 == 1)
        chars_to_read++;

    half_bytes = malloc(chars_to_read*sizeof(*half_bytes));
    assert(half_bytes);

    fread(half_bytes, sizeof(*half_bytes), chars_to_read, f);

    link->diff = half_bytes_to_ASCII(half_bytes, script_length);
    link->next = NULL;

    free(half_bytes);

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

    seq = cb_compressed_seq_init(id, h);

    if (h == NULL) {
        fprintf(stderr, "Could not get compressed sequence\n");
        return NULL;
    }
    free(h);

    fread(&l, sizeof(l), 1, f);
fprintf(stderr, "%ld\n", l);
    /*read_int_from_file(8, f);*/

    while (true) {
        struct cb_link_to_coarse *current_link = read_compressed_link(f);
        char c = 1;
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
fprintf(stderr, "*\n");
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
    if (h == NULL) {
        fprintf(stderr, "Could not get sequence length\n");
        return -1;
    }
    free(h);

    fread(&length, sizeof(length), 1, f);
    return length;
    /*return read_int_from_file(8, f);*/
}

/*Gets the lengths in bases for all sequences in the database*/
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

    fseek(links, 0, SEEK_SET);
    fseek(index, 0, SEEK_SET);
    return lengths;
}


/*Takes the compressed.cb generated by cablast-compress and parses it to get
  an array of compressed sequences.*/
struct cb_compressed_seq **read_compressed(FILE *f){
    int length = 0;

    struct cb_compressed_seq **compressed_seqs =
      malloc(1000*sizeof(*compressed_seqs));
    assert(compressed_seqs);

    /*Read each sequence*/
    while (true) {
        struct cb_link_to_coarse *links = NULL;
        char *header = get_compressed_header(f);

        if (header == NULL)
            break;

        /*Skip over the length of the original sequence*/
        read_int_from_file(8, f);

        /*Read each link in the sequence*/
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

        (compressed_seqs[length])->name = header;
        (compressed_seqs[length]->links) = links;
        length++;
    }

    compressed_seqs =
      realloc(compressed_seqs, (length+1)*sizeof(*compressed_seqs));
    assert(compressed_seqs);

    compressed_seqs[length] = NULL;
    return compressed_seqs;
}

/*Gets the offset in the compressed database's compressed.cb file for the
 *link whose index is passed into id.
 */
int64_t cb_compressed_link_offset(struct cb_compressed *comdb, int id){
    int64_t offset;
    int32_t try_off = id * 8;
    bool fseek_success = fseek(comdb->file_index, try_off, SEEK_SET) == 0;

    if (!fseek_success) {
        fprintf(stderr, "Error in seeking to offset %d", try_off);
        return (int64_t)(-1);
    }

    fread(&offset, sizeof(offset), 1, comdb->file_index);
    return offset;
}
