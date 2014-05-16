#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "ds.h"

#include "bitpack.h"
#include "compressed.h"
#include "edit_scripts.h"

struct cb_compressed *
cb_compressed_init(FILE *file_compressed, FILE *file_index)
{
    struct cb_compressed *com_db;

    com_db = malloc(sizeof(*com_db));
    assert(com_db);

    com_db->file_compressed = file_compressed;
    com_db->file_index = file_index;
    com_db->seqs = ds_vector_create_capacity(100);

    return com_db;
}

void
cb_compressed_free(struct cb_compressed *com_db)
{
    int i;

    fclose(com_db->file_compressed);
    fclose(com_db->file_index);

    for (i = 0; i < com_db->seqs->size; i++)
        cb_compressed_seq_free(cb_compressed_seq_at(com_db, i));

    ds_vector_free_no_data(com_db->seqs);
    free(com_db);
}

int32_t
cb_compressed_size(struct cb_compressed *com_db)
{
    return com_db->seqs->size;
}

void
cb_compressed_add(struct cb_compressed *com_db,
                   struct cb_compressed_seq *seq)
{
    ds_vector_append(com_db->seqs, (void*) seq);
}

/*Takes in as input the compressed database and converts its sequences to a
 *binary output format, which is printed to the file pointed to by
 *com_db->file_compressed.
 */
void
cb_compressed_save_binary(struct cb_compressed *com_db)
{
    int i;
    struct cb_compressed_seq *seq;
    struct cb_link_to_coarse *link;

    int16_t mask = (((int16_t)1)<<8)-1;
    uint64_t index = (uint64_t)0;

    for (i = 0; i < com_db->seqs->size; i++) {
        int j;
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

        struct cb_link_to_coarse *find_length = seq->links;
        uint64_t original_length = 0;
        for (; find_length; find_length = find_length->next)
            original_length = find_length->original_end + 1;

        output_int_to_file(original_length, 8, com_db->file_compressed);

        for (link = seq->links; link != NULL; link = link->next){
            /*Convert the start and end indices for the link to two
              characters.*/
            uint64_t org_start = link->original_start;
            uint64_t org_end   = link->original_end;
            int16_t cor_start = (int16_t)link->coarse_start;
            int16_t cor_end   = (int16_t)link->coarse_end;
            char cor_start_left  = (cor_start >> 8) & mask;
            char cor_start_right = cor_start & mask;
            char cor_end_left    = (cor_end >> 8) & mask;
            char cor_end_right   = cor_end & mask;

            char script_left, script_right;
            char *edit_script = link->diff;
            char *script;
            int16_t script_length = (int16_t)0;
            int j;

            /*Output the ID of the current chunk as 8 characters*/
            int coarse_seq_id = link->coarse_seq_id;
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

void
cb_compressed_save_plain(struct cb_compressed *com_db)
{
    int i;
    struct cb_compressed_seq *seq;
    struct cb_link_to_coarse *link;

    for (i = 0; i < com_db->seqs->size; i++) {
        seq = cb_compressed_seq_at(com_db, i);
        fprintf(com_db->file_compressed, "> %ld; %s\n", seq->id, seq->name);
        for (link = seq->links; link != NULL; link = link->next)
            fprintf(com_db->file_compressed,
                "reference sequence id: %d, reference range: (%d, %d), "
                  "original sequence range: (%ld %ld)\n%s\n",
                link->coarse_seq_id, link->coarse_start, link->coarse_end,
                link->original_start, link->original_end,
                link->diff);
    }
}

void
cb_compressed_write(struct cb_compressed *com_db,
                     struct cb_compressed_seq *seq)
{
    struct cb_link_to_coarse *link;

    fprintf(com_db->file_compressed, "> %ld; %s\n", seq->id, seq->name);
    for (link = seq->links; link != NULL; link = link->next) {
        fprintf(com_db->file_compressed,
            "reference sequence id: %d, reference range: (%d, %d)\n%s\n",
            link->coarse_seq_id, link->coarse_start, link->coarse_end,
            link->diff);
    }
}

/*Outputs a compressed sequence in the compressed database to the database's
  compressed file in binary format.*/
void
cb_compressed_write_binary(struct cb_compressed *com_db,
                            struct cb_compressed_seq *seq)
{
    int i;
    struct cb_link_to_coarse *link;
    int16_t mask = (((int16_t)1)<<8)-1;
    char *id_string = malloc(20*sizeof(*id_string));
    uint64_t index = ftell(com_db->file_compressed);
    output_int_to_file(index, 8, com_db->file_index);
    sprintf(id_string, "%ld", seq->id);

    /*Output the header for the sequence*/
    putc('>', com_db->file_compressed);
    putc(' ', com_db->file_compressed);
    for (i = 0; id_string[i] != '\0'; i++)
        putc(id_string[i], com_db->file_compressed);
    putc(';', com_db->file_compressed);
    putc(' ', com_db->file_compressed);
    for (i = 0; seq->name[i] != '\0'; i++)
        putc(seq->name[i], com_db->file_compressed);
    putc('\n', com_db->file_compressed);

    free(id_string);

    struct cb_link_to_coarse *find_length = seq->links;
    uint64_t original_length = 0;
    for (; find_length; find_length = find_length->next)
        original_length = find_length->original_end + 1;

    output_int_to_file(original_length, 8, com_db->file_compressed);

    for (link = seq->links; link != NULL; link = link->next){
        /*Convert the start and end indices for the link to two
          characters.*/
        uint64_t org_start = link->original_start;
        uint64_t org_end   = link->original_end;
        int16_t cor_start = (int16_t)link->coarse_start;
        int16_t cor_end   = (int16_t)link->coarse_end;
        char cor_start_left  = (cor_start >> 8) & mask;
        char cor_start_right = cor_start & mask;
        char cor_end_left    = (cor_end >> 8) & mask;
        char cor_end_right   = cor_end & mask;

        char script_left, script_right;
        char *edit_script = link->diff;
        char *script = edit_script_to_half_bytes(edit_script);
        int16_t script_length = (int16_t)0;
        int odd;

        /*Output the ID of the current chunk as 8 characters*/
        uint64_t coarse_seq_id = link->coarse_seq_id;
        for (i = 7; i >= 0; i--) {
            char b = (char)(shift_right(coarse_seq_id, i*8) & mask);
            putc(b, com_db->file_compressed);
        }

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
        output_int_to_file(org_start, 8, com_db->file_compressed);
        output_int_to_file(org_end, 8, com_db->file_compressed);
        putc(cor_start_left, com_db->file_compressed);
        putc(cor_start_right, com_db->file_compressed);
        putc(cor_end_left, com_db->file_compressed);
        putc(cor_end_right, com_db->file_compressed);
        putc(script_left, com_db->file_compressed);
        putc(script_right, com_db->file_compressed);

        /*Output all of the characters of the edit script as half-bytes*/
        for (i = 0; i < script_length/2+odd; i++)
            putc(script[i], com_db->file_compressed);

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

struct cb_compressed_seq *
cb_compressed_seq_init(int32_t id, char *name)
{
    struct cb_compressed_seq *seq;

    seq = malloc(sizeof(*seq));
    assert(seq);

    seq->id = id;
    seq->links = NULL;

    seq->name = malloc((1 + strlen(name)) * sizeof(*seq->name));
    assert(seq->name);
    strcpy(seq->name, name);

    return seq;
}

void
cb_compressed_seq_free(struct cb_compressed_seq *seq)
{
    struct cb_link_to_coarse *link1, *link2;

    for (link1 = seq->links; link1 != NULL; ) {
        link2 = link1->next;
        cb_link_to_coarse_free(link1);
        link1 = link2;
    }
    free(seq->name);
    free(seq);
}

void
cb_compressed_seq_addlink(struct cb_compressed_seq *seq,
                           struct cb_link_to_coarse *newlink)
{
    struct cb_link_to_coarse *link;

    assert(newlink->next == NULL);

    if (seq->links == NULL) {
        seq->links = newlink;
        return;
    }

    for (link = seq->links; link->next != NULL; link = link->next);
    link->next = newlink;
}


struct cb_link_to_coarse *
cb_link_to_coarse_init(int32_t coarse_seq_id,
                        uint64_t original_start, uint64_t original_end,
                        uint16_t coarse_start, uint16_t coarse_end,
                        struct cb_alignment alignment, bool dir){
    struct cb_link_to_coarse *link;

    link = malloc(sizeof(*link));
    assert(link);

    link->coarse_seq_id = coarse_seq_id;
    link->original_start = original_start;
    link->original_end = original_end;
    link->coarse_start = coarse_start;
    link->coarse_end = coarse_end;
    link->next = NULL;
    link->diff = make_edit_script(alignment.org, alignment.ref, dir,
                                                   alignment.length);
    assert(link->diff);

    return link;
}

struct cb_link_to_coarse *
cb_link_to_coarse_init_nodiff(int32_t coarse_seq_id,
                               uint64_t original_start, uint64_t original_end,
                               uint16_t coarse_start, uint16_t coarse_end,
                               bool dir){
    struct cb_link_to_coarse *link;

    link = malloc(sizeof(*link));
    assert(link);

    link->diff = malloc(2*sizeof(*(link->diff)));
    link->diff[0] = dir ? '0' : '1';
    link->diff[1] = '\0';
    link->coarse_seq_id = coarse_seq_id;
    link->original_start = original_start;
    link->original_end = original_end;
    link->coarse_start = coarse_start;
    link->coarse_end = coarse_end;
    link->next = NULL;

    return link;
}

void
cb_link_to_coarse_free(struct cb_link_to_coarse *link)
{
    free(link->diff);
    free(link);
}

struct cb_compressed_seq *
cb_compressed_seq_at(struct cb_compressed *com_db, int32_t i)
{
    return (struct cb_compressed_seq *) ds_vector_get(com_db->seqs, i);
}

/*A function for getting the header for an entry in the compressed links file.
  Returns NULL if EOF is found before a newline.*/
char *get_compressed_header(FILE *f){
    int c = 0;
    char *header = malloc(10000*sizeof(*header));
    int header_length = 10000;
    int i = 0;
    int chars_read = 0;
    while (c != EOF && (char)c != '\n') {
        chars_read++;
        c = getc(f);
        if (c != EOF && (char)c != '\n' && chars_read > 2) {
            header[i] = (char)c;
            i++;
            if (i == header_length-1) {
                header_length *= 2;
                header = realloc(header, header_length*sizeof(*header));
            }
        }
        if (c == EOF)
            return NULL;
    }
    header[i] = '\0';
    header = realloc(header, (i+1)*sizeof(*header));
    return header;
}

/*Reads one link from a file with the links to the coarse database and converts
  its data to a struct cb_link_to_coarse*/
struct cb_link_to_coarse *read_compressed_link(FILE *f){
    int i;
    unsigned int c = 0;
    uint16_t script_length = (uint16_t)0;
    struct cb_link_to_coarse *link = malloc(sizeof(*link));
    int chars_to_read;
    char *half_bytes;
    char *diff;

    link->coarse_seq_id = read_int_from_file(8, f);
    if (feof(f)) {
        free(link);
        return NULL;
    }

    link->original_start = read_int_from_file(8, f);
    if (feof(f)) {
        free(link);
        return NULL;
    }

    link->original_end = read_int_from_file(8, f);
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

    script_length = (uint16_t)read_int_from_file(2, f);

    chars_to_read = script_length / 2;
    if (script_length % 2 == 1)
        chars_to_read++;
    half_bytes = malloc(chars_to_read*sizeof(*half_bytes));
    for (i = 0; i < chars_to_read; i++) {
        c = getc(f);
        if (feof(f)) {
            free(link);
            return NULL;
        }
        half_bytes[i] = (char)c;
    }
    diff = half_bytes_to_ASCII(half_bytes, script_length);
    free(half_bytes);
    link->diff = diff;
    link->next = NULL;
    return link;
}

/*Takes in a pointer to the coarse.links file generated by cablast-compress and
 *the ID number of the sequence and returns the compressed sequence that the
 *data being read represents.  For this function to work properly the file
 *pointer must be pointing to the start of the header of a sequence entry in the
 *links file.
 */
struct cb_compressed_seq *get_compressed_seq(FILE *f, int id){
    char *h = get_compressed_header(f);
    struct cb_link_to_coarse *first_link = NULL;
    struct cb_link_to_coarse *last_link = NULL;
    struct cb_compressed_seq *seq = cb_compressed_seq_init(id, h);

    if (h == NULL) {
        fprintf(stderr, "Could not get compressed sequence\n");
        return NULL;
    }
    free(h);

    read_int_from_file(8, f);

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
struct cb_compressed_seq *cb_compressed_read_seq_at(
                                                 struct cb_compressed *comdb,
                                                 int32_t id){
    FILE *links = comdb->file_compressed;
    int64_t offset = cb_compressed_link_offset(comdb, id);
    if (offset < 0)
        return NULL;
    bool fseek_success = fseek(links, offset, SEEK_SET) == 0;
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
    char *h = get_compressed_header(f);
    if (h == NULL) {
        fprintf(stderr, "Could not get sequence length\n");
        return -1;
    }
    free(h);
    return read_int_from_file(8, f);
}

/*Gets the lengths in bases for all sequences in the database*/
int64_t *cb_compressed_get_lengths(struct cb_compressed *comdb){
    FILE *links = comdb->file_compressed;
    FILE *index = comdb->file_index;

    bool fseek_success;
    int64_t *lengths = NULL;
    fseek_success = fseek(index, 0, SEEK_END) == 0;
    if (!fseek_success) {
        fprintf(stderr, "Error in seeking to end of compressed.cb.index\n");
        return NULL;
    }

    int64_t num_sequences = ftell(index) / 8;
    fseek_success = fseek(index, 0, SEEK_SET) == 0;
    if (!fseek_success) {
        fprintf(stderr, "Error in seeking to start of compressed.cb.index\n");
        return NULL;
    }

    if (num_sequences > 0) {
        int i = 0;
        lengths = malloc(num_sequences*sizeof(*lengths));
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
            char c = 1;
            struct cb_link_to_coarse *current_link = read_compressed_link(f);
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
        (compressed_seqs[length])->name = header;
        (compressed_seqs[length]->links) = links;
        length++;
    }

    compressed_seqs = realloc(compressed_seqs,
                              (length+1)*sizeof(*compressed_seqs));
    compressed_seqs[length] = NULL;
    return compressed_seqs;
}

/*Gets the offset in the compressed database's compressed.cb file for the
 *link whose index is passed into id.
 */
int64_t cb_compressed_link_offset(struct cb_compressed *comdb, int id){
    int i;
    int try_off = id * 8;
    int64_t offset = (int64_t)0;
    bool fseek_success = fseek(comdb->file_index, try_off, SEEK_SET) == 0;
    int64_t mask = make_mask(8);
    if (!fseek_success) {
        fprintf(stderr, "Error in seeking to offset %d", try_off);
        return (int64_t)(-1);
    }
    return read_int_from_file(8,comdb->file_index);
}
