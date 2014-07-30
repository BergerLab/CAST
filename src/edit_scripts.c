#include <assert.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include "DNAutils.h"
#include "edit_scripts.h"

int minimum(int a, int b){return a<b?a:b;}
int maximum(int a, int b){return a>b?a:b;}

//Converts an integer to a string containing its octal representation.
char *to_octal_str(int i) {
    char *buf = malloc(16*sizeof(*buf));
    assert(buf);

    sprintf(buf, "%o", i);

    return buf;
}

/*Converts a character in an edit script to a half-byte representing the 
  character.*/
char to_half_byte(char c){
    switch (c) {
        case '0': return (char)0;
        case '1': return (char)1;
        case '2': return (char)2;
        case '3': return (char)3;
        case '4': return (char)4;
        case '5': return (char)5;
        case '6': return (char)6;
        case '7': return (char)7;
        case 'A': return (char)8;
        case 'C': return (char)9;
        case 'G': return (char)10;
        case 'T': return (char)11;
        case '-': return (char)12;
        case 'i': return (char)13;
        case 's': return (char)14;

        /*For direction bytes with a 1 as the first bit to indicate that
         *the script is from a match, move the 1 to the leftmost bit of the
         *half-byte.
         */
        case '0' | (char)0x80: return (char)8;
        case '9' | (char)0x80: return (char)9;

        default:  return (char)15; //'N' is represented as the half byte 1111.
    }
}

//Converts a half-byte to its corresponding edit script character
char half_byte_to_char(char h){
    if (h < 8) //h is an octal digit.
        return '0' + h;
    switch (h) {
        case (char)8:  return 'A';
        case (char)9:  return 'C';
        case (char)10: return 'G';
        case (char)11: return 'T';
        case (char)12: return '-';
        case (char)13: return 'i';
        case (char)14: return 's';
        default:       return 'N';
    }
}

/*Takes in an edit script in ASCII format and returns a new string containing
  the edit script in half-byte format.*/
char *edit_script_to_half_bytes(char *edit_script){
    int i = 0, length = 0, odd;
    char *half_bytes;

    while (edit_script[length] != '\0')
        length++;
    odd    = length % 2;
    length = length / 2 + odd;

    half_bytes = malloc(length*sizeof(*half_bytes));
    assert(half_bytes);
 
    while (edit_script[i] != '\0') {
        if (i%2 == 0)
            half_bytes[i/2] = (char)0;
        else
            half_bytes[i/2] <<= 4;

        //Insert the current half byte into the current byte.
        half_bytes[i/2] |= to_half_byte(edit_script[i]);
        i++;
    }
    /*If the length of the edit script is odd, to signify the end of the edit
     *script add a 1 to signify the end of the script since this is incorrect
     *edit script syntax.
     */
    if (odd) {
        half_bytes[i/2] <<= 4;
        half_bytes[i/2] |= (char)1;
    }

    return half_bytes;
}

/*Takes in an edit script in half-byte format and returns a new string
  containing the edit script in ASCII format.*/
char *half_bytes_to_ASCII(char *half_bytes, int length){
    char *edit_script = malloc((length+1)*sizeof(*edit_script));
    assert(edit_script);

    for (int i = 0; i < length; i++) {
        //Copy the left half-byte of the current byte.
        if (i % 2 == 0) {
            char left = half_bytes[i/2] & (((char)15) << 4);
            left >>= 4;
            left &= (char)15;
            edit_script[i] = half_byte_to_char(left);

            /*Handle the direction byte so that the match bit of the
             *half-byte is added to the direction byte separately from the
             *direction bit.
             */
            if (i == 0) {
                left &= (char)1;
                edit_script[i] = half_byte_to_char(left);
                if ((half_bytes[0] & (char)0x80) != 0)
                    edit_script[i] |= (char)0x80;
            }
        }
        else //Copy the right half-byte of the current byte
            edit_script[i] = half_byte_to_char(half_bytes[i/2] & (char)15);
    }
    edit_script[length] = '\0';

    return edit_script;
}

/*Takes in two strings, a bool representing whether or not they are in the same
 *direction, and the length of the strings and returns an edit script that can
 *convert the reference string to the original string.
 */
char *make_edit_script(char *str, char *ref, bool dir, int length){
    /*direction has its first bit set to 1 to indicate that the edit script
      was made from a match*/
    int last_edit = 0, current = 1;
    char *edit_script, *octal, direction = dir ? '0' : '1';
    bool insert_open = false, subdel_open = false;

    edit_script = malloc(3*length*sizeof(*edit_script));
    assert(edit_script);

    direction |= ((char)0x80);

    edit_script[0] = direction;
    for (int i = 0; i < length; i++) {
        if (str[i] == ref[i]) {
            insert_open = false;
            subdel_open = false;
        }
        else { //Mismatch
            //Insertion in str relative to ref (i.e., gap in ref)
            if (ref[i] == '-') {
                subdel_open = false;

                //Indicate start of insertion
                if (!insert_open) { 
                    insert_open = true;
                    octal = to_octal_str(i - last_edit);
                    edit_script[current++] = 'i';
                    for (int j = 0; octal[j] != '\0'; j++)
                        edit_script[current++] = octal[j];
                    last_edit = i;
                    free(octal);
                }
                edit_script[current++] = str[i];
            }
            /*Substitution or deletion in str (represented in script by '-')
              relative to ref*/
            else {
                insert_open = false;

                //Indicate start of subdel.
                if (!subdel_open) { 
                    subdel_open = true;
                    octal = to_octal_str(i - last_edit);
                    edit_script[current++] = 's';
                    for (int j = 0; octal[j] != '\0'; j++)
                        edit_script[current++] = octal[j];
                    last_edit = i;
                    free(octal);
                }
                edit_script[current++] = str[i];
            }
        }
    }

    edit_script = realloc(edit_script, (current+1)*sizeof(*edit_script));
    assert(edit_script);

    edit_script[current] = '\0';

    return edit_script;
}

/*Reads one edit worth of info (or fails if current decoded char is digit)
 *moves pos accordingly.  If there is an edit to read, next_edit reads in
 *the edit and returns true.  Otherwise, next_edit returns false.
 */
bool next_edit(char *edit_script, int *pos, struct edit_info *edit){
    int edit_length = 0, i = 0;

    if (isdigit(edit_script[(*pos)]) || edit_script[(*pos)] == '\0')
        return false;

    edit->is_subdel  = edit_script[(*pos)++] == 's';
    edit->last_dist  = 0;
    edit->str_length = 0;

    if (edit->str != NULL)
        free(edit->str);

    while (isdigit(edit_script[(*pos)])) {
        edit->last_dist *= 8; //Octal encoding
        edit->last_dist += edit_script[(*pos)++] - '0';
    }
    while (isupper(edit_script[(*pos)+edit_length]) ||
                   edit_script[(*pos)+edit_length] == '-')
        edit_length++;

    edit->str = malloc((edit_length+1)*sizeof(*edit_script));
    assert(edit->str);

    edit->str_length = edit_length;
    while (isupper(edit_script[(*pos)]) || edit_script[(*pos)] == '-')
        edit->str[i++] = edit_script[(*pos)++];

    return true;
}

/*Takes in as input an edit script in ASCII format, a DNA sequence to read, and
 *the length of the sequence and returns the sequence created when the edit
 *script is applied to the sequence passed in.
 */
char *read_edit_script(char *edit_script, char *orig, int length){
    struct edit_info edit;
    int orig_pos = 0, last_edit_str_len = 0, current = 0, script_pos = 1;

    char *str = malloc((2*length+1)*sizeof(*str));
    assert(str);

    edit.str = NULL;

    while (next_edit(edit_script, &script_pos, &edit)) {
        //Chunk after previous edit
        for (int i = 0; i < edit.last_dist - last_edit_str_len; i++)
            str[current++] = orig[orig_pos+i];

        //Update position in original string.
        orig_pos += edit.last_dist - last_edit_str_len;

        //Append replacement string in edit script; get rid of dashes.
        for (int i = 0; i < edit.str_length; i++)
            if (edit.str[i] != '-')
                str[current++] = edit.str[i];

        //Skip subdel along original string.
        if (edit.is_subdel) orig_pos += edit.str_length;

        last_edit_str_len = edit.str_length;
    }
    free(edit.str);

    while (orig_pos < length)
        str[current++] = orig[orig_pos++];

    str = realloc(str, (current+1)*sizeof(*str));
    assert(str);

    str[current] = '\0';
    if ((edit_script[0] & ((char)0x7f)) == '1') {
        char *str_fwd = str;

        str = string_revcomp(str_fwd, -1);
        free(str_fwd);
    }

    return str;
}

/*Takes in as input a string and returns a copy of the string with the '-'
  characters removed*/
char *no_dashes(char *sequence){
    int length, bases = 0, j = 0;
    char *n;

    //Get the length of the final string.
    for (length = 0; sequence[length] != '\0'; length++)
        if (sequence[length] != '-')
            bases++;

    n = malloc((bases+1)*sizeof(*n));
    assert(n);

    n[bases] = '\0';
    for (int i = 0; i < length; i++)
        if (sequence[i] != '-')
            n[j++] = sequence[i];

    return n;
}
