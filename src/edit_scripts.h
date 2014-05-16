#ifndef __CABLAST_EDITSCRIPTS_H__
#define __CABLAST_EDITSCRIPTS_H__

#include "coarse.h"
#include "link_to_coarse.h"

struct edit_info{
    bool is_subdel;
    int last_dist;
    char *str;
    int str_length;
};

char to_half_byte(char c);
char half_byte_to_char(char h);
char *edit_script_to_half_bytes(char *edit_script);
char *half_bytes_to_ASCII(char *half_bytes, int length);
char *to_octal_str(int i);
char *make_edit_script(char *str, char *ref, bool dir, int length);
char *read_edit_script(char *edit_script, char *orig, int length);
void decode_edit_script(char *orig, int dest_len, int dest0_coord,
                        struct cb_coarse *coarsedb,
                        struct cb_link_to_coarse *link);
bool next_edit(char *edit_script, int *pos, struct edit_info *edit);
char *no_dashes(char *sequence);
#endif
