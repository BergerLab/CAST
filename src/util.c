#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "util.h"

//Calls trim to remove leading and trailing whitespace
char *trim_space(char *s){
    return trim(s, " \r\n\t");
}

/*Takes in a string and an array of characters and removes all leading and
 *trailing characters from the string until the first character from the start
 *and end not in the "totrim" array.
 */
char *trim(char *s, const char *totrim){
    int32_t i, j, start, end, slen, totrimlen, newlen;
    char *news;
    bool trimmed;

    slen = strlen(s);
    totrimlen = strlen(totrim);

    start = 0;
    for (i = 0; i < slen; i++) {
        trimmed = false;
        for (j = 0; j < totrimlen; j++)
            if (totrim[j] == s[i]) {
                start++;
                trimmed = true;
                break;
            }
        if (!trimmed)
            break;
    }

    for (i = slen - 1; i >= 0; i--) {
        trimmed = false;
        for (j = 0; j < totrimlen; j++)
            if (totrim[j] == s[i]) {
                end--;
                trimmed = true;
            }
        if (!trimmed)
            break;
    }
    end = i+1;

    newlen = (end >= start) ? (end - start) : 0;

    news = malloc((newlen+1)*sizeof(*news));
    assert(news);

    strncpy(news, s + start, newlen);
    news[newlen] = '\0';
    free(s);

    return news;
}

/*Takes in a file pointer and a pointer to a string and reads in every
 *character in the current line of the file, stopping at a newline or the end of
 *the file.  The line is stored in "line" and the function returns the length of
 *the line.
 */
int32_t readline(FILE *f, char **line){
    int32_t allocated;
    char buf[1024];

    allocated = 1; //For the null terminator

    *line = malloc(allocated*sizeof(**line));
    assert(line);

    (*line)[0] = '\0';
    while (NULL != fgets(buf, 1024, f)) {
        allocated += strlen(buf);

        *line = realloc(*line, allocated*sizeof(**line));
        assert(*line);

        strcat(*line, buf);

        //If we have found a new line, quit.
        if ((*line)[allocated - 2] == '\n')
            break;
    }
    return allocated - 1;
}

//Gets the number of processors on the machine.
int32_t num_cpus(){
    int32_t cpus;

    cpus = (int32_t)sysconf(_SC_NPROCESSORS_ONLN);
    if (cpus <= 0)
        cpus = 1;
    return cpus;
}

/*Takes in a string and starting and ending indices and returns a new string
 *containing the characters of "str" from the starting index to the ending
 *index.
 */
char *str_slice(char *str, int32_t start, int32_t end){
    char *ret;

    assert(end > start);
    assert(start >= 0);
    assert(end <= strlen(str));

    ret = malloc((1+end-start)*sizeof(*ret));
    assert(ret);

    strncpy(ret, str + start, end - start);

    return ret;
}

/*Takes in two strings and checks if sub occurs in str*/
bool is_substring(char *sub, char *str){
    assert(sub != NULL && str != NULL);

    int len_sub = strlen(sub), len_str = strlen(str), i, j;
    for (i = 0; i < len_str-len_sub; i++) {
        for (j = 0; j < len_sub; j++)
            if (sub[j] != str[i+j])
                break;
        if (j == len_sub)
            return true;
    }

    return false;
}

/*Takes in a line and a character to split the line on and splits the line on
 *that character, returning a string for each section of the original string
 *separated by that character.
 */
char **split_char(char *line, char split){
    assert(line);

    int len = strlen(line), num_sections = 1;

    for (int i = 0; i < len; i++)
        if (line[i] == split)
            num_sections++;

    int *indices = malloc((num_sections+1)*sizeof(*indices));
    assert(indices);

    char **words = malloc((num_sections+1)*sizeof(*words));
    assert(words);
   
    indices[0] = -1;
    indices[num_sections] = len;
    
    for (int i = 0, j = 1; i < len; i++) {
        if (line[i] == split)
            indices[j++] = i;
    }

    words[num_sections] = NULL;

    for (int i = 0; i < num_sections; i++) {
        words[i] = malloc((indices[i+1]-indices[i]+1)*sizeof(*words));
        assert(words[i]);

        strncpy(words[i], line+indices[i]+1, indices[i+1]-indices[i]-1);
        words[i][indices[i+1]-indices[i]-1] = '\0';
    }

    free(indices);

    return words;
}
