#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "util.h"

char *trim_space(char *s){
    return trim(s, " \r\n\t");
}

char *trim(char *s, const char *totrim){
    int32_t i, j,
             start, end,
             slen, totrimlen, newlen;
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

    end = slen - start;
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

    newlen = (end >= start) ? (end - start) : 0;

    news = malloc((newlen + 1) * sizeof(*news));
    assert(news);

    strncpy(news, s + start, newlen);
    news[newlen] = '\0';
    free(s);

    return news;
}

int32_t readline(FILE *f, char **line){
    int32_t allocated;
    char buf[1024];

    allocated = 1; /* for \0 */

    *line = malloc(allocated * sizeof(**line));
    assert(line);

    (*line)[0] = '\0';
    while (NULL != fgets(buf, 1024, f)) {
        allocated += strlen(buf);

        *line = realloc(*line, allocated * sizeof(**line));
        assert(*line);

        strcat(*line, buf);

        /* if we have found a new line, quit */
        if ((*line)[allocated - 2] == '\n')
            break;
    }
    return allocated - 1;
}

int32_t num_cpus(){
    int32_t cpus;

    cpus = (int32_t) sysconf(_SC_NPROCESSORS_ONLN);
    if (cpus <= 0)
        cpus = 1;
    return cpus;
}

char *str_slice(char *str, int32_t start, int32_t end){
    int32_t len;
    char *ret;

    len = strlen(str);

    assert(end > start);
    assert(start >= 0);
    assert(end <= len);

    ret = malloc((1 + end - start) * sizeof(*ret));
    assert(ret);

    strncpy(ret, str + start, end - start);

    return ret;
}
