#include <assert.h>
#include <stdlib.h>

#include "progress-bar.h"

char *progress_bar(int current, int full){
    char *bar;
    int bars = (int)(((float)(current+1)/full)*50);
    int b = 0;

    bar = malloc(53*sizeof(*bar));
    assert(bar);

    bar[0] = '[';
    for(b = 0; b < 50; b++)
        bar[b+1] = b < bars ? '|' : ' ';
    bar[51] = ']';
    bar[52] = '\0';
    return bar;
}
