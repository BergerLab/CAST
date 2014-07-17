#include <stdbool.h>

#ifndef __CABLAST_DNAUTILS_H__
#define __CABLAST_DNAUTILS_H__

char *base_complement;

bool bases_match(char a, char b, int dir_prod);
void kmer_revcomp(char *revcomp, char *kmer, int k);
char *string_revcomp(char *sequence, int length);

#endif
