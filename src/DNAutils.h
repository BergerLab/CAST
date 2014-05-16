#ifndef __CABLAST_DNAUTILS_H__
#define __CABLAST_DNAUTILS_H__

int bases_match(char a, char b, int dir_prod);
char base_complement(char base);
char *get_kmer(char *DNA_string, int k);
char *kmer_revcomp(char *kmer, int k);
char *string_revcomp(char *sequence, int length);

#endif
