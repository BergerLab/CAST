#include <assert.h>
#include <stdbool.h>
#include <stdlib.h>
#include "DNAutils.h"

char *base_complement = "TNGNNNCNNNNNNNNNNNNANNNNNN";

/*Checks if two bases match.  If one or both bases is an N, it is an
  automatic mismatch*/
extern inline bool bases_match(char a, char b, int dir_prod){
    return a == (dir_prod > 0 ? b : base_complement[b-'A']) && a != 'N';
}

/*Takes in a pointer to the start of a k-mer in a DNA sequence and the length
  of the k-mer and returns a copy of the k-mer.*/
char *get_kmer(char *DNA_string, int k){
    int i = 0;

    char *kmer = malloc(k*sizeof(*kmer));
    assert(kmer);

    for (i = 0; i < k; i++)
        kmer[i] = DNA_string[i];
    kmer[k] = '\0';
    return kmer;
}

//Takes in a k-mer and its length and returns the k-mer's reverse complement.
char *kmer_revcomp(char *revcomp, char *kmer, int k){
    int i = 0;

    for (i = 0; i < k; i++)
        revcomp[i] = base_complement[kmer[k-i-1]-'A'];
    revcomp[k] = '\0';
    return revcomp;
}

/*Takes in a string representing a DNA sequence and its length and returns the
 *sequence's reverse complement.  If a length less than 0 is given, the length
 *used is the length of the whole string.
 */
char *string_revcomp(char *sequence, int length){
    char *revcomp;
    int i;

    /*Find the length of the sequence up to the null terminator if no length
      is given*/
    if (length < 0)
        for (length = 0; sequence[length] != '\0'; length++);

    revcomp = malloc((length+1)*sizeof(*revcomp));
    assert(revcomp);

    for (i = 0; i < length; i++)
        revcomp[i] = base_complement[sequence[length-i-1]-'A'];
    revcomp[length] = '\0';
    return revcomp;
}
