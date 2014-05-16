#include <assert.h>
#include <stdlib.h>
#include "DNAutils.h"

/*Returns the base complement for a DNA residue.  The complement for any base
  that is not A, C, G, or T in this function is N.*/
char base_complement(char base){
    switch (base) {
        case 'A':
            return 'T';
        case 'C':
            return 'G';
        case 'G':
            return 'C';
        case 'T':
            return 'A';
        default:
            return 'N';
    }
}

/*Checks if two bases match.  If one or both bases is an N, it is an
  automatic mismatch*/
int bases_match(char a, char b, int dir_prod){
    if (dir_prod > 0) /*Same-direction matches*/
        return (a == b && a != 'N') ? 1 : 0;
    else              /*Reverse-complement matches*/
        return a == base_complement(b) && a != 'N' ? 1 : 0;
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

/*Takes in a k-mer and its length and returns the k-mer's reverse complement.*/
char *kmer_revcomp(char *kmer, int k){
    int i = 0;

    char *revcomp = malloc(k*sizeof(*revcomp));
    assert(revcomp);

    for (i = 0; i < k; i++)
        revcomp[i] = base_complement(kmer[k-i-1]);
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
    revcomp = malloc((length+1) * sizeof(*revcomp));
    assert(revcomp);

    for (i = 0; i < length; i++)
        revcomp[i] = base_complement(sequence[length-i-1]);
    revcomp[length] = '\0';
    return revcomp;
}
