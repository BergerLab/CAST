#include <assert.h>
#include <stdbool.h>
#include <stdlib.h>

#include "DNAutils.h"

//Global string of base complements of all leters of the alphabet
char *base_complement = "TNGNNNCNNNNNNNNNNNNANNNNNN";

/*Checks if two bases match.  If one or both bases is an N, it is an
  automatic mismatch*/
extern inline bool bases_match(char a, char b, int dir_prod){
    return a == (dir_prod > 0 ? b : base_complement[b-'A']) && a != 'N';
}

/*Takes in a k-mer, its length, and a char * "revcomp" and puts the k-mer's
 *reverse complement in revcomp.
 *
 *Please note that this function does not null-terminate the reverse complement
 *since a null terminator is never needed in any of the code that uses the
 *reverse complement of a k-mer.
 */
extern inline void kmer_revcomp(char *revcomp, char *kmer, int k){
    for (int i = 0; i < k; i++)
        revcomp[i] = base_complement[kmer[k-i-1]-'A'];
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
