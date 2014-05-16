#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include "bitpack.h"

/*Takes in as input a number of bits and returns a 64-bit mask with that
  number of bits on the right of the mask set to 1 and all bits to the left
  set to 0.*/
uint64_t make_mask(int bits){
    uint64_t mask = (uint64_t)1;
    if (bits == 64)
        return ~(uint64_t)0;
    if (bits == 0)
        return (uint64_t)0;
    mask <<= bits;
    return mask-1;
}

/*Left shift for 64-bit integers*/
uint64_t shift_left(uint64_t x, int bits){
    return x<<bits;
}

/*Right shift for 64-bit integers*/
uint64_t shift_right(uint64_t x, int bits){
    uint64_t mask = make_mask(64-bits);
    x >>= bits;
    return x & mask;
}

/*Takes in an integer in uint64_t format and the number of bytes in that int
  and returns its representation as an array of characters.*/
char *read_int_to_bytes(uint64_t number, int length){
    int i;
    uint64_t mask = make_mask(8);

    char *bytes = malloc(length * sizeof(*bytes));
    assert(bytes);

    for (i = length-1; i >= 0; i--)
        bytes[length-i-1] = (char)(shift_right(number, 8*i) & mask);
    return bytes;
}

/*Takes in an integer in unint64_t format, the number of bytes in that int,
  and a file pointer and outputs the int to the file.*/
void output_int_to_file(uint64_t number, int length, FILE *f){
    int i;
    char *bytes = read_int_to_bytes(number, length);
    for (i = 0; i < length; i++)
        putc(bytes[i], f);
    free(bytes);
}

/*Takes in a number of bytes to read and a file pointer and reads that number
  of bytes from the file and copies them into a 64-bit integer.*/
uint64_t read_int_from_file(int length, FILE *f){
    int i;
    uint64_t bytes = (uint64_t)0;
    uint64_t mask = make_mask(8);
    for (i = 0; i < length; i++) {
        uint64_t current_byte = ((uint64_t)getc(f)) & mask;
        bytes <<= 8;
        bytes |= current_byte;
        if (feof(f))
            break;
    }
    return bytes;
}
