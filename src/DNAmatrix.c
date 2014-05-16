#include <stdint.h>

#include "DNAalphabet.h"

const int32_t DNA_MATRIX[DNA_SIZE][DNA_SIZE] = {
  { 1,-3,-3,-3,-3},
  {-3, 1,-3,-3,-3},
  {-3,-3, 1,-3,-3},
  {-3,-3,-3, 1,-3},
  {-3,-3,-3,-3,-3}
};
