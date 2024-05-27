#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "topotoolbox.h"

static float pcg4d(uint32_t a, uint32_t b, uint32_t c, uint32_t d) {
  uint32_t x = a * 1664525u + 1013904223u;
  uint32_t y = b * 1664525u + 1013904223u;
  uint32_t z = c * 1664525u + 1013904223u;
  uint32_t w = d * 1664525u + 1013904223u;

  x += y * w;
  y += z * x;
  z += x * y;
  w += y * z;

  x ^= x >> 16;
  y ^= y >> 16;
  z ^= z >> 16;
  w ^= w >> 16;

  x += y * w;
  y += z * x;
  z += x * y;
  w += y * z;

  return (float)(w >> 8) / (1 << 24);
}

int main(int argc, char *argv[]) {
  ptrdiff_t nrows = 5;
  ptrdiff_t ncols = 5;

  float *dem = malloc(nrows * ncols * sizeof(float));
  if (!dem) {
    return -1;
  }

  float *excess = malloc(nrows * ncols * sizeof(float));
  if (!excess) {
    return -1;
  }

  ptrdiff_t *heap = malloc(nrows * ncols * sizeof(ptrdiff_t));
  if (!heap) {
    return -1;
  }

  ptrdiff_t *back = malloc(nrows * ncols * sizeof(ptrdiff_t));
  if (!back) {
    return -1;
  }

  float *threshold = malloc(nrows * ncols * sizeof(float));
  if (!threshold) {
    return -1;
  }

  for (ptrdiff_t j = 0; j < ncols; j++) {
    for (ptrdiff_t i = 0; i < nrows; i++) {
      ptrdiff_t idx = j * nrows + i;
      dem[idx] = 100.0 * pcg4d(i, j, 0, 0);
      threshold[idx] = 1.0f / sqrtf(3.0);  // 30 degree threshold
    }
  }

  excesstopography(excess, heap, back, dem, threshold, 30.0, nrows, ncols);

  for (ptrdiff_t j = 0; j < ncols; j++) {
    for (ptrdiff_t i = 0; i < nrows; i++) {
      ptrdiff_t idx = j * nrows + i;
      if (excess[idx] > dem[idx]) {
        printf("(%ld,%ld) excess topography is greater than DEM",i,j);
        return -1;
      };
    }
  }
  return 0;
}
