#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>

extern "C" {
#include "topotoolbox.h"
}

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

int32_t random_dem_test(ptrdiff_t nrows, ptrdiff_t ncols, uint32_t seed) {
  float *dem = new float[nrows * ncols];

  for (ptrdiff_t col = 0; col < ncols; col++) {
    for (ptrdiff_t row = 0; row < nrows; row++) {
      dem[col * nrows + row] = 100.0f * pcg4d(row, col, seed, 1);
    }
  }
  
  float *filled_dem = new float[nrows * ncols];
  int32_t *flats = new int32_t[nrows * ncols];

  ptrdiff_t *conncomps = new ptrdiff_t[nrows * ncols];
  float *costs = new float[nrows * ncols];

  fillsinks(filled_dem, dem, nrows, ncols);
  identifyflats(flats, filled_dem, nrows, ncols);
  gwdt_computecosts(costs, conncomps, flats, dem, filled_dem, nrows, ncols);

  float *dist = new float[nrows * ncols];
  ptrdiff_t *heap = new ptrdiff_t[nrows * ncols];
  ptrdiff_t *back = new ptrdiff_t[nrows * ncols];
  
  gwdt(dist, costs, flats, heap, back, nrows, ncols);

  return 0;
}

int main(int argc, char *argv[]) {
  ptrdiff_t nrows = 5;
  ptrdiff_t ncols = 5;

  for (uint32_t test = 0; test < 100; test++) {
    int32_t result = random_dem_test(nrows, ncols, test);
    if (result < 0) {
      return result;
    }
  }

  return 0;
}
