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
  float *excess = new float[nrows * ncols];
  ptrdiff_t *heap = new ptrdiff_t[nrows * ncols];
  ptrdiff_t *back = new ptrdiff_t[nrows * ncols];
  float *threshold = new float[nrows * ncols];

  for (ptrdiff_t col = 0; col < ncols; col++) {
    for (ptrdiff_t row = 0; row < nrows; row++) {
      ptrdiff_t idx = col * nrows + row;
      dem[idx] = 100.0f * pcg4d(row, col, seed, 0);
      threshold[idx] = 1.0f / sqrtf(3.0);  // 30 degree threshold
    }
  }
  float cellsize = 30.0;
  excesstopography(excess, heap, back, dem, threshold, cellsize, nrows, ncols);

  for (ptrdiff_t col = 0; col < ncols; col++) {
    for (ptrdiff_t row = 0; row < nrows; row++) {
      ptrdiff_t idx = col * nrows + row;

      if (excess[idx] > dem[idx]) {
        std::cout << "Seed (" << seed << "): Pixel (" << row << ", " << col
                  << ") excess topography is greater than DEM" << std::endl;
        return -1;
      };

      // Test the upwind gradient
      if (col > 0 && col < ncols - 1 && row > 0 && row < nrows - 1) {
        float north_gradient =
            (excess[col * nrows + row] - excess[col * nrows + row - 1]) /
            cellsize;
        float south_gradient =
            (excess[col * nrows + row - 1] - excess[col * nrows + row]) /
            cellsize;

        float ns_gradient =
            std::fmaxf(0.0f, std::fmaxf(north_gradient, -south_gradient));

        float west_gradient =
            (excess[col * nrows + row] - excess[(col - 1) * nrows + row]) /
            cellsize;
        float east_gradient =
            (excess[(col + 1) * nrows + row] - excess[col * nrows + row]) /
            cellsize;

        float ew_gradient =
            std::fmaxf(0.0f, std::fmaxf(west_gradient, -east_gradient));
        float g = ns_gradient * ns_gradient + ew_gradient * ew_gradient;

        if (g > 1.0f / threshold[col * nrows + row]) {
          std::cout << "Seed (" << seed << "): Pixel (" << row << ", " << col
                    << ") discrete gradient " << g
                    << " is greater than threshold" << std::endl;
          return -1;
        };
      }
    }
  }

  return 0;
}

int main(int argc, char *argv[]) {
  ptrdiff_t nrows = 200;
  ptrdiff_t ncols = 100;

  for (uint32_t test = 0; test < 100; test++) {
    int32_t result = random_dem_test(nrows, ncols, test);
    if (result < 0) {
      return result;
    }
  }

  return 0;
}
