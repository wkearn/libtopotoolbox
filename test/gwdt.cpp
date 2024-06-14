#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>

extern "C" {
#include "topotoolbox.h"
}

#define SQRT2f 1.41421356237309504880f

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
  ptrdiff_t *prev = new ptrdiff_t[nrows * ncols];

  gwdt(dist, prev, costs, flats, heap, back, nrows, ncols);

  for (ptrdiff_t col = 0; col < ncols; col++) {
    for (ptrdiff_t row = 0; row < nrows; row++) {
      ptrdiff_t idx = col * nrows + row;
      if (dist[idx] < 0) {
        // Distances for should be non-negative for all pixels
        std::cout << "Pixel (" << row << ", " << col
                  << ") distance is negative." << std::endl;
        return -1;
      }
      if (flats[idx] & 4) {
        if (dist[idx] > 1.0) {
          // Distances should be 1.0 for presills
          std::cout << "Presill (" << row << ", " << col
                    << ") distance is greater than 1." << std::endl;
          return -1;
        }
      } else if (flats[idx] & 1) {
        if (dist[idx] <= 1.0) {
          // Distances should be > 1.0 for non-presill flats
          std::cout << "Flat pixel (" << row << ", " << col
                    << ") distance is less than 1." << std::endl;
          return -1;
        }

        ptrdiff_t parent = prev[idx];
        ptrdiff_t parent_col = parent / nrows;
        ptrdiff_t parent_row = parent % nrows;
        float chamfer =
            (parent_col != col) && (parent_row != row) ? SQRT2f : 1.0;
        float proposed_dist =
            dist[parent] + chamfer * (costs[idx] + costs[parent]) / 2;
        if (proposed_dist != dist[idx]) {
          // The distance of the current pixel should be the distance
          // to its parent pixel plus the geodesic distance between
          // the parent and the current pixel.
          std::cout << "Computed distance for pixel (" << row << ", " << col
                    << ") distance is incorrect: " << proposed_dist << ", "
                    << dist[idx] << std::endl;
          return -1;
        }

      } else {
        if (dist[idx] > 0) {
          // Distances should be == 0.0 for all non-flat pixels
          std::cout << "Nonflat pixel (" << row << ", " << col
                    << ") distance is nonzero." << std::endl;
          return -1;
        }
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
