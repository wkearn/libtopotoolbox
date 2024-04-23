#include <cstddef>
#include <cstdint>
#include <iostream>

extern "C" {
#include "topotoolbox.h"
}

/*
  PCG4D hash function


  Jarzynski, Mark and Olano, Marc. (2020). Hash functions for GPU
  rendering. Journal of Computer Graphics Techniques. Vol 9,
  No. 3. 21-38.
 */
float pcg4d(uint32_t a, uint32_t b, uint32_t c, uint32_t d) {
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
  // Initialize a random DEM
  float *dem = new float[nrows * ncols];

  for (uint32_t col = 0; col < ncols; col++) {
    for (uint32_t row = 0; row < nrows; row++) {
      dem[col * nrows + row] = 100.0f * pcg4d(row, col, seed, 1);
    }
  }

  // Allocate output for fillsinks
  float *filled_dem = new float[nrows * ncols];

  // Allocate output for identify flats
  int32_t *flats = new int32_t[nrows * ncols];

  fillsinks(filled_dem, dem, nrows, ncols);
  identifyflats(flats, filled_dem, nrows, ncols);

  // Test properties of filled DEM and identified flats
  ptrdiff_t col_offset[8] = {-1, -1, -1, 0, 1, 1, 1, 0};
  ptrdiff_t row_offset[8] = {1, 0, -1, -1, -1, 0, 1, 1};

  for (ptrdiff_t col = 1; col < ncols - 1; col++) {
    for (ptrdiff_t row = 1; row < nrows - 1; row++) {
      float z = filled_dem[col * nrows + row];
      int32_t flat = flats[col * nrows + row];

      // Each pixel of the filled raster should be >= the DEM
      if (z < dem[col * nrows + row]) {
        std::cout << "Pixel (" << row << ", " << col << ") is below the DEM"
                  << std::endl;
        std::cout << "Value: " << z << std::endl;
        std::cout << "DEM: " << dem[col * nrows + row] << std::endl;
        return -1;
      }

      // No pixel of the filled raster should be surrounded by
      // neighbors that are higher than it
      int32_t sink_neighbor_count = 0;

      // No flat pixel should have a neighbor that is lower than it
      int32_t flat_neighbor_count = 0;

      // Every sill pixel should have at least one neighbor lower than it
      int32_t sill_neighbor_count = 0;

      // Every sill pixel should border a flat
      int32_t sill_neighboring_flats = 0;

      for (ptrdiff_t neighbor = 0; neighbor < 8; neighbor++) {
        ptrdiff_t neighbor_row = row + row_offset[neighbor];
        ptrdiff_t neighbor_col = col + col_offset[neighbor];

        float neighbor_height = filled_dem[neighbor_col * nrows + neighbor_row];
        int32_t neighboring_flat = flats[neighbor_col * nrows + neighbor_row];

        if (z < neighbor_height) {
          sink_neighbor_count++;
        }

        if ((flat == 1) && (z > neighbor_height)) {
          flat_neighbor_count++;
        }

        if ((flat == 2) && (z > neighbor_height)) {
          sill_neighbor_count++;
        }

        if ((flat == 2) && (neighboring_flat == 1)) {
          sill_neighboring_flats++;
        }
      }

      if (sink_neighbor_count == 8) {
        std::cout << "Pixel (" << row << ", " << col << ") is a sink"
                  << std::endl;
        return -1;
      }

      if ((flat == 2) && (flat_neighbor_count > 0)) {
        std::cout << "Pixel (" << row << ", " << col
                  << ") is a flat but has a lower neighbor" << std::endl;
        return -1;
      }

      if ((flat == 2) && (sill_neighbor_count == 0)) {
        std::cout << "Pixel (" << row << ", " << col
                  << ") is a sill but has no lower neighbor" << std::endl;
        return -1;
      }

      if ((flat == 2) && (sill_neighboring_flats == 0)) {
        std::cout << "Pixel (" << row << ", " << col
                  << ") is a sill but does not border a flat" << std::endl;
        return -1;
      }
    }
  }

  return 0;
}

int main(int argc, char *argv[]) {
  ptrdiff_t nrows = 100;
  ptrdiff_t ncols = 200;

  for (uint32_t test = 0; test < 100; test++) {
    int32_t result = random_dem_test(nrows, ncols, test);
    if (result < 0) {
      return result;
    }
  }
}
