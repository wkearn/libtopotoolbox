#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>

extern "C" {
#include "topotoolbox.h"
#include "utils.h"
}

/*
  Each pixel of the filled DEM should be greater than or equal to
  the corresponding pixel in the original DEM.
 */
int32_t test_fillsinks_ge(float *original_dem, float *filled_dem,
                          ptrdiff_t nrows, ptrdiff_t ncols) {
  for (ptrdiff_t j = 0; j < ncols; j++) {
    for (ptrdiff_t i = 0; i < nrows; i++) {
      assert(filled_dem[i + nrows * j] >= original_dem[i + nrows * j]);
    }
  }
  return 0;
}
/*
  No pixel in the filled DEM should be completely surrounded by pixels higher
  than it.
 */
int32_t test_fillsinks_filled(float *filled_dem, ptrdiff_t nrows,
                              ptrdiff_t ncols) {
  ptrdiff_t col_offset[8] = {-1, -1, -1, 0, 1, 1, 1, 0};
  ptrdiff_t row_offset[8] = {1, 0, -1, -1, -1, 0, 1, 1};

  for (ptrdiff_t j = 0; j < ncols; j++) {
    for (ptrdiff_t i = 0; i < nrows; i++) {
      float z = filled_dem[i + nrows * j];
      ptrdiff_t up_neighbor_count = 0;
      for (int32_t neighbor = 0; neighbor < 8; neighbor++) {
        ptrdiff_t neighbor_i = i + row_offset[neighbor];
        ptrdiff_t neighbor_j = j + col_offset[neighbor];

        if (neighbor_i < 0 || neighbor_i >= nrows || neighbor_j < 0 ||
            neighbor_j >= ncols) {
          continue;
        }

        if (filled_dem[neighbor_i + nrows * neighbor_j] > z) {
          up_neighbor_count++;
        }
      }
      assert(up_neighbor_count != 8);
    }
  }
  return 0;
}

/*
  Every pixel not on the boundary with no lower neighbors and fewer than
  8 higher neighbors should be labeled a flat. Likewise, every flat
  should have no lower neighbors and fewer than 8 higher neighbors.
 */
int32_t test_identifyflats_flats(int32_t *flats, float *dem, ptrdiff_t nrows,
                                 ptrdiff_t ncols) {
  ptrdiff_t col_offset[8] = {-1, -1, -1, 0, 1, 1, 1, 0};
  ptrdiff_t row_offset[8] = {1, 0, -1, -1, -1, 0, 1, 1};

  for (ptrdiff_t j = 0; j < ncols; j++) {
    for (ptrdiff_t i = 0; i < nrows; i++) {
      float z = dem[i + nrows * j];
      int32_t flat = flats[i + nrows * j];

      int32_t up_neighbor_count = 0;
      int32_t down_neighbor_count = 0;
      for (int32_t neighbor = 0; neighbor < 8; neighbor++) {
        ptrdiff_t neighbor_i = i + row_offset[neighbor];
        ptrdiff_t neighbor_j = j + col_offset[neighbor];

        if (neighbor_i < 0 || neighbor_i >= nrows || neighbor_j < 0 ||
            neighbor_j >= ncols) {
          continue;
        }

        float neighbor_height = dem[neighbor_i + nrows * neighbor_j];

        if (neighbor_height > z) {
          up_neighbor_count++;
        } else if (neighbor_height < z) {
          down_neighbor_count++;
        }
      }
      int32_t current_pixel_on_boundary =
          i == 0 || i == nrows - 1 || j == 0 || j == ncols - 1;

      assert((!current_pixel_on_boundary && (down_neighbor_count == 0) &&
              (up_neighbor_count < 8)) == ((flat & 1) == 1));
    }
  }

  return 0;
}

/*
  Every pixel that has a lower neighbor, borders a flat, and has the
  same elevation as a flat that it touches should be labeled a sill.
*/
int32_t test_identifyflats_sills(int32_t *flats, float *dem, ptrdiff_t nrows,
                                 ptrdiff_t ncols) {
  ptrdiff_t col_offset[8] = {-1, -1, -1, 0, 1, 1, 1, 0};
  ptrdiff_t row_offset[8] = {1, 0, -1, -1, -1, 0, 1, 1};

  for (ptrdiff_t j = 0; j < ncols; j++) {
    for (ptrdiff_t i = 0; i < nrows; i++) {
      float z = dem[i + j * nrows];
      int32_t flat = flats[i + j * nrows];

      int32_t down_neighbor_count = 0;
      int32_t equal_neighbor_flats = 0;

      for (int32_t neighbor = 0; neighbor < 8; neighbor++) {
        ptrdiff_t neighbor_i = i + row_offset[neighbor];
        ptrdiff_t neighbor_j = j + col_offset[neighbor];

        if (neighbor_i < 0 || neighbor_i >= nrows || neighbor_j < 0 ||
            neighbor_j >= ncols) {
          // Count boundary pixels as down neighbors so sills can
          // drain off the map.
          down_neighbor_count++;
          continue;
        }

        float neighbor_height = dem[neighbor_i + nrows * neighbor_j];
        int32_t neighbor_flat = flats[neighbor_i + nrows * neighbor_j];

        if (neighbor_height < z) {
          down_neighbor_count++;
        }

        if ((neighbor_flat & 1) && (neighbor_height == z)) {
          equal_neighbor_flats++;
        }
      }
      assert(((flat & 2) > 0) ==
             ((down_neighbor_count > 0) && (equal_neighbor_flats > 0)));
    }
  }

  return 0;
}

/*
  Every pixel that is a flat and borders a sill of the same
  height should be labeled a presill
*/
int32_t test_identifyflats_presills(int32_t *flats, float *dem, ptrdiff_t nrows,
                                    ptrdiff_t ncols) {
  ptrdiff_t col_offset[8] = {-1, -1, -1, 0, 1, 1, 1, 0};
  ptrdiff_t row_offset[8] = {1, 0, -1, -1, -1, 0, 1, 1};

  for (ptrdiff_t j = 0; j < ncols; j++) {
    for (ptrdiff_t i = 0; i < nrows; i++) {
      float z = dem[i + j * nrows];
      int32_t flat = flats[i + j * nrows];

      int32_t equal_neighbor_sills = 0;

      for (int32_t neighbor = 0; neighbor < 8; neighbor++) {
        ptrdiff_t neighbor_i = i + row_offset[neighbor];
        ptrdiff_t neighbor_j = j + col_offset[neighbor];

        if (neighbor_i < 0 || neighbor_i >= nrows || neighbor_j < 0 ||
            neighbor_j >= ncols) {
          continue;
        }

        float neighbor_height = dem[neighbor_i + nrows * neighbor_j];
        int32_t neighbor_flat = flats[neighbor_i + nrows * neighbor_j];

        if (((neighbor_flat & 2) > 0) && (neighbor_height == z)) {
          equal_neighbor_sills++;
        }
      }
      assert(((flat & 1) > 0 && (equal_neighbor_sills > 0)) ==
             ((flat & 4) > 0));
    }
  }

  return 0;
}

/*
  Costs should be zero on nonflats and positive on flats.
 */
int32_t test_gwdt_costs(float *costs, int32_t *flats, ptrdiff_t nrows,
                        ptrdiff_t ncols) {
  for (ptrdiff_t j = 0; j < ncols; j++) {
    for (ptrdiff_t i = 0; i < nrows; i++) {
      float cost = costs[i + nrows * j];
      int32_t flat = flats[i + nrows * j];

      assert((cost >= 0) && (((flat & 1) > 0) == (cost > 0)));
    }
  }

  return 0;
}

/*
  Flats neighboring other flats should all have the same nonzero
  connected component label
 */
int32_t test_gwdt_conncomps(ptrdiff_t *conncomps, int32_t *flats,
                            ptrdiff_t nrows, ptrdiff_t ncols) {
  ptrdiff_t col_offset[8] = {-1, -1, -1, 0, 1, 1, 1, 0};
  ptrdiff_t row_offset[8] = {1, 0, -1, -1, -1, 0, 1, 1};

  for (ptrdiff_t j = 0; j < ncols; j++) {
    for (ptrdiff_t i = 0; i < nrows; i++) {
      ptrdiff_t label = conncomps[i + nrows * j];
      int32_t flat = flats[i + nrows * j];

      if (!(flat & 1)) {
        continue;
      }

      for (int32_t neighbor = 0; neighbor < 8; neighbor++) {
        ptrdiff_t neighbor_i = i + row_offset[neighbor];
        ptrdiff_t neighbor_j = j + col_offset[neighbor];

        if (neighbor_i < 0 || neighbor_i >= nrows || neighbor_j < 0 ||
            neighbor_j >= ncols) {
          continue;
        }

        int32_t neighbor_flat = flats[neighbor_i + nrows * neighbor_j];
        ptrdiff_t neighbor_label = conncomps[neighbor_i + nrows * neighbor_j];
        if (neighbor_flat & 1) {
          assert(neighbor_label == label);
        }
      }
    }
  }
  return 0;
}

int32_t random_dem_test(ptrdiff_t nrows, ptrdiff_t ncols, uint32_t seed) {
  // Allocate variables

  // Input DEM
  float *dem = new float[nrows * ncols];

  // Output for fillsinks
  float *filled_dem = new float[nrows * ncols];

  // Output for identifyflats
  int32_t *flats = new int32_t[nrows * ncols];

  // Outputs for compute_costs
  ptrdiff_t *conncomps = new ptrdiff_t[nrows * ncols];
  float *costs = new float[nrows * ncols];

  // Outputs and intermediate needs for gwdt
  float *dist = new float[nrows * ncols];
  ptrdiff_t *heap = new ptrdiff_t[nrows * ncols];
  ptrdiff_t *back = new ptrdiff_t[nrows * ncols];
  ptrdiff_t *prev = new ptrdiff_t[nrows * ncols];
  for (uint32_t col = 0; col < ncols; col++) {
    for (uint32_t row = 0; row < nrows; row++) {
      dem[col * nrows + row] = 100.0f * pcg4d(row, col, seed, 1);
    }
  }

  // Run flow routing algorithms
  fillsinks(filled_dem, dem, nrows, ncols);
  ptrdiff_t count_flats = identifyflats(flats, filled_dem, nrows, ncols);
  gwdt_computecosts(costs, conncomps, flats, dem, filled_dem, nrows, ncols);
  gwdt(dist, prev, costs, flats, heap, back, nrows, ncols);

  // Number of flats identified in the test
  test_fillsinks_ge(dem, filled_dem, nrows, ncols);
  test_fillsinks_filled(filled_dem, nrows, ncols);

  ptrdiff_t count_flats = identifyflats(flats, filled_dem, nrows, ncols);

  test_identifyflats_flats(flats, filled_dem, nrows, ncols);
  test_identifyflats_sills(flats, filled_dem, nrows, ncols);
  test_identifyflats_presills(flats, filled_dem, nrows, ncols);

  gwdt_computecosts(costs, conncomps, flats, dem, filled_dem, nrows, ncols);

  test_gwdt_costs(costs, flats, nrows, ncols);
  test_gwdt_conncomps(conncomps, flats, nrows, ncols);
  ptrdiff_t test_count_flats = 0;

  // Test properties of filled DEM and identified flats
  ptrdiff_t col_offset[8] = {-1, -1, -1, 0, 1, 1, 1, 0};
  ptrdiff_t row_offset[8] = {1, 0, -1, -1, -1, 0, 1, 1};

  for (ptrdiff_t col = 0; col < ncols; col++) {
    for (ptrdiff_t row = 0; row < nrows; row++) {
      float z = filled_dem[col * nrows + row];

      int32_t current_pixel_on_border =
          row == 0 || row == nrows - 1 || col == 0 || col == ncols - 1;

      // Number of neighbors lower than the current pixel
      int32_t down_neighbor_count = 0;

      // Number of neighbors that are flats
      int32_t neighboring_flats = 0;

      // Number of neighboring flats that have the same elevation as
      // the current pixel
      int32_t equal_neighboring_flats = 0;

      // Number of neighboring sills that have the same elevation as
      // the current pixel
      int32_t equal_neighboring_sills = 0;

      for (ptrdiff_t neighbor = 0; neighbor < 8; neighbor++) {
        ptrdiff_t neighbor_row = row + row_offset[neighbor];
        ptrdiff_t neighbor_col = col + col_offset[neighbor];

        if (neighbor_row < 0 || neighbor_row >= nrows || neighbor_col < 0 ||
            neighbor_col >= ncols) {
          continue;
        }

        float neighbor_height = filled_dem[neighbor_col * nrows + neighbor_row];
        int32_t neighboring_flat = flats[neighbor_col * nrows + neighbor_row];

        if (z > neighbor_height) {
          down_neighbor_count++;
        }

        if (neighboring_flat & 1) {
          neighboring_flats++;

          if (z == neighbor_height) {
            equal_neighboring_flats++;
          }
        }

        if (neighboring_flat & 2) {
          if (z == neighbor_height) {
            equal_neighboring_sills++;
          }
        }
      }

      if (!current_pixel_on_border && up_neighbor_count < 8 &&
          down_neighbor_count == 0) {
        // This pixel is a flat
        test_count_flats++;
      }
    }
  }

  if (test_count_flats != count_flats) {
    std::cout << "Number of flats identified in the test is not equal to the "
                 "number of flats computed by identifyflats"
              << std::endl;
    return -1;
  }

  delete[] dem;
  delete[] filled_dem;
  delete[] flats;
  delete[] conncomps;
  delete[] costs;

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
