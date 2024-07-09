#undef NDEBUG
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <iostream>

extern "C" {
#include "topotoolbox.h"
#include "utils.h"
}

#define SQRT2f 1.41421356237309504880f

/*
  Each pixel of the filled DEM should be greater than or equal to
  the corresponding pixel in the original DEM.
 */
int32_t test_fillsinks_ge(float *original_dem, float *filled_dem,
                          ptrdiff_t dims[2], ptrdiff_t strides[2]) {
  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      assert(filled_dem[i * strides[0] + j * strides[1]] >=
             original_dem[i * strides[0] + j * strides[1]]);
    }
  }
  return 0;
}

/*
  No pixel in the filled DEM should be completely surrounded by pixels higher
  than it.
 */
int32_t test_fillsinks_filled(float *filled_dem, ptrdiff_t dims[2],
                              ptrdiff_t strides[2]) {
  ptrdiff_t i_offset[8] = {1, 0, -1, -1, -1, 0, 1, 1};
  ptrdiff_t j_offset[8] = {-1, -1, -1, 0, 1, 1, 1, 0};

  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      float z = filled_dem[i * strides[0] + j * strides[1]];
      ptrdiff_t up_neighbor_count = 0;
      for (int32_t neighbor = 0; neighbor < 8; neighbor++) {
        ptrdiff_t neighbor_i = i + i_offset[neighbor];
        ptrdiff_t neighbor_j = j + j_offset[neighbor];

        if (neighbor_i < 0 || neighbor_i >= dims[0] || neighbor_j < 0 ||
            neighbor_j >= dims[1]) {
          continue;
        }

        if (filled_dem[neighbor_i * strides[0] + neighbor_j * strides[1]] > z) {
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

      if (!((!current_pixel_on_boundary && (down_neighbor_count == 0) &&
             (up_neighbor_count < 8)) == ((flat & 1) == 1))) {
        printf("Current pixel on boundary?: %d\n", current_pixel_on_boundary);
        printf("Down neighbor count: %d\n", down_neighbor_count);
        printf("Up neighbor count: %d\n", up_neighbor_count);
        printf("Flat?: %d", flat);
        assert(0);
      }
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

/*
  Gray-weighted distances should be

  1. non-negative for all pixels
  2. 1 for presills
  3. >1 for other flats
  4. equal to the distance to the parent pixel plus the geodesic distance
  between the parent and the current pixel.
 */
int32_t test_gwdt(float *dist, ptrdiff_t *prev, float *costs, int32_t *flats,
                  ptrdiff_t nrows, ptrdiff_t ncols) {
  for (ptrdiff_t j = 0; j < ncols; j++) {
    for (ptrdiff_t i = 0; i < nrows; i++) {
      float d = dist[i + j * nrows];
      int32_t flat = flats[i + j * nrows];

      assert(d >= 0);

      if (flat & 4) {
        assert(d == 1.0);
      } else if (flat & 1) {
        assert(d > 1.0);

        ptrdiff_t parent = prev[i + j * nrows];
        ptrdiff_t parent_j = parent / nrows;
        ptrdiff_t parent_i = parent % nrows;
        float chamfer = (parent_j != j) && (parent_i != i) ? SQRT2f : 1.0f;
        float proposed_dist =
            dist[parent] + chamfer * (costs[i + j * nrows] + costs[parent]) / 2;
        assert(proposed_dist == d);
      }
    }
  }

  return 0;
}

int32_t random_dem_test(ptrdiff_t dims[2], ptrdiff_t strides[2],
                        uint32_t seed) {
  ptrdiff_t nrows = dims[0];
  ptrdiff_t ncols = dims[1];
  // Allocate variables

  // Input DEM
  float *dem =
      new float[(dims[0] - 1) * strides[0] + (dims[1] - 1) * strides[1] + 1];

  // Output for fillsinks
  float *filled_dem =
      new float[(dims[0] - 1) * strides[0] + (dims[1] - 1) * strides[1] + 1];

  // Output for identifyflats
  int32_t *flats =
      new int32_t[(dims[0] - 1) * strides[0] + (dims[1] - 1) * strides[1] + 1];

  // Outputs for compute_costs
  ptrdiff_t *conncomps = new ptrdiff_t[(dims[0] - 1) * strides[0] +
                                       (dims[1] - 1) * strides[1] + 1];
  float *costs =
      new float[(dims[0] - 1) * strides[0] + (dims[1] - 1) * strides[1] + 1];

  // Outputs and intermediate needs for gwdt
  float *dist =
      new float[(dims[0] - 1) * strides[0] + (dims[1] - 1) * strides[1] + 1];
  ptrdiff_t *heap = new ptrdiff_t[(dims[0] - 1) * strides[0] +
                                  (dims[1] - 1) * strides[1] + 1];
  ptrdiff_t *back = new ptrdiff_t[(dims[0] - 1) * strides[0] +
                                  (dims[1] - 1) * strides[1] + 1];
  ptrdiff_t *prev = new ptrdiff_t[(dims[0] - 1) * strides[0] +
                                  (dims[1] - 1) * strides[1] + 1];

  for (uint32_t j = 0; j < dims[1]; j++) {
    for (uint32_t i = 0; i < dims[0]; i++) {
      dem[i * strides[0] + j * strides[1]] = 100.0f * pcg4d(i, j, seed, 1);
    }
  }

  // Run flow routing algorithms
  fillsinks(filled_dem, dem, dims, strides);

  test_fillsinks_ge(dem, filled_dem, dims, strides);
  test_fillsinks_filled(filled_dem, dims, strides);

  identifyflats(flats, filled_dem, nrows, ncols);

  test_identifyflats_flats(flats, filled_dem, nrows, ncols);
  test_identifyflats_sills(flats, filled_dem, nrows, ncols);
  test_identifyflats_presills(flats, filled_dem, nrows, ncols);

  gwdt_computecosts(costs, conncomps, flats, dem, filled_dem, nrows, ncols);

  test_gwdt_costs(costs, flats, nrows, ncols);
  test_gwdt_conncomps(conncomps, flats, nrows, ncols);

  gwdt(dist, prev, costs, flats, heap, back, nrows, ncols);
  test_gwdt(dist, prev, costs, flats, nrows, ncols);

  delete[] dem;
  delete[] filled_dem;
  delete[] flats;
  delete[] conncomps;
  delete[] costs;
  delete[] dist;
  delete[] heap;
  delete[] back;
  delete[] prev;

  return 0;
}

int main(int argc, char *argv[]) {
  ptrdiff_t cm_dims[2] = {100, 200};
  ptrdiff_t cm_strides[2] = {1, 100};

  ptrdiff_t rm_dims[2] = {200, 100};
  ptrdiff_t rm_strides[2] = {1, 200};

  ptrdiff_t strided_cm_dims[2] = {100, 200};
  ptrdiff_t strided_cm_strides[2] = {3, 300};

  for (uint32_t test = 0; test < 50; test++) {
    int32_t result = 0;
    result = random_dem_test(cm_dims, cm_strides, 3 * test);
    if (result < 0) {
      return result;
    }

    result = random_dem_test(rm_dims, rm_strides, 3 * test + 1);
    if (result < 0) {
      return result;
    }

    result = random_dem_test(strided_cm_dims, strided_cm_strides, 3 * test + 2);
    if (result < 0) {
      return result;
    }
  }
}
