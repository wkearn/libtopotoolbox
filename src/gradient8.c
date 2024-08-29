#define TOPOTOOLBOX_BUILD

#include <math.h>
#include <omp.h>
#include <stddef.h>
#include <stdint.h>

#include "topotoolbox.h"

#define PI 3.14159265358979323846f
#define SQRT2 1.41421356237309504880f
/*
  Compute the gradient of each cell in a digital elevation model.

  The arrays pointed to by `gradient` and `dem` should represent
  two-dimensional arrays of size specified in dims. The gradient is
  computed as the maximum difference in elevation between a cell
  and its 8 neighbors.

  int use_mp:
    0   --> don't use multiprocessing
    1   --> use multiprocessing

  char unit:
    't' --> tangent (default)
    'r' --> radian
    'd' --> degree
    's' --> sine
    'p' --> percent
 */

TOPOTOOLBOX_API
void gradient8(float *output, float *dem, float cellsize, char unit, int use_mp,
               ptrdiff_t dims[2]) {
  ptrdiff_t i_offset[8] = {-1, -1, -1, 0, 0, 1, 1, 1};
  ptrdiff_t j_offset[8] = {-1, 0, 1, -1, 1, -1, 0, 1};

#ifdef _MSC_VER
  ptrdiff_t j;
#pragma omp parallel for if (use_mp)
  for (j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
#else
  ptrdiff_t i, j;
#pragma omp parallel for collapse(2) if (use_mp)
  for (j = 0; j < dims[1]; j++) {
    for (i = 0; i < dims[0]; i++) {
#endif
      float max_gradient = 0;

      for (int k = 0; k < 8; k++) {
        ptrdiff_t neighbour_i = i + i_offset[k];
        ptrdiff_t neighbour_j = j + j_offset[k];

        // Check if cells in bounds
        if (neighbour_i >= 0 && neighbour_i < dims[0] && neighbour_j >= 0 &&
            neighbour_j < dims[1]) {
          float horizontal_dist;
          float vertical_dist;
          float local_gradient;

          if (neighbour_i != i && neighbour_j != j) {
            horizontal_dist = SQRT2 * cellsize;
          } else {
            horizontal_dist = cellsize;
          }
          vertical_dist = fabsf(dem[neighbour_j * dims[0] + neighbour_i] -
                                dem[j * dims[0] + i]);

          local_gradient = vertical_dist / horizontal_dist;
          if (local_gradient > max_gradient) {
            max_gradient = local_gradient;
          }
        }
      }
      // convert results to correct format and save them
      switch (unit) {
        case 'r':  // radian
          output[j * dims[0] + i] = atanf(max_gradient);
          break;
        case 'd':  // degree
          output[j * dims[0] + i] = atanf(max_gradient) * (180.0f / PI);
          break;
        case 's':  // sine
          output[j * dims[0] + i] = sinf(atanf(max_gradient));
          break;
        case 'p':  // percent
          output[j * dims[0] + i] = max_gradient * 100.0f;
          break;
        default:  // tangent (default)
          output[j * dims[0] + i] = max_gradient;
          break;
      }
    }
  }
}
