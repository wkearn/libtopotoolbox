#define TOPOTOOLBOX_BUILD

#include <stddef.h>
#include <stdint.h>
#include <math.h>

#include "topotoolbox.h"

/*
  Compute the gradient of each cell in a digital elevation model.

  The arrays pointed to by `gradient` and `dem` should represent
  two-dimensional arrays of size (nrows, ncols). The gradient is 
  computed as the maximum difference in elevation between a cell 
  and its 8 neighbors.


  char unit:
    't' --> tangent (default)
    'r' --> radian
    'd' --> degree
    's' --> sine
    'p' --> percent

 */

TOPOTOOLBOX_API
void gradient8(float *output, float *dem, ptrdiff_t nrows,
                     ptrdiff_t ncols, char unit) {
  ptrdiff_t col_offset[8] = {-1, -1, -1, 0, 0, 1, 1, 1};
  ptrdiff_t row_offset[8] = {-1, 0, 1, -1, 1, -1, 0, 1};

  // Loop through each cell in the DEM
  for (ptrdiff_t col = 0; col < ncols; col++) {
    for (ptrdiff_t row = 0; row < nrows; row++) {
      // Skip border pixels
      if (col == 0 || col == ncols - 1 || row == 0 || row == nrows - 1) {
        output[col * nrows + row] = 0;
        continue;
      }

      float dem_height = dem[col * nrows + row];
      float max_diff = 0;

      // Compute the maximum height difference in the neighborhood around the current pixel
      for (int32_t neighbor = 0; neighbor < 8; neighbor++) {
        ptrdiff_t neighbor_row = row + row_offset[neighbor];
        ptrdiff_t neighbor_col = col + col_offset[neighbor];

        // neighbor_row and neighbor_col are valid indices because we skipped border pixels above
        float neighbor_height = dem[neighbor_col * nrows + neighbor_row];
        float height_diff = fabs(dem_height - neighbor_height);
        if (height_diff > max_diff) {
          max_diff = height_diff;
        }
      }

      output[col * nrows + row] = max_diff;
    }
  }
}
