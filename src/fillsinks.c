#define TOPOTOOLBOX_BUILD

#include <math.h>
#include <stddef.h>
#include <stdint.h>

#include "morphology/reconstruct.h"
#include "topotoolbox.h"

/*
  Fill sinks in a digital elevation model.

  The arrays pointed to by `output` and `dem` should represent a two
  dimensional array of size (nrows, ncols). The filled DEM is written
  to the `output` array.

  `dem` is modified during the processing, but the modifications are
  reverted before the function returns. Be careful when accessing
  `dem` concurrently with `fillsinks`.

  Sinks are filled using grayscale morphological reconstruction.
*/
TOPOTOOLBOX_API
void fillsinks(float *output, float *dem, ptrdiff_t dims[2], ptrdiff_t strides[2]) { 
  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      // Invert the DEM
      dem[i * strides[0] + j * strides[1]] *= -1.0;

      // Set the boundary pixels of the output equal to the DEM and
      // the interior pixels equal to -INFINITY.
      if ((i == 0 || i == (dims[0] - 1)) ||
          (j == 0 || j == (dims[1] - 1))) {
        output[i * strides[0] + j * strides[1]] = dem[i * strides[0] + j * strides[1]];
      } else {
        output[i * strides[0] + j * strides[1]] = -INFINITY;
      }
    }
  }

  reconstruct(output, dem, dims, strides);

  // Revert the DEM and the output
  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      dem[i * strides[0] + j * strides[1]] *= -1.0;
      output[i * strides[0] + j * strides[1]] *= -1.0;
    }
  }
}
