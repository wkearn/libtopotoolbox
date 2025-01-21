#define TOPOTOOLBOX_BUILD

#include "topotoolbox.h"

TOPOTOOLBOX_API
void drainagebasins(ptrdiff_t *basins, ptrdiff_t *source, ptrdiff_t *target,
                    ptrdiff_t dims[2]) {
  // Initialize the basins array
  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      basins[j * dims[0] + i] = 0;
    }
  }

  ptrdiff_t basin_count = 1;  // Start the basin labels at one.

  ptrdiff_t edge_count = dims[0] * dims[1];

  for (ptrdiff_t e = edge_count - 1; e >= 0; e--) {
    ptrdiff_t src = source[e];
    ptrdiff_t tgt = target[e];

    if (tgt < 0) {
      // src is an outlet
      //
      // Don't initialize the basin here: we are refactoring to get
      // rid of these placeholder edges
      continue;
    }
    if (basins[tgt] == 0) {
      // If tgt does not yet belong to a drainage basin, create a new one
      basins[tgt] = basin_count++;
    }
    basins[src] = basins[tgt];
  }
}
