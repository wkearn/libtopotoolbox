#include <math.h>
#include <stddef.h>
#include <stdint.h>

#include "priority_queue.h"
#include "topotoolbox.h"

float eikonal_update(float *solution, float fi, ptrdiff_t i, ptrdiff_t j,
                     ptrdiff_t nrows, ptrdiff_t ncols) {
  float south = i < nrows - 1 ? solution[j * nrows + i + 1] : INFINITY;
  float north = i > 0 ? solution[j * nrows + i - 1] : INFINITY;

  float u1 = fminf(south, north);

  float east = j < ncols - 1 ? solution[(j + 1) * nrows + i] : INFINITY;
  float west = j > 0 ? solution[(j - 1) * nrows + i] : INFINITY;

  float u2 = fminf(east, west);

  if (fabsf(u1 - u2) < fi) {
    return (u1 + u2) / 2 +
           sqrtf((u1 + u2) * (u1 + u2) - 2 * (u1 * u1 + u2 * u2 - fi * fi)) / 2;
  } else {
    return fminf(u1 + fi, u2 + fi);
  }
}

TOPOTOOLBOX_API
void excesstopography(float *excess, ptrdiff_t *heap, ptrdiff_t *back,
                         float *dem, float *threshold, float cellsize,
                         ptrdiff_t nrows, ptrdiff_t ncols) {
  // Initialize the arrays
  // Pixels start with the elevation given by the DEM
  for (ptrdiff_t j = 0; j < ncols; j++) {
    for (ptrdiff_t i = 0; i < nrows; i++) {
      ptrdiff_t idx = j * nrows + i;
      back[idx] = idx;
      excess[idx] = dem[idx];
      heap[idx] = idx;
    }
  }
  ptrdiff_t count = ncols * nrows;
  PriorityQueue q = pq_create(count, heap, back, excess, 1);

  while (!pq_isempty(&q)) {
    ptrdiff_t trial = pq_deletemin(&q);
    float trial_elevation = pq_get_priority(&q, trial);

    ptrdiff_t col = trial / nrows;
    ptrdiff_t row = trial % nrows;

    // Neighbors are only visited if they are not boundary pixels

    // South neighbor
    if (row < nrows-1 && excess[col * nrows + row + 1] >= trial_elevation) {
      float f = cellsize * threshold[col * nrows + row + 1];
      float proposal =
          eikonal_update(q.priorities, f, row + 1, col, nrows, ncols);
      pq_decrease_key(&q, col * nrows + row + 1, proposal);
    }

    // North neighbor
    if (row > 0 && excess[col * nrows + row - 1] >= trial_elevation) {
      float f = cellsize * threshold[col * nrows + row - 1];
      float proposal =
          eikonal_update(q.priorities, f, row - 1, col, nrows, ncols);
      pq_decrease_key(&q, col * nrows + row - 1, proposal);
    }

    // East neighbor
    if (col < ncols - 1 && excess[(col + 1) * nrows + row] >= trial_elevation) {
      float f = cellsize * threshold[(col + 1) * nrows + row];
      float proposal =
          eikonal_update(q.priorities, f, row, col + 1, nrows, ncols);
      pq_decrease_key(&q, (col + 1) * nrows + row, proposal);
    }

    // West neighbor
    if (col > 0 && excess[(col - 1) * nrows + row] >= trial_elevation) {
      float f = cellsize * threshold[(col - 1) * nrows + row];
      float proposal =
          eikonal_update(q.priorities, f, row, col - 1, nrows, ncols);
      pq_decrease_key(&q, (col - 1) * nrows + row, proposal);
    }
  }
}
