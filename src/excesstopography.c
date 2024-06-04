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
void fsm_excesstopography(float *excess, float *dem, float *threshold,
                          float cellsize, ptrdiff_t nrows, ptrdiff_t ncols) {
  // Initialize excess == dem
  for (ptrdiff_t col = 0; col < ncols; col++) {
    for (ptrdiff_t row = 0; row < nrows; row++) {
      excess[col * nrows + row] = dem[col * nrows + row];
    }
  }
  ptrdiff_t count = nrows * ncols;

  while (count > 0) {
    count = 0;
    // Sweep 1
    for (ptrdiff_t col = 1; col < ncols - 1; col++) {
      for (ptrdiff_t row = 1; row < nrows - 1; row++) {
        float fi = cellsize * threshold[col * nrows + row];
        float proposal = eikonal_update(excess, fi, row, col, nrows, ncols);
        if (proposal < excess[col * nrows + row]) {
          excess[col * nrows + row] = proposal;
          count += 1;
        }
      }
    }

    for (ptrdiff_t col = ncols - 2; col > 0; col--) {
      for (ptrdiff_t row = 1; row < nrows - 1; row++) {
        float fi = cellsize * threshold[col * nrows + row];
        float proposal = eikonal_update(excess, fi, row, col, nrows, ncols);
        if (proposal < excess[col * nrows + row]) {
          excess[col * nrows + row] = proposal;
          count += 1;
        }
      }
    }

    // Sweep 3
    for (ptrdiff_t col = ncols - 2; col > 0; col--) {
      for (ptrdiff_t row = nrows - 2; row > 0; row--) {
        float fi = cellsize * threshold[col * nrows + row];
        float proposal = eikonal_update(excess, fi, row, col, nrows, ncols);
        if (proposal < excess[col * nrows + row]) {
          excess[col * nrows + row] = proposal;
          count += 1;
        }
      }
    }

    // Sweep 4
    for (ptrdiff_t col = 1; col < ncols - 1; col++) {
      for (ptrdiff_t row = nrows - 2; row > 0; row--) {
        float fi = cellsize * threshold[col * nrows + row];
        float proposal = eikonal_update(excess, fi, row, col, nrows, ncols);
        if (proposal < excess[col * nrows + row]) {
          excess[col * nrows + row] = proposal;
          count += 1;
        }
      }
    }
  }
}

TOPOTOOLBOX_API
void fmm_excesstopography(float *excess, ptrdiff_t *heap, ptrdiff_t *back,
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
    if (row < nrows - 1 && excess[col * nrows + row + 1] >= trial_elevation) {
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

TOPOTOOLBOX_API
void fmm_excesstopography3d(float *excess, ptrdiff_t *heap, ptrdiff_t *back,
                            float *dem, float *lithstack,
                            float *threshold_slopes, float cellsize,
                            ptrdiff_t nrows, ptrdiff_t ncols,
                            ptrdiff_t nlayers) {
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

    // South neighbor
    if (row < nrows - 1 &&
        pq_get_priority(&q, col * nrows + row + 1) >= trial_elevation) {
      float proposal = pq_get_priority(&q, col * nrows + row + 1);
      for (int layer = 0; layer < nlayers; layer++) {
        float fi = cellsize * threshold_slopes[layer];
        proposal = eikonal_update(q.priorities, fi, row + 1, col, nrows, ncols);
        if (proposal < lithstack[(col * nrows + row + 1) * nlayers + layer]) {
          pq_decrease_key(&q, col * nrows + row + 1, proposal);
          break;
        }
      }
    }

    // North neighbor
    if (row > 0 &&
        pq_get_priority(&q, col * nrows + row - 1) >= trial_elevation) {
      float proposal = pq_get_priority(&q, col * nrows + row - 1);
      for (int layer = 0; layer < nlayers; layer++) {
        float fi = cellsize * threshold_slopes[layer];
        proposal = eikonal_update(q.priorities, fi, row - 1, col, nrows, ncols);
        if (proposal < lithstack[(col * nrows + row - 1) * nlayers + layer]) {
          pq_decrease_key(&q, col * nrows + row - 1, proposal);
          break;
        }
      }
    }

    // East neighbor
    if (col < ncols - 1 &&
        pq_get_priority(&q, (col + 1) * nrows + row) >= trial_elevation) {
      float proposal = pq_get_priority(&q, (col + 1) * nrows + row);

      for (int layer = 0; layer < nlayers; layer++) {
        float fi = cellsize * threshold_slopes[layer];
        proposal = eikonal_update(q.priorities, fi, row, col + 1, nrows, ncols);
        if (proposal < lithstack[((col + 1) * nrows + row) * nlayers + layer]) {
          pq_decrease_key(&q, (col + 1) * nrows + row, proposal);
          break;
        }
      }
    }

    // West neighbor
    if (col > 0 &&
        pq_get_priority(&q, (col - 1) * nrows + row) >= trial_elevation) {
      float proposal = pq_get_priority(&q, (col - 1) * nrows + row);
      for (int layer = 0; layer < nlayers; layer++) {
        float fi = cellsize * threshold_slopes[layer];
        proposal = eikonal_update(q.priorities, fi, row, col - 1, nrows, ncols);
        if (proposal < lithstack[((col - 1) * nrows + row) * nlayers + layer]) {
          pq_decrease_key(&q, (col - 1) * nrows + row, proposal);
          break;
        }
      }
    }
  }
}
