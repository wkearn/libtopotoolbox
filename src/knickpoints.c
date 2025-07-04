#define TOPOTOOLBOX_BUILD

#include <math.h>

#include "topotoolbox.h"

ptrdiff_t nnz(ptrdiff_t *ix, ptrdiff_t node_count) {
  ptrdiff_t n = 0;
  for (ptrdiff_t i = 0; i < node_count; i++) {
    n += ix[i] ? 1 : 0;
  }
  return n;
}

TOPOTOOLBOX_API
void lowerenv(float *elevation, uint8_t *knickpoints, float *distance,
              ptrdiff_t *ix, uint8_t *onenvelope, ptrdiff_t *source,
              ptrdiff_t *target, ptrdiff_t edge_count, ptrdiff_t node_count) {
  // All nodes start out on the envelope
  for (ptrdiff_t i = 0; i < node_count; i++) {
    onenvelope[i] = 1;
  }

  for (ptrdiff_t e = edge_count - 1; e >= 0; e--) {
    ptrdiff_t u = source[e];
    ptrdiff_t v = target[e];

    if (onenvelope[u]) {
      // allpred
      for (ptrdiff_t i = 0; i < node_count; i++) {
        ix[i] = 0;
      }
      ix[u] = 1;
      for (ptrdiff_t e1 = edge_count - 1; e1 >= 0; e1--) {
        ptrdiff_t u1 = source[e1];
        ptrdiff_t v1 = target[e1];

        ix[u1] = ix[u1] | (ix[v1] & ~knickpoints[v1]);
      }
      // end allpred

      if (nnz(ix, node_count) == 0) {
        // Skip gradient computations if there are no upstream nodes
        // of u below a knickpoint.
        continue;
      }

      // Compute the minimum gradient between v and a point upstream
      // of v but below any knickpoints.
      float g = INFINITY;
      ptrdiff_t idx = -1;
      for (ptrdiff_t e1 = 0; e1 < edge_count; e1++) {
        ptrdiff_t u1 = source[e1];

        if (ix[u1]) {
          float gu =
              (elevation[u1] - elevation[v]) / (distance[u1] - distance[v]);
          if (gu < g) {
            g = gu;
            idx = u1;
          }
        }
      }

      // Set up fast indexing/tracking visited nodes
      //
      // MATLAB uses a separate ixcix array, but we will repurpose the
      // ix array from above, which is no longer required.
      //
      // Initialize ix to -1. If ix[u] = -1, u has been visited or u
      // has no downstream neighbors.
      for (ptrdiff_t i = 0; i < node_count; i++) {
        ix[i] = -1;
      }

      // ix[u] is the index in the edge lists of the unique edge that
      // starts at u.
      for (ptrdiff_t e1 = 0; e1 < edge_count; e1++) {
        ix[source[e1]] = e1;
      }

      // Visit node v
      ix[v] = -1;

      float z0 = elevation[idx];
      float d0 = distance[idx];
      // Loop until we hit a node that has been visited.
      while (ix[idx] >= 0) {
        ptrdiff_t idx2 = target[ix[idx]];

        // Adjust the distance downward by the minimum gradient we found earlier
        elevation[idx2] = z0 - g * (d0 - distance[idx2]);

        // Take idx2 off the envelope, because its elevation has been
        // updated. It will not be searched in later iterations of the
        // outer loop, though it may still be updated.
        onenvelope[idx2] = 0;

        // Mark idx as visited
        ix[idx] = -1;

        // Visit the downstream neighbor
        idx = idx2;
      }
    }
  }
}
