#define TOPOTOOLBOX_BUILD

#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>

#include "topotoolbox.h"

// Compute the steepest descent flow direction from the DEM with flow
// routed over flat regions by the auxiliary topography in dist.
uint8_t compute_flowdirection(ptrdiff_t i, ptrdiff_t j, float *dem, float *dist,
                              int32_t *flats, ptrdiff_t dims[2]) {
  int32_t is_flat = flats[j * dims[0] + i] & 1;
  float z = is_flat > 0 ? dist[j * dims[0] + i] : dem[j * dims[0] + i];
  uint8_t direction = 0;
  float max_gradient = 0.0;

  float chamfer[8] = {1.0, sqrtf(2.0), 1.0, sqrtf(2.0),
                      1.0, sqrtf(2.0), 1.0, sqrtf(2.0)};

  ptrdiff_t i_offset[8] = {0, 1, 1, 1, 0, -1, -1, -1};
  ptrdiff_t j_offset[8] = {1, 1, 0, -1, -1, -1, 0, 1};

  for (int32_t neighbor = 0; neighbor < 8; neighbor++) {
    uint8_t query_direction = 1 << neighbor;
    // step_downstream(idxs, query_direction, idx, dims);
    ptrdiff_t neighbor_i = i + i_offset[neighbor];
    ptrdiff_t neighbor_j = j + j_offset[neighbor];

    if (neighbor_i < 0 || neighbor_i >= dims[0] || neighbor_j < 0 ||
        neighbor_j >= dims[1]) {
      continue;
    }

    float neighbor_z = dem[neighbor_j * dims[0] + neighbor_i];
    if (neighbor_z > dem[j * dims[0] + i]) {
      // This will skip any surrounding higher neighbors of flats,
      // including sills that have higher elevations
      continue;
    }

    if (is_flat > 0) {
      neighbor_z = dist[neighbor_j * dims[0] + neighbor_i];
    }

    float g = (z - neighbor_z) / chamfer[neighbor];
    if (g > max_gradient) {
      max_gradient = g;
      direction = query_direction;
    }
  }
  return direction;
}

uint8_t compute_flowdirection_TT2(ptrdiff_t i, ptrdiff_t j, float *dem,
                                  float *dist, int32_t *flats,
                                  ptrdiff_t dims[2]) {
  uint8_t direction = 0;

  float z = dem[j * dims[0] + i];
  float max_gradient = 0.0;
  float Ggwdt = dist[j * dims[0] + i];
  float chamfer[8] = {1.0, sqrtf(2.0), 1.0, sqrtf(2.0),
                      1.0, sqrtf(2.0), 1.0, sqrtf(2.0)};

  int32_t is_flat = flats[j * dims[0] + i] & 1;

  ptrdiff_t i_offset[8] = {0, 1, 1, 1, 0, -1, -1, -1};
  ptrdiff_t j_offset[8] = {1, 1, 0, -1, -1, -1, 0, 1};

  for (int32_t neighbor = 0; neighbor < 8; neighbor++) {
    uint8_t query_direction = 1 << neighbor;

    ptrdiff_t neighbor_i = i + i_offset[neighbor];
    ptrdiff_t neighbor_j = j + j_offset[neighbor];

    if (neighbor_i < 0 || neighbor_i >= dims[0] || neighbor_j < 0 ||
        neighbor_j >= dims[1]) {
      continue;
    }
    float neighbor_z = dem[neighbor_j * dims[0] + neighbor_i];

    if (neighbor_z > z) {
      continue;
    }

    if (is_flat) {
      if (flats[neighbor_j * dims[0] + neighbor_i] & 2) {
        // neighbor is a sill
        direction = query_direction;
        break;
      } else {
        float ggwdt = dist[neighbor_j * dims[0] + neighbor_i];
        if (ggwdt < Ggwdt) {
          direction = query_direction;
          Ggwdt = ggwdt;
        }
      }
    } else {
      // Pixel is not flat
      float g = (z - neighbor_z) / chamfer[neighbor];
      if (g > max_gradient) {
        max_gradient = g;
        direction = query_direction;
      }
    }
  }
  return direction;
}

TOPOTOOLBOX_API
void flow_routing_d8_carve(ptrdiff_t *node, uint8_t *direction, float *dem,
                           float *dist, int32_t *flats, ptrdiff_t dims[2]) {
  // node contains an array of dims[0] * dims[1] linear pixel indices
  // into dem. These indices are sorted topologically, so that if
  // there is an edge from node u to node v, u comes before v in the
  // array.
  //
  // direction[i] is the bitfield-encoded flow direction for the
  // pixel at i
  //
  // 1<<5 1<<6  1<<7
  // 1<<4    0  1<<0
  // 1<<3 1<<2  1<<1
  //
  // To construct a topological ordering of the flow graph, we conduct
  // a depth first traversal of the flow graph. The topological order
  // is given by a reversed postorder of the nodes encountered during
  // the traversal. As a result, the node array fills up from the
  // bottom.
  //
  // Use node[--next] = u to append vertex u onto the node list.
  ptrdiff_t next = dims[0] * dims[1];

  // To implement the depth first traversal, we need a stack. The
  // stack can only hold up to as many vertices as have not yet been
  // assigned to the node list.  We can therefore using the top of the
  // node array to hold our stack.
  //
  // Use source[stack_top++] = u to push vertex u onto the stack
  // Use idx = source[--stack_top] top pop the stack
  ptrdiff_t stack_top = 0;

  // Initialize the directions with 0xff, which is our sentinel value
  // for unvisited nodes
  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      direction[j * dims[0] + i] = 0xff;
    }
  }

  ptrdiff_t i_offset[8] = {0, 1, 1, 1, 0, -1, -1, -1};
  ptrdiff_t j_offset[8] = {1, 1, 0, -1, -1, -1, 0, 1};

  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      if (direction[j * dims[0] + i] == 0xff) {
        // Pixel is unvisited, push it onto the stack and do a
        // depth-first search.
        // The stack should always be empty by the time we get here
        assert(stack_top == 0);
        node[stack_top++] = j * dims[0] + i;

        while (stack_top > 0) {
          // Pop a node from the stack
          ptrdiff_t u = node[--stack_top];

          if (direction[u] == 0xff) {
            // node has not been discovered yet

            // Compute its Cartesian indices
            ptrdiff_t u_i = u % dims[0];
            ptrdiff_t u_j = u / dims[0];

            uint8_t flowdir =
                compute_flowdirection(u_i, u_j, dem, dist, flats, dims);

            if (flowdir == 0) {
              // This node is a sink/outlet
              direction[u] = flowdir;
              // Prepend the node to the edge list.
              assert(next > 0);
              assert(next > stack_top);
              node[--next] = u;
              continue;
            }
            // flowdir is not an index, but 1<<index
            // Compute the index
            uint8_t d = flowdir;
            uint8_t r = 0;
            while (d >>= 1) {
              r++;
            }
            ptrdiff_t v_i = u_i + i_offset[r];
            ptrdiff_t v_j = u_j + j_offset[r];

            if (v_i < 0 || v_i >= dims[0] || v_j < 0 || v_j >= dims[1]) {
              continue;
            }

            ptrdiff_t v = v_j * dims[0] + v_i;

            if (direction[v] == 0xff) {
              // The neighbor has not been visited yet
              // Push the node back on the stack, then push the neighbor
              assert(stack_top < next);
              node[stack_top++] = u;

              assert(stack_top < next);
              node[stack_top++] = v;
            } else {
              // The downstream neighbor has been visited, visit the node
              direction[u] = flowdir;
              // Prepend the node to the edge list.
              assert(next > 0);
              assert(next > stack_top);
              node[--next] = u;
            }
          } else {
            // This is a cycle
            assert(0);
          }
        }
      }
    }
  }
}

TOPOTOOLBOX_API
ptrdiff_t flow_routing_d8_edgelist(ptrdiff_t *source, ptrdiff_t *target,
                                   ptrdiff_t *node, uint8_t *direction,
                                   ptrdiff_t dims[2]) {
  ptrdiff_t offsets[8] = {dims[0],  dims[0] + 1,  1,  -dims[0] + 1,
                          -dims[0], -dims[0] - 1, -1, dims[0] - 1};

  ptrdiff_t edge_count = 0;
  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      ptrdiff_t u = node[j * dims[0] + i];

      uint8_t flowdir = direction[u];

      if (flowdir != 0) {
        uint8_t v = flowdir;
        uint8_t r = 0;
        while (v >>= 1) {
          r++;
        }

        source[edge_count] = u;
        target[edge_count++] = u + offsets[r];
      }
    }
  }
  return edge_count;
}
