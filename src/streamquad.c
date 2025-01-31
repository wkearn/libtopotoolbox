#define TOPOTOOLBOX_BUILD

#include <math.h>
#include <stddef.h>
#include <stdint.h>

#include "topotoolbox.h"

TOPOTOOLBOX_API
void streamquad_trapz_f32(float *integral, float *integrand, ptrdiff_t *source,
                          ptrdiff_t *target, float *weight,
                          ptrdiff_t edge_count) {
  // Iterate over the edges in reverse topological order (upstream)
  for (ptrdiff_t e = edge_count - 1; e >= 0; e--) {
    // NOTE: source and target are 1-based indices into the node attribute
    // list, so we must subtract 1 here
    ptrdiff_t u = source[e];
    ptrdiff_t v = target[e];
    integral[u] = integral[v] + weight[e] * (integrand[u] + integrand[v]) / 2;
  }
}

TOPOTOOLBOX_API
void streamquad_trapz_f64(double *integral, double *integrand,
                          ptrdiff_t *source, ptrdiff_t *target, float *weight,
                          ptrdiff_t edge_count) {
  // Iterate over the edges in reverse topological order (upstream)
  for (ptrdiff_t e = edge_count - 1; e >= 0; e--) {
    // NOTE: source and target are 1-based indices into the node attribute
    // list, so we must subtract 1 here
    ptrdiff_t u = source[e];
    ptrdiff_t v = target[e];
    integral[u] = integral[v] + weight[e] * (integrand[u] + integrand[v]) / 2;
  }
}

TOPOTOOLBOX_API
void traverse_up_u32_and(uint32_t *output, uint32_t *input, ptrdiff_t *source,
                         ptrdiff_t *target, ptrdiff_t edge_count) {
  for (ptrdiff_t e = edge_count - 1; e >= 0; e--) {
    ptrdiff_t u = source[e];
    ptrdiff_t v = target[e];

    output[u] = output[v] & input[u];
  }
}

TOPOTOOLBOX_API
void traverse_down_f32_max_add(float *output, float *input, ptrdiff_t *source,
                               ptrdiff_t *target, ptrdiff_t edge_count) {
  for (ptrdiff_t e = 0; e < edge_count; e++) {
    ptrdiff_t u = source[e];
    ptrdiff_t v = target[e];

    output[v] = fmaxf(output[v], output[u] + input[e]);
  }
}
