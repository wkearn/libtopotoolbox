#include <math.h>

#include "topotoolbox.h"

#define PI 3.14159265358979323846f
#define PI_2 1.57079632679489661923f

TOPOTOOLBOX_API
void gradient_secondorder(float *restrict p0, float *restrict p1, float *dem,
                          float cellsize, ptrdiff_t dims[2]) {
  // Compute second order
  ptrdiff_t m = dims[0];
  ptrdiff_t n = dims[1];

  // Boundaries use a stencil of [-3 4 -1]/(2*cs) and [1 -4 3]/(2*cs)
  // Interior points use a stencil of [-1 0 1]/(2*cs)

  // Compute gradient in the first dimension
  for (ptrdiff_t j = 0; j < n; j++) {
    ptrdiff_t i = 0;

    p0[j * m + i] =
        -dem[j * m + i + 2] + 4 * dem[j * m + i + 1] - 3 * dem[j * m + i];
    p0[j * m + i] /= 2 * cellsize;

    for (i = 1; i < m - 1; i++) {
      p0[j * m + i] =
          (dem[j * m + i + 1] - dem[j * m + i - 1]) / (2 * cellsize);
    }

    i = m - 1;
    p0[j * m + i] =
        dem[j * m + i - 2] - 4 * dem[j * m + i - 1] + 3 * dem[j * m + i];
    p0[j * m + i] /= 2 * cellsize;
  }

  // Compute gradient in the second dimension
  ptrdiff_t j = 0;
  for (ptrdiff_t i = 0; i < m; i++) {
    p1[j * m + i] = -dem[(j + 2) * m + i] + 4 * dem[(j + 1) * m + i] -
                    3 * dem[(j + 0) * m + i];
    p1[j * m + i] /= 2 * cellsize;
  }

  for (j = 1; j < n - 1; j++) {
    for (ptrdiff_t i = 0; i < m; i++) {
      p1[j * m + i] =
          (dem[(j + 1) * m + i] - dem[(j - 1) * m + i]) / (2 * cellsize);
    }
  }

  j = n - 1;
  for (ptrdiff_t i = 0; i < m; i++) {
    p1[j * m + i] =
        dem[(j - 2) * m + i] - 4 * dem[(j - 1) * m + i] + 3 * dem[j * m + i];
    p1[j * m + i] /= 2 * cellsize;
  }
}

void normal_vectors(float *restrict nx, float *restrict ny, float *restrict nz,
                    float *dem, float cellsize, ptrdiff_t dims[2]) {
  gradient_secondorder(nx, ny, dem, cellsize, dims);

  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      float dx = nx[j * dims[0] + i];
      float dy = ny[j * dims[0] + i];

      float inorm = 1.0f / sqrtf(dx * dx + dy * dy + 1.0);

      nx[j * dims[0] + i] *= -inorm;
      ny[j * dims[0] + i] *= -inorm;
      nz[j * dims[0] + i] = inorm;
    }
  }
}

void hillshade(float *output, float *nx, float *ny, float *nz, float *dem,
               float azimuth, float altitude, float cellsize,
               ptrdiff_t dims[2]) {
  normal_vectors(nx, ny, nz, dem, cellsize, dims);

  float sx = sinf(PI_2 - altitude) * cosf(azimuth);
  float sy = sinf(PI_2 - altitude) * sinf(azimuth);
  float sz = cosf(PI_2 - altitude);

  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      float x = nx[j * dims[0] + i];
      float y = ny[j * dims[0] + i];
      float z = nz[j * dims[0] + i];
      output[j * dims[0] + i] = x * sx + y * sy + z * sz;
    }
  }
}
