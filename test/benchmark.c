#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>

#include "topotoolbox.h"

/*
  PCG4D hash function


  Jarzynski, Mark and Olano, Marc. (2020). Hash functions for GPU
  rendering. Journal of Computer Graphics Techniques. Vol 9,
  No. 3. 21-38.
 */
float pcg4d(uint32_t a, uint32_t b, uint32_t c, uint32_t d) {
  uint32_t x = a * 1664525u + 1013904223u;
  uint32_t y = b * 1664525u + 1013904223u;
  uint32_t z = c * 1664525u + 1013904223u;
  uint32_t w = d * 1664525u + 1013904223u;

  x += y * w;
  y += z * x;
  z += x * y;
  w += y * z;

  x ^= x >> 16;
  y ^= y >> 16;
  z ^= z >> 16;
  w ^= w >> 16;

  x += y * w;
  y += z * x;
  z += x * y;
  w += y * z;

  return (float)(w >> 8) / (1 << 24);
}

int main(int argc, char *argv[]) {
  ptrdiff_t nrows = 200;
  ptrdiff_t ncols = 400;

  float *original_dem = malloc(sizeof(float)*nrows*ncols);

  // Allocate output for fillsinks
  float *filled_dem = malloc(sizeof(float)*nrows*ncols);

  // Allocate output for identify flats
  int32_t *flats = malloc(sizeof(int32_t)*nrows*ncols);

  // Allocate output for compute_costs
  ptrdiff_t *labels = malloc(sizeof(ptrdiff_t)*nrows*ncols);
  float *costs = malloc(sizeof(float)*nrows * ncols);

  if (!original_dem || !filled_dem || !flats || !labels || !costs) {
    printf("Failed to allocate memory");
    return -1;
  }
  
  double duration = 0.0;
  for (uint32_t test = 0; test < 100; test++) {

    for (uint32_t col = 0; col < ncols; col++) {
      for (uint32_t row = 0; row < nrows; row++) {
        original_dem[col * nrows + row] = 100.0f * pcg4d(row, col, test, 1);
      }
    }
    
    // Run and time analysis
    clock_t t1 = clock();
    fillsinks(filled_dem, original_dem, nrows, ncols);
    identifyflats(flats,filled_dem,nrows,ncols);
    compute_costs(costs,labels,flats,original_dem,filled_dem,nrows,ncols);
    clock_t t2 = clock();
    duration += 1000.0 * (t2 - t1) / CLOCKS_PER_SEC;
  }
  printf("Duration: %f ms\n",duration);
  return 0; 
}
