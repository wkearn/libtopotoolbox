#include "reconstruct.h"

#include <stddef.h>
#include <stdint.h>

// Chris Wellons' queue
typedef struct {
  int32_t head;
  int32_t tail;
} Queue;

int32_t queue_push(Queue *q, int32_t exp) {
  int32_t mask = (1u << exp) - 1;
  int32_t head = q->head & mask;
  int32_t tail = q->tail & mask;
    
  int32_t next = (head + 1u) & mask;

  if (head & 0x80000000) { // Avoid overflow
    q->head &= ~0x80000000;
  }

  if (next == tail) {
    return -1;
  } else {
    q->head += 1;
    return head;
  }  
}

int32_t queue_pop(Queue *q, int32_t exp) {
  int32_t mask = (1u << exp) - 1;

  int32_t head = q->head & mask;
  int32_t tail = q->tail & mask;

  if (head == tail) {
    return -1;
  } else {
    q->tail += 1;
    return tail;
  }
}

/*
  Perform a partial reconstruction by scanning in the forward
  direction.

  Returns the number of pixels that were modified in the current scan.

  The forward scan replaces every pixel in `marker` with the maximum
  of a neighborhood consisting of the 4 neighbors already visited by
  the raster scan, denoted by 'x' in the diagram:

  x x .
   \|
  x-o .
   /
  x . .


  The new marker pixel value is constrained to lie below the
  corresponding pixel in `mask`.
 */
ptrdiff_t forward_scan(float *marker, float *mask, ptrdiff_t nrows,
                       ptrdiff_t ncols) {
  // Offsets for the four neighbors
  ptrdiff_t col_offset[4] = {-1, -1, -1, 0};
  ptrdiff_t row_offset[4] = {1, 0, -1, -1};

  ptrdiff_t count = 0;  // Number of modified pixels

  for (ptrdiff_t p = 0; p < nrows * ncols; p++) {
    ptrdiff_t col = p / nrows;
    ptrdiff_t row = p % nrows;

    // Compute the maximum of the marker at the current pixel and all
    // of its previously visisted neighbors
    float max_height = marker[p];
    for (ptrdiff_t neighbor = 0; neighbor < 4; neighbor++) {
      ptrdiff_t neighbor_row = row + row_offset[neighbor];
      ptrdiff_t neighbor_col = col + col_offset[neighbor];

      ptrdiff_t q = neighbor_col * nrows + neighbor_row;

      // Skip pixels outside the boundary
      if (neighbor_row < 0 || neighbor_row >= nrows || neighbor_col < 0 ||
          neighbor_col >= ncols) {
        continue;
      }

      max_height = max_height > marker[q] ? max_height : marker[q];
    }

    // Set the marker at the current pixel to the minimum of the
    // maximum height of the neighborhood and the mask at the current
    // pixel.

    float z = max_height < mask[p] ? max_height : mask[p];

    if (z != marker[p]) {
      // Increment count only if we change the current pixel
      count++;
      marker[p] = z;
    }
  }
  return count;
}

/*
  Perform a partial reconstruction by scanning in the backward
  direction.

  Returns the number of pixels that were modified in the current scan.

  The backward scan replaces every pixel in `marker` with the maximum
  of a neighborhood consisting of the 4 neighbors already visited by
  the raster scan, denoted by 'x' in the diagram:

  . . x
     /
  . o-x
    |\
  . x x


  The new marker pixel value is constrained to lie below the
  corresponding pixel in `mask`.
 */
ptrdiff_t backward_scan(float *marker, float *mask, ptrdiff_t nrows,
                        ptrdiff_t ncols) {
  // Offsets for the four neighbors
  ptrdiff_t col_offset[4] = {1, 1, 1, 0};
  ptrdiff_t row_offset[4] = {-1, 0, 1, 1};

  ptrdiff_t count = 0;  // Number of modified pixels

  // Note that the loop decreases. p must have a signed type for this
  // to work correctly.
  for (ptrdiff_t p = nrows * ncols - 1; p >= 0; p--) {
    ptrdiff_t col = p / nrows;
    ptrdiff_t row = p % nrows;

    // Compute the maximum of the marker at the current pixel and all
    // of its previously visited neighbors
    float max_height = marker[p];
    for (ptrdiff_t neighbor = 0; neighbor < 4; neighbor++) {
      ptrdiff_t neighbor_row = row + row_offset[neighbor];
      ptrdiff_t neighbor_col = col + col_offset[neighbor];
      ptrdiff_t q = neighbor_col * nrows + neighbor_row;

      // Skip pixels outside the boundary
      if (neighbor_row < 0 || neighbor_row >= nrows || neighbor_col < 0 ||
          neighbor_col >= ncols) {
        continue;
      }
      max_height = max_height > marker[q] ? max_height : marker[q];
    }

    // Set the marker at the current pixel to the minimum of the
    // maximum height of the neighborhood and the mask at the current
    // pixel.

    float z = max_height < mask[p] ? max_height : mask[p];
    if (z != marker[p]) {
      // Increment count only if we change the current pixel
      count++;
      marker[p] = z;
    }
  }
  return count;
}

/*
  Grayscale reconstruction

  Performs a grayscale reconstruction of the `mask` image by the
  `marker` image using the sequential reconstruction algorithm of Vincent
  (1993).

  Both `marker` and `mask` should point to two-dimensional arrays of
  size (nrows, ncols) with the first dimension (nrows) changing
  fastest. The `marker` array is updated with the result in-place.

  The algorithm alternately scans the marker image in the forward and
  reverse directions performing a partial reconstruction in each
  direction. It repeats these scans until no change is detected or
  until a maximum iteration threshold (currently 1000) is reached.

  Vincent, Luc. (1993). Morphological grayscale reconstruction in
  image analysis: applications and efficient algorithms. IEEE
  Transactions on Image Processing, Vol. 2, No. 2.
  https://doi.org/10.1109/83.217222
 */
void reconstruct(float *marker, float *mask, ptrdiff_t nrows, ptrdiff_t ncols) {
  ptrdiff_t n = ncols * nrows;

  const int32_t max_iterations = 1000;
  for (int32_t iteration = 0; iteration < max_iterations && n > 0;
       iteration++) {
    n = forward_scan(marker, mask, nrows, ncols);
    n += backward_scan(marker, mask, nrows, ncols);
  }
}
