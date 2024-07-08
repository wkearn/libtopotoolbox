#ifndef TOPOTOOLBOX_MORPHOLOGY_RECONSTRUCT_H
#define TOPOTOOLBOX_MORPHOLOGY_RECONSTRUCT_H

#include <stddef.h>

void reconstruct(float *marker, float *mask, ptrdiff_t dims[2],
                 ptrdiff_t strides[2]);

#endif  // TOPOTOOLBOX_MORPHOLOGY_RECONSTRUCT_H
