/**
   @file topotoolbox.h
   @author William Kearney <william.kearney@uni-potsdam.de>
   @version 3.0
 
   @brief Public API for libtopotoolbox
 */
#ifndef TOPOTOOLBOX_H
#define TOPOTOOLBOX_H

#ifdef TOPOTOOLBOX_STATICLIB
#define TOPOTOOLBOX_API
#elif defined(TOPOTOOLBOX_BUILD)
#if defined(_WIN32)
#define TOPOTOOLBOX_API __declspec(dllexport)
#else
#define TOPOTOOLBOX_API
#endif
#else
#if defined(_WIN32)
#define TOPOTOOLBOX_API __declspec(dllimport)
#else
#define TOPOTOOLBOX_API
#endif
#endif

#define TOPOTOOLBOX_VERSION_MAJOR 3
#define TOPOTOOLBOX_VERSION_MINOR 0
#define TOPOTOOLBOX_VERSION_PATCH 0

#include <stddef.h>
#include <stdint.h>

/**
   Used to ensure that topotoolbox is compiled and linked
   correctly. Always returns 1.
 */
TOPOTOOLBOX_API
int has_topotoolbox(void);

/**
   @brief Fills sinks in a digital elevation model

   Uses an algorithm based on grayscale morphological reconstruction.

   @param[out] output The filled DEM
   @param[in]  dem    The input DEM
   @param[in]  nrows  The size of both DEMs in the fastest changing dimension
   @param[in]  ncols  The size of both DEMs in the slowest changing dimension
 */
TOPOTOOLBOX_API
void fillsinks(float *output, float *dem, ptrdiff_t nrows, ptrdiff_t ncols);


/**
   @brief Labels flat, sill and presill pixels in the provided DEM

   A flat pixel is one surrounded by pixels with the same or higher
   elevations. A sill pixel has the same elevation as a neighboring
   flat pixel but borders a pixel with a lower elevation. A presill
   pixel is a flat pixel that borders a sill pixel.

   The pixels are labeled with a bit field:

   - Bit 0: Set if pixel is a flat
   - Bit 1: Set if pixel is a sill
   - Bit 2: Set if pixel is a presill

   Since all presill pixels are also flats, bits 0 and 2 are both set
   for presill pixels. In other words presill pixels have a value of 5
   in the output array. This allows one to test the identity of pixels
   with bitwise operations:

   ```
   if (output[pixel] & 1) {
     // Pixel is a flat
   }
   if (output[pixel] & 2) {
     // Pixel is a sill
   }
   if (output[pixel] & 4) {
     // Pixel is a presill
   }
   ```

   @param[out] output The integer-valued output array with pixel labels
   @param[in]  dem    The input DEM
   @param[in]  nrows  The size of both DEMs in the fastest changing dimension
   @param[in]  ncols  The size of both DEMs in the slowest changing dimension
 */
TOPOTOOLBOX_API
ptrdiff_t identifyflats(int32_t *output, float *dem, ptrdiff_t nrows,
                        ptrdiff_t ncols);

/**
  @brief Compute costs for the gray-weighted distance transform

  The costs used to route flow over flat regions with the
  gray-weighted distance transform are based on the difference between
  original and filled DEMs at each pixel. This difference is
  subtracted from the maximum difference over the flat region to which
  the pixel belongs. It is then squared and a small constant
  (currently set to 0.1) is added.

  This function identifies the connected components of flat pixels,
  provided in the `flats` array, computes the maximum difference
  between the `original_dem` and the `filled_dem` over each connected
  component and then computes the cost for each flat pixel. The
  `costs` are returned, as are the connected components labels
  (`conncomps`). These labels are the linear index of the pixel in the
  connected component with the maximum difference between the filled
  and the output DEMs.

  @param[out] costs        The costs for the gray-weighted distance transform
  @param[out] conncomps    Labeled connected components
  @param[in]  flats        Array identifying the flat pixels
  @param[in]  original_dem The DEM prior to sink filling
  @param[in]  filled_dem   The DEM after sink filling
  @param[in]  nrows  The size of both DEMs in the fastest changing dimension
  @param[in]  ncols  The size of both DEMs in the slowest changing dimension
 */
TOPOTOOLBOX_API
void gwdt_computecosts(float *costs, ptrdiff_t *conncomps, int32_t *flats,
                       float *original_dem, float *filled_dem, ptrdiff_t nrows,
                       ptrdiff_t ncols);

#endif  // TOPOTOOLBOX_H
