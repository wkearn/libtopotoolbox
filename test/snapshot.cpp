#include <gdal.h>
#undef NDEBUG

#include <gdal_priv.h>

#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <ctime>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

namespace tt {
extern "C" {
#include "topotoolbox.h"
}
}  // namespace tt

#include "profiler.h"
Profiler prof;

template <typename T, GDALDataType S>
void load_data_from_file(const std::string& filename, std::vector<T>& output,
                         std::array<ptrdiff_t, 2>& dims) {
  GDALDataset* dataset = GDALDataset::Open(filename.data(), GA_ReadOnly);

  GDALRasterBand* poBand = dataset->GetRasterBand(1);
  dims[0] = poBand->GetXSize();
  dims[1] = poBand->GetYSize();

  output.resize(dims[0] * dims[1]);

  assert(poBand->RasterIO(GF_Read, 0, 0, dims[0], dims[1], output.data(),
                          dims[0], dims[1], S, 0, 0) == CE_None);

  GDALClose(dataset);
}

struct SnapshotData {
  std::array<ptrdiff_t, 2> dims;
  float cellsize;

  std::vector<float> dem;
  std::vector<float> filled_dem;
  std::vector<int32_t> flats;
  std::vector<int32_t> sills;

  // Output arrays
  std::vector<float> test_dem;
  std::vector<float> test_filled_dem;
  std::vector<int32_t> test_flats;
  std::vector<int32_t> test_sills;

  // Intermediate arrays
  std::vector<uint8_t> bc;

  SnapshotData(const std::filesystem::path& snapshot_path) {
    GDALAllRegister();

    if (exists(snapshot_path / "dem.tif")) {
      load_data_from_file<float, GDT_Float32>(snapshot_path / "dem.tif", dem,
                                              dims);
      test_dem = dem;  // Create a copy of the DEM to modify in fillsinks
      assert(test_dem.size() == dem.size());
    }

    std::array<ptrdiff_t, 2> dims_check;
    if (exists(snapshot_path / "fillsinks.tif")) {
      load_data_from_file<float, GDT_Float32>(snapshot_path / "fillsinks.tif",
                                              filled_dem, dims_check);
      assert(dims[0] == dims_check[0] && dims[1] == dims_check[1]);
    }

    if (exists(snapshot_path / "identifyflats_flats.tif")) {
      load_data_from_file<int32_t, GDT_Int32>(
          snapshot_path / "identifyflats_flats.tif", flats, dims_check);
      assert(dims[0] == dims_check[0] && dims[1] == dims_check[1]);
    }

    if (exists(snapshot_path / "identifyflats_sills.tif")) {
      load_data_from_file<int32_t, GDT_Int32>(
          snapshot_path / "identifyflats_sills.tif", sills, dims_check);
      assert(dims[0] == dims_check[0] && dims[1] == dims_check[1]);
    }

    // Allocate and resize output and intermediate arrays
    if (dem.size() > 0) {
      if (filled_dem.size() > 0) {
        // Boundary conditions
        bc.resize(dims[0] * dims[1]);
        test_filled_dem.resize(dims[0] * dims[1]);
      }
      if (flats.size() > 0 && sills.size() > 0) {
        test_flats.resize(dims[0] * dims[1]);
      }
    }
  }

  void runtests() {
    // fillsinks
    if (test_filled_dem.size() > 0) {
      // Initialize bcs
      for (ptrdiff_t j = 0; j < dims[1]; j++) {
        for (ptrdiff_t i = 0; i < dims[0]; i++) {
          if (isnan(dem[j * dims[0] + i])) {
            bc[j * dims[0] + i] = 1;
            test_dem[j * dims[0] + i] = -INFINITY;
          } else {
            if (i == 0 || i == dims[0] - 1 || j == 0 || j == dims[1] - 1) {
              // 1 on the boundaries
              bc[j * dims[0] + i] = 1;
            } else {
              // 0 on the interior
              bc[j * dims[0] + i] = 0;
            }
          }
        }
      }

      tt::fillsinks(test_filled_dem.data(), test_dem.data(), bc.data(),
                    dims.data());

      for (ptrdiff_t j = 0; j < dims[1]; j++) {
        for (ptrdiff_t i = 0; i < dims[0]; i++) {
          if (!isnan(filled_dem[j * dims[0] + i])) {
            assert(test_filled_dem[j * dims[0] + i] ==
                   filled_dem[j * dims[0] + i]);
          }
        }
      }
    }
    if (flats.size() > 0 && sills.size() > 0) {
      // identifyflats
      //
      // Use the snapshot filled DEM rather than the generated one in
      // case fillsinks fails.
      tt::identifyflats(test_flats.data(), filled_dem.data(), dims.data());

      for (ptrdiff_t j = 0; j < dims[1]; j++) {
        for (ptrdiff_t i = 0; i < dims[0]; i++) {
          if (flats[j * dims[0] + i] == 1)
            assert(test_flats[j * dims[0] + i] & 1);

          if (sills[j * dims[0] + i] == 1)
            assert(test_flats[j * dims[0] + i] & 2);
        }
      }
    }
  }
};

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cout << "Usage: snapshot <snapshot_directory>" << std::endl;
    return 0;
  }

  // Get the snapshot directory path from the command line arguments
  std::filesystem::path snapshot_dirpath(argv[1]);

  try {
    std::filesystem::directory_iterator dir(snapshot_dirpath);
  } catch (const std::filesystem::filesystem_error& e) {
    // Do not fail the test if the snapshots/data directory is not
    // present
    std::cout << "Snapshot directory does not exist" << std::endl;
    return 0;
  }

  for (const auto& entry :
       std::filesystem::directory_iterator(snapshot_dirpath)) {
    std::cout << entry.path().filename() << std::endl;
    SnapshotData data(entry.path());

    data.runtests();
  }
}