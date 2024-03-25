# libtopotoolbox

A C++ library for the analysis of digital elevation models.

## Build instructions

From the top level directory of the repository, generate the project
buildsystem in the `build/` directory by running

```
> cmake -B build
```

and then build the project with

```
> cmake --build build
```

## Running the tests

Once the build has completed, run the tests by moving into the build
directory and running the CTest executable that comes with CMake.

```
> cd build
> ctest
```
