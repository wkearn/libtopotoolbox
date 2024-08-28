# Architecture

This document describes the design of libtopotoolbox and is intended
to help maintainers and developers understand the code base.

## High-level overview

libtopotoolbox is a library of basic algorithms for analyzing digital
elevation models. Each algorithm is implemented as a function that
accepts input data as in-memory arrays and writes its output data to
in-memory arrays. Reading and writing digital elevation models from
files is the concern of the program that calls into
libtopotoolbox. This significantly simplifies the architecture of
libtopotoolbox and allows libtopotoolbox to plug easily into existing
geospatial data ecosystems.

## Code structure

### `include/topotoolbox.h`

The public interface for libtopotoolbox is specified in this header
file. Every function that is exported by the library is forward
declared in this file with the `TOPOTOOLBOX_API` qualifier, which
expands to platform-specific attributes that ensure that the function
names are available in the compiled library.

Documentation for each function is provided as a Doxygen comment in
this file.

### `src/`

This directory contains C files that implement the API functions
declared in `topotoolbox.h`. Multiple related API functions are often
implemented in a single file. These functions are named with a prefix
corresponding to the implementation file name followed by an
underscore. For example, `excesstopography_fmm2d`,
`excesstopography_fsm2d`, and `excesstopography_fmm3d` are found in
`excesstopography.c`. This convention implements a kind of namespacing
and each file can be considered a module, though C does not have
actual modules.

Data structures and utility functions shared across modules can be
defined in their own C file. In this case, they should be accompanied
by a private header file in the `src/` directory (not `include/`!)
that forward declares the interface to the data structure or utility
functions. Modules that use these shared facilities can `#include` the
private header file. For an example, see the priority queue declared in
`priority_queue.h`, implemented in `priority_queue.c` and used in
`excesstopography.c` and `gwdt.c`.

Collections of related files can be included in a subdirectory of
`src/`. For example, mathematical morphology functions are implemented
as shared utility modules in the `src/morphology/` subdirectory. Using
subdirectories is another way to emulate a hierarchical module system.

### `test/`

Automated tests are implemented in the `test/` subdirectory. Automated
tests are compiled to executables that are linked to libtopotoolbox
and exercise some aspect of the libtopotoolbox code. The
[CTest](https://cmake.org/cmake/help/book/mastering-cmake/chapter/Testing%20With%20CMake%20and%20CTest.html)
framework included with CMake is used to run the test executables, but
each test executable contains many individual tests of libtopotoolbox
functionality. For example `test/random_dem.cpp` tests several
libtopotoolbox functions using randomly generated DEMs. It is compiled
into a single test executable that is run by CTest, but it implements
many tests.

At the moment the tests are implemented in C++ to ensure that our C
library can be successfully linked to C++ code, but tests can also be
implemented in C. More on the testing strategy can be found below.

### `docs/`

The documentation build system and additional documentation can be
found in the `docs/` subdirectory. The library itself is documented in
`topotoolbox.h` as described above. We use a combination of Doxygen to
generate documentation from the inline function documentation and
Sphinx to render the generated documentation. More in-depth narrative
documentation is written in reStructuredText files within this
subdirectory.

## Error handling

libtopotoolbox functions always succeed, so there is no need to handle
errors coming from libtopotoolbox.

libtopotoolbox functions generally do not perform IO or dynamic memory
allocation, eliminating two potential sources of errors. The only way
a libtopotoolbox function can "fail" is if it is provided data that
does not meet its assumptions. For example, most libtopotoolbox
functions will not perform as expected if the supplied digital
elevation model contains NaNs. In the case that they are provided
invalid data, they should successfully run to completion, but their
outputs will be invalid. The calling program is responsible for
ensuring that data passed to libtopotoolbox meets the expectations of
the called function and data returned from libtopotoolbox meets its
own expectations.

## Testing

Because it can be hard to specify a priori the result of an algorithm
on a given digital elevation model, we generally prefer property-based
tests. These tests generate random digital elevation models, pass them
to libtopotoolbox functions, and then verify that the output data
satisfies properties expected of those algorithms for every randomly
generated DEM.

For example, the properties that we test for `fillsinks` include

1. that each pixel of the filled digital elevation model should be
   higher than or at the same elevation as the same pixel of the
   original DEM and
2. that no pixel in the filled digital elevation model is completely
   surrounded by pixels higher than it.
   
These properties do not fully constrain the `fillsinks`
implementation, but they do give us some confidence that our
implementation behaves as expected. Additional properties can also be
easily added as they are developed or edge cases encountered.

This property-based testing strategy has the advantages of being fully
automatic and easy to implement. 

Other testing methods would be welcomed as a contribution.

## Build system

CMake is used to compile libtopotoolbox. The `CMakeLists.txt` file in
the libtopotoolbox root sets up the project and includes additional
`CMakeLists.txt` throughout the project.

The CMake build configuration should require only minor modifications
when changes are made to libtopotoolbox. Most importantly, new C
source files that need to be compiled must be added to the
`add_library` declaration at the top `src/CMakeLists.txt`.

New test executables must be declared in `test/CMakeLists.txt`. These
require some extra configuration to ensure they get compiled and
linked properly on different platforms. Unless a new test executable
requires special build configurations, it should be sufficient to copy
the declaration of an existing test to create a new test executable.
