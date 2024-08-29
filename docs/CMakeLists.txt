# Build the documentation using Doxygen and Sphinx
#
# This CMake file is largely derived from Sy Brand's 2019 article
# "Clear, Functional C++ Documentation with Sphinx + Breathe + Doxygen + CMake"
# on the Microsoft C++ Team Blog:
# https://web.archive.org/web/20240608084430/https://devblogs.microsoft.com/cppblog/clear-functional-c-documentation-with-sphinx-breathe-doxygen-cmake/

# Set up Doxygen
find_package(Doxygen REQUIRED)

# Doxygen-specific variables
set(DOXYGEN_INPUT_DIR "${PROJECT_SOURCE_DIR}/include")
set(DOXYGEN_OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}/docs/doxygen")
set(DOXYGEN_INDEX_FILE "${DOXYGEN_OUTPUT_DIR}/html/index.html")
set(DOXYFILE_IN "${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in")
set(DOXYFILE_OUT "${CMAKE_CURRENT_BINARY_DIR}/Doxyfile")

# Doxyfile.in is a template that is filled out and turned into the
# Doxyfile by configure_file
configure_file(${DOXYFILE_IN} ${DOXYFILE_OUT} @ONLY)

file(MAKE_DIRECTORY ${DOXYGEN_OUTPUT_DIR})

# Run Doxygen
#
# This will be run automatically because the Sphinx configuration
# below depends on the ${DOXYGEN_INDEX_FILE}
add_custom_command(
  OUTPUT ${DOXYGEN_INDEX_FILE}
  DEPENDS topotoolbox
  COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYFILE_OUT}
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  MAIN_DEPENDENCY ${DOXYFILE_OUT} ${DOXYFILE_IN}
  COMMENT "Generating API docs with Doxygen"
)

# Sphinx configuration

find_package(Sphinx REQUIRED)

# Sphinx-specific variables
set(SPHINX_SOURCE ${CMAKE_CURRENT_SOURCE_DIR})
set(SPHINX_BUILD "${CMAKE_CURRENT_BINARY_DIR}/docs/sphinx")
set(SPHINX_INDEX_FILE "${SPHINX_BUILD}/index.html")

# Find all of the reStructuredText files in the docs/ directory and
# its subdirectories. This is used to define the dependencies so that
# CMake knows when Sphinx needs to be rerun. This is only for the
# build system: Sphinx finds documentation source files itself.
file(GLOB_RECURSE TOPOTOOLBOX_DOCS_RST ${CMAKE_CURRENT_SOURCE_DIR}/*.rst)

# Run Sphinx
# 
# The Sphinx configuration is in conf.py
add_custom_command(
  OUTPUT ${SPHINX_INDEX_FILE}
  COMMAND
  ${SPHINX_EXECUTABLE} -b html
  -Dbreathe_projects.TopoToolbox="${DOXYGEN_OUTPUT_DIR}/xml"
  ${SPHINX_SOURCE} ${SPHINX_BUILD}
  DEPENDS ${TOPOTOOLBOX_DOCS_RST} ${DOXYGEN_INDEX_FILE}
  MAIN_DEPENDENCY ${CMAKE_CURRENT_SOURCE_DIR}/conf.py
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Generating documentation with Sphinx"
)

# The docs target can be built explicitly with
# cmake --build build --target docs
add_custom_target(docs ALL DEPENDS ${SPHINX_INDEX_FILE})

