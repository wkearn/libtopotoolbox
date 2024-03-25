cmake_minimum_required(VERSION 3.5)

project(
  topotoolbox
  VERSION 3.0.0
  LANGUAGES C CXX
)

enable_testing()

add_subdirectory(src)
add_subdirectory(test)
