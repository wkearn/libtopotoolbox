add_library(topotoolbox
  topotoolbox.c
  fillsinks.c
  identifyflats.c
  gwdt.c
  morphology/reconstruct.c
  excesstopography.c
)

target_include_directories(
  topotoolbox
  PUBLIC
  $<BUILD_INTERFACE:${topotoolbox_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}>
)

target_link_libraries(topotoolbox m)

if (BUILD_SHARED_LIBS)
else()
  target_compile_options(topotoolbox PUBLIC -DTOPOTOOLBOX_STATICLIB)
endif()

set_target_properties(topotoolbox PROPERTIES PUBLIC_HEADER ${topotoolbox_SOURCE_DIR}/include/topotoolbox.h)
