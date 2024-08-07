OPTION(TT_SANITIZE "Build with AddressSanitizer" OFF)

add_executable(versioninfo versioninfo.cpp utils.c)
target_link_libraries(versioninfo PRIVATE topotoolbox)
add_test(NAME versioninfo COMMAND versioninfo)
set_tests_properties(versioninfo PROPERTIES ENVIRONMENT_MODIFICATION
"PATH=path_list_prepend:$<$<BOOL:${WIN32}>:$<TARGET_FILE_DIR:topotoolbox>>")

add_executable(random_dem random_dem.cpp utils.c)
if (TT_SANITIZE AND NOT MSVC)
  target_compile_options(random_dem PRIVATE "$<$<CONFIG:DEBUG>:-fsanitize=address>")
  target_link_options(random_dem PRIVATE "$<$<CONFIG:DEBUG>:-fsanitize=address>")
endif()
target_link_libraries(random_dem PRIVATE topotoolbox)
add_test(NAME random_dem COMMAND random_dem)
set_tests_properties(random_dem PROPERTIES ENVIRONMENT_MODIFICATION
  "PATH=path_list_prepend:$<$<BOOL:${WIN32}>:$<TARGET_FILE_DIR:topotoolbox>>")

add_executable(excesstopography excesstopography.cpp utils.c)
if (TT_SANITIZE AND NOT MSVC)
  target_compile_options(random_dem PRIVATE "$<$<CONFIG:DEBUG>:-fsanitize=address>")
  target_link_options(random_dem PRIVATE "$<$<CONFIG:DEBUG>:-fsanitize=address>")
endif()
target_link_libraries(excesstopography PRIVATE topotoolbox)
add_test(NAME excesstopography COMMAND excesstopography)
set_tests_properties(excesstopography PROPERTIES ENVIRONMENT_MODIFICATION
  "PATH=path_list_prepend:$<$<BOOL:${WIN32}>:$<TARGET_FILE_DIR:topotoolbox>>")
