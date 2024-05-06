add_executable(versioninfo versioninfo.cpp)
target_link_libraries(versioninfo PRIVATE topotoolbox)
add_test(NAME versioninfo COMMAND versioninfo)
set_tests_properties(versioninfo PROPERTIES ENVIRONMENT_MODIFICATION
"PATH=path_list_prepend:$<$<BOOL:${WIN32}>:$<TARGET_FILE_DIR:topotoolbox>>")

add_executable(random_dem random_dem.cpp)
target_link_libraries(random_dem PRIVATE topotoolbox)
add_test(NAME random_dem COMMAND random_dem)
set_tests_properties(random_dem PROPERTIES ENVIRONMENT_MODIFICATION
  "PATH=path_list_prepend:$<$<BOOL:${WIN32}>:$<TARGET_FILE_DIR:topotoolbox>>")

add_executable(benchmark benchmark.cpp)
target_link_libraries(benchmark PRIVATE topotoolbox)
add_test(NAME benchmark COMMAND benchmark)
set_tests_properties(benchmark PROPERTIES ENVIRONMENT_MODIFICATION
"PATH=path_list_prepend:$<$<BOOL:${WIN32}>:$<TARGET_FILE_DIR:topotoolbox>>")
