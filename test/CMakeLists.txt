add_executable(versioninfo versioninfo.cpp)

target_link_libraries(versioninfo PRIVATE topotoolbox)

add_test(NAME versioninfo COMMAND versioninfo)

install(FILES $<TARGET_RUNTIME_DLLS:versioninfo> TYPE BIN)
