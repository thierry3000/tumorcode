
add_executable(test_belos test_belos.cpp)
target_link_libraries(test_belos ${myLibs})
install(TARGETS test_belos RUNTIME DESTINATION bin)

add_executable(PrecBiCGStabExFile PrecBiCGStabExFile.cpp)
target_link_libraries(PrecBiCGStabExFile ${myLibs})
