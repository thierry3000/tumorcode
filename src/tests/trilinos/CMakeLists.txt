# note Ifpack2 not provided on durga
string(COMPARE EQUAL "$ENV{CLUSTER}" "lusi-arch" string_compare_result)
#message("$ENV{CLUSTER}")
#message("string_compare_result: ${string_compare_result}")
if(string_compare_result)
  add_executable(test_belos test_belos.cpp)
  target_link_libraries(test_belos ${myLibs} )
  install(TARGETS test_belos RUNTIME DESTINATION bin)
endif()
add_executable(PrecBiCGStabExFile PrecBiCGStabExFile.cpp)
target_link_libraries(PrecBiCGStabExFile ${myLibs})
install(TARGETS PrecBiCGStabExFile RUNTIME DESTINATION bin)

add_executable(IfpackBelos2 IfpackBelos2.cpp)
#target_link_libraries(IfpackBelos2 ${myLibs} galeri-epetra)
target_link_libraries(IfpackBelos2 ${myLibs})
install(TARGETS IfpackBelos2 RUNTIME DESTINATION bin)
