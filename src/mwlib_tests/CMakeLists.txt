 cmake_minimum_required(VERSION 2.6)
#PROJECT(mwlib_tests)

include_directories(../common)

add_executable(test_array test_array.cpp)
# add_executable(test_hdf test_hdf.cpp)
# target_link_libraries(test_hdf mwlib ${testLibs} ${myLibs})
# install(TARGETS test_hdf RUNTIME DESTINATION bin)
add_executable(test_lookup_grid test_lookup_grid.cpp)

set(testLibs mwlib hdf5 png tbb)
target_link_libraries(test_array ${testLibs})

target_link_libraries(test_lookup_grid mwlib ${testLibs})

#set_target_properties(test_array test_hdf test_lookup_grid PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)

install(TARGETS test_array test_lookup_grid RUNTIME DESTINATION bin)