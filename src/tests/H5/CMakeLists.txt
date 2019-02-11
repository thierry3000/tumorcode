
# these examples are not working with older boost implementations

add_executable(testH5Cpp testH5Cpp.cpp)
find_package(Boost REQUIRED COMPONENTS "random")
target_link_libraries(testH5Cpp "common;${myLibs};${Boost_RANDOM_LIBRARY};${HDF5_CXX_LIBRARIES}")
#install(TARGETS testH5Cpp RUNTIME DESTINATION bin)

# 
add_executable(read testH5Cpp_read.cpp)
target_link_libraries(read "common;${myLibs}")
#install(TARGETS read RUNTIME DESTINATION bin)

add_executable(test_cluster_stackoverflow test_cluster_stackoverflow.cpp)
target_link_libraries(test_cluster_stackoverflow "common;${myLibs}")

add_executable(test_tumorcode_hdfio test_tumorcode_hdfio.cpp)
target_link_libraries(test_tumorcode_hdfio "common;${myLibs};${Boost_RANDOM_LIBRARY};${HDF5_CXX_LIBRARIES}")

