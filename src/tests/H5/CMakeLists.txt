
add_executable(testH5Cpp testH5Cpp.cpp)
find_package(Boost REQUIRED COMPONENTS "random")
target_link_libraries(testH5Cpp "common;${myLibs};${Boost_RANDOM_LIBRARY}")
install(TARGETS testH5Cpp RUNTIME DESTINATION bin)

add_executable(read testH5Cpp_read.cpp)
target_link_libraries(read "common;${myLibs}")
install(TARGETS read RUNTIME DESTINATION bin)

