cmake_minimum_required(VERSION 2.6)

# bash wrapper
message(STATUS "Installing bash wrapper")
#install_programs( "${CMAKE_INSTALL_PREFIX}/utils" bash_wrapper.sh)
add_executable(vessel_wrapper tumorcode_to_milotti_wrapper.cpp)
target_link_libraries(vessel_wrapper boost_program_options krebs_ milotti_mts gmp CGAL)
install_programs( "/utils" bash_wrapper.sh)
install(TARGETS vessel_wrapper RUNTIME DESTINATION bin)
# 	DESTINATION ${CMAKE_INSTALL_PREFIX}/utils)
