cmake_minimum_required(VERSION 2.6)

# bash wrapper
message(STATUS "Installing bash wrapper")
#install_programs( "${CMAKE_INSTALL_PREFIX}/utils" bash_wrapper.sh)
if(USE_MILOTTI_MTS)
  add_executable(vessel_wrapper tumorcode_to_milotti_wrapper.cpp)
  target_link_libraries(vessel_wrapper boost_program_options krebs_ vbl gmp CGAL)
  install(TARGETS vessel_wrapper RUNTIME DESTINATION bin)
endif()
install_programs( "/utils" bash_wrapper.sh)

# 	DESTINATION ${CMAKE_INSTALL_PREFIX}/utils)
