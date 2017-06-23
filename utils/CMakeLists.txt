cmake_minimum_required(VERSION 2.6)

# bash wrapper
message(STATUS "Installing bash wrapper")
#install_programs( "${CMAKE_INSTALL_PREFIX}/utils" bash_wrapper.sh)
install_programs( "/utils" bash_wrapper.sh)
# 	DESTINATION ${CMAKE_INSTALL_PREFIX}/utils)
