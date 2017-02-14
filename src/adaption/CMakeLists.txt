cmake_minimum_required(VERSION 2.6)
project(adaption)

#include_directories(../common)

file(GLOB AdaptionFiles "*.cpp")

add_library(adaption_ SHARED
            ${AdaptionFiles})
#target_link_libraries(adaption_ krebs_ pagmo ${MPI_CXX_LIBRARIES})-->already done by pagmo
target_link_libraries(adaption_ krebs_ pagmo )

#set_target_properties(adaption_ PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)

install(TARGETS adaption_
	LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
