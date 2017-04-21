cmake_minimum_required(VERSION 2.6)
project(adaption)

file(GLOB AdaptionFiles "*.cpp")

add_library(adaption_ SHARED
            ${AdaptionFiles})
#but not the headers
#by this flag we can use decide if pagmo uses MPI

if(USE_PAGMO)
  message("pagmo should be mpi")
  target_link_libraries(adaption_ krebs_ pagmoMPI ${MPI_CXX_LIBRARIES})
else()
  target_link_libraries(adaption_ krebs_)
endif()
#set_target_properties(adaption_ PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)

install(TARGETS adaption_
	LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
