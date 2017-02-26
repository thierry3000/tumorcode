cmake_minimum_required(VERSION 2.6)
project(adaption)

#include_directories(../common)

file(GLOB AdaptionFiles "*.cpp")

add_library(adaption_ SHARED
            ${AdaptionFiles})
#but not the headers
#by this flag we can use decide if pagmo uses MPI
option(PAGMO_MPI "is pagmo mpi?" ON)
if(PAGMO_MPI)
  ADD_DEFINITIONS(-DPAGMO_ENABLE_MPI)
  find_package(MPI)
  include_directories(${MPI_CXX_INCLUDE_PATH})
  target_link_libraries(adaption_ krebs_ pagmo ${MPI_CXX_LIBRARIES})
else()
  target_link_libraries(adaption_ krebs_ pagmo)
endif()
#set_target_properties(adaption_ PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)

install(TARGETS adaption_
	LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
