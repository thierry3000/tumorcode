cmake_minimum_required(VERSION 2.8.8)

#include_directories(../)
include_directories(../common)

if(USE_MILOTTI_MTS)
  file(GLOB Tumor_Files "*.cpp")
  set(myLibs "${myLibs};vbl")
  message("CMAKE_CXX_COMPILER:        ${CMAKE_CXX_COMPILER}")
  message("CMAKE_CXX_COMPILER_ID:     ${CMAKE_CXX_COMPILER_ID}")
  if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    message("using GNU build system")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -frounding-math")
  endif()
  ### Boost
  find_package(Boost REQUIRED COMPONENTS "thread")
  set(myLibs "${myLibs};${Boost_THREAD_LIBRARY}")

  ### OPENMP
  find_package(OpenMP REQUIRED)
  if (OPENMP_FOUND)
    message(STATUS "OpenMP_CXX_FLAGS      = ${OpenMP_CXX_FLAGS}")
    message(STATUS "OpenMP_C_FLAGS      = ${OpenMP_C_FLAGS}")
    set (CMAKE_C_FLAGS  "${OpenMP_C_FLAGS}   -DOPENMP_ENABLED=1 ${CMAKE_C_FLAGS}  ")
    set(CMAKE_CXX_FLAGS "${OpenMP_CXX_FLAGS} -DOPENMP_ENABLED=1 ${CMAKE_CXX_FLAGS}")
  endif()
  ### CGAL
  find_package(CGAL REQUIRED)
  if(CGAL_FOUND)
    message("Found CGAL_LIBRARY: ${CGAL_LIBRARY})")
    set(myLibs "${myLibs};${CGAL_LIBRARY}" )
  endif()

  ### GMP
  set(myLibs "${myLibs};gmp")
  
  ### ANN
  set(myLibs "${myLibs};ann")

  #message("set comiler flags: ${CMAKE_CXX_FLAGS}")
else()
  file(GLOB Tumor_Files bulktissue-no-vessels.cpp bulktissue-with-vessels.cpp faketum.cpp python-tumors.cpp)
endif()

add_library(tumors_ SHARED ${Tumor_Files})

set(myLibs "${myLibs};krebs_")
if(USE_ADAPTION)
  set(myLibs "${myLibs};adaption_")
endif()
# if(USE_MILOTTI_MTS)
#   set(myLibs "${myLibs};detailedo2_")
# endif()

# if(USE_ADAPTION)
#   set(link_tumor_to_this 
#   krebs_
#   adaption_
#   ${myLibs}
#   )
# else()
#   set(link_tumor_to_this 
#   krebs_
#   ${myLibs}
#   )
# endif()
message("link_tumor_to_this: ${myLibs} ")
target_link_libraries(tumors_ ${myLibs})
target_link_libraries(tumors_ krebs_)
target_link_libraries(tumors_ detailedo2_)
if(USE_MILOTTI_MTS)
  target_link_libraries(tumors_ vbl)
  #target_link_libraries(tumors_ gmp)
  #target_link_libraries(tumors_ CGAL)
  #target_link_libraries(tumors_ ANN)
endif()

set_target_properties(tumors_ PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)

install(TARGETS tumors_ 
	LIBRARY DESTINATION ${PYTHON_LIB_COMMON_INSTALL_DIR})
