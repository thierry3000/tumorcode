include_directories(../)
include_directories(../common)
file(GLOB Tumor_Files "*.cpp")
#file(GLOB Tumor_Files "bulktissue-with-vessels.cpp python_tumors.cpp")

#add_library(bulktissue_with STATIC bulktissue-with-vessels.cpp)

add_library(tumor_ SHARED ${Tumor_Files})
if(USE_MILOTTI_MTS)
  set(myLibs "${myLibs};milotti_mts")
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

  #message("set comiler flags: ${CMAKE_CXX_FLAGS}")
  
endif()
if(USE_ADAPTION)
  set(link_tumor_to_this 
  krebs_
  adaption_
  ${myLibs}
  )
else()
  set(link_tumor_to_this 
  krebs_
  ${myLibs}
  )
endif()
message("link_tumor_to_this: ${link_tumor_to_this} ")
target_link_libraries(tumor_ ${link_tumor_to_this})

set_target_properties(tumor_ PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)

install(TARGETS tumor_ 
	LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
