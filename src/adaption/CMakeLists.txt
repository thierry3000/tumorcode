cmake_minimum_required(VERSION 2.6)
#project(adaption)

# if(USE_PAGMO)
#   #set(AdaptionFiles "python-adaption.cpp adaption_model2.cpp adaption_as_pagmo_problem.cpp")
#   #file(GLOB AdaptionFiles "python-adaption.cpp adaption_model2.cpp adaption_as_pagmo_problem.cpp")
#   add_library(adaption_ SHARED python-adaption.cpp adaption_model2.cpp adaption_as_pagmo_problem.cpp)
#   add_executable(run_pagmo run_pagmo.cpp)
# #  target_link_libraries(run_pagmo adaption_ pagmoMPI ${MPI_CXX_LIBRARIES})
#   target_link_libraries(run_pagmo adaption_ pagmo)
#   
#   install(TARGETS run_pagmo RUNTIME DESTINATION bin)
# else()
#   #set(AdaptionFiles "python-adaption.cpp adaption_model2.cpp")
#   #file(GLOB AdaptionFiles "python-adaption.cpp adaption_model2.cpp")
#   add_library(adaption_ SHARED python-adaption.cpp adaption_model2.cpp)
# endif()
#file(GLOB ${AdaptionFiles})

#message("AdaptionFiles : ${AdaptionFiles}")
#add_library(adaption_ SHARED python-adaption.cpp adaption_model2.cpp)
#but not the headers
#by this flag we can use decide if pagmo uses MPI

if(USE_PAGMO_MPI)
  message("pagmo should be mpi")
  find_package(MPI)
  message("MPI_CXX_LIBRARIES: ${MPI_CXX_LIBRARIES} ")
  add_library(adaption_ SHARED python-adaption.cpp adaption_model2.cpp adaption_as_pagmo_problem.cpp)
  target_link_libraries(adaption_ krebs_ pagmo ${MPI_CXX_LIBRARIES})
  add_executable(run_pagmo run_pagmo.cpp)
  install(TARGETS run_pagmo RUNTIME DESTINATION bin)
else()
  add_library(adaption_ SHARED python-adaption.cpp adaption_model2.cpp adaption_as_pagmo_problem.cpp)
  target_link_libraries(adaption_ krebs_)
  add_executable(run_pagmo run_pagmo.cpp)
  target_link_libraries(run_pagmo adaption_ pagmo)
  install(TARGETS run_pagmo RUNTIME DESTINATION bin)
endif()
set_target_properties(adaption_ PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)

install(TARGETS adaption_
	LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
