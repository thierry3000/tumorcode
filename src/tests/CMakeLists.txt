cmake_minimum_required(VERSION 2.6)
#PROJECT(programs)

include_directories(../common)
include_directories(${ADDITIONAL_INCLUDE_DIRS})

add_executable(tests tests.cpp)
target_link_libraries(tests ${myLibs})

# add_executable(test_convection test_convection1.cpp)
# target_link_libraries(test_convection ${myLibs})

add_executable(test_stepper test_stepper.cpp)
target_link_libraries(test_stepper ${myLibs})

# add_executable(test_level_set test_level_set.cpp)
# target_link_libraries(test_level_set ${myLibs})

#add_executable(test_elliptic_solver test_elliptic_solver.cpp)
#target_link_libraries(test_elliptic_solver ${myLibs})

add_executable(test_pagmo test_pagmo.cpp)
#set_target_properties(test_pagmo PROPERTIES DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX} LINK_FLAGS "-rdynamic")
#but not the headers
#by this flag we can use decide if pagmo uses MPI

if(PAGMO_MPI)
  message("pagmo should be mpi")
  target_link_libraries(test_pagmo pagmoMPI ${myLibs})
else()
  target_link_libraries(test_pagmo krebs_ pagmo boost_system)
endif()

#target_link_libraries(test_pagmo pagmo boost_system ${myLibs} )

add_executable(test_pagmo_s test_pagmo_serialize.cpp)
#set_target_properties(test_pagmo PROPERTIES DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX} LINK_FLAGS "-rdynamic")
if(PAGMO_MPI)
  message("pagmo should be mpi")
  target_link_libraries(test_pagmo_s adaption_ common ${myLibs} pagmoMPI)
else()
  target_link_libraries(test_pagmo_s adaption_ common ${myLibs} pagmo)
endif()

#target_link_libraries(test_pagmo_s pagmo boost_system adaption_ common ${myLibs} )

add_executable(test_boost_serialize test_boost_serialize.cpp)
target_link_libraries(test_boost_serialize boost_serialization)



# add_executable(test_h5cpp test_hdf5_cpp.cpp)
# link_directories(${HDF5_LIBRARY_DIRS})
# target_link_libraries(test_h5cpp hdf5 hdf5_cpp)
# install(TARGETS test_h5cpp RUNTIME DESTINATION bin)

# add_executable(test_convection_diffusion test_convection_diffusion.cpp)
# target_link_libraries(test_convection_diffusion ${myLibs})


install(TARGETS tests test_pagmo test_stepper test_pagmo_s test_boost_serialize RUNTIME DESTINATION bin)
