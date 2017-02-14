cmake_minimum_required(VERSION 2.8.8)
project(python_krebsutils)

include_directories(../common)

add_library(krebs_ SHARED
            mainboost.cpp pylatticedata.cpp python-helpers.cpp
            povray_export.cpp samplevessels.cpp python_tumormodel_tests.cpp
            python_numericaltools_tests.cpp python_vesselgen.cpp
            python_elliptic_solver_test.cpp python_calcflow.cpp python_interpolate_between_vessels.cpp
            python_murray.cpp)
if (False)
  find_package(MPI REQUIRED)
  message(STATUS "Including MPI for libkrebs")
  message(STATUS "MPI libs         = ${MPI_CXX_LIBRARIES}")
  message(STATUS "MPI include dirs = ${MPI_CXX_INCLUDE_PATH}")
  message(STATUS "MPI_LIBRARY = ${MPI_LIBRARY}")
  message(STATUS "MPI_CXX_COMPILER_FLAGS = ${MPI_CXX_COMPILER_FLAGS}")
  message(STATUS "MPI_CXX_LINK_FLAGS = ${MPI_CXX_LINK_FLAGS}")
  target_link_libraries(krebs_ common vesselgenlib numpycpp ${myLibs} ${Boost_LIBRARIES} ${PYTHON_LIBRARIES} ${MPI_CXX_LIBRARIES})
else()
  target_link_libraries(krebs_ common vesselgenlib numpycpp ${myLibs} ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})
endif()
set_target_properties(krebs_ PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)

install(TARGETS krebs_
	LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)

install(DIRECTORY ${CMAKE_SOURCE_DIR}/py
	DESTINATION .
	FILES_MATCHING PATTERN "*.py"
	)
#experimental
#install(FILES_MATCHING PATTERN "submit*.py"
#        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ)

#changing permissions of sumit script, make them executable
install(DIRECTORY ${CMAKE_SOURCE_DIR}/py
	DESTINATION .
	FILES_MATCHING PATTERN "submit*.py"
	PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
	)
	
#note: crete symlinks after python_krebsutils is successfully installed
#symlinks to ease usage
install(DIRECTORY DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
#it turned out that post install stuff is not so easy in cmake
install(SCRIPT "${CMAKE_SOURCE_DIR}/cmake/PostInstall.cmake")