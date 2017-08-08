cmake_minimum_required(VERSION 3.3)

###############  set up spack stuff  ##############

# enable @rpath in the install name for any shared library being built
# note: it is planned that a future version of CMake will enable this by default
set(CMAKE_MACOSX_RPATH 1)

# Always use full RPATH
# http://www.cmake.org/Wiki/CMake_RPATH_handling
# http://www.kitware.com/blog/home/post/510

# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# the RPATH to be used when installing, but only if it's not a system directory
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
   #this one is needed for the project itself
   #SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/milotti_mts/lib")
ENDIF("${isSystemDir}" STREQUAL "-1")

# Include all the transitive dependencies determined by Spack.
# If we're not running with Spack, this does nothing...
include_directories($ENV{SPACK_TRANSITIVE_INCLUDE_PATH})

###############  end spack stuff  ##############

### build external stuff numpycpp is still needed, maybe we could get rid of this in future withe build in stuff
### adapted from:
#   https://crascit.com/2015/07/25/cmake-gtest/
# Download and unpack numpycpp at configure time
configure_file(externalProjects/numpycpp/CMakeLists.txt.in externalProjects/numpycpp/CMakeLists.txt)
message("CMAKE_COMMAND:   ${CMAKE_COMMAND}")
message("CMAKE_GENERATOR:   ${CMAKE_GENERATOR}")
# numpycpp has its own modules
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${CMAKE_BINARY_DIR}/externalProjects/numpycpp-src/.cmake/modules;")
message("CMAKE_MODULE_PATH at numpycpp build:   ${CMAKE_MODULE_PATH}")
execute_process(COMMAND "${CMAKE_COMMAND}" -G "${CMAKE_GENERATOR}" .
  WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/externalProjects/numpycpp" )
execute_process(COMMAND "${CMAKE_COMMAND}" --build .
  WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/externalProjects/numpycpp" )
include_directories("${CMAKE_BINARY_DIR}/externalProjects/numpycpp-src")
# Add numpycpp directly to our build. This adds
# the following targets: numpycpp
add_subdirectory("${CMAKE_BINARY_DIR}/externalProjects/numpycpp-src"
                 "${CMAKE_BINARY_DIR}/externalProjects/numpycpp-build")

## HDF5-cpp-wrapper does not need compilation
configure_file(externalProjects/HDF5-cpp-wrapper/CMakeLists.txt.in externalProjects/HDF5-cpp-wrapper/CMakeLists.txt)
execute_process(COMMAND "${CMAKE_COMMAND}" -G "${CMAKE_GENERATOR}" .
  WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/externalProjects/HDF5-cpp-wrapper" )
execute_process(COMMAND "${CMAKE_COMMAND}" --build .
  WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/externalProjects/HDF5-cpp-wrapper" )
include_directories("${CMAKE_BINARY_DIR}/externalProjects/HDF5-cpp-wrapper-src")
# add_subdirectory("${CMAKE_BINARY_DIR}/externalProjects/HDF5-cpp-wrapper-src"
#                  "${CMAKE_BINARY_DIR}/externalProjects/HDF5-cpp-wrapper-build")
                 
project(allmystuff C CXX)
list(APPEND CMAKE_CXX_SOURCE_FILE_EXTENSIONS c)
add_definitions(-DSILENT)
find_package(Threads REQUIRED)
# set a search path for additional custom cmake modules.
# according to recommendation in https://cmake.org/Wiki/CMake:How_To_Find_Libraries
message("CMAKE_ROOT: ${CMAKE_ROOT}")
message("CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}")
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules/;")
message("CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}")
message("CMAKE_EXE_LINKER_FLAGS: ${CMAKE_EXE_LINKER_FLAGS}")
message("TrilinosDIR: ${TrilinosDIR}")


include(${CMAKE_ROOT}/Modules/CMakeDependentOption.cmake)

option(USE_MILOTTI_MTS "decide to use libmilotti_mts.so" ON)
if(USE_MILOTTI_MTS)
  message("MILOTTI_MTS preprocessor defined")
  add_definitions(-DMILOTTI_MTS)
endif()
option(USE_ADAPTION "decide to build the adaption library" OFF)
CMAKE_DEPENDENT_OPTION(USE_PAGMO "Use Pagmo" OFF "USE_ADAPTION" ON)
CMAKE_DEPENDENT_OPTION(USE_PAGMO_MPI "MPI version of Pagmo" OFF
                       "USE_ADAPTION;USE_PAGMO" ON)
option(TOTAL_SILENCE "suppress all output, good for running on cluster" OFF)
if(TOTAL_SILENCE)
  message("TOTAL_SILENCE preprocessor defined")
  add_definitions(-DTOTAL_SILENCE)
endif()


#It would be nicer to see if pagmo was compilled 
#with mpi automatically, pagmo2 is supporting that
# option(USE_PAGMO_MPI OFF)
# if(USE_PAGMO_MPI)
#   set(USE_PAGMO ON)
# endif()
if(USE_PAGMO)
#   if(NOT USE_ADAPTION)
#     message("so far there is no usage for pagmo without adation, we set it to be true!!!")
#     set(USE_ADAPTION ON)
#   endif()
  message("USE_PAGMO preprocessor defined")
  add_definitions(-DUSE_PAGMO)
endif()

if(USE_ADAPTION)
  message(STATUS "ADAPTION: USE_ADAPTION preprocessor defined")
  add_definitions(-DUSE_ADAPTION)
endif()

set (MAX_OPTIMIZATION TRUE CACHE BOOL "compile with -march=native (important for VM) and some other optimization options")
# set (MPI_NEEDED FALSE CACHE BOOL "Must we link to the MPI libraries?")
set (ADDITIONAL_INCLUDE_DIRS CACHE PATH "Additional dirs to search for included files")
set (ADDITIONAL_LIBRARY_DIRS CACHE PATH "Additional dirs to search for linked libraries")

find_package(Boost REQUIRED COMPONENTS "program_options")
message("Boost_PROGRAM_OPTIONS_LIBRARY: ${Boost_PROGRAM_OPTIONS_LIBRARY}")
set(myLibs "${myLibs};${Boost_PROGRAM_OPTIONS_LIBRARY}")
find_package(Boost REQUIRED COMPONENTS "system")
message("Boost_SYSTEM_LIBRARY: ${Boost_SYSTEM_LIBRARY}")
set(myLibs "${myLibs};${Boost_SYSTEM_LIBRARY}")
find_package(Boost REQUIRED COMPONENTS "serialization")
set(myLibs "${myLibs};${Boost_SERIALIZATION_LIBRARY}")
find_package(Boost REQUIRED COMPONENTS "thread")
set(myLibs "${myLibs};${Boost_THREAD_LIBRARY}")
find_package(Boost REQUIRED COMPONENTS "filesystem")
set(myLibs "${myLibs};${Boost_FILESYSTEM_LIBRARY}")


#note there are system without eigen, so let the configure fail here
find_package(Eigen3 3.2 QUIET)
#find_package(HDF5 1.8.10 REQUIRED)
find_package(HDF5 REQUIRED)
if(HDF5_IS_PARALLEL)
  message(WARNING "your hdf5 is parallel, you probably need mpi")
endif()
find_package(LAPACK REQUIRED)
find_package(BLAS REQUIRED)

find_package(PNG REQUIRED)
find_package(TBB REQUIRED)
find_package(Trilinos REQUIRED)
message(STATUS "Your Trilions depends on: ${Trilinos_TPL_LIST}")
if (";${Trilinos_TPL_LIST};" MATCHES ";MPI;")
  message(STATUS "MPI is NEEDED")
  set(MPI_NEEDED TRUE)
endif()

if (MPI_NEEDED OR USE_PAGMO_MPI)
  #set(MPI_CXX_COMPILER "/usr/bin/mpic++.openmpi") --> system dependent, move to do configure
  find_package(MPI REQUIRED)
  add_definitions(-DMPI_NEEDED)
  message(STATUS "Including MPI")
  message(STATUS "MPI libs         = ${MPI_CXX_LIBRARIES}")
  message(STATUS "MPI include dirs = ${MPI_CXX_INCLUDE_PATH}")
  message(STATUS "MPI_LIBRARY = ${MPI_LIBRARY}")
  message(STATUS "MPI_CXX_COMPILER_FLAGS = ${MPI_CXX_COMPILER_FLAGS}")
  message(STATUS "MPI_CXX_LINK_FLAGS = ${MPI_CXX_LINK_FLAGS}")
  
  include_directories(${MPI_CXX_INCLUDE_PATH})
  #set(myLibs ${MPI_CXX_LIBRARIES})
  set(myLibs "${myLibs};${MPI_CXX_LIBRARY}")
  if (HDF5_IS_PARALLEL) #good here
  message(WARNING "good choice here")
  else()
  message(WARNING "You are using trilinos with mpi support and hdf5 without, do you know what you are doing?")
  endif()
endif()

#from cmake 2.8 module documentation:
#The Python_ADDITIONAL_VERSIONS variable can be used to specify a list of version numbers that should be # taken into account when searching for Python. You need to set this variable before calling find_package(PythonLibs).
#for manual configuration, set PYTHON_INCLUDE_DIR and PYTHON_LIBRARY
#see also do-configure-tumorcode-lusi.sh script
set(Python_ADDITIONAL_VERSIONS 2.7)
find_package(PythonInterp REQUIRED)
find_package(PythonLibs REQUIRED)
find_package(Numpy REQUIRED)
include_directories(${PYTHON_INCLUDE_DIR})
include_directories(${PYTHON_NUMPY_INCLUDE_DIR})#this is needed here so that the whole system uses the same numpy source, if this is let out here the compiler may fall back to default numpy
#this package has a different name on ubuntu14
#so either one of this is needed!

find_package(Boost REQUIRED COMPONENTS "python")
set(myLibs "${myLibs};${Boost_PYTHON_LIBRARY}")
#I have to deal with this on the new arch system
# if( NOT Boost_PYTHON-2.7_FOUND)
#   #we are on ubuntu14
#   find_package(Boost REQUIRED COMPONENTS "python-py27")
# endif()
# find_package(Boost COMPONENTS "python-py27")
# if( NOT Boost_PYTHON-PY27_FOUND)
#   #we are not on ubuntu
#   find_package(Boost REQUIRED COMPONENTS "python-2.7")
# endif()



message(STATUS "PYTHON_INCLUDE_DIRS 	=${PYTHON_INCLUDE_DIRS}")
message(STATUS "PYTHON_LIBRARIES         	=${PYTHON_LIBRARIES}")
IF(Boost_PYTHON-PY27_FOUND)
  message(STATUS "Boost_PYTHON-PY27_LIBRARY 	=${Boost_PYTHON-PY27_LIBRARY}")
  message(STATUS "Boost_PYTHON-PY27_LIBRARY_DEBUG = ${Boost_PYTHON-PY27_LIBRARY_DEBUG}")
  message(STATUS "Boost_PYTHON-PY27_LIBRARY_RELEASE = ${Boost_PYTHON-PY27_LIBRARY_RELEASE}")
  set(myLibs "${myLibs};${Boost_PYTHON-PY27_LIBRARY}")
ENDIF()
IF(Boost_PYTHON-2.7_FOUND)
  message(STATUS "Boost_PYTHON-2.7_LIBRARY 	=${Boost_PYTHON-2.7_LIBRARY}")
  message(STATUS "Boost_PYTHON-2.7_LIBRARY_DEBUG = ${Boost_PYTHON-2.7_LIBRARY_DEBUG}")
  message(STATUS "Boost_PYTHON-2.7_LIBRARY_RELEASE = ${Boost_PYTHON-2.7_LIBRARY_RELEASE}")
  set(myLibs "${myLibs};${Boost_PYTHON-2.7_LIBRARY}")
ENDIF()
message(STATUS "Numpy include dirs  	=${PYTHON_NUMPY_INCLUDE_DIR}")
message(STATUS "PYTHON_LIBRARIES 	=${PYTHON_LIBRARIES}")

message(STATUS "Boost_PROGRAM_OPTIONS_LIBRARY_RELEASE         	=${Boost_PROGRAM_OPTIONS_LIBRARY_RELEASE}")
message(STATUS "LAPACK lib 					=${LAPACK_LIBRARIES}")
message(STATUS "BLAS_LIBRARIES   	= ${BLAS_LIBRARIES}")
message(STATUS "HDF5_CXX_LIBRARIES      = ${HDF5_CXX_LIBRARIES}")
message(STATUS "HDF5_C_LIBRARIES    	= ${HDF5_C_LIBRARIES}")
message(STATUS "HDF5_LIBRARY_DIRS	= ${HDF5_LIBRARY_DIRS}")
message(STATUS "HDF5_INCLUDE_DIRS	= ${HDF5_INCLUDE_DIRS}")


message(STATUS "Boost_lib dirs     	= ${Boost_LIBRARY_DIRS}")
message(STATUS "Boost include dirs 	= ${Boost_INCLUDE_DIRS}")
message(STATUS "Boost_LIBRARIES 	= ${Boost_LIBRARIES}")

#TF: somehow the cmake is broken here, I overwrite with environmental variabels
message(STATUS "TBB_LIBRARIES          	= $ENV{TBB_LIBRARY_RELEASE}")
set(TBB_LIBRARIES_DIR $ENV{TBB_LIBRARY_RELEASE})
message(STATUS "TBB_INCLUDE_DIR 	      = $ENV{TBB_INCLUDE_DIR}")
set(TBB_INCLUDE_DIRS $ENV{TBB_INCLUDE_DIR})
message(STATUS "Trilinos_dir          	= ${Trilinos_DIR}")
message(STATUS "Trilinos libs         	= ${Trilinos_LIBRARIES}")
message(STATUS "Trilinos lib dirs     	= ${Trilinos_LIBRARY_DIRS}")
message(STATUS "Trilinos include dirs 	= ${Trilinos_INCLUDE_DIRS}")


include_directories(${PYTHON_INCLUDE_DIRS})
include_directories(${Trilinos_INCLUDE_DIRS})
link_directories(${Trilinos_LIBRARY_DIRS})
include_directories(${HDF5_INCLUDE_DIRS})
link_directories(${HDF5_LIBRARY_DIRS})
include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})
include_directories(${PNG_INCLUDE_DIRS})
include_directories(${TBB_INCLUDE_DIRS})
link_directories(${TBB_LIBRARIES_DIR})
include_directories(${ADDITIONAL_INCLUDE_DIRS})
link_directories(${ADDITIONAL_LIBRARY_DIRS})
include_directories("src/")

set(myLibs
  ${myLibs} #maybe something is already there
  common  # part of this software
  mwlib   # part of this software
  numpycpp #see https://github.com/DaWelter/numpycpp
  tbb
  #${TBB_LIBRARIES}
  ${Trilinos_LIBRARIES} # -> trilinos
  ${Boost_PROGRAM_OPTIONS_LIBRARY_RELEASE}
  ${Boost_PYTHON_LIBRARY}
  ${HDF5_LIBRARIES}
  ${LAPACK_LIBRARIES}
  ${BLAS_LIBRARIES}
  ${PNG_LIBRARIES}
  ${PYTHON_LIBRARIES}
)
#message(STATUS "myLibs = ${myLibs}")


### compiler stuff
# this must come first as it overwrites preceeding compiler flags
# need openmp
# all variables propagate down the tree so set up as much as possible here
option(USE_OPENMP CACHE ON)
if (USE_OPENMP)
  message(STATUS "openMP: USE_OPENMP preprocessor defined")
  find_package(OpenMP REQUIRED)
  #add_definitions(-D_OPENMP) #is this done by find_package
  if (NOT OPENMP_FOUND)
    message(FATAL_ERROR "OpenMP not found")
  endif()
  message(STATUS "OpenMP_CXX_FLAGS      = ${OpenMP_CXX_FLAGS}")
  message(STATUS "OpenMP_C_FLAGS      = ${OpenMP_C_FLAGS}")
  set (CMAKE_C_FLAGS  "${OpenMP_C_FLAGS}   -DOPENMP_ENABLED=1 ${CMAKE_C_FLAGS}  ")
  set(CMAKE_CXX_FLAGS "${OpenMP_CXX_FLAGS} -DOPENMP_ENABLED=1 ${CMAKE_CXX_FLAGS}")
  #to get rid of the no #pragma warnings
  #add_definitions(-DENABLE_OMP)
endif()


GET_FILENAME_COMPONENT(compiler_name ${CMAKE_C_COMPILER} NAME)
string(COMPARE EQUAL ${compiler_name} "icc" IS_INTEL_COMPILER)
if (NOT IS_INTEL_COMPILER)
  message("Set for gnu gxx")
  set(Warnings "-Wall -Wextra -Wno-sign-compare -Wno-reorder -Wno-deprecated -Wuninitialized -Wno-unused-parameter -Wno-unused-variable -Wno-unused-result -Wno-unused-function -Wno-unused-local-typedefs") # -Wno-unknown-pragmas   -Wpointer-arith -Wwrite-strings -Woverloaded-virtual -Wno-long-long  -Wuninitialized ")
  set(Warnings "${Warnings} -Wno-deprecated-declarations")  # -Wswitch -Wno-extra -Wsynth
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Warnings} -pthread -fPIC -ftemplate-depth-128")
else()
  #1478 is a deprecation warning for auto_ptr
  #set(Warnings "-gcc -cxxlib -wd1478") #gnu++0x
  #set(Warnings "-wd1478") #gnu++0x
  message("Set for intel using:")
  message("CMAKE_C_COMPILER: 	${CMAKE_C_COMPILER}")
  message("CMAKE_CXX_COMPILER: 	${CMAKE_CXX_COMPILER}")
endif()

#always use new c++11 standard
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")

if (IS_INTEL_COMPILER)
    #1478 is a deprecation warning for auto_ptr
    set(Warnings "${Warnings} -wd1478")
endif()

#-fstrict-aliasing
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Warnings} -pthread -fPIC -ftemplate-depth-128")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG -g -Wno-unknown-pragmas" )
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DNDEBUG -O2 -g") # -funroll-all-loops
if (MAX_OPTIMIZATION)
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -funroll-loops -march=native")
endif()


set(PYTHON_LIB_COMMON_INSTALL_DIR "lib")

#set(CMAKE_SKIP_RPATH TRUE)

set(LIBRARY_OUTPUT_PATH ${PYTHON_LIB_COMMON_INSTALL_DIR})

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")


# this makes that all executables are placed in the top level build dir
#if (CMAKE_BUILD_TYPE STREQUAL "Release")

#set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
#message("executables go to ${EXECUTABLE_OUTPUT_PATH}")

#endif()

message(STATUS "myLibs = ${myLibs}")

# defer to build scripts in the src dir
add_subdirectory(src)
add_subdirectory(utils)
