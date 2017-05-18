cmake_minimum_required(VERSION 2.6)

#PROJECT(openmptest)

set (TRILINOS_LIB "" CACHE PATH "Where to find trilinos libs")
set (TRILINOS_INCLUDE "" CACHE PATH "Where to find trilinos includes")

include_directories(
/home/usersHR/mwelter/tumorcode/src/mwlib
)

link_directories(
)

link_directories(
  ${TRILINOS_LIB}
)

include_directories(
  ${TRILINOS_INCLUDE}
)

#set(CMAKE_CXX_FLAGS "-fopenmp")

ADD_EXECUTABLE(openmptest main.cpp)

target_link_libraries(openmptest epetra teuchos ml galeri blas lapack)


