cmake_minimum_required(VERSION 2.6)

PROJECT(chombotest)

include_directories(
/home/mwelter/tumorcode/mylib2/include
/home/mwelter/local/chombo-3d-static/include
)

link_directories(
/home/mwelter/local/chombo-3d-static/lib
)


set(Chombovars "-DCH_SPACEDIM=3 -DCH_LANG_CC -DCH_USE_HDF5")
set(CMAKE_CXX_FLAGS "${Chombovars}")

ADD_EXECUTABLE(chombotest main.cpp)

target_link_libraries(chombotest amrelliptic3d amrtools3d boxtools3d basetools3d hdf5)
