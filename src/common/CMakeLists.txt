cmake_minimum_required(VERSION 2.6)

#include_directories(../adaption)
file(GLOB CommonFiles "*.cpp")
add_library(common STATIC ${CommonFiles})
# Das wuerde ich deswegen in das Haupt-CMakeLists skript verschieben.
target_link_libraries(common ${HDF5_LIBRARIES})
target_link_libraries(common mwlib)
if (MPI_NEEDED)
  target_link_libraries(common ${MPI_CXX_LIBRARIES})
endif()


file(GLOB CommonFiles2 "vesselgen/*.cpp")
add_library(vesselgenlib STATIC ${CommonFiles2})

install(DIRECTORY . DESTINATION "include/common" FILES_MATCHING PATTERN "*.h")
