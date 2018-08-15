cmake_minimum_required(VERSION 2.8)
#project(detailedo2)

include_directories(../common)

file(GLOB DetailedO2Files oxygen_model2.cpp python-oxygen2.cpp)
#oxygen_model2.cpp;python-oxygen2.cpp
# set(DetailedO2Files "${DetailedO2Files};../common/trilinos_linsys_construction.cpp")

add_library(detailedo2_ SHARED
            ${DetailedO2Files})

target_link_libraries(detailedo2_ krebs_ common)

#set_target_properties(detailedo2_ PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)

install(TARGETS detailedo2_
	LIBRARY DESTINATION ${PYTHON_LIB_COMMON_INSTALL_DIR})
