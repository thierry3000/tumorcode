cmake_minimum_required(VERSION 2.8)
project(detailedo2)

include_directories(../common)

file(GLOB DetailedO2Files "*.cpp")

add_library(detailedo2_ SHARED
            ${DetailedO2Files})

target_link_libraries(detailedo2_ krebs_)

#set_target_properties(detailedo2_ PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)

install(TARGETS detailedo2_
	LIBRARY DESTINATION ${PYTHON_LIB_COMMON_INSTALL_DIR})
