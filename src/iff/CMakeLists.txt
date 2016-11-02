cmake_minimum_required(VERSION 2.8.8)
project(iff)

include_directories(../common)

set(CMAKE_DEBUG_POSTFIX "d")

file(GLOB Iff_Files "*.cpp")

add_library(iff_ SHARED
            ${Iff_Files})

#set_target_properties(iff_ PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)

target_link_libraries(iff_ krebs_)

install(TARGETS iff_
	LIBRARY DESTINATION ${PYTHON_LIB_COMMON_INSTALL_DIR})
