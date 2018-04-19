cmake_minimum_required(VERSION 2.6)
PROJECT(mwlib)

#include_directories('.')

file(GLOB mwLibFiles "*.cpp")

add_library(mwlib STATIC ${mwLibFiles})
target_link_libraries(mwlib common)

install(DIRECTORY . DESTINATION "include/mwlib" FILES_MATCHING PATTERN "*.h")
