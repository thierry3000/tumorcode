cmake_minimum_required(VERSION 2.6)
PROJECT(mwlib)

#include_directories('.')

file(GLOB mwLibFiles any_tree.cpp;drawing.cpp;field.cpp;hdf_wrapper.cpp;hdf_wrapper_array3d.cpp;hdf_wrapper_ptree.cpp;helpers.cpp;histogram.cpp;lattice-data.cpp;mempool.cpp;ptree_ext.cpp)

add_library(mwlib STATIC ${mwLibFiles})
target_link_libraries(mwlib common)

install(DIRECTORY . DESTINATION "include/mwlib" FILES_MATCHING PATTERN "*.h")
