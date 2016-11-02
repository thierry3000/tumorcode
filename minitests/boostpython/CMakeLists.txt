cmake_minimum_required(VERSION 2.6)
PROJECT(programs)

find_package(PythonLibs)
if (NOT PYTHONLIBS_FOUND)
  message (FATAL_ERROR "Python Libs not found!")
else()
  include_directories(${PYTHON_INCLUDE_DIRS})
endif()

LINK_LIBRARIES(boost_python)
ADD_LIBRARY(hello_ext SHARED hello.cpp)

SET_TARGET_PROPERTIES(hello_ext PROPERTIES PREFIX "")