cmake_minimum_required(VERSION 2.8.8)
# we need quite a bit of python information for compiling the source
# declare that here, so that is it is recursive valid

set(CMAKE_DEBUG_POSTFIX "d")

add_subdirectory(detailedO2)
add_subdirectory(common)
add_subdirectory(mwlib)
add_subdirectory(mwlib_tests)
add_subdirectory(python_krebsutils)
add_subdirectory(tests)
add_subdirectory(adaption)
add_subdirectory(iff)
add_subdirectory(tumors)