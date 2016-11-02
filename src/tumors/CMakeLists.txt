include_directories(../)
include_directories(../common)
file(GLOB Tumor_Files "*.cpp")

add_library(bulktissue_with STATIC bulktissue-with-vessels.cpp)

add_library(tumor_ SHARED ${Tumor_Files})
target_link_libraries(tumor_
#vesselgenlib
bulktissue_with
adaption_ 
krebs_
numpycpp ${myLibs} ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})

set_target_properties(tumor_ PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)

install(TARGETS tumor_ 
	LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)