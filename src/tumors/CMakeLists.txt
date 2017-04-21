include_directories(../)
include_directories(../common)
file(GLOB Tumor_Files "*.cpp")
#file(GLOB Tumor_Files "bulktissue-with-vessels.cpp python_tumors.cpp")

#add_library(bulktissue_with STATIC bulktissue-with-vessels.cpp)

add_library(tumor_ SHARED ${Tumor_Files})

if(USE_ADAPTION)
  set(link_tumor_to_this 
  krebs_
  adaption_
  ${myLibs}
  )
else()
  set(link_tumor_to_this 
  krebs_
  ${myLibs}
  )
endif()
message("link_tumor_to_this: ${link_tumor_to_this} ")
target_link_libraries(tumor_ ${link_tumor_to_this})

set_target_properties(tumor_ PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)

install(TARGETS tumor_ 
	LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
