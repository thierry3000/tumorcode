#this add_custom_target stuff works only at build times!!!!!
#http://stackoverflow.com/questions/34876602/cmake-create-symlink-vs-ln
# add_custom_target(make_bin_dir ALL COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_INSTALL_PREFIX}/bin)
# add_custom_target(submitAdaption ALL COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_INSTALL_PREFIX}/py/krebsjobs/submitAdaption.py ${CMAKE_INSTALL_PREFIX}/bin/submitAdaption)
# add_custom_target(submitIff ALL COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_INSTALL_PREFIX}/py/krebsjobs/submitIff.py ${CMAKE_INSTALL_PREFIX}/bin/submitIff)
# add_custom_target(submitDetailedO2 ALL COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_INSTALL_PREFIX}/py/krebsjobs/submitDetailedO2.py ${CMAKE_INSTALL_PREFIX}/bin/submitDetailedO2)
# add_custom_target(submitVesselgeneration ALL COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_INSTALL_PREFIX}/py/krebsjobs/submitVesselgeneration.py ${CMAKE_INSTALL_PREFIX}/bin/submitVesselgeneration)
# add_custom_target(submitPreziosi ALL COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_INSTALL_PREFIX}/py/krebsjobs/submitPreziosi.py ${CMAKE_INSTALL_PREFIX}/bin/submitPreziosi)
# add_custom_target(submitFakeTum ALL COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_INSTALL_PREFIX}/py/krebsjobs/submitFakeTum.py ${CMAKE_INSTALL_PREFIX}/bin/submitFakeTum)
# add_custom_target(submitPovrayRender ALL COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_INSTALL_PREFIX}/py/krebsjobs/submitPovrayRender.py ${CMAKE_INSTALL_PREFIX}/bin/submitPovrayRender)

message(STATUS "Creating some symlink to make life easier!")
message(STATUS "see PostInstall.cmake")
execute_process(COMMAND ln -s ${CMAKE_INSTALL_PREFIX}/py/krebsjobs/submitAdaption.py ${CMAKE_INSTALL_PREFIX}/bin/submitAdaption)
execute_process(COMMAND ln -s ${CMAKE_INSTALL_PREFIX}/py/krebsjobs/submitIff.py ${CMAKE_INSTALL_PREFIX}/bin/submitIff)
execute_process(COMMAND ln -s ${CMAKE_INSTALL_PREFIX}/py/krebsjobs/submitDetailedO2.py ${CMAKE_INSTALL_PREFIX}/bin/submitDetailedO2)
execute_process(COMMAND ln -s ${CMAKE_INSTALL_PREFIX}/py/krebsjobs/submitVesselgeneration.py ${CMAKE_INSTALL_PREFIX}/bin/submitVesselgeneration)
execute_process(COMMAND ln -s ${CMAKE_INSTALL_PREFIX}/py/krebsjobs/submitBulkTissue.py ${CMAKE_INSTALL_PREFIX}/bin/submitBulkTissue)
execute_process(COMMAND ln -s ${CMAKE_INSTALL_PREFIX}/py/krebsjobs/submitFakeTum.py ${CMAKE_INSTALL_PREFIX}/bin/submitFakeTum)
execute_process(COMMAND ln -s ${CMAKE_INSTALL_PREFIX}/py/krebsjobs/submitPovrayRender.py ${CMAKE_INSTALL_PREFIX}/bin/submitPovrayRender)

