cmake_minimum_required(VERSION 2.6)

#include_directories(../adaption)
#file(GLOB CommonFiles "*.cpp")

add_library(common STATIC
bulktissuemodel1_new.cpp
calcflow.cpp
calcflow_linsys.cpp
cell_based_oxygen_update_model.cpp
common.cpp
continuum-flow.cpp
continuum-grid.cpp
continuum-surfacetension.cpp
continuum-utils.cpp
convection_diffusion_solver.cpp
distancemap.cpp
glucose_model.cpp
growthfactor_model.cpp
hdfio.cpp
lattice-data-polymorphic.cpp
levelset.cpp
shared-object.cpp
simple_oxygen_model.cpp
trilinos_helper.cpp
trilinos_linsys_construction.cpp
vesselmodel1.cpp
vessels3d_biconnected.cpp
vessels3d.cpp
)
# Das wuerde ich deswegen in das Haupt-CMakeLists skript verschieben.
target_link_libraries(common ${HDF5_LIBRARIES})
target_link_libraries(common mwlib)



file(GLOB CommonFiles2 "vesselgen/*.cpp")
add_library(vesselgenlib STATIC ${CommonFiles2})

install(DIRECTORY . DESTINATION "include/common" FILES_MATCHING PATTERN "*.h")
