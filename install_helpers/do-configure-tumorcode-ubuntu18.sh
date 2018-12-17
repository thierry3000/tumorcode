cmake \
    -DCMAKE_INSTALL_PREFIX=/localdisk/tc_install \
    -DCMAKE_BUILD_TYPE=Release \
    -DADDITIONAL_INCLUDE_DIRS="" \
    -DADDITIONAL_LIBRARY_DIRS="" \
    -DHDF5_INCLUDE_DIRS="/usr/include/hdf5/openmpi" \
    -DCMAKE_CXX_COMPILER="/usr/bin/mpic++" \
    $1
