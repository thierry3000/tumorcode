cmake \
    -DCMAKE_INSTALL_PREFIX=/localdisk/tc_install \
    -DCMAKE_BUILD_TYPE=Release \
    -DPYTHON_NUMPY_INCLUDE_DIR=/usr/include/python2.7/numpy \
    -DADDITIONAL_INCLUDE_DIRS="" \
    -DADDITIONAL_LIBRARY_DIRS="" \
    -DHDF5_INCLUDE_DIRS="/usr/include/hdf5/openmpi" \
    $1
