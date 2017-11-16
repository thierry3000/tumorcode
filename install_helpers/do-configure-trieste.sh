cmake \
    -DCMAKE_INSTALL_PREFIX=/home/thierry/tc_install \
    -DCMAKE_BUILD_TYPE=Release \
    -DPYTHON_NUMPY_INCLUDE_DIR=/usr/include/python2.7/numpy \
    -DADDITIONAL_INCLUDE_DIRS="/home/thierry/vbl_install/include;/home/thierry/CImg-2.0.5_pre100217" \
    -DADDITIONAL_LIBRARY_DIRS="/home/thierry/vbl_install/lib" \
    $1
