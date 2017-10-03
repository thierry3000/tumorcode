# Introduction

This is the software implementing the models used by the Group of Prof. Heiko Rieger
at Saarland University to study vascularized tumour growth of tumor multi-cellular
spheroids, interstitial fluid flow, drug delivery and oxygenation.

We publish it here to foster good scientific practices, i.e. to facilitate
reproduction of our work by interested parties.

For more information see
http://www.uni-saarland.de/fak7/rieger/homepage/research/tumor/tumor.html

# Building
## Requirements
- trilinos
- tbb
- cimg
- boost
- HDF5-cpp-wrapper
(https://github.com/DaWelter/HDF5-cpp-wrapper)
- numpycpp
(https://github.com/DaWelter/numpycpp)

## Successfully built on
- 4.1.11-gentoo-ARCH
- Ubuntu 16.04 4.4.0-45-generic #66-Ubuntu

## Build process

CMake is our build system. You most likely need to configure the build according
to your system specifics. To do so we recommended to create a shell script with
the call to cmake. We have several examples in the `install_helpers` directory.

Basic CMake usage is as follows
```
$mkdir build  # create your build directory
$cd build
$cmake -DCMAKE_INSTALL_PREFIX=your_install_dir path_to_tumorcode
$make
$make install
```

The install command is required since it copies generated python libraries to
their proper place within the `py` subdirectory. Note that you can set `CMAKE_INSTALL_PREFIX`
to your source directory, which is useful for development.

## Details
For a detailed list of Ubuntu packages you will need, see the wiki.

External projects needed:


Our HDF5-cpp-wrapper is a lightweight c++ library around HDF5. Get it by
`git clone https://github.com/DaWelter/HDF5-cpp-wrapper.git` As it is only a
single-header you can simply put a link to it in a directory in your include paths.

NumpyCPP is our C++ wrapper around the numpy C API. Get it by
`git clone https://github.com/DaWelter/numpycpp.git`. Compile and install by
```
$cd numpycpp
$mkdir build
$cmake -DCMAKE_INSTALL_PREFIX=your_install_dir_defaulting_to_usr/local_or_something_like_that ../
$make
$make install
```

## Issues

### Linking problems with HDF5 and h5py
IMPORTANT: We fetch object references pertaining to the HDF5 library linked with
h5py and use functions on them pertaining to the HDF5 lib linked to our code,
assuming both libraries are the same. THIS MIGHT, HOWEVER, NOT BE THE CASE.

If so, on Ubuntu 16 for instance, you have to compile h5py by yourself. Here is
how you can achieve this:

Install the requirements `sudo apt-get install cython python-mpi4py python-pkgconfig`.
Then get the h5py source, e.g. from GitHub `git clone https://github.com/h5py`.
Then build and install by
```
$cd h5py
$python2 setup.py configure --hdf5=/usr/include/hdf5/openmpi --mpi
$python2 setup.py build
$sudo python2 setup.py install
# or if you want to install it in a custom directory
$python2 setup.py install --prefix=target_dir
```
By `--hdf5=...` and `--mpi` we request to link against the hdf5 library
that is installed by respective debian packages, which happen to be compiled with
MPI Support.

### MPI
This software performs all computations locally, as opposed to distributed HPC.
Therefore MPI is NOT required. However, it should do no harm if trilinos and
other packages are linked to MPI.

You might have to link against `libmpi_cxx.so` or do
`$export LD_PRELOAD=/usr/lib/libmpi_cxx.so` prior to running our code, though.

### Compiling Trilinos without MPI

We usually completely eliminate MPI from your build to avoid any potential for trouble.
Therefore we compile the Trilinos suite by ourselves. To aid people
in doing the same we provide configure scripts located in `install_helpers`.
