This file is part of tumorcode project. For more information see
http://www.uni-saarland.de/fak7/rieger/homepage/research/tumor/tumor.html

Note: direct commands are quoted. We used KDevelop 4.7.2 as deloping framework.

*** code tested with 
- 4.1.11-gentoo-ARCH
- Ubuntu 16.04 4.4.0-45-generic #66-Ubuntu

*** Requirements
- trilinos
- tbb
- cimg
- boost
- HDF5-cpp-wrapper
(see https://github.com/DaWelter/HDF5-cpp-wrapper)
- numpycpp
(https://github.com/DaWelter/numpycpp)

If your sure your system fullfills all requirements
you can directly use the install script in:
"/install_helpers/do-configure-tumorcode-ubuntu16.sh"

Otherwise continue reading...

It is recommended to adapt this to your own system.

1) Adapt the scripts to your sytem e.g. create the installation directory and make 
   sure you have the access rights 

2) Usage of do-configure.sh

*) create build directory
	  */build   or  */buildopt
	  
*) cd into directory

*) run do-configure.sh
syntax:                   ". do-configure-tumorcode-lusi.sh ../"


****
****  Ubuntu 16.04
****
We set up a desktop PC with Ubuntu 16.04.
The system contains 2 partions. One for the system and one at /localdisk/ .
If you decide to use different locations adapt the scripts before running them.

*************Packages needed:
We assume you got this file by git. So this is already there ;-)

From Ubuntu source --> 
**general system
"sudo apt-get install \
cmake \
gfortran \
libblas-dev \
liblapack-dev \
python2.7-dev \
libboost-dev \
libboost-python-dev \
python-numpy \
python-matplotlib \
libboost-program-options-dev \
libpng16-dev \
libtbb-dev \
libeigen3-dev \
python-h5py \
python-scipy \
povray \
mpi-default-dev \
libmumps-dev \
libsuperlu-dev \
libhdf5-openmpi-dev \
libptscotch-dev \
binutils-dev \
libiberty-dev"

**trilinos
"sudo apt-get install libtrilinos-*"

****** HDF5-cpp-wrapper  (see https://github.com/DaWelter/HDF5-cpp-wrapper):
"git clone https://github.com/DaWelter/HDF5-cpp-wrapper.git"

****** numpycpp  (see https://github.com/DaWelter/numpycpp) :
"git clone https://github.com/DaWelter/numpycpp.git"
"cd numpycpp"
"mkdir build"
"cmake ../"
"make"

****** make h5py parallel
Unfortunatelly ubuntus python-h5py package is linked against the serial version of
libhdf5.so. Therefore we need to build, and install our own h5py with mpi support.

- get additional packages 
"sudo apt-get install cython python-mpi4py python-pkgconfig"
- get h5py source, e.g. 
"wget https://pypi.python.org/packages/22/82/64dada5382a60471f85f16eb7d01cc1a9620aea855cd665609adf6fdbb0d/h5py-2.6.0.tar.gz#md5=ec476211bd1de3f5ac150544189b0bf4"
-untar
"tar xvf h5py-2.6.0.tar.gz"
"cd h5py-2.6.0"
"python2 setup.py configure --hdf5=/usr/include/hdf5/openmpi --mpi"
"python2 setup.py build"
"sudo python2 setup.py install" 
or if you do not have administer priviledges
"python2 setup.py install --prefix="where/ever/you/have/rights""


****** tumorcode:
*) create install dir (e.g. mkdir /localdisk/tc_install )
*) create build dir (e.g. mkdir /localdisk/buildopt )
*) change to build dir (cd /localdisk/buildopt )
*) check paths in do-configure-tumorcode-ubuntu16.sh !!!!
*) run ". ../tumorcode/install_helpers/do-configure-tumorcode-ubuntu16.sh ../tumorcode"
*) run "make install"
*) enviroment variables for convenience "export PATH=$PATH:/localdisk/tc_install/bin"

*** BEFORE EXECUTE --- tell python to use mpi!!!!
"export LD_PRELOAD=/usr/lib/libmpi_cxx.so"



****** trilinos: (from source and WITHOUT MPI support)
Note: to get trilinos source you either need to register
on the website or use your github.com account and 
get it from there.

*) create install dir (e.g. "/localdisk/trilinos12")
*) get trilinos-12.6.4-Source.tar.gz from https://trilinos.org/download/previous-releases/download-12-6/
"tar xvf trilinos-12.6.4-Source.tar.gz"
"cd trilinos-12.6.4-Source"
"mkdir build"
"cd build"
". /tumorcode/do-configure-trilinos-ubuntu16.sh .. /localdisk/trilinos12"
"make install -j 8"







