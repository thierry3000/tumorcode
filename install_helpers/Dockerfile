# syntax=docker/dockerfile:1
FROM ubuntu:18.04
MAINTAINER sofija
ARG DEBIAN_FRONTEND=noninteractive
RUN mkdir /localdisk
RUN cd /localdisk
RUN apt-get update
RUN echo $SHELL
RUN yes | apt-get install git
RUN yes | apt-get update && yes | apt-get install build-essential cmake trilinos-dev libtbb-dev cimg-dev libboost-dev libhdf5-openmpi-dev
RUN yes | apt-get install cmake gfortran libblas-dev liblapack-dev python2.7-dev libboost-dev libboost-python-dev libboost-numpy-dev python-numpy python-matplotlib libboost-program-options-dev libtbb-dev libeigen3-dev python-h5py python-scipy povray mpi-default-dev libmumps-dev libsuperlu-dev libhdf5-openmpi-dev libptscotch-dev binutils-dev libiberty-dev libtrilinos-* python3-numpy libhdf5-dev libboost-date-time-dev libboost-system-dev libboost-serialization-dev cimg-dev libopenmpi-dev lsb-release libboost-thread-dev libboost-filesystem-dev libpng-dev libboost-random-dev
RUN yes | mkdir /root/.ssh
RUN eval "$(ssh-agent -s)"
RUN ssh-keygen -q -t rsa -N '' -f /root/.ssh/id_rsa
RUN exec ssh-agent bash
RUN echo '''[url https://github.com/]    insteadOf = <ssh://git@github.com>''' ~/.gitconfig
RUN yes | ssh-keyscan -t rsa github.com >> /root/.ssh/known_hosts
RUN yes | apt-get install cython python-mpi4py python-pkgconfig
RUN yes | git clone --branch dev_fredrich https://github.com/thierry3000/tumorcode.git
RUN mkdir /buildopt
RUN chmod +x /tumorcode/install_helpers/do-configure-tumorcode-ubuntu18.sh
RUN pwd
RUN cd /buildopt && ../tumorcode/install_helpers/do-configure-tumorcode-ubuntu18.sh /tumorcode/
RUN cd /buildopt && make && make install
CMD ["/localdisk/tc_install/bin/submitVesselgeneration","-h"]
CMD ["echo", "Hello World"]
