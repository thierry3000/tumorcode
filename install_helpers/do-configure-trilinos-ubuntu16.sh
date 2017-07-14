#!/bin/sh
DST=/usr/local/trilinos
cmake \
      -D Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=OFF \
      -D Trilinos_ENABLE_Amesos:BOOL=ON \
      -D Trilinos_ENABLE_Ifpack:BOOL=ON \
      -D Trilinos_ENABLE_Belos:BOOL=ON \
      -D Trilinos_ENABLE_Tpetra:BOOL=ON \
      -D Trilinos_ENABLE_OpenMP:BOOL=ON \
      -D Trilinos_ENABLE_Sacado:BOOL=ON \
      -D Trilinos_ENABLE_ML:BOOL=ON \
      -D Trilinos_ENABLE_Epetra:BOOL=ON \
      -D Trilinos_ENABLE_EpetraExt:BOOL=ON \
      -D Trilinos_ENABLE_AztecOO:BOOL=ON \
      -D Trilinos_ENABLE_DEBUG_SYMBOLS:BOOL=ON \
      -D ML_ENABLE_Zoltan:BOOL=OFF \
      -D ML_ENABLE_Aztec:BOOL=ON \
      -D ML_ENABLE_Galeri:BOOL=OFF \
      -D ML_ENABLE_Isorropia:BOOL=OFF \
      -D CMAKE_BUILD_TYPE:STRING=Release \
      -D Trilinos_ENABLE_OpenMP:BOOL=ON \
      -D TPL_ENABLE_MPI:BOOL=OFF \
      -D TPL_ENABLE_Pthread:BOOL=ON \
      -D CMAKE_INSTALL_PREFIX:PATH=$DST \
      -D BUILD_SHARED_LIBS:BOOL=ON \
      $1
