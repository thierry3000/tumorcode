/**
This file is part of tumorcode project.
(http://www.uni-saarland.de/fak7/rieger/homepage/research/tumor/tumor.html)

Copyright (C) 2016  Michael Welter and Thierry Fredrich

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef HDF_WRAPPER_ARRAY3D_H
#define HDF_WRAPPER_ARRAY3D_H


#include "hdf_wrapper_vec.h"
#include "field.h"
#include "lattice-data.h"


void WriteHdfLd( H5::Group f, const LatticeDataQuad3d &ld );
void ReadHdfLd( H5::Group f, LatticeDataQuad3d &ld );
void WriteHdfLd( H5::Group f, const LatticeDataFCC &ld );
void ReadHdfLd( H5::Group f, LatticeDataFCC &ld );


template<class T>
//H5::DataSet WriteArray3D(H5::Group file, const std::string &id, const ConstArray3d<T> &a, const H5::DataType &disktype = h5cpp::get_disktype<T>());
H5::DataSet WriteArray3D(H5::Group file, const std::string &id, const ConstArray3d<T> &a);

template<class Vector>
H5::DataSet WriteVectorArray3D(H5::Group  file,const std::string &id, const ConstArray3d<Vector> &a);

template<class T>
void ReadArray3D(H5::DataSet ds, Array3d<T> &a);

template<class Vector>
void ReadVectorArray3D(H5::DataSet ds, Array3d<Vector> &a);



#endif