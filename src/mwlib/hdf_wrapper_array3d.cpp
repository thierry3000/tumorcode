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
#include "hdf_wrapper_array3d.h"

template<class T>
H5::DataSet WriteArray3D(H5::Group file, const std::string &DATASET_NAME, const ConstArray3d<T> &a)
{
  
  const Int3 s = a.size();
  const int rank = 3;
  hsize_t     dims[rank];       // dataset dimensions
  dims[0] = s[0];
  dims[1] = s[1];
  dims[2] = s[2];
  H5::DataSpace dspace = H5::DataSpace( rank, dims);
  H5::DataSet dataset = file.createDataSet(DATASET_NAME, H5::PredType::NATIVE_FLOAT,dspace);
  // Write the data to the dataset using default memory space, file
	// space, and transfer properties.
  dataset.write(&a, H5::PredType::NATIVE_FLOAT);
  return dataset;
//     h5cpp::Dataspace dspace = h5cpp::Dataspace::simple_dims(s[0],s[1],s[2]);
//     Array3d<T> tmp(Int3(s[2],s[1],s[0]));
//     tmp.swapAxes(0,2);
//     tmp.fill(a);
//     return h5cpp::create_dataset<T>(file, id, dspace, tmp.getPtr(), h5cpp::CREATE_DS_COMPRESSED);
}


template<class Vector>
H5::DataSet WriteVectorArray3D(H5::Group  file,const std::string &id, const ConstArray3d<Vector> &a)
{
  const Int3 s = a.size();
//     h5cpp::Dataspace dspace = h5cpp::Dataspace::simple_dims(s[0],s[1],s[2], sizeof(Vector)/sizeof(typename Vector::value_type));
//     Array3d<Vector> tmp(Int3(s[2],s[1],s[0]));
//     tmp.swapAxes(0,2);
//     tmp.fill(a);
//     return h5cpp::create_dataset(file, id, dspace, (typename Vector::value_type const*)tmp.getPtr(), h5cpp::CREATE_DS_COMPRESSED);
}



template<class T>
void ReadArray3D(H5::DataSet ds, Array3d<T> &a)
{
//   h5cpp::Dataspace sp = ds.get_dataspace();
//   hsize_t dims[H5S_MAX_RANK];
//   int r = sp.get_dims(dims);
//   if (r != 3) throw h5cpp::Exception("ReadArray3d expected Dataset of Rank 3!");
//   Array3d<T> tmp(Int3(dims[2], dims[1], dims[0]));
//   tmp.swapAxes(0,2);
//   ds.read<>(tmp.getPtr());
//   a = Array3d<T>(tmp.size());
//   a.fill(tmp);
}


template<class Vector>
void ReadVectorArray3D(H5::DataSet ds, Array3d<Vector> &a)
{
//   h5cpp::Dataspace sp = ds.get_dataspace();
//   hsize_t dims[H5S_MAX_RANK];
//   int r = sp.get_dims(dims);
//   if (r != 4) throw h5cpp::Exception("ReadVectorArray3d expected Dataset of Rank 4.");
//   if (dims[3] != Vector::SizeAtCompileTime) throw h5cpp::Exception("ReadVectorArray3d expected Dataset dimension 3 to be the same size as the Vector type.");
//   Array3d<Vector> tmp(Int3(dims[2], dims[1], dims[0]));
//   tmp.swapAxes(0,2);
//   ds.read<>((typename Vector::value_type*)tmp.getPtr());
//   a = Array3d<Vector>(tmp.size());
//   a.fill(tmp);
}

/*
#define INSTANTIATE(T)\
  template h5cpp::Dataset WriteArray3D<T>(h5cpp::Group file, const std::string &id, const ConstArray3d<T> &a, const h5cpp::Datatype &disktype);\
  template void ReadArray3D<T>(h5cpp::Dataset ds, Array3d<T> &a);

#define INSTANTIATE_VEC(T)\
  template h5cpp::Dataset WriteVectorArray3D<Vec<T,3> >(h5cpp::Group  file,const std::string &id, const ConstArray3d<Vec<T,3> > &a);\
  template h5cpp::Dataset WriteVectorArray3D<Vec<T,2> >(h5cpp::Group  file,const std::string &id, const ConstArray3d<Vec<T,2> > &a);\
  template h5cpp::Dataset WriteVectorArray3D<Vec<T,1> >(h5cpp::Group  file,const std::string &id, const ConstArray3d<Vec<T,1> > &a);\
  template void ReadVectorArray3D<Vec<T,3> >(h5cpp::Dataset ds, Array3d<Vec<T,3> > &a);

INSTANTIATE(float)
INSTANTIATE(double)
INSTANTIATE(int)
INSTANTIATE(char)
INSTANTIATE_VEC(float)
INSTANTIATE_VEC(double)
*/

#define INSTANTIATE(T)\
  template H5::DataSet WriteArray3D<T>(H5::Group file, const std::string &id, const ConstArray3d<T> &a);\
  template void ReadArray3D<T>(H5::DataSet ds, Array3d<T> &a);

#define INSTANTIATE_VEC(T)\
  template H5::DataSet WriteVectorArray3D<Vec<T,3> >(H5::Group  file,const std::string &id, const ConstArray3d<Vec<T,3> > &a);\
  template H5::DataSet WriteVectorArray3D<Vec<T,2> >(H5::Group  file,const std::string &id, const ConstArray3d<Vec<T,2> > &a);\
  template H5::DataSet WriteVectorArray3D<Vec<T,1> >(H5::Group  file,const std::string &id, const ConstArray3d<Vec<T,1> > &a);\
  template void ReadVectorArray3D<Vec<T,3> >(H5::DataSet ds, Array3d<Vec<T,3> > &a);

INSTANTIATE(float)
INSTANTIATE(double)
INSTANTIATE(int)
INSTANTIATE(char)
INSTANTIATE_VEC(float)
INSTANTIATE_VEC(double)