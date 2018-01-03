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

/**
 * until 2018 tumorcode relied on some homebrewed hdfwrapper 
 * (https://github.com/DaWelter/HDF5-cpp-wrapper).
 * Because of numerous bugs, interface problem and inconvieniences,
 * I decided to change the api to the official supported HDF5 Cpp wrapper.
 * Probably this will take some time to behave well.
 */
#ifndef HDFIO_H
#define HDFIO_H


#include "mwlib/helpers-vec.h"

#include "mwlib/hdf_wrapper_array3d.h"
#include "mwlib/hdf_wrapper_ptree.h"
#include "mwlib/field.h"

#include <string>
#include <H5Cpp.h>
#include <H5Exception.h>



struct LatticeDataQuad3d;
class VesselList3d;
class ContinuumTumor;
struct SimParameters3d;
struct ContinuumTumorState;

void WriteHdfGraph( H5::Group f, const VesselList3d &vl );
void ReadHdfGraph( H5::Group f, VesselList3d *vl );


// template<class T>
// H5::DataSet WriteScalarField(H5::Group g, const string &name, ConstArray3d<T> arr, const LatticeDataQuad3d &ld, const h5cpp::Group ldgroup, const h5cpp::Datatype &disktype = h5cpp::get_disktype<T>());

template<class T>
H5::DataSet WriteScalarField(H5::Group g, const string &name, ConstArray3d<T> arr, const LatticeDataQuad3d &ld, const H5::Group ldgroup);

template<class Vector>
H5::DataSet WriteVectorField(H5::Group g, const string &name, ConstArray3d<Vector> arr, const LatticeDataQuad3d &ld, const H5::Group ldgroup);

template<class T>
H5::DataSet WriteAveragedFaceVariables(H5::Group file, const std::string &id, int dim, const ConstArray3d<T> *face_fields);

template<class T>
H5::DataSet WriteAveragedFaceVariableField(H5::Group file, const std::string &id, int dim, const ConstArray3d<T> *face_fields, const LatticeDataQuad3d &ld, const H5::Group ldgroup);

H5::Group RequireLatticeDataGroup(H5::H5File f, const string &name, const LatticeDataQuad3d &ld);
H5::Group RequireLatticeDataGroup(H5::Group g, const LatticeDataQuad3d &ld);

template<typename T>
T readAttrFromGroup(H5::Group g, string attr_name);

template<typename T>
T readAttrFromDataSet(H5::DataSet g, string attr_name);

/** @brief
 * write attributes to Group 
 */
void writeAttrToH5(H5::Group h, const string &attr_name,  const Int6 &value);
void writeAttrToH5(H5::Group h, const string &attr_name,  const Int3 &value);
void writeAttrToH5(H5::Group h, const string &attr_name,  const int &value);
void writeAttrToH5(H5::Group h, const string &attr_name,  const Float3 &value);
void writeAttrToH5(H5::Group h, const string &attr_name,  const float &value);
void writeAttrToH5(H5::Group h, const string &attr_name,  const Double3 &value);
void writeAttrToH5(H5::Group h, const string &attr_name,  const double &value);
void writeAttrToH5(H5::Group h, const string &attr_name,  const string &value);
void writeAttrToH5(H5::Group h, const string &attr_name,  const Bool3 &value);
void writeAttrToH5(H5::Group h, const string &attr_name,  const bool &value);
/** @brief
 * write attributes to DataSet 
 */
void writeAttrToH5(H5::DataSet h, const string &attr_name,  const Int6 &value);
void writeAttrToH5(H5::DataSet h, const string &attr_name,  const Int3 &value);
void writeAttrToH5(H5::DataSet h, const string &attr_name,  const int &value);
void writeAttrToH5(H5::DataSet h, const string &attr_name,  const Float3 &value);
void writeAttrToH5(H5::DataSet h, const string &attr_name,  const float &value);
void writeAttrToH5(H5::DataSet h, const string &attr_name,  const Double3 &value);
void writeAttrToH5(H5::DataSet h, const string &attr_name,  const double &value);
void writeAttrToH5(H5::DataSet h, const string &attr_name,  const string &value);
void writeAttrToH5(H5::DataSet h, const string &attr_name,  const Bool3 &value);
void writeAttrToH5(H5::DataSet h, const string &attr_name,  const bool &value);



template<class T>
H5::DataSet writeDataSetToGroup(H5::Group g, string attr_name, T value);

template<class T>
void readDataSetFromGroup(H5::Group, string attr_name, T *placeToStore);

template<class T>
string getH5Name(T g);
#endif
