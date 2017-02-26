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
#ifndef HDFIO_H
#define HDFIO_H

#include "mwlib/helpers-vec.h"
#include "hdf_wrapper.h"
#include "mwlib/hdf_wrapper_array3d.h"
#include "mwlib/hdf_wrapper_ptree.h"
#include "mwlib/field.h"

struct LatticeDataQuad3d;
class VesselList3d;
class ContinuumTumor;
struct SimParameters3d;
struct ContinuumTumorState;

void WriteHdfGraph( h5cpp::Group f, const VesselList3d &vl );
void ReadHdfGraph( h5cpp::Group f, VesselList3d *vl );


template<class T>
h5cpp::Dataset WriteScalarField(h5cpp::Group g, const string &name, ConstArray3d<T> arr, const LatticeDataQuad3d &ld, const h5cpp::Group ldgroup, const h5cpp::Datatype &disktype = h5cpp::get_disktype<T>());

template<class Vector>
h5cpp::Dataset WriteVectorField(h5cpp::Group g, const string &name, ConstArray3d<Vector> arr, const LatticeDataQuad3d &ld, const h5cpp::Group ldgroup);

template<class T>
h5cpp::Dataset WriteAveragedFaceVariables(h5cpp::Group file, const std::string &id, int dim, const ConstArray3d<T> *face_fields);

template<class T>
h5cpp::Dataset WriteAveragedFaceVariableField(h5cpp::Group file, const std::string &id, int dim, const ConstArray3d<T> *face_fields, const LatticeDataQuad3d &ld, const h5cpp::Group ldgroup);

h5cpp::Group RequireLatticeDataGroup(h5cpp::Group g, const string &name, const LatticeDataQuad3d &ld);

#endif
