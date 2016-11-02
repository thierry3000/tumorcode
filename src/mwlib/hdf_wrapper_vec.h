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
#ifndef HDF_WRAPPER_VEC_H
#define HDF_WRAPPER_VEC_H

#include "helpers-vec.h"
#include "hdf_wrapper.h"

#if 1  //multiple declaration is always bad!!! before this was h5=h5cpp
namespace h5cpp_intern {

template<class T, int dim>
inline void set_array(h5cpp::Attributes attrs, const std::string &name, const Vec<T,dim> &vec)
{
  attrs.set(name, h5cpp::Dataspace::simple_dims(dim), vec.data());
}

template<class T, int dim>
inline Vec<T,dim> get_array(h5cpp::Attributes attrs, const std::string &name)
{
  Vec<T,dim> r;
  auto a = attrs.open(name);
  assert (a.get_dataspace().get_npoints() == dim);
  a.read<T>(r.data());
  return r;
}

}
#endif

#endif