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

#include "field.h"


template<class T>
static DynArray<T> ComputeConvolutionValues(ConstArray3d<T> &stencil, ConstArray3d<bool> mask)
{
  const BBox3 stencil_box = stencil.getBox();
  int n = Volume(stencil_box);
  DynArray<T> values(n);
  int k = 0;
  FOR_BBOX3(p, stencil_box)
  {
    if (!mask.empty() && !mask(p)) continue;
    values[k] = stencil(p);
    ++k;
  }
  values.resize(k);
  return values;
}


template<class T>
void PrintArray(const ConstArray3d<T> a, std::ostream &os)
{
  using std::endl;
  os << "<<size=" << a.size() <<  ", box=" << a.getBox() << "," << endl;
  const BBox3 b = a.getBox();
  Int3 p;
  os << "[" << endl;
  for(p.z()=b.min[2]; p.z()<=b.max[2]; ++p.z())
  {
    os << " [" << endl;
    for(p.y()=b.min[1]; p.y()<=b.max[1]; ++p.y())
    {
      os << "  [ ";
      for(p.x()=b.min[0]; p.x()<=b.max[0]; ++p.x())
      {
        os << a(p) << " ";
      }
      os << "]" << endl;
    }
    os << " ]" << endl;
  }
  os << "]>>" << endl;
}


