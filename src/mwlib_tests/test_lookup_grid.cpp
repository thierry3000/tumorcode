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

#include <iostream>

#include "mwlib/lookup_grid.h"

template<class T>
boost::array<T, 2> make_array(T x, T y)
{
  boost::array<T, 2> p = {{x, y}};
  return p;
}

template<class T, std::size_t dim>
inline std::ostream& operator<<(std::ostream &os, const boost::array<T, dim> &a)
{
  os << "[";
  for (int i=0; i<dim; ++i){
    os << a[i];
    if (i < dim-1) os << ",";
  }
  os << "]";
  return os;
}

 
int main(int argc, char **argv)
{
  //std::cout << make_array<int>(1,2) << std::endl;
  typedef my::LookupGrid<int, 2>::Point Point;
  my::LookupGrid<int, 2> g;
  g.initialize(make_array<double>(-1., -1.), make_array<double>(1., 1.), 0.5);
  g(make_array<double>(0.,0.)) = 5;
  g.print(std::cout);
  std::cout << "point at 0., 0. = " << g(make_array<double>(0.,0.)) << std::endl;
  std::vector<my::LookupGrid<int, 2>::const_iterator> search_dat;
  
  g.region_search(make_array<double>(0., 0.), make_array<double>(10., 10.), search_dat);
  for (int i=0; i<search_dat.size(); ++i)
  {
    std::cout << "search " << g.cell_center(search_dat[i]) << " = " << *search_dat[i] << std::endl;
  }
} 
