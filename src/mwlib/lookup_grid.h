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
#ifndef LOOKUP_GRID_H
#define LOOKUP_GRID_H

#include <vector>

#include <boost/range/begin.hpp>
#include <boost/range/end.hpp>
#include <boost/array.hpp>
#include <boost/mpl/int.hpp>

#include "math_ext.h"

namespace my
{


template<class T, int dim>
class LookupGrid
{
  public:
    typedef boost::array<double, dim> Point;
    typedef boost::array<int, dim> Index;

    void initialize(const Point &a_, const Point &b_, double cell_size_)
    {
      for(int i=0; i<dim; ++i)
      {
        cellsize[i] = cell_size_;
        ncells[i] = my::iceil((b_[i] - a_[i]) / cellsize[i]);
        double l = ncells[i] * cellsize[i];
        double f = l - (b_[i] - a_[i]);
        a[i] = a_[i] - 0.5 * f;
        b[i] = a[i] + l;
        scale[i] = 1./cellsize[i];
      }

      unsigned int ncells_total = 1;
      for (int i=0; i<dim; ++i) ncells_total *= ncells[i];

      data.clear();
      data.resize(ncells_total);
    }

    T& operator()(const Point &p)
    {
      return data[index_to_site(point_to_index(p))];
    }

    const T& operator()(const Point &p) const
    {
      return data[index_to_site(point_to_index(p))];
    }

    void print(std::ostream &os) const
    {
      os << "LookupGrid ";
      for (int i=0; i<dim; ++i)
      {
        os << ncells[i];
        if (i < dim-1) os << "x";
      }
      os << std::endl;
      for (int i=0; i<dim; ++i)
        os << " " << a[i] << " - " << b[i] << std::endl;
      Index coord;
      for (int i=0; i<dim; ++i) coord[i] = 0;
      print_rec(os, dim-1, coord);
    }

    typedef typename std::vector<T>::size_type size_type;
    typedef typename std::vector<T>::const_iterator const_iterator;
    typedef typename std::vector<T>::iterator iterator;
    const_iterator begin() const { return data.begin(); }
    const_iterator end() const { return data.end(); }
    iterator begin() { return data.begin(); }
    iterator end() { return data.end(); }
    size_type size() const { return data.size(); }

    Point cell_center(const_iterator it) const
    {
      return index_to_cell_center(iterator_to_index(it));
    }

    void region_search(const Point &a, const Point &b, std::vector<const_iterator> &search_data) const
    {
      Index ia = point_to_index(a);
      Index ib = point_to_index(b);
      for (int i=0; i<dim; ++i)
      {
        assert(ia[i] <= ib[i]);
        ia[i] = std::max(ia[i], 0);
        ib[i] = std::min(ib[i], ncells[i]-1);
      }
      Index c;
      rec_search(c, ia, ib, boost::mpl::int_<dim-1>(), search_data);
    }

    void region_search(const Point &a, const Point &b, std::vector<iterator> &search_data)
    {
      region_search(a, b, reinterpret_cast<std::vector<const_iterator>&>(search_data));
    }

  private:
    Point a, b, scale, cellsize;
    Index ncells;
    std::vector<T> data;

    Point index_to_cell_center(const Index &coord) const
    {
      Point p;
      for (int i=0; i<dim; ++i)
        p[i] = a[i] + (coord[i]+0.5) * cellsize[i];
      return p;
    }

    Index iterator_to_index(const_iterator it) const
    {
      size_type s = it - data.begin();
      Index c;
      for(int i=0; i<dim; ++i)
      {
        c[i] = s % ncells[i];
        s = s / ncells[i];
      }
      return c;
    }

    Index point_to_index(const Point &p) const
    {
      Index coord;
      for(int i=0; i<dim; ++i)
        coord[i] = int((p[i] - a[i]) * scale[i]);
      return coord;
    }

    int index_to_site(const Index &coord) const
    {
      int s = 0;
      for(int i=dim-1; i>=0; --i)
        s = s * ncells[i] + coord[i];
      return s;
    }

    void print_rec(std::ostream &os, int axis, Index &coord) const
    {
      if (axis >= 0)
      {
        for (coord[axis]=0; coord[axis]<ncells[axis]; ++coord[axis])
          print_rec(os, axis-1, coord);
      }
      else
      {
        for (int i=0; i<dim; ++i)
        {
          os << coord[i];
          if (i < dim-1) os << " ";
        }
        os << ":" << data[index_to_site(coord)] << std::endl;
      }
    }

    void rec_search(Index &c, const Index &a, const Index &b, boost::mpl::int_<-1>, std::vector<const_iterator> &search_data) const
    {
      search_data.push_back(data.begin() + index_to_site(c));
    }

    template<int axis>
    void rec_search(Index &c, const Index &a, const Index &b, boost::mpl::int_<axis>, std::vector<const_iterator> &search_data) const
    {
      for(c[axis] = a[axis]; c[axis] <= b[axis]; ++c[axis])
      {
        rec_search(c, a, b, boost::mpl::int_<axis-1>(), search_data);
      }
    }
};

}

#endif
