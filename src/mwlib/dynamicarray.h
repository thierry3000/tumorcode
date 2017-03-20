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
#ifndef _DYNAMIC_ARRAY_H
#define _DYNAMIC_ARRAY_H

#include "vector_ext.h"
#include "myAssert.h"
#include "boost/serialization/access.hpp"
#include "boost/serialization/base_object.hpp"

namespace ConsTags
{
  constexpr struct RESERVE_t {} RESERVE;
  constexpr struct SHARED_t {} SHARED;
  constexpr struct COPY_t {} COPY;
  constexpr struct LIKE_t {} LIKE;
}

template <class T>
class DynArray : public std::vector<T>
{
  public:
    typedef std::vector<T> Base;
    DynArray() : Base() {}
    template<class U>
    DynArray(const DynArray<U> &other) : Base(other) {}
    DynArray(std::size_t size, const T &filler = T()) : Base(size, filler) {}
    DynArray(std::size_t size, ConsTags::RESERVE_t) : Base() { Base::reserve(size); }

    using Base::iterator;
    using Base::const_iterator;

    //const T& atmod(int i) const { myAssert(size()>0); return at(Mod(i,size())); }
    //T& atmod(int i)             { myAssert(size()>0); return at(Mod(i,size())); }

    typename Base::const_reference operator[](std::size_t i) const { myAssert(i<Base::size()); return Base::operator[](i); }
    typename Base::reference operator[](std::size_t i) { myAssert(i<Base::size()); return Base::operator[](i); }

    void remove_all() { Base::clear(); } // relies upon std::vector not freeing memory on clear
    void fill(const T &x) { std::fill(Base::begin(), Base::end(), x); }
};


template <class T>
inline std::size_t estimateMemoryUsage(const std::vector<T> &a)
{
  if (a.size()<=0) return sizeof(std::vector<T>);
  return sizeof(std::vector<T>) + a.size() * estimateMemoryUsage(get_ptr(a)[0]);
}



#endif
