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
#ifndef DISTANCEMAP_H
#define DISTANCEMAP_H

#include "python_krebsutils/python_helpers.h"

#include "mwlib/dynamicarray.h"
#include "mwlib/helpers-containers.h" // for the heap
#include "mwlib/helpers-vec.h"
#include "mwlib/field.h"


struct DistanceFieldComputer
{
  typedef LatticeDataQuad3d LD;
  LD ld;
  int dim;
  const int FLAG_OUTER;
  const int FLAG_INNER;
  const float DIST_MAX;
  
  struct HeapLessFunc {
    Array3d<float> f;
	  inline bool operator()( const Int3 &a, const Int3 &b ) { return f(a)<f(b); }
  };
  struct HeapSwapFunc {
    Array3d<int> flags;
	  inline void operator()( Int3 &a, Int3 &b ) { Int3 c(a); a=b; b=c;  int x=flags(a); flags(a)=flags(b); flags(b)=x;  myAssert(flags(a)>=0 && flags(b)>=0); }
  };    
  
  Heap<Int3,HeapLessFunc,HeapSwapFunc>  heap;
  Array3d<float> f;
  Array3d<int> flags;
    
  inline void InsertActive( const Int3 &p );
  inline void RemoveActive( const Int3 &p );
  inline float SolveDistance(float q, float x, float y, float z, bool bx, bool by, bool bz) const;
  inline void CheckOkDir(const Int3 &pp, int dir, bool &ok, float &val) const;
  float ComputeCenterGradient( const Int3 &p );
  
public:    
  DistanceFieldComputer();    
  void Do(const LatticeDataQuad3d &ld_, Array3d<float> &f );
#if BOOST_VERSION>106300
  void Do(const LatticeDataQuad3d &ld_, np::ndarray &f );
#endif
};


#endif
