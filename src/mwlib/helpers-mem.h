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
#ifndef HELPERS_MEM_H
#define HELPERS_MEM_H

#include "helpers-defs.h"
#include "myAssert.h"

#include <malloc.h>

// calling the constructor on an uninitialized piece of mem
#ifdef WIN32
#ifndef __PLACEMENT_NEW_INLINE
#  define __PLACEMENT_NEW_INLINE
  inline void* operator new(size_t, void* __p) throw() { return __p; }
  inline void  operator delete  (void*, void*) throw() { }
#endif

#ifndef __PLACEMENT_VEC_NEW_INLINE
#  define __PLACEMENT_VEC_NEW_INLINE
  inline void* operator new[](size_t, void* __p) throw() { return __p; }
  inline void  operator delete[](void*, void*) throw() { }
#endif
#endif


void _SetMem( void *ptr, int c, size_t n );
void _CopyMem( const void *src, void *dst, size_t n );


template< class T >
inline void ClearMem( T* ptr, size_t n )
{
  _SetMem( ptr, 0, n*sizeof(T) );
}

template< class T >
inline void CopyMem( const T* src, T* dst, size_t n )
{
  _CopyMem( src, dst, n*sizeof(T) );
}

inline const void* MovePtr( const void* p, ptrdiff_t bytes_offset )
{
  return ((char*)p)+bytes_offset;
}

inline void* MovePtr( void* p, ptrdiff_t bytes_offset )
{
  return ((char*)p)+bytes_offset;
}

/*----------------------------------------------------------
  Warpers for Placement New / Delete
  ---------------------------------------------------------*/

template<class T> inline void Construct( void *p ) { new(p) T(); }
template<class T, class U> inline void Construct( void *p, const U &init ) { new(p) T(init); }
template<class T> inline void ConstructVec( void *p, size_t n )
{
  T *end,*pp = (T*)(p);
  end = pp + n;
  while( pp < end )
  {
    Construct<T>( pp++ );
  }
}
template<class T> inline T* ConstructVec( void *p, size_t n, const T &init )
{
  T *end,*pp = (T*)(p);
  end = pp + n;
  while( pp < end )
  {
    Construct<T>( pp++, init );
  }
  return (T*)p;
}
template<class U, class T> inline U* ConstructVecCopy( void *p, size_t n, const T* init )
{
  U *end,*pp = (U*)(p);
  end = pp + n;
  while( pp < end )
  {
    Construct<U,T>( pp++, *(init++) );
  }
  return (U*)p;
}
template<class T> inline void Destruct( T &x )
{
	x.~T();
}
template<class T> inline void DestructVec( T *p, size_t n )
{
	size_t i;
	for (i=0;i<n;i++)
	{
		Destruct(*(p++));
  }
}



#endif
