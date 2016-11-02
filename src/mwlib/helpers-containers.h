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
#ifndef HELPERS_CONTAINERS_H
#define HELPERS_CONTAINERS_H

#include "helpers-defs.h"
#include "helpers-mem.h"
#include "dynamicarray.h"

#include <vector>
#include <limits>

using std::vector;

#include <stdlib.h>



template< class V >
inline void SwapRemove( V &v, size_t i ) {
  v[i] = v.back();
  v.pop_back();
}

template< class V >
inline typename V::value_type GetSwapRemoved( V &v, size_t i ) {
  typename V::value_type tmp( v[i] );
  v[i] = v.back();
  v.pop_back();
  return tmp;
}


/*---------------------------------------*/
// helpers for very small arrays
/*---------------------------------------*/

namespace StatArray
{
  template<class T, class C>
  inline void remove( T* _data, C &_cnt, int index ) {
    myAssert( index>=0 );
    myAssert( index<_cnt );
    for( int i=index+1; i<_cnt; ++i ) _data[i-1]=_data[i];
    --_cnt;
    _data[_cnt] = T();
  }

  template<class T, class C>
  inline int find( const T* const _data, const C _cnt, const T& x)
  {
    int i;
    for( i=0; i<_cnt; ++i ) if( _data[i]==x ) return i;
    return -1;
  }

  template<class T, class C>
  inline void push_back( T* _data, C& _cnt, const T& x)
  {
    _data[_cnt++] = x;
  }

  template<class T, class C>
  std::ostream& printArray(std::ostream &os, T* _data, C _cnt )
  {
    os << "{ " << std::endl;
    for( int i=0; i<_cnt; ++i )
    {
      os << "["<<i<<"] = " << _data[i] << std::endl;
    }
    os << "}" << std::endl;
    return os;
  }
}

namespace SmArray
{
  template<class T, class C>
  inline void resize( T* &_data, C& _cnt, size_t newcnt, bool clearNew = true )
  {
    myAssert( newcnt>=0 && newcnt<std::numeric_limits<C>::max() );
    if(_cnt==newcnt) return;

    if( _data==NULL ) {
      _data = (T*)malloc( newcnt * sizeof(T) );
    } else {
      _data = (T*)realloc( _data, newcnt * sizeof(T) );
    }
    if( !( _data || newcnt==0 ) ) { fprintf(stderr,"Error: Out of Mem\n"); exit(0); }

    if( clearNew && newcnt>_cnt )
      memset( _data+_cnt, 0, (newcnt-_cnt)*sizeof(T) );

    _cnt  = newcnt;
  }

  template<class T, class C>
  inline void remove( T* &_data, C& _cnt, int index )
  {
    StatArray::remove( _data, _cnt, index );
    resize<T,C>( _data, _cnt+1, _cnt );
  }

  template<class T, class C>
  inline void push_back( T* &_data, C& _cnt, const T& x)
  {
    resize( _data, _cnt, _cnt+1 );
    _data[_cnt-1] = x;
  }

  template<class T, class C>
  inline void free( T* &_data, C& _cnt )
  {
    free( _data );
    _cnt = 0;
  }

  template< class T, class C>
  inline void copyFrom( T* &_data, C& _cnt, const T* src, C srccnt )
  {
    if( _data ) free( _data );

    _cnt = srccnt;
    _data = (T*)malloc( _cnt * sizeof(T) );
    if( !( _data || _cnt==0 ) ) { fprintf(stderr,"Error: Out of Mem\n"); exit(0); }

    memcpy( _data, src, _cnt * sizeof(T) );
  }
}


/*---------------------------------------*/
// link list helper
/*---------------------------------------*/

// append b after a, updating begin and end also. a may be NULL, b must not be NULL
#define DBLLL_APPEND(T, _a, _b, next, prev, _begin, _end)\
{\
  T* a = _a;\
  T* b = _b;\
  T* &begin = _begin;\
  T* &end = _end;\
  T* c = a ? a->next : NULL;\
  if(a) a->next = b;\
  b->next = c;\
  if(c) c->prev = b;\
  b->prev = a;\
  if(end == a) end = b;\
  if(begin==NULL) begin = b;\
}

#define DBLLL_REMOVE(T, b, next, prev, begin, end)\
{\
  T* b = _b;\
  T* &begin = _begin;\
  T* &end = _end;\
  T* a = b->prev;\
  T* c = b->next;\
  if(a) a->next = c;\
  if(c) c->prev = a;\
  if(begin == b) begin=c;\
  if(end == b) end=a;\
}



template<class ThisType>
class DblLLNode
{
public:
    ThisType *next,*prev;
    DblLLNode() : next(NULL),prev(NULL) {}

    inline const ThisType* Next() const { return next; }
    inline const ThisType* Prev() const { return prev; }
    inline ThisType* Next() { return next; }
    inline ThisType* Prev() { return prev; }
    inline ThisType* This() { return static_cast<ThisType*>(this); }

	  inline void Append(ThisType *e)
	  {
	    e->next = next;
	    e->prev = This();
	    if(next) next->prev=e;
	    next = e;
	  }
	  inline void Append(ThisType *e, ThisType* &end)
	  {
	    Append( e );
	    if( end == This() ) end = e;
	  }

	  inline void Prepend(ThisType *e)
	  {
	    e->prev = prev;
	    e->next = This();
	    if(prev) prev->next=e;
	    prev = e;
	  }
	  inline void Prepend(ThisType *e, ThisType* &begin)
	  {
	    Prepend( e );
	    if( begin == This() ) begin = e;
	  }

    inline void InsertFirst(ThisType* &begin)
    {
      if( begin ) (begin)->prev = This();
      This()->prev = NULL;
      This()->next = begin;
      begin = This();
    }
    inline void InsertFirst(ThisType* &begin, ThisType* &end)
    {
      InsertFirst( begin );
      if( !end ) end = This();
    }

    inline void InsertLast(ThisType* &begin, ThisType* &end)
    {
      if( end ) end->next = This();
      This()->next = NULL;
      This()->prev = end;
      end = This();
      if( !begin ) begin = This();
    }

	  inline void Remove()
	  {
	    ThisType *n=next,*p=prev;
	    if(n) n->prev = prev;
	    if(p) p->next = next;
	  }
	  inline void Remove(ThisType* &begin)
	  {
	    Remove();
	    if( begin == This() )  begin = next;
	  }
	  inline void Remove(ThisType* &begin,ThisType* &end)
	  {
	    Remove( begin );
	    if( end == This() ) end = prev;
	  }
};



/*---------------------------------------*/
// a heap container
/*---------------------------------------*/


template<class T> struct HeapLessFunc {
	inline bool operator()( const T &a, const T &b ) { return a<b; };
};

template<class T> struct HeapSwapFunc {
	inline void operator()( T &a, T &b ) { T c(a); a=b; b=c; };
};


template<class T, class LessFunc = HeapLessFunc<T>, class SwapFunc = HeapSwapFunc<T> >
class Heap : DynArray<T>
{
public:
  typedef typename DynArray<T>::size_type size_type;
	LessFunc lessfunc;
	SwapFunc swapfunc;
	inline void Swap( T &a, T &b ) { swapfunc(a,b); }
	inline bool CmpLess( const T &a, const T &b ) { return lessfunc(a,b); }
	inline bool CmpGreater( const T &a, const T &b ) { return lessfunc(b,a); }

	Heap(int size, ConsTags::RESERVE_t) : DynArray<T>(size, ConsTags::RESERVE)
	{
	}

	Heap(int size, Cons::DONT_TYPE) : DynArray<T>(size)
	{
	}

  inline T& extremum()
  {
    return this->front();
  }

  inline const T& extremum() const
  {
    return this->front();
  }

  using DynArray<T>::operator[];
  using DynArray<T>::size;

  inline const T pop_extremum()
  {
    const T r(extremum());
    remove_extremum();
    return r;
  }

	void remove_extremum()
	{
		this->Swap(this->at(0),this->at(this->size()-1));
		this->pop_back();
		if(this->size()<=1) return;
		this->ShiftDown(0);
	}

	int insert(const T &x)
	{
		this->push_back(x);
		return this->ShiftUp(this->size()-1);
	}

	void remove(int k)
	{
		myAssert(k>=0 && k<this->size());
		swapfunc(this->at(k),this->at(this->size()-1));
		this->pop_back();
		if(this->size()<=1 || k==this->size()) return;
		if(this->CmpLess(this->at(k),this->at(this->PARENT(k))))
			this->ShiftUp(k);
		else if(
      (this->LEFT(k) <this->size() && this->CmpGreater(this->at(k),this->at(this->LEFT(k) ))) ||
			(this->RIGHT(k)<this->size() && this->CmpGreater(this->at(k),this->at(this->RIGHT(k))))
    )
      this->ShiftDown(k);
	}

	bool check()
	{
		printf("Checking Heap:");
		int visited=0;
		if(!this->CheckHeapRec(0,visited)) {
      myAssert(false);
			return false;
		}
		if(visited<this->size()) {
			printf("not all nodes visited");
      myAssert(false);
			return false;
		}
		printf("all clear");
		return true;
	}
private:
	inline int LEFT(int x) { return x*2+1; }
	inline int RIGHT(int x) { return x*2+2; }
	inline int PARENT(int x) { return (x-1)/2; }

	int ShiftUp(int k)
	{
		myAssert(k>=0 && k<this->size());
		int kk = this->PARENT(k);
		while(kk>=0 && this->CmpLess(this->at(k),this->at(kk))) {
			this->Swap(this->at(k),this->at(kk));
			k = kk;
			kk = this->PARENT(k);
		}
		return k;
	}

	void ShiftDown(int k)
	{
		myAssert(k>=0 && k<this->size());
		while(true)
    {
      int lk = this->LEFT(k);
		  int rk = this->RIGHT(k);
      // test the smaller child first, this is important
      if( lk<this->size() && rk<this->size() && this->CmpLess(this->at(rk),this->at(lk)) ) std::swap( lk,rk );

			if(lk<this->size() && this->CmpGreater(this->at(k),this->at(lk))) {
				this->Swap(this->at(k),this->at(lk));
				k = lk;
			}
			else if(rk<this->size() && this->CmpGreater(this->at(k),this->at(rk))) {
				this->Swap(this->at(k),this->at(rk));
				k = rk;
			}
			else break;
		}
	}

	bool CheckHeapRec(int k,int &visited)
	{
		++visited;
		int lk=this->LEFT(k);
		int rk=this->RIGHT(k);
		if(lk<this->size() && this->CmpGreater(this->at(k),this->at(lk))) {
			printf("heap violated at %i -> %i",k,lk);
			return false;
		}
		if(rk<this->size() && this->CmpGreater(this->at(k),this->at(rk))) {
			printf("heap violated at %i -> %i",k,rk);
			return false;
		}
		if(lk<this->size())
			if(!CheckHeapRec(lk,visited)) return false;
		if(rk<this->size())
			if(!this->CheckHeapRec(rk,visited)) return false;
		return true;
	}
};

#endif
