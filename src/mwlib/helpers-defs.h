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
#ifndef HELPERS_DEFS_H
#define HELPERS_DEFS_H

#include <stddef.h>
#include <string>
#include <ostream>

/*---------------------------------------*/
//  usefull typedefs
/*---------------------------------------*/
typedef unsigned int uint;
typedef unsigned char uchar;
typedef unsigned char ubyte;
typedef unsigned short uword;
typedef short word;
typedef unsigned long ulong;
typedef unsigned long uint64;
typedef long int64;
typedef unsigned int uint32;
typedef int int32;
static_assert(sizeof(int64)==8, "wrong type definition");
static_assert(sizeof(int32)==4, "wrong type definition");
// __GUNC__ is defined on linux gcc
// __LP64__ is defined on 64 bit linux gcc

using std::string;

typedef double PdeReal;

#define FOR_EACH(ITEM, A, CNT)\
  for(int __i=0,__t=1; __i<CNT; ++__i,__t=1)\
    for(ITEM=A[__i];__t;__t=!__t)

#define FOR_EACH2(ITEM1, A1, ITEM2, A2, CNT)\
  for(int __i=0,__t1=1,__t2=1; __i<CNT; ++__i,__t1=1,__t2=1)\
    for(ITEM2=A2[__i];__t2;__t2=!__t2)\
      for(ITEM1=A1[__i];__t1;__t1=!__t1)


#ifndef __GNUC__
  template <bool x> struct STATIC_ASSERTION_FAILURE;
  template <> struct STATIC_ASSERTION_FAILURE<true> { enum { value = 1 }; };
  #define STATIC_ASSERT( B ) \
    enum { _static_assert_test \
        = sizeof(STATIC_ASSERTION_FAILURE<(bool)(B)>) }
#else
  #define STATIC_ASSERT( B ) 
#endif

#ifdef DEBUG
#  define DEBUG_CODE(x) x
#  define IS_DEBUG true
#else
#  define DEBUG_CODE(x)
#  define IS_DEBUG false
#endif


#ifndef PRETTY_FUNCTION
#if defined _MSC_VER
#  define PRETTY_FUNCTION __FUNCTION__
#elif defined __GNUC__
#  define PRETTY_FUNCTION __func__
#else
#  define PRETTY_FUNCTION "no pretty function"
#endif
#endif

#define PING fprintf(stderr,PRETTY_FUNCTION)

#if DOTIMING
// my timer class
#include "timer.h"
// create time variable with some name
#define NAMED_TIMING_START(name) \
  my::Time my_timer_##name;
// output time between calls to output streams os
#define NAMED_TIMING_END_OS(os, name) \
  os << #name ": " << (my::Time() - my_timer_##name) << std::endl;
// timing shortcuts for function names
#define FUNC_TIMING_START \
  my::Time my_timer_##PRETTY_FUNCTION;
#define FUNC_TIMING_END_OS(os) \
  os << PRETTY_FUNCTION << ": " << (my::Time() - my_timer_##PRETTY_FUNCTION) << std::endl;
#else
#define NAMED_TIMING_START(name) (void)0;
#define NAMED_TIMING_END_OS(os, name) (void)0;
#define FUNC_TIMING_START (void)0;
#define FUNC_TIMING_END_OS(os) (void)0;
#endif


#define CALL_MEMBER_FN(object,ptrToMember)  ((object).*(ptrToMember))

template <bool flag, typename T, typename U>
struct Select
{
   typedef T Result;
};
template <typename T, typename U>
struct Select<false, T, U>
{
   typedef U Result;
};

template <int X>
struct Int2Type
{
  enum Value { value = X };
};

template <class T>
struct IdentityFunctor
{
  T operator()(const T &x) const { return x; }
};

template <class V>
struct ConstValueFunctor
{ 
  V v;
  ConstValueFunctor( const V &v ) : v(v) {}
  ConstValueFunctor() {}
  template<class A>
  const V operator()(const A &) const { return v; }
};

template<class T> struct RemoveConst {
  typedef T result;
};
template<class T> struct RemoveConst<const T> {
  typedef T result;
};

template<class T> struct IsConst {
  enum { result = false };
};
template<class T> struct IsConst<const T> {
  enum { result = true };
};


template<class T>
class SafeBool
{
  typedef void (SafeBool::*bool_type)() const;
  void this_type_does_not_support_comparisons() const {}

public:
  operator bool_type() const { return ((const T*)(this))->isTrue() ? &SafeBool::this_type_does_not_support_comparisons : NULL; }
};


namespace Cons
{
  enum DONT_TYPE { DONT };
  enum COPY_TYPE { COPY };
  enum DEEP_COPY_TYPE { DEEP_COPY };
  enum SHARED_TYPE { SHARED };
  enum FROM_GENERATOR_TYPE { FROM_GENERATOR };
  enum FROM_RANGE_TYPE { FROM_RANGE };
}

enum ConsMode
{
  AS_COPY,
  CLEAN,
};


namespace my
{

template<class T>
struct eqpair : public std::pair<T,T>
{
  eqpair(const T &a, const T &b) : std::pair<T,T>(a,b) {}
  eqpair() : std::pair<T,T>() {}

  using std::pair<T,T>::first;
  using std::pair<T,T>::second;

  inline T operator[](int i) const { return *((&first)+i); }
  inline T& operator[](int i) { return *((&first)+i); }
};

template<class T>
eqpair<T> make_eqpair(const T &a, const T &b) { return eqpair<T>(a,b); }

}

namespace std
{

template<class A, class B>
std::ostream& operator<<(std::ostream &os, const std::pair<A, B> &p)
{
  os << "(" << p.first << ", " << p.second << ")";
  return os;
}

}

template<class Container>
void clear_and_free_memory(Container &v) { Container().swap(v); }


template<class Container>
void print_container(std::ostream &os, const Container &v)
{
  typedef typename Container::const_iterator Iter;
  os << "[";
  Iter it = v.begin();
  if (it != v.end())
  {
    while(true)
    {
      os << (*it);
      ++it;
      if (it == v.end()) break;
      os << ", ";
    }
  }
  os << "]";
}

template<class T>
inline T* get_ptr(T *p) { return p; }

template<class T>
inline const T* get_ptr(const T* p) { return p; }


  /*
   * this struct is used to return the edges and nodes neighboring a
   * certain node which is queried for its neighbors
   */
template< class NodeT, class EdgeT >
struct NodeNeighbor
{
  NodeNeighbor( NodeT n, EdgeT e ) : node(n),edge(e) {}
  NodeT node;
  EdgeT edge;
};

template<class N, class E>
inline const NodeNeighbor<N, E> make_node_neighbor(const N &n, const E &e) { return NodeNeighbor<N,E>(n, e); }


template<class T> inline std::size_t estimateMemoryUsage(const T &) { return sizeof(T); }
// inline std::size_t estimateMemoryUsage(const double &) { return sizeof(double); }
// inline std::size_t estimateMemoryUsage(const long &) { return sizeof(long); }
// inline std::size_t estimateMemoryUsage(const float &) { return sizeof(float); }
// inline std::size_t estimateMemoryUsage(const uint &) { return sizeof(uint); }
// inline std::size_t estimateMemoryUsage(const char &) { return sizeof(char); }



#endif
