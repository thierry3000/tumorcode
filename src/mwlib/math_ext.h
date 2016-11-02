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
#ifndef MY_MATH_H
#define MY_MATH_H

#include <cmath>
#include <limits>
#include <ostream>
#include <assert.h>

namespace my
{
  
namespace mconst
{
#define DECL_MATH_CONST(name, type, val) inline type name() { return val; }
DECL_MATH_CONST(pi,double,M_PI)
DECL_MATH_CONST(pi05,double,M_PI*0.5)
DECL_MATH_CONST(pi2,double,M_PI*2)
DECL_MATH_CONST(fpi,float,(float)M_PI)
DECL_MATH_CONST(fpi05,float,(float)M_PI*0.5)
DECL_MATH_CONST(fpi2,float,(float)M_PI*2)
DECL_MATH_CONST(fln2,float,(float)M_LN2);
DECL_MATH_CONST(ln2,double,M_LN2);
#undef DECL_MATH_CONST

const int primes[] = {
  2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997
};

}

template<class T>
T clamp(const T &x, const T &a, const T &b) 
{
    return std::max(a, std::min(b, x));
}

template<class T>
T cut(const T &x, const T &a, const T &b) { return clamp(x, a, b); }

template<class T>
inline int iceil(T x)
{
  T intpart;
  T fracpart = std::modf(x, &intpart);
  if (intpart < x)
    intpart += T(1.0);
  return int(intpart);
}

template<class T>
inline int ifloor(T x)
{
  T intpart;
  T fracpart = std::modf(x, &intpart);
  if (intpart > x)
    intpart -= T(1.0);
  return int(intpart);
}


inline int ipow(int a, int b)
{
  int r = 1;
  for (int i=0; i<b; ++i)
    r *= a;
  return r;
}


inline float round( float x, float decimal_place = 1. ) { return float(int(x/decimal_place+0.5f)*decimal_place); }

template<class T>
inline T sqr(const T &x) { return x*x; }

template<class T>
inline T cubed(const T &x) { return x*x*x; }

template< class T, class U >
inline T lerp( U x, const T &a, const T &b ) {  return T((U(1)-x)*a + x*b); }

template< class T >
inline T normalize( const T x, const T min, const T max )
{
  return (max != min) ? (x-min)/(max-min) : x;
}

template< class T >
inline T normalize_cut( const T x, const T min, const T max ) {
  return my::clamp( normalize(x,min,max), T(0), T(1) );
}

template< class T >
inline T normalize_log( T x, T min, T max, T eps ) { return normalize( std::log10(x+eps), std::log10(min+eps), std::log10(max+eps) ); }

template<class T>
inline void sincos( T x, T &s, T &c )
{
  s = std::sin(x);
  c = std::cos(x);
}

template<class T>
inline T smooth_heaviside(T x, double width)
{
  assert(width > 0.);
  T q = (x*10)/(M_PI*width);
  if (q < -20.) return T(0);
  if (q > 20.) return T(1);
  return (tanh(q)+1)*0.5;
}


/*
 * finite width smooth delta function approximation, see [peskin_2002]
 * has bounded support on  (-2*radius, 2*radius)
 */
template<class T, class S>
inline T smooth_delta_cos(T x, S radius)
{
  x /= radius;
  if (x < S(2.) && x > S(-2.))
    return S(0.25)*(S(1.)+cos(my::mconst::pi05()*x)) / radius;
  else
    return T();
}

/*
 * integral over smooth_delta_cos
 */
template<class T, class S>
inline T smooth_heaviside_sin(T x, S radius)
{
  x /= radius;
  if (x <= S(-2))
    return T(0);
  else if(x >= S(2))
    return T(1);
  else
    return S(0.25)*(S(2) + x + S(1./my::mconst::pi05()) * sin(my::mconst::pi05()*x));
}


template<class T>
T smooth_lowclip( T x, const T w )
{
  assert(w > 0.);
  const T w2(0.5*w);
  if(x<-w2) return T(0);
  if(x>w2) return x;
  return (0.5/w)*sqr<T>(x+w2);
}

template<class T>
T smooth_max( T x, T y, const T w )
{
  assert(w > 0.);
  return x + smooth_lowclip<T>(y-x,w);
}

template<class T>
T sign( T x )
{
  return x >= T(0) ? T(1) : T(-1);
}

inline int mod( int a, int n ) { const int m=a%n; return m<0 ? n+m : m; }

inline float mod( float a, float n ) { return a>=0 ? fmodf(a,n) : n-fmodf(-a,n); }

using std::min;

template<class T>
T min(const T &a, const T &b, const T &c)
{
  return std::min(std::min(a,b), c);
}

template<class T>
T min(const T &a, const T &b, const T &c, const T &d)
{
  return std::min(my::min(a,b,c), d);
}



using std::max;

template<class T>
T max(const T &a, const T &b, const T &c)
{
  return std::max(std::max(a,b), c);
}

template<class T>
T max(const T &a, const T &b, const T &c, const T &d)
{
  return std::max(my::max(a,b,c), d);
}


namespace averaging
{
  template< class T, class S >
  inline T avg( const T w, const S cnt ) { return w/cnt; }

  template< class T, class S >
  inline T var( const T w, const T w2, const S cnt ) {
    return (w2-w*w/cnt)/cnt;
  }
//   template< class T, class S >
//   inline T var( const T avg, const T w2, const S cnt ) {
//     return (w2/cnt-avg*avg);
//   }

  template< class T, class S >
  inline T stddev( const T w, const T w2, const S cnt ) {
    const T v = var<T,S>(w,w2,cnt);
    assert(v > -(cnt*std::numeric_limits<T>::epsilon()));
    return v<T(0) ? T(0) : std::sqrt(v);
  }

//   template< class T, class S >
//   inline T stddev( const T avg, const T w2, const S cnt ) {
//     return std::sqrt(var(avg,w2,cnt));
//   }
}



template< class T >
struct MinMax
{
  T min,max;
  MinMax() { clear(); }
  MinMax( const T &mm ) : min(mm),max(mm) {}
  MinMax( const T& min, const T& max ) : min(min),max(max) {}
  MinMax& add( const T &x ) { if( x<min ) min=x; if(x>max) max=x; return *this; }
  MinMax& add( const MinMax<T> &x ) { add(x.min); add(x.max); return *this; }
  void clear() { min=std::numeric_limits<T>::max(); max=-std::numeric_limits<T>::max(); }
  inline friend bool operator==( const MinMax<T> &a, const MinMax<T> &b ) { return a.min==b.min && a.max==b.max; }
  inline friend const MinMax<T> operator&( const MinMax& a, const MinMax& b )
  {
    if(a.max<b.min || b.max<a.min) return MinMax<T>();
    else return MinMax(std::max(a.min,b.min),std::min(a.max,b.max));
  }
  inline friend const MinMax<T> operator|( const MinMax& a, const MinMax& b )
  {
    return MinMax(std::min(a.min,b.min),std::max(a.max,b.max));
  }
};


template<class T>
class Averaged
{
  T w,w2,min,max;
  uint n;
public:
  explicit Averaged() { Clear(); }
  explicit Averaged( const T &x ) { Clear(); Add(x); }
  Averaged( const Averaged &a ) { Clear(); Add(a.Avg()); }

  Averaged& operator=( const Averaged &a )
  {
    w2=a.w2;
    w=a.w;
    n=a.n;
    min=a.min;
    max=a.max;
  }
  
  void Clear()
  {
    w2=w=T();
    n=0;
    min=max=T();
  }

  void Add( const T &x )
  {
    w  += x;
    w2 += x*x;
    if( n<=0 ) { min=max=x; }
    else {
      if(x<min) min=x;
      if(x>max) max=x;
    }
    ++n;
  }

  void Add( const Averaged<T> &x )
  {
    w += x.w;
    w2 += x.w2;
    n += x.n;
    if( x.min<min ) min=x.min;
    if( x.max>max ) max=x.max;
  }

  const T Avg() const { return averaging::avg(w,n); }
  const T RmsE() const { return averaging::stddev(w,w2,n); }
  const T Min() const { return min; }
  const T Max() const { return max; }
  const T Sum() const { return w; }
  uint Count() const { return n; }
  const T MaxAbs() const { return std::max(std::abs(min), std::abs(max)); }

  MinMax<T> getMinMax() const { return MinMax<T>(min, max); }
};

template<class T>
inline std::ostream& operator<<(std::ostream& os,const Averaged<T> &a)
{
  os <<a.Avg()<<"+/-"<<a.RmsE()<<"\\in["<<a.Min()<<","<<a.Max()<<"]";
  return os;
}



template<class T>
struct Bitfield
{
  T bits;
  inline Bitfield& AddBits( T b ) { bits |= b; return *this; }
  inline Bitfield& AddBits( T b, const bool bOn ) { bits = bOn ? (bits|b)  : bits&(~b);  return *this; }
  inline Bitfield& AddBits( T b, T mask ) { bits = (bits&(~mask)) | (b&mask);  return *this;  }
  inline const T GetBits( T mask ) const { return bits&mask; }
  inline bool AllBitsSet( T mask ) const { return (bits&mask)==mask; }
  inline Bitfield& DelBits( T b ) { bits &= (~b);  return *this;  }
  inline Bitfield& ToggleBits( T b ) { bits = (bits^b);  return *this;  }

  inline Bitfield& operator |= (const T x) { AddBits(x); return *this; }

  operator const T() const { return bits; }
  const T operator=( const T x ) { bits=x; return bits; }
  Bitfield( T x ) : bits(x) {}
  Bitfield() : bits(T()) {}
  Bitfield( const Bitfield<T> &x ) : bits(x.bits) {}
};



template<class T>
inline std::ostream& operator<<(std::ostream& os, const MinMax<T> &mm)
{
  os << "[" << mm.min << ", " << mm.max << "]";
  return os;
}

}

namespace my
{

// old stuff
template< class T >
inline T Cut( const T x, const T a, const T b ) { return (x<a) ? a : ((x>b) ? b : x); }

template< class T >
inline T Normalize( const T x, const T min, const T max ) { return (x-min)/(max-min); }

template< class T >
inline T NormalizeCut( const T x, const T min, const T max ) {
  return Cut( Normalize(x,min,max), T(0), T(1) );
}

template< class T, class U >
inline T Lerp( U x, const T &a, const T &b ) {  return T((1.0-x)*a + x*b); }

template< class T, class U >
inline T Lerp2D( U x, U y, const T &a00, const T &a01, const T &a10, const T &a11 )
{
/*
  y 10--11
    |    |
    00--01
         x
*/
  return Lerp<T,U>(y,Lerp<T,U>(x,a00,a01),Lerp<T,U>(x,a10,a11));
}

template<class T>
inline T Min(const T &a, const T &b) { return std::min(a,b); }

template<class T>
inline T Max(const T &a, const T &b) { return std::max(a,b); }


//this class is meant to protect from accientially doing a = b instead of a = min(a,b)
template<class T>
class Smallest
{
  T x;
  Smallest& operator=(T); // don't want that
public:
  Smallest() : x() {}
  explicit Smallest(T x_) : x(x_) {}
  operator T() const { return x; }  // can be used where double is required
  template<class U>
  Smallest& min(const U &y) { x = std::min<T>(x, y); return *this; }
};

template<class T>
inline Smallest<T> MakeSmallestMax()
{
  return Smallest<T>(std::numeric_limits<T>::max());
}

template<class T>
static std::ostream& operator<<(std::ostream &os, const Smallest<T> &x)
{
  if (x >= MakeSmallestMax<T>())
    os << "max";
  else
    os << static_cast<T>(x);
  return os;
}


inline double compare_approx_less(double a, double b)
{
  double diff = a-b;
  double sum  = std::max(std::abs(a), std::abs(b));
  // return diff < 0.
  return diff <= 16. * std::numeric_limits<double>::epsilon() * sum;
}

inline double compare_approx_equal(double a, double b)
{
  double diff = std::abs(a-b);
  double sum  = std::max(std::abs(a), std::abs(b));
  return diff <= 16. * std::numeric_limits<double>::epsilon() * sum;
}


}

#endif

class T;

class T;
