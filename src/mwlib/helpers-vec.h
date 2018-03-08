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
#ifndef HELPERS_VEC_H
#define HELPERS_VEC_H

#include "helpers-defs.h"
#include "myAssert.h"

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>
#include <boost/format.hpp>
#include <boost/functional/hash.hpp>

/* teach boost how to serialize eigen 
 * http://stackoverflow.com/questions/12851126/serializing-eigens-matrix-using-boost-serialization
 */
// namespace boost
// {
//     template<class Archive, typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
//     inline void serialize(
//         Archive & ar, 
//         Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> & t, 
//         const unsigned int file_version
//     ) 
//     {
//         ar & boost::serialization::make_array(t.data(), t.size());
//     }
// }
template<class T >
class DynArray;
template<class T, int d>
class Vec : public Eigen::Matrix<T, d, 1>
{
  public:
    typedef Eigen::Matrix<T, d, 1> Base;
    typedef T value_type;
    // constructors
    Vec() : Base(Base::Constant(0.)) {}
    explicit Vec(const T &v) : Base() { this->setConstant(d, v); }
    Vec(const T &a, const T &b) : Base(a, b) {}
    Vec(const T &a, const T &b, const T &c) : Base(a, b, c) {}
    template<typename OtherDerived >
    Vec (const Eigen::MatrixBase< OtherDerived > &other) : Base(other) {}

    template<class FUNC>
    static inline Vec<T, d> mapIndex(FUNC f)
    {
      Vec<T, d> ret;
      for (int i=0; i<d; ++i)
        ret[i] = f(i);
      return ret;
    }
// private:
//   friend class boost::serialization::access;
//   template <typename Archive>
//   void serialize(Archive &ar, const unsigned int version)
//   {
//     //ar & ;
//     ar & boost::serialization::base_object<Eigen::Matrix<T, d, 1>>(*this);
//   }
};


#if 0
typedef Eigen::Matrix<int, dim, 1> Intd;
typedef Eigen::Matrix<double, dim, 1> Doubled;
typedef Eigen::Matrix<int, 3, 1> Int3;
typedef Eigen::Matrix<int, 2, 1> Int2;
typedef Eigen::Matrix<double, 3, 1> Double3;
typedef Eigen::Matrix<double, 2, 1> Double2;
typedef Eigen::Matrix<float, 3, 1> Float3;
typedef Eigen::Matrix<float, 2, 1> Float2;
#define VECDARG class T,int d
#define VECD Eigen::Matrix<T, d, 1>
#else
//typedef Vec<int, dim> Intd;
//typedef Vec<double, dim> Doubled;
typedef Vec<int, 6> Int6;
typedef Vec<int, 3> Int3;
typedef Vec<int, 2> Int2;
typedef Vec<bool, 3> Bool3;
typedef Vec<bool, 2> Bool2;
//typedef Vec<double, 6> Double6;
typedef Vec<double, 3> Double3;
typedef Vec<double, 2> Double2;
//typedef Vec<float, 6> Float6;
typedef Vec<float, 3> Float3;
typedef Vec<float, 2> Float2;
typedef Vec<long, 2> Long2;
typedef Vec<long, 2> Long3;
#define VECDARG class T,int d
#define VECD Vec<T, d>
#endif

namespace std
{
  template<class T, int dim>
  struct hash<Vec<T,dim> >
  {
    std::size_t operator()(const Vec<T, dim> &key) const
    {
      std::size_t seed = 0;
      for (int i=0; i<dim; ++i)
        boost::hash_combine(seed, boost::hash_value(key[i]));
      return seed;
    }
  };
}

//namespace Eigen
//{
template<VECDARG>
inline std::ostream& operator<<(std::ostream &os, const VECD &v)
{
  os << "<";
  for (int i=0; i<d; ++i)
  {
    os << v[i];
    if (i < d -1) os << ",";
  }
  os << ">";
  return os;
}

template<VECDARG>
inline std::istream& operator>>(std::istream &is, VECD &v)
{
  char c;
  is >> c;
  if (c != '<') { is.setstate(std::ios_base::failbit); return is; }
  for (int i=0; i<d; ++i)
  {
    is >> v[i];
    if (i < d -1) is >> c;
    if (c != ',') { is.setstate(std::ios_base::failbit); return is; }
  }
  is >> c;
  if (c != '>') { is.setstate(std::ios_base::failbit); return is; }
  return is;
}

template <class U, class T, int dim>
inline Vec<U, dim> Cast(const Vec<T, dim> &v)
{
  return v.template cast<U>().eval();
}



//}

#if 0
template<class T, int d>
inline VECD make_vec(T a)
{
  VECD q;
  for(int i=0; i<d; ++i) q[i] = a;
  return q;
}

template<class T, int d>
inline VECD make_vec(T a, T b)
{
  VECD q;
  q[0] = a;
  for(int i=1; i<d; ++i) q[i] = b;
  return q;
}

template<class T, int d>
inline VECD make_vec(T a, T b, T c)
{
  VECD q;
  q[0] = a;
  q[1] = b;
  for(int i=2; i<d; ++i) q[i] = c;
  return q;
}
#endif

template<class T>
inline T cross( const Vec<T,2>& u, const Vec<T,2> &v )
{
  return u.x()*v.y()-u.y()*v.x();
}

template<class T>
inline Vec<T,3> cross( const Vec<T,3>& u, const Vec<T,3> &v )
{
  return u.cross(v);
}

// template<class Derived>
// inline typename Eigen::internal::traits<Derived>::Scalar product_of_components(const Eigen::MatrixBase<Derived> &a) {
//   return a.prod();
// }
// 
// template<class Derived>
// inline typename Eigen::internal::traits<Derived>::Scalar sum_of_components(const Eigen::MatrixBase<Derived> &a) {
//   return a.sum();
// }


template<class Derived>
inline typename Eigen::internal::plain_matrix_type<Derived>::type safe_normalized(const Eigen::MatrixBase<Derived> &a) {
  typedef typename Eigen::internal::traits<Derived>::Scalar S;
  S n = a.norm();
  return (n > std::numeric_limits<S>::epsilon()) ? (a / n).eval() : (0 * a).eval();
}

#define DECLARE_EIGEN_COMPARE_FUNCTION_REDUCTION(op, name) \
  template<class Derived, class OtherDerived>\
  inline bool name( const Eigen::MatrixBase<Derived> &a, const Eigen::MatrixBase<OtherDerived> &b ) {\
    return (a.array() op b.array()).all();\
  }\
  template<class Derived>\
  inline bool name( const Eigen::MatrixBase<Derived> &a, const typename Derived::Scalar &b) {\
    return (a.array() op b).all();\
  }

DECLARE_EIGEN_COMPARE_FUNCTION_REDUCTION(==,allEq)
DECLARE_EIGEN_COMPARE_FUNCTION_REDUCTION(<,allL)
DECLARE_EIGEN_COMPARE_FUNCTION_REDUCTION(>,allG)
DECLARE_EIGEN_COMPARE_FUNCTION_REDUCTION(<=,allLe)
DECLARE_EIGEN_COMPARE_FUNCTION_REDUCTION(>=,allGe)
  
#undef DECLARE_EIGEN_COMPARE_FUNCTION_REDUCTION

// Note that the return type is an actual matrix object!
// It is unfortunately not allowed to return an expression object! Doing so apparently results in dangling references to temporaries ...
#define DECLARE_EIGEN_COMPARE_FUNCTION(op, name) \
  template<class DerivedA, class DerivedB> \
  inline Eigen::Matrix<bool, DerivedA::RowsAtCompileTime, DerivedA::ColsAtCompileTime> name(const Eigen::MatrixBase<DerivedA> &a, const Eigen::MatrixBase<DerivedB> &b) \
  { \
    return (a.array() op b.array()).matrix(); \
  }\
  template<class DerivedA> \
  inline Eigen::Matrix<bool, DerivedA::RowsAtCompileTime, DerivedA::ColsAtCompileTime> name(const Eigen::MatrixBase<DerivedA> &a, const typename DerivedA::Scalar &b) \
  { \
    return (a.array() op b).matrix(); \
  }
DECLARE_EIGEN_COMPARE_FUNCTION(<, cwiseL)
DECLARE_EIGEN_COMPARE_FUNCTION(>, cwiseG)
DECLARE_EIGEN_COMPARE_FUNCTION(<=, cwiseLe)
DECLARE_EIGEN_COMPARE_FUNCTION(>=, cwiseGe)

#undef DECLARE_EIGEN_COMPARE_FUNCTION

#if 0
// other operators
template<class Derived, class OtherDerived>
inline typename Eigen::internal::traits<Derived>::Scalar DiffLenSqr( const Eigen::MatrixBase<Derived>& a, const Eigen::MatrixBase<OtherDerived> &b ) {
  return (a - b).squaredNorm();
}
#endif

template<class Derived>
inline double norm( const Eigen::MatrixBase<Derived>& v )
{
  return v.norm();
}

template<class Derived>
inline typename Eigen::internal::traits<Derived>::Scalar squaredNorm( const Eigen::MatrixBase<Derived>& v )
{
  return v.squaredNorm();
}


template<class Derived>
inline typename Eigen::internal::traits<Derived>::Scalar maxCoeff(const Eigen::MatrixBase<Derived>& a) {
  return a.maxCoeff();
}

template<class Derived>
inline typename Eigen::internal::traits<Derived>::Scalar minCoeff(const Eigen::MatrixBase<Derived>& a) {
  return a.minCoeff();
}


template<class Derived, class DerivedB, class DerivedC>
inline typename Eigen::internal::plain_matrix_type<Derived>::type clamp( const Eigen::MatrixBase<Derived>& v, const Eigen::MatrixBase<DerivedB>& a, const Eigen::MatrixBase<DerivedC>& b )
{
  return typename Eigen::internal::plain_matrix_type<Derived>::type(v.cwiseMax(a).cwiseMin(b));
}

template<class T>
inline Eigen::Matrix<T, 2, 1> Rot90( const Eigen::Matrix<T, 2, 1> &v )
{
  Eigen::Matrix<T, 2, 1> r;
  r[1] = v[0];
  r[0] = -v[1];
  return r;
}


template<class Derived>
inline int MajorAxis(const Eigen::MatrixBase<Derived> &v ) {
  typedef typename Eigen::internal::traits<Derived>::Scalar Scalar;
  Scalar m = Scalar();
  int a = -1;
  for (int i=0; i<v.size(); ++i)
    if (std::abs(v[i]) > m)
    {
      m = std::abs(v[i]);
      a = i;
    }
  return a;
}

template<class Derived>
inline typename Eigen::internal::plain_matrix_type<Derived>::type add_to_axis(const Eigen::MatrixBase<Derived>& p, int axis, int offset)
{
  typename Eigen::internal::plain_matrix_type<Derived>::type q(p); q[axis]+=offset;
  return q;
}


template<class Derived>
inline typename Eigen::internal::plain_matrix_type<Derived>::type replace(const Eigen::MatrixBase<Derived>& p, int axis, typename Eigen::internal::traits<Derived>::Scalar value)
{
  typename Eigen::internal::plain_matrix_type<Derived>::type q(p); q[axis]=value;
  return q;
}



#undef VECD
#undef VECDARG




template<class T, int d>
struct BBoxd
{
  typedef ::Vec<T, d> Vec;
  enum {
    MIN = 0,
    MAX = 1,
    AXISX = 0,
    AXISY = 1,
    AXISZ = 2
  };

  Vec min,max;

  BBoxd(T xmin,T ymin, T zmin,T xmax,T ymax, T zmax) :
    min(Vec(xmin,ymin,zmin)),max(Vec(xmax,ymax,zmax)) {}
  BBoxd(T xmin, T ymin, T xmax, T ymax) :
    min(Vec(xmin, ymin)), max(Vec(xmax, ymax)) {}
    
  BBoxd(const Vec &min, const Vec &max) : min(min), max(max)
  {
    myAssert(allLe(min, max));
  }

  BBoxd() :
    min(Vec::Constant(std::numeric_limits<T>::max())), max(Vec::Constant(-std::numeric_limits<T>::max())) {}

  explicit BBoxd(Cons::DONT_TYPE) : min(Cons::DONT), max(Cons::DONT) {}

  int size() const { return d; }

  friend bool operator==(const BBoxd &a, const BBoxd &b)
  {
    return a.min==b.min && a.max==b.max;
  }

  friend bool operator!=(const BBoxd &a, const BBoxd &b)
  {
    return !(a == b);
  }

  inline bool IsEmpty() const
  {
    return !allGe(max,min);
  }

  inline BBoxd& Clear(){
    *this = BBoxd();
    return *this;
  }

  inline BBoxd& Add( const BBoxd<T,d>::Vec &p ) {
    min = min.cwiseMin(p);
    max = max.cwiseMax(p);
    return *this;
  }

  inline BBoxd& Add( const BBoxd &bbox ) {
    if(bbox.IsEmpty()) return *this;
    Add( bbox.min );
    Add( bbox.max );
    return *this;
  }

  inline BBoxd& Extend( const Vec &q ) {
    max += q;
    min -= q;
    return *this;
  }

  inline BBoxd& Extend( const T q ) {
    return Extend(Vec::Constant(q));
  }

  inline BBoxd& ExtendFactor( const T q ) {
    Vec mp = 0.5f*(max+min);
    Vec ex = 0.5f*(max-min);
    min = mp-q*ex;
    max = mp+q*ex;
    return *this;
  }

  inline BBoxd& Intersection( const BBoxd &clip ) {
    max = max.cwiseMin(clip.max);
    min = min.cwiseMax(clip.min);
    return *this;
  }

  inline BBoxd& Set( int axis, const T q ) {
    min[axis] = max[axis] = q;
    return *this;
  }
  
  inline BBoxd& Set(int side, int axis, const T q) {
    if(side<=0) min[axis]=q;
    else max[axis]=q;
    return *this;
  };

  inline const Vec operator[](int side) const {
    return side<=0 ? min : max;
  }
   
  inline Vec& operator[](int side) {
    return side<=0 ? min : max;
  }

  inline BBoxd& Move(const Vec &p) {
    max += p;
    min += p;
    return *this;
  }

  inline BBoxd& Move(int axis, const T q)
  {
    min[axis] += q;
    max[axis] += q;
    return *this;
  }

  inline BBoxd& Move(int side, int axis, const T q)
  {
    (*this)[side][axis] += q;
    return *this;
  }

  inline BBoxd& MoveMax(const T q)
  {
    max.array() += q;
    return *this;
  }

  inline BBoxd& Multiply(const Vec &p) {
    max.array() *= p.array();
    min.array() *= p.array();
    return *this;
  }

  inline BBoxd& Divide(const Vec &p) {
    max.array() /= p.array();
    min.array() /= p.array();
    return *this;
  }

  inline BBoxd& Multiply(const T p) {
    return Multiply(Vec(p));
  };

  inline BBoxd& Divide(const T p) {
    return Divide(Vec(p));
  }

  inline BBoxd& operator+=(const Vec &p) { return Move(p); }
  inline BBoxd& operator*=(const Vec &p) { return Multiply(p); }
  inline BBoxd& operator*=(const T p) { return Multiply(p); }
  inline BBoxd& operator/=(const Vec &p) { return Divide(p); }
  inline BBoxd& operator/=(const T p) { return Divide(p); }

  inline bool Overlaps( const Vec &p ) const {
    return allGe(p, min) && allLe(p, max);
  }
  inline bool Overlaps(const BBoxd &b) const {
    // check for overlap, including the boundary points
    // comparison operators return true if comparison statement is true for all components
    return !(!allLe(b.min, max) || !allGe(b.max, min));
  }

  inline bool Includes(const BBoxd &b) const {
    return allGe(b.min, min) && allLe(b.max, max);
  }

  friend const BBoxd operator+(const Vec &p, const BBoxd &b) { return BBoxd(b).Move(p); }
  friend const BBoxd operator*( const Vec &p, const BBoxd &b) { return BBoxd(b).Multiply(p); }
  friend const BBoxd operator/(const BBoxd &b, const Vec &p) { return BBoxd(b).Divide(p); }
  friend const BBoxd operator*(const T &p, const BBoxd &b) { return BBoxd(b).Multiply(p); }
  friend const BBoxd operator/(const BBoxd &b, const T &p) { return BBoxd(b).Divide(p); }
  friend const BBoxd Extend(const BBoxd &b, const Vec &p) { return BBoxd(b).Extend(p); }
  friend const BBoxd Extend(const BBoxd &b, const T &p)  { return BBoxd(b).Extend(p); }
  friend const BBoxd Intersection(const BBoxd &b, const BBoxd &c) { return BBoxd(b).Intersection(c); }
  friend const BBoxd Set(const BBoxd &b, int side, int axis, const T &p)  { return BBoxd(b).Set(side, axis, p); }
  friend const BBoxd Set(const BBoxd &b, int axis, const T &p)  { return BBoxd(b).Set(axis, p); }
  friend const BBoxd Move(const BBoxd &b, const Vec &p) { return BBoxd(b).Move(p); }
  friend const BBoxd Move(const BBoxd &b, int side, int axis, const T &p)  { return BBoxd(b).Move(side, axis, p); }
  friend const BBoxd Move(const BBoxd &b, int axis, const T &p)  { return BBoxd(b).Move(axis, p); }
};


template<int d>
inline Vec<int, d> Size(const BBoxd<int, d> &bb)
{
  // the min/max values always include the indices of the boundary points
  // so the number of points is one more than the difference
  if(bb.IsEmpty()) return Vec<int, d>(0);
  else return bb.max - bb.min + Vec<int, d>(1);
}


template<int d>
inline int Volume(const BBoxd<int, d> &bb)
{
  // the total number of points
  const Vec<int, d> s = Size(bb);
  return s.prod();
}



template<class T>
inline void write_gnuplot(std::ostream &os, const BBoxd<T,3> &bb, double w)
{
  boost::format fmt("%f %f %f %f");
  for (int x=0; x<2; ++x)
    for (int y=0; y<2; ++y)
    {
      for (int z=0; z<2; ++z)
        os << fmt % bb[x][0] % bb[y][1] % bb[z][2] % w << std::endl;
      os << std::endl << std::endl;
      for (int z=0; z<2; ++z)
        os << fmt % bb[z][0] % bb[y][1] % bb[x][2] % w << std::endl;
      os << std::endl << std::endl;
      for (int z=0; z<2; ++z)
        os << fmt % bb[y][0] % bb[z][1] % bb[x][2] % w << std::endl;
      os << std::endl << std::endl;
    }
}


template<class T>
inline void write_gnuplot(std::ostream &os, const BBoxd<T,2> &bb, double w)
{
  boost::format fmt("%f %f %f");
  for (int x=0; x<2; ++x)
  {
    for (int z=0; z<2; ++z)
      os << fmt % bb[x][0] % bb[z][1] % w << std::endl;
    os << std::endl << std::endl;
    for (int z=0; z<2; ++z)
      os << fmt % bb[z][0] % bb[x][1] % w << std::endl;
    os << std::endl << std::endl;
  }
}



typedef BBoxd<int, 3> BBox3;
typedef BBoxd<float, 3> FloatBBox3;
typedef BBoxd<int, 2> BBox2;
typedef BBoxd<float, 3> FloatBBox3;

inline const Float3 Size(const FloatBBox3 &bb) { return bb.max-bb.min; }
inline const Float3 MP(const FloatBBox3 &bb) { return 0.5f*(bb.max+bb.min); }
inline BBox3 MkBBox3( const Int3 &size, const Int3 &offset=Vec<int,3>(0) ) { return BBox3().Add(offset).Add(size+offset-Vec<int,3>(1)); }

inline BBox3 ExtendForDim(const BBox3 &bb, int dim, int e)
{
  Int3 ee;
  for (int i=0; i<dim; ++i)
    ee[i] = e;
  return BBox3(bb).Extend(ee);
}

#define FOR_BBOX2(p,bb) for(Int2 p=bb.min; p[1]<=bb.max[1]; ++p[1]) for(p[0]=bb.min[0]; p[0]<=bb.max[0]; ++p[0])
#define FOR_BBOX3(p,bb)\
  for(std::pair<const BBox3, bool> __forbb(bb,true); __forbb.second; __forbb.second=!__forbb.second)\
    for(Int3 p=__forbb.first.min; p[2]<=__forbb.first.max[2]; ++p[2]) for(p[1]=__forbb.first.min[1]; p[1]<=__forbb.first.max[1]; ++p[1]) for(p[0]=__forbb.first.min[0]; p[0]<=__forbb.first.max[0]; ++p[0])
#define FOR_REG3(T,x,x0,x1,y,y0,y1,z,z0,z1) for(T z=z0; z<=z1; ++z) for(T y=y0; y<=y1; ++y) for(T x=x0; x<=x1; ++x)
#define FOR_REG3V(p,x0,x1,y0,y1,z0,z1) for(Int3 p = Vec<int,3>(x0,y0,z0); p[2]<=z1; ++p[2]) for(p[1]=y0; p[1]<=y1; ++p[1]) for(p[0]=x0; p[0]<=x1; ++p[0])
#define FOR_REG2V(p,x0,x1,y0,y1) for(Int2 p = Vec<int,2>(x0,y0); p[1]<=y1; ++p[1]) for(p[0]=x0; p[0]<=x1; ++p[0])
#define FOR_REG3V2(p,s) for(Int3 p = Vec<int,3>(0),_s(s); p[2]<_s[2]; ++p[2]) for(p[1]=0; p[1]<_s[1]; ++p[1]) for(p[0]=0; p[0]<_s[0]; ++p[0])
#define FOR_REG3V3(p, s, incx, incy, incz) for(Int3 p(0),_s(s); p[2]<_s[2]; p[2]+=incz) for(p[1]=0; p[1]<_s[1]; p[1]+=incy) for(p[0]=0; p[0]<_s[0]; p[0]+=incx)
#define FOR_REG2V2(p,s) for(Int2 p = Vec<int,2>(0),_s(s); p[1]<_s[1]; ++p[1]) for(p[0]=0; p[0]<_s[0]; ++p[0])
#define FOR_BBOX3_BORDER(T,p,bb)\
  for(T inc,p=bb.min  ; p[2]<=bb.max[2]; ++p[2])\
  for(       p[1]=bb.min[1]; p[1]<=bb.max[1]; ++p[1])\
  for(inc[0]=(p[1]==bb.min[1] || p[1]==bb.max[1] || p[2]==bb.min[2] || p[2]==bb.max[2])?1:(bb.max[0]-bb.min[0]),p[0]=bb.min[0]; p[0]<=bb.max[0]; p[0]+=inc[0])


template<class T, int d>
inline std::ostream& operator<<(std::ostream& out,const BBoxd<T, d> &r )
{
  out << "[";
  for (int i=0; i<r.size(); ++i) {
    out << r.min[i] << "," << r.max[i];
    if (i < r.size()-1)
      out << "]x[";
  }
  out << "]";
  return out;
}

template<class T> class DynArray;
enum {
  BBOXA = 1,
  BBOXB = 2,
  INCLUDE_UNION = 4
};
int SplitBBoxByBox(const BBox3 &bb, const BBox3 &splitbox, DynArray<BBox3> &res, bool includeUnion, int boundarySeparation);

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
namespace calc_strides
{

void first_dim_varies_fastest(int rank, const int *dims, int *strides);
void last_dim_varies_fastest(int rank, const int *dims, int *strides);

template<int rank>
static Eigen::Matrix<int, rank, 1> first_dim_varies_fastest(const Eigen::Matrix<int, rank, 1> &dims)
{
  Eigen::Matrix<int, rank, 1> strides;
  first_dim_varies_fastest(rank, dims.data(), strides.data());
  return strides;
}


template<int rank>
static Eigen::Matrix<int, rank, 1> last_dim_varies_fastest(const Eigen::Matrix<int, rank, 1> &dims)
{
  Eigen::Matrix<int, rank, 1> strides;
  last_dim_varies_fastest(rank, dims.data(), strides.data());
  return strides;
}


}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

struct Float3x4
{
  Float3x4() : X(0),Y(0),Z(0),T(0) {}
  explicit Float3x4(float s) : X(s,0,0),Y(0,s,0),Z(0,0,s),T(0) { }
  Float3 X,Y,Z,T;
  float Det() const;
  Float3& operator[](int i) { myAssert(i>=0 && i<3); return (&X)[i]; }
  const Float3& operator[](int i) const { myAssert(i>=0 && i<3); return (&X)[i]; }
  inline float& operator()( int i, int j ) { return ((float*)this)[j*3+i]; }
  void ScaleAxes(const Float3 &s) { X*=s.x(); Y*=s.y(); Z*=s.z(); }
};

const Float3x4 Inverse(const Float3x4 &M);
const Float3x4 InverseTranspose( const Float3x4 &M );
Float3 operator^( const Float3x4 &M, const Float3 &x );
inline Float3 operator*( const Float3x4 &M, const Float3 &x ) { return (M^x)+M.T; };
void MultMatrixToGL( const Float3x4 &M );
void GetMatrixFromGL( Float3x4 &M );
void HPBToMatrix(const Float3 &w, Float3x4 &M);
void Diagonal( float s, Float3x4 &D );
void Inverse( const Float3x4 &M, Float3x4 &R );
void Transpose( const Float3x4 &M, Float3x4 &R );
void InverseTranspose( const Float3x4 &M, Float3x4 &R );
void Eigenvalues( const Float3x4 &m, Float3 &ev, Float3x4 &o );
const Float3x4 operator*(const Float3x4 &a, const Float3x4 &b);
Float3x4 OrthogonalSystem(const Float3 &Z);

#endif
