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
#include <Eigen/Core>
//#include <Eigen/Array>
#include <iostream>

// import most common Eigen types
//USING_PART_OF_NAMESPACE_EIGEN

template<class T, int d>
class Vec : public Eigen::Matrix<T, d, 1>
{
  public:
    typedef Eigen::Matrix<T, d, 1> Base;
    typedef T value_type;
    // constructors
    Vec() : Base(Base::Constant(0.)) {}
    explicit Vec(const T &v) : Base() { setConstant(d, v); }
    Vec(const T &a, const T &b) : Base(a, b) {}
    Vec(const T &a, const T &b, const T &c) : Base(a, b, c) {}
    template<typename OtherDerived >
    Vec (const Eigen::MatrixBase< OtherDerived > &other) : Base(other) {}
};


#define VECDARG class T,int d
#define VECD Vec<T, d>

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


template<class Derived>
inline typename Eigen::internal::plain_matrix_type<Derived>::type foobar(const Eigen::MatrixBase<Derived>& v, int x)
{
}

template<class Derived, class DerivedB>
inline auto cwiseL(const Eigen::MatrixBase<Derived> &a, const Eigen::MatrixBase<DerivedB> &b) -> decltype((a.array() < b.array()).matrix().eval())
{
  return (a.array() < b.array()).matrix();
}


template<class DerivedA, class DerivedB>
Eigen::Matrix<bool, DerivedA::RowsAtCompileTime, DerivedA::ColsAtCompileTime> cwiseLV(const Eigen::MatrixBase<DerivedA> &a, const Eigen::MatrixBase<DerivedB> &b)
{
  return (a.array() < b.array()).matrix();
}


// template<class Derived, class DerivedB>
// inline const Eigen::MatrixWrapper<
//   const Eigen::CwiseBinaryOp<
//     std::less<float>,  
//     const Eigen::ArrayWrapper< const Derived >, 
//     const Eigen::ArrayWrapper< const DerivedB > 
//   >
// > cwiseL(const Eigen::MatrixBase<Derived> &a, const Eigen::MatrixBase<DerivedB> &b)
// {
//   return (a.array() < b.array()).matrix();
// }

template<class Derived>
inline const Eigen::ArrayWrapper<const Derived> forwardArray(const Eigen::MatrixBase<Derived> &a)
{
  return a.array();
}


int main(int, char *[])
{
  using namespace Eigen;
  using namespace std;
  
  Matrix3f m3;
  m3 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  Matrix4f m4 = Matrix4f::Identity();
  Vector4i v4(1, 2, 3, 4);

  Vec<float, 3> q(4,5,100);
  foobar(q, 5);
 
  Vec<float, 3> r(10,10,3);
  Vec<float,3> s(5.);
  //Vec<float, 3> comp = 
  //cwiseL(q, r + Vec<float,3>(5.));
  //cwiseL(q, r);
  forwardArray(r+q);
    
  cout << "m3\n" << m3 << "\nm4:\n"
    << m4 << "\nv4:\n" << v4 << endl;
  
  //decltype((q.array() < q.array()).matrix()) x;
  //what_type<Vec<float ,3> >::type();
  
  r.array() >= 5.f;
    
  cout << "----" << q << " " << Vec<float,3>(r+s) << endl;
  Vec<bool,3> comp(cwiseL(q, r+s));
  cout << (cwiseL(q, r+s)).eval() << endl;
  cout << comp << endl;
  cout << (q.array() < (r).array()).matrix() << endl;
  cout << cwiseLV(q, r) << endl;
}
