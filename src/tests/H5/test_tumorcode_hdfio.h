#ifndef TEST_TUMORCODE_HDFIO_H
#define TEST_TUMORCODE_HDFIO_H

# include <iostream>
# include <string>
# include <boost/random.hpp>
# include <boost/random/random_device.hpp>  // creates real random number from system functions
# include <boost/multi_array.hpp>


#include <H5Cpp.h>
#include <eigen3/Eigen/Core>
#include "../../common/hdfio.h"

//store fake data
#include "../../mwlib/dynamicarray.h"

// template<class T, int d>
// class Vec : public Eigen::Matrix<T, d, 1>
// {
//   public:
//     typedef Eigen::Matrix<T, d, 1> Base;
//     typedef T value_type;
//     // constructors
//     Vec() : Base(Base::Constant(0.)) {}
//     explicit Vec(const T &v) : Base() { this->setConstant(d, v); }
//     Vec(const T &a, const T &b) : Base(a, b) {}
//     Vec(const T &a, const T &b, const T &c) : Base(a, b, c) {}
//     template<typename OtherDerived >
//     Vec (const Eigen::MatrixBase< OtherDerived > &other) : Base(other) {}
// 
//     template<class FUNC>
//     static inline Vec<T, d> mapIndex(FUNC f)
//     {
//       Vec<T, d> ret;
//       for (int i=0; i<d; ++i)
//         ret[i] = f(i);
//       return ret;
//     }
// };
// typedef Vec<float, 3> Float3;
// typedef Vec<float, 2> Float2;
// typedef Vec<int,2> Int2;
// typedef Vec<int,3> Int3;
// typedef Vec<int,6> Int6;
// typedef Vec<double,6> Double6;
// typedef Vec<double,3> Double3;
// typedef Vec<double,2> Double2;
// typedef Vec<bool,3> Bool3;
// typedef u_char uchar;
// typedef int64_t int64;

using namespace std;
void function_to_create_h5_file();
#endif
