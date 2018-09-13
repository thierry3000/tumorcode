/**
This file is part of tumorcode project.
(http://www.uni-saarland.de/fak7/rieger/homepage/research/tumor/tumor.html)

Copyright (C) 2018  Thierry Fredrich

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

/**
 * this file is dedicated to check out the behaviour of the new 
 * hdf5 cpp interface. I intend to transfer the result to the 
 * tumor code in order to get rid of ther HDF5cpp wrapper
 */

//#define _GLIBCXX_USE_CXX11_ABI 0

# include <iostream>
# include <string>
# include <boost/random.hpp>
# include <boost/random/random_device.hpp>  // creates real random number from system functions
# include <boost/multi_array.hpp>


#include <H5Cpp.h>
#include <eigen3/Eigen/Core>
#include "../../common/hdfio.h"

//create fake data
boost::random::random_device dev;         // produces randomness out of thin air
boost::random::uniform_01<> uni_f;                 // see pseudo-random number generators
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

int main() 
{
    string fn = "writeFile.h5";
    std::cout << "creating h5 file " << fn << std::endl;
    
    H5::H5File f = H5::H5File(fn, H5F_ACC_TRUNC);
    H5::Group g = f.createGroup("bla");
    int myInt = 42;
    float myFloat = 44.0;
    double myDouble = 42.0;
    string myString = string("a test string to");
    Float3 myFloat3= Float3(2.0,4.0,6.6);
    Int2 myInt2 = Int2(33,33);
    
    
    
    Float2 myFloat2= Float2(42.0,43.0);
    //Vec<float, 6> myFloat6=Vec<float, 6>(float(42),float(43),float(42),float(43),float(42),float(43));
    //Double6 myDouble6 =Double6(42.0,43.0,42.0,43.0,42.0,43.0);
    Double2 myDouble2= Double2(42.0,43.0);
    Int6 myInt6;
    
    Bool3 myBool3 = Bool3(false,true,false);
    bool myBool = false;
    //Eigen::Matrix<int, 6, 1> myInt6;
    myInt6[0] =42;
    myInt6[3] = 42;
    
    Eigen::Matrix<float, 2, 1> myBase;
    
    writeAttrToH5(g, string("test int"), myInt);
    writeAttrToH5(g, string("test float"), myFloat);
    writeAttrToH5(g, string("test double"), myDouble);
    writeAttrToH5(g, string("test string"), myString);
    writeAttrToH5(g, string("test Float3"), myFloat3);
    //writeAttrToH5(g, string("test Float2"), myFloat2);
    //writeAttrToH5(g, string("test Double2"), myDouble2);
    //writeAttrToH5(g, string("test Int2"), myInt2);
    writeAttrToH5(g, string("test myInt6"), myInt6);
    writeAttrToH5(g, string("test myBool3"), myBool3);
    writeAttrToH5(g, string("test myBool"), myBool);
    
    
    int result;
    readAttrFromH5(g, string("test int"), result);
    cout << result << endl;
    float result2;
    readAttrFromH5(g, string("test float"), result2);
    cout << result2 << endl;
    Float3 result3;
    readAttrFromH5(g, string("test Float3"), result3);
    cout << "1: " << result3[0] << endl;
    cout << "2: " << result3[1] << endl;
    cout << "3: " << result3[2] << endl;
    Bool3 result4;
    readAttrFromH5(g, string("test myBool3"), result4);
    cout << "1: " << result4[0] << endl;
    cout << "2: " << result4[1] << endl;
    cout << "3: " << result4[2] << endl;
    string readOutString;
    readAttrFromH5(g, string("test string"), readOutString);
    cout << "test string: " << readOutString <<endl;
    
    H5::Group myDataGroup = f.createGroup("data");
    
    const int n=20;
    DynArray<float> fakeData;
    DynArray<Float3> fakeData3D;
    DynArray<int> fakeDataInt;
    std::vector<bool> fakeDataBool;
    fakeData.resize(n);
    fakeData3D.resize(n);
    fakeDataInt.resize(n);
    fakeDataBool.resize(n);
    float avg =0;
    for (int i=0; i<n; ++i)
    {
      float b = uni_f(dev);
      bool c;
      if (b<0.5)
        c = true;
      else
        c = false;
      cout << c << endl;
      fakeDataBool[i] = c;
      fakeData[i]= b;
      fakeData3D[i] = Float3(uni_f(dev),uni_f(dev),uni_f(dev));
      fakeDataInt[i] = (int) 300*uni_f(dev);
      avg = avg+b;
      cout << b << endl;
    }
    cout << avg/n << endl;
    fakeData3D[10] = Float3(33,22,11);
    
    writeDataSetToGroup(myDataGroup, string("fakeData"), fakeData);
    writeDataSetToGroup(myDataGroup, string("fakeData3D"), fakeData3D);
    writeDataSetToGroup(myDataGroup, string("fakeDataInt"), fakeDataInt);
    writeDataSetToGroup(myDataGroup, string("fakeDataBool"), fakeDataBool);
    
    DynArray<float> readBack;
    readDataSetFromGroup(myDataGroup, string("fakeData"), readBack);
    cout << readBack[5] << endl;
    
    DynArray<Float3> readBack3D;
    readDataSetFromGroup(myDataGroup, string("fakeData3D"), readBack3D);
    cout << readBack3D[5] << endl;
    
    DynArray<int> readBackInt;
    readDataSetFromGroup(myDataGroup, string("fakeDataInt"), readBackInt);
    cout << readBackInt[5] << endl;
    
    f.close();
    return 0;
}

