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

/**
 * this file is dedicated to check out the behaviour of the new 
 * hdf5 cpp interface. I intend to transfer the result to the 
 * tumor code in order to get rid of ther HDF5cpp wrapper
 */

//#define _GLIBCXX_USE_CXX11_ABI 0

#include <iostream>
#include <string>
# include <boost/random.hpp>
# include <boost/random/random_device.hpp>  // creates real random number from system functions
using namespace std;

#include <H5Cpp.h>
#include <eigen3/Eigen/Core>

//create fake data
boost::random::random_device dev;         // produces randomness out of thin air
boost::random::uniform_01<> uni_f;                 // see pseudo-random number generators
//store fake data
#include "../../mwlib/dynamicarray.h"

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
typedef Vec<float, 3> Float3;
typedef Vec<float, 2> Float2;
typedef Vec<int,2> Int2;
typedef Vec<int,2> Int3;
typedef Vec<int,6> Int6;
typedef Vec<double,6> Double6;
typedef Vec<double,2> Double2;
typedef Vec<bool,3> Bool3;


void writeAttrToH5(H5::Group g, const string &attr_name,  const int &value)
{ 
  { 
    // Create new dataspace for attribute
    H5::DataSpace attr_dataspace = H5::DataSpace(H5S_SCALAR);
    H5::Attribute attr_out = g.createAttribute(attr_name, H5::PredType::NATIVE_INT, attr_dataspace);
    attr_out.write(H5::PredType::NATIVE_INT, &value);
  }
};
void writeAttrToH5(H5::Group g, const string &attr_name,  const double &value)
{ 
  { 
    // Create new dataspace for attribute
    H5::DataSpace attr_dataspace = H5::DataSpace(H5S_SCALAR);
    H5::Attribute attr_out = g.createAttribute(attr_name, H5::PredType::NATIVE_DOUBLE, attr_dataspace);
    attr_out.write(H5::PredType::NATIVE_DOUBLE, &value);
  }
};
void writeAttrToH5(H5::Group g, const string &attr_name,  const float &value)
{ 
  { 
    // Create new dataspace for attribute
    H5::DataSpace attr_dataspace = H5::DataSpace(H5S_SCALAR);
    H5::Attribute attr_out = g.createAttribute(attr_name, H5::PredType::NATIVE_FLOAT, attr_dataspace);
    attr_out.write(H5::PredType::NATIVE_FLOAT, &value);
  }
};
void writeAttrToH5(H5::Group g, const string &attr_name, const  bool &value)
{ 
  { 
    // Create new dataspace for attribute
    H5::DataSpace attr_dataspace = H5::DataSpace(H5S_SCALAR);
    H5::Attribute attr_out = g.createAttribute(attr_name, H5::PredType::NATIVE_UCHAR, attr_dataspace);
    attr_out.write(H5::PredType::NATIVE_UCHAR, &value);
  }
};
void writeAttrToH5(H5::Group g, const string &attr_name, const  string &value)
{ 
  { 
    // Create new dataspace for attribute
    H5::DataSpace attr_dataspace = H5::DataSpace(H5S_SCALAR);
    
    H5::StrType strdatatype(H5::PredType::C_S1, H5T_VARIABLE); // of length 256 characters
    
    const H5std_string strwritebuf (value);
    
    try{
      H5::Attribute myatt_in = g.createAttribute(attr_name, strdatatype, attr_dataspace);
      myatt_in.write(strdatatype, strwritebuf);
    }
    catch(H5::Exception error)
    {
      error.printErrorStack();
    }
   
  }
};

void writeAttrToH5(H5::Group h, const string &attr_name,  const Float3 &value)
{ 
  {
    const int rank = 2;
    hsize_t dims[rank] = {1,3};
    H5::DataSpace mspace( rank, dims);
    H5::Attribute attr_out = h.createAttribute(attr_name, H5::PredType::NATIVE_FLOAT, mspace);
    
    attr_out.write(H5::PredType::NATIVE_FLOAT, &value);
  }
};

void writeAttrToH5(H5::Group h, const string &attr_name,  const Float2 &value)
{ 
  {
    const int rank = 2;
    hsize_t dims[rank] = {1,2};
    H5::DataSpace mspace( rank, dims);
    H5::Attribute attr_out = h.createAttribute(attr_name, H5::PredType::NATIVE_FLOAT, mspace);
    
    attr_out.write(H5::PredType::NATIVE_FLOAT, &value);
  }
};
void writeAttrToH5(H5::Group h, const string &attr_name,  const Double2 &value)
{ 
  {
    const int rank = 2;
    hsize_t dims[rank] = {1,2};
    H5::DataSpace mspace( rank, dims);
    H5::Attribute attr_out = h.createAttribute(attr_name, H5::PredType::NATIVE_DOUBLE, mspace);
    
    attr_out.write(H5::PredType::NATIVE_DOUBLE, &value);
  }
};

void writeAttrToH5(H5::Group h, string attr_name,  Int6 &value)
{ 
  {
    const int rank = 2;
    hsize_t dims[rank] = {1,6};
    
    H5::DataSpace mspace( rank, dims);
    H5::Attribute attr_out = h.createAttribute(attr_name, H5::PredType::NATIVE_INT, mspace);
    
    attr_out.write(H5::PredType::NATIVE_INT, &value);
  }
};

void writeAttrToH5(H5::Group h, string attr_name,  Int2 &value)
{ 
  {
    const int rank = 2;
    hsize_t dims[rank] = {1,2};
    
    H5::DataSpace mspace( rank, dims);
    H5::Attribute attr_out = h.createAttribute(attr_name, H5::PredType::NATIVE_INT, mspace);
    
    attr_out.write(H5::PredType::NATIVE_INT, &value);
  }
};

void writeAttrToH5(H5::Group h, string attr_name,  Bool3 &value)
{ 
  {
    const int rank = 2;
    hsize_t dims[rank] = {1,3};
    H5::DataSpace mspace( rank, dims);
    H5::Attribute attr_out = h.createAttribute(attr_name, H5::PredType::NATIVE_CHAR, mspace);
    
    attr_out.write(H5::PredType::NATIVE_CHAR, &value);
  }
};



// template <class T>
// void writeAttrToH5(H5::Group g, const string &attr_name, const  T &value)
// { 
//   { 
//     int noOfEntries = value.size();
//     std::cout << "found size: " << noOfEntries << std::endl;
// //     if (typeid(int()) == value.value_type)
// //     {
// //       std::cout << "found int" << endl;
// //     }
//     //typeid(a).name()
//     std::cout << "found value type" << typeid(value[0]).name() << std::endl;
//     std::cout << "found value type" << typeid(value).name() << std::endl;
//     int rank=2;
//     const hsize_t dims[rank] = {1, noOfEntries};
//     
//     auto firstValue = value[0];
//     std::cout << "fist value found by auto type: " << firstValue << endl;
//     //std::cout << "fist value type id: " << typeid(firstValue) << endl;
//     std::cout << "fist value type id name: " << typeid(firstValue).name() << endl;
//     
//     H5::DataSpace mspace( rank, dims);
//     // Create new dataspace for attribute
//     H5::DataSpace attr_dataspace = H5::DataSpace(rank, dims);
//     
//     //compare results 0, for equal strings!
//     if(!string("f").compare(string(typeid(firstValue).name())))
//     {
//       H5::Attribute attr_out = g.createAttribute(attr_name, H5::PredType::NATIVE_FLOAT, attr_dataspace);
//       attr_out.write(H5::PredType::NATIVE_FLOAT, &value);
//     }
//     if(!string("d").compare(string(typeid(firstValue).name())))
//     {
//       H5::Attribute attr_out = g.createAttribute(attr_name, H5::PredType::NATIVE_DOUBLE, attr_dataspace);
//       attr_out.write(H5::PredType::NATIVE_DOUBLE, &value);
//     }
//     if(!string("i").compare(string(typeid(firstValue).name())))
//     {
//       H5::Attribute attr_out = g.createAttribute(attr_name, H5::PredType::NATIVE_INT, attr_dataspace);
//       attr_out.write(H5::PredType::NATIVE_INT, &value);
//     }
//     if(!string("b").compare(string(typeid(firstValue).name())))
//     {
//       H5::Attribute attr_out = g.createAttribute(attr_name, H5::PredType::NATIVE_CHAR, attr_dataspace);
//       attr_out.write(H5::PredType::NATIVE_CHAR, &value);
//     }
//   }
// };

// template<class T>
// T readAttrFromH5(H5::Group g, string attr_name)
// {
//   H5::Attribute att_to_read = g.openAttribute(attr_name);
//   T buffer =0;
//   H5::DataType type = att_to_read.getDataType();
//   att_to_read.read(type, &buffer);
//   return buffer;
// }

/* this uses c++14 features! */
// auto readAttrFromH5(H5::Group g, string attr_name)
// {
//   H5::Attribute att_to_read = g.openAttribute(attr_name);
//   H5::DataType type = att_to_read.getDataType();
//   if(type == H5::PredType::NATIVE_INT)
//   {
//     cout << "reading an integer" <<endl;
//     int buffer =0;
//     att_to_read.read(type, &buffer);
//     return buffer;
//   }
//   else if(type == H5::PredType::NATIVE_FLOAT)
//   {
//     cout << "reading a float" <<endl;
//     float buffer =0;
//     att_to_read.read(type, &buffer);
//     return buffer;
//   }
//   else
//   {
//     cout << "warning: unsupported hdf read out" << endl;
//     return 0;
//   }
//   
// }

template<class T>
void readAttrFromH5(H5::Group g, const string &attr_name, T &output_buffer)
{
  H5::Attribute att_to_read = g.openAttribute(attr_name);
  H5::DataType type = att_to_read.getDataType();
  try
  {
    att_to_read.read(type, &output_buffer);
  }
  catch(H5::Exception error)
  {
    error.printErrorStack();
  }
}
template<>
void readAttrFromH5<string>(H5::Group g, const string &attr_name, string &output_buffer)
{ 
  H5::Attribute att_to_read = g.openAttribute(attr_name);
  //H5::DataType type = att_to_read.getDataType();
  
  // Create new dataspace for attribute
  H5::DataSpace attr_dataspace = H5::DataSpace(H5S_SCALAR);

  // Create new string datatype for attributeH5T_VARIABLE
  H5::StrType strdatatype(H5::PredType::C_S1, H5T_VARIABLE); // of length 256 characters

  // Set up read buffer for attribute
  H5std_string strreadbuf ("");

  // Create attribute and write to it
  att_to_read.read(strdatatype, strreadbuf); 
  output_buffer = strreadbuf;
  cout << "read string works" << endl;
}
template<class T>
H5::DataType getH5TypeFromCpp()
{
  H5::DataType thisWritingType;
  if(typeid(T) == typeid(float))
  {
    thisWritingType = H5::PredType::NATIVE_FLOAT;
  }
  else if(typeid(T) == typeid(Float3))
  {
    thisWritingType = H5::PredType::NATIVE_FLOAT;
  }
  else if(typeid(T) == typeid(int))
  {
    thisWritingType = H5::PredType::NATIVE_INT;
  }
  else if(typeid(T) == typeid(Int3))
  {
    thisWritingType = H5::PredType::NATIVE_INT;
  }
  else if(typeid(T) == typeid(bool))
  {
    thisWritingType = H5::PredType::NATIVE_CHAR;
  }
  else
  {
    cout << "unsupported Template type in writeDataSetToGroup!" << endl;
    exit(1);
  }
  return thisWritingType;
}
template<class T>
H5::DataSet writeDataSetToGroup(H5::Group g, const string &dataset_name, DynArray<T> &value)
{
  int sizeOfDynArray = value.size();
  T continousMemoryArrary[sizeOfDynArray];
  for( int i=0; i<sizeOfDynArray;++i)
  {
    continousMemoryArrary[i] = T(value[i]);
  }
  int rank = 2;
  hsize_t dims[rank];
  dims[0]=sizeOfDynArray;
  H5::DataType thisWritingType = getH5TypeFromCpp<T>();
  if(typeid(T) == typeid(Float3) or typeid(T) == typeid(Int3) or typeid(T) == typeid(Bool3))
  {
    dims[1] = 3;
  }
  else
  {
    dims[1] = 1;
  }
#ifdef DEBUG
  cout << "we are writting data of size (" << dims[0] << ", " << dims[1] << "!" <<endl;
#endif
  H5::DataSpace dataspace( rank, dims);
  H5::DataSet ds = g.createDataSet(dataset_name, thisWritingType, dataspace);
  ds.write(&continousMemoryArrary, thisWritingType, dataspace);
  return ds;
}

template <class T>
void readDataSetFromGroup(H5::Group g, const string &dataset_name, DynArray<T> &readIn)
{
  H5::DataSet dset = g.openDataSet(dataset_name);
  H5::DataSpace dataspace = dset.getSpace();
  hsize_t dims[dataspace.getSimpleExtentNdims()];
  cout << dataspace.getSimpleExtentDims(dims) << endl;
  cout << dims[0] << endl;
  cout << dims[1] << endl;
  T arr[dims[0]];
  H5::DataType thisWritingType = getH5TypeFromCpp<T>();
  dset.read(&arr, thisWritingType);
  readIn.resize(dims[0]);
  for( int i=0;i<dims[0];++i)
  {
    readIn[i] = arr[i];
  }
}

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
    writeAttrToH5(g, string("test Float2"), myFloat2);
    writeAttrToH5(g, string("test Double2"), myDouble2);
    writeAttrToH5(g, string("test Int2"), myInt2);
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
    fakeData.resize(n);
    fakeData3D.resize(n);
    fakeDataInt.resize(n);
    float avg =0;
    for (int i=0; i<n; ++i)
    {
      float b = uni_f(dev);
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
