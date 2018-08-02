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
# include <boost/multi_array.hpp>
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
typedef Vec<double,3> Double3;
typedef Vec<double,2> Double2;
typedef Vec<bool,3> Bool3;
typedef u_char uchar;
typedef int64_t int64;

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
  else if(typeid(T) == typeid(Float2))
  {
    thisWritingType = H5::PredType::NATIVE_FLOAT;
  }
  else if(typeid(T) == typeid(double))
  {
    thisWritingType = H5::PredType::NATIVE_DOUBLE;
  }
  else if(typeid(T) == typeid(Double3))
  {
    thisWritingType = H5::PredType::NATIVE_DOUBLE;
  }
  else if(typeid(T) == typeid(Bool3))
  {
    thisWritingType = H5::PredType::NATIVE_CHAR;
  }
  else if(typeid(T) == typeid(int))
  {
    thisWritingType = H5::PredType::NATIVE_INT;
  }
  else if(typeid(T) == typeid(Int3))
  {
    thisWritingType = H5::PredType::NATIVE_INT;
  }
  else if(typeid(T) == typeid(Int6))
  {
    thisWritingType = H5::PredType::NATIVE_INT;
  }
  else if(typeid(T) == typeid(bool))
  {
    thisWritingType = H5::PredType::NATIVE_CHAR;
  }
  else if(typeid(T) == typeid(char))
  {
    thisWritingType = H5::PredType::NATIVE_CHAR;
  }
  else if(typeid(T) == typeid(uchar))
  {
    thisWritingType = H5::PredType::NATIVE_UCHAR;
  }
  else if(typeid(T) == typeid(long))
  {
    thisWritingType = H5::PredType::NATIVE_LONG;
  }
  else if(typeid(T) == typeid(int64))
  {
    thisWritingType = H5::PredType::NATIVE_INT64;
  }
  else if(typeid(T) == typeid(unsigned long))
  {
    thisWritingType = H5::PredType::NATIVE_UINT64;
  }
  else
  {
    cout << "unsupported Template type in writeDataSetToGroup!" << endl;
    exit(1);
  }
  return thisWritingType;
}


/**
 * one attribute could be an int3 or a point in 3D space
 */
#if H5_VERS_MINOR > 9
template <class T>
void writeAttrToH5(H5::H5Object &h, const string &attr_name,  const T &value)
{ 
  H5::DataType thisType = getH5TypeFromCpp<T>();
  const int rank = 2;
  hsize_t dims[rank];
  dims[0] = 1;
  if(typeid(T) == typeid(Float3) or typeid(T) == typeid(Int3) or typeid(T) == typeid(Bool3))
  {
    dims[1] = 3;
  }
  else
  {
    dims[1] = 1;
  }
  if(typeid(T) == typeid(Int6))
  {
    dims[1] = 6;
  }
  H5::DataSpace mspace = H5::DataSpace( rank, dims);
  H5::Attribute attr_out;
  try{
    attr_out = h.createAttribute(attr_name, thisType, mspace);
  }
  catch(H5::Exception &e)
  {
    std::cout << "unable for write: " << attr_name << std::endl;
    e.printErrorStack();
  }
  attr_out.write(thisType, &value);
};
#else //#if H5_VERS_MINOR > 9
template <class T>
void writeAttrToH5(H5::H5Location &h, const string &attr_name,  const T &value)
{ 
  H5::DataType thisType = getH5TypeFromCpp<T>();
  const int rank = 2;
  hsize_t dims[rank];
  dims[0] = 1;
  if(typeid(T) == typeid(Float3) or typeid(T) == typeid(Int3) or typeid(T) == typeid(Bool3))
  {
    dims[1] = 3;
  }
  else
  {
    dims[1] = 1;
  }
  if(typeid(T) == typeid(Int6))
  {
    dims[1] = 6;
  }
  
  H5::DataSpace mspace = H5::DataSpace( rank, dims);
  H5::Attribute attr_out;
  try
  {
    attr_out = h.createAttribute(attr_name, thisType, mspace);
  }
  catch(H5::Exception &e)
  {
    std::cout << "unable for write: " << attr_name << std::endl;
    e.printErrorStack();
  }
  attr_out.write(thisType, &value);
};
#endif //#if H5_VERS_MINOR > 9

#if H5_VERS_MINOR > 9
template<>
void writeAttrToH5<string>(H5::H5Object &h, const string &attr_name, const string &value)
{ 
  // Create new dataspace for attribute
  H5::DataSpace attr_dataspace = H5::DataSpace(H5S_SCALAR);
  // Create new string datatype for attribute
  H5::StrType strdatatype(H5::PredType::C_S1, H5T_VARIABLE); // of length 256 characters
  // Set up write buffer for attribute
  const H5std_string strwritebuf (value);
  try
  {
    H5::Attribute myatt_in = h.createAttribute(attr_name, strdatatype, attr_dataspace);
    myatt_in.write(strdatatype, strwritebuf);
  }
  catch(H5::Exception &e)
  {
    std::cout << "unable for write: " << attr_name << std::endl;
    e.printErrorStack();
  }
};
#else //#if H5_VERS_MINOR > 9
template<>
void writeAttrToH5<string>(H5::H5Location &h, const string &attr_name, const string &value)
{ 
  // Create new dataspace for attribute
  H5::DataSpace attr_dataspace = H5::DataSpace(H5S_SCALAR);
  // Create new string datatype for attribute
  H5::StrType strdatatype(H5::PredType::C_S1, H5T_VARIABLE); // of length 256 characters
  // Set up write buffer for attribute
  const H5std_string strwritebuf (value);
  
  try
  {
    H5::Attribute myatt_in = h.createAttribute(attr_name, strdatatype, attr_dataspace);
    myatt_in.write(strdatatype, strwritebuf);
  }
  catch(H5::Exception &e)
  {
    std::cout << "unable for write: " << attr_name << std::endl;
    e.printErrorStack();
  }
};
#endif //#if H5_VERS_MINOR > 9



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

#if H5_VERS_MINOR > 9
template<class T>
void readAttrFromH5(H5::H5Object &g, const string &attr_name, T &output_buffer)
{
  H5::Attribute att_to_read;
  H5::DataType type; 
  try
  {
    att_to_read = g.openAttribute(attr_name);
    type = att_to_read.getDataType();
    att_to_read.read(type, &output_buffer);
  }
  catch(H5::Exception &e)
  {
    std::cout << "unable for read: " << attr_name << std::endl;
    e.printErrorStack();
  }
}
#else //#if H5_VERS_MINOR > 9
template<class T>
void readAttrFromH5(H5::H5Location &g, const string &attr_name, T &output_buffer)
{
  H5::Attribute att_to_read;
  H5::DataType type; 
  try
  {
    att_to_read = g.openAttribute(attr_name);
    type = att_to_read.getDataType();
    att_to_read.read(type, &output_buffer);
  }
  catch(H5::Exception &e)
  {
    std::cout << "unable for read: " << attr_name << std::endl;
    e.printErrorStack();
  }
}
#endif //#if H5_VERS_MINOR > 9



#if H5_VERS_MINOR > 9
template<>
void readAttrFromH5<string>(H5::H5Object &g, const string &attr_name, string &output_buffer)
{ 
  H5::Attribute att_to_read = g.openAttribute(attr_name);
  
  // Create new dataspace for attribute
  H5::DataSpace attr_dataspace = H5::DataSpace(H5S_SCALAR);

  // Create new string datatype for attribute
  H5::StrType strdatatype(H5::PredType::C_S1, H5T_VARIABLE); // of length 256 characters

  // Set up read buffer for attribute
  H5std_string strreadbuf ("");

  // Create attribute and write to it
  try
  {
    att_to_read.read(strdatatype, strreadbuf);
  }
  catch(H5::Exception &e)
  {
    std::cout << "unable for read: " << attr_name << std::endl;
    e.printErrorStack();
  }
  output_buffer = strreadbuf;
}
#else //#if H5_VERS_MINOR > 9
template<>
void readAttrFromH5<string>(H5::H5Location &g, const string &attr_name, string &output_buffer)
{ 
  H5::Attribute att_to_read = g.openAttribute(attr_name);
  
  // Create new dataspace for attribute
  H5::DataSpace attr_dataspace = H5::DataSpace(H5S_SCALAR);

  // Create new string datatype for attribute
  H5::StrType strdatatype(H5::PredType::C_S1, H5T_VARIABLE); // of length 256 characters

  // Set up read buffer for attribute
  H5std_string strreadbuf ("");

  // Create attribute and write to it
  try
  {
    att_to_read.read(strdatatype, strreadbuf);
  }
  catch(H5::Exception &e)
  {
    std::cout << "unable for read: " << attr_name << std::endl;
    e.printErrorStack();
  }
  output_buffer = strreadbuf;
}
#endif //#if H5_VERS_MINOR > 9





template<class T>
H5::DataSet writeDataSetToGroup(H5::Group &g, const string &dataset_name, DynArray<T> &value)
{
  int sizeOfDynArray = value.size();
  boost::multi_array<T,1> continousMemoryArrary(boost::extents[sizeOfDynArray]);
  //T *continousMemoryArrary = new T[sizeOfDynArray];
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
  else if( typeid(T) == typeid(Float2))
  {
    dims[1] = 2;
  }
  else
  {
    dims[1] = 1;
  }
#ifndef NDEBUG
  cout << "we are writting data of size (" << dims[0] << ", " << dims[1] << "!" <<endl;
#endif
  H5::DataSet ds;
  try 
  {
    H5::DataSpace dataspace( rank, dims);
    ds = g.createDataSet(dataset_name, thisWritingType, dataspace);
    ds.write(continousMemoryArrary.data(), thisWritingType, dataspace);
  }
  catch( H5::Exception &e)
  {
    cout<< "failed to write " << dataset_name << endl;
    e.printErrorStack();
  }
  return ds;
}



template<class T>
void writeDataSetToGroup(H5::Group &g, const string &dataset_name, const std::vector<T> &value)
{
  int sizeOfvector = value.size();
  int rank = 2;
  hsize_t dims[rank];
  dims[0]=sizeOfvector;
  H5::DataType thisWritingType = getH5TypeFromCpp<T>();
  dims[1] = 1;
#ifndef NDEBUG
  cout << "value[0]: " << value[0] << endl;
  cout<< "writing std::vector: " << dataset_name << " to hdf5" << endl;
  cout << "we are writting data of size (" << dims[0] << ", " << dims[1] << "!" <<endl;
#endif
  H5::DataSet ds;
  try 
  {
    H5::DataSpace dataspace( rank, dims);
    ds = g.createDataSet(dataset_name, thisWritingType, dataspace);
    ds.write(value.data(), thisWritingType, dataspace);
  }
  catch( H5::Exception &e)
  {
    cout<< "failed to write " << dataset_name << endl;
    e.printErrorStack();
  }
}

/** bool s work differently
 * 
 * untested!!!
 * 
 */
template<>
void writeDataSetToGroup(H5::Group &g, const string &dataset_name, const std::vector<bool> &value)
{
  int sizeOfvector = value.size();
  int rank = 2;
  hsize_t dims[rank];
  dims[0]=sizeOfvector;
  H5::DataType thisWritingType = getH5TypeFromCpp<bool>();
  dims[1] = 1;
#ifndef NDEBUG
  cout<< "writing std::vector: " << dataset_name << " to hdf5" << endl;
  cout << "we are writting data of size (" << dims[0] << ", " << dims[1] << "!" <<endl;
  cout << "value[0]: " << value[0] << endl;
#endif
  H5::DataSet ds;
  try 
  {
    H5::DataSpace dataspace( rank, dims);
    ds = g.createDataSet(dataset_name, thisWritingType, dataspace);
    ds.write(&value, thisWritingType, dataspace);
  }
  catch( H5::Exception &e)
  {
    cout<< "failed to write " << dataset_name << endl;
    e.printErrorStack();
  }
}

template <class T>
void readDataSetFromGroup(H5::Group &g, const string &dataset_name, DynArray<T> &readIn)
{
  H5::DataSet dset;
  H5::DataSpace dataspace;
  
  try{
    dset = g.openDataSet(dataset_name);
    dataspace = dset.getSpace();
    hsize_t rank = dataspace.getSimpleExtentNdims();
    hsize_t dims[rank];
    hsize_t max_dims[rank];
    hsize_t again = dataspace.getSimpleExtentDims(dims,max_dims);
#ifndef NDEBUG
    cout << "rank: " << rank << endl;
    for(int i=0; i<rank; i++)
    {
      cout << "dims[" << i << "]" << " = " << dims[i] << endl;
    }
#endif
  
    boost::multi_array<T,1> arr_data(boost::extents[dims[0]]);
    H5::DataType thisWritingType = getH5TypeFromCpp<T>();

    dset.read(arr_data.data(), thisWritingType);
    readIn.resize(dims[0]);
    
    for( int i=0;i<dims[0];++i)
    {
      readIn[i] = arr_data[i];
    }
  }
  catch(H5::Exception &e)
  {
    cout << "failed to read: " << dataset_name << endl;
    e.printErrorStack();
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
