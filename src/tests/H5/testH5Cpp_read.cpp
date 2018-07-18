//#define _GLIBCXX_USE_CXX11_ABI 0

#include <iostream>
#include <string>
using namespace std;

//#include <hdfio.h>
//#include <H5CompType.h>
// #include <hdf5.h>
// #include <hdf5_hl.h>

// #define H5std_string std::string

#include <H5Cpp.h>
#include <eigen3/Eigen/Core>

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
typedef Vec<int,2> Int2;
typedef Vec<double,6> Double6;


void writeAttrToH5(H5::Group g, const string &attr_name,  int &value)
{ 
  { 
    // Create new dataspace for attribute
    H5::DataSpace attr_dataspace = H5::DataSpace(H5S_SCALAR);
    H5::Attribute attr_out = g.createAttribute(attr_name, H5::PredType::NATIVE_INT, attr_dataspace);
    attr_out.write(H5::PredType::NATIVE_INT, &value);
  }
};
void writeAttrToH5(H5::Group g, const string &attr_name,  double &value)
{ 
  { 
    // Create new dataspace for attribute
    H5::DataSpace attr_dataspace = H5::DataSpace(H5S_SCALAR);
    H5::Attribute attr_out = g.createAttribute(attr_name, H5::PredType::NATIVE_DOUBLE, attr_dataspace);
    attr_out.write(H5::PredType::NATIVE_DOUBLE, &value);
  }
};
void writeAttrToH5(H5::Group g, const string &attr_name,  float &value)
{ 
  { 
    // Create new dataspace for attribute
    H5::DataSpace attr_dataspace = H5::DataSpace(H5S_SCALAR);
    H5::Attribute attr_out = g.createAttribute(attr_name, H5::PredType::NATIVE_FLOAT, attr_dataspace);
    attr_out.write(H5::PredType::NATIVE_FLOAT, &value);
  }
};
void writeAttrToH5(H5::Group g, const string &attr_name,  string &value)
{ 
  { 
    // Create new dataspace for attribute
    H5::DataSpace attr_dataspace = H5::DataSpace(H5S_SCALAR);

    // Create new string datatype for attribute
    H5::StrType strdatatype(H5::PredType::C_S1, 256); // of length 256 characters

    // Set up write buffer for attribute
    const H5std_string strwritebuf (value);

    // Create attribute and write to it
    H5::Attribute myatt_in = g.createAttribute(attr_name, strdatatype, attr_dataspace);
    myatt_in.write(strdatatype, strwritebuf); 
  }
};
template <class T>
void writeAttrToH5(H5::Group g, const string &attr_name,  T &value)
{ 
  { 
    int noOfEntries = value.size();
    std::cout << "found size: " << noOfEntries << std::endl;
//     if (typeid(int()) == value.value_type)
//     {
//       std::cout << "found int" << endl;
//     }
    //typeid(a).name()
    std::cout << "found value type" << typeid(value[0]).name() << std::endl;
    std::cout << "found value type" << typeid(value).name() << std::endl;
    int rank=2;
    const hsize_t dims[rank] = {1, noOfEntries};
    
    auto firstValue = value[0];
    std::cout << "fist value found by auto type: " << firstValue << endl;
    //std::cout << "fist value type id: " << typeid(firstValue) << endl;
    std::cout << "fist value type id name: " << typeid(firstValue).name() << endl;
    
    H5::DataSpace mspace( rank, dims);
    // Create new dataspace for attribute
    H5::DataSpace attr_dataspace = H5::DataSpace(rank, dims);
    
    //compare results 0, for equal strings!
    if(!string("f").compare(string(typeid(firstValue).name())))
    {
      H5::Attribute attr_out = g.createAttribute(attr_name, H5::PredType::NATIVE_FLOAT, attr_dataspace);
      attr_out.write(H5::PredType::NATIVE_FLOAT, &value);
    }
    if(!string("d").compare(string(typeid(firstValue).name())))
    {
      H5::Attribute attr_out = g.createAttribute(attr_name, H5::PredType::NATIVE_DOUBLE, attr_dataspace);
      attr_out.write(H5::PredType::NATIVE_DOUBLE, &value);
    }
    if(!string("i").compare(string(typeid(firstValue).name())))
    {
      H5::Attribute attr_out = g.createAttribute(attr_name, H5::PredType::NATIVE_INT, attr_dataspace);
      attr_out.write(H5::PredType::NATIVE_INT, &value);
    }
  }
};

template <class T>
T readAttrFromGroup(H5::Group g, string attr_name)
{
  try
  {
    // Turn off the auto-printing when failure occurs so that we can
    // handle the errors appropriately
    H5::Exception::dontPrint();
    
    H5::Attribute att_to_read = g.openAttribute(attr_name);
    T buffer =0;
    H5::DataType type = att_to_read.getDataType();
    att_to_read.read(type, &buffer);
    return buffer;
  }
  catch()
  {
    cout << "attribute " << attr_name << " not found in file: " << g.getFileName() << endl;
  }
}


int main() 
{
    string fn_read = "writeFile.h5";
    H5::H5File f(fn_read, H5F_ACC_RDONLY);
    int readIntFromFile = readAttrFromGroup<int>(f.openGroup("bla"), string("test int"));
    float readFloatFromFile = readAttrFromGroup<float>(f.openGroup("bla"), string("test float"));
    
    cout << "read int: " << readIntFromFile << endl;
    cout << "read float: " << readFloatFromFile << endl;
    
    f.close();
    return 0;
}