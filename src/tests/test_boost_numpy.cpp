//compile line
//g++ ~/tumorcode/src/tests/test_boost_numpy.cpp -I/usr/include/python2.7/ -lboost_numpy -lboost_python -lpython2.7

#include <boost/python/numpy.hpp>
#include <iostream>

namespace p = boost::python;
namespace np = boost::python::numpy;

int main(int argc, char **argv)
{
  Py_Initialize();
  np::initialize();
  
  np::ndarray myTest = np::zeros(p::make_tuple(6), np::dtype::get_builtin<float>());
  std::cout << std::endl
            << "Python ndarray :" << p::extract<char const *>(p::str(myTest)) << std::endl;
  //int arr[] = {1,2,3,4,5};
  //int arr[2][5] = {{1,8,12,20,25}, {5,9,13,24,26}};
  int arr[5][2] = {{1,2},{3,4},{5,6},{7,8},{9,10}};
  std::cout << "C++ array created" << std::endl;
//   np::ndarray py_array = np::from_data(arr, np::dtype::get_builtin<int>(),
//                                      p::make_tuple(2,5),
//                                      p::make_tuple(sizeof(int),sizeof(int)),
//                                      p::object());
  np::ndarray py_array = np::from_data(arr, np::dtype::get_builtin<int>(),
                                     p::make_tuple(5,2),
                                     p::make_tuple(2*sizeof(int),sizeof(int)),
                                     p::object());
  
  std::cout << "C++ array :" << std::endl;
  for (int j=0;j<5;j++)
  {
    std::cout << arr[j][0] << arr[j][1] << ' ';
  }
  std::cout << std::endl
            << "Python ndarray :" << p::extract<char const *>(p::str(py_array)) << std::endl;
  //py_array.get_data()[4] = 5;//gives second element
  //py_array.get_data()[8] = 1;//gives second element
  
  std::cout << "t.f. test" << std::endl;
  std::cout << "shape()" << py_array.get_shape() << std::endl;
  std::cout << "strides()" << py_array.get_strides() << std::endl;
  std::cout << "get_nd()" << py_array.get_nd() << std::endl;
  
  std::cout << "strides dim 0: " << py_array.strides(0) << std::endl;
  std::cout << "strides dim 1: " << py_array.strides(1) << std::endl;
  std::cout << "strides dim 2: " << py_array.strides(2) << std::endl;
  //std::cout << p::extract<int>(py_array[1,2]) << std:: endl;
  std::cout << "Is the change reflected in the C++ array used to create the ndarray ? " << std::endl;
  for (int j = 0; j < 5; j++)
  {
    std::cout << arr[j][0]<< arr[j][1] << ' ';
  }
  std::cout<<std::endl;
  for (int j = 0; j < 5; j++)
  {
    std::cout << "j: " << j << ":" << std::endl << "py_array["<<j<<"][0]: " << p::extract<int>(py_array[j][0]) << "py_array["<<j<<"][1]: " << p::extract<int>(py_array[j][1])<< std::endl;
  }
//   std::cout << "try strides loop" << std::endl;
//   for (int j = 0; j < py_array.get_shape()[0]; j++)
//   {
//     std::cout << (int)py_array.get_data()[j*py_array.get_strides()[0]] << std::endl;
//     for(int jj = 0; jj< py_array.get_shape()[1]; jj++)
//     {
//       std::cout << (int)py_array.get_data()[(jj)*py_array.get_strides()[1]] << std::endl;
//       //std::cout << "j: " << j << " jj: " << jj << std::endl;
//     }
//     std::cout << std::endl;
//   }
  std::cout << std::endl
            << "Is the change reflected in the Python ndarray ?" << std::endl
            << p::extract<char const *>(p::str(py_array)) << std::endl;
}
