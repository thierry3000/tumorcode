#include <boost/python.hpp>
#include "test_tumorcode_hdfio.h"
#include <mpi.h>

char const* greet()
{
   return "hello, world";
}
bool is_mpi_on()
{
  int mpi_is_initialized = 0;
  MPI_Initialized(&mpi_is_initialized);
  if (mpi_is_initialized>0)
    return true;
  else
    return false;
}

BOOST_PYTHON_MODULE(py_h5_tumorcode_io)
{
    boost::python::def("greet", greet);
    boost::python::def("run_file_test", function_to_create_h5_file);
    boost::python::def("is_mpi_on", is_mpi_on);
}
