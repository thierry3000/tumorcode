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
#include "python_helpers.h"
#include <boost/property_tree/info_parser.hpp>


/*
h5cpp::Dataset PythonToCppDataset(const py::object &op_)
{
  py::object id_obj1 = py::getattr(op_, "id");
  py::object id_obj2 = py::getattr(id_obj1, "id");
  ssize_t id = py::extract<ssize_t>(id_obj2);
  return h5cpp::Dataset(id);
}
*/

using boost::property_tree::ptree;

ptree convertInfoStr(const py::str &param_info_str, const ptree &defaults)
{
  string param_string = py::extract<string>(param_info_str);
  std::istringstream param_stream(param_string);
  ptree pt_params = defaults, pt_tmp;
  boost::property_tree::read_info(param_stream, pt_tmp);
  boost::property_tree::update(pt_params, pt_tmp);
  return pt_params;
}



namespace mw_py_impl
{
// 
// template<class Details>
// struct H5FromPy : public Details
// {
//   static void Register()
//   {
//     boost::python::converter::registry::push_back(
//       &convertible,
//       &construct,
//       boost::python::type_id<typename Details::ResultType>());
//   }
// 
//   static void* convertible(PyObject* obj_ptr)
//   {
//     py::object h5py = boost::python::import("h5py"); // not very efficient ...
//     py::object cls = py::getattr(h5py, Details::ResultPyClass());
//     py::object o(py::handle<>(py::borrowed(obj_ptr)));
//     if (PyObject_IsInstance(o.ptr(), cls.ptr()))
//       return obj_ptr;
//     else
//       return NULL;
//   }
// 
//   static void construct(
//     PyObject* obj_ptr,
//     boost::python::converter::rvalue_from_python_stage1_data* data)
//   {
//     py::object o(py::handle<>(py::borrowed(obj_ptr)));
// 
//     py::object id_obj1 = py::getattr(o, "id");
//     py::object id_obj2 = py::getattr(id_obj1, "id");
//     int id = py::extract<int>(id_obj2);
// 
//     typedef typename Details::ResultType Result;
//     // Grab pointer to memory into which to construct the new QString
//     void* storage = ((boost::python::converter::rvalue_from_python_storage<Result>*)data)->storage.bytes;
//     //new (storage) Result(id);
//     Details::Construct(storage, id);
//     // Stash the memory chunk pointer for later use by boost.python
//     data->convertible = storage;
//   }
// };
// 
// 
// // this is a bit more general than needed, but still ...
// struct H5DetailFile
// {
//   typedef H5::H5File ResultType;
//   static const char* ResultPyClass() { return "File"; }
//   static void Construct(void* storage, int id) { new (storage) H5::H5File(id); }
// };
// 
// struct H5DetailGroup
// {
//   typedef H5::Group ResultType;
//   static const char* ResultPyClass() { return "Group"; }
//   static void Construct(void* storage, int id) { new (storage) H5::Group(id); }
// };
// 
// struct H5DetailDataset
// {
//   typedef H5::DataSet ResultType;
//   static const char* ResultPyClass() { return "Dataset"; }
//   static void Construct(void* storage, int id) { new (storage) H5::DataSet(id); }
// };



  

template<class T, int dim>
struct VecFromPy
{
    typedef Vec<T,dim> V;
    static void Register()
    {
      boost::python::converter::registry::push_back(
        &convertible,
        &construct,
        boost::python::type_id<V>());
    }

    static void* convertible(PyObject* obj_ptr)
    {
      if (!PySequence_Check(obj_ptr) || PySequence_Size(obj_ptr)!=dim)
        return 0;
      py::object o(py::handle<>(py::borrowed(obj_ptr)));
      for (int i=0; i<dim; ++i)
      {
        py::extract<T> ex(o[i]);
        if (!ex.check()) return 0;
      }
      return obj_ptr;
    }

    // Convert obj_ptr into a QString
    static void construct(
      PyObject* obj_ptr,
      boost::python::converter::rvalue_from_python_stage1_data* data)
    {
      py::object o(py::handle<>(py::borrowed(obj_ptr)));
      V value;
      for (int i=0; i<dim; ++i)
      {
        value[i] = py::extract<T>(o[i]);
      }

      // Grab pointer to memory into which to construct the new QString
      void* storage = ((boost::python::converter::rvalue_from_python_storage<V>*)data)->storage.bytes;
      new (storage) V(value);
      // Stash the memory chunk pointer for later use by boost.python
      data->convertible = storage;
    }
};

// template<class T, int dim>
// struct VecToPy
// {
//   typedef Vec<T,dim> V;
// 
//   static PyObject* convert(const V& p)
//   {
//     Py_ssize_t dims[1] = { dim };
//     
//     auto r = np::arrayt<T>(np::zeros(1, dims, np::getItemtype<T>()));
//     for (int i=0; i<dim; ++i)
//       r[i] = p[i];
//     return py::incref(r.getObject().ptr());
//   }
// 
//   static void Register()
//   {
//     py::to_python_converter<V, VecToPy<T,dim> >();
//   }
// };


// void exportH5Converters()
// {
//   mw_py_impl::H5FromPy<H5DetailDataset>::Register();
//   mw_py_impl::H5FromPy<H5DetailFile>::Register();
//   mw_py_impl::H5FromPy<H5DetailGroup>::Register();
// }


void exportVectorClassConverters()
{
  // register automatic conversion routines.
  mw_py_impl::VecFromPy<int, 3>::Register();
  mw_py_impl::VecFromPy<float, 3>::Register();
  mw_py_impl::VecFromPy<double, 3>::Register();
  mw_py_impl::VecFromPy<bool, 3>::Register();
//   mw_py_impl::VecToPy<int, 3>::Register();
//   mw_py_impl::VecToPy<float, 3>::Register();
//   mw_py_impl::VecToPy<double, 3>::Register();
//   mw_py_impl::VecToPy<bool, 3>::Register();
}

template<class T>
double checkedExtractFromDict(const py::dict &d, const char* name)
{
  try
  {
    return py::extract<double>(d.get(name));
  }
  catch (const py::error_already_set &e) 
  {
    std::cerr << format("unable to extract parameter '%s': ") % name;
    throw e; // don't every try to handle this!
  }
}


}//namespace mw_py_impl

bool PyCheckAbort()
{
  const int res = PyErr_CheckSignals();
  //std::cout << "py check abort called " << res << endl;
  return res != 0;
}

// static void PyMPI_OPENMPI_dlopen_libmpi(void)
// {
// void handle = 0;
// int mode = RTLD_NOW | RTLD_GLOBAL;
// /* GNU/Linux and others */
// #ifdef RTLD_NOLOAD
// mode |= RTLD_NOLOAD;
// #endif
// if (!handle) handle = dlopen("libmpi.so.20", mode);
// if (!handle) handle = dlopen("libmpi.so.12", mode);
// if (!handle) handle = dlopen("libmpi.so.1", mode);
// if (!handle) handle = dlopen("libmpi.so.0", mode);
// if (!handle) handle = dlopen("libmpi.so", mode);
// }
// 
// static int PyMPI_OPENMPI_MPI_Init(int *argc, char ***argv)
// {
// PyMPI_OPENMPI_dlopen_libmpi();
// return MPI_Init(argc, argv);
// }
// #undef MPI_Init
// #define MPI_Init PyMPI_OPENMPI_MPI_Init
// 
// static int PyMPI_OPENMPI_MPI_Init_thread(int *argc, char ***argv,
// int required, int *provided)
// {
// PyMPI_OPENMPI_dlopen_libmpi();
// return MPI_Init_thread(argc, argv, required, provided);
// }
// #undef MPI_Init_thread
// #define MPI_Init_thread PyMPI_OPENMPI_MPI_Init_thread
