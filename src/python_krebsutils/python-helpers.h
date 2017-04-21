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
#ifndef PYTHON_HELPERS_H
#define PYTHON_HELPERS_H

#include <boost/python.hpp>
#include <boost/python/tuple.hpp>
#include <boost/format.hpp>
#include <fenv.h>

#include "mwlib/helpers-vec.h"
#include "hdf_wrapper.h"
#include "mwlib/math_ext.h"
#include "mwlib/ptree_ext.h"

#include "numpy.hpp"

namespace py = boost::python;
namespace np = boost::python::numpy;
namespace nm = boost::python::numeric;
using boost::format;
using boost::str;



template<class T, int dim>
inline BBoxd<T, dim> BBoxFromPy(const py::object &p)
{
  BBoxd<T, dim> r;
  for (int i=0; i<dim; ++i)
  {
    r.min[i] = py::extract<T>(p[i*2]);
    r.max[i] = py::extract<T>(p[i*2+1]);
  }
  return r;
}

template<class T, int dim>
inline py::object BBoxToPy(const BBoxd<T, dim> &bb)
{
  np::ssize_t dims[1] = { dim*2 };
  np::arrayt<T> r = np::empty(1, dims, np::getItemtype<T>());
  for (int i=0; i<dim; ++i)
  {
    r[i*2  ] = bb.min[i];
    r[i*2+1] = bb.max[i];
  }
  return r.getObject();
}


// this is crazy shit! python generates a string representation of the entire
// configuration (.info file), then it is passed to the python function where
// it is parsed into a ptree object which is used to initialize the model.
boost::property_tree::ptree convertInfoStr(const py::str &param_info_str, const boost::property_tree::ptree &defaults);



static void _CheckArray(const py::object &obj, int itemtype, int rank, const int *dims, const char* file, const int line )
{
  np::arraytbase a(obj);
  if(rank>0 && a.rank()!=rank) throw std::invalid_argument(str(format("array rank %i expected, got %i @ %s(%i)") % rank % a.rank() % file % line));
  //if(itemtype>0 && a.itemtype()!=itemtype) throw std::invalid_argument(str(format("array itemtype %i expected, got %i @ %s(%i)") % itemtype % a.itemtype() % file % line));
  //TODO: fix itemtypes. For some reason a python array of type np.float64 gives a different itemtype than NPY_DOUBLE ...!!!!!
  if(dims)
  {
    const auto* adims = a.shape();
    for(int i=0; i<a.rank(); ++i)
    {
      if(dims[i]>0 && adims[i]!=dims[i])
      {
        throw std::invalid_argument(str(format("array dims[%i] %i expected, got %i @ %s(%i)") %  i % dims[i] % adims[i] % file % line));
      }
    }
  }
}

#define CheckArray(a,itemtype,rank,dims) _CheckArray(a,itemtype,rank,dims,__FILE__,__LINE__)

template<class T>
static void _CheckArray1(const py::object &obj, int dim1, const char* file, const int line)
{
  const int dims[] = { dim1 };
  _CheckArray(obj, np::getItemtype<T>(), 1, dims, file, line);
}

template<class T>
static void _CheckArray2(const py::object &obj, int dim1, int dim2, const char* file, const int line)
{
  const int dims[] = { dim1, dim2 };
  _CheckArray(obj, np::getItemtype<T>(), 2, dims, file, line);
}

template<class T>
static void _CheckArray3(const py::object &obj, int dim1, int dim2, int dim3, const char* file, const int line)
{
  const int dims[] = { dim1, dim2, dim3 };
  _CheckArray(obj, np::getItemtype<T>(), 3, dims, file, line);
}

#define CheckArray1(T, a,dim1) _CheckArray1<T>(a, dim1, __FILE__,__LINE__)
#define CheckArray2(T, a,dim1,dim2) _CheckArray2<T>(a, dim1, dim2, __FILE__,__LINE__)
#define CheckArray3(T, a,dim1,dim2,dim3) _CheckArray3<T>(a, dim1, dim2, dim3, __FILE__,__LINE__)

// #define THROW_ERROR(msg)        throw std::invalid_argument(str(format(msg + " @ %s(%i)") % __FILE__ % __LINE__))
// #define THROW_ERROR1(msg,a1)    throw std::invalid_argument(str(format(msg + " @ %s(%i)") % a1 % __FILE__ % __LINE__))
// #define THROW_ERROR2(msg,a1,a2) throw std::invalid_argument(str(format(msg + " @ %s(%i)") % a1 % a2 % __FILE__ % __LINE__))

/**brief With these functions we can take hdf5 entities (group and datasets) as arguments to python functions and use them on the c++ side.
   
   TODO: This approach is however flawed because the default h5py install might use different hdf5 libraries than we link tumorcode against. It might therefore be a good idea to rewrite regarding functions. They could take filename-path pairs or memory buffers as arguments instead.   
*/
/** even worse:
 * this function does not seem to be hdf5 1.10 ready. got some strange errors with it.
 * solved:
 * this was due to extraction of 32bit int as type, now size_t is used
 */
h5cpp::Group PythonToCppGroup(const py::object &op_);
h5cpp::Dataset PythonToCppDataset(const py::object &op_);

#include "mwlib/field.h"

// wrap data from a numpy array in a Array3d<T>
// the numpy array must be kept alive while this array exists
template<class T>
inline Array3d<T> Array3dFromPy(np::arrayt<T> &a)
{
  assert(a.rank() <= 3);
  Int3 size(1);
  for (int i=0; i<a.rank(); ++i) size[i] = a.shape()[i];
  Array3d<T> arr3d(size,
                   Int3(a.strides()[0],a.strides()[1],a.strides()[2])/sizeof(T),
                       (T*)a.bytes(), size.prod(),false);
  return arr3d;
}


struct FpExceptionStateGuard
{
  fenv_t fp_state;
  FpExceptionStateGuard(int excepts = FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW)
  {
    fegetenv(&fp_state);
    feenableexcept(excepts);
  }
  ~FpExceptionStateGuard()
  {
    fesetenv(&fp_state);
  }
};

// Allow the python interpreter to do things (multithreading) while c stuff is running
// Shameless copy & paste from http://stackoverflow.com/questions/8009613/boost-python-not-supporting-parallelism
// BEWARE: It is VERY important to note that you MUST NOT touch any python code, or python data
// or call in to the interpreter while not holding the GIL. This will cause your interpreter to crash. 
class ReleaseGIL{
public:
    inline ReleaseGIL(){
        save_state = PyEval_SaveThread();
    }

    inline ~ReleaseGIL() {
        PyEval_RestoreThread(save_state);
    }
private:
    PyThreadState *save_state;
};


bool PyCheckAbort();

#endif
