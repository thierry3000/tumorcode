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
#ifndef COMMON_H
#define COMMON_H
//#define mwOMP
#include <omp.h>
#include "mwlib/helpers-defs.h"
#include "mwlib/myAssert.h"
#include "mwlib/helpers-vec.h"
#include "mwlib/math_ext.h"
#include "mwlib/vector_ext.h"
#include "mwlib/dynamicarray.h"
#include "mwlib/field.h"
#include "mwlib/timer.h"

#include "mwlib/ptree_ext.h"
#include "trilinos_helper.h"
#include <boost/format.hpp>
#include <boost/optional.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/asio/ip/host_name.hpp> // get host name

#include <cstdlib> // std::getenv


using boost::optional;
using boost::tie;
using boost::str;
using boost::format;
using boost::property_tree::ptree;
using boost::property_tree::make_ptree;

namespace my
{
  // ininitalize mpi and multithreading
//   class MultiprocessingInitializer
//   {
//     // usage: construct this at the beginning of the program and keep it alive.
//     // Call Get/SetNumThreads somewhere in your program.
//     public:
//       MultiprocessingInitializer(int argc, char **argv, int num_threads = 1);
//       MultiprocessingInitializer(int num_threads = 1);
//       ~MultiprocessingInitializer();
//   };
//   void initMultithreading(int argc, char **argv, int num_threads = 1);
//   int GetNumThreads();
//   void SetNumThreads(int n);
//   int OmpGetCurrentThread();
//   int OmpGetMaxThreadCount();
//   bool MultiprocessingInitializer_exists();
}


namespace my
{
  namespace my_log_impl {
    class Log;
  }
  my_log_impl::Log& log();

  /* This is used to abort lenghy computations on key pressed.
     In a normal program it goes automatically without calling
     a special check function but not so in python modules. So to
     make keyboard interrupts work for python first the module init
     code must set checkAbort to something suitable and the computation
     code must call it to see if it should stop. */
  typedef bool (*AbortFunction)();
  extern AbortFunction checkAbort;
}

struct SystemParameters
{
  int num_threads;
  string cluster;
  string computing_node;
  int num_threads_queuing;
  double mem_in_GB;
  
  bool isRerun;
  int reRunNumber;
  
  SystemParameters();
  void assign(const ptree &pt);
  ptree as_ptree() const;
};

void readSystemParameters(SystemParameters &sysParamsToFill);

#endif
