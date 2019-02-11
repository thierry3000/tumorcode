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
#include "common.h"
#include <tbb/task_scheduler_init.h>

#if defined( _OPENMP )
  #include <omp.h>
#endif

// //it is good to know the epetra settings before we set up mpi and multithreading
#include <Epetra_ConfigDefs.h>

#include "mwlib/log.h"


namespace my
{

// struct MPPrivate
// {
//   // constructor
//   MPPrivate(int argc, char **argv, int num_threads_) : num_threads(0)
//   {
// #ifdef EPETRA_MPI
//     {
// #ifndef TOTAL_SILENCE
//       printf("\n\nEPETRA_MPI flag is set!\n");
// #endif
//       int mpi_is_initialized = 0;
//       int prov;
//       MPI_Initialized(&mpi_is_initialized);
//       if (!mpi_is_initialized)
// 	MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE,&prov);
// 	//kind of works with 
// 	//MPI_Init(&argc, &argv);
//     }
// #endif
//     SetNumThreads(num_threads_);
//   }
//   MPPrivate(int num_threads_) : num_threads(0)
//   {
//     SetNumThreads(num_threads_);
//   }
//   // destructor
// //   ~MPPrivate()
// //   {
// // #ifdef EPETRA_MPI
// //     int mpi_is_finalized = 0;
// //     MPI_Finalized(&mpi_is_finalized);
// //     if (!mpi_is_finalized)
// //       MPI_Finalize();
// // #endif
// //   }
// 
//   void SetNumThreads(int n)
//   {
//     if (n==0)
//       n=1;
//     else if (n <0) 
//       n = tbb::task_scheduler_init::default_num_threads();
//     num_threads = n;
//     tbbinit.terminate();
//     tbbinit.initialize(num_threads);
// #if defined( _OPENMP )
// #ifdef DEBUG
//     printf("omp num threads <- %i\n", num_threads);
// #endif
//     omp_set_num_threads(num_threads);
// #endif
//   }
// 
//   //members
//   int num_threads;
//   tbb::task_scheduler_init tbbinit;
// };
// 
// static std::auto_ptr<MPPrivate> mp;
// 
// bool MultiprocessingInitializer_exists()
// {
//   return mp.get();
// }
// 
// int GetNumThreads()
// {
//   if (!mp.get()) 
//     throw std::runtime_error("Multiprocessing not initialized! Construct a MultiprocessingInitializer!");
//   return mp->num_threads;
// }
// 
// int OmpGetCurrentThread()
// {
// #if defined( _OPENMP )
//     assert(omp_in_parallel() || (omp_get_num_threads() == 1 && omp_get_thread_num()==0));
//     return omp_get_thread_num();
// #else
//     return 0;
// #endif
// }
// 
// int OmpGetMaxThreadCount()
// {
// #if defined( _OPENMP )
//     return omp_get_max_threads();
// #else
//     return 1;
// #endif  
// };
// 
// 
// void SetNumThreads(int n)
// {
//   if (!mp.get()) 
//     throw std::runtime_error("Multiprocessing not initialized! Construct a MultiprocessingInitializer! from SetNumThreads");
//   mp->SetNumThreads(n);
// }
// 
// MultiprocessingInitializer::MultiprocessingInitializer(int argc, char **argv, int num_threads)
// {
//   if (mp.get()) 
//     throw std::runtime_error("A MultiprocessingInitializer lives already!");
//   mp.reset(new MPPrivate(argc, argv, num_threads));
// }
// MultiprocessingInitializer::MultiprocessingInitializer(int num_threads)
// {
//   if (mp.get()) 
//     throw std::runtime_error("A MultiprocessingInitializer lives already!");
//   mp.reset(new MPPrivate(num_threads));
// }
// 
// MultiprocessingInitializer::~MultiprocessingInitializer()
// {
//   mp.reset(NULL);
// }
// 
// void initMultithreading(int argc, char** argv, int num_threads)
// {
//   if (mp.get()) {
//     myAssert(false);
//     throw std::runtime_error("A MultiprocessingInitializer lives already!"); 
//   }
//   mp.reset(new MPPrivate(argc, argv, num_threads));
// }





Log& log()
{
  static Log l(std::cout);
  static bool init = false;
  if (!init)
  {
    std::streambuf *b_cout = std::cout.rdbuf();
    std::streambuf *b_log = l.rdbuf();
    std::cout.rdbuf(b_log);
    init = true;
  }
  return l;
}


static bool DefaultCheckAbort()
{
  return false;
};
AbortFunction checkAbort = DefaultCheckAbort;

}//end namespace my

SystemParameters::SystemParameters()
{
  num_threads = 1;
  cluster = "local";
  computing_node = "local";
  num_threads_queuing = 1;
  mem_in_GB = 1;
  
  isRerun = false;
  reRunNumber = 0;
  // Get the current time in local Timezone   
  boost::posix_time::ptime timeLocal = boost::posix_time::second_clock::local_time();
  date_of_run = boost::posix_time::to_simple_string(timeLocal);
}
void SystemParameters::assign(const boost::property_tree::ptree& pt)
{
  #define DOPT(name) boost::property_tree::get(name, #name, pt)
  DOPT(num_threads);
  DOPT(cluster);
  DOPT(computing_node);
  DOPT(num_threads_queuing);
  DOPT(mem_in_GB);
  
  DOPT(isRerun);
  DOPT(reRunNumber);
  DOPT(date_of_run);
  #undef DOPT
}
boost::property_tree::ptree SystemParameters::as_ptree() const
{
  boost::property_tree::ptree pt;
  #define DOPT(name) pt.put(#name, name)
  DOPT(cluster);
  DOPT(computing_node);
  DOPT(num_threads);
  DOPT(num_threads_queuing);
  DOPT(mem_in_GB);
  
  DOPT(isRerun);
  DOPT(reRunNumber);
  DOPT(date_of_run);
  DOPT(JobID);
  #undef DOPT
  return pt;
}

void readSystemParameters(SystemParameters &sysParamsToFill)
{
  std::cout << "reading System parameters" << std::endl;
  // update cluster information, if we are on a cluster
  if( std::getenv("SLURM_CLUSTER_NAME") )
  {
    //systemSettings.put("cluster", std::getenv("SLURM_CLUSTER_NAME"));
    sysParamsToFill.cluster = std::getenv("SLURM_CLUSTER_NAME");
  }
  if( std::getenv("SLURMD_NODENAME") )
  {
    //systemSettings.put("computing_node", std::getenv("SLURMD_NODENAME"));
    sysParamsToFill.computing_node = std::getenv("SLURMD_NODENAME");
  }
  else
  {
    sysParamsToFill.computing_node = boost::asio::ip::host_name();
//     if( std::getenv("HOSTNAME") ) //HOSTNAME is not cross plattform an not supported under arch linux
//     {
//       sysParamsToFill.computing_node = std::getenv("HOSTNAME");
//     }
  }
  if( std::getenv("SLURM_CPUS_ON_NODE") )
  {
    //systemSettings.put("num_threads_queuing", std::getenv("SLURM_CPUS_ON_NODE"));
    sysParamsToFill.num_threads_queuing = atoi( std::getenv("SLURM_CPUS_ON_NODE"));
  }
  //systemSettings.put("num_threads", omp_get_max_threads());
  sysParamsToFill.num_threads = omp_get_max_threads();
}
