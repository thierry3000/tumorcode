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
#include "../python_krebsutils/python_helpers.h"
#include "vesselgen/vesselgen.h"
#include "calcflow.h"
#include <boost/property_tree/info_parser.hpp>

#include <fenv.h>
#include <cstdlib> // std::getenv()
#include "mwlib/hdf_wrapper_ptree.h"
#include "mwlib/log.h"
#include "shared-objects.h"
#include "common.h"

#include <boost/filesystem.hpp>
//#include <boost/python/def.hpp>
#include "faketum.h"

#include "bulktissue-no-vessels.h"
#include "bulktissue-with-vessels.h"

#define BOOST_RESULT_OF_USE_DECLTYPE 1

#ifdef MILOTTI_MTS
  #include "faketum_mts.h"
  #include "../detailedO2/oxygen_model2.h"
#endif

#ifdef USE_ADAPTION
  #include "../adaption/adaption_model2.h"
#endif

#include <stdio.h>
#include <execinfo.h> //backtrace()
#include <signal.h>
#include <stdlib.h>
#include <unistd.h>
//https://stackoverflow.com/questions/77005/how-to-automatically-generate-a-stacktrace-when-my-gcc-c-program-crashes
void handler(int sig) {
  void *array[42];
  size_t size;

  // get void*'s for all entries on the stack
  size = backtrace(array, 42);

  // print out all the frames to stderr
  fprintf(stderr, "Error: signal %d:\n", sig);
  backtrace_symbols_fd(array, size, STDERR_FILENO);
  exit(sig);
}
void baz() {
 int *foo = (int*)-1; // make a bad pointer
  printf("%d\n", *foo);       // causes segfault
}

void bar() { baz(); }
void foo() { bar(); }
// #endif

namespace Tumors{
  
  
#ifdef MILOTTI_MTS
/** @brief fake tumor with multicellular tumor spheroids
 * here a growing sphere of tumor cells is assumed, no tumor model is used for
 * that. One needs a growing speed
 */
void run_fakeTumor_mts(const py::str &param_info_str)
{
  /* reading parameter string from info file */
  ptree pt_params = convertInfoStr(param_info_str, ptree());
  std::cout << "run_fakeTumor_mts on c++ called" << std::endl;
#ifdef DEBUG
  std::cout << "with params: " << std::endl;
  printPtree(pt_params);
#endif
  // enable standard exception handling
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
  
  /* construct default simulation */
  FakeTumMTS::FakeTumorSimMTS s;
  
  /* 
   * create ptree with default settings!!!
   */
  ptree detailedO2Settings = s.o2_params.as_ptree();
  ptree bfSettings = s.o2_sim.bfparams.as_ptree();
  ptree fakeTumMTSSettings = s.params.as_ptree();
  ptree systemSettings = s.mySystemParameters.as_ptree();
  ptree vesselSettings = s.vessel_model.params.as_ptree();

  #ifndef NDEBUG
  std::cout << "with detailed params: " << std::endl;
  printPtree(detailedO2Settings);
  std::cout << "with calcflow params: " << std::endl;
  printPtree(bfSettings);
  #endif
  // update settings with the read in data
  boost::property_tree::update(detailedO2Settings, pt_params.get_child("detailedo2"));
  boost::property_tree::update(bfSettings, detailedO2Settings.get_child("calcflow"));
  boost::property_tree::update(fakeTumMTSSettings, pt_params);
  boost::property_tree::update(vesselSettings, pt_params.get_child("vessels"));
  
  // update cluster information, if we are on a cluster
  if( std::getenv("SLURM_CLUSTER_NAME") )
  {
    systemSettings.put("cluster", std::getenv("SLURM_CLUSTER_NAME"));
  }
  if( std::getenv("SLURMD_NODENAME") )
  {
    systemSettings.put("computing_node", std::getenv("SLURMD_NODENAME"));
  }
  if( std::getenv("SLURM_CPUS_ON_NODE") )
  {
    systemSettings.put("num_threads_queuing", std::getenv("SLURM_CPUS_ON_NODE"));
  }
  systemSettings.put("num_threads", omp_get_max_threads());
  
  #ifdef DEBUG
  std::cout << "detailed params after update: " << std::endl;
  printPtree(detailedO2Settings);
  std::cout << "calcflow params after update: " << std::endl;
  printPtree(bfSettings);
  #endif
  
  // assign o2 parameters to the simulation
  s.o2_params.assign(detailedO2Settings);
  // assign bfparams parameters to the simulation
  s.o2_sim.bfparams.assign(bfSettings);
  s.bfparams.assign(bfSettings);
  // assign vessel parameters to the simulation
  s.vessel_model.params.assign(vesselSettings);
  s.params.assign(fakeTumMTSSettings);
  s.mySystemParameters.assign(systemSettings);
  /* 
   * if we are on a cluster, we expect multiple runs
   * and create a directory for each run 
   */
  boost::filesystem::path P = boost::filesystem::path(s.params.fn_out);
//   for(auto& part : boost::filesystem::path(s.params.fn_out))
//         std::cout << part << "\n";
  
  if( std::getenv("SLURM_JOB_ID") )// environmental variable present, we are on a slurm cluster queue!
  {
    boost::filesystem::path pathOfNewFolder = P.parent_path()/boost::filesystem::path(std::getenv("SLURM_JOB_ID"));
    boost::filesystem::create_directory(pathOfNewFolder);
    std::cout << pathOfNewFolder << std::endl;
    boost::filesystem::path newPath = pathOfNewFolder / boost::filesystem::path(P.stem().string());
    std::cout << newPath << std::endl;
    s.params.fn_out = newPath.string();
  }
  
    //ReleaseGIL unlock(); // allow the python interpreter to do things while this is running
    /*
     * run major simulation, hopefully all parameters are set correct at this stage
     */
#ifdef EPETRA_MPI
  std::cout << "EPETRA_MPI flag is set!\n" << std::endl;
  int mpi_is_initialized = 0;
  int prov;
  MPI_Initialized(&mpi_is_initialized);
  if (!mpi_is_initialized)
    //MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE,&prov);
    MPI_Init_thread(0, NULL, 1,&prov);
#endif
  try
  {

// #ifndef NDEBUG
  signal(SIGSEGV, handler);
  signal(SIGFPE, handler);
// #endif
  int returnCode = s.run();
  }
  catch(std::exception &ex)
  {
    std::cout << ex.what();
  }

  if (PyErr_Occurred() != NULL) return; // don't save stuff
}
void export_faketum_mts()
{
  py::def("run_faketum_mts_", run_fakeTumor_mts);
}
#endif

/** @brief fake tumor
 * here a growing sphere of tumor cells is assumed, no tumor model is used for
 * that. One needs a growing speed
 */
void run_fakeTumor(const py::str &param_info_str)
{
  ptree pt_params = convertInfoStr(param_info_str, ptree());
  std::cout << "run_fakeTumor on c++ called" << std::endl;
#ifdef DEBUG
  std::cout << "with params: " << std::endl;
  printPtree(pt_params);
#endif
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
  
  FakeTum::FakeTumorSim s;
  //construct default parameters
  //get default params
  ptree bfSettings = s.params.bfparams.as_ptree();
  ptree vesselSettings = s.vessel_model.params.as_ptree();
#ifdef USE_ADAPTION
  ptree adaptionSettings = s.params.adap_params.as_ptree();
  if(pt_params.count("adaption")>0)
  {
    boost::property_tree::update(adaptionSettings, pt_params.get_child("adaption"));
  }
  s.params.adap_params.assign(adaptionSettings);
#endif
  ptree fakeTumSettings = s.params.as_ptree();
  //update with read in params
  boost::property_tree::update(vesselSettings, pt_params.get_child("vessels"));
  boost::property_tree::update(bfSettings, pt_params.get_child("calcflow"));
  boost::property_tree::update(fakeTumSettings, pt_params);
  
  s.vessel_model.params.assign(vesselSettings);
  s.params.bfparams.assign(bfSettings);
  s.params.assign(fakeTumSettings);
  
  try
  {
#ifdef EPETRA_MPI
    std::cout << "EPETRA_MPI flag is set!\n" << std::endl;
    int mpi_is_initialized = 0;
    int prov;
    MPI_Initialized(&mpi_is_initialized);
    if (!mpi_is_initialized)
      //MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE,&prov);
      MPI_Init_thread(0, NULL, 1,&prov);
#endif
    int returnCode = s.run();
  }
  catch(std::exception &ex)
  {
    std::cout << ex.what();
  }
}

// void printObjectName(H5::H5Location&, H5std_string, void*)
// {
//   
// }
//typedef void(* 	attr_operator_t )(H5Object &loc, const H5std_string attr_name, void *operator_data)

//H5::attr_operator_t printObjectNameH5;
// void printObjectName( H5::H5Location &g, const H5std_string attr_name, void* operator_data)
// {
//   try
//   {
//     std::cout << attr_name << std::endl;
//   }
//   catch(H5::Exception e)
//   {
//     e.printError();
//   }
// }


void rerun_fakeTumor(const py::str &filename_of_previous_run)
{
  std::cout << "rerun_fakeTumor on c++ called" << std::endl;
  char const* fn_of_previous_sim_c_str = py::extract<char const*>(filename_of_previous_run);
#ifdef DEBUG
  std::cout << "with previous simulation file: " << std::endl;
  //std::printf("%s\n",filename_of_previous_run );
  
   // Print it using printf
  std::printf("%s\n", fn_of_previous_sim_c_str);
#endif
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
  FakeTum::FakeTumorSim s;
  //construct default parameters
  //get default params
  ptree bfSettings = s.params.bfparams.as_ptree();
  ptree vesselSettings = s.vessel_model.params.as_ptree();
  
#ifdef USE_ADAPTION
  ptree adaptionSettings = s.params.adap_params.as_ptree();
  if(pt_params.count("adaption")>0)
  {
    boost::property_tree::update(adaptionSettings, pt_params.get_child("adaption"));
  }
  s.params.adap_params.assign(adaptionSettings);
#endif
  ptree fakeTumSettings = s.params.as_ptree();
  //update with read in params
  H5::H5File file;
  H5::Group h5_params_of_previous_run;
  H5::Group h5_vessel_params_of_previous_run;
  H5::Group h5_calcflow_of_previous_run;
  try{
    file = H5::H5File(fn_of_previous_sim_c_str, H5F_ACC_RDONLY);
    h5_params_of_previous_run = file.openGroup("/parameters");
    h5_vessel_params_of_previous_run = h5_params_of_previous_run.openGroup("vessels");
    h5_calcflow_of_previous_run = h5_params_of_previous_run.openGroup("calcflow");
  }
  catch(H5::Exception e)
  {
    e.printErrorStack();
  }
  ReadHdfPtree(vesselSettings, h5_vessel_params_of_previous_run);
  ReadHdfPtree(bfSettings, h5_calcflow_of_previous_run);
  ReadHdfPtree(fakeTumSettings, h5_params_of_previous_run);
  
  //boost::property_tree::update(vesselSettings, pt_params.get_child("vessels"));
  //boost::property_tree::update(bfSettings, pt_params.get_child("calcflow"));
  //boost::property_tree::update(fakeTumSettings, pt_params);
  
  s.vessel_model.params.assign(vesselSettings);
  s.params.bfparams.assign(bfSettings);
  s.params.assign(fakeTumSettings);
  
  try
  {
#ifdef EPETRA_MPI
    std::cout << "EPETRA_MPI flag is set!\n" << std::endl;
    int mpi_is_initialized = 0;
    int prov;
    MPI_Initialized(&mpi_is_initialized);
    if (!mpi_is_initialized)
      //MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE,&prov);
      MPI_Init_thread(0, NULL, 1,&prov);
#endif
    int returnCode = s.run();
  }
  catch(std::exception &ex)
  {
    std::cout << ex.what();
  }
}
void export_faketum()
{
  py::def("run_faketum_", run_fakeTumor);
  py::def("rerun_faketum_", rerun_fakeTumor);
}

/** @brief BulkTissue no vessels
 */
void run_bulktissue_no_vessels(const py::str &param_info_str)
{
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
  /* Prameter Handling */
  // construct default parameters
  BulkTissueWithoutVessels::SimulationParameters sparams;
  NewBulkTissueModel::Params pparams;
  BloodFlowParameters bfparams;
  VesselModel1::Params vessel_params;
  O2Model::SimpleO2Params prezO2params;
  //Adaption::Parameters adaption_params;
  
  ptree all_pt_params;
  all_pt_params = sparams.as_ptree();
  all_pt_params.put_child("tumor", pparams.as_ptree());
  all_pt_params.put_child("calcflow", bfparams.as_ptree());
  all_pt_params.put_child("vessels", vessel_params.as_ptree());
  all_pt_params.put_child("simple_o2", prezO2params.as_ptree());
  //all_pt_params.put_child("adaption", adaption_params.as_ptree());
  cout.rdbuf(my::log().rdbuf());
  {
  //boost::optional<ptree> read_params = pt_params;
  //boost::optional<ptree> read_params = HandleSimulationProgramArguments(all_pt_params, argc, argv);
  //if (!read_params) 
  //  return 0;
  /** get the read params*/
  ptree pt_params = convertInfoStr(param_info_str, all_pt_params);
  BulkTissueWithoutVessels::SimulationParameters::update_ptree(all_pt_params, pt_params);
  
  all_pt_params.put<Int3>("lattice_size", Int3(200,1,1));
  sparams.assign(all_pt_params);
  //boost::property_tree::update(params, BulkTissueWithoutVessels::Params().as_ptree());
  
  all_pt_params.put_child("tumor", NewBulkTissueModel::Params().as_ptree());
  { 
#ifdef DEBUG
    cout << "read params in main are: ";
    boost::property_tree::write_info(cout, all_pt_params);
    cout << endl;
#endif
  }

  }//end cout.buffer
  /* start */
  try
  {
#ifdef EPETRA_MPI
    std::cout << "EPETRA_MPI flag is set!\n" << std::endl;
    int mpi_is_initialized = 0;
    int prov;
    MPI_Initialized(&mpi_is_initialized);
    if (!mpi_is_initialized)
      //MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE,&prov);
      MPI_Init_thread(0, NULL, 1,&prov);
#endif
    signal(SIGSEGV, handler);
    signal(SIGFPE, handler);
    BulkTissueWithoutVessels::run(all_pt_params);
  }
  catch(std::exception &ex)
  {
    std::cout << ex.what();
  }
}
void export_bulktissue_no_vessels()
{
  py::def("run_bulktissue_no_vessels_", run_bulktissue_no_vessels);
}

/** @brief BulkTissue with vessels
 * the complete stuff
 */
//static py::object
void run_bulktissue_with_vessels(const py::str &param_info_str)
{
  BulkTissue::Params someDefaults;
  const ptree neededDefaults = someDefaults.as_ptree();
  const ptree pt_params = convertInfoStr(param_info_str, neededDefaults);
  
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
   
  BulkTissue::NewTumorSim theBulkTissueSim;
  try{
#ifdef EPETRA_MPI
    std::cout << "EPETRA_MPI flag is set!\n" << std::endl;
    int mpi_is_initialized = 0;
    int prov;
    MPI_Initialized(&mpi_is_initialized);
    if (!mpi_is_initialized)
      //MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE,&prov);
      MPI_Init_thread(0, NULL, 1,&prov);
#endif
    signal(SIGSEGV, handler);
    signal(SIGFPE, handler);
    theBulkTissueSim.run(pt_params);
  }
  catch(std::exception &ex)
  {
    std::cout << ex.what();
  }
  catch(H5::Exception e)
  {
    e.printErrorStack();
  }
  //return 0;
}

void export_bulktissue_with_vessels()
{
  py::def("run_bulktissue_with_vessels_", run_bulktissue_with_vessels);
}

}//end namespace Tumors

#ifdef DEBUG
BOOST_PYTHON_MODULE(libtumors_d)
#else
BOOST_PYTHON_MODULE(libtumors_)
#endif
{ 
  Py_Initialize();
#if BOOST_VERSION>106300
  np::initialize();
#endif
  PyEval_InitThreads(); // need for release of the GIL (http://stackoverflow.com/questions/8009613/boost-python-not-supporting-parallelism)
//HACK2018
  //   if (my::MultiprocessingInitializer_exists())
//   {
//   }
//   else
//   {
//     my::initMultithreading(0, NULL, 1);
//   }
  my::checkAbort = PyCheckAbort; // since this is the python module, this is set to use the python signal check function
  Tumors::export_faketum();
  Tumors::export_bulktissue_no_vessels();
  Tumors::export_bulktissue_with_vessels();
  
#ifdef MILOTTI_MTS
  Tumors::export_faketum_mts();
#endif
}




template<class T>
inline boost::optional<T> getOptional(const char* name, py::dict &d)
{
  py::object o = d.get(name);
  if (!o.is_none())
    return boost::optional<T>(py::extract<T>(o));
  else
    return boost::optional<T>();
}

