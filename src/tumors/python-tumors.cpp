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
#if BOOST_VERSION>106300
  #include <boost/convert/lexical_cast.hpp>
#else
  #include <boost/lexical_cast.hpp>
#endif
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
//#ifdef MILOTTI_MTS
//H5::H5File file;
//#endif
//https://stackoverflow.com/questions/77005/how-to-automatically-generate-a-stacktrace-when-my-gcc-c-program-crashes
void handler(int sig) 
{
  //printf("%s", file.getFileName());
  //std::cout << "handler close of file: " << std::endl;
  //file.flush(H5F_SCOPE_GLOBAL);
  //file.close();
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
//H5::H5File file;
void run_fakeTumor_mts(const py::str &param_info_str_or_filename_of_pr, bool isRerun)
{
  //******GLOBALS
  ptree pt_params;
  char const* fn_of_previous_sim_c_str;
  if( !isRerun )
  {
    /* reading parameter string from info file */
    pt_params = convertInfoStr(param_info_str_or_filename_of_pr, ptree());
    std::cout << "run_fakeTumor_mts without rerun called on c++ side" << std::endl;
    #ifdef DEBUG
      std::cout << "with params: " << std::endl;
      printPtree(pt_params);
    #endif
  }
  else
  {
    std::cout << "run_fakeTumor_mts rerun called on c++ side" << std::endl;
    fn_of_previous_sim_c_str = py::extract<char const*>(param_info_str_or_filename_of_pr);
    //#ifdef DEBUG
      std::cout << "with filename: " << std::endl;
      std::printf("%s\n", fn_of_previous_sim_c_str);
    //#endif
  }
  // enable standard exception handling
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
  
  /* construct default simulation */
  FakeTumMTS::FakeTumorSimMTS s;
  //initialize cell system 
  //use memory on heap to not mess up allocation
  s.tumorcode_pointer_to_currentCellsSystem = new vbl::CellsSystem();
  
  /* 
   * create ptree with default settings!!!
   */
  ptree detailedO2Settings = s.o2_sim.params.as_ptree();
  ptree bfSettings = s.o2_sim.bfparams.as_ptree();
  ptree fakeTumMTSSettings = s.params.as_ptree();
  ptree systemSettings = s.mySystemParameters.as_ptree();
  ptree vesselSettings = s.vessel_model.params.as_ptree();
  ptree vblSettings;

  #ifndef NDEBUG
  std::cout << "with detailed params: " << std::endl;
  printPtree(detailedO2Settings);
  std::cout << "with calcflow params: " << std::endl;
  printPtree(bfSettings);
  std::cout << "with fakeTumMTS params: " << std::endl;
  printPtree(fakeTumMTSSettings);
  #endif
  
  if(const char* slurmJobID = std::getenv("SLURM_JOB_ID"))
  {
    s.mySystemParameters.JobID = boost::lexical_cast<int>(std::getenv("SLURM_JOB_ID"));
  }
  else 
  {
    s.mySystemParameters.JobID = 0;
  }
  if( !isRerun )
  {
    // update settings with the read in data
    boost::property_tree::update(detailedO2Settings, pt_params.get_child("detailedo2"));
    boost::property_tree::update(vesselSettings, pt_params.get_child("vessels"));
    boost::property_tree::update(bfSettings, detailedO2Settings.get_child("calcflow"));
    boost::property_tree::update(fakeTumMTSSettings, pt_params);
    fakeTumMTSSettings.put("vessel_path", "vessels");
    fakeTumMTSSettings.put("isRerun", 0);
    s.mySystemParameters.reRunNumber = 0;
    s.mySystemParameters.isRerun = false;
    
    readSystemParameters(s.mySystemParameters);
    
    boost::property_tree::update(systemSettings, s.mySystemParameters.as_ptree());
    #ifndef NDEBUG
      std::cout << "detailed params after update: " << std::endl;
      printPtree(detailedO2Settings);
      std::cout << "calcflow params after update: " << std::endl;
      printPtree(bfSettings);
    #endif
  }
  else
  {
    //update with params of previous run
    H5::Group h5_params_of_previous_run;
    H5::Group h5_vessel_params_of_previous_run;
    H5::Group h5_calcflow_of_previous_run;
    H5::Group h5_o2_params_of_previous_run;
    H5::Group h5_system_of_first_run;
    H5::Group h5_system_of_current_run;
    H5::Group last_state;
    H5::Group h5_vbl_param;
    H5::Group h5_vbl_Environment;
    H5::Group h5_vbl_Environment_0;
    H5::Group h5_vbl_dose_rateSignal;
    H5::Group h5_vbl_flowSignal;
    
    H5::Group h5_params_of_permanent;
    H5::Group h5_system_of_permanent;
    
    int reRunNumber;
    try
    {
      //file = H5::H5File(fn_of_previous_sim_c_str, H5F_ACC_RDONLY);
      //storing current system information needs write permissions
      std::unique_ptr<H5::H5File> file = std::unique_ptr<H5::H5File>(new H5::H5File(fn_of_previous_sim_c_str, H5F_ACC_RDWR));
      std::unique_ptr<H5::H5File> last_state_file = std::unique_ptr<H5::H5File>(new H5::H5File("last_state.h5", H5F_ACC_RDONLY));
    
      //file = H5::H5File(fn_of_previous_sim_c_str, H5F_ACC_RDWR);
      last_state = last_state_file->openGroup("/last_state");
      h5_params_of_previous_run = last_state_file->openGroup("/parameters");
      h5_params_of_permanent = file->openGroup("/parameters");
      
      h5_vessel_params_of_previous_run = h5_params_of_previous_run.openGroup("vessels");
      h5_calcflow_of_previous_run = h5_params_of_previous_run.openGroup("calcflow");
      h5_system_of_first_run = h5_params_of_previous_run.openGroup("system");
      h5_o2_params_of_previous_run = h5_params_of_previous_run.openGroup("o2_sim");
      h5_vbl_param = last_state.openGroup("vbl");
      h5_vbl_Environment =h5_vbl_param.openGroup("Environment");
      h5_vbl_Environment_0=h5_vbl_param.openGroup("Environment_0");
      h5_vbl_dose_rateSignal=h5_vbl_param.openGroup("dose_rateSignal");
      h5_vbl_flowSignal=h5_vbl_param.openGroup("flowSignal");
    }
    catch(H5::Exception &e)
    {
      cout << "Error opening the parameters" << endl;
      e.printErrorStack();
      cout << "emergency close called I hope unique_ptr close the file handler" << endl;
    }
    catch(std::exception &ex)
    {
      cout << "emergency close called 2" << endl;
    }
    ReadHdfPtree(vesselSettings, h5_vessel_params_of_previous_run);
    ReadHdfPtree(bfSettings, h5_calcflow_of_previous_run);
    ReadHdfPtree(fakeTumMTSSettings, h5_params_of_previous_run);
    //I do not need the system settings of the previous run, do I?
    //ReadHdfPtree(systemSettings, h5_system_of_first_run);
    
    ReadHdfPtree(detailedO2Settings, h5_o2_params_of_previous_run);
    
    ReadHdfPtree(vblSettings, h5_vbl_param);
    ptree vblSettings_Environment;
    ReadHdfPtree(vblSettings_Environment, h5_vbl_Environment);
    vblSettings.put_child("Environment", vblSettings_Environment);
#ifndef NDEBUG
    printPtree(vblSettings);
#endif
    ptree vblSettings_Environment_0;
    ReadHdfPtree(vblSettings_Environment_0, h5_vbl_Environment_0);
    vblSettings.put_child("Environment_0", vblSettings_Environment_0);
    ptree vblSettings_dose_rateSignal;
    ReadHdfPtree(vblSettings_dose_rateSignal, h5_vbl_dose_rateSignal);
    vblSettings.put_child("dose_rateSignal", vblSettings_dose_rateSignal);
    ptree vblSettings_flowSignal;
    ReadHdfPtree(vblSettings_flowSignal, h5_vbl_flowSignal);
    vblSettings.put_child("flowSignal", vblSettings_flowSignal);
    
    //ready types
    int ntypes = vblSettings.get<int>("ntypes");
    for(int i=0;i<ntypes; i++)
    {
      boost::format my_string_template("type_%i");
      ptree readOutType;
      H5::Group h5_vbl_type;
      string current_type_name = boost::str(my_string_template %  i);
      h5_vbl_type = h5_vbl_param.openGroup(current_type_name);
      ReadHdfPtree(readOutType, h5_vbl_type);
      h5_vbl_type.close();
      vblSettings.put_child(current_type_name, readOutType);
    }
    s.all_pt_params.put_child("vbl", vblSettings);
    try 
    {
      #ifndef NDEBUG
      cout << "**************** VBL Settings ***************"<< endl;
    printPtree(vblSettings);
#endif
      //s.tumorcode_pointer_to_currentCellsSystem->get_params_pointer()->assign(vblSettings);
      s.tumorcode_pointer_to_currentCellsSystem->assign(vblSettings);
    }
    catch(std::runtime_error &e)
    {
      cout << "error assigning vblSettings" << endl;
      cout << e.what() << endl;
      cout << "emergency close called 3" << endl;
    }
    //s.tumorcode_pointer_to_currentCellsSystem->assign(vblSettings);
    //s.set_CellTypeFromIndexVector();
    //s.get_CellTypeIndexVector();
    
    s.readVBLDataFromHDF(h5_vbl_param);
    
#ifndef NDEBUG
      std::cout << "found vbl Settings: " << std::endl;
      printPtree(vblSettings);
#endif
 
    //override read in
    fakeTumMTSSettings.put("vessel_path", "last_state/vessels");
    fakeTumMTSSettings.put("fn_vessel", "last_state.h5");
    systemSettings.put("isRerun", true);
    
    readAttrFromH5(last_state, "CURRENT_RERUN_NUMBER", reRunNumber);
    reRunNumber++;
    systemSettings.put("reRunNumber", reRunNumber);
    try
    {
      h5_system_of_permanent = h5_params_of_permanent.createGroup(str(format("system_rerun_%02i") % reRunNumber));
      //h5_system_of_current_run = h5_params_of_previous_run.createGroup(str(format("system_rerun_%02i") % reRunNumber));
    }
    catch(H5::Exception &e)
    {
      cout << "error: creating system_rerun parameters group" << endl;
      e.printErrorStack();
    }
    //WriteHdfPtree(h5_system_of_current_run, systemSettings);
    WriteHdfPtree(h5_system_of_permanent, systemSettings);
    //run time variables NO parameters
    readAttrFromH5(last_state, string("OUTPUT_NUM"), s.output_num);
    readAttrFromH5(last_state, string("NUM_ITERATION"), s.num_iteration);
    s.num_iteration++;
    readAttrFromH5(last_state, string("time"), s.time);
    readAttrFromH5(last_state, string("NEXT_OUTPUT_TIME"), s.next_output_time);
    readAttrFromH5(last_state, string("NEXT_ADAPTION_TIME"), s.next_adaption_time);
    
    h5_params_of_previous_run.close();
    h5_vessel_params_of_previous_run.close();
    h5_calcflow_of_previous_run.close();
    h5_system_of_first_run.close();
    h5_system_of_current_run.close();
    last_state.close();
    h5_vbl_param.close();
    h5_vbl_Environment.close();
    h5_vbl_Environment_0.close();
    h5_vbl_dose_rateSignal.close();
    h5_vbl_flowSignal.close();
    h5_system_of_permanent.close();
    h5_params_of_permanent.close();
  }
  
  // assign o2 parameters to the simulation
  
  try 
  {
    s.vessel_model.params.assign(vesselSettings);
    s.bfparams.assign(bfSettings);
    s.mySystemParameters.assign(systemSettings);
    s.params.assign(fakeTumMTSSettings);
  }
  catch(std::runtime_error &e)
  {
    cout << "error assigning fakeTumMTSSettings" << endl;
    cout << e.what() << endl;
  }
  
  s.o2_sim.bfparams.assign(bfSettings);
  s.o2_sim.params.assign(detailedO2Settings);
  

  
  if( !isRerun )
  {
    /* 
    * if we are on a cluster, we expect multiple runs
    * and create a directory for each run 
    */
    /*
     * this was not so nice for the reruns and the storage of last_state.h5
     */
    /*
    boost::filesystem::path P = boost::filesystem::path(s.params.fn_out);
    if( std::getenv("SLURM_JOB_ID") and !s.mySystemParameters.isRerun )// environmental variable present, we are on a slurm cluster queue!
    {
      boost::filesystem::path pathOfNewFolder = P.parent_path()/boost::filesystem::path(std::getenv("SLURM_JOB_ID"));
      boost::filesystem::create_directory(pathOfNewFolder);
      std::cout << pathOfNewFolder << std::endl;
      boost::filesystem::path newPath = pathOfNewFolder / boost::filesystem::path(P.stem().string()+".h5");
      std::cout << newPath << std::endl;
      s.params.fn_out = newPath.string();
    }
    */
  }

#ifdef EPETRA_MPI
  std::cout << "EPETRA_MPI flag is set!\n" << std::endl;
  int mpi_is_initialized = 0;
  int prov;
  MPI_Initialized(&mpi_is_initialized);
  if (!mpi_is_initialized)
    //MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE,&prov);
    MPI_Init_thread(0, NULL, 1,&prov);
#endif
  /**
   * run major simulation, hopefully all parameters are set correct at this stage
   */
  try
  {
    /**
     * output error stack on segfault and floating point error
     */
// #ifndef NDEBUG
    signal(SIGSEGV, handler);
    signal(SIGFPE, handler);
    signal(SIGTERM, handler);
    signal(SIGKILL, handler);
// #endif
    int returnCode = s.run();
    delete s.tumorcode_pointer_to_currentCellsSystem;
  }
  catch(std::exception &ex)
  {
    std::cout << "either SIGFPE or SIGFPE occured" << std::endl;
    std::cout << ex.what();
    std::cout << "emergency close of hdf5 file!!!!" << std::endl;
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
void run_fakeTumor(const py::str &param_info_str_or_filename_of_pr, bool isRerun)
{
  ptree pt_params;
  char const* fn_of_previous_sim_c_str;
  if( !isRerun )
  {
    /* reading parameter string from info file */
    pt_params = convertInfoStr(param_info_str_or_filename_of_pr, ptree());
    std::cout << "run_fakeTumor without rerun called on c++ side" << std::endl;
    #ifdef DEBUG
      std::cout << "with params: " << std::endl;
      printPtree(pt_params);
    #endif
  }
  else
  {
    std::cout << "run_fakeTumor rerun called on c++ side" << std::endl;
    fn_of_previous_sim_c_str = py::extract<char const*>(param_info_str_or_filename_of_pr);
    #ifdef DEBUG
      std::cout << "with filename: " << std::endl;
      std::printf("%s\n", fn_of_previous_sim_c_str);
    #endif
  }
  // enable standard exception handling
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
  
  FakeTum::FakeTumorSim s;
  
  /* 
   * create ptree with default settings!!!
   */
  ptree bfSettings = s.params.bfparams.as_ptree();
  ptree fakeTumSettings = s.params.as_ptree();
  ptree systemSettings = s.mySystemParameters.as_ptree();
  ptree vesselSettings = s.vessel_model.params.as_ptree();
#ifdef USE_ADAPTION
  ptree adaptionSettings = s.params.adap_params.as_ptree();
  if(pt_params.count("adaption")>0)
  {
    boost::property_tree::update(adaptionSettings, pt_params.get_child("adaption"));
  }
  s.params.adap_params.assign(adaptionSettings);
#endif
  
  if( !isRerun )
  {
    //update with read in params
    boost::property_tree::update(vesselSettings, pt_params.get_child("vessels"));
    boost::property_tree::update(bfSettings, pt_params.get_child("calcflow"));
    boost::property_tree::update(fakeTumSettings, pt_params);
    fakeTumSettings.put("vessel_path", "vessels");
    fakeTumSettings.put("isRerun", 0);
    #ifndef NDEBUG
      std::cout << "with calcflow params: " << std::endl;
      printPtree(bfSettings);
    #endif
    s.mySystemParameters.reRunNumber = 0;
    s.mySystemParameters.isRerun = false;
    
    readSystemParameters(s.mySystemParameters);
    
    boost::property_tree::update(systemSettings, s.mySystemParameters.as_ptree());
  }
  else
  {
    //update with params of previous run
    //H5::H5File file;
    std::unique_ptr<H5::H5File> file;
    std::unique_ptr<H5::H5File> last_state_file;
    H5::Group h5_params_of_permanent;
    H5::Group h5_system_of_permanent;
    
    H5::Group h5_params_of_previous_run;
    H5::Group h5_vessel_params_of_previous_run;
    H5::Group h5_calcflow_of_previous_run;
    H5::Group h5_system_of_first_run;
    H5::Group h5_system_of_current_run;
    H5::Group last_state;
    int reRunNumber;
    try{
      //file = H5::H5File(fn_of_previous_sim_c_str, H5F_ACC_RDONLY);
      //storing current system information needs write permissions
      file = std::unique_ptr<H5::H5File>(new H5::H5File(fn_of_previous_sim_c_str, H5F_ACC_RDWR));
      last_state_file = std::unique_ptr<H5::H5File>(new H5::H5File(fn_of_previous_sim_c_str, H5F_ACC_RDONLY));
      h5_params_of_permanent = file->openGroup("/parameters");
      //file = H5::H5File(fn_of_previous_sim_c_str, H5F_ACC_RDWR);
      h5_params_of_previous_run = last_state_file->openGroup("/parameters");
      h5_vessel_params_of_previous_run = h5_params_of_previous_run.openGroup("vessels");
      h5_calcflow_of_previous_run = h5_params_of_previous_run.openGroup("calcflow");
      h5_system_of_first_run = h5_params_of_previous_run.openGroup("system");
    }
    catch(H5::Exception &e)
    {
      e.printErrorStack();
    }
    ReadHdfPtree(vesselSettings, h5_vessel_params_of_previous_run);
    ReadHdfPtree(bfSettings, h5_calcflow_of_previous_run);
    ReadHdfPtree(fakeTumSettings, h5_params_of_previous_run);
    ReadHdfPtree(systemSettings, h5_system_of_first_run);
 
    //override read in
    fakeTumSettings.put("vessel_path", "last_state/vessels");
    fakeTumSettings.put("fn_vessel", fn_of_previous_sim_c_str);
    systemSettings.put("isRerun", true);
    try
    {
      last_state = last_state_file->openGroup("/last_state");
      readAttrFromH5(last_state, "CURRENT_RERUN_NUMBER", reRunNumber);
      reRunNumber++;
      systemSettings.put("reRunNumber", reRunNumber);
      
      h5_system_of_current_run = h5_params_of_previous_run.createGroup(str(format("system_rerun_%02i") % reRunNumber));
      h5_system_of_permanent = h5_params_of_permanent.createGroup(str(format("system_rerun_%02i") % reRunNumber)); 
      WriteHdfPtree(h5_system_of_current_run, systemSettings);
      WriteHdfPtree(h5_system_of_permanent, systemSettings);
    }
    catch(H5::Exception &e)
    {
      e.printErrorStack();
    }
    //run time variables NO parameters
    readAttrFromH5(last_state, string("OUTPUT_NUM"), s.output_num);
    readAttrFromH5(last_state, string("NUM_ITERATION"), s.num_iteration);
    readAttrFromH5(last_state, string("time"), s.time);
    readAttrFromH5(last_state, string("NEXT_OUTPUT_TIME"), s.next_output_time);
    readAttrFromH5(last_state, string("NEXT_ADAPTION_TIME"), s.next_adaption_time);
    
    h5_params_of_previous_run.close();
    h5_vessel_params_of_previous_run.close();
    h5_calcflow_of_previous_run.close();
    h5_system_of_first_run.close();
    h5_system_of_current_run.close();
    last_state.close();
  }
  
  // assign o2 parameters to the simulation
  s.vessel_model.params.assign(vesselSettings);
  s.params.bfparams.assign(bfSettings);
  s.params.assign(fakeTumSettings);
  s.mySystemParameters.assign(systemSettings);
  
  if( !isRerun )
  {
    /* 
    * if we are on a cluster, we expect multiple runs
    * and create a directory for each run 
    */
    boost::filesystem::path P = boost::filesystem::path(s.params.fn_out);
    if( std::getenv("SLURM_JOB_ID") and !s.mySystemParameters.isRerun )// environmental variable present, we are on a slurm cluster queue!
    {
      boost::filesystem::path pathOfNewFolder = P.parent_path()/boost::filesystem::path(std::getenv("SLURM_JOB_ID"));
      boost::filesystem::create_directory(pathOfNewFolder);
      std::cout << pathOfNewFolder << std::endl;
      boost::filesystem::path newPath = pathOfNewFolder / boost::filesystem::path(P.stem().string()+".h5");
      std::cout << newPath << std::endl;
      s.params.fn_out = newPath.string();
    }
  }
  
#ifdef EPETRA_MPI
    std::cout << "EPETRA_MPI flag is set!\n" << std::endl;
    int mpi_is_initialized = 0;
    int prov;
    MPI_Initialized(&mpi_is_initialized);
    if (!mpi_is_initialized)
      //MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE,&prov);
      MPI_Init_thread(0, NULL, 1,&prov);
#endif
  /**
   * run major simulation, hopefully all parameters are set correct at this stage
   */
  try
  {
    /**
     * output error stack on segfault and floating point error
     */
// #ifndef NDEBUG
    signal(SIGSEGV, handler);
    signal(SIGFPE, handler);
    signal(SIGKILL, handler);
    signal(SIGTERM, handler);
// #endif
    int returnCode = s.run();
  }
  catch(std::exception &ex)
  {
    std::cout << ex.what();
  }
  
  if (PyErr_Occurred() != NULL) return; // don't save stuff
}

void export_faketum()
{
  py::def("run_faketum_", run_fakeTumor);
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
//   catch()
//   {
//     e.printError();
//   }
// }


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
  catch(H5::Exception &e)
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

