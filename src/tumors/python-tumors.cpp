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
#include "../python_krebsutils/python-helpers.h"
#include "vesselgen/vesselgen.h"
#include "calcflow.h"
#include <boost/property_tree/info_parser.hpp>

#include <fenv.h>
#include "mwlib/hdf_wrapper_ptree.h"
#include "mwlib/log.h"
#include "shared-objects.h"
#include "common.h"

//#include <boost/python/module.hpp>
//#include <boost/python/def.hpp>
#include "faketum.h"

#include "bulktissue-no-vessels.h"
#include "bulktissue-with-vessels.h"

#define BOOST_RESULT_OF_USE_DECLTYPE 1

#ifdef MILOTTI_MTS
  #include "faketum_mts.h"
#endif

#ifdef USE_ADAPTION
  #include "../adaption/adaption_model2.h"
#endif

namespace py = boost::python;

namespace Tumors{
#ifdef MILOTTI_MTS
/** @brief fake tumor with multicellular tumor spheroids
 * here a growing sphere of tumor cells is assumed, no tumor model is used for
 * that. One needs a growing speed
 */
void run_fakeTumor_mts(const py::str &param_info_str)
{
  std::cout << "run_fakeTumor_mts on c++ called" << std::endl;
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
  //construct default parameters
  FakeTum::Parameters defaultParams;
  printPtree(defaultParams.as_ptree());
  ptree pt_params = convertInfoStr(param_info_str, defaultParams.as_ptree());
  printPtree(pt_params);
  FakeTumMTS::FakeTumorSimMTS s;
  int returnCode = s.run(pt_params);
  //return returnCode;
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
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
  //construct default parameters
  FakeTum::Parameters defaultParams;
  printPtree(defaultParams.as_ptree());
  ptree pt_params = convertInfoStr(param_info_str, defaultParams.as_ptree());
  printPtree(pt_params);
  FakeTum::FakeTumorSim s;
  int returnCode = s.run(pt_params);
  //return returnCode;
}
void export_faketum()
{
  py::def("run_faketum_", run_fakeTumor);
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
    O2Model::PrezO2Params prezO2params;
    //Adaption::Parameters adaption_params;
    
    ptree all_pt_params;
    all_pt_params = sparams.as_ptree();
    all_pt_params.put_child("tumor", pparams.as_ptree());
    all_pt_params.put_child("calcflow", bfparams.as_ptree());
    all_pt_params.put_child("vessels", vessel_params.as_ptree());
    all_pt_params.put_child("prez_o2", prezO2params.as_ptree());
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
    BulkTissueWithoutVessels::run(all_pt_params);
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
  theBulkTissueSim.run(pt_params);
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
  PyEval_InitThreads(); // need for release of the GIL (http://stackoverflow.com/questions/8009613/boost-python-not-supporting-parallelism)
//   if (my::MultiprocessingInitializer_exists())
//   {
//   }
//   else
//   {
//     my::initMultithreading(0, NULL, 1);
//   }
  my::checkAbort = PyCheckAbort; // since this is the python module, this is set to use the python signal check function
  Tumors::export_faketum();
  //Tumors::export_faketum_mts();
  //Tumors::export_bulktissue_no_vessels();
  //bla
  
// #ifdef MILOTTI_MTS
//   Tumors::export_bulktissue_with_vessels();
// #endif
  //bla
}
