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
  #ifdef USE_DETAILED_O2
    #include "../detailedO2/oxygen_model2.h"
  #endif
#endif

#ifdef USE_ADAPTION
  #include "../adaption/adaption_model2.h"
#endif


namespace Tumors{
  
void InitSimpleO2ParametersFromPy(O2Model::SimpleO2Params &params, py::dict py_parameters)
{
//   #define GET_SimpleO2_PARAM_FROM_DICT(TYPE, NAME) checkedExtractFromDict<TYPE>(py_parameters, NAME);
//   #define GET_SimpleO2_PARAM_IF_NONNONE(TARGET, TYPE, NAME) { py::object o(py_parameters.get(NAME)); if (!o.is_none()) TARGET=py::extract<TYPE>(o); }
//  params.hematocrit_init = GET_SimpleO2_PARAM_FROM_DICT(int, "hematocrit_init");
  params.hematocrit_init = checkedExtractFromDict<double>(py_parameters,"hematocrit_init");
  cout<<params.hematocrit_init<<endl;
//   params.sat_curve_exponent = GET_PO2_PARAM_FROM_DICT(double, "S_n");
//   params.sat_curve_p50 = GET_PO2_PARAM_FROM_DICT(double, "S_p50");
//   GET_PO2_PARAM_IF_NONNONE(params.axial_integration_step_factor, double, "axial_integration_step_factor");
//   GET_PO2_PARAM_IF_NONNONE(params.debug_zero_o2field, bool, "debug_zero_o2field");
//   GET_PO2_PARAM_IF_NONNONE(params.debug_fn, string, "debug_fn");
//   GET_PO2_PARAM_IF_NONNONE(params.transvascular_ring_size, double, "transvascular_ring_size");
//   double kd, rd_norm = 100., rd_tum = 100., rd_necro = 100.;
//   kd = GET_PO2_PARAM_FROM_DICT(double, "D_tissue");
//   params.tissue_solubility = GET_PO2_PARAM_FROM_DICT(double, "solubility_tissue");
//   params.plasma_solubility = GET_PO2_PARAM_FROM_DICT(double, "solubility_plasma");
//   
//   GET_PO2_PARAM_IF_NONNONE(rd_norm, double, "rd_norm");
//   GET_PO2_PARAM_IF_NONNONE(rd_tum, double, "rd_tum");
//   GET_PO2_PARAM_IF_NONNONE(rd_necro, double, "rd_necro");
//   //rd_tum = GET_PO2_PARAM_FROM_DICT(double, "rd_tum");
//   //rd_necro = GET_PO2_PARAM_FROM_DICT(double, "rd_necro");
//   params.SetTissueParamsByDiffusionRadius(kd, params.tissue_solubility, rd_norm, rd_tum, rd_necro);
//   params.po2init_r0 = GET_PO2_PARAM_FROM_DICT(double, "po2init_r0");
//   params.po2init_dr = GET_PO2_PARAM_FROM_DICT(double, "po2init_dr");
//   params.po2init_cutoff = GET_PO2_PARAM_FROM_DICT(double, "po2init_cutoff");
//   GET_PO2_PARAM_IF_NONNONE(params.michaelis_menten_uptake, bool, "michaelis_menten_uptake");
//   GET_PO2_PARAM_IF_NONNONE(params.po2_mmcons_k[DetailedPO2::NORMAL], double, "mmcons_k_norm");
//   GET_PO2_PARAM_IF_NONNONE(params.po2_mmcons_k[DetailedPO2::TUMOR], double, "mmcons_k_tum");
//   GET_PO2_PARAM_IF_NONNONE(params.po2_mmcons_k[DetailedPO2::NECRO], double, "mmcons_k_necro");
//   GET_PO2_PARAM_IF_NONNONE(params.po2_mmcons_m0[DetailedPO2::NORMAL], double, "mmcons_m0_norm");
//   GET_PO2_PARAM_IF_NONNONE(params.po2_mmcons_m0[DetailedPO2::TUMOR], double, "mmcons_m0_tum");
//   GET_PO2_PARAM_IF_NONNONE(params.po2_mmcons_m0[DetailedPO2::NECRO], double, "mmcons_m0_necro");
//   GET_PO2_PARAM_IF_NONNONE(params.haemoglobin_binding_capacity, double, "haemoglobin_binding_capacity");
//   string bc_type("neumann");
//   GET_PO2_PARAM_IF_NONNONE(bc_type, string, "tissue_po2_boundary_condition");
//   if (bc_type == "dirichlet_x") 
//     params.tissue_boundary_condition_flags = 1; //FiniteVolumeMatrixBuilder;
//   else if (bc_type == "dirichlet_yz") 
//     params.tissue_boundary_condition_flags = 2;
//   else if (bc_type == "dirichlet") 
//     params.tissue_boundary_condition_flags = 3;
//   else if (bc_type == "neumann") 
//     params.tissue_boundary_condition_flags = 0;
//   else 
//     throw std::invalid_argument("tissue_po2_boundary_condition must be 'dirichlet','dirichlet_x', 'dirichlet_yz' or 'neumann'");
//   GET_PO2_PARAM_IF_NONNONE(params.tissue_boundary_value, double, "tissue_boundary_value");
//   GET_PO2_PARAM_IF_NONNONE(params.extra_tissue_source_linear, double, "extra_tissue_source_linear");
//   GET_PO2_PARAM_IF_NONNONE(params.extra_tissue_source_const, double, "extra_tissue_source_const");
//   // required parameters for transvascular transport
//   params.conductivity_coeff1 = GET_PO2_PARAM_FROM_DICT(double, "conductivity_coeff1"); 
//   params.conductivity_coeff2 = GET_PO2_PARAM_FROM_DICT(double, "conductivity_coeff2");
//   params.conductivity_coeff3 = GET_PO2_PARAM_FROM_DICT(double, "conductivity_coeff3");
//   GET_PO2_PARAM_IF_NONNONE(params.convergence_tolerance, double, "convergence_tolerance");
//   GET_PO2_PARAM_IF_NONNONE(params.loglevel, int, "loglevel");
//   GET_PO2_PARAM_IF_NONNONE(params.approximateInsignificantTransvascularFlux, bool, "approximateInsignificantTransvascularFlux");
//   GET_PO2_PARAM_IF_NONNONE(params.D_plasma, double, "D_plasma");
//   GET_PO2_PARAM_IF_NONNONE(params.massTransferCoefficientModelNumber, int, "massTransferCoefficientModelNumber");
#undef GET_SimpleO2_PARAM_FROM_DICT
#undef GET_SimpleO2_PARAM_IF_NONNONE
//   params.UpdateInternalValues();
}  
  
  
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
  
#ifdef USE_DETAILED_O2
  /* 
   * detailedO2 case
   */
  //create ptree with default settings!!!
  ptree detailedO2Settings = s.o2_params.as_ptree();
  ptree bfSettings = s.o2_sim.bfparams.as_ptree();
  ptree fakeTumMTSSettings = s.params.as_ptree();
  #ifdef DEBUG
  std::cout << "with detailed params: " << std::endl;
  printPtree(detailedO2Settings);
  std::cout << "with calcflow params: " << std::endl;
  printPtree(bfSettings);
  #endif
  // update settings with the read in data
  boost::property_tree::update(detailedO2Settings, pt_params.get_child("detailedo2"));
  boost::property_tree::update(bfSettings, detailedO2Settings.get_child("calcflow"));
  boost::property_tree::update(fakeTumMTSSettings, pt_params);
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
  s.params.assign(fakeTumMTSSettings);
  std::string slurm_id = std::getenv("SLURM_JOB_ID");
  boost::filesystem::path P = boost::filesystem::path(s.params.fn_out);
  for(auto& part : boost::filesystem::path(s.params.fn_out))
        std::cout << part << "\n";
  

  const std::string rndString = "quz";
  boost::filesystem::path pathOfNewFolder = P.parent_path()/boost::filesystem::path(slurm_id);
  boost::filesystem::create_directory(pathOfNewFolder);
  std::cout << pathOfNewFolder << std::endl;
  boost::filesystem::path newPath = pathOfNewFolder / boost::filesystem::path(P.stem().string());
  std::cout << newPath << std::endl;

//   boost::filesystem::path p(s.params.fn_out);
//   std::string new_filename = p.leaf() + ".foo";
//   p.remove_leaf() /= new_filename;
//   std::cout << p << '\n';
  s.params.fn_out = newPath.string();
#else
  /*
   * simple or no o2 case
   * 
   */
  //get default params
  ptree detailedO2Settings = s.o2_params.as_ptree();
  ptree bfSettings = s.bfparams.as_ptree();
  ptree fakeTumMTSSettings = s.params.as_ptree();
  //update with read in params
  boost::property_tree::update(detailedO2Settings, pt_params.get_child("detailedo2"));
  boost::property_tree::update(bfSettings, pt_params.get_child("calcflow"));
  boost::property_tree::update(fakeTumMTSSettings, pt_params);
  //update read default parameters with read in parameters
  //boost::property_tree::update(destination, source);
  //O2Model::SimpleO2Params simpleO2params;
  ptree simpleO2Settings = s.o2_params.as_ptree();
  //boost::property_tree::update(simpleO2Settings, pt_params.get_child("simple_o2"));
  // assign params
  s.o2_params.assign(simpleO2Settings);
  s.bfparams.assign(bfSettings);
  s.params.assign(fakeTumMTSSettings);
  
#endif
  { 
    //ReleaseGIL unlock(); // allow the python interpreter to do things while this is running
    /*
     * run major simulation, hopefully all parameters are set correct at this stage
     */
    int returnCode = s.run();
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
  boost::property_tree::update(bfSettings, pt_params.get_child("calcflow"));
  boost::property_tree::update(fakeTumSettings, pt_params);
  
  s.params.bfparams.assign(bfSettings);
  s.params.assign(fakeTumSettings);
  
//   FakeTum::Parameters defaultParams;  /// default parameters
//   ptree faketumSettings = defaultParams.as_ptree();
//   boost::property_tree::update(faketumSettings,pt_params);
  //printPtree(defaultParams.as_ptree());
  
  //printPtree(faketumSettings);
  
  int returnCode = s.run();
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
  Py_Initialize();
#if BOOST_VERSION>106300
  np::initialize();
#endif
  PyEval_InitThreads(); // need for release of the GIL (http://stackoverflow.com/questions/8009613/boost-python-not-supporting-parallelism)
  if (my::MultiprocessingInitializer_exists())
  {
  }
  else
  {
    my::initMultithreading(0, NULL, 1);
  }
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

