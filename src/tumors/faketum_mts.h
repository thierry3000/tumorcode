/**
This file is part of tumorcode project.
(http://www.uni-saarland.de/fak7/rieger/homepage/research/tumor/tumor.html)

Copyright (C) 2017 Thierry Fredrich

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
#ifndef _FAKETUMMTS_H_
#define _FAKETUMMTS_H_

#include "common/common.h"
#include "common/calcflow.h"
#include "common/shared-objects.h"
#include <math.h> // calculte 3rd root
// to check if parameter files are present
#include <boost/filesystem.hpp>
#include <memory> //std::shared_ptr
#include "common/vesselmodel1.h"
#include "common/cell_based_oxygen_update_model.h"
#include "common/growthfactor_model.h"
#include "common/hdfio.h"

#include "mwlib/log.h"


#include "../common/glucose_model.h"

#include "trilinos_linsys_construction.h"

#include <ANN/ANN.h>
#include <vbl.h>
#include "../detailedO2/oxygen_model2.h"


#ifdef W_timing
  #include <chrono>
#endif


namespace FakeTumMTS{
  
struct State : boost::noncopyable
{
  State() :
    tumor_checksum(0), vessels_checksum(0), chem_checksum(0)
    {}
    
  //NewBulkTissueModel::State tumor;
  int tumor_checksum;
  //boost::scoped_ptr<VesselList3d> vessels;
  int vessels_checksum;

  boost::optional<Array3df> previous_po2field;
  boost::optional<DetailedPO2::VesselPO2Storage> previous_po2vessels;

  //Array3d<float> o2field;
  //Array3d<float> glucoseField;
  Array3d<float> gffield;
  Array3d<float> cell_O2_consumption;
  int chem_checksum;
};  


enum TissuePressureDistribution
{
  TISSUE_PRESSURE_SPHERE = 0,
  TISSUE_PRESSURE_SHELL = 1,
};
struct nearest
{
  double distance = std::numeric_limits<double>::max();
  uint indexOfVessel = 0;
};

struct Parameters
{
  double out_intervall, tend;
  double dt;
  double latest_executed_timepoint;
  
  double apply_adaption_intervall;
  string message;
  string fn_out, fn_vessel, vessel_path, vesselfile_message;
  string paramset_name;
  int vesselfile_ensemble_index;
  double rGf, gf_production_threshold;
  double rO2Consumtion;
  double tumor_radius, tumor_speed;
  double stopping_radius_fraction; // actual-stopping-radius = 0.5*system-length-along-the-longest-dimension * stopping_radius_fraction
  int tissuePressureDistribution;
  double tissuePressureWidth;
  double tissuePressureCenterFraction;
  //for continuum oxygen calculation
  Int3 lattice_size = {20,20,20};
  double lattice_scale = 10;
  bool useConstO2;
  bool useTumorcodeVessels;
  int output_num;
  double time;
  
  Parameters();
  void assign(const ptree &pt);
  ptree as_ptree() const;
  static void update_ptree(ptree &dst, const ptree &src);
};

//#ifdef W_timing
//everything in mu seconds
struct Timing
{
  std::chrono::steady_clock::time_point begin_init_o2;
  std::chrono::steady_clock::time_point end_init_o2;
  std::chrono::steady_clock::time_point begin_o2;
  std::chrono::steady_clock::time_point end_o2;
  std::chrono::steady_clock::time_point begin_ann;
  std::chrono::steady_clock::time_point end_ann;
  std::chrono::steady_clock::time_point begin_doStep; 
  std::chrono::steady_clock::time_point end_doStep;
  std::chrono::steady_clock::time_point begin_doMilottiStep; 
  std::chrono::steady_clock::time_point end_doMilottiStep;
  std::chrono::steady_clock::time_point begin_findNearestVessel;
  std::chrono::steady_clock::time_point end_findNearestVessel;
  std::chrono::steady_clock::time_point begin_mts_main_loop;
  std::chrono::steady_clock::time_point end_mts_main_loop;
  std::chrono::duration<double> time_diff;
  double run_init_o2 = 0;
  double run_o2 = 0;
  double run_ann = 0;
  double run_doStep = 0;
  double run_doMilottiStep = 0;
  double run_findNearestVessel = 0;
  double run_mts_main_loop = 0;

  void reset()
  {
    run_init_o2 = 0;
    run_o2 = 0;
    run_ann = 0;
    run_doStep =0;
    run_doMilottiStep =0 ;
    run_findNearestVessel=0;
    run_mts_main_loop=0;
  };
};
//#endif

struct FakeTumorSimMTS : public boost::noncopyable
{
  Timing currentTiming;
  std::shared_ptr<VesselList3d> vl;
  
  /** global ANN stuff
   */
const int ANN_dim = 3;            // space dimension
  const int ANN_maxPts = 25000;   // maximum number of data points --> to limit memory allocation
  const double ANN_eps = 0.0;     // error bound
  const int	   ANN_k= 5;          // number of nearest neighbors
  
  DynArray<nearest> vectorOfnearestVessels;
	
	
  VesselModel1::Model vessel_model;
  DetailedPO2::DetailedPO2Sim o2_sim;
  
  void calcChemFields();
  // lattice definition of the continuum field lattice
  ContinuumGrid grid;
  DomainDecomposition mtboxes;
  // the state
  State state;
  
//   void insertGlucoseCoefficients(int box_index, const BBox3& bbox, const State &state, FiniteVolumeMatrixBuilder& mb);
  //void insertGFCoefficients(int box_index, const BBox3& bbox, const State &state, FiniteVolumeMatrixBuilder& mb);
  Array3d<float> vessel_volume_fraction;
//   Array3d<float> vessel_o2src_clin, vessel_o2src_crhs;
//   Array3d<float> vessel_glucosesrc_clin, vessel_glucosesrc_crhs;
  Array3d<float> cell_GFsrc;
  Array3d<float> cell_O2src;
  
  Array3dOps<float> oxyops;
  Array3dOps<float> glucoseOps;
  double last_chem_update;
  int last_vessels_checksum;
  void UpdateVesselVolumeFraction();
  
//   GfModel gf_model;
  TissuePhases phases;//Declaration
  
  FakeTumMTS::Parameters params;
  SystemParameters mySystemParameters;
  //GlucoseModel::GlucoseParams glucoseParams;
  GfModel_Cell gf_model;
  CellBasedO2Uptake o2_uptake_model;
  
  BloodFlowParameters bfparams;

  ptree all_pt_params;

  double tumor_radius;
  double time;
  double next_output_time;
  double next_adaption_time;
  int num_iteration;
  int output_num;

  /**
   * interface functions for VBL
   */
  float estimateTumorRadiusFromCells();
  
  //growth factor
  float getGf(const Float3 &pos);
  float getPress(const Float3 &pos) ;
  float getTumorDens(const Float3 &pos);
  Float3 getGfGrad(const Float3 &pos) const;
  
  //oxygen
/** @brief core interface for vbl BloodVesselVector
 * calculates the distances to the nearest vessels by means of the ANN library.
 */
  void findNearestVessel( DetailedPO2::VesselPO2Storage &po2Store);

  
  /** 
   * main functions
   */
  int run();
  void doStep(double dt);
  
/** @brief
 * writes the vbl based stuff to the HDF5 file
 * returns the H5 Group name 
 */
  std::string writeOutput(bool doPermanentSafe);
  void writeVBLDataToHDF(H5::Group &g);
  void WriteCellsSystemHDF_with_nearest_vessel_index(H5::Group &out_cell_group);
  
  vbl::CellsSystem *tumorcode_pointer_to_currentCellsSystem;
  void doMilottiStep();
  void initMilotti();
};
}//end FakeTum

#endif //_FAKETUM_H_
