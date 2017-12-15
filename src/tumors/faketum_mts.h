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
#ifdef USE_ADAPTION
  #include "adaption/adaption_model2.h"
#endif
#include "common/vesselmodel1.h"
#include "common/growthfactor_model.h"

#include "mwlib/log.h"


#include "../common/glucose_model.h"

#include "trilinos_linsys_construction.h"

/** from milotti
 */

#ifdef MILOTTI_MTS
#define ANN
  #ifdef ANN
    #include <ANN/ANN.h>
  #endif
  #include <vbl.h>
#endif

#define USE_DETAILED_O2

#ifdef USE_DETAILED_O2
  #include "../detailedO2/oxygen_model2.h"
#else
  #include "../common/simple_oxygen_model.h"
#endif

#ifdef W_timing
  #include <chrono>
#endif
namespace h5 = h5cpp;

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
  double apply_adaption_intervall;
  string message;
  string fn_out, fn_vessel, vesselfile_message;
  string paramset_name;
  int vesselfile_ensemble_index;
  double rGf, gf_production_threshold;
  double tumor_radius, tumor_speed;
  double stopping_radius_fraction; // actual-stopping-radius = 0.5*system-length-along-the-longest-dimension * stopping_radius_fraction
  int num_threads;
  int tissuePressureDistribution;
  double tissuePressureWidth;
  double tissuePressureCenterFraction;
  
#ifdef USE_ADAPTION
  Adaption::Parameters adap_params;
#endif
  //for continuum oxygen calculation
  Int3 lattice_size = {20,20,20};
  double lattice_scale = 10;
  
  Parameters();
  void assign(const ptree &pt);
  ptree as_ptree() const;
  static void update_ptree(ptree &dst, const ptree &src);
};
#ifdef W_timing
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
//   void calculate_timings();
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
#endif
struct FakeTumorSimMTS : public boost::noncopyable
{
  Timing currentTiming;
  std::auto_ptr<VesselList3d> vl;
  // ANN stuff
//   ANNkd_tree* kd_tree_of_vl;        // ann kd tree structurs
  const int ANN_dim = 3;            // space dimension
  const int ANN_maxPts = 25000;      // maximum number of data points --> to limit memory allocation
  const double ANN_eps = 0.0;       // error bound
//   ANNpointArray    dataPts;         // data points
//   ANNpoint         queryPt;         // query point
//   ANNidxArray      ANN_nnIdx;       // near neighbor indices   --> will be filled during search
//   ANNdistArray     ANN_dists;       // near neighbor distances --> will be filled during search
  const int	   ANN_k= 5;        // number of nearest neighbors
  
  DynArray<nearest> vectorOfnearestVessels;
	
	
  VesselModel1::Model vessel_model;
#ifdef USE_DETAILED_O2
  DetailedPO2::DetailedP02Sim o2_sim;
#else
  void insertO2Coefficients(int box_index, const BBox3& bbox, const State &state, FiniteVolumeMatrixBuilder& mb);
#endif
  
  void calcChemFields();
  // lattice definition of the continuum field lattice
  ContinuumGrid grid;
  DomainDecomposition mtboxes;
  // the state
  State state;
  
//   void insertGlucoseCoefficients(int box_index, const BBox3& bbox, const State &state, FiniteVolumeMatrixBuilder& mb);
  void insertGFCoefficients(int box_index, const BBox3& bbox, const State &state, FiniteVolumeMatrixBuilder& mb);
  Array3d<float> vessel_volume_fraction;
//   Array3d<float> vessel_o2src_clin, vessel_o2src_crhs;
//   Array3d<float> vessel_glucosesrc_clin, vessel_glucosesrc_crhs;
  Array3d<float> cell_GFsrc;
  Array3dOps<float> oxyops;
  Array3dOps<float> glucoseOps;
  double last_chem_update;
  int last_vessels_checksum;
  void UpdateVesselVolumeFraction();
  
//   GfModel gf_model;
  TissuePhases phases;//Declaration
  
  FakeTumMTS::Parameters params;
  //GlucoseModel::GlucoseParams glucoseParams;
  GfModel_Cell gf_model;
  
  BloodFlowParameters bfparams;
#ifdef USE_ADAPTION
  Adaption::Parameters adap_params;
#endif
#ifdef USE_DETAILED_O2
  DetailedPO2::Parameters o2_params;
#else
  O2Model::SimpleO2Params o2_params;
#endif
  ptree all_pt_params;

  double tumor_radius;
  double time;
  int num_iteration;
  int output_num;

  // interface functions
  //float getGf(const Float3 &pos) const;
  float getGf(const Float3 &pos);
  //float getPress(const Float3 &pos) const;
  float getPress(const Float3 &pos) ;
  //float getTumorDens(const Float3 &pos) const;
  float getTumorDens(const Float3 &pos);
  Float3 getGfGrad(const Float3 &pos) const;
  
  // main functions
  //int run(int argc, char **argv);
  int run();
  void doStep(double dt);
  std::string writeOutput();
  
  //milotti mts
  vbl::CellsSystem *currentCellsSystem;
  void doMilottiStep();
  //void update_milotti_vessels(vbl::CellsSystem &currentCellSys, VesselList3d &vl, DetailedPO2::VesselPO2Storage &po2Store);
  //void WriteCellsSystemHDF(h5cpp::Group &out_cell_group);
  void WriteCellsSystemHDF_with_nearest_vessel_index(h5cpp::Group &out_cell_group);
  //void fillKdTreeFromVl();
  void findNearestVessel( DetailedPO2::VesselPO2Storage &po2Store);
  float estimateTumorRadiusFromCells();
};
}//end FakeTum

#endif //_FAKETUM_H_
