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
#ifdef USE_ADAPTION
  #include "adaption/adaption_model2.h"
#endif
#include "common/vesselmodel1.h"
#include "common/growthfactor_model.h"

#include "mwlib/log.h"
#include "../detailedO2/oxygen_model2.h"
#include "../common/simple_oxygen_model.h"
#include "../common/glucose_model.h"

#include "trilinos_linsys_construction.h"

/** from milotti
 */
//#define undo

#ifdef MILOTTI_MTS
#ifndef undo
  #include "sim.h"
  #include "InputFromFile.h"
  #include "CellType.h"
  #include "Environment.h"
  #include "EnvironmentalSignals.h"
  #include "geom-2.h"
  #include "CellsSystem.h"
#endif
#endif

#define USE_DETAILED_O2

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
  Array3d<float> o2field;
  Array3d<float> glucoseField;
  //Array3d<float> gffield;
  int chem_checksum;
};  
enum TissuePressureDistribution
{
  TISSUE_PRESSURE_SPHERE = 0,
  TISSUE_PRESSURE_SHELL = 1,
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
  BloodFlowParameters bfparams;
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
struct FakeTumorSimMTS : public boost::noncopyable
{
  std::auto_ptr<VesselList3d> vl;
  VesselModel1::Model model;
#ifdef USE_DETAILED_O2
  DetailedPO2::DetailedP02Sim o2_sim;
#else
  void insertO2Coefficients(int box_index, const BBox3& bbox, const State &state, FiniteVolumeMatrixBuilder& mb);
#endif
  // simple o2
  void calcChemFields();
  // lattice definition of the continuum field lattice
  ContinuumGrid grid;
  DomainDecomposition mtboxes;
  // the state
  State state;
  
  void insertGlucoseCoefficients(int box_index, const BBox3& bbox, const State &state, FiniteVolumeMatrixBuilder& mb);
  Array3d<float> vessel_volume_fraction;
  Array3d<float> vessel_o2src_clin, vessel_o2src_crhs;
  Array3d<float> vessel_glucosesrc_clin, vessel_glucosesrc_crhs;
  Array3dOps<float> oxyops;
  Array3dOps<float> glucoseOps;
  double last_chem_update;
  int last_vessels_checksum;
  void UpdateVesselVolumeFraction();
  
//   GfModel gf_model;
  TissuePhases phases;//Declaration
  
  FakeTumMTS::Parameters params;
  GlucoseModel::GlucoseParams glucoseParams;
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
  float getGf(const Float3 &pos) const;
  float getPress(const Float3 &pos) const;
  float getTumorDens(const Float3 &pos) const;
  Float3 getGfGrad(const Float3 &pos) const;

  // main functions
  //int run(int argc, char **argv);
  int run(const ptree &params);
  void doStep(double dt);
  std::string writeOutput();
  
  //milotti mts
#ifndef undo
  CellsSystem currentCellsSystem;
  void doMilottiStep();
#endif
};
}//end FakeTum
#endif //_FAKETUM_H_
