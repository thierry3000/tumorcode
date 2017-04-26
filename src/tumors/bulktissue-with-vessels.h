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
#pragma once // include this file only once per compilation unit (see https://en.wikipedia.org/wiki/Pragma_once)

#include <fenv.h>
#include <algorithm>
#include <stdexcept>

#include "vesselmodel1.h"
#include "bulktissuemodel1_new.h"
#include "growthfactor_model.h"
#include "oxygen_model.h"
#include "mwlib/log.h"

#include "shared-objects.h"
#include "mwlib/hdf_wrapper_ptree.h"
#include "mwlib/ptree_ext.h"
#include "convection_diffusion_solver.h"
#include "calcflow.h"
#include "time_stepper_utils_new.h"

namespace BulkTissue
{
enum TissueId {
  TCS = 0,
  DEAD = 1,
  TISSUE = 2,
};

static const string tissue_name[] = {
  "tumor",
  "necro",
  "normal"
};


struct Params
{
  double override_scale;  // used when loading vessel trees
  double rGf, gf_production_threshold;
  double out_intervall, tend;
  double lattice_scale;
  string message;
  string fn_out, fn_vessel, vesselfile_message;
  string paramset_name;
  int vesselfile_ensemble_index;
  bool create_single_output_file;
  Int3 lattice_size;
  bool vessel_volume_exclusion;
  BloodFlowParameters bfparams;
  double
         reference_intercapillary_distance,
         o2_range[3],
         o2_cons_coeff[3],
         capillary_wall_permeability,
         o2_level_normal,
         hematocrit_init;
  
  Params();
  static void update_ptree(ptree &dst, const ptree &src);
  void assign(const ptree &pt);
  ptree as_ptree() const;
};

struct State : boost::noncopyable
{
  State() :
    tumor_checksum(0), vessels_checksum(0), chem_checksum(0)
    {}
    
  NewBulkTissueModel::State tumor;
  int tumor_checksum;
  boost::scoped_ptr<VesselList3d> vessels;
  int vessels_checksum;
  Array3d<float> o2field;
  Array3d<float> gffield;
  int chem_checksum;
};

struct NewTumorSim : public  boost::noncopyable
{
  // the state
  State state;
  
  // vessel system
  VesselModel1::Model vessel_model; // dynamics   
  Array3d<float> vessel_volume_fraction, vessel_o2src_clin, vessel_o2src_crhs; // cached data, could also be read from hdf
  int last_vessels_checksum;
  NewSteppers::StepControl vessel_step_ctrl;

  // lattice definition of the continuum field lattice
  ContinuumGrid grid;
  DomainDecomposition mtboxes;
  
  // the tumor
  NewBulkTissueModel::Model tumor_model;
  Array3d<float> cell_pressure, tumor_phase; // cached data
  int last_tumor_checksum;
  double tumor_radius_estimate;
  NewSteppers::StepControl tumor_step_ctrl;

  // chemical field data
  GfModel gf_model;
  Array3dOps<float> oxyops;
  double last_chem_update;
  
  // paramters
  ptree all_pt_params;
  Params params;

  // iteration & timing data
  my::Time real_start_time;
  int iteration_num;
  int output_num;
  bool failFlag, checkStop;

  //FieldInterpolator<float> fieldinterp_extrapolate;
  //FieldInterpolator<float> fieldinterp_const;
  
  // interface functions for the vessel model
  float getGf(const Float3 &pos) const;
  float getPress(const Float3 &pos) const;
  float getTumorDens(const Float3 &pos) const;
  Float3 getGfGrad(const Float3 &pos) const;
  void UpdateVesselVolumeFraction();
  void insertO2Coefficients(int box_index, const BBox3& bbox, const State &state, FiniteVolumeMatrixBuilder& mb);

  // main functions
  //int run(int argc, char **argv);
  int run(const ptree &params);
  bool doStep(NewSteppers::StepControl &ctrl);
  void calcChemFields();
  void advanceTumorState(double dt);
  void advanceVesselState();
  void writeOutput(double time);
};
}//end namespace BulkTissue

