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
#ifndef _FAKETUMMTS_H_
#define _FAKETUMMTS_H_

#include "common/common.h"
#ifdef USE_ADAPTION
  #include "adaption/adaption_model2.h"
#endif
#include "common/vesselmodel1.h"

#include "mwlib/log.h"

/** from milotti
 */
#ifdef MILOTTI_MTS
  #include "sim.h"
  #include "InputFromFile.h"
  #include "CellType.h"
  #include "Environment.h"
  #include "EnvironmentalSignals.h"
  #include "geom-2.h"
  #include "CellsSystem.h"
#endif

namespace FakeTumMTS{
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
  double rGf;
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
  
  Parameters();
  void assign(const ptree &pt);
  ptree as_ptree() const;
  static void update_ptree(ptree &dst, const ptree &src);
};
struct FakeTumorSimMTS : public boost::noncopyable
{
  std::auto_ptr<VesselList3d> vl;
  VesselModel1::Model model;

  FakeTumMTS::Parameters params;
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
  void writeOutput();
  
  //milotti mts

  CellsSystem currentCellsSystem;
  void doMilottiStep();

};
}//end FakeTum
#endif //_FAKETUM_H_
