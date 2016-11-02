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
#ifndef VESSEL_MODEL1_H
#define VESSEL_MODEL1_H

#include "shared-objects.h"
#include "vessels3d.h"
#include <adaption/adaption_model2.h>

#include <boost/function.hpp>
#include <tbb/spin_mutex.h>

namespace VesselModel1
{


struct Params
{
  Params();
  void assign(const ptree &pt);
  static ptree pt_defaults();
  ptree as_ptree() const;
  
  int seed;
  double pressMax,pressMin;

  double gfVessProl;
  int    timeProlEcSprout;
  int    timeProlEcSproutLifetime;
  int    timeEcRadiusInflate;
  int    timeEcRadiusDeflate;

  bool   bRelativeShearforceCollapse;
  bool   onlyReduceMaturationIfUnderperfused;
  bool   isShearForceEffectLinear;
  
  double probCollapse;
  double forceCollapse;
  double forceEnlarge;

  double sproutMaxVesselWallThickness;
  double sproutMaxVesselWallThicknessArterial;

  int    distSproutMin;
  bool   bSproutModelSimple;
  bool   badaption_on_off;
  int    sproutDelay;

  float  radMax;
 
  float  dematuration_rate;
  float  maturation_crit;
  float  radInit;

  float  vesselCompressionFactor;

  bool   bShearStressControlledDilatationAndRadiusDependentCollapse;
  //double radiusShrinkingTendency;
  double radMin;
  bool   bRadiusSmoothing;
  double radiusSmoothingDiffusionCoeff;
};

struct Callbacks
{
  // oh well, the vessel system is very dependent on the other parts of the system.
  // by these functions i want to make the dependencies as loose as possible
  boost::function1<float, Float3> getGf, getPress, getTumorDens;
  boost::function1<Float3, Float3> getGfGrad;
};


struct Model
{
  typedef VesselList3d::LatticeData LatticeData;
  Params params;
  const LatticeData *m_ld;
  Float3 lattice_axis_directions[32];
  VesselList3d* vl;
  HemodynamicBounds hemodynBounds;
  int num_iteration;
  double time;

  const HemodynamicBounds* GetHdBounds() const { return &hemodynBounds; }
  void Init(VesselList3d *vl_, const Params &params_, Callbacks &callbacks_);
  //void DoStep(double dt);
  void DoStep(double dt, const Adaption::Parameters *adap_params, const BloodFlowParameters *bfparams);//for using with adaption

private:
  double dt;
  Callbacks callbacks;
  //interface functions, using the vesselfield lattice, get values from fields
  float  GetGfConc( const Float3 &pos ) const { return callbacks.getGf(pos); }
  float  GetGfConc( const Int3 &pos ) const { return GetGfConc(Ld().LatticeToWorld(pos)); }
  float  GetTissPress( const Int3 &pos ) const { return callbacks.getPress(Ld().LatticeToWorld(pos)); }
  Float3 GetGfGrad( const Int3 &pos ) const { return callbacks.getGfGrad(Ld().LatticeToWorld(pos)); }
  float  GetTumorDens( const Int3 &pos ) const { return callbacks.getTumorDens(Ld().LatticeToWorld(pos)); }
  
  // private stuff
  tbb::spin_mutex main_mutex;
  Random main_rnd;
  // this is thread safe and should be used to seed thread local random generators
  uint GetThreadRandomSeed(); 
  
  void ClassifyTumorVessels();
  void GenerateSprouts();
  Vessel* GenerateSprout(Random& rnd, const Int3& pos, const Float3& forward_dir);
  Vessel* AddSproutVessel(VesselNode* vc, VesselNode* dst_vc);
  bool CheckCanSprout(Random &rnd, const Int3 &lpos, const VesselNode* nd, const Vessel* v);
  
  void CollapseVessels();
  void EnlargeVessels();
  void SmoothVessels();
  void MaturateVessel();
  const LatticeData& Ld() const { return *m_ld; }

  //int FindInitialSproutDir( const Int3 &pos, const Vessel* v, const VesselNode* vc );
  //bool CheckSproutTarget( int dir, Int3 &pos, int &len, Vessel* &vDst, VesselNode* &vcDst );
  //Vessel* AddSproutVessel( VesselNode* vc, const Int3 &lposDst, Vessel* vDst, VesselNode* vcDst );
  int DistanceToJunction( const Int3 &pos, int maxdist );

  float GetAverageGrowthfactorOverSegment(const Vessel* v);
  void EnlargeRadiusSubstep(Random &rnd, Vessel* v, float fr);
  double GetProbCollapse( const Vessel* v, const Int3 &pos  );
  double StabilityAmountDueToShearStress(const Vessel* v);
};

void myprint(boost::property_tree::ptree const& pt);

}

#endif