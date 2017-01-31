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
#ifndef CALCFLOW_H
#define CALCFLOW_H

#pragma once // include this file only once per compilation unit (see https://en.wikipedia.org/wiki/Pragma_once)

#include "common.h"
#include "mwlib/math_ext.h"

#include <boost/unordered_map.hpp>

typedef double FlReal;
class VesselList3d;

/*
  units:
    time - sec
    space - micron
    pressure - kPa
*/


enum Rheology
{
  RheologyForHuman = 0,
  RheologyForRats = 1,
  /* according to Pries and Secomb 2005
   * American Journal of Physiology - Heart and Circulatory Physiology
   * Microvascular blood viscosity in vivo and the endothelial surface layer.
   */
  RheologySecomb2005 = 2
};
ostream& operator<<(ostream &os, Rheology rheology);
istream& operator>>(istream &os, Rheology &rheology);

struct BloodFlowParameters
{
  BloodFlowParameters();
  void assign(const ptree &pt);
  ptree as_ptree() const;
  
  double viscosityPlasma; // [kPa s]
  Rheology rheology;
  double inletHematocrit;
  bool includePhaseSeparationEffect; // non-homogeneous intravascular hematocrit distribution, Warning: computationally very expensive
};

ostream& operator<<(ostream &os, const BloodFlowParameters &bfparams);



typedef double FlReal;
FlReal CalcRelViscosity( FlReal r, FlReal h, Rheology rheology);
//FlReal CalcRelViscosityByTable( FlReal r, FlReal h, Rheology rheology);
FlReal CalcFahraeusEffect(double h, double r, Rheology rheology);
FlReal PressureRadiusRelation( FlReal rad, bool bArtery );

inline FlReal CalcFlowCoeff( FlReal viscos, FlReal rad, FlReal len )
{
  myAssert( viscos>0 );
  myAssert( len>0);
  FlReal r = rad;
  r = r*r;
  r = r*r;
  const FlReal l = len;
#if 0
#ifdef DEBUG
  printf("r**4: %f, viscos: %f, %f, all: %f\n", r, viscos, l, (my::mconst::pi()/8.0)*r/(viscos*l));
#endif
#endif
  return (my::mconst::pi()/8.0)*r/(viscos*l);
}

inline FlReal CalcShearFromFlow( FlReal q, FlReal rad, FlReal visc )
{
  return FlReal(4.0/my::mconst::pi())*visc/std::pow(rad,3.)*q;
}

inline FlReal CalcShearFromFlowAndCond( FlReal q, FlReal rad, FlReal len, FlReal cond)
{
#if 0
#ifdef DEBUG
  printf("rad: %f, cond: %f, len: %f, q: %f \n", rad, cond,len,q);
#endif
#endif
  return 0.5 * rad / (cond * len) * q;
}

inline FlReal CalcVelocity( FlReal q, FlReal r )
{
  return q/(r*r*my::mconst::pi());
}

struct FlowBC
{
  enum Pin { PIN = 1 };  // set the node to a fixed blood pressure
  enum Current { CURRENT = 2 }; // a fixed blood flow rate into the node
  enum Resist { RESIST = 3 }; // a resistor between the node and a fixed blood pressure
  FlowBC() : type(PIN),w(0.),val(0.) {}
  FlowBC(Current, double current) : w(0),val(current),type(CURRENT) {}
  FlowBC(Pin, double val) : w(0),val(val),type(PIN) {}
  FlowBC(Resist, double w, double val) : w(w),val(val),type(RESIST) {}
  uint type; // PIN, CURRENT or RESIST
  double w,val; // w = flow conuctivity (?) of series "resistor", val = either blood pressure or flow rate
};


void CalcFlowSimple(VesselList3d &vl, const BloodFlowParameters &params, bool keepTheVesselHematocrit); // Either override the vessel hematocrit with a constant or leave the vessel hematocrit values alone and just compute pressure, flow and force using the segments hematocrit
void CalcFlow(VesselList3d &vl, const BloodFlowParameters &params);
void ChangeBoundaryConditions(VesselList3d &vl, uint boundary_handling  );
//void ChangeBoundaryConditionsFix(VesselList3d &vl);
//void ChangeBoundaryConditionsLARGE_2D(VesselList3d &vl);
bool MarkFailingVesselHidden(VesselList3d &vl);
#endif // CALCFLOW_H
