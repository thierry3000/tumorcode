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

#include "calcflow_linsys.h"
#include "common.h"
#include "mwlib/compressed_row_undirected_graph.h"
#include "mwlib/dynamicarray.h"
#include <boost/unordered_map.hpp>
#include <boost/property_tree/info_parser.hpp>


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
  RheologySecomb2005 = 2,
  RheologyEnd = 3
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
  friend class boost::serialization::access;
  template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
      ar & viscosityPlasma;
      ar & rheology;
      ar & inletHematocrit;
      ar & includePhaseSeparationEffect;
    }
};

ostream& operator<<(ostream &os, const BloodFlowParameters &bfparams);

//DA FUCK!!!! When did i ever think that global variables are a good thing ...
//extern BloodFlowParameters bloodFlowParameters;

// template<class T>
// inline T getNAN() { return std::numeric_limits<T>::quiet_NaN(); }

template<class T>
inline bool isFinite(T x) { return std::isfinite(x); }

struct CompressedFlowNetwork
{
  FlArray press;
  FlArray len, rad, hema, flow;
  FlBoundaryList bcs;

  DynArray<my::eqpair<int> > edges;
  /*
   * This map is needed in order to get rid of the unperfused vessels 
   * while mantaining the network structure.
   * NB Matrix solver do not like gaps in indeces!
   */
  IntegerMap org2new_vertex;

  int num_vertices() const { return org2new_vertex.num_added(); }
  int num_edges() const { return edges.size(); }
};


template<class Map1, class Map2, class Keymap>
static void remap_keys(const Map1 &src, Map2 &dst, const Keymap &keymap)
{
  typedef typename Map1::const_iterator Iter;
  for (Iter it = src.begin(); it != src.end(); ++it)
  {
#if 0
#ifdef DEBUG
    cout<<format("first: %i, second: %i, keymap[first]: %i\n") % it->first->Index() % it->second.val % keymap(it->first);
#endif
#endif
    dst[keymap(it->first)] = it->second;
#if 0
#ifdef DEBUG
    cout<<format("dst[keymap(it->first)] %f\n") % dst[keymap(it->first)].val;
#endif
#endif
  }
}

void CalcViscosities( const FlArray &rad, const FlArray &hema, const BloodFlowParameters &bloodFlowParameters, FlArray &visc);
void CalcConductivities(const FlArray &rad, const FlArray &len, const FlArray &visc, FlArray &cond);

void SetFlowValues( VesselList3d &vl,
                    const CompressedFlowNetwork &fl,
                    const FlArray &cond,
                    const FlArray &press,
                    const FlArray &hema
		    );
uint GetFlowNetwork(CompressedFlowNetwork &fl,
                    const VesselList3d &vl,
                    const BloodFlowParameters &bfparams,
                    bool keepTheVesselHematocrit
		    );

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




void CalcFlowSimple(VesselList3d &vl, const BloodFlowParameters &params, bool keepTheVesselHematocrit); // Either override the vessel hematocrit with a constant or leave the vessel hematocrit values alone and just compute pressure, flow and force using the segments hematocrit
void CalcFlow(VesselList3d &vl, const BloodFlowParameters &params);
bool MarkFailingVesselHidden(VesselList3d &vl);

#endif
