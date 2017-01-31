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
#ifndef CALCFLOW_COMMON_H
#define CALCFLOW_COMMON_H

#include "calcflow.h"

#include "mwlib/compressed_row_undirected_graph.h"
#include "mwlib/dynamicarray.h"

typedef DynArray<FlReal> FlArray;
typedef DynArray<my::eqpair<int> > FeArray;
typedef DynArray<my::Bitfield<uchar> > FbArray;
typedef boost::unordered_map<int, FlowBC> FlBoundaryList;

//DA FUCK!!!! When did i ever think that global variables are a good thing ...
//extern BloodFlowParameters bloodFlowParameters;

template<class T>
inline T getNAN() { return std::numeric_limits<T>::quiet_NaN(); }

template<class T>
inline bool isFinite(T x) { return std::isfinite(x); }

enum BoundaryHandling
  {
    KEEP = 0,
    VEIN_AS_FLOW_ARTERY_PRESSURE  = 1,
    LARGE_2D = 2,
    LARGE_2D_2 = 3,
    LARGE_2D_like_paper = 4,
  };

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

void SetFlowValues( VesselList3d* vl,
                    const CompressedFlowNetwork &fl,
                    const FlArray &cond,
                    const FlArray &press,
                    const FlArray &hema
		    );
void GetFlowNetwork(CompressedFlowNetwork &fl,
                    const VesselList3d* vl,
                    const BloodFlowParameters &bfparams,
                    bool keepTheVesselHematocrit
		    );

#endif