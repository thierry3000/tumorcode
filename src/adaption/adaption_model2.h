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
#ifndef _ADAPTION_MODEL_H_
#define _ADAPTION_MODEL_H_
#include <boost/concept_check.hpp>
#include <boost/graph/graph_concepts.hpp>

//#include "hdfcppwrapper/hdf_wrapper.h"
#include "mwlib/compressed_row_undirected_graph.h"
#include "mwlib/dynamicarray.h"

#include "common/shared-objects.h"
#include "common/continuum-grid.h"
#include "common/vessels3d.h"

#include "common/calcflow_common.h" //inhere from the CompressedFlowNetwork

namespace Adaption
{
  enum TissueIDs
  {
    NORMAL = 0,
    TUMOR  = 1,
    NECRO  = 2
  };
  enum BoundaryHandling
  {
    KEEP = 0,
    VEIN_AS_FLOW_ARTERY_PRESSURE  = 1,
    LARGE_2D = 2,
    LARGE_2D_2 = 3,
    LARGE_2D_3 = 4,
    LARGE_2D_like_paper = 5,
    VALUE = 42,
  };
  

  typedef DynArray<FlReal> FlArray;
  typedef DynArray<my::eqpair<int> > FeArray;
  typedef DynArray<my::Bitfield<uchar> > FbArray;
  typedef boost::unordered_map<int, FlowBC> FlBoundaryList;

  template<class T>
  inline T getNAN() { return std::numeric_limits<T>::quiet_NaN(); }

  template<class T>
  inline bool isFinite(T x) { return std::isfinite(x); }


  class Parameters
  {
    
    Parameters& operator=(const Parameters&);
    
  public:
    Parameters();
    double k_c;
    double k_m;
    double k_s;
    double S_0;
    double Q_refdot;
    int max_nun_iterations;
    double qdev;
    double starting_radii;
    double delta_t;
    int no_of_roots;
    double max_flow;
    double min_flow;
    double avg_flow;
    double avgRootNodeConductivity;
    double cond_length;
    bool tum_manitulate_s1;
    bool tum_manitulate_s2;
    bool tum_manitulate_s3;
    bool tum_manitulate_s4;
    bool tum_manitulate_s5;
    bool write2File;
    
    double radMin_for_kill;
    uint boundary_Condition_handling;
    double a_pressure;
    double a_flow;
    
    void assign(const ptree &pt);
    ptree as_ptree() const;
  };

  struct CompressedAdaptionNetwork : CompressedFlowNetwork
  {
    FlArray shearforce, metabolic, conductive, condSignal;
    FlArray S_tot_array;
    FlReal max_delta_r;
    FlReal max_stot;
    FlReal delta_r_square;
    DynArray<bool> isVein;
    DynArray<bool> isArtery;
    DynArray<bool> isCapillary;
    DynArray<bool> inTumor;
    DynArray<bool> isBoundaryVessel;
    //AdaptionParameters myParameters;
    /*
    * This map is needed in order to get rid of the unperfused vessels 
    * while mantaining the network structure.
    * NB Matrix solver do not like gaps in indeces!
    */
    IntegerMap newEdge2original;//other direction already defined in mother struct
    IntegerMap org2newEdge;
  };
  
  struct ExtremeFlows{
    double max_flow=0.0;
    double min_flow=std::numeric_limits<double>::max();
    double avg_flow=0.0;
  };

  void SetAdaptionValues( VesselList3d* vl,
		      CompressedAdaptionNetwork &fl,
		      double delta_t,
		      double max_r
		      );
  void GetAdaptionNetwork(CompressedAdaptionNetwork &fl,
		      const VesselList3d* vl
			);

  void TestAdaption();
  //double CalcRadiiChange(const Parameters &params, VesselList3d &vl);
  //void CalcRadiiChange2(const Parameters &params, VesselList3d &vl);
  std::tuple<FlReal,FlReal,FlReal> CalcRadiiChange_mw(const Adaption::Parameters &params, VesselList3d &vl, float delta_t_calc);
  uint runAdaption_Loop(const Adaption::Parameters *params, const BloodFlowParameters *bfparams, VesselList3d *vl, h5cpp::Group *vessels_after_adaption, bool doDebugOutput);
  void GetExtremFlows(const VesselList3d *vl, Adaption::ExtremeFlows *myExtrems);
  //void UpdateBoundaryConditions(VesselList3d &vl);
  void ChangeBoundaryConditions(VesselList3d &vl, const Adaption::Parameters &params);
};//end namespace
#endif