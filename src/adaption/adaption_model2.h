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

#define ADAPTION_OUTPUT 1
#include <boost/concept_check.hpp>
#include <boost/graph/graph_concepts.hpp>
#include <boost/archive/archive_exception.hpp>
#include <boost/optional/optional.hpp>

//#include "hdfcppwrapper/hdf_wrapper.h"
#include "H5Cpp.h"
// #ifndef H5_NO_NAMESPACE
//   using namespace H5;
// #endif
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
    BOTH_PRESSURE = 43,
  };
  

  typedef DynArray<FlReal> FlArray;
  typedef DynArray<my::eqpair<int> > FeArray;
  typedef DynArray<my::Bitfield<uchar> > FbArray;
  typedef boost::unordered_map<int, FlowBC> FlBoundaryList;

//   template<class T>
//   inline T getNAN() { return std::numeric_limits<T>::quiet_NaN(); }

  template<class T>
  inline bool isFinite(T x) { return std::isfinite(x); }


  struct Parameters
  {
    //Parameters& operator=(const Parameters&);
    //Parameters(const Parameters &obj);
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
    double tum_manitulate_s1,tum_manitulate_s2,tum_manitulate_s3,tum_manitulate_s4,tum_manitulate_s5;
    bool write2File;
    string outputFileName;
    string parameterSetName;
    double radMin_for_kill;
    uint boundary_Condition_handling;
    double a_pressure;
    double a_flow;
    
    int pop;
    int individuals;
    int opt_iter;
    string vesselFileName;
    string vesselGroupName;
    
    void assign(const ptree &pt);
    ptree as_ptree() const;
    
    template<class Archive>
      void serialize(Archive &ar, const unsigned int version)
      {
  ar & parameterSetName;
	ar & k_c;
	ar & k_m;
	ar & k_s;
	ar & S_0;
	ar & Q_refdot;
	ar & max_nun_iterations;
	ar & qdev;
	ar & starting_radii;
        ar & delta_t;
	ar & no_of_roots;
	ar & max_flow;
	ar & min_flow;
	ar & avg_flow;
	ar & avgRootNodeConductivity;
	ar & cond_length;
	ar & tum_manitulate_s1;
	ar & tum_manitulate_s2;
	ar & tum_manitulate_s3;
	ar & tum_manitulate_s4;
	ar & tum_manitulate_s5;
	ar & write2File;
  ar & outputFileName;
  
	ar & radMin_for_kill;
	ar & boundary_Condition_handling;
	ar & a_pressure;
	ar & a_flow;
	
	ar & pop;
	ar & individuals;
	ar & opt_iter;
  ar & vesselFileName;
  ar & vesselGroupName;
      }
  };
  ostream& operator<<(ostream &os, const Parameters &params);

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

  void SetAdaptionValues( VesselList3d &vl,
		      CompressedAdaptionNetwork &fl,
		      double delta_t,
		      double max_r
		      );
  uint GetAdaptionNetwork(CompressedAdaptionNetwork &fl,
		      const VesselList3d *vl
			);

  void TestAdaption();
  //double CalcRadiiChange(const Parameters &params, VesselList3d &vl);
  //void CalcRadiiChange2(const Parameters &params, VesselList3d &vl);
  std::tuple<FlReal,FlReal,FlReal> CalcRadiiChange_mw(const Adaption::Parameters &params, VesselList3d &vl, float delta_t_calc);
  /// return state saying if convergent and mean of capillary flow 
  //std::tuple<uint,FlReal> runAdaption_Loop(Adaption::Parameters params, BloodFlowParameters bfparams, VesselList3d &vl, bool doDebugOutput);
  std::tuple<uint,FlReal,FlReal,FlReal> runAdaption_Loop(Adaption::Parameters params, BloodFlowParameters bfparams, bool doDebugOutput);
  std::tuple<uint,FlReal,FlReal,FlReal> runAdaption_Loop(VesselList3d &vl, Adaption::Parameters params, BloodFlowParameters bfparams, bool doDebugOutput);
  //std::tuple<uint, FlReal> runAdaption_Loop(Adaption::Parameters params, BloodFlowParameters bfparams, std::auto_ptr<VesselList3d> vl);
  //uint runAdaption_Loop(boost::shared_ptr<Adaption::Parameters> params, boost::shared_ptr<BloodFlowParameters> bfparams, boost::shared_ptr<VesselList3d> vl, bool doDebugOutput);
//  uint run_optimization(Adaption::Parameters params, BloodFlowParameters bfparams, std::auto_ptr<VesselList3d> vl);
  Adaption::ExtremeFlows GetExtremFlows(VesselList3d *vl);
  //void UpdateBoundaryConditions(VesselList3d &vl);
  void ChangeBoundaryConditions(VesselList3d &vl, const Adaption::Parameters &params);
};//end namespace
#endif
