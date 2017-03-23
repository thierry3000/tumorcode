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
#include "common/shared-objects.h"
#include "common/continuum-flow.h"
#include "common/time_stepper_utils_new.h"
#include "common/trilinos_linsys_construction.h"
#include "common/vessels3d.h"
#include "common/calcflow_common.h" //at least for remap_keys

#include "mwlib/math_ext.h"
//#include "hdf_wrapper.h"
//#include "hdf_wrapper_stl.h"

#include <boost/math/tools/roots.hpp>
#include <boost/math/tools/tuple.hpp>
#include <boost/foreach.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
//#include <boost/container/vector.hpp> //used to be threadsafe

#include <unordered_set>
#define OUTPUT_PO2INTEGRATION(x)
#define OUTPUT_PO2MEASURECOMP(x)

#include "adaption_model2.h"
#include <fenv.h>
//#include <atomic>
//#include "adaption_as_pagmo_problem.h"

namespace Adaption
{
// int test_2()
// {
//   std::cout<<"running pagmo test 2"<<std::endl;
//   // Initialise the MPI environment.
//   int mc_steps=0;
//   int dim=0;
// #ifdef PAGMO_ENABLE_MPI
//   pagmo::mpi_environment env;
//   mc_steps = 10000000;
//   dim = 400;
// #else
//   mc_steps = 10;
//   dim = 4;
// #endif
//   // Create a problem and an algorithm.
//   pagmo::problem::dejong prob(dim);
//   pagmo::algorithm::monte_carlo algo(mc_steps);
//   // Create an archipelago of 10 MPI islands.
//   pagmo::archipelago a;
//   a.set_topology(pagmo::topology::ring());
//   for (int i = 0; i < 2; ++i) {
// #ifdef PAGMO_ENABLE_MPI
// 	  a.push_back(pagmo::mpi_island(algo,prob,1));
// #else
// 	  a.push_back(pagmo::island(algo,prob,1));
// #endif
//   }
//   // Evolve the archipelago 10 times.
//   a.evolve(10);
//   a.join();
//   return 0;
// }  
// int test_1()
// {
//   std::cout<<"running pagmo test 1"<<std::endl;
//   // Initialise the MPI environment.
//   int mc_steps=0;
//   int dim=0;
// #ifdef PAGMO_ENABLE_MPI
//   pagmo::mpi_environment env;
//   mc_steps = 10000000;
//   dim = 400;
// #else
//   mc_steps = 10;
//   dim = 4;
// #endif
//   // Create a problem and an algorithm.
//   pagmo::problem::dejong prob(dim);
//   pagmo::algorithm::monte_carlo algo(mc_steps);
//   // Create an archipelago of 10 MPI islands.
//   pagmo::archipelago a;
//   a.set_topology(pagmo::topology::ring());
//   for (int i = 0; i < 10; ++i) {
// #ifdef PAGMO_ENABLE_MPI
// 	  a.push_back(pagmo::mpi_island(algo,prob,1));
// #else
// 	  a.push_back(pagmo::island(algo,prob,1));
// #endif
//   }
//   // Evolve the archipelago 10 times.
//   a.evolve(10);
//   a.join();
//   return 0;
// }  
void ChangeBoundaryConditions(VesselList3d &vl, const Adaption::Parameters &params)
{
#ifdef DEBUG
  printf("entered ChangeBoundaryConditions\n");
#endif
  uint handling=params.boundary_Condition_handling;
  int ncnt = vl.GetNCount();
  double min_boundary_pressure = std::numeric_limits< double >::max();
  double max_boundary_flow = std::numeric_limits< double >::min();
  switch(handling){
    case KEEP:
#ifdef DEBUG
      std::printf("no boundary changes in");
#endif
      break;
    case VEIN_AS_FLOW_ARTERY_PRESSURE:
      
      for( int i=0;i<ncnt; ++i)
      {
	const VesselNode* vc= vl.GetNode(i);
	if( vc->Count() >0 )
	{
	  if( vc->IsBoundary() and vc->GetEdge(0)->IsCirculated())
	  {
	    if( vc->press < min_boundary_pressure)
	    {
	      min_boundary_pressure = vc->press;
	    }
	    if( vc->GetEdge(0)->q > max_boundary_flow)
	    {
	      max_boundary_flow = vc->GetEdge(0)->q;
	    }
	  }
	}
      }
      for(int i=0; i<ncnt; ++i)
      {
	const VesselNode* vc= vl.GetNode(i);
	// insert boundary nodes into the bcs array
	if (vc->Count() > 0 and vc->IsBoundary() and vc->GetEdge(0)->IsCirculated())
	{
    #if 0  //keep max flow as flow boundary, everything else to pressure
	  if(vc->GetEdge(0)->q== max_boundary_flow)
	  {
	    if(vc->GetEdge(0)->IsArtery())
	    {
	      vl.SetBC(vc,FlowBC(FlowBC::CURRENT, -vc->GetEdge(0)->q));
	    }
	    if(vc->GetEdge(0)->IsVein())
	    {
	      vl.SetBC(vc,FlowBC(FlowBC::CURRENT, vc->GetEdge(0)->q));
	    }
	  }
	  else
	  {
	    vl.SetBC(vc,FlowBC(FlowBC::PIN, vc->press));
	  }
    #endif
    #if 0
	  //all arteries flow, all veins pressure
	  if(vc->GetEdge(0)->IsArtery())
	  {
	    vl.SetBC(vc,FlowBC(FlowBC::CURRENT, -vc->GetEdge(0)->q));
	  }
	  if(vc->GetEdge(0)->IsVein())
	  {
	    vl.SetBC(vc,FlowBC(FlowBC::PIN, vc->press));
	  }
    #endif
    #if 1
	  //minmal boundary pressure becomes PIN condition everthing else flow
	  if(vc->press == min_boundary_pressure)
	  {
	    vl.SetBC(vc,FlowBC(FlowBC::PIN, vc->press));
    #ifdef DEBUG
	    cout << format("min pressure boundary: %f") % vc->press <<endl;
    #endif
	  }
	  else
	  {
	    if( vc->GetEdge(0)->IsArtery() )
	    {
	      vl.SetBC(vc,FlowBC(FlowBC::CURRENT, -vc->GetEdge(0)->q));
    #ifdef DEBUG
	      cout << format("flow boundary: %f") % -vc->GetEdge(0)->q <<endl;
    #endif
	    }
	    if( vc->GetEdge(0)->IsVein())
	    {
	      vl.SetBC(vc,FlowBC(FlowBC::CURRENT, vc->GetEdge(0)->q));
    #ifdef DEBUG
	      cout << format("flow boundary: %f") % vc->GetEdge(0)->q <<endl;
    #endif
	    }
	  }
    #endif
	  //cout<<format("press : %f\n") % vc->press;
	  //cout << format("flow boundary node: %i") % id << endl;
    #ifdef DEBUG
	  cout<<"Changed bcs map!"<<endl;
    #endif
	}
      }
      break;
    case LARGE_2D:
      for(int i=0; i<ncnt; ++i)
      {
	const VesselNode* vc= vl.GetNode(i);
	// insert boundary nodes into the bcs array
	if (vc->Count() > 0 and vc->IsBoundary() and vc->GetEdge(0)->IsCirculated() and vc->GetEdge(0)->IsVein())
	{
	  vl.SetBC(vc,FlowBC(FlowBC::PIN, 3.7));
	}
	if (vc->Count() > 0 and vc->IsBoundary() and vc->GetEdge(0)->IsCirculated() and vc->GetEdge(0)->IsArtery())
	{
	  vl.SetBC(vc,FlowBC(FlowBC::CURRENT, -500000.0));
	}
      }
      break;
    case LARGE_2D_2:
      for(int i=0; i<ncnt; ++i)
      {
	const VesselNode* vc= vl.GetNode(i);
	// insert boundary nodes into the bcs array
	if (vc->Count() > 0 and vc->IsBoundary() and vc->GetEdge(0)->IsCirculated() and vc->GetEdge(0)->IsVein())
	{
	  vl.SetBC(vc,FlowBC(FlowBC::PIN, 1.8));
	}
	if (vc->Count() > 0 and vc->IsBoundary() and vc->GetEdge(0)->IsCirculated() and vc->GetEdge(0)->IsArtery())
	{
	  vl.SetBC(vc,FlowBC(FlowBC::CURRENT, -800000.0));
	}
      }
      break;
    case LARGE_2D_3:
      for(int i=0; i<ncnt; ++i)
      {
	const VesselNode* vc= vl.GetNode(i);
	// insert boundary nodes into the bcs array
	if (vc->Count() > 0 and vc->IsBoundary() and vc->GetEdge(0)->IsCirculated() and vc->GetEdge(0)->IsVein())
	{
	  vl.SetBC(vc,FlowBC(FlowBC::PIN, 3.5));
	}
	if (vc->Count() > 0 and vc->IsBoundary() and vc->GetEdge(0)->IsCirculated() and vc->GetEdge(0)->IsArtery())
	{
	  vl.SetBC(vc,FlowBC(FlowBC::CURRENT, -800000.0));
	}
      }
      break;
    case LARGE_2D_like_paper:
      for(int i=0; i<ncnt; ++i)
      {
	const VesselNode* vc= vl.GetNode(i);
	// insert boundary nodes into the bcs array
	if (vc->Count() > 0 and vc->IsBoundary() and vc->GetEdge(0)->IsCirculated() and vc->GetEdge(0)->IsVein())
	{
	  vl.SetBC(vc,FlowBC(FlowBC::PIN, 1.8));
	}
	if (vc->Count() > 0 and vc->IsBoundary() and vc->GetEdge(0)->IsCirculated() and vc->GetEdge(0)->IsArtery())
	{
	  vl.SetBC(vc,FlowBC(FlowBC::CURRENT, -16666666.0));
	}
      }
      break;
    case VALUE:
      for(int i=0; i<ncnt; ++i)
      {
	const VesselNode* vc= vl.GetNode(i);
	// insert boundary nodes into the bcs array
	if (vc->Count() > 0 and vc->IsBoundary() and vc->GetEdge(0)->IsCirculated() and vc->GetEdge(0)->IsVein())
	{
	  vl.SetBC(vc,FlowBC(FlowBC::PIN, params.a_pressure));
	}
	if (vc->Count() > 0 and vc->IsBoundary() and vc->GetEdge(0)->IsCirculated() and vc->GetEdge(0)->IsArtery())
	{
	  vl.SetBC(vc,FlowBC(FlowBC::CURRENT, -1*params.a_flow));
	}
      }
      break;
    default:
      cout<<"no proper boundary handling found!"<<endl;
  }//switch
  
  
#ifdef DEBUG
  printf("leave ChangeBoundaryConditions\n");
#endif
}

void TestAdaption()
{
  std::printf("hello world\n");
}
/* @brief default parameters*/
Parameters::Parameters()
{
  k_c = 1.42;
  k_m = 1.42;
  k_s = 1.42;
  Q_refdot = 42;
  max_nun_iterations =42;
  qdev = 42;
  starting_radii = 42;
  delta_t =42;
  no_of_roots = 42;
  max_flow = 42;
  min_flow = 42;
  avgRootNodeConductivity = 0;
  S_0=42;
  cond_length=4242;
  tum_manitulate_s1=tum_manitulate_s2=tum_manitulate_s3=tum_manitulate_s4=tum_manitulate_s5=false;
  /*
   * important value!
   * if adaption without tumor this goes to this default value
   */
  radMin_for_kill = 2.5;
  write2File = true;
  boundary_Condition_handling = KEEP;
  a_pressure = 1.8;
  a_flow = 200000.;
  pop = 5;
  individuals = 42;
  opt_iter = 1;
}

void Parameters::assign(const ptree& pt)
{
  k_c = pt.get<double>("k_c", 1.42);
  k_m = pt.get<double>("k_m", 1.42);
  k_s = pt.get<double>("k_s", 1.42);
  Q_refdot = pt.get<double>("Q_refdot", 42);
  max_nun_iterations = pt.get<int>("max_nun_iterations", 42);
  qdev = pt.get<double>("qdev", 42);
  starting_radii = pt.get<double>("starting_radii", 42);
  delta_t = pt.get<double>("delta_t", 42);
  no_of_roots = pt.get<int>("no_of_roots", 42);
  max_flow = pt.get<double>("max_flow", 42);
  min_flow = pt.get<double>("min_flow", 42);
  avgRootNodeConductivity = pt.get<double>("avgRootNodeConductivity", 42);
  S_0 = pt.get<double>("S_0", 42);
  cond_length = pt.get<double>("cond_length", 42);
  tum_manitulate_s1 = pt.get<bool>("tum_manitulate_s1", false);
  tum_manitulate_s2 = pt.get<bool>("tum_manitulate_s2", false);
  tum_manitulate_s3 = pt.get<bool>("tum_manitulate_s3", false);
  tum_manitulate_s4 = pt.get<bool>("tum_manitulate_s4", false);
  tum_manitulate_s5 = pt.get<bool>("tum_manitulate_s5", false);
  radMin_for_kill = pt.get<double>("radMin_for_kill", 2.5);
  boundary_Condition_handling = pt.get<uint>("boundary_Condition_handling", KEEP);
  a_pressure = pt.get<double>("a_pressure", 1.8);
  a_flow = pt.get<double>("a_flow", 1.8);
  write2File = pt.get<bool>("write2File", true);
  pop = pt.get<int>("pop",5);
  individuals = pt.get<int>("individuals",42);
  opt_iter = pt.get<int>("opt_iter",42);
}

ptree Parameters::as_ptree() const
{
  return make_ptree("k_c", k_c)
                   ("k_m", k_m)
                   ("k_s", k_s)
                   ("Q_refdot", Q_refdot)
		   ("max_nun_iterations", max_nun_iterations)
		   ("qdev", qdev)
		   ("starting_radii", starting_radii)
		   ("delta_t", delta_t)
		   ("no_of_roots", no_of_roots)
		   ("max_flow", max_flow)
		   ("min_flow", min_flow)
		   ("S_0", S_0)
		   ("cond_length", cond_length)
		   ("avgRootNodeConductivity", avgRootNodeConductivity)
		   ("tum_manitulate_s1", tum_manitulate_s1)
		   ("tum_manitulate_s2", tum_manitulate_s2)
		   ("tum_manitulate_s3", tum_manitulate_s3)
		   ("tum_manitulate_s4", tum_manitulate_s4)
		   ("tum_manitulate_s5", tum_manitulate_s5)
		   ("radMin_for_kill", radMin_for_kill)
		   ("boundary_Condition_handling", boundary_Condition_handling)
		   ("a_pressure", a_pressure)
		   ("a_flow", a_flow)
		   ("write2File", write2File)
		   ("pop", pop)
		   ("individuals", individuals)
		   ("opt_iter", opt_iter)
		   ;
}

/* @brief propergates conductive condSignal
 * 
 * Here this is now done like in the HematocritCalculator
 * in calcflow.cpp
 * 
 */
class ConductiveTransport
{
  CompressedAdaptionNetwork &nw; // stores edges and properties
  CompressedRows nbs;  // stores graph neighbor relations in sparse compressed column format
  const Parameters adaptionParams;
  DynArray<uchar> markers; // marks nodes and edges which have been visited
  uchar             iteration_id; //current iteration
  // trivial getters and setters
  FlReal            GetFlow(int i) const                { return nw.flow[i]; }
  FlReal            GetLen(int i) const                 { return nw.len[i]; }
  FlReal            GetPress(int i) const               { return nw.press[i]; }
  FlReal            GetHema(int i) const		{ return nw.hema[i]; }
  my::eqpair<int>   GetEdge(int i) const                { return nw.edges[i]; }
  bool              IsBoundary(int i) const             { return nw.bcs.find(i) != nw.bcs.end(); }
  FlReal	    GetCondSignal(int i) const		{ return nw.condSignal[i]; }
  FlReal	    GetRadius(int i) const		{ return nw.rad[i];}
  void              SetConductive(int i,FlReal x) const;
  void		    SetCondSignal(int i,FlReal x) const;
  void		    SetFlow(int i,FlReal x) const;
  void		    SetHema(int i,FlReal x) const;
  void		    SetMetabolic(int i,FlReal x) const;
  void		    SetS_tot(int i,FlReal x) const;
  void		    set_S_tot_sq(FlReal x) const;
  // not trivial getter an setters
  FlReal GetConductive(int i) const{
    //if( nw.inTumor[i] and adaptionParams.tum_manitulate_s4)
    //{
    //  return 0.0;
    //}
    return nw.conductive[i];
  }
  FlReal GetMetabolic(int i) const{
    //if( nw.inTumor[i] and adaptionParams.tum_manitulate_s3)
    //{
    //  return 0.0;
    //}
    return nw.metabolic[i]; 
  }
  // get neighboring info i.e. node -> #neighbors, neighbor[k] -> edge, node
  int GetNeighborCount(int i) const { return nbs.row_length(i); }
  NodeNeighbor<int,int> GetNeighbor(int i, int k) const
  {
    int id_of_edge = nbs.column_number(i, k);
    const my::eqpair<int> &edge = nw.edges[id_of_edge];
    return make_node_neighbor(edge.first==i ? edge.second : edge.first, id_of_edge);
  }

  void UpdateConductiveRecursiveVenousTree(int iv, int ind_up);
  void UpdateConductiveAtNodeVenous(int ind, int iv);
  void UpdateConductiveRecursiveArterialTree(int iv, int ind_up);
  void UpdateConductiveAtNodeArterial(int ind, int iv);

public:
  bool check_conductive_range;
  ConductiveTransport(CompressedAdaptionNetwork &fl, const Adaption::Parameters& adaptionParams);
  void UpdateConductive();
  void Calculate_S_tot();
  //Parameters getParameters(){ return this->adaptionParams; };
};



ConductiveTransport::ConductiveTransport(CompressedAdaptionNetwork &fl_,const Adaption::Parameters& adaptionParams_) :
  nw(fl_),iteration_id(0), 
  check_conductive_range(false),
  adaptionParams(adaptionParams_)
{
  int ecnt = fl_.edges.size();
  int ncnt = fl_.num_vertices();
  myAssert(nw.shearforce.size()==ecnt);
  myAssert(nw.press.size()==ncnt);
  myAssert(nw.flow.size()==ecnt);
  myAssert(nw.hema.size()==ecnt);
  nw.metabolic.resize(ecnt);
  nw.conductive.resize(ecnt);
  nw.condSignal.resize(ecnt);
  markers.resize(ecnt);
  std::vector<unsigned int> col_counts(ncnt);
  for (int i=0; i<ecnt; ++i)
  {
    double aflow = fl_.flow[i];
    double ahema = fl_.hema[i];
    SetFlow(i, aflow);
    SetHema(i, ahema);
    //failsafe --> after using my own set- functions the values are in reasonable ranges!!!!
    aflow = GetFlow(i);
    ahema = GetHema(i);
    double ametabolic = adaptionParams.k_m*log10(adaptionParams.Q_refdot/(aflow*ahema)+1);
    SetMetabolic(i, ametabolic);
    const my::eqpair<int> &e = nw.edges[i];
    col_counts[e.first]++;
    col_counts[e.second]++;
  }
  nbs.reinit(ncnt, col_counts);
  nbs.allow_identical_edges = false;
  for (int i=0; i<ecnt; ++i)
  {
    const my::eqpair<int> &e = nw.edges[i];
    nbs.add(e.first, i);
    nbs.add(e.second, i);
  }
  markers.resize(std::max(ecnt, ncnt));
}
void ConductiveTransport::SetCondSignal(int i, FlReal x) const
{
  myAssert(std::isfinite(x));
  if (!check_conductive_range) x = my::cut(x, 0., 9999999.);
  else if (x<0. || x>99.) 
  { 
    printf("warning conductive signal %f way to high at vessel %i\n", x, i);
    myAssert(x<0. || x>99.); 
  }
  nw.condSignal[i] = x;
}


void ConductiveTransport::SetConductive(int i,FlReal x) const
{
  myAssert(std::isfinite(x));
  if (!check_conductive_range) x = my::cut(x, 0., 9999999.);
  else if (x<0. || x>99.) 
  { 
    printf("warning conductive signal %f way to high at vessel %i\n", x, i);
    myAssert(x<0. || x>99.); 
  }
  nw.conductive[i]=x;
}
void ConductiveTransport::SetFlow(int i, FlReal x) const
{
  myAssert(std::isfinite(x));
  if ( x<1e-6 )//1e-5 worked for type D,E,F
  {
    printf("Warning: Low flow segment at %i\n", i);
    //myAssert(x > 0.01);
    x = 1e-6;
  }
  nw.flow[i] = x;
}
void ConductiveTransport::SetHema(int i, FlReal x) const
{
  if( x < 0.1 )
  {
    //printf("Warning: Low hema segment at %i\n", i);
    //myAssert(x > 0.001);
    x = 0.1;
  }
  if( x > 0.9 )
  {
    //printf("Warning: High hema segment at %i\n", i);
    //myAssert(x > 0.001);
    x = 0.9;
  }
  nw.hema[i] = x;
}
void ConductiveTransport::SetMetabolic(int i,FlReal x) const
{
  myAssert(std::isfinite(x));
  if (!check_conductive_range) x = my::cut(x, 0., 99.);
  else if (x<0. || x>99.) 
  { 
    printf("warning conductive signal %f way to high at vessel %i\n", x, i);
    myAssert(x<0. || x>99.); 
  }
  this->nw.metabolic[i]=x;
}
void ConductiveTransport::SetS_tot(int i, FlReal x) const
{
  myAssert(std::isfinite(x));
  myAssert(i<nw.num_edges());
  nw.S_tot_array[i] = x;
}

// void ConductiveTransport::set_S_tot_sq(FlReal x) const
// {
//   nw.S_tot_squarre = x;
// }

void ConductiveTransport::UpdateConductive()
{
  ++iteration_id;
  /*march arterial tree*/
  nw.conductive.fill(0.);
  for (int i=0; i<nw.num_edges(); ++i)
  {
    if (nw.isVein[i]) continue;
    if (markers[i] == iteration_id) continue;
    const my::eqpair<int> &e = GetEdge(i);
    int downstream_node = (GetPress(e.first) < GetPress(e.second)) ? e.first : e.second;
    UpdateConductiveRecursiveArterialTree(i, downstream_node);
  }
  /*march venous tree*/
  for (int i=0; i<nw.num_edges(); ++i)
  {
    if (nw.isArtery[i]) continue;
    if (markers[i] == iteration_id) continue;
    const my::eqpair<int> &e = GetEdge(i);
    int upstream_node = (GetPress(e.first) > GetPress(e.second)) ? e.first : e.second;
    UpdateConductiveRecursiveVenousTree(i, upstream_node);
  }
}



void ConductiveTransport::UpdateConductiveRecursiveArterialTree( int edge_id, int downstream_node )
{
  // downstream_node is the node of edge_id which is more downstream than the other
  // we march through the vessel graph with the stream
  // at the capillary nodes the conductive will be set to the metabolic.
  // when the function call stack is rewound all the conductive values are computed
  // in the correct order in upstream direction
  if( markers[edge_id]==iteration_id ) return; //if we were already here

  if( !nw.isCapillary[edge_id] ) // continue downstream segments
  {
    for( int k=0; k<GetNeighborCount(downstream_node); ++k )
    {      
      NodeNeighbor<int,int> nb = GetNeighbor(downstream_node,k);
      if (GetPress(nb.node)>=GetPress(downstream_node)) continue; // do not go upstream the the conductive signal is not yet calculated
      if (nw.isCapillary[nb.edge]) continue; // in this branch we only march arteries
      UpdateConductiveRecursiveArterialTree( nb.edge, nb.node ); // this is an downstream edge, go further downstream
    }
  }
  if( markers[edge_id]==iteration_id ) 
  {
    return;
  }
  UpdateConductiveAtNodeArterial(downstream_node, edge_id);
}

void ConductiveTransport::UpdateConductiveAtNodeArterial(int downstream_node, int edge_id)
{
  // assume downstream segments have already their conductive updated
  // that is segments adjacent to upstream_node, which are not edge_id
  if( nw.isCapillary[edge_id] ) //this means we are at the capillary
  {
    SetConductive(edge_id, GetMetabolic(edge_id));
    markers[edge_id] = iteration_id;
    return;
  }
  int nbcnt = GetNeighborCount(downstream_node);
  if( nbcnt==1 )
  {
    SetConductive(edge_id,0.0);
  }
  else if(nbcnt<=3)
  {
    int numIn=0,numOut=0;
    FlReal flowIn=0;
    FlReal flowOut=0;
    int vin[3];
    int vout[3];
    FlReal conducGoingUp =0;
    ClearMem(vin, 3);
    ClearMem(vout, 3);
    for( int i=0; i<nbcnt; ++i ) 
    {
      const int  iw = GetNeighbor(downstream_node,i).edge;
      const my::eqpair<int> ee = GetEdge(iw);
      const int  imd = ee[0]==downstream_node ? ee[1] : ee[0];
      myAssert(imd == GetNeighbor(downstream_node, i).node);
      if( GetPress(imd) < GetPress(downstream_node) )
      { // outflow
        flowOut += GetFlow(iw);
        vout[numOut++] = iw;
	conducGoingUp += GetConductive(iw) * exp(-GetLen(iw)/adaptionParams.cond_length) + GetMetabolic(iw);
      } 
      else 
      {
        flowIn += GetFlow(iw);
        vin[numIn++] = iw;
      }
    }
    if( numOut==2 && numIn==1 )
    {
      const FlReal S_c = conducGoingUp;
      SetConductive(edge_id, S_c);
      markers[vin[0]] = iteration_id;
    }
    else if( numOut==1 && numIn==2 )
    {
      const int v1 = vin[0];
      const int v2 = vin[1];
      const int w = vout[0];
      SetConductive(v1,conducGoingUp);
      SetConductive(v2,conducGoingUp);
      markers[v1] = iteration_id;
      markers[v2] = iteration_id;
    }
    else if(numOut==1 && numIn==1)
    {
      SetConductive(vin[0], conducGoingUp);
      markers[vin[0]] = iteration_id;
    }
  }
  else // if num outgoing vessels > 2 or num ingoing vessel > 1
  {
    printf("BAD ERROR at arterial Tree, \n");
    printf("Original vessel no.: %i \n", nw.newEdge2original[edge_id]);
    throw std::exception();
    myAssert(nbcnt<4);
  }
  markers[edge_id] = iteration_id;
}

void ConductiveTransport::UpdateConductiveRecursiveVenousTree( int edge_id, int upstream_node )
{
  // upstream_node is the node of edge_id which is more upstream than the other
  // we march through the vessel graph with the opposite of the stream
  // at the capillary nodes the conductive will be set to the metabolic.
  // when the function call stack is rewound all the conductive values are computed
  // in the correct order in downstream direction
  if( markers[edge_id]==iteration_id ) return; //if we were already here

  if( !nw.isCapillary[edge_id] ) // continue downstream segments
  {
    for( int k=0; k<GetNeighborCount(upstream_node); ++k )
    {      
      NodeNeighbor<int,int> nb = GetNeighbor(upstream_node,k);
      if (GetPress(nb.node)<=GetPress(upstream_node)) continue; // do not go downstream the the conductive signal is not yet calculated
      if (nw.isCapillary[nb.edge]) continue; // in this branch we only march veins
      UpdateConductiveRecursiveVenousTree( nb.edge, nb.node ); // this is an downstream edge, go further downstream
    }
  }
  if( markers[edge_id]==iteration_id ) 
  {
    return;
  }
  UpdateConductiveAtNodeVenous(upstream_node, edge_id);
}

void ConductiveTransport::UpdateConductiveAtNodeVenous(int upstream_node, int edge_id)
{

#if 1
  // assume upstream segments have already their conductive updated
  // that is segments adjacent to upstream_node, which are not edge_id
  
  if( nw.isCapillary[edge_id] ) //this means we are at the capillary
  {
    return;
  }
  int nbcnt = GetNeighborCount(upstream_node);
  if( nbcnt==1 )// dead end
  {
    SetConductive(edge_id,0.0);
  }
  else if(nbcnt<=16)
  {
    int numIn=0,numOut=0;
    FlReal flowIn=0;
    FlReal flowOut=0;
    int vin[3];
    int vout[3];
    FlReal conducGoingDown =0;
    ClearMem(vin, 3);
    ClearMem(vout, 3);
    for( int i=0; i<nbcnt; ++i ) 
    {
      const int  iw = GetNeighbor(upstream_node,i).edge;
      const my::eqpair<int> ee = GetEdge(iw);
      const int  imd = ee[0]==upstream_node ? ee[1] : ee[0];
      myAssert(imd == GetNeighbor(upstream_node, i).node);
      if( GetPress(imd) < GetPress(upstream_node) )
      { // outflow
        flowOut += GetFlow(iw);
        vout[numOut++] = iw;     
      } 
      else 
      {
	conducGoingDown += GetConductive(iw) * exp(-GetLen(iw)/adaptionParams.cond_length) + GetMetabolic(iw);
        flowIn += GetFlow(iw);
        vin[numIn++] = iw;
      }
    }
    if( numOut==1 && numIn==2 )
    {
      //myAssert( vin[0]==edge_id );
      const FlReal S_c = conducGoingDown;
      SetConductive(edge_id, S_c);
      markers[vout[0]] = iteration_id;
    }
    else if( numOut==2 && numIn==1 )
    {
      const int v1 = vout[0];
      const int v2 = vout[1];
      const int w = vin[0];
      SetConductive(v1,conducGoingDown);
      SetConductive(v2,conducGoingDown);
      markers[v1] = iteration_id;
      markers[v2] = iteration_id;
    }
    else if(numOut==1 && numIn==1)
    {
      SetConductive(vout[0], conducGoingDown);
      //SetHematocrit(vout[0], GetHematocrit(vin[0]));
      markers[vout[0]] = iteration_id;
    }
  }
  else // if num outgoing vessels > 2 or num ingoing vessel > 1
  {
    printf("BAD ERROR at venous Tree, \n");
    printf("Original vessel no.: %i \n", nw.newEdge2original[edge_id]);
    throw std::exception();
    myAssert(nbcnt<4);
  }
  markers[edge_id] = iteration_id;
#endif
}

void ConductiveTransport::Calculate_S_tot()
{
 /*
  * Pries, Secomb
  * Structural adaption and stability of microvascular networks:
  * theory and simulations
  * H354 equation 4
  */
  //FlReal s_tot_sq=0.0;
  //nw.max_delta_r = std::numeric_limits< FlReal >::min();
  //nw.max_stot = std::numeric_limits< FlReal >::min();
  this->nw.S_tot_array.resize(this->nw.num_edges());
  FlReal max_t1 = std::numeric_limits< FlReal >::min();
  FlReal min_t1 = std::numeric_limits< FlReal >::max();
  FlReal max_t2 = std::numeric_limits< FlReal >::min();
  FlReal min_t2 = std::numeric_limits< FlReal >::max();
  FlReal max_t3 = 0.0;
  FlReal max_t4 = 0.0;
  for(int i=0;i<this->nw.num_edges();++i)
  {    
    //v->f is in kpa
    double t1=-4.0;
//    t1 = log10(nw.shearforce[i]*10000.0+0.1); //done like this in the McDougall Anderson paper
    if( nw.shearforce[i]*10000.0 <0.0001)
    {
      //t1 = -2.0; // =log10(0.01);
      //t1 = -4.0; // =log10(0.0001);
      //t1=0.0;
    }
    else
    {
      t1 = log10(nw.shearforce[i]*10000.0);// 1kpa = 10 000 dyne/cm^2
    }
    // Attention: this is assuming that all pressures are positve >0
    const my::eqpair<int> ee = GetEdge(i);
    double t2 = (nw.press[ee[0]]+nw.press[ee[1]])/2.;
    t2 = t2 *7.50061683;// from kpa to mmhg
    myAssert(t2>=10.);
  #ifdef DEBUG
    //printf("P upstream: %f, P downstream: %f \n", GetUpstreamNode(v)->press, GetDownstreamNode(v)->press);
    printf("vessel id: %i, P: %f, log10(P): %f, log10(log10(P)): %f \n", i, t2, log10(t2), log10(log10(t2)));
  #endif
    //formula works only for pressures >10mmHg
    
    if(t2>10)
    {
      //logarithm of this will be calculated
      t2=log10(log10(t2));
      t2=-5000.0 * pow(t2,5.4);
      t2=exp(t2);
      t2=100-86.0*t2;
      t2=-log10(t2);
    }
    else
    {
      //t2=0.;
      t2=-log10(14.);
    }
    double t3 = GetMetabolic(i);
    double t4 = adaptionParams.k_c*(GetConductive(i)/(GetConductive(i)+adaptionParams.S_0));
    SetCondSignal(i,t4);
    double t5 = -adaptionParams.k_s;
    
    myAssert(!std::isinf(t1));myAssert(!std::isnan(t1));
    myAssert(!std::isinf(t2));myAssert(!std::isnan(t2));
    myAssert(!std::isinf(t3));myAssert(!std::isnan(t3));
    myAssert(!std::isinf(t4));myAssert(!std::isnan(t4));
    myAssert(!std::isinf(t5));myAssert(!std::isnan(t5));
    
    if(adaptionParams.tum_manitulate_s1 and nw.inTumor[i])
    {
      t1=1.3*t1;
    }
    if(adaptionParams.tum_manitulate_s2 and nw.inTumor[i])
    {
      t2=t2;
    }
    if(adaptionParams.tum_manitulate_s3 and nw.inTumor[i])
    {
      t3=0.7*t3;
    }
    if(adaptionParams.tum_manitulate_s4 and nw.inTumor[i])
    {
      t4=0.;
    }
    if(adaptionParams.tum_manitulate_s5 and nw.inTumor[i])
    {
      t5=t5*1.3;
    }
    if( t1 > max_t1) max_t1 = t1;
    if( t1 < min_t1) min_t1 = t1;
    if( t2 > max_t2) max_t2 = t2;
    if( t2 < min_t2) min_t2 = t2;
    if( t3 > max_t3) max_t3 = t3;
    if( t4 > max_t4) max_t4 = t4;
    double s_tot_this_edge = t1+t2+t3+t4+t5;
    
    myAssert(std::isfinite(s_tot_this_edge));
#ifdef DEBUG
    if(true){
      printf("t1: %f, t2: %f, t3: %f, t4: %f, t5: %f, S_tot %f\n",t1,t2,t3,t4,t5, s_tot_this_edge);
    }
#endif
    SetS_tot(i,s_tot_this_edge);
//    double delta_r=s_tot_this_edge*GetRadius(i);
//    delta_r = fabs(delta_r);
//    s_tot_this_edge = fabs(s_tot_this_edge);
    
//    s_tot_sq = s_tot_sq + s_tot_this_edge*s_tot_this_edge;
//     if( delta_r > nw.max_delta_r)
//     {
//       nw.max_delta_r = delta_r;
//     }
//     if( s_tot_this_edge > nw.max_stot )
//     {
//       nw.max_stot = s_tot_this_edge;
//     }
  }
#ifndef SILENT
  printf("->t1 min: %f, max: %f\n", min_t1,max_t1);
  printf("->t2 min: %f, max: %f\n", min_t2,max_t2);
  printf("max_t3: %f, max_t4: %f\n",max_t3,max_t4);
#endif
//  set_S_tot_sq(s_tot_sq);
}


/*
 * Adaption According to  Pries, Secomb, Gaehtgens 1998
 * "Structural adaption and stability of microvascular networks: theory and simulations"
 */
std::tuple<FlReal,FlReal,FlReal> CalcRadiiChange_mw(const Adaption::Parameters &params, VesselList3d &vl, float delta_t_calc)
{
  // generate a CompressedAdaptionNetwork structure
  CompressedAdaptionNetwork adaption_network;
  // read in the needed network stuff
  GetAdaptionNetwork(adaption_network, vl);
  
#ifdef DEBUG 
  //if debug, we can look at this data
  for (int i=0;i<adaption_network.edges.size();++i)
  {
    int i_a = adaption_network.edges[i].first;
    int i_b = adaption_network.edges[i].second;
#ifndef SILENT
    printf("no: %i, a: %i, b: %i\n",i,i_a,i_b);
    printf("press at %i: %f, press at %i: %f\n",i_a, adaption_network.press[i_a], i_b, adaption_network.press[i_b]);
#endif
  }
#endif
  //set up the conductiveTransport
  ConductiveTransport conductiveTransport(adaption_network,params);
  conductiveTransport.UpdateConductive();
  
  //here the work is done!!!
  conductiveTransport.Calculate_S_tot();
  
  //set the results back to the vessellist!
  SetAdaptionValues(vl, adaption_network, delta_t_calc, params.radMin_for_kill);
  //return minimization value
  std::tuple<double,double,double> myreturn(adaption_network.delta_r_square, adaption_network.max_stot, adaption_network.max_delta_r);
  return myreturn;
}
void KillSmallVessels(VesselList3d &vl, double rad_min)
{
  int ecnt = vl.GetECount();
  #define VESSEL_THREAD_CHUNK_SIZE 1024
  
  DynArray<Vessel*> toKill;
  tbb::spin_mutex mutex;
  #pragma omp parallel
  {
    DynArray<Vessel*> th_toKill(4024, ConsTags::RESERVE);
    #pragma omp for schedule(dynamic, VESSEL_THREAD_CHUNK_SIZE)
    for( int i=0;i<ecnt; ++i )
    {
      Vessel* v = vl.GetEdge( i );
      bool bKill = false;
      if( v->timeSprout>=0 ) continue;
      if(v->r < rad_min || !v->IsCirculated())
      {
	bKill = true;
      }
      if(bKill) 
      {
	th_toKill.push_back(v);//I think this is the crucial part, were threadsafety of boost is required!!!
	//note wrong: each thread get his own vector and then they are added
      }
    }
    mutex.lock();
    toKill.insert(toKill.end(), th_toKill.begin(), th_toKill.end());
    mutex.unlock();
    th_toKill.remove_all();
  }//end omp parallel
  #pragma omp barrier

  for( int i=0; i<toKill.size(); ++i )
  {
    vl.DeleteVessel( toKill[i] );
  }
}
void SetAdaptionValues(VesselList3d &vl, CompressedAdaptionNetwork& fl, double delta_t, double rad_min)
{
  
  int ecnt = vl.GetECount();
  int ncnt = vl.GetNCount();
  
  int negativeRadiusSuggested =0;
  int tooSmallRadiusSuggested =0;
  int tooBigRadiusSuggested =0;
  #define VESSEL_THREAD_CHUNK_SIZE 1024
  
  DynArray<Vessel*> toKill;
  tbb::spin_mutex mutex;

  FlReal max_delta_r_again = 0.0;
  FlReal max_stot_again = 0.0;
  FlReal delta_r_square = 0.0;
  
//  int k =0; //note definition of kk is even more powerful than a threadsafe counter!!!!
//#pragma omp parallel private(k) // this is needed because otherwise k++ not threadsafe

  #pragma omp parallel reduction(+:negativeRadiusSuggested,tooSmallRadiusSuggested,tooBigRadiusSuggested,delta_r_square),reduction(max:max_delta_r_again,max_stot_again)
  {
//     DynArray<Vessel*> th_toKill(4024, RESERVE);
//     DynArray<FlReal> th_max_delta_r(64,RESERVE);
//     DynArray<FlReal> th_max_stot(64,RESERVE);
    #pragma omp for schedule(dynamic, VESSEL_THREAD_CHUNK_SIZE)
    for( int i=0;i<vl.GetECount(); ++i )
    {
      Vessel* v = vl.GetEdge( i );
      int kk = fl.org2newEdge[v->Index()];
//       bool bKill = false;
//       if(!v->IsCirculated())  // moved to KillSmallVessels()
//       {
// // 	bKill = true;
// 	continue;
//       }
//       if(v->timeSprout>0)//do not kill the new created sprouts
//       {
// // 	bKill = false;
// 	continue;//ignore sprouts for adaptation
//       }
#if 0
      if(v->IsCirculated() and !v->flags.GetBits(BOUNDARY)
	or v->IsCirculated() and v->flags.GetBits(BOUNDARY) and v->IsVein()
      ) //is it important to include boundaries?
#endif
#if 1
      if( v->IsCirculated() )
#endif
      {
#ifdef DEBUG
	if(kk<0)
	  printf("kk: %i\n", kk);
#endif
	//node map
	int a = fl.org2new_vertex[v->NodeA()->Index()];
	int b = fl.org2new_vertex[v->NodeB()->Index()];
	myAssert(a != IntegerMap::unused() && b != IntegerMap::unused());
	myAssert(kk<fl.num_edges());
	v->S_total = fl.S_tot_array[kk];
	double aconductiveSignal = fl.condSignal[kk];
	myAssert( std::isfinite(aconductiveSignal));
	v->conductivitySignal = aconductiveSignal;
	v->metabolicSignal = fl.metabolic[kk];
  #if 1 //activate this in order to leafe out the change in radii
	double buff = v->r * fl.S_tot_array[kk];
	v->r = v->r + buff*delta_t;
	delta_r_square = delta_r_square + buff*buff;
	{
	  //if we change here, we store the amount
	  if(fabs(buff)>max_delta_r_again)
	  {
	    //max_delta_r_again = fabs(v->r * fl.S_tot_array[kk]);
	    max_delta_r_again = fabs(buff);
	  }
	  if(fabs(fl.S_tot_array[kk])>max_stot_again)
	  {
	    max_stot_again = fabs(fl.S_tot_array[kk]);
	  }
	}
	if( v->r <=rad_min )//minimal value
	{
	  tooSmallRadiusSuggested++;
	  if( v->r <=0 )
	  {
	    negativeRadiusSuggested++;
	  }
	  v->r = rad_min;
	  //uncomment this for not killing
	  //bKill = true;
	}
	if( v->r>142.)//maximun value
	{
	  v->r=142.;
	  tooBigRadiusSuggested++;
	}

  #endif
      }//end if circulated
    }//end parallel for all the vessels

    //parallel kill ;-)
//     mutex.lock();
//     //ensures that a vessel is killed only once!
//     toKill.insert(toKill.end(), th_toKill.begin(), th_toKill.end());
//     mutex.unlock();
//  printf("max_stot_again: %f\n", max_stot_again);
  }//end omp parallel
  #pragma omp barrier
#if 0
  printf("########## max delta r: %f \n", max_delta_r_again);
  printf("########## max s_tot : %f \n", max_stot_again);
#endif
  fl.max_delta_r = max_delta_r_again;
  fl.max_stot = max_stot_again;
  fl.delta_r_square = delta_r_square;
//  printf("toKill.size(): %i, vl.GetECount(): %i\n", toKill.size(),vl->GetECount());
  
//   for( int i=0; i<toKill.size(); ++i )
//   {
//     vl->DeleteVessel( toKill[i] );
//     //printf("toKill[i]: %s, toKill[i]->Index(): %i\n", toKill[i],toKill[i]->Index());
//   }
#ifndef SILENT
  if(negativeRadiusSuggested>0)
  {
    printf("WARNING: %i negative radii suggested!\n", negativeRadiusSuggested);
  }
  if(tooSmallRadiusSuggested>0)
  {
    printf("WARNING: %i too small radii suggested!\n", tooSmallRadiusSuggested);
  }
  if(tooBigRadiusSuggested>0)
  {
    printf("WARNING: %i too big radii suggested!\n", tooBigRadiusSuggested);
  }
#endif
}

#if 1
void GetAdaptionNetwork(
  CompressedAdaptionNetwork& fl, 
  VesselList3d& vl
 		      )
{
  int ecnt = vl.GetECount();
  int ncnt = vl.GetNCount();

  fl.edges.reserve(ecnt);
  fl.newEdge2original.resize(ecnt);
  fl.org2newEdge.resize(ecnt);
  fl.org2new_vertex.resize(ncnt);

  // copy only perfused edges
  for(int i=0; i<ecnt; i++)
  {
    const Vessel* v = vl.GetEdge(i);
    if (!v->IsCirculated()) continue;
    //edges
    int kk = fl.newEdge2original.add(i);
    int c = fl.org2newEdge.add(v->Index());
    //nodes
    int a = fl.org2new_vertex.add(v->NodeA()->Index());
    int b = fl.org2new_vertex.add(v->NodeB()->Index());
    myAssert(a!=b);
    myAssert(v->NodeA()->press!=0.0);
    myAssert(v->NodeB()->press!=0.0);
    //myAssert(v->NodeA()->press!=v->NodeB()->press);
    fl.edges.push_back(my::make_eqpair(a, b));
  }

  myAssert(fl.edges.size() > 0);

  // copy bcs and map index numbers
  auto mapKey = [&](const VesselNode* nd) -> int 
  { 
    return fl.org2new_vertex[nd->Index()];
  };
  remap_keys(vl.GetBCMap(), fl.bcs, mapKey);

#if 0 //useless vl->BCMap is used for calcflow
  for(int i=0; i<ncnt; ++i)
  {
    const VesselNode* vc= vl->GetNode(i);
    // insert boundary nodes into the bcs array
    if (vc->IsBoundary())
    {
      int id = fl.org2new_vertex[vc->Index()];
      if (id != IntegerMap::unused() && fl.bcs.find(id) == fl.bcs.end())
      {
	if(vc->GetEdge(0)->IsArtery())
	{
	  fl.bcs[id] = FlowBC(FlowBC::CURRENT, vc->GetEdge(0)->q);
	}
	if(vc->GetEdge(0)->IsVein())
	{
	  fl.bcs[id] = FlowBC(FlowBC::PIN, vc->press);
	}
	//cout<<format("press : %f\n") % vc->press;
        //cout << format("flow boundary node: %i") % id << endl;
      }
    }
  }
#endif
  {
    // copy edge properties
    fl.len.resize(fl.edges.size());
    fl.hema.resize(fl.edges.size());
    fl.flow.resize(fl.edges.size());
    fl.shearforce.resize(fl.edges.size());
    fl.isArtery.resize(fl.edges.size());
    fl.isVein.resize(fl.edges.size());
    fl.isCapillary.resize(fl.edges.size());
    fl.inTumor.resize(fl.edges.size());
    fl.isBoundaryVessel.resize(fl.edges.size());
    fl.rad.resize(fl.edges.size());
    /*
     * vertices loop
     */
    fl.press.resize(fl.num_vertices());
#ifdef DEBUG
#ifndef SILENT
    for( int i=0; i<vl.GetNCount();++i)
    {
      printf("vl press: %f\n", vl.GetNode(i)->press);
    }
#endif
#endif
    for(int i=0; i<ncnt;++i)
    {
      if (fl.org2new_vertex[i] != IntegerMap::unused())
      {
	const VesselNode *nd = vl.GetNode(i);
	fl.press[fl.org2new_vertex[i]] = nd->press;
      }
    }
    /*
     * edges loop
     */
    for(int i=0, k=0; i<ecnt; ++i)
    {
      const Vessel* v = vl.GetEdge(i);
      if (!v->IsCirculated()) continue;
      fl.hema[k] = v->hematocrit;
      fl.flow[k] = v->q*60./1e6;
      fl.shearforce[k] = v->f;
      fl.isArtery[k] = v->IsArtery();
      fl.isVein[k] = v->IsVein();
      fl.isCapillary[k] = v->IsCapillary();
      fl.inTumor[k] = v->flags.GetBits(WITHIN_TUMOR);
      fl.isBoundaryVessel[k] = v->flags.GetBits(BOUNDARY);
      fl.rad[k] = v->r;
      double wl=0;
      if (!vl.HasLattice())
      {
        wl = (v->NodeA()->worldpos.transpose()-v->NodeB()->worldpos.transpose()).norm();
	myAssert(wl>0);
#ifdef DEBUG
#if ADAPTION_OUTPUT
cout<<"vessel no.: "<<v->Index()<<" length: "<<wl<<endl;
cout<<"posA : "<<v->NodeA()->worldpos<<endl;
cout<<"posB : "<<v->NodeB()->worldpos<<endl;
#endif
#endif
      }
      else
      {
	const VesselList3d::LatticeData &ld = vl.Ld();
	wl = v->WorldLength(ld);
      }
      fl.len[k] = wl;
      ++k;
    }
  }
}
#endif
bool check_while_break( int no_vessels, double min_error, double nqdev, double max_delta_r, double qdev_from_params, double max_stot)
{
  double stopp = min_error/sqrt(no_vessels);
  //global conditions
  if (qdev_from_params >0)
  {
    if( nqdev < stopp)
      return false;
  }
  //local conditions
  else
  {
#if 0
    if (nqdev < stopp*0.42)
    {
      return false;
    }
    else
    {
    return max_delta_r>1.0;
    }
#endif
#if 0
    if( max_delta_r < 1.0)
      return false;
    if( max_stot < 0.1)
      return false;
#endif
#if 1
    //printf("max_stot: %f\n", max_stot);
    if( max_stot < 0.01)
      return false;
#endif
  }
  // if not in one case from above, continue while loop
  return true;
}


std::tuple<uint,FlReal> runAdaption_Loop( Parameters params, BloodFlowParameters bfparams, VesselList3d &vl,bool doDebugOutput)
{
  /*
   * first, change the boundary Conditions
   * by default vesselgen only creates pressure boundary conditions,
   * it turns out that they are not so nice for adaption
   */
  //is good to have well definde pressures and flow, even if this may be redundant
  //VesselList3d *p_vl = &vl;
  //std::auto_ptr<VesselList3d> p_vl(&vl);
#ifdef DEBUG
  printf("running first CalcFlow in adaptionLoop\n");
  CalcFlow(vl, bfparams);
  printf("done!\n");
#endif
#if 1
  //this seems to be necessary to get a stable result for the apj stuff
  //used together with change of radii of veins
  // is this threadsafe when multiple instances of p_vl exist?
  #pragma omp single
  {
    ChangeBoundaryConditions(vl,params);
  }
  
#endif
  CalcFlow(vl, bfparams);
  int no_Vessels_before_adaption = vl.GetECount();
  //adaption loop local variables
  FlReal qdev=std::numeric_limits<FlReal>::max();
  FlReal nqdev=std::numeric_limits<FlReal>::max();
  FlReal max_delta_r = std::numeric_limits<FlReal>::max();
  FlReal max_stot = std::numeric_limits<FlReal>::min();
  FlReal delta_t = params.delta_t;
  bool use_dynamic_t = false;
  if (delta_t == 0)
  {
    use_dynamic_t = true;
    delta_t = 0.1; //as initial guess
  }
    
  int level=0;
  std::string temp;
  /*
   * some data we can look at after convergence, container for local variables
   */
  int how_often=0;//counts iterations untils convergence and can cause failsafe stop
  DynArray<FlReal> last_vessel_radii(vl.GetECount());
  std::vector<FlReal> qdevs;
  std::vector<FlReal> nqdevs;
  std::vector<FlReal> max_delta_r_s;
  std::vector<FlReal> delta_t_s;
  
  //we label all small vessels as capillaries, maybe this this helps for convergence
  //no, this makes it even worse
//   for(int i=0;i<vl->GetECount();++i)
//   {
//     Vessel *v = vl->GetEdge(i);
//     if (v->r < 4.2 ) //note: why 42???
//     {
//       v->flags.AddBits(CAPILLARY);
//     }
//   }
    
   
  /*
   * this is the MAIN loop!!!!!!!!!!
   */
  //lot of crazy attemps
  //max_delta_r=1.;
  //while(max_delta_r>params->qdev)
  //while(nqdev > params->qdev)
  //while(nqdev>params->qdev and max_delta_r >0.05)
  //while(max_delta_r >0.042)
  //double stopp = 1./vl->GetECount();
  //while(max_stot >0.042 or how_often<1)
  //while(max_delta_r > params->qdev)
  
  
  
  /*
   * stopping condition is dynamical now,
   * therefore print it before we start.
   */
  double limit_error = 0.001;//in mu m
  double stopp = 0;
#ifndef SILENT
  printf("limit_error: %f, stopp: %f\n", limit_error,stopp);
  printf("MAX IT: %i, break at: %f\n", params->max_nun_iterations,stopp);
#endif
  bool mybreakcondition= true;
  while(mybreakcondition)
  //while(nqdev > stopp)
  //while(max_delta_r>2.0)
  //while(max_delta_r>2.5 or max_stot >0.2)
  {
    how_often++;
    //break if we reach max_nun_iterations or for failsafe more than 200K steps
    if(params.max_nun_iterations>0)
    {
      if(how_often>30000 or how_often> params.max_nun_iterations-1)
      {
      break;
      }
    }
    
    //ExtremeFlows myExtremFlows = GetExtremFlows(p_vl);
    //GetExtremFlows(vl, &myExtremFlows);
    //params.Q_refdot = (params.max_flow * params.min_flow)/params.no_of_roots;
#ifndef SILENT
    //cout << "no_of_roots: " << params.no_of_roots << endl;
    cout << "max_flow: " << myExtremFlows.max_flow << endl;
    cout << "min_flow: " << myExtremFlows.min_flow << endl;
    cout << "min_flow*max_flow: " << myExtremFlows.min_flow*myExtremFlows.max_flow  << endl;
    cout << "arithmetica avg: " << (myExtremFlows.min_flow + myExtremFlows.max_flow)/2 <<endl;
    cout << "avg_flow: " << myExtremFlows.avg_flow <<endl;
#endif

    if( false )//this is a little like simulated annealing
    {
      if(params.max_nun_iterations>0 and (how_often/params.max_nun_iterations>0.5 and level == 0))
      {
	delta_t = 0.1*delta_t;
	level++;
      }
      if(params.max_nun_iterations>0 and (how_often/params.max_nun_iterations>0.8 and level == 1))
      {
	delta_t = 0.1*delta_t;
	level++;
      }
    }
    
    // **************MAIN Routine ******************** for single iteration
    std::tie(qdev,max_stot,max_delta_r) = Adaption::CalcRadiiChange_mw(params, vl,delta_t);
#if 1
    // am I serious? this would cause a never ending story.
    //UpdateBoundaryConditions(*vl);//new boundary flows, can show up after adaption
    /* 
     * segfault at trilinos? AZ_manage_memory?
     */
    CalcFlow(vl, bfparams);
    //flow infos will be passed as return values, so better use new ones
    //vice versa, network can expire different flows because of change in BC
    //an change in radii
#endif
    qdev = sqrt(qdev);
    nqdev = qdev/vl.GetECount();
    //qdev = Adaption::CalcRadiiChange_RK4(params,*vl);
    qdevs.push_back(qdev);
    nqdevs.push_back(nqdev);
    max_delta_r_s.push_back(max_delta_r);
    
    //update stopping condition
    stopp = limit_error/sqrt(vl.GetECount());
    
    //printf("hack delta t!");
    delta_t_s.push_back(delta_t);
    /*
     * dynamic change of stepsize, approx 1/largest eigenvalue
     */
    if (use_dynamic_t)
    {
      delta_t = 1/(4*max_stot);
    }
    
    
    //ComputeCirculatedComponents(vl); //done also in calcflow!!!!
#if ADAPTION_OUTPUT
#ifndef SILENT
    printf("Iterated adaption #%i, nqdev: %f, max_stot: %f, max_delta_r: %f\n", how_often, nqdev, max_stot, max_delta_r);
#endif
#ifdef SILENT
    if( how_often%200 == 0)
    {
      printf("                                                        Iterated adaption #%i, nqdev: %f/%f, max_stot: %f, max_delta_r: %f\n", how_often, nqdev,stopp, max_stot,max_delta_r);
      if (use_dynamic_t) printf("dynamic t change: t= %f!\n",delta_t);
      printf("limit_error: %f, stopp: %f\n", limit_error,stopp);
    }
#endif
#endif

    
    //debug output everystep
#if 0
#ifdef DEBUG
  if(doDebugOutput)
  {
    if( vessels_after_adaption != nullptr )
    {
      temp = "vessels_after_adaption_" + std::to_string(how_often);
      h5cpp::Group grp_temp = vessels_after_adaption->create_group(temp);
      //h5cpp::Group grp_temp = the_debug_file.root().create_group(temp);
      //h5cpp::Group grp_temp = f.create_group(temp);
      ptree getEverytingPossible = make_ptree("w_all", true);
      WriteVesselList3d(*vl, grp_temp, getEverytingPossible);
      h5cpp::Dataset dsqdevs = h5cpp::create_dataset(grp_temp, "qdev", h5cpp::Dataspace::simple_dims(qdevs.size()), &qdevs[0]);
      h5cpp::Dataset dsnqdevs = h5cpp::create_dataset(grp_temp, "nqdev", h5cpp::Dataspace::simple_dims(nqdevs.size()), &nqdevs[0]);
      }
    }
#endif //DEBUG
#endif
  mybreakcondition = check_while_break(vl.GetECount(), limit_error,nqdev,max_delta_r, params.qdev, max_stot);
  }
  //end main adaption loop, hopefully convergent here
  
  
  KillSmallVessels(vl,params.radMin_for_kill);
  
  
  int no_Vessels_after_adaption = vl.GetECount();
  int no_killed = no_Vessels_before_adaption-no_Vessels_after_adaption;
  if( no_killed >0)
    printf("Killed %i of %i Vessels below %f\n", no_Vessels_before_adaption-no_Vessels_after_adaption, no_Vessels_before_adaption,params.radMin_for_kill);

#if 0
#ifndef DEBUG //release output, only last convergent part
  if(vessels_after_adaption != nullptr)
  {
    if(params->write2File)
    {
      std::string temp = "vessels_after_adaption";
      //temp = "vessels_after_adaption";
      h5cpp::Group grp_temp = vessels_after_adaption->create_group(temp);
      ptree getEverytingPossible = make_ptree("w_all", true);
      WriteVesselList3d(*vl, grp_temp, getEverytingPossible);
      h5cpp::Dataset dsqdevs = h5cpp::create_dataset(grp_temp, "qdev", h5cpp::Dataspace::simple_dims(qdevs.size()), &qdevs[0]);
      h5cpp::Dataset dsnqdevs = h5cpp::create_dataset(grp_temp, "nqdev", h5cpp::Dataspace::simple_dims(nqdevs.size()), &nqdevs[0]);
      h5cpp::Dataset dsmaxdelta_r = h5cpp::create_dataset(grp_temp, "max_delta_r", h5cpp::Dataspace::simple_dims(max_delta_r_s.size()), &max_delta_r_s[0]);
      h5cpp::Dataset dsdelta_t = h5cpp::create_dataset(grp_temp, "delta_t", h5cpp::Dataspace::simple_dims(delta_t_s.size()), &delta_t_s[0]);
    }
  }
#endif
#endif
  
 
  //write the input and output radii to console  
  for(int w=0;w<vl.GetECount();++w)
  {
    const Vessel *v = vl.GetEdge(w);
#if ADAPTION_OUTPUT
#ifdef DEBUG
    printf("change at #%i: %f and %f \n",w, last_vessel_radii[w],v->r);
#endif
#endif
  }
#if ADAPTION_OUTPUT
  cout<<format(" nqdev: %f, max_stot: %f, max_delta_r: %f # %i iterations\n ") % nqdev %max_stot %max_delta_r %how_often;
#endif
  /*
   * return statement   0: convergent
   * 			1: not convergent */
  FlReal mean_value; 
  FlReal mean_std;
  if(how_often >= params.max_nun_iterations)
  {
#if ADAPTION_OUTPUT
    cout<<"NOT Convergent!"<<endl;
#endif
    return std::make_tuple(1, 1000000000.);
  }
  else
  {
#if ADAPTION_OUTPUT
    cout<<"Convergent!"<<endl;
#endif
    using namespace boost::accumulators;
    accumulator_set<double, features<tag::mean, tag::variance>> acc;
#pragma omp parallel for
  for(int i =0;i<vl.GetECount();++i)
  {
    Vessel* v = vl.GetEdge(i);
    acc(v->q);
  }
  mean_value = mean(acc);
  mean_std = sqrt(variance(acc));
#pragma omp barrier
    
    return std::make_tuple(0, mean_value);
  }
}

Adaption::ExtremeFlows GetExtremFlows(VesselList3d *vl_)
{
  Adaption::ExtremeFlows *theExtrems = new ExtremeFlows();
  int n = vl_->GetECount();
  double max_flow=0.0;
  double min_flow=std::numeric_limits<double>::max();
  double avg_flow=0.0;
  for(int i=0;i<vl_->GetECount();++i)
  {
    const Vessel* v = vl_->GetEdge(i);
    if(v->q>max_flow and v->IsCirculated())max_flow=v->q;
    if(v->q<min_flow and v->IsCirculated())min_flow=v->q;
    if(v->IsCirculated()) avg_flow = avg_flow + v->q;
  }
#ifdef DEBUG
  myAssert(n>0);
#endif
  theExtrems->avg_flow = avg_flow / n;
  theExtrems->max_flow = avg_flow*60/1000000.;
  theExtrems->max_flow = max_flow*60/1000000.;
  theExtrems->min_flow = min_flow*60/1000000.;
  
  return *theExtrems;
}
#if 0
void UpdateBoundaryConditions(VesselList3d &vl)//after radii are changed!!!
{
  int ncnt = vl.GetNCount();
  double min_boundary_pressure = std::numeric_limits< double >::max();
  double max_boundary_flow = std::numeric_limits< double >::min();
//   boost::container::vector<VesselNode*> th_b_nd_container;
//   boost::container::vector<double> th_min_b_p_container;
  DynArray<VesselNode*> th_b_nd_container;
  DynArray<VesselNode*>th_min_b_p_container;
  for(auto itr = vl.GetBCMap().begin();itr != vl.GetBCMap().end();++itr)
  {
    const VesselNode* nd = itr->first;
    FlowBC aflowbc = itr->second;
    if(aflowbc.type == FlowBC::PIN)
    {
      vl.SetBC(nd,FlowBC(FlowBC::PIN, vl.GetNode(nd->Index())->press));
    }
    if(aflowbc.type == FlowBC::CURRENT and nd->GetEdge(0)->IsVein())
    {
      vl.SetBC(nd,FlowBC(FlowBC::CURRENT, nd->GetEdge(0)->q));
    }
    if(aflowbc.type == FlowBC::CURRENT and nd->GetEdge(0)->IsArtery())
    {
      vl.SetBC(nd,FlowBC(FlowBC::CURRENT, -nd->GetEdge(0)->q));
    }

  }
#if 0
#pragma omp for
  for( int i=0;i<ncnt; ++i)
  {
    const VesselNode* vc= vl.GetNode(i);
    if( vc->IsBoundary())
    {
      th_b_nd_container.push_back(vc);
      th_min_b_p_container.push_back(vc->press);
    }
  }
#pragma omp barrier
  for(int i=0; i<ncnt; ++i)
  {
    const VesselNode* vc= vl.GetNode(i);
    // insert boundary nodes into the bcs array
    if (vc->IsBoundary())
    {
#if 0  //keep max flow as flow boundary, everything else to pressure
      if(vc->GetEdge(0)->q== max_boundary_flow)
      {
	if(vc->GetEdge(0)->IsArtery())
	{
	  vl.SetBC(vc,FlowBC(FlowBC::CURRENT, -vc->GetEdge(0)->q));
	}
	if(vc->GetEdge(0)->IsVein())
	{
	  vl.SetBC(vc,FlowBC(FlowBC::CURRENT, vc->GetEdge(0)->q));
	}
      }
      else
      {
	vl.SetBC(vc,FlowBC(FlowBC::PIN, vc->press));
      }
#endif
#if 0
      //all arteries flow, all veins pressure
      if(vc->GetEdge(0)->IsArtery())
      {
	vl.SetBC(vc,FlowBC(FlowBC::CURRENT, -vc->GetEdge(0)->q));
      }
      if(vc->GetEdge(0)->IsVein())
      {
	vl.SetBC(vc,FlowBC(FlowBC::PIN, vc->press));
      }
#endif
#if 1
      //minmal boundary pressure becomes PIN condition everthing else flow
      if(vc->press == min_boundary_pressure)
      {
	vl.SetBC(vc,FlowBC(FlowBC::PIN, vc->press));
	//cout << format("min pressure boundary: %f") % vc->press <<endl;
      }
      else
      {
	if( vc->GetEdge(0)->IsArtery() )
	{
	  vl.SetBC(vc,FlowBC(FlowBC::CURRENT, -vc->GetEdge(0)->q));
	  //cout << format("flow boundary: %f") % -vc->GetEdge(0)->q <<endl;
	}
	if( vc->GetEdge(0)->IsVein())
	{
	  vl.SetBC(vc,FlowBC(FlowBC::CURRENT, vc->GetEdge(0)->q));
	  //cout << format("flow boundary: %f") % vc->GetEdge(0)->q <<endl;
	}
      }
#endif
      //cout<<format("press : %f\n") % vc->press;
      //cout << format("flow boundary node: %i") % id << endl;
      //cout<<"Changed bcs map!"<<endl;
    }
  }
#endif
}
#endif

}//namespace
