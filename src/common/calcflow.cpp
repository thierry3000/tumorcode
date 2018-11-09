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
//#define DOTIMING 1

#include "calcflow.h"
#include "vessels3d.h"
#include "shared-objects.h"

#include <boost/foreach.hpp>
#include <limits>
#include <fstream>
#include "Ifpack_Utils.h"
#include <mutex>

std::mutex calcflow_mutex;


// void ChangeBoundaryConditionsFix(VesselList3d &vl)
// {
// #ifdef DEBUG
//   printf("entered ChangeBoundaryConditionsFix\n");
// #endif
//   int ncnt = vl.GetNCount();
//   
//   for(int i=0; i<ncnt; ++i)
//   {
//     const VesselNode* vc= vl.GetNode(i);
//     // insert boundary nodes into the bcs array
//     if (vc->Count() > 0 and vc->IsBoundary() and vc->GetEdge(0)->IsCirculated() and vc->GetEdge(0)->IsVein())
//     {
//       vl.SetBC(vc,FlowBC(FlowBC::PIN, 2.0));
//     }
//     if (vc->Count() > 0 and vc->IsBoundary() and vc->GetEdge(0)->IsCirculated() and vc->GetEdge(0)->IsArtery())
//     {
//       vl.SetBC(vc,FlowBC(FlowBC::CURRENT, -5000000.0));
//     }
//       //cout<<format("press : %f\n") % vc->press;
//       //cout << format("flow boundary node: %i") % id << endl;
// #ifdef DEBUG
//       cout<<"Changed bcs map!"<<endl;
// #endif
//   }
// #ifdef DEBUG
//   printf("leave ChangeBoundaryConditionsFix\n");
// #endif
// }


template<class Flux>
double CalcFluxResidual(const CompressedFlowNetwork &flownet, Flux flux)
{
  FlArray absflux(flownet.num_vertices());
  FlArray nodeflux(flownet.num_vertices());
  for (int i=0; i<flownet.num_edges(); ++i)
  {
    int a, b; boost::tie(a,b) = flownet.edges[i];
    double q = flux(i, a, b);
    myAssert(q >= 0.);
    if (flownet.press[a] > flownet.press[b])
    { // a -> b
      nodeflux[a] -= q;
      nodeflux[b] += q;
    }
    else
    {
      nodeflux[a] += q;
      nodeflux[b] -= q;
    }
    absflux[a] += q;
    absflux[b] += q;
  }
  double res = 0.;
  for (int i=0; i<nodeflux.size(); ++i)
  {
    if (absflux[i]!=0. && flownet.bcs.find(i)==flownet.bcs.end())
      res = std::max(res, std::abs(nodeflux[i])/absflux[i]);
  }
  return res;
}


double CalcFlowResidual(const CompressedFlowNetwork &flownet, const FlArray &cond)
{
  auto flux = [&](int i, int a, int b) -> double
  {
    return cond[i]*std::abs(flownet.press[a]-flownet.press[b]);
  };
  return CalcFluxResidual(flownet, flux);
}


double CalcHemaResidual(const CompressedFlowNetwork &flownet, const FlArray &cond, const FlArray &hema)
{
  auto flux = [&](int i, int a, int b) -> double
  {
    return hema[i]*cond[i]*std::abs(flownet.press[a]-flownet.press[b]);
  };
  return CalcFluxResidual(flownet, flux);
}


bool MarkFailingVesselHidden(VesselList3d &vl)
{
  double max_res = 0.;
  bool ok = true;
  for (int i=0; i<vl.GetNCount(); ++i)
  {
    VesselNode* nd = vl.GetNode(i);
    nd->residual = 0.;
    
    if (nd->IsBoundary())
    {
      continue;
    }
    double qsum = 0.;
    double qabs = 0.;
    int n = 0;
    for (int i=0; i<nd->Count(); ++i)
    {
      auto nb = nd->GetNode(i);
      if (!nb.edge->IsCirculated())
      {
        continue;
      }
      if (nb.node->press < nd->press)
      {
        qsum -= nb.edge->q;
      }
      else if (nb.node->press > nd->press)
      {
        qsum += nb.edge->q;
      }
      else
      {
        myAssert(nb.edge->q == 0.);
      }
      qabs += nb.edge->q;
      n++;
    }
    if (qabs!=0)
    {
      double r = std::abs(qsum)/qabs;
      max_res = std::max(r, max_res);
      nd->residual = r;
      if (r >= 1.)
      {
        nd->flags |= HIDDEN;
        ok = false;
      }
    }
  }
#ifndef NDEBUG
#ifndef SILENT
  cout << "max relative fluid loss: " << max_res << endl;
#endif
#endif
  return ok;
}


template<class T>
T Median(const std::vector<T> &x)
{
  std::vector<T> v(x); // copy
  int n = v.size()/2;
  std::nth_element(v.begin(), v.begin()+n, v.end());
  return v[n];
}


/*------------------------------------------------------------------
--------------------------------------------------------------------*/
void CalcFlowSimple(VesselList3d &vl, const BloodFlowParameters &bloodFlowParameters, bool keepTheVesselHematocrit)
{
  {
#ifndef SILENT
  cout << "calcflow (no hematocrit)" << endl;
#endif
  CompressedFlowNetwork flownet;
  uint returnOfGetFlowNetwork = GetFlowNetwork(flownet, vl, bloodFlowParameters, keepTheVesselHematocrit);
  if(returnOfGetFlowNetwork>0)
  {
    return;
  }

  FlArray visc(flownet.num_edges());
  CalcViscosities(flownet.rad, flownet.hema,bloodFlowParameters, visc);
  //from header
  //FlArray press;
  //FlArray len, rad, hema, flow;
  //clear_and_free_memory(flownet.hema);
  
  FlArray cond(flownet.num_edges());

  CalcConductivities(flownet.rad, flownet.len, visc, cond);
  clear_and_free_memory(flownet.len);
  clear_and_free_memory(flownet.rad);
  
  
  Linsys flowsys;
  flowsys.scaling_const = CalcFlowCoeff(bloodFlowParameters.viscosityPlasma, 4., 100.); //Median(cond);
  flowsys.initialize_pattern(flownet.num_vertices(), flownet.edges);
  flowsys.fill_values(flownet.edges, cond, flownet.bcs);

  //sparse_matrix_print(*flowsys.sys);

  flowsys.solve();
  
  //vector_print(*flowsys.lhs);

  flownet.press.resize(flownet.num_vertices());
  //flownet.len.resize(flownet.num_edges());
  //flownet.rad.resize(flownet.num_edges());
  //flownet.hema.resize(flownet.num_edges());
  //flownet.flow.resize(flownet.num_edges());
  
  for (int i=0; i<flownet.num_vertices(); ++i) flownet.press[i] = flowsys.lhs_get(i);
  SetFlowValues(vl, flownet, cond, flownet.press, flownet.hema);
    
    
#if 0
    {
      Epetra_Vector residuals(flowsys.lhs->Map(), true);
      flowsys.sys->Multiply(false, *flowsys.lhs, residuals);
      residuals.Update(-1., *flowsys.rhs, 1.);
    
      std::ofstream f("failsol.txt");
      f << format("%i x %i matrix)") % flowsys.sys->NumGlobalRows() % flowsys.sys->NumGlobalCols() << endl;
      for (int i=0; i<flownet.num_vertices(); ++i)
      {
        int cnt;
        double *values;
        int *indices;
        flowsys.sys->ExtractMyRowView(i, cnt, values, indices);
        f << format("m[%i,...] = ") % i;
        for (int k=0; k<cnt; ++k) f << format("%10i ") % indices[k];
        f << "\n";
        f << "             ";
        for (int k=0; k<cnt; k++) f << format("%.10lf ") % values[k];
        f << "\n";
        //f << format(" lhs = %.20lf rhs = %.20lf\n") %  flowsys.lhs_get(i) % residuals[i];
      }
      f.flush();
      std::exit(0);
    }
#endif

#ifdef DEBUG
  for (auto iter = flownet.bcs.begin(); iter != flownet.bcs.end(); ++iter)
  {
    cout << format("boundary node id: %i, value = %f") % iter->first % iter->second.val << endl;
  }
#endif
  //cout << "max relative mass loss: " << CalcFlowResidual(flownet, cond) << endl;

  //H5::H5File f("debugvessels.h5","w");
  //WriteVesselList3d(*vl, f.root().create_group("vessels"), make_ptree("w_all",true));

  }
  bool ok = MarkFailingVesselHidden(vl);
  if (!ok)
  {
    ComputeCirculatedComponents(&vl);
    CalcFlowSimple(vl, bloodFlowParameters, keepTheVesselHematocrit);
  }
}



#if 1 // hematocrit calculations
/*------------------------------------------------------
// hematocrit calculations
------------------------------------------------------*/


class HematocritCalculator
{  
  CompressedFlowNetwork &nw; // stores edges and properties
  CompressedRows nbs;  // stores graph neighbor relations in sparse compressed column format
  const FlArray   &flow,&press; // more properties
  FlReal          h0; // initial hematocirt (from inflow)
  DynArray<uchar> markers; // marks nodes and edges which have been visited
  uchar             iteration_id; //current iteration

  // trivial getters and setters
  FlReal            GetFlow(int i) const                
  {
    myAssert(flow[i]>0.0);
    return flow[i]; 
  }
  FlReal            GetPress(int i) const               { return press[i]; }
  my::eqpair<int>   GetEdge(int i) const                { return nw.edges[i]; }
  bool              IsBoundary(int i) const             { return nw.bcs.find(i) != nw.bcs.end(); }
  FlReal            GetRad(int i) const                 { return nw.rad[i]; }
  FlReal            GetHematocrit(int i) const          { return nw.hema[i]; }
  void              SetHematocrit(int i,FlReal x) const;

  // get neighboring info i.e. node -> #neighbors, neighbor[k] -> edge, node
  int GetNeighborCount(int i) const { return nbs.row_length(i); }
  NodeNeighbor<int,int> GetNeighbor(int i, int k) const
  {
    int id_of_edge = nbs.column_number(i, k);
    const my::eqpair<int> &edge = nw.edges[id_of_edge];
    return make_node_neighbor(edge.first==i ? edge.second : edge.first, id_of_edge);
  }

  void UpdateHematocritRecursive(int iv, int ind_up);
  void UpdateHematocritAtNode(int ind, int iv);

public:
  bool check_hematocrit_range;
  HematocritCalculator(CompressedFlowNetwork &fl, const FlArray &flow, const FlArray &press, const FlReal h0);
  void UpdateHematocrit();
};

HematocritCalculator::HematocritCalculator(CompressedFlowNetwork &fl_, const FlArray &flow_, const FlArray &press_, const FlReal h0_) :
  nw(fl_), flow(flow_), press(press_), h0(h0_), iteration_id(0), check_hematocrit_range(false)
{
  int ecnt = nw.edges.size();
  int ncnt = nw.num_vertices();
  myAssert(flow.size()==ecnt);
  myAssert(press.size()==ncnt);
  myAssert(nw.rad.size()==ecnt);
  std::vector<unsigned int> col_counts(ncnt);
  for (int i=0; i<ecnt; ++i)
  {
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


void HematocritCalculator::SetHematocrit(int i,FlReal x) const
{
  myAssert(std::isfinite(x));
  if (!check_hematocrit_range) x = my::cut(x, 0., 0.9999999);
  else if (x<0. || x>1.01) 
  {
#ifndef TOTAL_SILENCE
    printf("warning hematocrit %f out of bounds at vessel %i\n", x, i);
#endif
    myAssert(x<0. || x>1.01); 
  }
  nw.hema[i]=x;
}


void HematocritCalculator::UpdateHematocrit()
{
  ++iteration_id;
  nw.hema.fill(-1.);
  for (int i=0; i<nw.num_edges(); ++i)
  {
    if (markers[i] == iteration_id) continue;
    const my::eqpair<int> &e = GetEdge(i);
    int upstream_node = (GetPress(e.first) > GetPress(e.second)) ? e.first : e.second;
    UpdateHematocritRecursive(i, upstream_node);
  }
}


void HematocritCalculator::UpdateHematocritRecursive( int edge_id, int upstream_node )
{
  // upstream_node is the node of edge_id which is more upstream than the other
  // we march through the vessel graph against the stream
  // at the upmost nodes the hematocrit is set to the initial value
  // when the function call stack is rewound all the hematocrit values are computed
  // in the correct order in downstream direction
  if( markers[edge_id]==iteration_id ) return;

  if( !IsBoundary(upstream_node) ) // compute upstream segments
  {
    for( int k=0; k<GetNeighborCount(upstream_node); ++k )
    {      
      NodeNeighbor<int,int> nb = GetNeighbor(upstream_node,k);
      if (GetPress(nb.node)<=GetPress(upstream_node)) continue; // do not go downstream
      UpdateHematocritRecursive( nb.edge, nb.node ); // this is an upstream edge, go further upstream
    }
  }
  if( markers[edge_id]==iteration_id ) return;
  UpdateHematocritAtNode(upstream_node, edge_id);
}


void HematocritCalculator::UpdateHematocritAtNode( int upstream_node, int edge_id )
{
  // assume upstream segments have already their hematocrit updated
  // that is segments adjacent to upstream_node, which are not edge_id
  
  if( IsBoundary(upstream_node) )
  {
    SetHematocrit(edge_id, h0);
    markers[edge_id] = iteration_id;
    return;
  }
  int nbcnt = GetNeighborCount(upstream_node);
  if( nbcnt==1 )
  {
    SetHematocrit(edge_id,0.0);
  }
#if 0
  else if( nbcnt==2 )
  {
    NodeNeighbor<int, int> parent = GetNeighbor(upstream_node,0);
    if( parent.edge==edge_id ) parent = GetNeighbor(upstream_node,1);
    myAssert(GetPress(parent.node)>=GetPress(upstream_node));
    myAssert( markers[parent.edge]==iteration_id );
    SetHematocrit(edge_id,GetHematocrit(parent.edge));
  }
#endif
  else if(nbcnt<=16)
  {
    //if( vc->Count()>2 ) {
    //  printf("WARNING: UpdateHematocritAtNode: #vessels at node (%i,%i) = %i\n",vc->lpos.x,vc->lpos.y,vc->Count());
    //}
    int numIn=0,numOut=0;
    FlReal flowIn=0;
    FlReal flowOut=0;
    FlReal rbcIn=0;
    int vin[16];
    int vout[16];
    FlReal rbcOut[16];
    ClearMem(vin, 16);
    ClearMem(vout, 16);
    ClearMem(rbcOut, 16);
    for( int i=0; i<nbcnt; ++i ) 
    {
      const int  iw = GetNeighbor(upstream_node,i).edge;
      const my::eqpair<int> ee = GetEdge(iw);
      const int  imd = ee[0]==upstream_node ? ee[1] : ee[0];
      myAssert(imd == GetNeighbor(upstream_node, i).node);
      if( GetPress(imd) > GetPress(upstream_node) )
      { // inflow
        flowIn += GetFlow(iw);
        rbcIn  += GetFlow(iw) * GetHematocrit(iw);
        vin[numIn++] = iw;
      } 
      else 
      {
        flowOut += GetFlow(iw);
        rbcOut[numOut] = 0.0;
        vout[numOut++] = iw;
      }
    }
    
    if( numOut==1 )
    {
      myAssert( vout[0]==edge_id );
      const FlReal h = rbcIn/GetFlow(edge_id);
      SetHematocrit(edge_id, h);
    }
    else if( numOut==2 && numIn==1 )
    {
      const int w= vin[0];
      const int v1 = vout[0];
      const int v2 = vout[1];      
      const FlReal wr = GetRad(w);
      const FlReal v1r = GetRad(v1);
      const FlReal v2r = GetRad(v2);
      const FlReal wh = GetHematocrit(w);
      double v1h = GetHematocrit(v1);
      double v2h = GetHematocrit(v2);
      const FlReal wq = GetFlow(w);
      const FlReal v1q = GetFlow(v1);
      const FlReal v2q = GetFlow(v2);

      myAssert( v1==edge_id || v2==edge_id );
      FlReal x0 = (0.4f*0.5f)/wr;
      FlReal A = (-6.96f*0.5f)*std::log( v1r/v2r )/wr;
      FlReal B = 1.0f+(6.98f*0.5f)*(1.0-wh)/wr;
      FlReal FQ1 = v1q/wq;
      FlReal x = (FQ1-x0)/(1.f-2.f*x0);
      if( x<=0.0f || x>=1.0f )
      {
        if( x>=1.0 ) {
          v1h = std::min(1., rbcIn/v1q);
          v2h = 0.0f;
        }
        else {
          v2h = std::min(1., rbcIn/v2q);
          v1h = 0.0f;
        }
      }
      else
      {
        FlReal xx = (1.0f-x)/x;
        FlReal rbc1 = rbcIn/( std::exp(-A)*std::pow(xx,B) + 1.0f );
        FlReal rbc2 = rbcIn-rbc1;
        const FlReal v1q = GetFlow(v1);
        const FlReal v2q = GetFlow(v2);
        v1h = rbc1/(v1q);
        v2h = rbc2/(v2q);
        //myAssert( v1h<=1.0f || v2h<=1.0f );
        if (v1h>1.0f && v2h>1.0f)
        {
          v1h = 1.f;
          v2h = 1.f;
        }
        else if( v1h>1.0f )
        {   // one branch is saturated with rbc
          v1h = 1.0f;
          v2h = (rbcIn-v1q)/v2q;
        }
        else if( v2h>1.0f ) { // or the other
          v2h = 1.0f;
          v1h = (rbcIn-v2q)/v1q;
        }
      }
      SetHematocrit(v1,v1h);
      SetHematocrit(v2,v2h);
      markers[v1] = iteration_id;
      markers[v2] = iteration_id;
    }
    else if(numOut==1 && numIn==1)
    {
      SetHematocrit(vout[0], GetHematocrit(vin[0]));
      markers[vout[0]] = iteration_id;
    }
    else // if num outgoing vessels > 2 or num ingoing vessel > 1
    {
      myAssert((numOut>0 && numIn==0) || numIn>1 || numOut>2);
      for( int i=0; i<numOut; ++i )
      {
        const FlReal h = rbcIn/flowOut;
        SetHematocrit(vout[i],h);
        markers[vout[i]] = iteration_id;
      }
    }

    if( !IsBoundary(upstream_node) && check_hematocrit_range)
    {
      FlReal rbcout = 0.f;
      for( int i=0; i<numOut; ++i )
      {
        rbcout += GetHematocrit(vout[i])*GetFlow(vout[i]);
      }
      const FlReal diff = std::abs(rbcout-rbcIn)/(rbcout+rbcIn+1.e-20);
      if( diff>1.0e-3f ) 
      {
#ifndef SILENCE
        printf("WARNING: UpdateHematocritAtNode: RBC conservation violated %f %% (out: %e, in: %e) at node %i\n",diff,rbcout,rbcIn,upstream_node);
#endif
        /*if( diff>0.5f ) 
        {
          for( int i=0; i<vc->Count(); ++i )
          {
            const Vessel *v = vc->GetEdge(i);
            const VesselNode* wc = v->GetOther( vc );
            printf("  neighbors: v=%i, flags=%s, q=%lf, h=%lf   %s\n",v->Index(), PrintVessFlags(v->flags).c_str(), v->q, v->hematocrit, vc->p<wc->p ? "IN" : "OUT");
          }
        }*/    
      }
    }
  }
  else 
  {
    printf("ERROR: UpdateHematocritAtNode: Too many branches at junction %i\n",upstream_node);
    myAssert(nbcnt<16);
  }
  markers[edge_id] = iteration_id;
}



int CalcFlowWithPhaseSeparation(VesselList3d &vl, const BloodFlowParameters &bloodFlowParameters)
{
  /* I am looking for a memory leak and create additional output with that flag */
//#define trilinos_bug_output
  
#ifndef NDEBUG
  cout << "CalcFlowWithPhaseSeparation called" << endl;
  cout.flush();
#endif

  bool ok = true;
  CompressedFlowNetwork flownet;
  uint returnOfGetFlowNetwork = GetFlowNetwork(flownet, vl, bloodFlowParameters, false);
  if(returnOfGetFlowNetwork>0)
  {
    return 42;
  }

  
//   
  flownet.flow.resize(flownet.num_edges());
  flownet.press.resize(flownet.num_vertices());
  flownet.hema.resize(flownet.num_edges());
  
//   flownet.rad.resize(flownet.num_edges());
  //WARNING
  //flownet.hema.fill(bloodFlowParameters.inletHematocrit); // well defined values for the first compuation of viscosities // Dafug??? this was hardcoded to 0.45 but we would rather actually like to have it filled in GetFlowNetwork i suppose
  //if read in could be zero, which is bad for oxygen calculations
//   cout << "resized " << endl;
//   cout.flush();
  FlArray visc(flownet.num_edges());
  FlArray cond(flownet.num_edges());
  FlArray hema_last(flownet.num_edges(), std::numeric_limits<FlReal>::max());
  FlArray flow_last(flownet.num_edges(), std::numeric_limits<FlReal>::max());
  
  HematocritCalculator hematocritCalculator(flownet, flownet.flow, flownet.press, bloodFlowParameters.inletHematocrit);

  
  Linsys flowsys;

  flowsys.initialize_pattern(flownet.num_vertices(), flownet.edges);
  flowsys.scaling_const = CalcFlowCoeff(bloodFlowParameters.viscosityPlasma, 4., 100.);

  
  //EllipticEquationSolver *solver= new EllipticEquationSolver();
  
  ptree solver_params = make_ptree("output", 1)("preconditioner","multigrid")("use_smoothed_aggregation", false)("max_iter", 200)("throw",false)("conv","rhs")("max_resid",1.e-10);
  
  
  int last_solver_iterations = std::numeric_limits<int>::max();
  /** note:
    * earlier we initialized with nan.
    * This worked for the gxx but not for icc
    * hematocritCalculator.check_hematocrit_range = delta_h<1.e-3 && delta_q<1.e-3 && iteration >= 5;
    * gxx: delta_h<1.e-3 with delta_h as NAN -> false, where icc runtine error
    */
  FlReal delta_h = std::numeric_limits<FlReal>::max(), delta_q = std::numeric_limits<FlReal>::max();
  const int max_iter = 20;
  
  int returnCode = 0; // 1 linsys error
  flowsys.initialize_pattern(flownet.num_vertices(), flownet.edges);
  
  for( int iteration=0; iteration<max_iter && !my::checkAbort(); ++iteration )
  {
    
    //visc is local variable
    CalcViscosities(flownet.rad, flownet.hema, bloodFlowParameters, visc);
    CalcConductivities(flownet.rad, flownet.len, visc, cond);

    
    flowsys.fill_values(flownet.edges, cond, flownet.bcs);
    
#if 1
    /** TF
     * probably not a good idea, put I will go on and NOT consider the 
     * preconditioner at the moment. For the small systems this should work.
     * 
     * I am also not sure, why this additional EllipticEquationSolver stuff is here.
     */
    solver_params.put("keep_preconditioner", iteration>3 && last_solver_iterations<15); 
    //solver->init(flowsys.sys, flowsys.rhs, solver_params);
    //EllipticEquationSolver &&solver = EllipticEquationSolver(flowsys.sys, flowsys.rhs, solver_params);
    //EllipticEquationSolver solver;
    //solver.init(*flowsys.sys, *flowsys.rhs, solver_params);
    //int return_of_solver = solver.solve(*flowsys.lhs);
    flowsys.solve();
#else
    solver_params.put("keep_preconditioner", iteration>3); // && solver.time_iteration<solver.time_precondition);
    int return_of_solver = SolveEllipticEquation(flowsys.sys, flowsys.rhs, flowsys.lhs, solver_params);
#endif
    
//     if( return_of_solver > 0 )
//     {
//       cout << "solver failed to solve!!!!" << endl;
//       returnCode=1;
//       break;
//     }
//     cout << "solver solve called" << endl;
//     cout.flush();
//     last_solver_iterations = solver.iteration_count;
//     cout << "iterations read" << endl;
//     cout.flush();

    flowsys.clear_values();//set rhs and sys to zero, not deleting


    for (int i=0; i<flownet.num_vertices(); ++i)
    {
      flownet.press[i] = flowsys.lhs_get(i);
    }
    SetFlowValues(vl, flownet, cond, flownet.press, flownet.hema);
    
    ok = MarkFailingVesselHidden(vl);
#ifdef trilinos_bug_output
    cout << "MarkFailingVesselHidden" << endl;
    cout.flush();
#endif
    
    if (!ok)
    {
      /* this is where due to numerical inaccuracies (???) small random flows occur which dont respect mass conservation
       * after breaking the computation is restarted with the offending vessels hidden.
       */
      cout << "offending vessels hidden --> restart" << std::endl;
      cout.flush();
      break;
    }

    for (int i=0; i<flownet.num_edges(); ++i)
    {
      flownet.flow[i] = cond[i] * std::abs(flownet.press[flownet.edges[i][0]]-flownet.press[flownet.edges[i][1]]);
    }
#ifdef trilinos_bug_output
    cout << "read flow" << endl;
    cout.flush();
#endif
    //cout << "max relative mass loss: " << CalcFlowResidual(flownet, cond) << endl;

    hematocritCalculator.check_hematocrit_range = delta_h<1.e-3 && delta_q<1.e-3 && iteration >= 5;
    hematocritCalculator.UpdateHematocrit();

#ifndef NDEBUG    
#ifndef SILENT
    cout << "max relative rbc loss: " << CalcHemaResidual(flownet, cond, flownet.hema) << endl;
#endif
#endif

    delta_h = delta_q = 0.;
    FlReal dampen = 0.5;
    for (int i=0; i<flownet.num_edges(); ++i)
    {
      FlReal h = flownet.hema[i], q = flownet.flow[i];
      if (iteration>0)
      {
        FlReal h_last = hema_last[i], q_last = flow_last[i];
        h = h_last*dampen + h*(1.-dampen);
        q = q_last*dampen + q*(1.-dampen);
        #if 0
        delta_h = std::max(delta_h, std::abs(h-h_last)/(h+1.e-20));
        delta_q = std::max(delta_q, std::abs(q-q_last)/(q+1.e-20));
        #endif
        delta_h += my::sqr((h-h_last));//sqare
        delta_q += my::sqr((q-q_last)/(q+1.e-20));
      }
      hema_last[i] = flownet.hema[i] = h;
      flow_last[i] = flownet.flow[i] = q;
    }
    delta_h = sqrt(delta_h)/flownet.num_edges();
    delta_q = sqrt(delta_q)/flownet.num_edges();
    
#ifndef NDEBUG
#ifndef SILENT
    cout << format("dh = %e, dq = %e") % delta_h % delta_q << endl;
#endif
#endif
    if (delta_h < 1.e-6 && delta_q < 1.e-6 && iteration >= 10)
    {
#ifdef trilinos_bug_output
      cout << "Manually call break" << endl;
      cout.flush();
#endif
      returnCode = 2;
      break;
    }
#ifndef NDEBUG
    cout << "next in hematoric calculating loop" << std::endl;cout.flush();
#endif
  }
  if (ok)
  {
    SetFlowValues(vl, flownet, cond, flownet.press, flownet.hema);
#ifndef NDEBUG
    cout << "OK here" << endl;cout.flush();
#endif
  }
  
  if (!ok)
  {
    cout << "before cc" << endl;
    cout.flush();
    ComputeCirculatedComponents(&vl);
    cout << "cc again called" << endl;
    cout.flush();
    returnCode = CalcFlowWithPhaseSeparation(vl, bloodFlowParameters);
  }
  //delete solver;
#ifdef trilinos_bug_output
  cout << "solver deleted" << endl;
  cout.flush();
#endif
  return returnCode;
}

/* NB: this is completely NOT threadsafe!!!!*/
void CalcFlow(VesselList3d &vl, const BloodFlowParameters &params)
{
  calcflow_mutex.lock();
    ComputeCirculatedComponents(&vl);
#ifndef NDEBUG
    std::cout << "ComputeCirculatedComponents called" << std::endl;
    std::cout.flush();
#endif
    if (params.includePhaseSeparationEffect)
    {
      try
      {
        int calcflowReturn = CalcFlowWithPhaseSeparation(vl, params);
        if(calcflowReturn == 1 )
        {
          throw 1;
        }
      }
      catch( int e)
      {
        if( e>0)
        {
          std::cout << "hematocrit calculator Error" << std::endl;
        }
        if( e==1)
        {
          std::cout << "linsys solver did not converge" << std::endl;
        }
        if( e==2)
        {
          std::cout << "Hematocrit error" << std::endl;
        }
      }
    }
    else
      CalcFlowSimple(vl, params, false);
#ifndef NDEBUG
    std::cout << "CalcFlow done" << std::endl;
    std::cout.flush();
#endif
  calcflow_mutex.unlock();
}

#endif //#if 1 // hematocrit calculations


/** 
 * this was in common before
 */
//Clin Hemorheol Microcirc. 2008;39(1-4):243-6.
//Plasma viscosity: a forgotten variable.
// they give the viscosity as 1.1 - 1.3 mPa at 37 deg C as normal value!!!!

//Clin Chem. 1996 Aug;42(8 Pt 1):1189-95.
//Distribution of blood viscosity values and biochemical correlates in healthy adults.
// gives 1.39 +/- 0.08

BloodFlowParameters::BloodFlowParameters()
{
  viscosityPlasma = 4.0e-6; // [kPa s] // i think this is actually the blood viscosity measured in-vitro perhaps?!
  rheology = RheologyForHuman;
  inletHematocrit = 0.45;
  includePhaseSeparationEffect = false;
}

void BloodFlowParameters::assign(const ptree &pt)
{
  viscosityPlasma = pt.get<double>("viscosityPlasma", 4.0e-6); // this is the WRONG!!! value, for backward compatibility!!!!!
  rheology        = pt.get<Rheology>("rheology", RheologyForRats);
  inletHematocrit = pt.get<double>("inletHematocrit", 0.45);
  includePhaseSeparationEffect = pt.get<bool>("includePhaseSeparationEffect", false);
  //includePhaseSeparationEffect = true;
}

ptree BloodFlowParameters::as_ptree() const
{
  return make_ptree("viscosityPlasma", viscosityPlasma)
                   ("rheology", rheology)
                   ("inletHematocrit", inletHematocrit)
                   ("includePhaseSeparationEffect", includePhaseSeparationEffect);
}

ostream& operator<<(ostream &os, Rheology rheology)
{
  switch (rheology)
  {
    case RheologyForHuman: os << "RheologyForHuman"; break;
    case RheologyForRats: os << "RheologyForRats"; break;
    case RheologySecomb2005: os << "RheologySecomb2005"; break;
    default:
      throw std::runtime_error("not implemented");
  }
  return os;
}

istream& operator>>(istream &is, Rheology &rheology)
{
  string s;
  is >> s;
  if (s == "RheologyForHuman")     rheology = RheologyForHuman;
  else if (s == "RheologyForRats") rheology = RheologyForRats;
  else if (s == "RheologySecomb2005") rheology = RheologySecomb2005;
  else throw std::runtime_error("reading enum Rheology error: unknown "+s);
  return is;
}


ostream& operator<<(ostream &os, const BloodFlowParameters &bfparams)
{
  const ptree pt = bfparams.as_ptree();
  boost::property_tree::write_info(os, pt);
  return os;
}

/*------------------------------------------------------------------
--------------------------------------------------------------------*/

double CalcFahraeusEffect(double h, double r, Rheology rheology)
{
  double d = r * 2.;
  if (rheology == RheologyForRats) d /= 0.84;
  const double e1 = 0.415;
  const double e2 = 0.011;
  double t1 = 1 + 1.7 * std::exp(-e1 * d) - 0.6 * std::exp(-e2 * d);
  return h*(h + (1.-h)*t1);
}

namespace ImprovedRelativeViscosityInternal 
{
// Code by Thierry; implements relative viscosity formula from Secomb et al. (2005)
// The main change is reduced viscosity for thinnest possible capillaries (r < 4 mu)
// where the previous formula produced a sharp rise in viscosity.
typedef  double myacc;
inline myacc etha_45(myacc d)
{
  return 220.*std::exp(-1.3*d)+3.2-2.44*std::exp(-0.06*std::pow(d,0.645));
}
inline myacc ca(myacc d)
{
  return (0.8+std::exp(-0.075*d))*(-1+1/(1+std::pow(10,-11)*std::pow(d,12)))+1/(1+std::pow(10,-11)*std::pow(d,12));
}

inline myacc etha_vitro(myacc d, myacc h)
{
  return 1+(etha_45(d)-1)*(std::pow(1-h,ca(d))-1)/(std::pow(1-0.45,ca(d))-1);
}
inline myacc was(myacc d)
{
  if( d<2.4)
  {
    return 0.;
  }
  else
  {
    return (d-2.4)/(d+100.-2.*2.4)*2.6;
  }
}
inline myacc wpeak(myacc d)
{
  if(d<2.4)
  {
    return 0.;
  }
  else if ((2.4<d) and (d<10.5))
  {
    return 1.1*(d-2.4)/(10.5-2.4);
  }
  else
  {
    return 1.1*std::exp(-0.03*(d-10.5));
  }
}
inline myacc weff(myacc d, myacc h)
{
  return was(d)+wpeak(d)*(1+h*1.18);
}
inline myacc Deff(myacc d, myacc h)
{
  return d-2*weff(d,h);
}
inline myacc EtaRel(myacc r, myacc h)
{
  return etha_vitro(2*r,h)*std::pow(2*r/Deff(2*r,h),4);
}
} //ImprovedRelativeViscosityInternal

// "In Vivo" Law from Secomb et al. (1994) Resistance to blood flow in microvessels in vivo
double CalcRelViscosity( double r, double h, Rheology rheology)
{
  h = my::cut(h, 0., 0.999); // hematocrit could be slightly beyond allowed limits [0,1]
  r = my::cut(r, 0.1, 2000.); // protect against division by zero due to r=0; flow resistance of smaller vessels (r<0.1) hence depends on hematocrit and r^4, however eta ist constant in this case.
  if (rheology == RheologyForHuman || rheology == RheologyForRats)
  {
    if (rheology == RheologyForRats) r /= 0.84;  // increases the radius that the viscosity formula "sees" because RBC of rats are smaller
    const double u = 6.0*std::exp( -0.17*r ) + 3.2 - 2.44*std::exp(-0.06*std::pow(2.0*r,0.645) );
    const double rr1 = r/5.0;
    const double rr2 = rr1*rr1;
    const double rr4 = rr2*rr2;
    const double rr8 = rr4*rr4;
    const double rr12 = rr4*rr8;
    const double rc = 1.0/(1.0+10.0*rr12);
    const double C = (0.8+std::exp(-0.15*r))*(-1.0+rc)+rc;
    const double f = (2.0*r)/(2.0*r-1.1);
    const double ff = f*f;
    return (1.0 + (u-1.0)*ff*(std::pow(1.0-h,C)-1.0)/(std::pow(1.0-0.45,C)-1.0))*ff;
  }
  else
  {
    return ImprovedRelativeViscosityInternal::EtaRel(r, h);
  }
}

/*------------------------------------------------------------------
  lookup table for viscosity
--------------------------------------------------------------------*/
#if 0
class ViscosityLUT : private LookupTable2D<FlReal,FlReal>
{
  typedef LookupTable2D<FlReal,FlReal> SUPER;
public:
  ViscosityLUT()
  {
    SUPER::Init( 1.0, 400.0, 200, true,
                 0.01, 0.99, 100, false );

    for( uint y=0; y<ny; ++y )
    {
      const float h = yAxis.GetBucketBegin(y);
      for( uint x=0; x<nx; ++x )
      {
        const float r = xAxis.GetBucketBegin(x);
        SUPER::operator()(x,y) = CalcRelViscosity( r, h, RheologyForRats);
      }
    }
  }
  FlReal DoLookup(FlReal r, FlReal h) const { return SUPER::DoLookup(r,h); }
};


FlReal CalcRelViscosityByTable( FlReal r, FlReal h, Rheology rheology)
{
  myAssert(r > 1.e-3 && h > -1.e-3);
  static ViscosityLUT viscosLut;
  if (rheology == RheologyForRats) r /= 0.84;
  return viscosLut.DoLookup( r, h );
}

void CalcViscosities( const FlArray &rad, const FlArray &hema, FlArray &visc )
{
  // init static table
  static ViscosityLUT viscosLut;

  int ecnt = (int)rad.size();
  if(visc.size()!=ecnt) visc.resize(ecnt);
  
  #pragma omp parallel for
  for(int i=0; i<ecnt; ++i)
  {
    float x = bloodFlowParameters.viscosityPlasma;
    //x *= viscosLut.DoLookup(rad[i], hema[i]);
    
    myAssert(x > 0. && isFinite(x));
    visc[i] = x;
  }
}
#endif

void CalcConductivities(const FlArray &rad, const FlArray &len, const FlArray &visc, FlArray &cond)
{
  int ecnt = (int)rad.size();
  if(cond.size()!=ecnt) cond.resize(ecnt);

  #pragma omp parallel for
  for(int i=0; i<ecnt; ++i)
  {
    if(len[i] == 0)
    {
      cout << " len[" << i << "] = 0" << endl;
    }
    FlReal coeff = CalcFlowCoeff(visc[i],rad[i],len[i]);
    myAssert(coeff > 0.);
    myAssert(isFinite(coeff));
    cond[i] = coeff;
  }
}


void CalcViscosities( const FlArray &rad, const FlArray &hema, const BloodFlowParameters &bloodFlowParameters, FlArray &visc)
{
  int ecnt = (int)rad.size();
  if(visc.size()!=ecnt) visc.resize(ecnt);
  
  #pragma omp parallel for
  for(int i=0; i<ecnt; ++i)
  {
    float x = bloodFlowParameters.viscosityPlasma;
    x *= CalcRelViscosity(rad[i], hema[i], bloodFlowParameters.rheology);
    myAssert(x > 0.);
    myAssert(isFinite(x));
    visc[i] = x;
  }
}


/*------------------------------------------------------------------
--------------------------------------------------------------------*/
/*
 * pries_secomb_1995: design principles of vascular beds; fig5.
 * fit to the plotted data
 * fits quite well and goes indeed from 2.5kPa to 10kPa
 * from 100 micron diameter veins to arteries
 */
#if 1
FlReal PressureRadiusRelation( FlReal rad, bool bArtery )
{
  if( bArtery ) rad = -rad;
  const FlReal A2 = 18.;
  const FlReal A1 = 89.;
  const FlReal x0 = -21.;
  const FlReal dx = 16;
  FlReal p = A2 + (A1-A2)/(1.+std::exp((rad-x0)/dx));
  p *= 0.1f*1.33f;
  return p;
}
#else


/* from Marieb Human Anatomy & Physiology,
 * blood pressure [mmHg]
 * arterioles: 80 - 35 
 * capillaries: 35 - 15
 * venules:     15 - 10
 * veins:       10 - 0
 */

FlReal PressureRadiusRelation( FlReal rad, bool bArtery )
{
  // see testpressureradiusrelation.py
  if( bArtery ) rad = -rad;
  FlReal A2 = 0.;
  FlReal A1 = 89.;
  FlReal x0 = -17.;
  FlReal dx = 10.;
  FlReal p1 = A2 + (A1-A2)/(1.+std::exp((rad-x0)/dx));
  x0 = -21.;
  dx = 2000;
  FlReal p2 = A2 + (A1-A2)/(1.+std::exp((rad-x0)/dx));
  FlReal p = 0.2*p2+0.8*p1;
  p *= 0.1*1.33;
  return p;
}
#endif

/*------------------------------------------------------------------
--------------------------------------------------------------------*/
uint GetFlowNetwork(CompressedFlowNetwork &fl,
                    const VesselList3d &vl,
                    const BloodFlowParameters &bfparams,
                    bool keepTheVesselHematocrit
 		  )
{
  int ecnt = vl.GetECount();
  int ncnt = vl.GetNCount();

  fl.edges.reserve(ecnt);
  fl.org2new_vertex.resize(ncnt);

  // copy only perfused edges
  for(int i=0; i<ecnt; ++i)
  {
    const Vessel* v = vl.GetEdge(i);
    if (!v->IsCirculated()) 
      continue;
    if (!vl.HasLattice())
    {
      if (v->NodeA()->worldpos == v->NodeB()->worldpos) 
	continue; //well obviously this happens :-(
    }
    int a = fl.org2new_vertex.add(v->NodeA()->Index());
    int b = fl.org2new_vertex.add(v->NodeB()->Index());
    myAssert(a!=b);
    fl.edges.push_back(my::make_eqpair(a, b));
  }
  /*sometimes this happens for to small system sizes, 
   *when there are independant arteries and veins
   * ComputeCirculatedComponents fails an all vessels
   * are labeled as uncirculated. For the next hierachical 
   * iteration this fail here!
   */
  if(fl.org2new_vertex.num_added()== 0 or fl.edges.size() == 0)
  {
    //that does not make sence, but could happen during optimization
    return 1;
  }
  myAssert(fl.org2new_vertex.num_added()>0);
  myAssert(fl.edges.size() > 0);

#if 0//this was micheal way-->freak
  // copy bcs and map index numbers
  auto mapKey = [&](const VesselNode* nd) -> int 
  { 
    return fl.org2new_vertex[nd->Index()];
  };
  remap_keys(vl->GetBCMap(), fl.bcs, mapKey);
#endif
  fl.bcs.clear();
  for(int i=0; i<ncnt; ++i)
  {
    const VesselNode* vc= vl.GetNode(i);
    // insert boundary nodes into the bcs array
    if (vc->IsBoundary())
    {
      int id = fl.org2new_vertex[vc->Index()];
      //myAssert(id<=fl.org2new_vertex.num_added());
      /*
       * fancy new stuff
       * use the BCMap if node is in
       */
      if (id != IntegerMap::unused())//and fl.bcs.find(vc->Index()) == fl.bcs.end())
      {
	if( vl.GetBCMap().find(vc) == vl.GetBCMap().end())//not present
	{
	  fl.bcs[id] = FlowBC(FlowBC::PIN, vc->press);
	}
	else //is present
	{
	  fl.bcs[id] = vl.GetBCMap().at(vc);
	}
      }//if not unused
    }//if boundary
  }//for all nodes
#ifdef DEBUG
#ifndef TOTAL_SILENCE
  for(auto bc: fl.bcs)
  {
    printf("first: %i, second: %f\n", bc.first, bc.second.val);
  }
#endif
#endif
  //if (!(flags & FLOWNET_NOPROPERTIES))
  {
    // copy edge properties
    
    fl.len.resize(fl.edges.size());
    fl.rad.resize(fl.edges.size());
    fl.hema.resize(fl.edges.size());
    fl.flow.resize(fl.edges.size());
    
    fl.press.resize(fl.num_vertices());
    
    const VesselList3d::LatticeData &ld = vl.Ld();
    for(int i=0, k=0; i<ecnt; ++i)
    {
      const Vessel* v = vl.GetEdge(i);
      if (!v->IsCirculated()) continue;
      fl.hema[k] = keepTheVesselHematocrit ? v->hematocrit : bfparams.inletHematocrit;
      fl.rad[k] = v->r;
      double wl=0;
      if (!vl.HasLattice())
      {
        wl = (v->NodeA()->worldpos.transpose()-v->NodeB()->worldpos.transpose()).norm();
	myAssert(wl>0);
#ifdef DEBUG
cout<<"vessel no.: "<<v->Index()<<" length: "<<wl<<endl;
cout<<"posA : "<<v->NodeA()->worldpos<<endl;
cout<<"posB : "<<v->NodeB()->worldpos<<endl;
#endif
      }
      else
      {
	wl = v->WorldLength(ld);
        myAssert(wl>0);
      }
      fl.len[k] = wl;
      ++k;
    }
  }
  return 0;
}


void SetFlowValues( VesselList3d &vl,
                    const CompressedFlowNetwork &fl,
                    const FlArray &cond,
                    const FlArray &press,
                    const FlArray &hema
		    )
{
  int ecnt = vl.GetECount();
  int ncnt = vl.GetNCount();
  
  //set node stuff
  for(int i=0; i<ncnt; ++i)
  {
    VesselNode *nd= vl.GetNode(i);
    if (fl.org2new_vertex[nd->Index()] != IntegerMap::unused())
    {
      nd->press = press[fl.org2new_vertex[nd->Index()]];
    }
  }
  //set edge stuff
  for(int i=0, k=0; i<ecnt; ++i)
  {
    Vessel* v = vl.GetEdge(i);
    if(v->IsCirculated())
    {
      int a = fl.org2new_vertex[v->NodeA()->Index()];
      int b = fl.org2new_vertex[v->NodeB()->Index()];
      myAssert(a != IntegerMap::unused() && b != IntegerMap::unused());
      if(hema[k]>0)
        v->hematocrit = hema[k];
#ifndef NDEBUG
      else
        cout<< "bad hema at " << k << "... " << hema[k] <<endl;
#endif
      double wl;
      if( !vl.HasLattice() )
      {
	wl = (v->NodeA()->worldpos.transpose()-v->NodeB()->worldpos.transpose()).norm();
      }
      else
      {
	const VesselList3d::LatticeData &ld = vl.Ld();
	wl = v->WorldLength(ld);
      }
      myAssert(wl>0);
      v->q = cond[k] * std::abs(press[a]-press[b]);
      v->f = CalcShearFromFlowAndCond(v->q, v->r, wl, cond[k]);
      //cout<<v->f*10000.<<endl;
      ++k;
    }
    else
    {
      v->q = 0;
      v->f = 0;
      v->hematocrit = 0.;
    }
  }
}


