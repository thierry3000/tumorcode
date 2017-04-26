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

#include "calcflow_linsys.h"
#include "calcflow_optnet.h"
#include "calcflow_common.h"
#include "vessels3d.h"
#include "shared-objects.h"

#include "calcflow.h"
#include <boost/foreach.hpp>
#include <limits>
#include <fstream>
#include "Ifpack_Utils.h"



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
    
    if (nd->IsBoundary()) continue;
    double qsum = 0.;
    double qabs = 0.;
    int n = 0;
    for (int i=0; i<nd->Count(); ++i)
    {
      auto nb = nd->GetNode(i);
      if (!nb.edge->IsCirculated()) continue;
      if (nb.node->press < nd->press)
      {
        qsum -= nb.edge->q;
      }
      else if (nb.node->press > nd->press)
      {
        qsum += nb.edge->q;
      }
      else
        myAssert(nb.edge->q == 0.);
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
#ifndef SILENT
  cout << "max relative fluid loss: " << max_res << endl;
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
  GetFlowNetwork(flownet, &vl, bloodFlowParameters, keepTheVesselHematocrit);

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
  
  {
    Linsys flowsys;
    flowsys.scaling_const = CalcFlowCoeff(bloodFlowParameters.viscosityPlasma, 4., 100.); //Median(cond);
    flowsys.initialize_pattern(flownet.num_vertices(), flownet.edges);
    flowsys.fill_values(flownet.edges, cond, flownet.bcs);

#if 0   // thierrys sparse output
    time_t rawtime;
    struct tm*timeinfo;
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    std::ofstream myMatrixFile;
    std::string myString ("Epetra_Matrix");
    myString = myString + asctime(timeinfo) +".txt";
    myMatrixFile.open(myString.c_str());
    flowsys.sys->Print(myMatrixFile);
    std::string myString2 ("Epetra_Sparsity_Pattern");
    myString2 = myString2 + asctime(timeinfo) + ".ps";
    Ifpack_PrintSparsity(*flowsys.sys,myString2.c_str());
    myMatrixFile.close();
#endif

    flowsys.solve();
  
    flownet.press.resize(flownet.num_vertices());
    //flownet.len.resize(flownet.num_edges());
    //flownet.rad.resize(flownet.num_edges());
    //flownet.hema.resize(flownet.num_edges());
    //flownet.flow.resize(flownet.num_edges());
    
    for (int i=0; i<flownet.num_vertices(); ++i) flownet.press[i] = flowsys.lhs_get(i);
    SetFlowValues(&vl, flownet, cond, flownet.press, flownet.hema);
    
    
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
  }

#ifdef DEBUG
  for (auto iter = flownet.bcs.begin(); iter != flownet.bcs.end(); ++iter)
  {
    cout << format("boundary node id: %i, value = %f") % iter->first % iter->second.val << endl;
  }
#endif
  //cout << "max relative mass loss: " << CalcFlowResidual(flownet, cond) << endl;

  //h5cpp::File f("debugvessels.h5","w");
  //WriteVesselList3d(*vl, f.root().create_group("vessels"), make_ptree("w_all",true));

  }
  bool ok = MarkFailingVesselHidden(vl);
  if (!ok)
  {
    ComputeCirculatedComponents(&vl);
    CalcFlowSimple(vl, bloodFlowParameters, keepTheVesselHematocrit);
  }
}



#if 1
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
    printf("warning hematocrit %f out of bounds at vessel %i\n", x, i); 
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
        printf("WARNING: UpdateHematocritAtNode: RBC conservation violated %f %% (out: %e, in: %e) at node %i\n",diff,rbcout,rbcIn,upstream_node);
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



void CalcFlowWithPhaseSeparation(VesselList3d &vl, const BloodFlowParameters &bloodFlowParameters)
{
#ifndef SILENT
  cout << "calcflow (with hematocrit)" << endl;
#endif
  bool ok = true;
  {
    CompressedFlowNetwork flownet;
    GetFlowNetwork(flownet, &vl, bloodFlowParameters, false);
    flownet.flow.resize(flownet.num_edges());
    flownet.press.resize(flownet.num_vertices());
    //flownet.hema.fill(params.inletHematocrit); // well defined values for the first compuation of viscosities // Dafug??? this was hardcoded to 0.45 but we would rather actually like to have it filled in GetFlowNetwork i suppose
    FlArray visc(flownet.num_edges());
    FlArray cond(flownet.num_edges());
    FlArray hema_last(flownet.num_edges(), getNAN<FlReal>()), flow_last(flownet.num_edges(), getNAN<FlReal>());

    HematocritCalculator hematocritCalculator(flownet, flownet.flow, flownet.press, bloodFlowParameters.inletHematocrit);
    Linsys flowsys;
    flowsys.initialize_pattern(flownet.num_vertices(), flownet.edges);
    flowsys.scaling_const = CalcFlowCoeff(bloodFlowParameters.viscosityPlasma, 4., 100.);
    EllipticEquationSolver solver;
    ptree solver_params = make_ptree("output", 1)("preconditioner","multigrid")("use_smoothed_aggregation", false)("max_iter", 200)("throw",false)("conv","rhs")("max_resid",1.e-10);

    int last_solver_iterations = std::numeric_limits<int>::max();
    /* note:
     * earlier we initialized with nan.
     * This worked for the gxx but not for icc
     * hematocritCalculator.check_hematocrit_range = delta_h<1.e-3 && delta_q<1.e-3 && iteration >= 5;
     * gxx: delta_h<1.e-3 with delta_h as NAN -> false, where icc runtine error
     */
    FlReal delta_h = std::numeric_limits<FlReal>::max(), delta_q = std::numeric_limits<FlReal>::max();

    const int max_iter = 20;
    for( int iteration=0; iteration<max_iter && !my::checkAbort(); ++iteration )
    {
      CalcViscosities(flownet.rad, flownet.hema, bloodFlowParameters, visc);
      CalcConductivities(flownet.rad, flownet.len, visc, cond);

      flowsys.fill_values(flownet.edges, cond, flownet.bcs);

      solver_params.put("keep_preconditioner", iteration>3 && last_solver_iterations<15); // && solver.time_iteration<solver.time_precondition);
      solver.init(*flowsys.sys, *flowsys.rhs, solver_params);
      solver.solve(*flowsys.lhs);
      last_solver_iterations = solver.iteration_count;

      flowsys.clear_values();

      for (int i=0; i<flownet.num_vertices(); ++i) flownet.press[i] = flowsys.lhs_get(i);

      SetFlowValues(&vl, flownet, cond, flownet.press, flownet.hema);
      ok = MarkFailingVesselHidden(vl);
      if (!ok) break; // this is where due to numerical inaccuracies (???) small random flows occur which dont respect mass conservation
      // after breaking the computation is restarted with the offending vessels hidden.

      for (int i=0; i<flownet.num_edges(); ++i)
      {
        flownet.flow[i] = cond[i] * std::abs(flownet.press[flownet.edges[i][0]]-flownet.press[flownet.edges[i][1]]);
      }

      //cout << "max relative mass loss: " << CalcFlowResidual(flownet, cond) << endl;

      hematocritCalculator.check_hematocrit_range = delta_h<1.e-3 && delta_q<1.e-3 && iteration >= 5;
      hematocritCalculator.UpdateHematocrit();

#ifndef SILENT
      cout << "max relative rbc loss: " << CalcHemaResidual(flownet, cond, flownet.hema) << endl;
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
          delta_h += my::sqr((h-h_last));
          delta_q += my::sqr((q-q_last)/(q+1.e-20));
        }
        hema_last[i] = flownet.hema[i] = h;
        flow_last[i] = flownet.flow[i] = q;
      }
      delta_h = sqrt(delta_h)/flownet.num_edges();
      delta_q = sqrt(delta_q)/flownet.num_edges();
#ifndef SILENT
      cout << format("dh = %e, dq = %e") % delta_h % delta_q << endl;
#endif
      if (delta_h < 1.e-6 && delta_q < 1.e-6 && iteration >= 10) break;
    }
    if (ok)
      SetFlowValues(&vl, flownet, cond, flownet.press, flownet.hema);
  }
  if (!ok)
  {
    ComputeCirculatedComponents(&vl);
    CalcFlowWithPhaseSeparation(vl, bloodFlowParameters);
  }
}

/* NB: this is completely NOT threadsafe!!!!*/
void CalcFlow(VesselList3d &vl, const BloodFlowParameters &params)
{
#pragma omp single
  {
  ComputeCirculatedComponents(&vl);
  if (params.includePhaseSeparationEffect)
    CalcFlowWithPhaseSeparation(vl, params);
  else
    CalcFlowSimple(vl, params, true);
  }//#pragma omp single
}
#endif
