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
#ifndef CALCFLOW_OPTNET_H
#define CALCFLOW_OPTNET_H

#include "calcflow_common.h"
#include "mwlib/compressed_row_undirected_graph.h"

struct CNEdgeSeg
{
  my::eqpair<int> nodes;
  FlReal w;  // effective conductivity of all contained vessels
};

struct CNBackLink
{
  float w_vess;  // conductivity of the vessel segment
  float w_accum; // from seg.a to seg.b, the successively computed effective conductivity
  int node_id;
};

class NetworkLinearChainRemover
{
public:
  const FeArray& edges() const { return new_edges; }
  const FlArray& cond() const { return new_cond; }

  void initialize(const CompressedGraph &old_graph, const FlArray &old_cond, const DynArray<bool> &boundary_points_)
  {
    boundary_points = &boundary_points_;
    _old_edges = &old_graph.edges();
    _old_nbs   = &old_graph;
    _old_cond  = &old_cond;
    _old_press = NULL;
    ecnt = (int)_old_edges->size();
    ncnt = old_graph.num_nodes();
    scnt = 0;
    new_edges.clear();
    new_cond.clear();

    new_edges.reserve(GetECount()*0.75f);
    new_cond.reserve(GetECount()*0.75f);

    segs_contained.reserve(GetECount());
    segs_contained_begin.reserve(GetECount()*0.75);
    DynArray<uchar> markers(GetECount());
    for( int iv=0; iv<GetECount(); ++iv )
    {
      if( markers[iv] ) continue;
      const my::eqpair<int> e = GetEdge(iv);

      for( int k=0; k<2; ++k )
      { // look at both endnodes, to start a new segment
        const int ind = e[k];
        if( FindPipeContinuation(ind,iv)<0 )
        {
          CNEdgeSeg seg;
          seg.nodes=my::make_eqpair(ind,-1);
          seg.w = 0;

          segs_contained_begin.push_back( segs_contained.size() );
          RecurseBuildEdge( ind, iv, seg, markers );
          k = 2; // done
          AddNewEdge(seg.nodes,seg.w);
          ++scnt;
        }
      }
    }
    segs_contained_begin.push_back( segs_contained.size() );

    _old_edges = NULL;
    _old_cond = NULL;
    _old_nbs = NULL;
  }

  void distribute_pressures( FlArray& press )
  {
    // assume pressures at actually used nodes have been set correctly
    // put the proper pressures at the nodes of linear chains
    _old_press = &press; // setup internal vars

    for( int is=0; is<scnt; ++is )
    {
      my::eqpair<int> e; FlReal w;
      GetNewEdge(is,e,w);
      const CNBackLink* const begin = GetBackLinkBegin(is);
      const CNBackLink* const end   = GetBackLinkEnd(is);
      float press_a = GetPress(e[0]);
      float press_b = GetPress(e[1]);
      myAssert( press_a>=0 && press_b>=0 );
      for( const CNBackLink* lnk=end; lnk!=begin; --lnk )
      {
        const CNBackLink* prev=lnk-1;
        float p = (prev->w_accum*press_a + lnk->w_vess*press_b)/(prev->w_accum+lnk->w_vess);
        SetPress(lnk->node_id,p);
        press_b = p;
      }
    }
  }

private:
  DynArray<CNBackLink> segs_contained;
  DynArray<int> segs_contained_begin;
  int ecnt,ncnt,scnt;

  const DynArray<bool> *boundary_points;
  const FlArray *_old_cond;
  const FeArray *_old_edges;
  const CompressedGraph *_old_nbs;
  FlArray *_old_press;

  FeArray new_edges;
  FlArray new_cond;

  int GetSCount() const { return scnt; }
  int GetECount() const { return ecnt; }
  int GetNCount() const { return ncnt; }
  FlReal GetConductivity(int i) const { return (*_old_cond)[i]; }

  my::eqpair<int> GetEdge(int i) const { return (*_old_edges)[i]; }
  int GetEdgesOfNodeCount(int i) const { return (*_old_nbs).degree(i); }
  NodeNeighbor<int,int> GetEdgeOfNode(int i, int k) const { return (*_old_nbs).adjacent(i,k); }

  FlReal  GetPress(int i) const         { return (*_old_press)[i]; }
  void    SetPress(int i, FlReal p)     { (*_old_press)[i]=p; }

  void AddNewEdge(my::eqpair<int> e, FlReal cond)
  {
    new_edges.push_back(e);
    new_cond.push_back(cond);
  }
  void GetNewEdge(int i, my::eqpair<int> &e, FlReal &cond) const
  {
    e=new_edges[i];
    cond=new_cond[i];
  }
  const CNBackLink* GetBackLinkBegin( int seg ) const { return get_ptr(segs_contained)+segs_contained_begin[seg]; }
  const CNBackLink* GetBackLinkEnd( int seg ) const { return get_ptr(segs_contained)+segs_contained_begin[seg+1]-1; }

  int FindPipeContinuation( const int ind, const int iv)
  {
    if ((*boundary_points)[ind]) return -1;
    uint vc_num_circ = 0;
    int res_nbe = -1;
    int cnt = GetEdgesOfNodeCount(ind);
    for( int i=0; i<cnt; ++i )
    {
      int e = GetEdgeOfNode(ind,i).edge;
      if( e==iv ) continue;
      ++vc_num_circ;
      res_nbe = e;
    }
    return vc_num_circ==1 ? res_nbe : -1;
  }

  void RecurseBuildEdge( int ind, int iv, CNEdgeSeg &seg, std::vector<uchar> &markers )
  {
    // check adjacent vessels
    bool bfirst = true;
    while( iv>=0 )
    {
      // init backlink entry
      CNBackLink lnk;
      lnk.w_vess = GetConductivity(iv);
      lnk.node_id = ind;
      seg.w = bfirst ? lnk.w_vess : (seg.w*lnk.w_vess)/(seg.w+lnk.w_vess);

      lnk.w_accum = seg.w;
      segs_contained.push_back(lnk);

      myAssert(markers[iv]==0);
      markers[iv] = 1;

      const my::eqpair<int> e = GetEdge(iv);
      bfirst = false;
      ind = e[0]==ind ? e[1] : e[0];
      iv = FindPipeContinuation(ind,iv);
    }
    seg.nodes[1] = ind;
  }
};

#endif