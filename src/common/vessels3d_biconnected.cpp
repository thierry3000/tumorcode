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

#include "vessels3d.h"
#include "common.h"
#include "hdfio.h"

#include <boost/range/begin.hpp>

enum {
  BICONCOMP_INIT = -1,
  BICONCOMP_IGNORE = -2,
  BICONCOMP_NONE = -3,
};


#if 0
template<class GRAPH, class NODE_ITER, class IntV>
static void FindBiComponents(
                             GRAPH &g,
                             const NODE_ITER begin,
                             const NODE_ITER end,
                             IntV &compBi, int& numCompBi, IntV &compCon, int& numCompCon
                             )
{
  typedef typename GRAPH::N N;
  typedef typename GRAPH::E E;
  const N NODE_NULL = GRAPH::null_node();
  DynArray<N> vcStack;
  DynArray<N> parent;
  DynArray<word> visited;
  DynArray<int> num;
  DynArray<int> low;
  DynArray<E> edgeStack;

  const size_t vcnt  = g.GetECount();
  const size_t vccnt = g.GetNCount();
  int counter = 1;
  numCompBi = 0;
  numCompCon = 0;

  edgeStack.reserve( vcnt/10 );
  vcStack.reserve( vccnt/10 );
  parent.resize( vccnt, NODE_NULL );
  visited.resize( vccnt, -1 );  // not visited if < 0, otherwise the last visited neighbor
  num.resize( vccnt, 0 );
  low.resize( vccnt, -1 );

  for(NODE_ITER node_iter=begin; node_iter!=end; ++node_iter)
  {
    N vc = node_iter;
    const int i = g.nindex(vc);
    if( visited[i]>=0 ) continue;
    if( compCon[i]==BICONCOMP_IGNORE ) continue;

    vcStack.push_back(vc);

    while( !vcStack.empty() )
    {
      N vc = vcStack.back();
      int vci = g.nindex(vc);
      if( visited[vci]<0 ) // if previously unvisited
      {
        low[g.nindex(vc)] = num[vci] = counter++;
        visited[vci] = 0;
        compCon[vci] = numCompCon;
      }
      else // vc was visited before, now returns from processing some child
      {
        NodeNeighbor<N,E> nb = g.neighbor_of_node(vc,visited[vci]);
        if( low[g.nindex(nb.node)] >= num[vci] )  // biconnected component found
        {
          E vbcc = NULL;
          do
          {
            vbcc = edgeStack.back();
            myAssert(compBi[g.eindex(vbcc)] != BICONCOMP_IGNORE);
            compBi[g.eindex(vbcc)] = numCompBi;
            edgeStack.pop_back();
          }
          while( vbcc!=nb.edge );
          ++numCompBi;
        }
        low[vci] = std::min( low[vci], low[g.nindex(nb.node)] );
        ++visited[vci];
      }
      int cnt = g.neighbor_count(vc);
      while( visited[vci]<cnt ) // loop through unprocessed children
      {
        NodeNeighbor<N,E> nb = g.neighbor_of_node(vc,visited[vci]);
        if (compBi[g.eindex(nb.edge)] != BICONCOMP_IGNORE && compCon[g.nindex(nb.node)] != BICONCOMP_IGNORE)
        {
          int nbni = g.nindex(nb.node);
          if(!(( num[vci]<num[nbni] && visited[nbni]>=0 ) ||
              ( num[vci]>num[nbni] && nb.node==parent[vci] )))
          {
            edgeStack.push_back( nb.edge );
          }
          if( visited[nbni]<0 ) // unvisited neighbor
          {
            parent[nbni] = vc;
            vcStack.push_back( nb.node );  // put neighbor on stack, so it will be visited next
            break; // break from the neighbor loop, and recover later when all children are popped from the stack
          }
          else if( parent[vci]!=nb.node ) // visited neighbor, which is not parent -> backedge
          {
            low[vci] = std::min( low[vci], num[nbni] );
          }
        }
        ++visited[vci];
      }

      if( visited[vci]>=cnt )// if all children visited -> return from "recursive" bccomps call
      {
        vcStack.pop_back();
      }
    } // while stack not empty

    ++numCompCon;
  } // loop over all vessels nodes
#if 0
  for (int i=0; i<vcnt; ++i)
    myAssert(compBi[i] >= 0);
  for (int i=0; i<vccnt; ++i)
    myAssert(compCon[i] >= 0);
#endif
}
#endif


template<class GRAPH, class IntV, class EdgeIndex, class NodeIndex>
static void FindBiComponents(const GRAPH &g,
                             NodeIndex nindex,
                             EdgeIndex eindex,
                             IntV &compBi, int& numCompBi, IntV &compCon, int& numCompCon
                             )
{
  typedef typename GRAPH::node_type N;
  typedef typename GRAPH::edge_type E;
  DynArray<N> vcStack;
  DynArray<int> parent; // references node by index
  DynArray<word> visited;
  DynArray<int> num;
  DynArray<int> low;
  DynArray<E> edgeStack;

  const size_t vcnt  = g.num_edges();
  const size_t vccnt = g.num_nodes();
  int counter = 1;
  numCompBi = 0;
  numCompCon = 0;

  edgeStack.reserve( vcnt/10 );
  vcStack.reserve( vccnt/10 );
  parent.resize( vccnt, -1 );
  visited.resize( vccnt, -1 );  // not visited if < 0, otherwise the last visited neighbor
  num.resize( vccnt, 0 );
  low.resize( vccnt, -1 );

  auto node_range_end = boost::end(g.node_range());
  for(auto node_iter=boost::begin(g.node_range()); node_iter!=node_range_end; ++node_iter)
  {
    N vc = *node_iter;
    const int i = nindex(vc);
    if( visited[i]>=0 ) continue;
    if( compCon[i]==BICONCOMP_IGNORE ) continue;

    vcStack.push_back(vc);

    while( !vcStack.empty() )
    {
      N vc = vcStack.back();
      int vci = nindex(vc);
      if( visited[vci]<0 ) // if previously unvisited
      {
        low[nindex(vc)] = num[vci] = counter++;
        visited[vci] = 0;
        compCon[vci] = numCompCon;
      }
      else // vc was visited before, now returns from processing some child
      {
        auto nb = g.adjacent(vc,visited[vci]);
        if( low[nindex(nb.node)] >= num[vci] )  // biconnected component found
        {
          E vbcc;
          do
          {
            vbcc = edgeStack.back();
            myAssert(compBi[eindex(vbcc)] != BICONCOMP_IGNORE);
            compBi[eindex(vbcc)] = numCompBi;
            edgeStack.pop_back();
          }
          while( vbcc!=nb.edge );
          ++numCompBi;
        }
        low[vci] = std::min( low[vci], low[nindex(nb.node)] );
        ++visited[vci];
      }
      int cnt = g.degree(vc);
      while( visited[vci]<cnt ) // loop through unprocessed children
      {
        auto nb = g.adjacent(vc,visited[vci]);
        if (compBi[eindex(nb.edge)] != BICONCOMP_IGNORE && compCon[nindex(nb.node)] != BICONCOMP_IGNORE)
        {
          int nbni = nindex(nb.node);
          if(!(( num[vci]<num[nbni] && visited[nbni]>=0 ) ||
              ( num[vci]>num[nbni] && nindex(nb.node)==parent[vci] )))
          {
            edgeStack.push_back( nb.edge );
          }
          if( visited[nbni]<0 ) // unvisited neighbor
          {
            parent[nbni] = nindex(vc);
            vcStack.push_back( nb.node );  // put neighbor on stack, so it will be visited next
            break; // break from the neighbor loop, and recover later when all children are popped from the stack
          }
          else if( parent[vci]!=nindex(nb.node)) // visited neighbor, which is not parent -> backedge
          {
            low[vci] = std::min( low[vci], num[nbni] );
          }
        }
        ++visited[vci];
      }

      if( visited[vci]>=cnt )// if all children visited -> return from "recursive" bccomps call
      {
        vcStack.pop_back();
      }
    } // while stack not empty

    ++numCompCon;
  } // loop over all vessels nodes
#if 0
  for (int i=0; i<vcnt; ++i)
    myAssert(compBi[i] >= 0);
  for (int i=0; i<vccnt; ++i)
    myAssert(compCon[i] >= 0);
#endif
}



void ComputeCirculatedComponents(VesselList3d* list)
{
  //printf("We are in ComputeCirculatedComponents!");
  vector<Vessel*> vboundary;
  vboundary.reserve( 32 );
  vector<VesselNode*> vcboundary;
  vcboundary.reserve( 32 );
  for( uint i=0; i<list->GetNCount(); ++i )
  {
    VesselNode *vc = list->GetNode(i);
    vc->flags.DelBits(CIRCULATED);
    if( !vc->flags.GetBits(BOUNDARY) ) continue;
    vcboundary.push_back( vc );
  }
  myAssert( vcboundary.size()>1 );

  uint nboundary = (uint)vcboundary.size();
  if (nboundary <= 2)  --nboundary; // (if only 2 boundary points one cannot add a cycle, because it would introduce a doublicate edge)
  for( uint i=0; i<nboundary; ++i )
  {
    VesselNode* vc1 = vcboundary[i];
    VesselNode* vc2 = vcboundary[(i+1)%vcboundary.size()];
    Vessel* v = list->InsertTempEdge( vc1, vc2 );
    vboundary.push_back( v  );
  }

  vector<int> compBi,compCon;
  compBi.resize(list->GetECount(), BICONCOMP_INIT);
  compCon.resize(list->GetNCount(), BICONCOMP_INIT);
  for (int i=0; i<list->GetECount(); ++i)
  {
    const Vessel *v = list->GetEdge(i);
    const VesselNode *a = v->NodeA(),
                     *b = v->NodeB();
    if (v->flags&HIDDEN)
      compBi[v->Index()] = BICONCOMP_IGNORE;
    if (a->flags&HIDDEN)
      compCon[a->Index()] = BICONCOMP_IGNORE;
    if (b->flags&HIDDEN)
      compCon[b->Index()] = BICONCOMP_IGNORE;
  }

  auto nindex = [=](const VesselNode* nd) { return nd->Index(); };
  auto eindex = [=](const Vessel* v) { return v->Index(); };
  
  int cntBi,cntCon;
  FindBiComponents(
                    list->Graph(),
                    nindex,
                    eindex,
                    compBi, cntBi, compCon, cntCon);

  const int maincompBi  = compBi[vboundary[0]->Index()];
  const int maincompCon = compCon[vboundary[0]->NodeA()->Index()];
  
  // cleanup and check
  for( int i=0; i<vboundary.size(); ++i )
  {
    myAssert( compBi[vboundary[i]->Index()]==maincompBi );
    myAssert( compCon[vboundary[i]->NodeA()->Index()]==maincompCon );
    myAssert( compCon[vboundary[i]->NodeB()->Index()]==maincompCon );
    list->DeleteTempEdge( vboundary[i] );
  }

  // update flags
  for( int i=0; i<list->GetECount(); ++i )
  {
    Vessel* v = list->GetEdge(i);
    //myAssert( !(v->flags&BOUNDARY) );
    if( compBi[v->Index()]==maincompBi ) {
      v->flags.AddBits(CIRCULATED);
      v->NodeA()->flags.AddBits(CIRCULATED);
      v->NodeB()->flags.AddBits(CIRCULATED);
    }
    else {
      v->flags.DelBits(CIRCULATED);
    }

    myAssert((v->flags&HIDDEN || v->NodeA()->flags&HIDDEN || v->NodeB()->flags&HIDDEN) || (compCon[v->NodeB()->Index()]==maincompCon) == (compCon[v->NodeA()->Index()]==maincompCon) );
    if( compCon[v->NodeA()->Index()]==maincompCon ) {
      v->flags.AddBits(CONNECTED);
      v->NodeA()->flags.AddBits(CONNECTED);
      v->NodeB()->flags.AddBits(CONNECTED);
    }
    else {
      v->flags.DelBits(CONNECTED);
      v->NodeA()->flags.DelBits(CONNECTED);
      v->NodeB()->flags.DelBits(CONNECTED);
    }
  }

#if 0
  void WriteVesselList3d(const VesselList3d& vl, h5cpp::Group vesselgroup, const ptree& params = ptree());
  {// maximum debug
    static int counter = 0;
    compBi.resize(list->GetECount());
    h5cpp::File f("biconnecteddebug.h5", counter==0 ? "w" : "r+");
    h5cpp::Group g = f.root().create_group(str(format("no%i") % (counter++))).create_group("vessels");
    WriteVesselList3d(*list, g, ptree());
    h5cpp::create_dataset_range(g,"edges/compbi", compBi);
    h5cpp::create_dataset_range(g,"nodes/connecti", compCon);
  }
#endif
}


#include "mwlib/compressed_row_undirected_graph.h"
#include <boost/unordered_set.hpp>


typedef DynArray<my::eqpair<int> > EdgeArray;
typedef DynArray<uint> FlagArray;



struct UndirectedEdge : public my::eqpair<int>
{
  UndirectedEdge(const my::eqpair<int> &e) : my::eqpair<int>(e) {}
  UndirectedEdge(int a, int b) : my::eqpair<int>(a, b) {}
};

inline bool operator==(const UndirectedEdge &e1, const UndirectedEdge &e2)
{
    return (e1.first==e2.first && e1.second==e2.second) ||
           (e1.first==e2.second && e1.second==e2.first);
}

std::size_t hash_value(UndirectedEdge const& e)
{
    std::size_t seed = 0;
    bool t = e.first<e.second;
    boost::hash_combine(seed, t ? e.first : e.second);
    boost::hash_combine(seed, t ? e.second : e.first);
    return seed;
}

#if 0
void ComputeCirculatedComponents(const DynArray<my::eqpair<int> > &org_edges, const int num_nodes, DynArray<uint> &eflags, DynArray<uint> &nflags)
{
  const int num_edges = org_edges.size();

  boost::unordered_set<UndirectedEdge> boundary_edges;
  boost::unordered_set<int> boundary_nodes;

  CompressedGraph graph;
  graph.beginConstruct();
  
  for (int i=0; i<num_edges; ++i)
  {
    int a, b; tie(a,b) = org_edges[i];
    const bool is_bd[2] = {
      bool(nflags[a] & BOUNDARY),
      bool(nflags[b] & BOUNDARY)
    };
    if (is_bd[0]) boundary_nodes.insert(a);
    if (is_bd[1]) boundary_nodes.insert(b);
    if (is_bd[0] && is_bd[1])
      boundary_edges.insert(UndirectedEdge(org_edges[i]));
    graph.add_edge(a, b);

    eflags[i] &= ~(CIRCULATED|CONNECTED);
  }
  for (int i=0; i<num_nodes; ++i)
    nflags[i] &= ~(CIRCULATED|CONNECTED);
  

  DynArray<int> compBi,compCon;
  compBi.resize(graph.num_edges(), BICONCOMP_INIT);
  compCon.resize(graph.num_nodes(), BICONCOMP_INIT);
  int maincompBi = BICONCOMP_NONE, maincompCon = BICONCOMP_NONE;
  
  if (boundary_nodes.size() > 1)
  {
    for (auto it = boundary_nodes.begin(); it != boundary_nodes.end(); ++it)
    {
      auto itnext = it; ++itnext;
      if (itnext == boundary_nodes.end()) itnext = boundary_nodes.begin();

      UndirectedEdge new_edge(*it, *itnext);
      if (boundary_edges.find(new_edge) == boundary_edges.end())
      {
        boundary_edges.insert(new_edge); // no duplicate edge allowed
        graph.add_edge(new_edge[0], new_edge[1]);
      }
    }

    graph.endConstruct();
    /*
    * graph is constructed. It should have the same topology as the original.
    * In addition it has edges at indices >= num_edges which connect boundary
    * nodes in order to form a kind of ring. The nodes should be 100% identical
    * to the ones given by the original edges.
    */

    int cntBi,cntCon;
    auto nindex = [=](int i) { return i; }; // just identity (we need this only for VesselList3d
    auto eindex = [=](int i) { return i; };
    FindBiComponents<CompressedGraph>(graph,
                    nindex,
                    eindex,
                    compBi, cntBi, compCon, cntCon);

    maincompBi  = compBi[num_edges]; // first inserted edge
    maincompCon = compCon[graph.edge_nodes(num_edges)[0]];
  }


  for( int i=0; i<num_edges; ++i ) // for the original edges
  {
    int a, b; tie(a,b) = org_edges[i];
    if( compBi[i]==maincompBi )
    {
      eflags[i] |= CIRCULATED;
      nflags[a] |= CIRCULATED;
      nflags[b] |= CIRCULATED;
    }

    //myAssert((v->flags&HIDDEN || v->NodeA()->flags&HIDDEN || v->NodeB()->flags&HIDDEN) || (compCon[v->NodeB()->Index()]==maincompCon) == (compCon[v->NodeA()->Index()]==maincompCon) );
    if (compCon[a]==maincompCon )
    {
      nflags[a] |= CONNECTED;
      nflags[b] |= CONNECTED;
      eflags[i] |= CONNECTED;
    }
  }
}
#endif
