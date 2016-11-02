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
#ifndef LISTGRAPH_H
#define LISTGRAPH_H

#include "dynamicarray.h"
#include <algorithm>
#include "helpers-containers.h"
#include "helpers-vec.h"

#include <boost/range/iterator_range.hpp>

// forward declare from lin-system-helpers.h
void CalcMatrixColumnMemoryOffsets( int n, int *ja );

/*----------
  a vector holding points to objects that have a public index member
  deletion from the vector is done by swapping with the last element
  the index is such that cont[i]->index == i, 
  This is usefull when you don't know where in the array the object is
  and you want to delete.
  ----------*/
template<class T>
struct comp_addr
{
  bool operator()(const T* a, const T* b) { return a<b; }
};

template<class T> 
class IndexObjVector : public DynArray<T*>
{  
public:
  typedef DynArray<T*> SUPER;
  inline T* &indexed_push_back( T* p ) { p->index=(int)SUPER::size(); SUPER::push_back(p); return SUPER::back(); }
  inline void indexed_remove( T* p ) 
  { 
    myAssert( SUPER::operator[](p->index)==p );
    const int index=p->index;
    if(SUPER::size()>0) 
    {      
      p->index = SUPER::back()->index;
      SUPER::back()->index = index;
      SUPER::operator[](index) = SUPER::back();
    }
    SUPER::pop_back();
  }

  void sort_by_addr()
  {
    std::sort(SUPER::begin(),SUPER::end(),comp_addr<T>());
    for( uint i=0; i<SUPER::size(); ++i )
    {
      SUPER::operator[](i)->index = i;
    }
  }
  
  std::size_t estimateMemoryUsage() const
  {
    return sizeof(*this) + (SUPER::size()>0 ? SUPER::size() * estimateMemoryUsage((*this)[0]) : 0);
  }
};


template< class NodeT, class EdgeT, size_t NBCNT >
class ListTreeNode 
{
  EdgeT* nb[NBCNT];
  uint cnt;
  int index;
  friend class IndexObjVector<NodeT>;
  
public:
  
  ListTreeNode() : cnt(0),index(-1)  { ClearMem(nb,NBCNT); }

  void AddEdge( EdgeT* v ) { myAssert( cnt<NBCNT );  StatArray::push_back<EdgeT*,uint>( nb, cnt, v ); }
  void RemoveEdge( EdgeT* v ) { StatArray::remove<EdgeT*,uint>( nb, cnt, StatArray::find<EdgeT*,uint>( nb, cnt, v ) ); }

  int FindEdgeIndex( const EdgeT* v ) const { return StatArray::find<EdgeT*,uint>( nb, cnt, const_cast<EdgeT*>(v) ); }
  uint Count() const { return cnt; }
  int Index() const { return index; }
  
  EdgeT* GetEdge( int i ) { myAssert(i>=0 && i<cnt ); return nb[i]; }
  const EdgeT* GetEdge( int i ) const { myAssert(i>=0 && i<cnt ); return nb[i]; }
  NodeNeighbor<NodeT*,EdgeT*> GetNode( int i )
  {
    myAssert( i>=0 && i<cnt && nb[i] && (nb[i]->NodeA()==this || nb[i]->NodeB()==this) );
    return NodeNeighbor<NodeT*,EdgeT*>( (nb[i]->NodeA()==this) ? nb[i]->NodeB() : nb[i]->NodeA(), nb[i] );
  }
  NodeNeighbor<const NodeT*,const EdgeT*> GetNode( int i ) const
  {
    myAssert( i>=0 && i<cnt && nb[i] && (nb[i]->NodeA()==this || nb[i]->NodeB()==this) );
    return NodeNeighbor<const NodeT*,const EdgeT*>( (nb[i]->NodeA()==this) ? nb[i]->NodeB() : nb[i]->NodeA(), nb[i] );
  }
};




template< class NodeT, class EdgeT >
class ListTreeEdge
{
  NodeT *a,*b;
  int index;
  friend class IndexObjVector< EdgeT >;

public:

  ListTreeEdge() : a(NULL),b(NULL) {}  

  NodeT* NodeA() { return a; }
  NodeT* NodeB() { return b; }
  const NodeT* NodeA() const { return a; }
  const NodeT* NodeB() const { return b; }
  NodeT* GetNode( int indx ) { return indx==0 ? a : b; }
  const NodeT* GetNode( int indx ) const { return indx==0 ? a : b; }
  int Index() const { return index; }

  void AttachNodes( NodeT* vc1, NodeT* vc2 )
  {
    myAssert( this->a==NULL || this->a->FindEdgeIndex(static_cast<const EdgeT*>(this))<0 );
    myAssert( this->b==NULL || this->b->FindEdgeIndex(static_cast<const EdgeT*>(this))<0 );
    myAssert( vc1!=NULL && vc2!=NULL );
    this->a = vc1;
    this->b = vc2;
    vc1->AddEdge( static_cast<EdgeT*>(this) );
    vc2->AddEdge( static_cast<EdgeT*>(this) );
  }

  bool HasNode( const NodeT* vc ) const { return a==vc || b==vc; }
  
  NodeT* GetOther( const NodeT* vc )
  { 
    if(a==vc) return b; 
    else if(b==vc) return a; 
    myAssert(false);
    return NULL;
  }

  const NodeT* GetOther( const NodeT* vc ) const 
  { 
    return const_cast<ListTreeEdge*>(this)->GetOther( vc );
  }
};


#define ASSERT_AT_DELETE_UNKNOWN_ELEMENT

/* A Graph class.
 * Random Note: Const Correctness is not properly implemented. You could obtain an edge or node pointer from a const ListGraph and modify its neighbor information.
 */
template< class NodeT, class EdgeT>
class ListGraph : boost::noncopyable
{
public:
  typedef NodeT* node_type;
  typedef EdgeT* edge_type;
  typedef NodeNeighbor<NodeT*,EdgeT*> NeighborType;

  ListGraph() {}
  ~ListGraph() { Flush(); }

public:
  inline uint num_nodes() const { return (uint)nodelist.size(); }
  inline uint num_edges() const { return (uint)edgelist.size(); }

  inline const my::eqpair<node_type> edge_nodes(edge_type e) const { return my::eqpair<node_type>(e->NodeA(),e->NodeB()); }
  inline const NeighborType adjacent(node_type n, int i) const { return n->GetNode(i); }
  inline int degree(node_type n) const { return n->Count(); }

  inline EdgeT*             edge( size_t i ) const { return edgelist[i]; }
  //inline const EdgeT*       edge( size_t i ) const { return edgelist[i]; }
  inline NodeT*             node( size_t i ) const { return nodelist[i]; }
  //inline const NodeT*       node( size_t i ) const { return nodelist[i]; }

  typedef typename IndexObjVector<NodeT>::iterator node_iterator;
  typedef typename IndexObjVector<NodeT>::const_iterator const_node_iterator;
  typedef typename IndexObjVector<EdgeT>::iterator edge_iterator;
  typedef typename IndexObjVector<EdgeT>::const_iterator const_edge_iterator;

  boost::iterator_range<node_iterator> node_range() { return boost::make_iterator_range(nodelist.begin(), nodelist.end()); }
  boost::iterator_range<const_node_iterator> node_range() const { return boost::make_iterator_range(nodelist.begin(), nodelist.end()); }

  boost::iterator_range<edge_iterator> edge_range() { return boost::make_iterator_range(edgelist.begin(), edgelist.end()); }
  boost::iterator_range<const_edge_iterator> edge_range() const { return boost::make_iterator_range(edgelist.begin(), edgelist.end()); }

  bool isInGraph( const EdgeT* v ) const { return v->Index()<edgelist.size() && edgelist[v->Index()]==v; }
  bool isInGraph( const NodeT* n ) const { return n->Index()<nodelist.size() && nodelist[n->Index()]==n; }

  std::size_t estimateMemoryUsage() const { return nodelist.estimateMemoryUsage()+edgelist.estimateMemoryUsage(); }

public: // manipulation functions
  void Reserve( size_t s )
  {
    edgelist.reserve( s );
    nodelist.reserve( s );
  }

  void Flush()
  {
    for( int i=0; i<edgelist.size(); ++i ) { FreeEdge(edgelist[i]); }
    for( int i=0; i<nodelist.size(); ++i ) { FreeNode(nodelist[i]); }
    edgelist.clear();
    nodelist.clear();
  }

  void Swap( ListGraph& g )
  {
    nodelist.swap( g.nodelist );
    edgelist.swap( g.edgelist );
  }

	EdgeT* InsertEdge( NodeT* &vc1, NodeT* &vc2 )
  {
    if( !vc1 ) vc1 = InsertNode();
    if( !vc2 ) vc2 = InsertNode();
    EdgeT* v = AllocEdge();  
    myAssert( v && vc1 && vc2 );
    edgelist.indexed_push_back( v );
    v->AttachNodes( vc1, vc2 );
    return v;
  }

	NodeT* InsertNode()
  {
    NodeT* vc = AllocNode();
    myAssert( vc );
    nodelist.indexed_push_back( vc );
    return vc;
  }

  void DeleteEdgeButDoNotTouchNodePointers(EdgeT* &v)
  {
    myAssert( this->isInGraph(v) );
    edgelist.indexed_remove( v );
    FreeEdge( v );
  }
  
  void DeleteEdge( EdgeT* &v, bool bDeleteNodes = true )
  {
#ifdef ASSERT_AT_DELETE_UNKNOWN_ELEMENT  
    myAssert( this->isInGraph(v) );
#else
    if( this->IsInGraph(v) )
#endif        
    {
      v->NodeA()->RemoveEdge( v );
      if( bDeleteNodes && v->NodeA()->Count()<=0 ) {
        DeleteUnusedNode( v->NodeA() );
      }
      v->NodeB()->RemoveEdge( v );
      if( bDeleteNodes && v->NodeB()->Count()<=0 ) {
        DeleteUnusedNode( v->NodeB() );
      }  
      DeleteEdgeButDoNotTouchNodePointers(v);
    }
  }

  void DeleteUnusedNode( NodeT* vc )
  {
    myAssert( vc->Count()==0 );
#ifdef ASSERT_AT_DELETE_UNKNOWN_ELEMENT
    myAssert( this->isInGraph(vc) );
#else    
    if( this->isInGraph(vc) )
#endif        
    {
      nodelist.indexed_remove( vc );
      FreeNode( vc );
    }
  }

  void Optimize()
  {
    edgelist.sort_by_addr();
    nodelist.sort_by_addr();
  }

protected:
  EdgeT* AllocEdge()              { return new EdgeT; }
  void   FreeEdge( EdgeT* v )     { delete v; }
  NodeT* AllocNode()              { return new NodeT; }
  void   FreeNode( NodeT* n )     { delete n; }

  IndexObjVector<NodeT> nodelist;
  IndexObjVector<EdgeT> edgelist;
};


//-----------------------------------
//-----------------------------------


#endif

