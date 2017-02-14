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
#ifndef COMPRESSED_ROW_UNDIRECTED_GRAPH
#define COMPRESSED_ROW_UNDIRECTED_GRAPH

#include <stdio.h>
#include "dynamicarray.h"
#include "../common/common.h"
#include <boost/iterator/counting_iterator.hpp>
#include <boost/range/iterator_range.hpp>

class CompressedRows
{
    DynArray<int> offsets;
    DynArray<int> col_id;

    void initialize(int num_rows_) { offsets.clear(); offsets.resize(num_rows_+1); }
    void inc_count(int i, int n) { offsets[i+1] += n; }
    void end_count()
    {
      for (int i=1; i<offsets.size(); ++i)
      {
        offsets[i] += offsets[i-1];
      }
      col_id.clear();
      col_id.resize(offsets[offsets.size()-1], -1);
    }

  public:
    CompressedRows() : allow_identical_edges(false) {}

    bool allow_identical_edges;
    
    void reinit(int num_rows_, const std::vector<unsigned int> &col_counts)
    {
      initialize(num_rows_);
      std::copy(col_counts.begin(), col_counts.end(), offsets.begin()+1);
      end_count();
    }

    int add(int row, int col)
    {
      const int a = offsets[row];
      const int b = offsets[row+1];
      for(int p=a; p<b; ++p)
      {
        if (col_id[p] == -1) { col_id[p]=col; return p; }
        if (col_id[p] == col && allow_identical_edges) { return p; }
      }
      myAssert(false);
      return -1;
    }

    void compress()
    {
#ifdef DEBUG
      for(int i=0; i<col_id.size(); ++i)
        myAssert(col_id[i] != -1);
#endif
    }

    int n_rows() const { return offsets.size()-1; }
    int operator()(int row, int col) const
    {
      myAssert(allow_identical_edges == false);
      const int a = offsets[row];
      const int b = offsets[row+1];
      for(int p=a; p<b; ++p)
      {
        if (col_id[p] == col) return p;
      }
      return -1;
    }    
    // Vorsicht. Wenn ich die klasse hier durch dealii SparsityPatterns ersetze sind noch Diagonaleinträge mit drin.
    int n_nonzero_elements() const { return col_id.size(); }
    int row_length(int row) const { return offsets[row+1]-offsets[row]; }
    int nz_number(int row, int index) const { return offsets[row]+index; }
    int column_number(int row, int index) const { return col_id[nz_number(row, index)]; }

    void print(std::ostream &os) const
    {
      os << "nrows = " << n_rows() << std::endl;
      for (int i=0; i<n_rows(); ++i)
      {
        os << "row " << i <<"-col: ";
        for (int j=0; j<row_length(i); ++j)
          os << column_number(i, j) << ", ";
        os << std::endl;
        if (!allow_identical_edges)
        {
          os << "    " << i <<"-id: ";
          for (int j=0; j<row_length(i); ++j)
            os << (*this)(i, column_number(i,j)) << ", ";
          os << std::endl;
        }
      }
    }
};


class CompressedGraph
{
  CompressedRows sp;
  DynArray<int> nz_to_edge;
  DynArray<my::eqpair<int> > m_edges;
  // for construction
  int nv;
public:
  typedef int node_type;
  typedef int edge_type;

  void beginConstruct()
  {
    sp.allow_identical_edges = true;
    nv = 0;
    m_edges.reserve(128);
  }

  void add_edge(int va, int vb)
  {
    myAssert(va >= 0 && vb >= 0);
    nv = std::max(nv, std::max(va, vb));
    m_edges.push_back(my::make_eqpair(va, vb));
  }

  template<class Range>
  void add_edges(const Range &r)
  {
    m_edges.reserve(m_edges.size() + (r.end() - r.begin()));
    for (typename Range::const_iterator it = r.begin(); it != r.end(); ++it)
    {
      add_edge(it->first, it->second);
    }
  }

  void endConstruct()
  {
    nv += 1;
    std::vector<unsigned int> rowlengths(nv, 0); // initially 0 entries in the rows. each edge only adds off-diagonals now
    for(int i=0; i<m_edges.size(); ++i)
    {
      int a, b; std::tie(a,b) = m_edges[i];
      rowlengths[a]+=1;
      rowlengths[b]+=1;
    }
    sp.reinit(nv, rowlengths);

    nz_to_edge.resize(m_edges.size() * 2);
    
    for (int i=0; i<m_edges.size(); ++i)
    {
      int a, b; std::tie(a,b) = m_edges[i];
      int gid = sp.add(a, b);
      nz_to_edge[gid] = i;
      gid = sp.add(b, a);
      nz_to_edge[gid] = i;
    }
    sp.compress();
  }

  const DynArray<my::eqpair<int> >& edges() const { return m_edges; }
  
  int num_edges() const { return m_edges.size(); }
  int num_nodes() const { return sp.n_rows(); }

  int degree(int v) const { return sp.row_length(v); }
  int adjacent_edge(int v, int i) const { assert(v>=0 && v<num_nodes()); assert(i < degree(v)); return nz_to_edge[sp.nz_number(v, i)]; }
  int adjacent_node(int v, int i) const { assert(v>=0 && v<num_nodes()); assert(i < degree(v)); return sp.column_number(v, i); }
  NodeNeighbor<int, int> adjacent(int v, int i) const { return make_node_neighbor(adjacent_node(v, i), adjacent_edge(v, i)); }

  my::eqpair<int> edge_nodes(int e) const { return m_edges[e]; }

  void print(std::ostream &os) const
  {
    sp.print(os);
    os << "v: " << num_nodes() << ", e: " << num_edges() << std::endl;
    for (int i=0; i<num_nodes(); ++i)
    {
      os << "v " << i << "-(v,e): ";
      for (int j=0; j<degree(i); ++j)
      {
        NodeNeighbor<int, int> nb = adjacent(i, j);
        os << boost::format("(%i, %i), ") % nb.node % nb.edge;
      }
      os << std::endl;
    }
  }

  //typedef DynArray<my::eqpair<int> >::iterator edge_iterator;
  //typedef DynArray<my::eqpair<int> >::const_iterator const_edge_iterator;

  typedef boost::counting_iterator<int> const_edge_iterator;
  boost::iterator_range<const_edge_iterator> edge_range() const { return boost::make_iterator_range(const_edge_iterator(0), const_edge_iterator(num_edges())); }

  //boost::iterator_range<edge_iterator> edge_range() { return boost::make_iterator_range(m_edges.begin(), m_edges.end()); }
  //boost::iterator_range<const_edge_iterator> edge_range() const { return boost::make_iterator_range(m_edges.begin(), m_edges.end()); }

  typedef boost::counting_iterator<int> const_node_iterator;
  boost::iterator_range<const_node_iterator> node_range() const { return boost::make_iterator_range(const_node_iterator(0), const_node_iterator(num_nodes())); }
};


class IntegerMap
{
    /* maps integers from an intervall to other integers */
    DynArray<int> vmap;
    int n;
  public:
    IntegerMap() : n(0) {}
    
    void resize(int vertex_count) { vmap.resize(vertex_count, unused()); }

    static int unused() { return -std::numeric_limits<int>::max(); }

    int add(int x)
    {
      myAssert(x < vmap.size());
      if (vmap[x] != unused()) return vmap[x];
      myAssert(vmap[x] == unused());
      vmap[x] = n++;
      return n-1;
    }

    int operator[](int x) const { myAssert(x < vmap.size()); return vmap[x]; }
    int num_added() const { return n; }
};



#endif
