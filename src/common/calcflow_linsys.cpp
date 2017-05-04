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

#include "calcflow_common.h"
#include "trilinos_linsys_construction.h"
#include "mwlib/ptree_ext.h"

// trilinos stuff
//see trilinos_linsys_construction.h
//#include <Epetra_SerialComm.h>
// is include in trilinos_linsys_construction.h
#include <Epetra_CrsMatrix.h>
#include <Epetra_LinearProblem.h>
//#include <AztecOO.h>
#include <ml_epetra_preconditioner.h>
#include <EpetraExt_RowMatrixOut.h>
#include <EpetraExt_VectorOut.h>
#include <mwlib/timer.h>
#include "calcflow_linsys.h"

using boost::property_tree::ptree;
using boost::property_tree::make_ptree;

void Linsys::initialize_pattern(int num_vertices, const std::vector< my::eqpair< int > >& edges)
{
    std::vector<int> column_cnt(num_vertices, 1);
    for (EdgeIter it = edges.begin(); it != edges.end(); ++it)
    {
      column_cnt[it->first]++;
      column_cnt[it->second]++;
    }
#ifdef EPETRA_MPI
#ifdef DEBUG //some mpi error left in ubuntu 16
    std::printf("is the bad thing happening here?");
#endif
    Epetra_MpiComm epetra_comm(MPI_COMM_SELF);
#endif
#ifndef EPETRA_MPI
    Epetra_SerialComm epetra_comm;
#endif
    //epetra_comm.Print(cout); cout << endl;
    Epetra_Map epetra_map(num_vertices, 0, epetra_comm);
    SparsityPattern sp(Copy, epetra_map, get_ptr(column_cnt), true);

    for (int i=0; i<num_vertices; ++i)
    {
      sp.InsertGlobalIndices(i, 1, &i);
    }
    for (int i=0; i<edges.size(); ++i)
    {
      int v[2] = { edges[i].first, edges[i].second };
      sp.InsertGlobalIndices(v[0], 1, &v[1]);
      sp.InsertGlobalIndices(v[1], 1, &v[0]);
    }
    sp.FillComplete();
    sp.OptimizeStorage();

    sys.reset(new Epetra_CrsMatrix(Copy, sp));
    rhs.reset(new Epetra_Vector(epetra_map));
    lhs.reset(new Epetra_Vector(epetra_map));
    //scaling_const = 1.;
}
void Linsys::sys_add(int a, int b, double val)
{
  sys->SumIntoMyValues(a, 1, &val, &b);
}
void Linsys::sys_set(int a, int b, double val)
{
  sys->ReplaceMyValues(a, 1, &val, &b);
}
void Linsys::rhs_add(int a, double val)
{
  (*rhs)[a] += val;
}
void Linsys::rhs_set(int a, double val)
{
  (*rhs)[a] = val;
}
double Linsys::lhs_get(int a)
{
  return (*lhs)[a]; 
}
void Linsys::clear_values()
{
  sys->PutScalar(0.);
  rhs->PutScalar(0.);
}
void Linsys::fill_edge_values(int a, int b, FlReal cond, const Linsys::BcsMap& bcs)
{
  fill_edge_values(a, b, -cond, cond, 0, 0, bcs);
}
void Linsys::fill_edge_values(int a, int b, FlReal coeff_diag, FlReal coeff_offdiag, FlReal rhsa, FlReal rhsb, const Linsys::BcsMap& bcs)
{
  rhs_add(a, -rhsa);
  rhs_add(b, -rhsb);
  {
    sys_add(a, a, coeff_diag);
    sys_add(b, b, coeff_diag);
    sys_add(a, b, coeff_offdiag);
    sys_add(b, a, coeff_offdiag);
  }
}
void Linsys::fill_values(std::vector< my::eqpair< int > >& edges, const FlArray& cond, const Linsys::BcsMap& bcs)
{
  // sum_i w_i (p_i - p_0) = in_flussW
  for (int i=0; i<edges.size(); ++i)
  {
    fill_edge_values(edges[i].first, edges[i].second, cond[i], bcs);
  }
  end_filling(bcs);
}
void Linsys::end_filling(const Linsys::BcsMap& bcs)
{
    //scaling_const /= sys->NumGlobalRows();
#if 1
#ifdef DEBUG
    printf("bcs.size(): %i\n",bcs.size());
    for(BcsMap::const_iterator bc = bcs.begin(); bc != bcs.end(); ++bc)
    {
      printf("first: %i, second: %f\n", bc->first, bc->second.val);
    }
#endif
#endif
    for(BcsMap::const_iterator bc = bcs.begin(); bc != bcs.end(); ++bc)
    {
      // Add elements to the matrix at the place of the boundary node. So boundary nodes get their value also during the matrix solve.
      int id = bc->first;
      const FlowBC &bcx = bc->second;
      if (bcx.type == FlowBC::PIN)
      {
        int cnt;
        #if 0
        // leave the matrix unsymmetric -> probably bad
        // note: it IS bad
        double *values;
        sys->ExtractGlobalRowView(id, cnt, values);
        for (int i=0; i<cnt; ++i) values[i] = 0.;  // clear all row values
        sys_set(id, id, 1. * -scaling_const); // set diagonal, must be negative, otherwise the matrix is not pos. definite
        rhs_set(id, bcx.val * -scaling_const);
        #else
        double *values;
        int *indices;
        sys->ExtractMyRowView(id, cnt, values, indices);
        for (int i=0; i<cnt; ++i)
        {
          values[i] = 0.;  // clear all row values
          // remove entry from other rows
          int othercnt, otherid = indices[i];
          double *other_values;
          int *other_indices;
          sys->ExtractMyRowView(otherid, othercnt, other_values, other_indices);
          for (int j=0; j<othercnt; ++j)
          {
            if (other_indices[j] == id)
            {
              double coeff = other_values[j];
              other_values[j] = 0.;  // no contribution from pinned value
              rhs_add(otherid, -coeff * bcx.val);
            }
          }
        }
        sys_set(id, id, 1. * -scaling_const); // set diagonal, must be negative, otherwise the matrix is not pos. definite
        rhs_set(id, bcx.val * -scaling_const);
        #endif
      }
      else if (bcx.type == FlowBC::CURRENT)
      {
        rhs_add(id, bcx.val);
      }
      else if (bcx.type == FlowBC::RESIST)
      {
        sys_add(id, id, -bcx.w);
        rhs_add(id, -bcx.w * bcx.val);
      }
    }
    sys->Scale(-1. / scaling_const); // reversed sign for pos definiteness
    rhs->Scale(-1. / scaling_const);
}
void Linsys::zero_out()
{
  sys->PutScalar(0.);
  rhs->PutScalar(0.);
  lhs->PutScalar(0.);
}
void Linsys::solve()
{
  /*
    * some major trilinos troubles appeared here with version 12.2 currently (14.04.2016)
    * running on lusi cluster.
    * mw suggested to switch solver from cg to bicgstab
    * we installed trilinos 11.0.3 and it seems to work again
    * NOTE: probably everythin was due to MKL incompability
    * solution: only use intel compiler in future
    */
  ptree pt = make_ptree("output", true)("max_iter", 1000)("throw",false)("conv","rhs")("max_resid",1.e-14)("solver","cg")("preconditioner","multigrid")("use_smoothed_aggregation", false);
  int returnSolve = SolveEllipticEquation(sys, rhs, lhs, pt);
}


// struct Linsys
// {
//   
// 
//   Linsys() : scaling_const(0.) {}
// 
// };
