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
#ifndef CALCFLOW_LINSYS_H
#define CALCFLOW_LINSYS_H

#include "calcflow_common.h"
#include "trilinos_linsys_construction.h"
#include "mwlib/ptree_ext.h"

#include <Epetra_CrsMatrix.h>
#include <Epetra_LinearProblem.h>
#include <ml_epetra_preconditioner.h>
#include <EpetraExt_RowMatrixOut.h>
#include <EpetraExt_VectorOut.h>
#include <mwlib/timer.h>
#include <boost/unordered_map.hpp>

using boost::property_tree::ptree;
using boost::property_tree::make_ptree;

struct FlowBC
{
  enum Type{ PIN = 1, CURRENT = 2, RESIST = 3};
  //enum Pin { PIN = 1 };  // set the node to a fixed blood pressure
  //enum Current { CURRENT = 2 }; // a fixed blood flow rate into the node
  //enum Resist { RESIST = 3 }; // a resistor between the node and a fixed blood pressure
  FlowBC() : typeOfInstance(Type::PIN),w(0.),val(0.) {}
  FlowBC(Type theType, double val_) : w(0),val(val_),typeOfInstance(theType) {}
  FlowBC(Type theType, double cond_, double val_) : w(cond_), val(val_), typeOfInstance(theType) {}
  //FlowBC_pin(double val) : w(0),val(val),typeOfInstance(Type::PIN) {}
  //FlowBC_resit(double w, double val) : w(w),val(val),typeOfInstance(Type::RESIST) {}
  Type typeOfInstance; // PIN, CURRENT or RESIST
  double w,val; // w = flow conuctivity (?) of series "resistor", val = either blood pressure or flow rate
//   friend class boost::serialization::access;
//   template<class Archive>
//   void serialize(Archive &ar, const unsigned int version)
//   {
//       // save/load base class information
//       ar & type & w & val;
//   }
};

typedef DynArray<FlReal> FlArray;
typedef DynArray<my::eqpair<int> > FeArray;
typedef DynArray<my::Bitfield<uchar> > FbArray;
typedef boost::unordered_map<int, FlowBC> FlBoundaryList;

struct Linsys
{
  
  typedef Epetra_CrsGraph SparsityPattern;
  //typedef Epetra_CrsMatrix Matrix;
  //typedef Epetra_Vector Vector;

  typedef std::vector<my::eqpair<int> > EdgeRange;
  typedef EdgeRange::const_iterator EdgeIter;
  typedef boost::unordered_map<int, FlowBC> BcsMap;

  Teuchos::RCP<Epetra_CrsMatrix> sys;
  Teuchos::RCP<Epetra_Vector> rhs,lhs;
  //Epetra_CrsMatrix *sys;
  //Epetra_Vector *rhs, *lhs;
  double scaling_const;
  Teuchos::RCP<SparsityPattern> sp;
  Teuchos::RCP<Epetra_Map> epetra_map;
  //SparsityPattern *sp;
  //Epetra_Map *epetra_map;
  std::vector<int> column_cnt;

  Linsys();
  ~Linsys();
  //Linsys(){} // void constructor make sure you initialize everything in initialize_pattern
  
  void initialize_pattern(int num_vertices, const std::vector<my::eqpair<int> > &edges);

  void sys_add(int a, int b, double val);
  void sys_set(int a, int b, double val);
  void rhs_add(int a, double val);
  void rhs_set(int a, double val);
  double lhs_get(int a);

  void clear_values();
  void fill_edge_values(int a, int b, FlReal coeff_diag, FlReal coeff_offdiag, FlReal rhsa, FlReal rhsb, const BcsMap &bcs);
  void fill_edge_values(int a, int b, FlReal cond, const BcsMap &bcs);
  void fill_values(std::vector<my::eqpair<int> > &edges, const FlArray &cond, const BcsMap &bcs);

  void end_filling(const BcsMap &bcs);
  void zero_out();
  void solve();
};

#endif //CALCFLOW_LINSYS_H
