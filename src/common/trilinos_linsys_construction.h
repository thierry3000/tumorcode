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
#ifndef TRILINOS_SYSTEM_CONSTRUCTION_H
#define TRILINOS_SYSTEM_CONSTRUCTION_H

#include "mwlib/helpers-vec.h"
#include "mwlib/lattice-data.h"
#include "mwlib/dynamicarray.h"
#include "mwlib/field.h"
#include "continuum-flow.h"

#if (defined __GNUC__) && !(defined __INTEL_COMPILER)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wextra"
#endif
#include <Epetra_ConfigDefs.h> //knows if EPETRA_MPI is there!
#ifdef EPETRA_MPI
    #include <Epetra_MpiComm.h> //import mpi on its own!
#endif
#ifndef EPETRA_MPI
  #include <Epetra_SerialComm.h>
#endif
#include <Epetra_CrsMatrix.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_MultiVector.h>
#include <Epetra_BasicRowMatrix.h>
#include <ml_epetra_preconditioner.h>
#include <Ifpack.h>

#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_Version.hpp>
#include <BelosTpetraAdapter.hpp>

#include <BelosSolverFactory.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosEpetraAdapter.hpp>

typedef double BelosScalarType;
typedef Epetra_MultiVector BelosMultiVector;
typedef Epetra_Operator BelosOperator;
typedef double                            ST;
typedef Teuchos::ScalarTraits<ST>        SCT;
typedef SCT::magnitudeType                MT;
typedef Epetra_MultiVector                MV;
typedef Epetra_Operator                   OP;
typedef Belos::MultiVecTraits<ST,MV>     MVT;
typedef Belos::OperatorTraits<ST,MV,OP>  OPT;


#if (defined __GNUC__) && !(defined __INTEL_COMPILER)
#pragma GCC diagnostic pop
#endif
#include <boost/property_tree/ptree.hpp>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/scoped_ptr.hpp>

#if 1
#include <boost/foreach.hpp>


inline Array3d<double> Array3dFromEpetraMultiVector(Epetra_MultiVector &mv, int vec_idx, const LatticeDataQuad3d &ld)
{
  // works only when number of sites or array entries is less than 2**31!
  return Array3d<double>(ld.Box(), ld.Strides().cast<int>(), mv[vec_idx], ld.NumSites(), false);
}

inline ConstArray3d<double> ConstArray3dFromEpetraMultiVector(const Epetra_MultiVector &mv, int vec_idx, const LatticeDataQuad3d &ld)
{
  // works only when number of sites or array entries is less than 2**31!
  // otherwise the stride type in Array3d must be changed.
  return ConstArray3d<double>(ld.Box(), ld.Strides().cast<int>(), const_cast<Epetra_MultiVector&>(mv)[vec_idx], ld.NumSites(), false);
}

template<class T>
inline Array3d<T> SetupBufferWithNeumannBC(T* storage, std::size_t storage_size, ConstArray3d<T> source_data, const BBox3 &buffer_bb, int dim, const BBox3 &full_bb)
{
  BBox3 full_buffer_box = ExtendForDim(buffer_bb, dim, 1);
  Array3d<double> buffer(full_buffer_box, Int3(-1), storage, storage_size, false);
  BBox3 storage_box_wo_boundary(full_buffer_box);
  storage_box_wo_boundary.Intersection(full_bb);
  buffer.setBox(storage_box_wo_boundary);
  buffer[storage_box_wo_boundary].fill(source_data[storage_box_wo_boundary]);
  //CopyBorder(buffer, true, dim>1, dim>2, full_buffer_box);
  CopyBorder(buffer, buffer_bb, full_bb, dim, 1);
  return buffer;
}



template<class T>
struct MapNeZero_
{
  typedef bool result_type;
  inline void inplace(bool &b, const T &a) const
  {
    b = Abs(a) > T(1.e-18);
  }
};

template<class T>
inline Array3d<bool> WhereNonZero(const ConstArray3d<T> a)
{
  Array3d<bool> b; b.initFromBox(a.Box());
  _ArrayNdInternal::apply(b, a, MapNeZero_<T>());
  return b;
}
#endif


struct FiniteVolumeMatrixBuilder : public boost::noncopyable
{
  enum BoundaryCondition {
    NEUMANN = 0,
    DIRICHLET_X = 1,
    DIRICHLET_YZ = 2,
  };
 
//   std::auto_ptr<Epetra_CrsMatrix> m;
//   std::auto_ptr<Epetra_Vector> rhs;
  Teuchos::RCP<Epetra_CrsMatrix> m;
  Teuchos::RCP<Epetra_Vector> rhs;
  LatticeDataQuad3d ld;
  //BBox3 bb, bb_inner, boxes[9];
  //DynArray<Int3> stencil_points;
  int dim;

  void Init7Point(const LatticeDataQuad3d& ld_, int dim_); // in 3d, 5 in 2d
  void Init27Point(const LatticeDataQuad3d& ld_, int dim_); // 3*3*3 in 3d, 9 in 2d

  template<class ArrayType>
  void AddDiffusion(const BBox3 &bb, const ArrayType &conductivities, double prefactor=1.);
  template<class FaceVarFunctor>
  void AddDiffusion2(const BBox3& bb, int dim, const FaceVarFunctor& conductivities, double prefactor=1., bool isotropic=false, const Array3d<float> local_prefactors = Array3d<float>()); // uses face centered conductivities
  void AddLocally(const Int3 &p, double diagonal, double rhs);
  void ZeroOut();
  void Add(int site_row, int site_col, double value);
  void AddRhs(int site, double value);
  void SetDirichletBoundaryConditions(const BBox3 &bb, int bc, double boundary_value);
#if 1
  void InitFromStencil(const LatticeDataQuad3d &ld, int dim, const ConstArray3d<bool> &pattern);
  //void AddStencil(const Int3 &pos, const ConstArray3d<double> &stencil, const T rhs);
#endif
};



template<class ArrayType>
inline void FiniteVolumeMatrixBuilder::AddDiffusion(const BBox3& bb, const ArrayType& conductivities, double prefactor)
{
  /*this is the 7 point approximation of the diffusion operator nabla * D * nabla, with variable D.*/
  /*warning: neumann boundary condition with non-zero gradient not implemented */
  const int dim = 3; // dimension
  prefactor /= (ld.Scale() * ld.Scale());
  DynArray<int> columns(2 * dim + 1);
  DynArray<double> values(2 * dim + 1);
  FOR_BBOX3(p,bb)
  {
    int site = ld.LatticeToSite(p);
    columns[0] = site;
    values[0] = 0.;
    int i=1;
    double k0 = conductivities(p);
    for(int d=0; d<2*dim; ++d)
    {
      Int3 nbp = ld.NbLattice(p,d);
      if(!ld.IsInsideLattice(nbp))continue;

      int nbsite = ld.LatticeToSite(nbp);
      double k1 = conductivities(nbp);
      double c = CalcInterfaceCondCoeff(k0,k1); //2.0*(k0*k1)/(k0+k1+std::numeric_limits<double>::epsilon());
      c *= prefactor;
      columns[i] = nbsite;
      values[0] -= c;
      values[i]  = c;
      ++i;
    }
    m->SumIntoMyValues(site, i, get_ptr(values), get_ptr(columns));
  }
}

/**
 * \brief face centered conductivities
 */
template<class FaceVarFunctor>
inline void FiniteVolumeMatrixBuilder::AddDiffusion2(const BBox3& bb, int dim, const FaceVarFunctor& conductivities, double prefactor, bool isotropic, const Array3d< float > local_prefactors)
{
  /*this is the 7 and 9 point approximation of the diffusion operator nabla * D * nabla, with variable D.*/
  prefactor = prefactor/(ld.Scale() * ld.Scale());
  const BBox3 ldbox = ld.Box();
  
  if (isotropic)
  {
    if (!local_prefactors.empty()) throw std::runtime_error("not implemented");
    
    Int3 ss(0); 
    Array3d<float> diff_stencils[3][2];
    for (int axis=0; axis<dim; ++axis)
    {
      Array3d<float> s = Stencils::Diff1<float>(axis, dim, 1., Stencils::W3_ISOTROPIC, Stencils::FWD);
      s.move(-(s.size()-Int3(1))/2);
      diff_stencils[axis][1].initCopy(s);
      Int3 off(0); off[axis]=-1;
      s.move(off);
      s *= -1.;
      diff_stencils[axis][0] = s;
      ss[axis] = 1;
    }
    LatticeIndexing<3, int> stencil_lattice;
    stencil_lattice.Init(BBox3(-ss[0], -ss[1], -ss[2], ss[0], ss[1] ,ss[2]));
    stencil_lattice.MoveSiteRange(-ss);

//     cout << prefactor << endl;
//     cout << stencil_lattice.LatticeToSite(-ss) << endl;
//     cout << stencil_lattice.LatticeToSite(ss) << endl;
//     for (int i=0; i<dim; ++i)
//     {
//       PrintArray(diff_stencils[i][0]);
//       PrintArray(diff_stencils[i][1]);
//     }
    Array3d<float> stencil_buffer(stencil_lattice.Box());
    
    DynArray<int> columns(my::ipow(3, dim));
    DynArray<double> values(my::ipow(3, dim));
    FOR_BBOX3(p,bb)
    {
      columns.fill(-1);
      values.fill(0);
      double stencil_rhs = 0.;
      stencil_buffer.fill(0);
      for (int axis=0; axis<dim; ++axis)
      {
        for (int side=0; side<=1; side+=1)
        {
          const Int3 nbp = ld.NbLattice(p, axis*2+side);
          if(!ld.IsInsideLattice(nbp)) continue;
          //Int3 fp(p); fp[axis]+=side;
          double c = conductivities(p, nbp, axis, side);
          FOR_BBOX3(q, diff_stencils[axis][side].getBox())
          {
            stencil_buffer(q) += c * diff_stencils[axis][side](q);
          }
        }
      }

      int num_entries = 0;
      FOR_BBOX3(q, stencil_lattice.Box())
      {
        Int3 pp = q+p;
        for (int axis=0; axis<dim; ++axis)
        {
          if (pp[axis] > ldbox.max[axis]) { --pp[axis];  }
          if (pp[axis] < ldbox.min[axis]) { ++pp[axis];  }
        }
        int index = stencil_lattice.LatticeToSite(pp-p);
        //if (is_bc) // zero gradient/flux bc
        //{
          //stencil_buffer(pp-p) += stencil_buffer(q);
          //stencil_buffer(q) = 0;
        num_entries += columns[index]<0 ? 1 : 0;
          values[index] += prefactor * stencil_buffer(q);
          columns[index] = ld.LatticeToSite(pp);
        //}
      }
//       cout << "point " << p << " -> " << ld.LatticeToSite(p) << endl;
//       for (int k=0; k<columns.size(); ++k)
//       {
//         cout << format("(%i, %f), ") % columns[k] % values[k];
//       }
//      cout << endl;
      if (num_entries != columns.size())
      {
        num_entries=0;
        for (int k=0; k<columns.size(); ++k)
        {
          if (columns[k]<0) continue;
          columns[num_entries] = columns[k];
          values[num_entries] = values[k];
          ++num_entries;
        }
      }
//       cout << "point " << p << " -> " << ld.LatticeToSite(p) << endl;
//       for (int k=0; k<num_entries; ++k)
//       {
//         cout << format("(%i, %f), ") % columns[k] % values[k];
//         myAssert(columns[k]>=0 && columns[k]<m->NumMyCols());
//       }
//      cout << endl;
      m->SumIntoMyValues(ld.LatticeToSite(p), num_entries, get_ptr(values), get_ptr(columns));
    }
  }
  else
  {
    bool use_local_prefactor = !local_prefactors.empty();
    
    DynArray<int> columns(2 * dim + 1);
    DynArray<double> values(2 * dim + 1);
    FOR_BBOX3(p,bb)
    {
      int site = ld.LatticeToSite(p);
      columns[0] = site;
      values[0] = 0.;
      int i=1;
      double loc_prefactor0 = prefactor*(use_local_prefactor ? local_prefactors(p) : 1.);
      
      for (int axis=0; axis<dim; ++axis)
      {
        for (int side=0; side<=1; side+=1)
        {
          int dir = axis*2 + side;
          Int3 nbp = ld.NbLattice(p,dir);
          if(!ld.IsInsideLattice(nbp)) continue;
          int nbsite = ld.LatticeToSite(nbp);
          //Int3 fp(p); fp[axis]+=side;
          double c = conductivities(p, nbp, axis, side);
          double loc_prefactor1 = prefactor*(use_local_prefactor ? local_prefactors(nbp) : 1.);
          values[0] -= c*loc_prefactor0;
          values[i]  = c*loc_prefactor1;
          columns[i] = nbsite;
          ++i;
        }
      }
      m->SumIntoMyValues(site, i, get_ptr(values), get_ptr(columns));
    }
  }
}


inline void FiniteVolumeMatrixBuilder::ZeroOut()
{
  m->PutScalar(0.);
  rhs->PutScalar(0.);
}

inline void FiniteVolumeMatrixBuilder::AddLocally(const Int3 &p, double diagonal, double rhsvalue)
{
  const int site = ld.LatticeToSite(p);
  m->SumIntoMyValues(site, 1, &diagonal, &site);
  (*rhs)[site] += rhsvalue;
}

inline void FiniteVolumeMatrixBuilder::Add(int site_row, int site_col, double value)
{
  m->SumIntoMyValues(site_row, 1, &value, &site_col);
}

inline void FiniteVolumeMatrixBuilder::AddRhs(int site, double value)
{
  (*rhs)[site] += value;
}


class ConvergenceFailureException : public std::runtime_error
{
public:
  ConvergenceFailureException(const std::string &msg, int reason_) : std::runtime_error(msg), reason(reason_) {}
  enum Reason
  {
    OTHER = 0,
    MAX_ITERATIONS = 1
  };
  int reason;
};



/*
  parameters:
  verbosity:
    "verbose", "none", "normal"
  output (bool)
  use_multigrid true
  max_levels 10
  max_iter 50
  max_resid 1.e-9
  keep_preconditioner false
*/
#if 1
class EllipticEquationSolver : boost::noncopyable
{
  typedef double scalar_type;
  typedef int local_ordinal_type;
  typedef long global_ordinal_type;
  typedef KokkosClassic::DefaultNode::DefaultNodeType node_type;
  
  Teuchos::RCP<Epetra_CrsMatrix> sys_matrix;
  Teuchos::RCP<const Epetra_Vector> rhs;
  
  typedef Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type> vec_type;
  typedef Tpetra::CrsMatrix<scalar_type, local_ordinal_type, global_ordinal_type, node_type> matrix_type;
  typedef Tpetra::Operator<scalar_type, local_ordinal_type, global_ordinal_type, node_type> op_type;
    
  typedef Belos::LinearProblem<BelosScalarType, BelosMultiVector, BelosOperator> BelosLinearProblem;
  Teuchos::RCP<BelosLinearProblem> problem;
  //Teuchos::RCP<Ifpack_Preconditioner> ifpack_Prec;
  Teuchos::RCP<Belos::EpetraPrecOp> belos_Prec;
  Teuchos::RCP<ML_Epetra::MultiLevelPreconditioner> ml_prec;
  Teuchos::RCP<Teuchos::ParameterList> belosList;
  
  ptree params;
  bool keep_preconditioner;
public:
  int solve(const Teuchos::RCP<Epetra_Vector> &lhs);
  double time_precondition, time_iteration;
  int iteration_count;
  bool convergent;
  double residual;
  int init(const Teuchos::RCP<Epetra_CrsMatrix> &_matrix, const Teuchos::RCP<const Epetra_Vector> &_rhs, const boost::property_tree::ptree &_params);
};
#endif

int SolveEllipticEquation(Teuchos::RCP <Epetra_CrsMatrix> &matrix, Teuchos::RCP<Epetra_Vector> &rhs, Teuchos::RCP<Epetra_Vector> &lhs, const boost::property_tree::ptree &params);

/*------------------------------------------------------
------------------------------------------------------*/

typedef boost::function3<void, int, const BBox3&, FiniteVolumeMatrixBuilder&> DiffSolveBuildFuncType;
template<class T>
void StationaryDiffusionSolve(const LatticeDataQuad3d &ld,
                              const DynArray<BBox3> &mtboxes,
                              int dim,
                              Array3d<T> result,
                              DiffSolveBuildFuncType buildfunc,
                              const ptree &pt_params = ptree());

template<class T>
void StationaryDiffusionSolve(const ContinuumGrid &grid,
                              const DomainDecomposition &mtboxes,
                              Array3d<T> result,
                              DiffSolveBuildFuncType buildfunc,
                              const ptree &pt_params = ptree());


#endif
