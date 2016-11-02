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
#include "trilinos_linsys_construction.h"
#include <boost/scoped_array.hpp>


#if 1
template<class FaceVarFunctor, class DiagonalFunctor>
class TrilinosMatrixFreeEllipticOperator : boost::noncopyable, public Epetra_BasicRowMatrix
{
  Epetra_Map epetra_map;
  LatticeDataQuad3d ld;
  int dim;
  const DomainDecomposition &mtboxes;
  int max_mt_box_size;
public:
  const FaceVarFunctor face_coefficient_functor;
  const DiagonalFunctor diagonal_functor;
  const double prefactor;

  TrilinosMatrixFreeEllipticOperator(const LatticeDataQuad3d &ld_, int dim_, const DomainDecomposition &mtboxes_, double prefactor_, const FaceVarFunctor &face_coefficient_functor_, const DiagonalFunctor &diagonal_functor_, Epetra_Map &map_)
    : epetra_map(map_), face_coefficient_functor(face_coefficient_functor_), ld(ld_), dim(dim_), mtboxes(mtboxes_), prefactor(prefactor_), diagonal_functor(diagonal_functor_), Epetra_BasicRowMatrix(map_.Comm())
  {
    Epetra_BasicRowMatrix::SetMaps(map_, map_);
    max_mt_box_size = 0;
    for (int i=0; i<mtboxes.size(); ++i)
    {
      max_mt_box_size = std::max(max_mt_box_size, Volume(ExtendForDim(mtboxes[i], dim, 1)));
    }
  }

  virtual int   ExtractMyRowCopy (int MyRow, int Length, int &NumEntries, double *Values, int *Indices) const
  {
    if (Length < 1 + 2 * dim)
    {
      NumEntries = 1+2*dim;
      return -1;
    }
    
    const BBox3 box = ld.Box();
    const Int3 p = ld.SiteToLattice(MyRow);
    const LatticeDataQuad3d::SVec strides = ld.Strides();
    const double f = prefactor/my::sqr(ld.Scale());
    const int p_offset = MyRow;
    
    NumEntries = 1;
    double diag_value = prefactor*diagonal_functor(p);
    for (int axis=0; axis<dim; ++axis)
    {
      if (p[axis]>box.min[axis])
      {
        double fco = f * face_coefficient_functor(axis, p);
        int q_offset = p_offset-strides[axis];
        diag_value -= fco;
        myAssert(NumEntries < Length);
        Values[NumEntries] = fco;
        Indices[NumEntries] = q_offset;
        ++NumEntries;
      }
      if (p[axis]<box.max[axis])
      {
        Int3 q(p); ++q[axis];
        double fco = f * face_coefficient_functor(axis, q);
        int q_offset = p_offset+strides[axis];
        diag_value -= fco;
        myAssert(NumEntries < Length);
        Values[NumEntries] = fco;
        Indices[NumEntries] = q_offset;
        ++NumEntries;
      }      
    }
    Values[0] = diag_value;
    Indices[0] = p_offset;
    return 0;
  }



  virtual int Multiply (bool TransA, const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
  {
    /* laplacian with no-flux boundary conditions
     */
    #pragma omp parallel
    {
      int storage_buffer_size = max_mt_box_size;
      boost::scoped_array<double> storage_array(new double[storage_buffer_size]);
      const BBox3 box = ld.Box();
      int length = X.MyLength();
      double h2inv = 1./my::sqr(ld.Scale());
      const int number_of_vectors = X.NumVectors();
      myAssert(length == ld.NumSites());

      BOOST_FOREACH(const BBox3 &mt_box, mtboxes.getCurrentThreadRange())
      {
        for (int vec_idx=0; vec_idx < number_of_vectors; ++vec_idx)
        {
          ConstArray3d<double> viewX = ConstArray3dFromEpetraMultiVector(X, vec_idx, ld);
          Array3d<double> viewY = Array3dFromEpetraMultiVector(Y, vec_idx, ld);
          Array3d<double> buffer = SetupBufferWithNeumannBC(storage_array.get(), storage_buffer_size, viewX, mt_box, dim, box);

          const Int3 strides = buffer.strides();
          
          FOR_BBOX3(p, mt_box)
          {
            int p_offset = buffer.offset(p);
            double x_at_p = buffer.dereference(p_offset);
            double result = diagonal_functor(p) * x_at_p;
            for (int axis=0; axis<dim; ++axis)
            {
              {
                double fco = h2inv * face_coefficient_functor(axis, p);
                result += fco * (buffer.dereference(p_offset-strides[axis]) - x_at_p);
              }
              {
                Int3 q(p); ++q[axis];
                double fco = h2inv * face_coefficient_functor(axis, q);
                result += fco * (buffer.dereference(p_offset+strides[axis]) - x_at_p);
              }
            }
            viewY(p) = prefactor * result;
          }
        }
      }
    }
    return 0;
  }


  virtual int   NumMyRowEntries (int MyRow, int &NumEntries) const
  {
    const BBox3 box = ld.Box();
    const Int3 p = ld.SiteToLattice(MyRow);
    int cnt = 1;
    for (int axis=0; axis<dim; ++axis)
    {
      if (p[axis]>box.min[axis]) cnt++;
      if (p[axis]<box.max[axis]) cnt++;
    }
    NumEntries = cnt;
    return 1;
  }

  virtual int   ExtractMyEntryView (int CurEntry, double *&Value, int &RowIndex, int &ColIndex)
  {
    return -1;
  }

  virtual int   ExtractMyEntryView (int CurEntry, double const *&Value, int &RowIndex, int &ColIndex) const
  {
    return -1;
  }

  virtual int   LeftScale (const Epetra_Vector &x)
  {
    return -1;
  }

  virtual int   RightScale (const Epetra_Vector &x)
  {
    return -1;
  }

  virtual const char *  Label () const
  {
    return "TrilinosMatrixFreeEllipticOperator";
  }
};
#endif