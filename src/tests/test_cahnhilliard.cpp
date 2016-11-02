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

#include <algorithm>
#include <fstream>

#include <boost/range.hpp>
#include <boost/range/numeric.hpp>
#include <boost/range/algorithm/transform.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>


#include "field.h"
#include "globals.h"
#include "drawing.h"
#include "multigrid.h"

#include "amreqgen.h"
#include "continuumnumerics.h"
#include "matrixconstruction.h"
#include "amrheaders.h"
#include "algebraicmultigrid.h"
#include "testamrutils.h"


namespace boost
{
    namespace range_detail
    {
        template<class InputIterator1, class InputIterator2,  class InputIterator3, class Fn3>
        inline Fn3 for_each_impl(InputIterator1 first1, InputIterator1 last1,
                                 InputIterator2 first2, InputIterator2 last2,
                                 InputIterator3 first3, InputIterator3 last3,
                                 Fn3 fn)
        {
            for (; first1 != last1 && first2 != last2 && first3 != last3; ++first1, ++first2, ++first3)
            {
                fn(*first1, *first2, *first3);
            }
            return fn;
        }
    }

    namespace range
    {
        template<class SinglePassRange1, class SinglePassRange2,  class SinglePassRange3, class Fn3>
        inline Fn3 for_each(const SinglePassRange1& rng1, const SinglePassRange2& rng2,  const SinglePassRange2& rng3, Fn3 fn)
        {
            BOOST_RANGE_CONCEPT_ASSERT(( SinglePassRangeConcept<const SinglePassRange1> ));
            BOOST_RANGE_CONCEPT_ASSERT(( SinglePassRangeConcept<const SinglePassRange2> ));
            BOOST_RANGE_CONCEPT_ASSERT(( SinglePassRangeConcept<const SinglePassRange3> ));

            return ::boost::range_detail::for_each_impl(
                ::boost::begin(rng1), ::boost::end(rng1),
                ::boost::begin(rng2), ::boost::end(rng2),
                ::boost::begin(rng3), ::boost::end(rng3), fn);
        }

        template<class SinglePassRange1, class SinglePassRange2,  class SinglePassRange3, class Fn3>
        inline Fn3 for_each(SinglePassRange1& rng1, const SinglePassRange2& rng2,  const SinglePassRange3& rng3, Fn3 fn)
        {
            BOOST_RANGE_CONCEPT_ASSERT(( SinglePassRangeConcept<SinglePassRange1> ));
            BOOST_RANGE_CONCEPT_ASSERT(( SinglePassRangeConcept<const SinglePassRange2> ));
            BOOST_RANGE_CONCEPT_ASSERT(( SinglePassRangeConcept<const SinglePassRange3> ));

            return ::boost::range_detail::for_each_impl(
                ::boost::begin(rng1), ::boost::end(rng1),
                ::boost::begin(rng2), ::boost::end(rng2),
                ::boost::begin(rng3), ::boost::end(rng3), fn);
        }

    }
    using range::for_each;
};



void TestArrayMapping()
{
  Array3d<float> laplace = Stencils::Laplace<float>(3, 1., Stencils::W3_SPARSE);
  Array3d<bool> pattern = Amr::WhereNonZero(laplace);
  PrintArray(laplace);
  PrintArray(pattern);
}


#if 0
void Testfunc()
{
  int a = 5;
  struct F
  { // does not work :(
    void print() { cout << a << endl; }
  };
  F().print();
}
#endif

template<class T>
struct DivGradStencil
{
  const Amr::GridVarT<T> mobility;
  const Amr::Grid &grid;
  const int ndim;
  ConstArray3d<T> diff13;
  mutable Array3d<T> buffer;
  mutable Array3d<T> buffer2;
  Int3 m_size;
  const T scale;

  DivGradStencil(const Amr::GridVarT<T> mobility, const Amr::Grid &grid, int ndim) :
    mobility(mobility), grid(grid), ndim(ndim), scale(1./grid.spacing())
    {
      m_size = Int3(3, ndim>1 ? 3 : 1, ndim>2 ? 3 : 1);
      diff13 = Stencils::Diff1<T>(0, ndim, grid.spacing(), Stencils::W3_SPARSE, Stencils::FWD);
      buffer.init(Int3(2,1,1));
      buffer2.init(m_size);
    }

  void operator()(const Int3 &center, Array3d<T> out) const
  {
    T mobbuff[3];
    Assert(out.size() == m_size);
    for(int axis=0; axis<ndim; ++axis)
    {
      // 1st order forward differences (part of divergence operator), multiplied by average mobility between cells
      Int3 p(center);
      for(int i=0; i<3; ++i)
      {
        mobbuff[i] = mobility(p.replace(axis,center[axis]-1+i));
      }
      buffer(0) = -0.5*scale*(mobbuff[0]+mobbuff[1]);
      buffer(1) = 0.5*scale*(mobbuff[1]+mobbuff[2]);
      Stencils::Convolve<T>(buffer, diff13, buffer2);      
      if(axis>0)
      {
        buffer2.swapAxes(0, axis);
        out += buffer2;
        buffer2.swapAxes(0, axis);
      }
      else
      {
        out.fill(buffer2);
      }
    }
  }

  Array3d<T> operator()(const Int3 &center) const { 
    Array3d<T> out(m_size, Cons::DONT);
    (*this)(center, out);
    return out;
  }

  const Int3 size() const { return m_size; }

  typedef T value_type;
};


template<class T>
void CahnHilliardImplicitStep2(const Amr::Grid &grid, const Amr::GridVarT<T> &ucurr, const Amr::GridVarT<T> &umid, const Amr::BoundaryConditions<T> &bcs, const T eps, const T dt, Amr::GridVarT<T> &unext)
{ 
  using namespace boost;
  using namespace boost::lambda; 
  // building the linear time stepping scheme following "unconditionally gradient stable time marching the Cahn-Hilliard equation.ps"
  // the equation to solve is: (I + k eps^2 A_h^2 - k A_h D(Un)^2) Un+1 = (I - k A_h) Un,
  // where k is the time step
  const int ndim = 2;

  Amr::GridVarT<T> mobility(grid, 1, umid.cellCentered());
#if 1
  FOR_BBOX3(p, mobility.extendedIndexRange())
  {
    T val = umid(p);
    val = (T(1.)+val)*T(0.5);
    mobility(p) = Cut(val*(1.-val), T(0.), T(1.));
  }
#else
  boost::transform(umid[mobility.extendedIndexRange()], mobility[mobility.extendedIndexRange()].begin(),
    ret<T>(bind(Func::cut<T>(0.,1.), bind(unlambda(_1 * (1.-_1)), (1.+_1)*0.5)))
    );
#endif

  Amr::GridVarT<T> qvar = unext;
  struct StencilGen
  {
#if 0  
    struct AssignScaled
    {
      const T s;
      AssignScaled(const T s) : s(s) {}
      void inplace(T &dst,const T src) const { dst = s *src; }
    };

    struct AddComp2
    {
      const T s;
      AddComp2(const T s) : s(s) {}
      void inplace(T &dst, const T u, const T f) const { dst += s * u * u * f; }
    };

#endif  
    ConstArray3d<T> stencilBiLaplace, stencilLaplace;
    const T eps, spacing ,dt;
    const Amr::GridVarT<T> ucurr;
    Int3 center;
    BBox3 inner;
    DivGradStencil<T> divgrad;

    StencilGen(int ndim, const T eps, const Amr::Grid &grid, const T dt, const Amr::GridVarT<T> &ucurr, const Amr::GridVarT<T> &mobility) :
      eps(eps), spacing(grid.spacing()), dt(dt), ucurr(ucurr), divgrad(mobility, grid, ndim)
    {
      stencilLaplace = Stencils::Laplace<T>(ndim, spacing, Stencils::W3_SPARSE);
      stencilBiLaplace = Stencils::Convolve(stencilLaplace, stencilLaplace);
      center = stencilBiLaplace.size()/2;
      inner = stencilLaplace.getBox().Move((stencilBiLaplace.size()-stencilLaplace.size())/2);
    }
    const ConstArray3d<bool> getPattern() const
    {
      //PrintArray(stencilBiLaplace);
      ConstArray3d<bool> a(Amr::WhereNonZero(stencilBiLaplace));
      //PrintArray(a);
      return a;
    }    
    void getValues(const Int3 &pos, Array3d<T> &a) const
    { 
      // rather slow variants here
#if 0    
      const ConstArray3d<T> uslice = ucurr[Move(inner,pos-center)];
      a = (dt*eps*eps) * stencilBiLaplace;
      a(center) += 1.;
      a[inner] += (-dt)*((uslice*uslice)*stencilLaplace);
#endif
#if 0
      const ConstArray3d<T> uslice(ucurr[Move(inner,pos-center)]);
      _ArrayNdInternal::apply(a, stencilBiLaplace, AssignScaled(dt*eps*eps));
      a(center) += 1.;
      _ArrayNdInternal::apply(a[inner], uslice, stencilLaplace, AddComp2(-dt));
#endif
#if 0
      FOR_REG3V2(p, stencilBiLaplace.size())
      {
        a(p) = (dt * eps * eps) * stencilBiLaplace(p);
      }
      a(center) += 1.;
      
      const ConstArray3d<T> uslice(ucurr[Move(inner,pos-center)]);
      Array3d<T> aslice(a[inner]);
      FOR_REG3V2(p, stencilLaplace.size())
      {
        aslice(p) -= dt * (uslice(p)*uslice(p)) * stencilLaplace(p);
      }
#endif      
#if 0
      const ConstArray3d<T> dg = divgrad(pos);
      const ConstArray3d<T> blp = Stencils::Convolve(dg, stencilLaplace);
      _ArrayNdInternal::apply(a, blp, AssignScaled(dt*eps*eps));
      a(center) += 1.;
      const ConstArray3d<T> uslice(ucurr[Move(inner,pos-center)]);
       _ArrayNdInternal::apply(a[inner], uslice, dg, AddComp2(-dt));
       //PrintArray(a);
#endif
#if 1
      const ConstArray3d<T> dg = divgrad(pos);
      const ConstArray3d<T> blp = Stencils::Convolve(dg, stencilLaplace);
      boost::transform(blp, a.begin(), constant(dt*eps*eps) * _1);
      a(center) += 1.;
      const ConstArray3d<T> uslice(ucurr[Move(inner,pos-center)]);
      boost::for_each(a[inner].range(), uslice, dg, _1 += constant(-dt) * _2 * _2 * _3);
#endif
    }
  };
  StencilGen gen(ndim, eps, grid, dt, ucurr, mobility);
  
  // fill right hand side
  convolution_pos_dependent(ucurr[ucurr.indexRange()], gen.divgrad, -gen.divgrad.size()/2, qvar[qvar.indexRange()], Ops::Mult<T>(), Ops::Assign<T>());
  FOR_BBOX3(p, qvar.indexRange())
  {
    qvar(p) = ucurr(p) - dt * qvar(p);
  }

  Amr::SolveLinearSystem<T, StencilGen>(gen, bcs, qvar);
 }




void TestCahnHilliard2()
{  
  typedef double Real;
  const Bool3 cc(true, true, false);
  Array3d<Real> a = ReadArrayFromImage<Real>("chtest/structure1.png",0,1);
  a += -0.5;
  a *= 2.0;
  
  SmartStatic<Amr::Grid> pgrid(Amr::Grid(Amr::ToVertexCentered(a.getBox(), cc), 1.));
  Amr::GridVarT<Real> conc(*pgrid, 1, cc), concnext(*pgrid, 1, cc), concmid(*pgrid, 1, cc);
  conc[conc.indexRange()].fill(a);
  
  Amr::BoundaryConditions<Real> bcs;
  FOR_EACH(int i, Int2(LatticeDataQuad3d::DIR_MZ, LatticeDataQuad3d::DIR_PZ), 2)
  {
    bcs.setDiffCondition(0, i, 1, 0.);
    bcs.setDiffCondition(1, i, 3, 0.);
  }
  bcs.setPeriodicCondition(LatticeDataQuad3d::DIR_PX);
  bcs.setPeriodicCondition(LatticeDataQuad3d::DIR_PY);
  bcs.init(*pgrid, 2, cc);

  const Real dt = 20.;
  uint _t = GetTimer();
  for(int i=0; i<100; ++i)
  {
    printf("iter %i\n",i);
#if 0
    bcs.applyTo(conc);
    CahnHilliardImplicitStep2<Real>(*pgrid, conc, conc, bcs, pgrid->spacing(), 10, concnext);
    conc.fill(concnext);
#endif
#if 1
    bcs.applyTo(conc);
    CahnHilliardImplicitStep2<Real>(*pgrid, conc, conc, bcs, pgrid->spacing(), dt/2., concmid);
    bcs.applyTo(concmid);
    CahnHilliardImplicitStep2<Real>(*pgrid, conc, concmid, bcs, pgrid->spacing(), dt, concnext);
    conc.fill(concnext);
#endif
    Real totalmass = boost::accumulate(conc[conc.indexRange()].range(), 0.);
    printf("mass = %lf\n", totalmass);

    Array3d<double> out(conc[conc.indexRange()], Cons::COPY);
    out *= 0.5;
    out += 0.5;
    WriteArrayAsImage(strprintf("chtest/res_with_mobility_midpointint_%02i.png",i),out);

  }
  printf("total time: %u\n",GetTimer()-_t);
}






void TestNewConvolution()
{
  typedef double Real;
  const Bool3 cc(true, true, false);
  Array3d<Real> a = ReadArrayFromImage<Real>("chtest/structure1.png",0,1);
  a += -0.5;
  a *= 2.0;

  SmartStatic<Amr::Grid> pgrid(Amr::Grid(Amr::ToVertexCentered(a.getBox(), cc), 1.));
  Amr::GridVarT<Real> conc(*pgrid, 1, cc);
  Amr::GridVarT<Real> concnext(*pgrid, 1, cc);
  conc[conc.indexRange()].fill(a);

  Amr::GridVarT<Real> mobility(*pgrid, 1, conc.cellCentered());
  mobility.fill(1.);
  BBox3 bb = mobility.extendedIndexRange();
  bb.min[0] += 20;
  bb.max[0] -= 20;
  mobility[bb].fill(0.5);


  DivGradStencil<Real> divgrad(mobility, *pgrid, 2);
  convolution_pos_dependent(conc[conc.indexRange()], divgrad, -divgrad.size()/2, concnext[concnext.indexRange()], Ops::Mult<Real>(), Ops::Assign<Real>());

  //convolution(conc[conc.indexRange()], Stencils::Laplace(2, pgrid->spacing(), Stencils::W3_SPARSE), -divgrad.size()/2, concnext[concnext.indexRange()], Ops::Mult<Real>(), Ops::Assign<Real>());

  WriteArrayAsImage("testconvolution.png",concnext[concnext.indexRange()], 0.1, 0.5);
}








static void CahnHilliardImplicitStep(const LatticeDataQuad3d &ld, const double eps, DynArray<double> currvalue, double dt, int iter)
{
  // building the linear time stepping scheme following "unconditionally gradient stable time marching the Cahn-Hilliard equation.ps"
  // the equation to solve is: (I + k eps^2 A_h^2 - k A_h D(Un)^2) Un+1 = (I - k A_h) Un,
  // where k is the time step
  SparseMatrix<double> laplace;
  Init7PointFiniteVolumeMatrix(ld, laplace);
  AddLaplacian(ld, laplace, 1.);

  DynArray<double> rhs(currvalue.size());

  SparseMatrix<double> m;
  SparseMultiply(laplace, laplace, m);

  m *= -eps*eps;

  for(int ri=0; ri<laplace.numRows(); ++ri)
  {
    double laplace_x_currval = 0.;

    const SparseMatrix<double>::Row row = laplace.getRow(ri);
    for(int i=0; i<row.count(); ++i)
    {
      int ci = row.col(i);
      double l = row.val(i);    
      laplace_x_currval += l*currvalue[ci];
      m.Add(ri,ci, l*Sqr(currvalue[ci]));
    }

    rhs[ri] = -laplace_x_currval;
  }
  laplace.clear();

  // done building the rhs of the system (u^{n+1}-u^{n})/dt = rhs(u^{n+1},u^{n})
  // now reorder to solve for the u^{n+1}
  m *= -dt;
  rhs *= dt;
  for(int i=0; i<m.numRows(); ++i) { m.Add(i,i,1.); }
  rhs += currvalue;    

  //{
  //  BasicHistogram1D<double> h;
  //  MinMax<double> mm;
  //  for(int i=0; i<m.NumNz(); ++i) mm.Add(m.Val(i));
  //  HistogramInfo<double> hi(mm, 300, false);
  //  h.Init(hi, HistoBase::NoError);
  //  for(int i=0; i<m.NumRows(); ++i)
  //  {
  //    for(int j=m.RowBegin(i); j<m.RowEnd(i); ++j)
  //    {
  //      double val = m.Val(j);
  //      int c = m.Index(j);
  //      //if(c!=i) {
  //        h.Fill(val, 1);
  //      //}
  //    }
  //  }
  //  std::ofstream os(strprintf("mstats%02i.txt",iter).c_str());
  //  PrintHistograms1D(os, &h, 1);
  //}

  //{
  //Image img;
  //SparsityPlot(m.NumRows(),m.GetRowPtr(),m.GetIndexPtr(),m.GetDataPtr(),img);
  //img.WritePng(strprintf("m%02i.png",iter));
  //}
  //CheckMatrix(m.NumRows(),m.GetRowPtr(),m.GetIndexPtr(),m.GetDataPtr());

  //PrintMatrix(strprintf("testchmatrix%02i_ok.txt",iter).c_str(), m, rhs);

  //currvalue.fill(0.0);
  DirectSolve(m, rhs.get_ptr(), currvalue.get_ptr(), true, true);

  //Multigrid<double,AlgebraicInterpolation> mg(m, AlgebraicInterpolation());  
  //mg.numSmootherIterations = 4;
  //mg.numVCycleIterations = 4;
  //mg.FMG(currvalue, rhs);

  double r = ResidualNorm(m, rhs, currvalue);
  printf("ch step r=%f\n",r);
}


void TestCahnHilliard()
{
  Array3d<double> conc = ReadArrayFromImage<double>("chtest/noise2.png",0,1);
  LatticeDataQuad3d ld; ld.Init(conc.size()[0],conc.size()[1],1);
  conc += -0.5;
  conc *= 2.0;

  for(int i=0; i<1; ++i)
  {
    CahnHilliardImplicitStep(ld, 2., conc.dynArray(), 100.0, i);

    Array3d<double> out(conc, Cons::COPY);
    out *= 0.5;
    out += 0.5;
    WriteArrayAsImage(strprintf("chtest/res%02i.png",i),out);
  }  
}
