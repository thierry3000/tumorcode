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

#include "continuum-flow.h"
#include "mwlib/drawing.h"
#include <stdio.h>

// #include "calcflow_common.h"
// #include "calcflow_linsys.h"

using namespace my;

/*------------------------------------------------------
 * copies a boundary shell worth of one site from boundary offset 'level'-1 to offset level 'level'
------------------------------------------------------*/
template<class T>
static void CopyBorder_(Array3d<T> f, const BBox3 &box_, const BBox3& full_box, int ndims, int level)
{
  BBox3 box(box_);
  bool at_edge[3][2];
  // extend the box where it hits the boundaries of the full_box
  // have to extend the box first in all directions in order to capture boundary points
  for(int dim=0; dim<ndims; ++dim)
  {
    at_edge[dim][0] = box.min[dim] <= full_box.min[dim];
    at_edge[dim][1] = box.max[dim] >= full_box.max[dim];
    if (at_edge[dim][0]) box.min[dim] -= level;
    if (at_edge[dim][1]) box.max[dim] += level;
  }
  // set up slices and copy
  for (int dim=0; dim<ndims; ++dim)
  {
    for (int side=0; side<2; ++side)
    {
      if (!at_edge[dim][side]) continue;
      BBox3 bdst(box);
      bdst.Set(dim, box[side][dim]);
      BBox3 bsrc(bdst);
      bsrc.Set(dim, bdst[side][dim]+(side==0 ? 1 : -1));
      f[bdst].fill(f[bsrc]);
    }
  }
//     Int3 bxyz(bx, by, bz);
//     const BBox3 fbox = f.getBox();
//     BBox3 bbext(fbox); bbext.Extend(bxyz);
//
//     for(int dim=0; dim<3; ++dim)
//     {
//       if(!bxyz[dim]) continue;
//       for(int side=0; side<2; ++side)
//       {
//         BBox3 bb_src = bbext, bb_dst = bbext;
//         bb_src.Set(dim, bb_src[side][dim] + (side==0 ? 1 : -1));
//         bb_dst.Set(dim, bb_dst[side][dim]);
//         if (!clip_box.IsEmpty())
//         {
//           bb_src.Clip(clip_box);
//           bb_dst.Clip(clip_box);
//           if (bb_dst.IsEmpty()) continue;
//         }
//         Array3d<T> dst(f.slice(bb_dst));
//         Array3d<T> src(f.slice(bb_src));
//         dst.fill(src);
//       }
//     }
}


template<class T>
void CopyBorder(Array3d<T> f, const BBox3 &bbox, const BBox3 &full_box, int dim, int border)
{
  myAssert(border>=0);
  for (int i=0; i<border; ++i)
  {
    CopyBorder_<T>(f, bbox, full_box, dim, i+1);
  }
}


template<class T>
void CopyBorder(Array3d<T> f, int dim, int border)
{
  myAssert(border>=0);
  const BBox3 bbox = f.getBox();
  for (int i=0; i<border; ++i)
  {
    CopyBorder_<T>(f, bbox, bbox, dim, i+1);
  }
}

#define INSTANTIATE_CopyBorder(T) \
  template void CopyBorder<T>(Array3d<T> f, const BBox3 &bbox, const BBox3 &full_box, int dim, int border); \
  template void CopyBorder<T>(Array3d<T> f, int dim, int border); 

INSTANTIATE_CopyBorder(double)
INSTANTIATE_CopyBorder(float)
INSTANTIATE_CopyBorder(int)

  
/*------------------------------------------------------
------------------------------------------------------*/

template<class T>
void FctStep(const BBox3& bb, Array3d< T > val, const ConstArray3d< T > u, T dt, T sourcedt, T spacing, int axis, ConstArray3d< T > *sources, ConvectionMode mode)
{
  const BBox3 outer_bb = BBox3(bb).Move(0, axis, -2).Move(1, axis, 2); // 2 site border in axis dir
  Array3d<T> res(outer_bb);
  // normal convection step
  FOR_BBOX3(p, bb)
  {
    Int3 pa(p),pb(p); pa[axis]--; pb[axis]++;
    const T q_a = val(pa);
    const T q = val(p);
    const T q_b = val(pb);
    const T e_a = u(p) * dt/spacing;
    const T e_b = u(pb) * dt/spacing;
    T nu_a = 1./6. + 1./3.*e_a*e_a;
    T nu_b = 1./6. + 1./3.*e_b*e_b;
    //T mu_a = 1./6. - 1./6.*e_a*e_a;
    //T mu_b = 1./6. - 1./6.*e_b*e_b;
    res(p) = q - 0.5 * (e_b * (q_b + q) - e_a * (q + q_a)) + nu_b * (q_b - q) - nu_a * (q - q_a);
    if (mode == NON_CONSERVATIVE) // experimental
    {
      // add dt*val*div(v) = val_i*dt/h * (u_i+1 - u_i) + other dimension which are treated by operator splitting
      res(p) += q*(e_b - e_a); 
    }
  }
  if (sources)
  {
    res[bb].addScaled(sourcedt, (*sources)[bb]);
  }
  FOR_BBOX3(p, bb)
  {
    Int3 pp[5] = { p, p, p, p, p };
    pp[0][axis]-=2;
    pp[1][axis]--;
    pp[3][axis]++;
    pp[4][axis]+=2;
    T r[5];
    for(int i=0; i<5; ++i) r[i] = res(pp[i]);
    T fc[2];
    for(int i=1; i<=2; ++i)
    {
      T s = my::sign(r[i+1]-r[i]);
      const T eps = u(pp[i+1]) * dt/spacing;
      T mu = 1./6. + 1./3.*eps*eps;
      T t1 = std::abs(mu * (r[i+1]-r[i]));
      T t0 = s * (r[i+2]-r[i+1]);
      T t2 = s * (r[i]-r[i-1]);
      fc[i-1] = s * std::max<T>(0., std::min(t0, std::min(t1, t2)));
    }
    val(p) = res(p) - fc[1] + fc[0];
  }
}


template<class T>
void FctStep(const BBox3 &bb, Array3d<T> val, const ConstArray3d<T>* const u, int ndim, T dt, T spacing, ConstArray3d<T> sources)
{
  for(int axis=0; axis<ndim; ++axis)
  {
    FctStep<T>(bb, val, u[axis], dt, dt/ndim, spacing, axis, &sources, CONSERVATIVE);
  }
}

template<class T>
void FctStep(const BBox3 &bb, Array3d<T> val, const ConstArray3d<T>* const u, int ndim, T dt, T spacing, ConvectionMode mode)
{
  for(int axis=0; axis<ndim; ++axis)
  {
    FctStep<T>(bb, val, u[axis], dt, dt/ndim, spacing, axis, NULL, mode);
  }
}


#define INSTANTIATE_FCT(TYPE)\
  template void FctStep<TYPE>(const BBox3 &bb, Array3d<TYPE> val, const ConstArray3d<TYPE> u, TYPE dt, TYPE sourcedt, TYPE spacing, int axis, ConstArray3d<TYPE> *sources, ConvectionMode mode); \
  template void FctStep<TYPE>(const BBox3 &bb, Array3d<TYPE> val, const ConstArray3d<TYPE>* const u, int ndim, TYPE dt, TYPE spacing, ConstArray3d<TYPE> sources); \
  template void FctStep<TYPE>(const BBox3 &bb, Array3d<TYPE> val, const ConstArray3d<TYPE>* const u, int ndim, TYPE dt, TYPE spacing, ConvectionMode mode); \

INSTANTIATE_FCT(float)
INSTANTIATE_FCT(double)




/*------------------------------------------------------
------------------------------------------------------*/

template<class T>
T inline minmod(T a, T b, T c)
{
  T abcmin = std::min(a, std::min(b, c));
  T abcmax = std::max(a, std::max(b, c));
  if (abcmin > 0.)
    return abcmin;
  else if (abcmax < 0.)
    return abcmax;
  else
    return 0;
}


template<class T>
T inline minmodgrad(T *v, T one_over_dx)
{
  // van leer MC-limiter
  T theta = 2.;
  return one_over_dx * minmod<T>(
    theta*(v[1]-v[0]),
    0.5*(v[2]-v[0]),
    theta*(v[2]-v[1])
  );
}


// gugst du
// kurganov_tadmor_new high res central schemes for nonlinear conservation laws and convection diffusion eq

template<class T, class S>
void AddKTAdvection(const BBox3& bb, const ConstArray3d< T >& val, const ConstArray3d< S >*const velocities, int ndim, T spacing, Array3d< T > res, ConvectionMode mode)
{
  FOR_BBOX3(p, bb)
  {
    T fluxes = 0.;
    for (int axis=0; axis<ndim; ++axis)
    {
      Int3 px = p;
      T v[5];
      T u[2];
      // copy stuff to local variables
      for (int i=0; i<2; ++i) {
        px[axis] = p[axis]+i;
        u[i] = velocities[axis](px);
      }
      for (int i=0; i<5; ++i) {
        px[axis] = p[axis]-2+i;
        v[i] = val(px);
      }
      // derivatives
      T vx[5];
      for (int i=0; i<3; ++i) {
        vx[i+1] = minmodgrad<T>(v+i, 1./spacing);
      }
      // compute fluxes
      for (int side=0; side<2; ++side)
      {
        T v_minus = v[side+1] + spacing*0.5* vx[side+1];
        T v_plus = v[side+2] - spacing*0.5* vx[side+2];
        T flux = 0.5*((v_minus+v_plus)*u[side] - std::abs(u[side])*(v_plus-v_minus));
        fluxes += (side==0) ? -flux : flux;
        if (mode == NON_CONSERVATIVE)
        {
          fluxes -= ((side==0) ? -u[side] : u[side]) * v[2];
        }
      }
    }
    res(p) -= fluxes / spacing;
  }
}





// gugst du
// sussman_an improved level set method for incompressible two phase flows
// shu_osher_1989_efficient implementation of ENO
// an assessment of semi-discrete central schemes for hyper bolic conservation laws

template<class S, class T>
void ENOReconstructFaceValues(const LatticeIndexRanges &ir, int dim, double spacing, int axis, const ConstArray3d<S> &values, const ConstArray3d<T> &velocities, Array3d<S> facevalues, int order)
{
  S factor = 1./spacing;
  S factor2 = 1./my::sqr(spacing);

  //LatticeIndexRanges ir = LatticeIndexRangesFromCellRange(bb, dim);
  //for (int axis=0; axis<dim; ++axis)
  //{
    FOR_BBOX3(p, ir.faces[axis])
    {
      Int3 q(p);
      if (velocities(p) > 0) q[axis]--; // upwind
      int k1 = q[axis];
      S fk1 = values(q);
      S f1 = fk1;
      if (order > 1)
      {
        S fk1_m_1 = values(add_to_axis(q, axis,-1));
        S fk1_p_1 = values(add_to_axis(q, axis,1));
        S a = (fk1-fk1_m_1)*factor;
        S b = (fk1_p_1 - fk1)*factor;
        bool t = std::abs(a)<std::abs(b);
        S c = t ? a : b;
        int k2 = t ? k1-1 : k1;
        q[axis] = k2;
        S f2 = f1 + spacing*0.5*c*(1-2*(k1-p[axis]+1));
        if (order > 2)
        {
          S fk2 = values(q);
          S fk2_m_1 = values(add_to_axis(q, axis,-1));
          S fk2_p_1 = values(add_to_axis(q, axis, 1));
          S fk2_p_2 = values(add_to_axis(q, axis, 2));
          a = (fk2_m_1 - 2.*fk2 + fk2_p_1)*factor2;
          b = (fk2 - 2.*fk2_p_1 + fk2_p_2)*factor2;
          c = std::abs(a)<std::abs(b) ? a : b;
          S f3 = f2 + my::sqr(spacing)/6. * c * (3*my::sqr(k2-p[axis]+1)-1);
          facevalues(p) = f3;
        }
        else facevalues(p) = f2;
      }
      else facevalues(p) = f1;
    }
  //}
}


template<class T, class S>
void AddENOAdvection(const BBox3 &bb, const ConstArray3d<T> &val, const ConstArray3d<S>*const velocities, int dim, T spacing, Array3d<T> res, ConvectionMode mode, int order)
{
  LatticeIndexRanges ir = LatticeIndexRanges::FromCellRange(bb, dim);
  
  T factor = 1./spacing;
  for (int axis=0; axis<dim; ++axis)
  {
    Array3d<T> facevalues(ir.faces[axis]);
    ENOReconstructFaceValues<T, S>(ir, dim, spacing, axis, val, velocities[axis], facevalues, order);
    FOR_BBOX3(p, ir.cells)
    {
      T r_p = 0.;
      Int3 q = add_to_axis(p, axis, 1);
      S u1 = velocities[axis](q);
      S u0 = velocities[axis](p);
      S phi0 = facevalues(p);
      S phi1 = facevalues(q);
      if (mode == CONSERVATIVE)
      {
        r_p = factor*(u1*phi1 - u0*phi0);
      }
      else
      {
        r_p = 0.5*(u0+u1)*(phi1-phi0)*factor;
      }
      res(p) -= r_p;
    }
  }
}


template<class S, class T>
void WENOReconstructFaceValues(const LatticeIndexRanges &ir, int dim, double spacing, int axis, const ConstArray3d<S> &values, Array3d<S> face_minus, Array3d<S> face_plus, int order)
{
  S sigma0, sigma1, sigma2, omega0, omega1, omega2, acc, omega3, omega4, omega5, fr0, fr1, fr2, fr3, fr4, fr5, fs0, fs1;
  const int fsi = 1;

  BBox3 bb = ir.cells;
  bb.max[axis]++;
  bb.min[axis]--;
  FOR_BBOX3(p, bb)
  {
    if (order == 5)
    {
      S f[5];
      f[0] = values(add_to_axis(p, axis, -2));
      f[1] = values(add_to_axis(p, axis, -1));
      f[2] = values(p);
      f[3] = values(add_to_axis(p, axis, 1));
      f[4] = values(add_to_axis(p, axis, 2));
      const int i = 2;
      if (std::abs(f[0])>200 || std::abs(f[1])>200 || std::abs(f[2])>200 || std::abs(f[3])>200 || std::abs(f[4])>200)
        fprintf(stderr, "%f %f %f %f %f", f[0], f[1], f[2], f[3], f[4]);
      // generated with pyweno!
      //--- smoothness ---
      sigma0 = 3.3333333333333333333333333333333333*f[(i+0)*fsi]*f[(i+0)*fsi] - 10.333333333333333333333333333333333*f[(i+0)*fsi]*f[(i+1)*fsi] + 3.6666666666666666666666666666666667*f[(i+0)*fsi]*f[(i+2)*fsi] + 8.3333333333333333333333333333333333*f[(i+1)*fsi]*f[(i+1)*fsi] - 6.3333333333333333333333333333333333*f[(i+1)*fsi]*f[(i+2)*fsi] + 1.3333333333333333333333333333333333*f[(i+2)*fsi]*f[(i+2)*fsi];
      sigma1 = 4.3333333333333333333333333333333333*f[(i+0)*fsi]*f[(i+0)*fsi] - 4.3333333333333333333333333333333333*f[(i+0)*fsi]*f[(i+1)*fsi] - 4.3333333333333333333333333333333333*f[(i+0)*fsi]*f[(i-1)*fsi] + 1.3333333333333333333333333333333333*f[(i+1)*fsi]*f[(i+1)*fsi] + 1.6666666666666666666666666666666667*f[(i+1)*fsi]*f[(i-1)*fsi] + 1.3333333333333333333333333333333333*f[(i-1)*fsi]*f[(i-1)*fsi];
      sigma2 = 3.3333333333333333333333333333333333*f[(i+0)*fsi]*f[(i+0)*fsi] - 10.333333333333333333333333333333333*f[(i+0)*fsi]*f[(i-1)*fsi] + 3.6666666666666666666666666666666667*f[(i+0)*fsi]*f[(i-2)*fsi] + 8.3333333333333333333333333333333333*f[(i-1)*fsi]*f[(i-1)*fsi] - 6.3333333333333333333333333333333333*f[(i-1)*fsi]*f[(i-2)*fsi] + 1.3333333333333333333333333333333333*f[(i-2)*fsi]*f[(i-2)*fsi];
      //--- weights ---
      omega0 = 0.1/(1.0e-6*1.0e-6 + 2*1.0e-6*sigma0 + sigma0*sigma0);
      omega1 = 0.6/(1.0e-6*1.0e-6 + 2*1.0e-6*sigma1 + sigma1*sigma1);
      omega2 = 0.3/(1.0e-6*1.0e-6 + 2*1.0e-6*sigma2 + sigma2*sigma2);
      acc = omega0 + omega1 + omega2;
      omega0 = omega0/acc;
      omega1 = omega1/acc;
      omega2 = omega2/acc;
      omega3 = 0.29999999999999971134201359745929949/(1.0e-6*1.0e-6 + 2*1.0e-6*sigma0 + sigma0*sigma0);
      omega4 = 0.59999999999999986677323704498121515/(1.0e-6*1.0e-6 + 2*1.0e-6*sigma1 + sigma1*sigma1);
      omega5 = 0.1000000000000000055511151231257827/(1.0e-6*1.0e-6 + 2*1.0e-6*sigma2 + sigma2*sigma2);
      acc = omega3 + omega4 + omega5;
      omega3 = omega3/acc;
      omega4 = omega4/acc;
      omega5 = omega5/acc;
      //--- reconstructions ---
      fr0 = 1.8333333333333333333333333333333333*f[(i+0)*fsi] - 1.1666666666666666666666666666666667*f[(i+1)*fsi] + 0.33333333333333333333333333333333333*f[(i+2)*fsi];
      fr1 = 0.83333333333333333333333333333333333*f[(i+0)*fsi] - 0.16666666666666666666666666666666667*f[(i+1)*fsi] + 0.33333333333333333333333333333333333*f[(i-1)*fsi];
      fr2 = 0.33333333333333333333333333333333333*f[(i+0)*fsi] + 0.83333333333333333333333333333333333*f[(i-1)*fsi] - 0.16666666666666666666666666666666667*f[(i-2)*fsi];
      fr3 = 0.33333333333333325931846502498956397*f[(i+0)*fsi] + 0.83333333333333303727386009995825589*f[(i+1)*fsi] - 0.16666666666666674068153497501043603*f[(i+2)*fsi];
      fr4 = 0.83333333333333348136306995002087206*f[(i+0)*fsi] + 0.33333333333333337034076748750521801*f[(i+1)*fsi] - 0.16666666666666674068153497501043603*f[(i-1)*fsi];
      fr5 = 1.833333333333333259318465024989564*f[(i+0)*fsi] - 1.1666666666666665186369300499791279*f[(i-1)*fsi] + 0.33333333333333337034076748750521801*f[(i-2)*fsi];
      fs0 = fr0*omega0 + fr1*omega1 + fr2*omega2;
      fs1 = fr3*omega3 + fr4*omega4 + fr5*omega5;
      // assign values
      /*
      * ---------------------------------------------------------------------
      * |         p-1        |         p          |            p+1          |
      * | p-1^+          p^- | p^+          p+1^- | p+1^+                   |
      */
    }
    else
    {
      S f[3];
      f[0] = values(add_to_axis(p, axis, -1));
      f[1] = values(p);
      f[2] = values(add_to_axis(p, axis, 1));
      const int i = 1;
      //--- smoothness ---
      sigma0 = f[(i+0)*fsi]*f[(i+0)*fsi] - 2.0*f[(i+0)*fsi]*f[(i+1)*fsi] + f[(i+1)*fsi]*f[(i+1)*fsi];
      sigma1 = f[(i+0)*fsi]*f[(i+0)*fsi] - 2.0*f[(i+0)*fsi]*f[(i-1)*fsi] + f[(i-1)*fsi]*f[(i-1)*fsi];
      //--- weights ---
      omega0 = 0.33333333333333333333333333333333333/(1.0e-6*1.0e-6 + 2*1.0e-6*sigma0 + sigma0*sigma0);
      omega1 = 0.66666666666666666666666666666666667/(1.0e-6*1.0e-6 + 2*1.0e-6*sigma1 + sigma1*sigma1);
      acc = omega0 + omega1;
      omega0 = omega0/acc;
      omega1 = omega1/acc;
      omega2 = 0.66666666666666666666666666666666667/(1.0e-6*1.0e-6 + 2*1.0e-6*sigma0 + sigma0*sigma0);
      omega3 = 0.33333333333333333333333333333333333/(1.0e-6*1.0e-6 + 2*1.0e-6*sigma1 + sigma1*sigma1);
      acc = omega2 + omega3;
      omega2 = omega2/acc;
      omega3 = omega3/acc;
      //--- reconstructions ---
      fr0 = 1.5*f[(i+0)*fsi] - 0.5*f[(i+1)*fsi];
      fr1 = 0.5*f[(i+0)*fsi] + 0.5*f[(i-1)*fsi];
      fr2 = 0.5*f[(i+0)*fsi] + 0.5*f[(i+1)*fsi];
      fr3 = 1.5*f[(i+0)*fsi] - 0.5*f[(i-1)*fsi];
      fs0 = fr0*omega0 + fr1*omega1;
      fs1 = fr2*omega2 + fr3*omega3;
    }
    if (p[axis]>bb.min[axis])
      face_plus(p) = fs0;
    if (p[axis]<bb.max[axis])
      face_minus(add_to_axis(p,axis,1)) = fs1;
  }
}



template<class T, class S>
void AddWENOAdvection(const BBox3 &bb, const ConstArray3d<T> &val, const ConstArray3d<S>*const velocities, int dim, T spacing, Array3d<T> res, ConvectionMode mode, int order)
{
  LatticeIndexRanges ir = LatticeIndexRanges::FromCellRange(bb, dim);

  T factor = 1./spacing;
  for (int axis=0; axis<dim; ++axis)
  {
    Array3d<T> face_m(ir.faces[axis]), face_p(ir.faces[axis]);
    
    WENOReconstructFaceValues<T, S>(ir, dim, spacing, axis, val, face_m, face_p, order);
    FOR_BBOX3(p, ir.cells)
    {
      T r_p = 0.;
      Int3 q = add_to_axis(p, axis, 1);
      S u1 = velocities[axis](q);
      S u0 = velocities[axis](p);
      S phi0 = u0 > 0 ? face_m(p) : face_p(p); // upwind flux
      S phi1 = u1 > 0 ? face_m(q) : face_p(q);
      if (mode == CONSERVATIVE)
      {
        r_p = factor*(u1*phi1 - u0*phi0);
      }
      else
      {
        r_p = 0.5*(u0+u1)*(phi1-phi0)*factor;
      }
      res(p) -= r_p;
    }
  }
}



#define INSTANTIATE_ADVECTION_SCHEMES(T, S)\
  template void AddENOAdvection(const BBox3 &bb, const ConstArray3d<T> &val, const ConstArray3d<S>*const velocities, int dim, T spacing, Array3d<T> res, ConvectionMode mode, int order); \
  template void AddWENOAdvection(const BBox3 &bb, const ConstArray3d<T> &val, const ConstArray3d<S>*const velocities, int dim, T spacing, Array3d<T> res, ConvectionMode mode, int order); \
  template void AddNaiveAdvection<T,S>(const BBox3 &bb, const ConstArray3d<T> values, ConstArray3d<S> velocity, int axis, int dim, double spacing , Array3d<T> res, ConvectionMode mode); \
  template void AddKTAdvection<T,S>(const BBox3 &bb, const ConstArray3d<T> &val, const ConstArray3d<S>*const u, int ndim, T spacing, Array3d<T> res, ConvectionMode mode); \

INSTANTIATE_ADVECTION_SCHEMES(double, double)
INSTANTIATE_ADVECTION_SCHEMES(double, float)
INSTANTIATE_ADVECTION_SCHEMES(float, float)

/*------------------------------------------------------
------------------------------------------------------*/

template<class T, class S>
void CalcFlowField( const ConstArray3d<T> &pfield, const ConstArray3d<S> &cond, const LatticeDataQuad3d &ld, Array3d<Float3> &ifvel )
{
  typedef LatticeDataQuad3d LD;
  myAssert(pfield.size()==ld.Size());
  myAssert(cond.size()==ld.Size());
  if(ifvel.empty()) ifvel.init(ld.Size());
  myAssert(ifvel.size()==ld.Size());

  Array3d<float> diff[3];
  for(int i=0; i<3; ++i)
  {
    BuildDiff1FilterQuad3d(diff[i],i);
  }

  FOR_REG3V2(p,pfield.size())
  {
    const S k0 = cond(p);
    Float3 v;
    for(int i=0; i<3; ++i) v[i] = Convolve3x3x3( diff[i], pfield, p );
    ifvel(p) = (k0/ld.Scale())*v;
  }
}


template<class T, class S>
void CalcDivergence( const ConstArray3d<T> &press, const ConstArray3d<S> &cond, const LatticeDataQuad3d &ld, Array3d<T> &div )
{
  typedef LatticeDataQuad3d LD;
  myAssert(press.size()==ld.Size());
  myAssert(cond.size()==ld.Size());
  if(div.empty()) div.init(ld.Size());
  myAssert(div.size()==ld.Size());

  T factor = (1.)/(ld.Scale()*ld.Scale());
  const BBox3 bb = press.getBox();
  FOR_BBOX3(p,bb)
  {
    T p0 = press(p);
    S k0 = cond(p);
    // total flow through a face of the voronoi cell:  unit_cell_length*ld.scale*velocity
    // total flow per area: sum over flows per face / (unit_cell_area*ld.scale*ld.scale)
    T srcval = 0;
    for(int k=0; k<LD::DIR_CNT; ++k)
    {
      const Int3 nbp = ld.NbLattice(p,k);
      S k1 = cond(nbp);
      S kk = CalcInterfaceCondCoeff(k0,k1);
      T p1 = press(nbp);
      T vk = factor*kk*(p0-p1);
      srcval += vk;
    }
    div(p) = srcval;
  }
}


void CalcDivergence( const ConstArray3d<Float3> &ifvel, const LatticeDataQuad3d &ld, Array3d<float> &div )
{
  typedef LatticeDataQuad3d LD;
  myAssert(ifvel.size()==ld.Size());
  if(div.empty()) div.init(ld.Size());
  myAssert(div.size()==ld.Size());

  Float3 wdelta[LD::DIR_CNT];
  GetWorldDirectionAxes(ld, wdelta);

  PdeReal factor = (1.)/(ld.Scale()*ld.Scale());
  const BBox3 bb = ifvel.getBox();
  FOR_BBOX3(p,bb)
  {
    Float3 vel = ifvel(p);
    PdeReal srcval = 0;

    for(int k=0; k<LD::DIR_CNT; ++k)
    {
      const Int3 nbp = ld.NbLattice(p,k);
      PdeReal vk = factor*0.5f*(ifvel(nbp)+vel).dot(wdelta[k]);
      srcval += vk;
    }
    div(p) = srcval;
  }
}



template<class T, class S>
void AddVelocitiesFromPotential(const BBox3 &bb, int axis, int dim, double spacing, ConstArray3d<S> vcoeff, ConstArray3d<T> vpot, Array3d<S> &res, ptree &params)
{
  LatticeIndexRanges ir = LatticeIndexRanges::FromCellRange(bb, dim);

  optional<double> max_vel = params.get_optional<double>("max_vel");
  bool geometric_coeff_mean = params.get<bool>("geometric_coeff_mean", true);
  T prefactor = params.get<double>("prefactor", 1.) / spacing;
  bool isotropic = params.get<string>("stencil", "")=="isotropic";

  BBox3 bf = ir.faces[axis];
  if (params.get<bool>("exclude_right_boundary", false))
    bf.max[axis]--;

  Convolutor<double, T> diff;
  diff.Init(Stencils::Diff1(axis, dim, 1., isotropic ? Stencils::W3_ISOTROPIC : Stencils::W1, Stencils::FWD), vpot, true);
  
  FOR_BBOX3(p, bf)
  {
    Int3 pleft(p); --pleft[axis];
    T cl = vcoeff(pleft);
    T cr = vcoeff(p);
    T loc_coeff = geometric_coeff_mean ? CalcInterfaceCondCoeff(cl, cr) : (0.5*(cl+cr));
#if 0
    T vel = -prefactor * loc_coeff * (vpot(p) - vpot(pleft));
#else
    T vel = -prefactor * loc_coeff * diff.point_convolve(pleft);
#endif
    res(p) += vel;
    if (max_vel) *max_vel = std::max<double>(*max_vel, std::abs(vel));
  }

  if (max_vel)
    params.put<double>("max_vel", *max_vel);
}


template<class T, class S>
void AddNaiveAdvection(const BBox3 &bb, const ConstArray3d<T> values, ConstArray3d<S> velocity, int axis, int dim, double spacing , Array3d<T> res, ConvectionMode mode)
{
  myAssert(mode == NON_CONSERVATIVE);
  FOR_BBOX3(p, bb)
  {
    Int3 face_l(p), face_r(p); face_r[axis]++;
    Int3 cell_l(p), cell_r(p); cell_l[axis]--; cell_r[axis]++;
    float ls_grad = 0.5 * (values(cell_r) - values(cell_l)) / spacing;
    float loc_vel = 0.5 * (velocity(face_l) + velocity(face_r));
    res(p) -= ls_grad * loc_vel;
  }
}




#define INSTANTIATE_FLOW_STUFF(T, S) \
  template void CalcFlowField<T,S>( const ConstArray3d<T> &pfield, const ConstArray3d<S> &cond, const LatticeDataQuad3d &ld, Array3d<Float3> &ifvel ); \
  template void CalcDivergence<T,S>( const ConstArray3d<T> &press, const ConstArray3d<S> &cond, const LatticeDataQuad3d &ld, Array3d<T> &div ); \
  template void AddVelocitiesFromPotential<T,S>(const BBox3 &bb, int axis, int dim, double spacing, ConstArray3d<S> vcoeff, ConstArray3d<T> vpot, Array3d<S> &res, ptree &params);

INSTANTIATE_FLOW_STUFF(float, float)
INSTANTIATE_FLOW_STUFF(double, float)
INSTANTIATE_FLOW_STUFF(double, double)


template<class T, class S>
void AddLaplacian(const BBox3 &bb, ConstArray3d<T> val, ConstArray3d<S> kdiff, int dim, T spacing, Array3d<T> res_arr, ptree &params)
{
  optional<double> max_kdiff = params.get_optional<double>("max_kdiff");
  bool geometric_coeff_mean = params.get<bool>("geometric_coeff_mean", true);
  T prefactor = params.get<double>("prefactor", 1.) / spacing;
  
  Convolutor<T, T> diff[3];
  for (int a=0; a<dim; ++a) {
    diff[a].Init(Stencils::Diff1(a, dim, spacing, Stencils::W1, Stencils::FWD), val);
  }
  FOR_BBOX3(p, bb)
  {
    T res = 0.;
    S d1 = kdiff(p);
    if (max_kdiff && d1>*max_kdiff) *max_kdiff = d1;
    for (int a=0; a<dim; ++a)
    {
      for (int side=0; side<2; ++side)
      {
        Int3 pp(p);
        if (side==0) --pp[a];
        T gnc = diff[a].point_convolve(pp); // gradient
        if (side==1) ++pp[a];
        S d2 = kdiff(pp); // diffusion const
        S dd = geometric_coeff_mean ? CalcInterfaceCondCoeff(d1, d2) : (0.5*(d1+d2)); // averaged face diffusion const
        res += (side==0 ? -1. : 1.) * gnc * dd;
      }
    }
    res_arr(p) += prefactor * res;
  }

  if (max_kdiff)
    params.put<double>("max_kdiff", *max_kdiff);
}


template<class T, class S>
void AddLaplacian(const BBox3 &bb, ConstArray3d<T> val, ConstArray3d<S> kdiff, int dim, T spacing, Array3d<T> res_arr, double *max_kdiff)
{
  ptree pt;
  if (max_kdiff) pt.put<double>("max_kdiff", *max_kdiff);
  AddLaplacian(bb, val, kdiff, dim, spacing, res_arr, pt);
  if (max_kdiff) *max_kdiff = pt.get<double>("max_kdiff");
}



#define INSTANTIATE_DIFFUSION(T, S) \
  template void AddLaplacian<T,S>(const BBox3 &bb, ConstArray3d<T> val, ConstArray3d<S> kdiff, int dim, T spacing, Array3d<T> res_arr, ptree &params);\
  template void AddLaplacian<T,S>(const BBox3 &bb, ConstArray3d<T> val, ConstArray3d<S> kdiff, int dim, T spacing, Array3d<T> res, double *max_kdiff);

INSTANTIATE_DIFFUSION(float, float)
INSTANTIATE_DIFFUSION(double, float)

/*------------------------------------------------------
------------------------------------------------------*/

namespace EulerDt
{

double ByEigenvalue(double max_ev)
{
  myAssert(max_ev>=0.);
  if (max_ev>0.)
    return 0.5/max_ev;
  else
    return std::numeric_limits<double>::max();
}
  
double ByOrder(double max_ev_factor, int dim, double spacing, double order)
{
  myAssert(max_ev_factor>=0. && dim>0 && spacing>0.);
  double denom = max_ev_factor * dim;
  if (denom <= 0.) return std::numeric_limits<double>::max();
  double nom = 0.5 * std::pow(spacing, order);
  return  nom/denom;
}

}
