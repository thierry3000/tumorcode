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
#ifndef CONTINUUM_FLOW_H
#define CONTINUUM_FLOW_H

#include "common.h"
//#include "mwlib/drawing.h"
#include "continuum-utils.h"
#include "continuum-grid.h"


/*------------------------------------------------------
------------------------------------------------------*/

enum ConvectionMode {
  CONSERVATIVE,
  NON_CONSERVATIVE
};

/*------------------------------------------------------
 * These functions use the ENO or WENO scheme to compute
 * the advection operator:
 * (grad u * v) (if NON_CONSERVATIVE)
 * div (u * v) (if CONSERVATIVE)
 * u is defined on cell centers, v is defined on the edges.
 * The algorithms approximate u at the edges and then compute
 * above quantities with these edge values at the left and right
 * cell boundaries.
 * The bounding box bb corresponds to cell center indices which are
 * included. The velocity array should include the left and right
 * boundaries. The val array must include ghost points around the boundary.
 * For ENO two layers (?) and 3 layers for WENO 5th order.
------------------------------------------------------*/

template<class T, class S>
void AddENOAdvection(const BBox3 &bb, const ConstArray3d<T> &val, const ConstArray3d<S>*const velocities, int dim, T spacing, Array3d<T> res, ConvectionMode mode, int order);

template<class T, class S>
void AddWENOAdvection(const BBox3 &bb, const ConstArray3d<T> &val, const ConstArray3d<S>*const velocities, int dim, T spacing, Array3d<T> res, ConvectionMode mode, int order);

// Kurganov and Tadmor scheme. Does the same thing as ENO and WENO but differently
template<class T, class S>
void AddKTAdvection(const BBox3& bb, const ConstArray3d< T >& val, const ConstArray3d< S >*const velocities, int ndim, T spacing, Array3d< T > res, ConvectionMode mode = CONSERVATIVE);

/*------------------------------------------------------
 * This is a naive advection scheme
 * res += vel * grad values, using centered differences for the gradient, UNSTABLE with euler scheme and no diffusion
------------------------------------------------------*/
template<class T, class S>
void AddNaiveAdvection(const BBox3 &bb, const ConstArray3d<T> values, ConstArray3d<S> velocity, int axis, int dim, double spacing , Array3d<T> res, ConvectionMode mode = NON_CONSERVATIVE);

/*------------------------------------------------------
 * This is fairly old so called flux corrected transport.
 * It directly computes an update to the advected variable
 * for a fixed time step. So it is not suited for modern
 * Runge-Kutta time stepping schemes.
------------------------------------------------------*/
// this is for advection in one direction
template<class T>
void FctStep(const BBox3 &bb, Array3d<T> val, const ConstArray3d<T> u, T dt, T sourcedt, T spacing, int axis, ConstArray3d<T> *sources, ConvectionMode mode);
// 3d advection
template<class T>
void FctStep(const BBox3& bb, Array3d< T > val, const ConstArray3d< T >*const u, int ndim, T dt, T spacing, ConstArray3d< T > sources);
template<class T>
void FctStep(const BBox3& bb, Array3d< T > val, const ConstArray3d< T >*const u, int ndim, T dt, T spacing, ConvectionMode mode = CONSERVATIVE);

#if 0
template<class T>
static void AddFctStepSlope(const BBox3& bb, ConstArray3d< T > val, const ConstArray3d< T >*const u, int ndim, T dt, T spacing, Array3d<T> res, ConvectionMode mode = CONSERVATIVE)
{
  BBox3 extbb = bb;
  for (int i=0; i<ndim; ++i) {
    extbb.Move(0, i, -2).Move(1, i, 2);
  }
  Array3d<T> tmp = Array3d<T>::fromBox(extbb);
  tmp[extbb].fill(val[extbb]);
  FctStep(bb, tmp, u, ndim, dt, spacing, mode);
  tmp[bb] -= val[bb];
  tmp[bb] *= (1./dt);
  res[bb] += tmp[bb];
}
#endif

/*------------------------------------------------------
------------------------------------------------------*/

template<class T>
inline double ComputeMaximalExplicitAdvectionEulerTimeStep(const ConstArray3d<T>* const u, int dim, T spacing)
{
  double max_vel = 0;
  for (int i=0; i<dim; ++i)
    max_vel = std::max(max_vel, u[i].valueStatistics().MaxAbs());
  return 0.5 * 0.98 * spacing / (max_vel * dim + 1.e-16);
}

/*------------------------------------------------------
------------------------------------------------------*/

namespace EulerDt
{

double ByOrder(double max_ev_factor, int dim, double spacing, double order);

double ByEigenvalue(double max_ev);

inline double Convection(double max_vel, int dim, double spacing)
{
  return ByOrder(max_vel, dim, spacing, 1.);
}

inline double Diffusion(double max_kdiff, int dim, double spacing)
{
  return ByOrder(max_kdiff, dim, spacing, 2.);
}

}

/*------------------------------------------------------
------------------------------------------------------*/

template<class T>
void CopyBorder(Array3d<T> f, int dim, int border);

template<class T>
void CopyBorder(Array3d<T> f, const BBox3 &bbox, const BBox3 &full_box, int dim, int border);


/*------------------------------------------------------
------------------------------------------------------*/

Float2 LeastSquaresSolveNx2( const float* A, const float* b, const int N );
Float3 LeastSquaresSolveNx3( const float* A, const float* b, const int N );

/*------------------------------------------------------
------------------------------------------------------*/

template<class T, class S>
void CalcFlowField( const ConstArray3d<T> &pfield, const ConstArray3d<S> &cond, const LatticeDataQuad3d &ld, Array3d<Float3> &ifvel );

template<class T, class S>
void CalcDivergence( const ConstArray3d<T> &press, const ConstArray3d<S> &cond, const LatticeDataQuad3d &ld, Array3d<T> &div );

void CalcDivergence( const ConstArray3d<Float3> &ifvel, const LatticeDataQuad3d &ld, Array3d<float> &div );

/*------------------------------------------------------
------------------------------------------------------*/

/*
* res += (vcoeff grad vpot)_axis
* bbox -> cell centered box, resulting face box includes left and right borders, coefficient array should include a valid ghost boundary
* params:
*   max_vel (double) - in/out maximum value
*   geometric_coeff_mean (bool) - compute vcoeff face value by geometric mean, otherwise by normal mean
*   prefactor (double) - multiplies the velocity by prefactor before it is added to res
*   exclude_right_boundary (bool) - does not compute the velocity at the right most cell of bb
*/
template<class T, class S>
void AddVelocitiesFromPotential(const BBox3 &bb, int axis, int dim, double spacing, ConstArray3d<S> vcoeff, ConstArray3d<T> vpot, Array3d<S> &res, ptree &params);


template<class T, class S, class FaceVarFunctor>
static void AddVelocitiesFromPotential2(const BBox3 &bb, int axis, int dim, double spacing, const FaceVarFunctor vcoeff, ConstArray3d<T> vpot, Array3d<S> &res, ptree &params)
{
  LatticeIndexRanges ir = LatticeIndexRanges::FromCellRange(bb, dim);

  optional<double> max_vel = params.get_optional<double>("max_vel");
  T prefactor = params.get<double>("prefactor", 1.) / spacing;

  BBox3 bf = ir.faces[axis];
  if (params.get<bool>("exclude_right_boundary", false))
    bf.max[axis]--;

  FOR_BBOX3(p, bf)
  {
    Int3 pleft(p); --pleft[axis];
    T loc_coeff = vcoeff(p, pleft, axis, 0);
    T vel = -prefactor * loc_coeff * (vpot(p) - vpot(pleft));
    res(p) += vel;
    if (max_vel) *max_vel = std::max<double>(*max_vel, std::abs(vel));
  }

  if (max_vel)
    params.put<double>("max_vel", *max_vel);
}



// explicit diffusion term: res += div (kdiff * grad val)
template<class T, class S>
void AddLaplacian(const BBox3 &bb, ConstArray3d<T> val, ConstArray3d<S> kdiff, int dim, T spacing, Array3d<T> res, double* max_kdiff=NULL);

/*
 * // explicit diffusion term: res += prefactor * div (kdiff * grad val)
 * parameters:
 *  max_kdiff -> records the maximal diffusion coefficient (no ptree entry means no recording, otherwise the ptree value is changed)
 *  geometric_coeff_mean -> use geometric mean for the diffusion coefficient values for the cells faces
 *  prefactor -> default is 1
 */
template<class T, class S>
void AddLaplacian(const BBox3 &bb, ConstArray3d<T> val, ConstArray3d<S> kdiff, int dim, T spacing, Array3d<T> res_arr, ptree &params);


template<class T, class S, class FaceVarFunctor>
static void AddLaplacian2(const BBox3 &bb, ConstArray3d<T> val, const FaceVarFunctor kdiff, int dim, T spacing, Array3d<T> res_arr, ptree &params)
{
  optional<double> max_kdiff = params.get_optional<double>("max_kdiff");
  T prefactor = params.get<double>("prefactor", 1.) / spacing;
  bool isotropic = params.get<string>("stencil", "")=="isotropic";
  
  Convolutor<T, S> diff[3];
  for (int a=0; a<dim; ++a) {
    diff[a].Init(Stencils::Diff1(a, dim, spacing, isotropic ? Stencils::W3_ISOTROPIC : Stencils::W1, Stencils::FWD), val, true);
    //diff[a].Init(Stencils::Diff1(a, dim, spacing, Stencils::W1, Stencils::FWD), val);
  }

  for (int axis=0; axis<dim; ++axis)
  {
    BBox3 bbext(bb);
    bbext.min[axis]--;
    
    FOR_BBOX3(p, bbext)
    {
      Int3 q(p); ++q[axis];
      T grad = diff[axis].point_convolve(p); // gradient
      T kd   = kdiff(p, q, axis, 1);
      T flux = -prefactor * kd * grad;
      if (p[axis]>=bb.min[axis])
        res_arr(p) -= flux;
      if (q[axis]<=bb.max[axis])
        res_arr(q) += flux;
      if (max_kdiff) *max_kdiff = std::max<double>(*max_kdiff, kd);
    }
  }

  if (max_kdiff)
    params.put<double>("max_kdiff", *max_kdiff);
}


#endif
