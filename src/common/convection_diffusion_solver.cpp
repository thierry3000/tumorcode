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
#include "convection_diffusion_solver.h"





/*------------------------------------------------------
------------------------------------------------------*/

template<class T_>
void ConvectionReactionDiffusionModel<T_>::calcSlope( const State &state,
                                                  State &slope,
                                                  Steppers::StepControl &ctrl)
{
  State tmp;
  calcSlopesIMEX(state, slope, tmp, ctrl);
  slope += tmp;

  ctrl.euler_dt = my::min<double>(euler_dt_kdiff, euler_dt_srcimpl, ctrl.euler_dt);
}




template<class T_>
void ConvectionReactionDiffusionModel<T_>::calcSlopesIMEX(  const State &state,
                                                        State &slope_expl,
                                                        State &slope_impl,
                                                        Steppers::StepControl &ctrl)
{
  if (slope_expl.empty())
    slope_expl.initFromBox(grid->ld.Box());
  else
    slope_expl.fill(0);
  if (slope_impl.empty())
    slope_impl.initFromBox(grid->ld.Box());
  else
    slope_impl.fill(0);

  CopyBorder(const_cast<State&>(state)[grid->ld.Box()], grid->dim, 3);
  
  #pragma omp parallel
  {
    //VirtualGridFunctions<2>::List arr2;
    //VirtualGridFunctions<3>::List arr3;
    Array3d<T> src_val, src_clin;
    Array3d<float> kdiff;
    T thread_max_src_impl=0., thread_max_src_expl=0., thread_max_kdiff=0., thread_max_vel=0.;

    #pragma omp for nowait schedule(dynamic, 1)
    for (int boxindex=0; boxindex<mtboxes->size(); ++boxindex)
    {
      const BBox3 &bbox = (*mtboxes)[boxindex];
      // init local vars
      src_val.initFromBox(bbox);
      src_clin.initFromBox(bbox);
      // explicit source terms
      if (source_explicit_f)
      {
        source_explicit_f(boxindex, bbox, state, src_val, thread_max_src_expl);
        slope_expl[bbox] += src_val;
      }
      if (source_implicit_f)
      {
        source_implicit_f(boxindex, bbox, state, src_clin, src_val);
        FOR_BBOX3(p, bbox)
        {
          slope_impl(p) += src_clin(p) * state(p) + src_val(p);
          thread_max_src_impl = std::max(thread_max_src_impl, std::abs(src_clin(p)));
        }
      }
      if (diffusivity_f)
      {
        const BBox3 bbext = ExtendForDim(bbox, grid->dim, 1);
        kdiff.initFromBox(bbext);
        diffusivity_f(boxindex, bbext, state, kdiff);
        AddLaplacian<T>(bbox, state, kdiff, grid->dim, grid->ld.Scale(), slope_impl);
        thread_max_kdiff = std::max<T>(thread_max_kdiff, kdiff.maxAbs());
      }
      if (velocity_f)
      {
        const BBox3 bbext = ExtendForDim(bbox, grid->dim, 3);
        Array3d<T> conc(bbext); conc.fill(state[bbext]);
        FaceVarArrays velocities(bbox, grid->dim);
        velocity_f(boxindex, bbext, bbox, conc, velocities, max_vel);
        //AddWENOAdvection<T>(bbox, conc, velocities.data(), grid->dim, grid->ld.Scale(), slope_expl, CONSERVATIVE, 5);
        AddKTAdvection<T>(bbox, conc, velocities.data(), grid->dim, grid->ld.Scale(), slope_expl, CONSERVATIVE);
      }
    }

    #pragma omp critical
    {
      max_src_expl = std::max(max_src_expl, thread_max_src_expl);
      max_src_impl = std::max(max_src_expl, thread_max_src_impl);
      max_kdiff = std::max(max_kdiff, thread_max_kdiff);
      max_vel = std::max(max_vel, thread_max_vel);
    }
  }

  euler_dt_srcexpl = 0.5/(1.e-13+max_src_expl);
  euler_dt_vel  = 0.8 * 0.5*grid->ld.Scale()/(grid->dim*max_vel+1.e-13);
  euler_dt_kdiff = 0.8 * 0.5*my::sqr(grid->ld.Scale())/(grid->dim*max_kdiff+1.e-13);
  euler_dt_srcimpl = 0.5/(1.e-13+max_src_impl);
  
  ctrl.euler_dt = my::min(euler_dt_srcexpl, euler_dt_vel);
}



template<class T_>
void ConvectionReactionDiffusionModel<T_>::insertLinearSystemCoefficients(int boxindex, const BBox3 &bbox, const State &rhs, const State& extrapolation, double identity_coeff, double operator_coeff, FiniteVolumeMatrixBuilder &mb)
{
  if (diffusivity_f)
  {
    const BBox3 bbext = ExtendForDim(bbox, grid->dim, 1);
    Array3d<float> kdiff(bbext);
    diffusivity_f(boxindex, bbext, extrapolation, kdiff);
    mb.AddDiffusion(bbox, kdiff, operator_coeff);
  }

  if (source_implicit_f)
  {
    Array3d<T> src_clin(bbox), src_val(bbox);
    source_implicit_f(boxindex, bbox, extrapolation, src_clin, src_val);
    FOR_BBOX3(p, bbox)
      mb.AddLocally(p, operator_coeff*src_clin(p), -operator_coeff*src_val(p));
  }

  FOR_BBOX3(p, bbox)
    mb.AddLocally(p, identity_coeff, rhs(p));
}



template<class T_>
void ConvectionReactionDiffusionModel<T_>::invertImplicitOperator(State &lhs,
                            const State &rhs,
                            double identity_coeff,
                            double operator_coeff,
                            Steppers::StepControl &ctrl,
                            State &extrapolation)
{
  /*
  * solve: (a I + b g) lhs = rhs
  * ->
  * {a I + b * (div (D*grad .) + clin)} lhs = rhs - b * cconst
  */

  ptree pt_params = make_ptree("output", 1)("max_iter", 100);

  StationaryDiffusionSolve(grid->ld, *mtboxes, grid->dim, lhs, boost::bind(&ConvectionReactionDiffusionModel::insertLinearSystemCoefficients, this, _1, _2, rhs, extrapolation, identity_coeff, operator_coeff, _3) , pt_params);
}



template class ConvectionReactionDiffusionModel<float>;
template class ConvectionReactionDiffusionModel<double>;

