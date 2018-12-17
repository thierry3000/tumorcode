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
#ifndef CONVECTION_DIFFUSION_SOLVER_H
#define CONVECTION_DIFFUSION_SOLVER_H

#include "continuum-flow.h"
#include "shared-objects.h"

#include "trilinos_linsys_construction.h"

#include "mwlib/time_stepper.h"
#include "mwlib/ptree_ext.h"

#include <boost/function.hpp>

/*------------------------------------------------------
------------------------------------------------------*/
/*
* solve
*  dc/dt = - div (c * v) + div (D * grad c) + clin * c + cconst
*
*  explicit: F[c] = -div (c*v)
*  implicit: G[c] = div (D*grad c) + clin * c + cconst
*/
template<class T_>
class ConvectionReactionDiffusionModel
{
  public:
    typedef T_ T;
    typedef ConvectionReactionDiffusionModel ThisType;
    typedef Array3d<T> State;
    enum { stepperFlags = Steppers::CAN_SOLVE_IMPLICIT };
    // read only plx
    T max_kdiff, max_src_expl, max_src_impl, max_vel, euler_dt_kdiff, euler_dt_srcexpl, euler_dt_srcimpl, euler_dt_vel;
  private:
    const ContinuumGrid *grid;
    const DynArray<BBox3> *mtboxes;

    typedef Steppers::Stepper<ThisType> Stepper;
    std::unique_ptr<Stepper> stepper;
    
  public:
    void calcSlope(const State &state,
                  State &slope,
                  Steppers::StepControl &ctrl);
    void calcSlopesIMEX(const State &state,
                        State &slope_expl,
                        State &slope_impl,
                        Steppers::StepControl &ctrl);
    void insertLinearSystemCoefficients(int boxindex,
                                        const BBox3 &bbox,
                                        const State &rhs,
                                        const State& extrapolation,
                                        double identity_coeff,
                                        double operator_coeff,
                                        FiniteVolumeMatrixBuilder &mb);
    void invertImplicitOperator(State &lhs,
                                const State &rhs,
                                double identity_coeff,
                                double operator_coeff,
                                Steppers::StepControl &ctrl,
                                State &extrapolation);

  public:
    ConvectionReactionDiffusionModel() {}
    ConvectionReactionDiffusionModel(const ContinuumGrid &grid_, const DynArray<BBox3> &mtboxes_, const ptree &params_) { init(grid_, mtboxes_, params_); }

    // set these callbacksz
    boost::function4<void, int, BBox3, ConstArray3d<T>, Array3d<float> > diffusivity_f;
    boost::function5<void, int, BBox3, ConstArray3d<T>, Array3d<T>, T &> source_explicit_f;
    boost::function5<void, int, BBox3, ConstArray3d<T>, Array3d<T>, Array3d<T> > source_implicit_f;
    boost::function6<void, int, BBox3, BBox3, Array3d<T>, FaceVarArrays &, T &> velocity_f;

    void init(const ContinuumGrid &grid_, const DynArray<BBox3> &mtboxes_, const ptree &params)
    {
      grid = &grid_;
      mtboxes = &mtboxes_;
      stepper = Stepper::create(params.get<string>("stepper"));
      stepper->init(*this);
    }

    Steppers::StepControl doStep(State &state, Steppers::StepControl ctrl)
    {
      max_kdiff = max_src_expl = max_src_impl = max_vel = 0.;
      euler_dt_kdiff = euler_dt_srcexpl = euler_dt_srcimpl = euler_dt_vel = 0;
      CopyBorder(state, grid->dim, 2);
      ctrl = stepper->doStep(state, ctrl);
      //cout << format("conv-diff-react-step: t=%f, euler_dt=%f, max_kdiff=%f, max_src_expl=%f, max_src_impl=%f, max_vel=%f") % ctrl.t % ctrl.euler_dt % max_kdiff % max_src_expl % max_src_impl % max_vel << endl;
      //cout << format("conv-diff-react-step: t=%f, euler_dt=%f, dt_kdiff=%f, dt_src_expl=%f, dt_src_impl=%f,dt_vel=%f") % ctrl.t % ctrl.euler_dt % euler_dt_kdiff % euler_dt_srcexpl % euler_dt_srcimpl % euler_dt_vel << endl;
      return ctrl;
    }
};



#endif
