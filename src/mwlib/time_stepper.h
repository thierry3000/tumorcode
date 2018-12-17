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
#ifndef TIME_STEPPER_H
#define TIME_STEPPER_H

#include "myAssert.h"
#include "refcounted.h"

#include <cmath>
#include <iostream>
#include <boost/format.hpp>
#include <boost/assign/std/vector.hpp>

#include <boost/function.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/mpl/if.hpp>

using boost::format;
using boost::str;

namespace Steppers
{

// definition of mathematical operations for different types

typedef long long llong;

template<class State>
struct Operations
{
  static void addScaled(double fa, State &a, double fb, const State &b) { a.addScaled(fb, b, fa); }
  static void initCloned(State &a, const State &b) { a.initCloned(b); }
  static void initLike(State &a, const State &b) { a.initLike(b); }
  static void setId(State &a, llong id) { a.setId(id); }
  static llong  getId(const State &a) { return a.getId(); }
  static double errorMeasure(double atol, double rtol, const State &y0, const State &y1, const State &err) {
    return State::errorMeasure(atol, rtol, y0, y1, err);
  }
};


template<>
struct Operations<double>
{
  typedef double State;
  static void addScaled(double fa, State &a, double fb, const State &b) {  a = fa*a + fb*b;  }
  static void initCloned(State &a, const State &b)  { a = b; }
  static void initLike(State &a, const State &b) { a = 0; }
  static void setId(State &a, llong id) {  }
  static llong  getId(const State &a) { return 0; }
  static double errorMeasure(double atol, double rtol, const State &y0, const State &y1, const State &err)
  {
    double z = 1./(rtol*std::max(std::abs(y0), std::abs(y1))+atol) * err;
    return std::abs(z);
  }
};



enum StoragePolicy
{
  STORAGE_REF = 1,
  STORAGE_COPY = 2,
};


#define STEPPER_TEMPLATE_ARGS  class Model_, int STORAGE = STORAGE_REF
#define STEPPER_TEMPLATE_ARGS3 class Model_, int STORAGE
#define STEPPER_TEMPLATE_ARGS2 Model_, STORAGE

enum StepperFlags
{
  HAS_PRESTEP = 1,
  HAS_POSTSTEP = 2,
  CAN_SOLVE_IMPLICIT = 4,
  HAS_CUSTOM_STATE_OPERATORS = 8
};



#define STEPPER_MODELARGS_PRESTEP State&, StepControl&

#if 0
template<class State_>
class ModelBase
{
  public:
    typedef State_ State;
    enum { stepperFlags = HAS_PRESTEP | HAS_POSTSTEP | CAN_SOLVE_IMPLICIT };
    virtual void preStep(STEPPER_MODELARGS_PRESTEP) {}
    virtual void postStep(State&, const StepControl&) {}
    virtual void calcSlope(const State&, State&, StepControl&) { assert(false); }
    virtual void calcSlopesIMEX(const State&, State&, State&, StepControl&) { assert(false); }
    virtual void invertImplicitOperator(State&, const State&, double, double, StepControl &, State &) {}
};
#endif


struct StepControl
{
  double dt, t, euler_dt;
  bool stopflag, // Set this at any point. actually ignored by
       is_rewind_step; // The previous step was invalid and the solution was rolled back to the state before that. This is required for BDF steppers

  StepControl() : dt(0), t(0), euler_dt(std::numeric_limits<double>::max()),
                  stopflag(false), is_rewind_step(false) {}

  void update()
  {
    t += dt;
    is_rewind_step = false;
  }
};


#if 1
template<class State_>
class Callback
{
  public:
    typedef State_ State;
    enum { stepperFlags = HAS_PRESTEP | HAS_POSTSTEP | CAN_SOLVE_IMPLICIT };
    boost::function2<void, STEPPER_MODELARGS_PRESTEP> preStep;
    boost::function2<void, State&, const StepControl&> postStep;
    boost::function3<void, const State&, State&, StepControl&> calcSlope;
    boost::function4<void, const State&, State&, State&, StepControl&> calcSlopesIMEX;
    boost::function6<void, State&, const State&, double, double, StepControl &, State &> invertImplicitOperator;

    Callback()
    {
      preStep = preStepDefault;
      postStep = postStepDefault;
      calcSlope = calcSlopeDefault;
      calcSlopesIMEX = calcSlopesIMEXDefault;
      invertImplicitOperator = invertImplicitOperatorDefault;
    }

    static void preStepDefault(STEPPER_MODELARGS_PRESTEP) {}
    static void postStepDefault(State&, const StepControl&) {}
    static void calcSlopeDefault(const State&, State&, StepControl&) { assert(false); }
    static void calcSlopesIMEXDefault(const State&, State&, State&, StepControl&) { assert(false); }
    static void invertImplicitOperatorDefault(State&, const State&, double, double, StepControl &, State &) { assert(false); }
};
#endif

#undef STEPPER_MODELARGS_PRESTEP

/**
 *\brief Abstract stepper class
 * later we use this to inhere differen time stepping method
 * and explicitly implement them 
 */
template<STEPPER_TEMPLATE_ARGS>
class Stepper
{
public:
  typedef Model_ Model;
  typedef typename Model::State State;
  typedef Operations<State> Ops;
  typedef std::vector<boost::shared_ptr<State> > BdfStateArray;
  typedef typename boost::mpl::if_c<STORAGE==STORAGE_REF, StorageRef<Model>, StorageCopy<Model> >::type ModelStorage;
  
protected:
  ModelStorage m_model;
  State *u;
  StepControl m_ctrl;

  void incrId(State &dst, const State &src, llong val) {  Ops::setId(dst, Ops::getId(src)+val);  }
  StepControl& ctrl() { return m_ctrl; }

  /**
    this function has to do the actual work of advancing the solution .
    ctrl.dt has been set with the requested step size
    ctrl.dt can be decreased if nessesary
    ctrl.euler_step will be set in the first call to calcSlope/calcSlopesIMEX
    ctrl.t must not be increased
  */
  virtual void doActualStep() = 0;

  virtual void internalInit() {}
  
  void preStep(boost::mpl::true_) { model().preStep(*u, ctrl()); }
  void preStep(boost::mpl::false_) {}
  void postStep(boost::mpl::true_) { model().postStep(*u, const_cast<const StepControl&>(ctrl())); }
  void postStep(boost::mpl::false_) {}
  
public:
  Stepper() : u(NULL), extrapolate(true) {}
  virtual ~Stepper() {}

  bool extrapolate;
  
  void init(Model &model_)
  {
    m_model = ModelStorage(model_);
    init();
  }

  void init()
  {
    internalInit();
  }
  /** calculate serveral callbacks */
  StepControl doStep(State &state, const StepControl &ctrl_)
  {
    u = &state;
    m_ctrl = ctrl_;
    m_ctrl.euler_dt = std::numeric_limits<double>::max();
    //ctrl.dt = dt;
    
    preStep(boost::mpl::bool_<Model::stepperFlags & HAS_PRESTEP>());
    doActualStep();
    postStep(boost::mpl::bool_<Model::stepperFlags & HAS_POSTSTEP>());

    incrId(*u, *u, 1000);
    //ctrl.t += ctrl.dt;
    //ctrl.iter_num += 1;

    u = NULL;
    m_ctrl.update();
    return m_ctrl;
    //return ctrl().dt;
  };

  Model& model() { return m_model; }
  typedef Stepper<STEPPER_TEMPLATE_ARGS2> Base;
  typedef std::unique_ptr<Base> BasePointer;
  static BasePointer create(const std::string &name)
  {
    BasePointer p = create(name, boost::mpl::bool_<Model::stepperFlags & CAN_SOLVE_IMPLICIT>());
    if (!p.get()) throw std::runtime_error(str(format("unable to create stepper %s") % name));
    return p;
  }
  
private:  
  static BasePointer create(const std::string &name, boost::mpl::true_);
  static BasePointer create(const std::string &name, boost::mpl::false_);
};


#define STEPPER_DERIVED_SETUP \
  public:\
  typedef Stepper<STEPPER_TEMPLATE_ARGS2> Super; \
  typedef typename Super::Model Model; \
  typedef typename Super::State State; \
  typedef typename Super::Ops Ops; \
  private:\
  using Super::ctrl; \
  using Super::u; \
  using Super::model; \
  using Super::incrId; 

  //using Super::nextStateId;
  
/**
 * \brief Explitic Euler timestepping
 * The most simple one.
 * k1 is the slope. State u is updated inplace.
 */
template<STEPPER_TEMPLATE_ARGS>
class ExplicitEuler : public Stepper<STEPPER_TEMPLATE_ARGS2>
{
  STEPPER_DERIVED_SETUP
public:
  void doActualStep()
  {
    State k1;
    model().calcSlope(*u, k1, ctrl());
    ctrl().dt = std::min(ctrl().dt, 0.99*ctrl().euler_dt);
    Ops::addScaled(1.,*u, ctrl().dt, k1);// Ops is operations!
  }
};

/**
 * ImprovedExplicitEuler method
 * I think we have some intermitant states.
 */
template<STEPPER_TEMPLATE_ARGS>
class ImprovedExplicitEuler : public Stepper<STEPPER_TEMPLATE_ARGS2>
{
  STEPPER_DERIVED_SETUP
public:
  void doActualStep()
  {
    State k1, k2, u1;

    model().calcSlope(*u, k1, ctrl());
    assert(ctrl().euler_dt > 0.);

    ctrl().dt = std::min(ctrl().dt, 0.99*ctrl().euler_dt);
    double first_step_size = ctrl().dt;

    Ops::initCloned(u1, *u);
    Ops::addScaled(1., u1, ctrl().dt, k1);
    Super::incrId(u1, *u, 1);
    
    model().calcSlope(u1, k2, ctrl());

    ctrl().dt = first_step_size;

    Ops::addScaled(1.,*u, 0.5*ctrl().dt, k1);
    Ops::addScaled(1.,*u, 0.5*ctrl().dt, k2);
  }
};

/**
 * Ruge Kutta 4 stepper
 */
template<STEPPER_TEMPLATE_ARGS>
class Rk4 : public Stepper<STEPPER_TEMPLATE_ARGS2>
{
  STEPPER_DERIVED_SETUP
public:
  void initState(State &state, const State &u, const State &slope, double dt, int id)
  {
    Ops::initCloned(state, u);
    Ops::addScaled(1., state, dt, slope);
    incrId(state, u, id);
  }

  void calcSlope(State &slope_v, const State &u, const State &slope_u, double dt, int id)
  {
    State v;
    initState(v, u, slope_u, dt, id);
    StepControl ctrl = this->ctrl();
    ctrl.t += dt;
    model().calcSlope(v, slope_v, ctrl);
  }

  void doActualStep()
  {
    State k[4];

    model().calcSlope(*u, k[0], ctrl());
    assert(ctrl().euler_dt > 0.);
    
    ctrl().dt = std::min(ctrl().dt, 0.99*ctrl().euler_dt * 2.);
    double dt = ctrl().dt;

    calcSlope(k[1], *u, k[0], 0.5 * dt, 1);
    calcSlope(k[2], *u, k[1], 0.5 * dt, 2);
    calcSlope(k[3], *u, k[2], dt, 3);
    double weights[4] = { 1., 2., 2., 1. };

    myAssert(dt == ctrl().dt);
    for (int i=0; i<4; ++i)
      Ops::addScaled(1., *u, weights[i]*dt/6., k[i]);
  }
};
/**
 * Ruge Kutta 3 stepper
 */
template<STEPPER_TEMPLATE_ARGS>
class Rk3 : public Stepper<STEPPER_TEMPLATE_ARGS2>
{
  STEPPER_DERIVED_SETUP
public:
  void correctionStep(State &ul, const State &un, bool first, double eta)
  {
    double dt_backup = ctrl().dt;
    State s;
    model().calcSlope(ul, s, ctrl());
    assert(ctrl().euler_dt > 0.);
    
    if (first) {
      ctrl().dt = std::min(ctrl().dt, 0.99*ctrl().euler_dt);
    } else {
      ctrl().dt = dt_backup;
    }
    Ops::addScaled(1., ul, ctrl().dt, s);
    Ops::addScaled(1.-eta, ul, eta, un);
    incrId(ul, ul, 1);
  }

  void doActualStep()
  {
    /*
      from New High-Resolution Central Schemes for Nonlinear Conservation Laws and Convection-Diffusion Equation (Kurganov 2000)
      u(1) = un + dt * C[un]
      u(2) = eta1 * un + (1-eta1)(u(1) + dt C[u(1)])
      u(3) = eta2 * un + (1-eta2)(u(2) + dt C[u(2)])
      un+1 = u(3)
    */
    State ul; Ops::initCloned(ul, *u);
    correctionStep(ul, *u, true, 0.);
    correctionStep(ul, *u, false, 3./4.);
    correctionStep(ul, *u, false, 1./3.);
    Ops::addScaled(0., *u, 1., ul); // it's a copy
  }
};
/**
 * Implizit Explicit Backward Differentiation Formulas ????
 * ????VSImExBdf2 stepper
 * where does it come from
 * obviously MW like it
 */
template<STEPPER_TEMPLATE_ARGS>
class VSImExBdf2 : public Stepper<STEPPER_TEMPLATE_ARGS2>
{
  STEPPER_DERIVED_SETUP
  bool needs_init;

  std::unique_ptr<State> var_u[2];
  std::unique_ptr<State> var_f[2];
  std::unique_ptr<State> var_g[2];
  double step_sizes[2];
  double alpha[3], beta[3], gamma[3];

protected:
  void internalInit()
  {
    needs_init = true;
    for (int i=0; i<2; ++i) {
      var_u[i].reset(); var_f[i].reset(); var_g[i].reset();
    }
  }
public:
  void doActualStep()
  {
    if (ctrl().is_rewind_step) {
      std::cout << "Warning: vsimexbdf2 rewind!" << std::endl;
      internalInit();
    }
    var_u[1].reset(new State());
    var_g[1].reset(new State());
    var_f[1].reset(new State());
    Ops::initCloned(*var_u[1], *u);  
    /* arguments are:
        current state,
        evaluation of (implicit) operator g 
        evaluation of (explicit) operator f
        control struct
      todo:
        evaluate g and f,
        also limit the timestep dt in ctrl
    */
    model().calcSlopesIMEX(*var_u[1], *var_f[1], *var_g[1], ctrl());
    assert(ctrl().euler_dt > 0.);
    
    ctrl().dt = std::min(ctrl().dt, 0.99*ctrl().euler_dt); // numerical experiments indicate that the stability region is the same as with euler and backward euler
    bool is_first = needs_init;

    calcCoefficients();
    State rhs;
    Ops::initLike(rhs, *var_u[1]);

    for (int i=0; i<2; ++i)
    {
      myAssert(var_f[i].get() || beta[i] == 0.);
      myAssert(var_g[i].get() || gamma[i] == 0.);
      if (var_u[i].get()) Ops::addScaled(1., rhs, -alpha[i], *var_u[i]);
      if (var_f[i].get()) Ops::addScaled(1., rhs, step_sizes[1]*beta[i], *var_f[i]);
      if (var_g[i].get()) Ops::addScaled(1., rhs, step_sizes[1]*gamma[i], *var_g[i]);
    }

    // free some memory, these correspond to the last step, so they are obsolete now that the rhs has been computed
    var_f[0].reset();
    var_g[0].reset();
    
    /* arguments are:
        lhs aka state at next time (to be computed),
        rhs,
        coefficient before identity; a,
        coefficient before operator; b,
        current state,
        control struct
      todo:
        solve: (a I + b g) lhs = rhs
    */
    State* extrapolation = NULL; 
    if (Super::extrapolate)
    {
      if (is_first)
      {
         // copy current state
        var_u[0].reset(new State());
        Ops::initCloned(*var_u[0], *var_u[1]);
      }
      else
      {
        double w = step_sizes[1]/step_sizes[0]; // curr / old step
        Ops::addScaled(-w, *var_u[0], 1.+w, *var_u[1]); // last_state = extrapolation = 2x current state - 1x last state
      }
      extrapolation = var_u[0].get();
      incrId(*extrapolation, *u, 1);
    }
    else
    {
      // default is extrapolation = current state, model is not allowed to free extrapolation
      extrapolation = var_u[1].get();
      var_u[0].reset();
    }

    model().invertImplicitOperator(*u, rhs, alpha[2], -step_sizes[1]*gamma[2], ctrl(), *extrapolation);

    var_u[0] = std::move(var_u[1]);
    var_f[0] = std::move(var_f[1]);
    var_g[0] = std::move(var_g[1]);
    step_sizes[0] = step_sizes[1];

    //Ops::setId(*u, last_state_id);
  }

private:
  void calcCoefficients()
  {
    step_sizes[1] = ctrl().dt;
    if (needs_init)
    {
      // for gamma = 1/2, first order
      alpha[0] = beta[0] = gamma[0] = 0.;
      alpha[1] = -1.;
      alpha[2] = 1.;
      beta[1] = 1.;
      beta[2] = 0.;
      gamma[1] = 0.5;
      gamma[2] = 0.5;
//       gamma[1] = 0.;
//       gamma[2] = 1.;
      needs_init = false;
    }
    else
    {
      double w = step_sizes[1]/step_sizes[0];
      // for gamma, c = 1/2, 1/8, second order
      alpha[0] = 0.;
      alpha[1] = -1.;
      alpha[2] = 1.;
      beta[0] = -0.5*w;
      beta[1] = 0.5*(2.+w);
      beta[2] = 0.;
      gamma[0] = 1./16.;
      gamma[1] = (7.*w - 1.)/(16.*w);
      gamma[2] = (8.*w + 1.)/(16.*w);
    }
  }
};


template<STEPPER_TEMPLATE_ARGS3>
typename Stepper<STEPPER_TEMPLATE_ARGS2>::BasePointer Stepper<STEPPER_TEMPLATE_ARGS2>::create(const std::string &name, boost::mpl::false_)
{
  BasePointer p;
  if (name == "euler")
    p.reset(new ExplicitEuler<STEPPER_TEMPLATE_ARGS2>());
  else if (name == "impeuler")
    p.reset(new ImprovedExplicitEuler<STEPPER_TEMPLATE_ARGS2>());
  else if (name == "rk4")
    p.reset(new Rk4<STEPPER_TEMPLATE_ARGS2>());
  else if (name == "rk3")
    p.reset(new Rk3<STEPPER_TEMPLATE_ARGS2>());
  return p;
}

template<STEPPER_TEMPLATE_ARGS3>
typename Stepper<STEPPER_TEMPLATE_ARGS2>::BasePointer Stepper<STEPPER_TEMPLATE_ARGS2>::create(const std::string &name, boost::mpl::true_)
{
  BasePointer p;
  if (name == "vsimexbdf2")
    p.reset(new VSImExBdf2<STEPPER_TEMPLATE_ARGS2>());
  else
    p = create(name, boost::mpl::false_());
  return p;
}


template<class Model>
static typename Stepper<Model>::BasePointer createRefBased(const std::string &name)
{
  return Stepper<Model>::create(name);
}

template<class State>
static typename Stepper<Callback<State>, STORAGE_COPY>::BasePointer createCallbackBased(const std::string &name)
{
  typename Stepper<Callback<State>, STORAGE_COPY>::BasePointer p = Stepper< Callback<State>, STORAGE_COPY>::create(name);
  p->init();
  return p;
}



}

#undef STEPPER_TEMPLATE_ARGS3
#undef STEPPER_TEMPLATE_ARGS2
#undef STEPPER_TEMPLATE_ARGS
#undef STEPPER_DERIVED_SETUP

#endif
