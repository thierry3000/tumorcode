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
#ifndef TIME_STEPPER_V2_H
#define TIME_STEPPER_V2_H

#include "myAssert.h"
#include "helpers-mem.h"
#include "mwlib/math_ext.h"

#include <iostream>
#include <boost/format.hpp>


#include <boost/function.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/mpl/if.hpp>
#include <boost/ref.hpp>

#include <memory>

/**
@brief Time stepping schemes, aka. numerical integration schemes but we assume that the independent variable is time.

NOTE: Examples can be found in src/tests/test_stepper.cpp, and wrapped in python we have py/krebsutils/pdetests.py which comprises tests of almost the full differential equation solver suite.

Provided are Euler, Implicit Euler, 2nd order Runge-Kutter (improved Euler), and less common schemes designed for advection equations (RK3) and possibly other things.

These stepper objects are templated by the type of state that is propagated in time.
*/
namespace NewSteppers
{

// definition of mathematical operations for different types

typedef long long llong;


/**
@brief Abstration of mathematical operations on state variables

The default is to assume that the state type has corresponding members addScaled and initFrom. If it does not, this template may be specialized.
*/
template<class State>
struct Operations
{
  // @brief a = fa * a + fb * b
  void addScaled(double fa, State &a, double fb, const State &b) { a.addScaled(fb, b, fa); }
  // @brief Init a from b. Can be used to initialize a to zeros using b as prototype for dimensions. Or can make a copy of b. 
  void initFrom(State &a, const State &b, ConsMode mode) { a.initFrom(b, mode); }
};

/// specialization for double numbers
template<>
struct Operations<double>
{
  typedef double State;
  void addScaled(double fa, State &a, double fb, const State &b) {  a = fa*a + fb*b;  }
  void initFrom(State &a, const State &b, ConsMode mode)  { a = mode==ConsMode::AS_COPY ? b : 0.; }
};

/**
@brief messy structure to coordinate actions and for communication with callback functions.

For instance, when the F in du/dt = F(u,t) is evaluated, t is supplied via this struct. On the other hand
euler_dt should be set in the evaluation of F - normally in calcSlope.
**/
struct StepControl
{
  double dt, t; /// Filled in by the stepper algorithm.
  my::Smallest<double> euler_dt;  /// You need to fill in this value! It defines the time step width that we would have to use with the first order explicit Euler method.
  bool is_rewind_step; /// The previous step was invalid and the solution was rolled back to the state before that. This is required for BDF steppers

  StepControl() : dt(0), t(0), euler_dt(std::numeric_limits<double>::max()),
                  is_rewind_step(false) {}

  void update()
  {
    t += dt;
    is_rewind_step = false;
  }
};


#if 1
template<class State_>
struct Callback
{
  typedef State_ State;
  boost::function<void (const State&, State&, StepControl&)> calcSlope; /// Compute application of rhs operator; F in du/dt = F(u)
  boost::function<void (const State&, State&, State&, StepControl&)> calcSlopesIMEX;  /// Compute first and second operators (F and G). See paper on imex methods.
  /* Let our equation be du/dt = F(u_t, t), and the arguments to this function be
     lhs - left hand side,
     rhs - right hand side,
     alpha - some factor,
     beta - some factor,
     stepctrl - StepControl 
     extrapolation - used in vsimex method (see below), 
     
     Then the function should solve
     (alpha I + beta F)[lhs] = rhs, for lhs.  
     Return result in first argument. */
  boost::function<void (State&, const State&, double, double, StepControl &, State &)> invertImplicitOperator; 
};
#endif


template<class State_>
struct OpsCallback
{
  typedef State_ State;
  boost::function<void (double, State &, double, const State &)> addScaled;
  boost::function<void (State &, const State &, ConsMode)> initFrom;
};


class StepsizeTooLarge : public std::exception
{
public:
  virtual const char* what() const throw()
  {
    return "most recent euler_dt < previously set euler_dt";
  }
};

/**
 * @brief Used to remove the pointer quality from a pointer type. 
 * Unwrapped::type gives the type that is pointed to, if M is a pointer type. Else Unwrapped::type is equal to M.
 */
template<class M>
struct Unwrapped
{
  typedef M type;
  typedef const M& argument_type;
  inline static M& get(M &x) { return x; }
  inline static const M& get(const M &x) { return x; }
};

template<class M>
struct Unwrapped<M*>
{
  typedef M type;
  typedef M* argument_type;
  inline static M& get(M *x) { return *x; }
  inline static const M& get(const M* x) { return *x; }
};


#define STEPPER_TEMPLATE_ARGS  class Model_, class Ops_
#define STEPPER_TEMPLATE_ARGS2 Model_, Ops_


/**
@brief The base class that defines a common interface.

It takes template arguments for a model (Model) and operations (Ops). Model must contain methods
to compute the rate of change of the state at the current state at the current time (e.g. see above, 
the class Callbacks). Model should also contain a type called State, representing the system state under
consideration. 
Ops stands for Operations! must contain mathematical operations on the state type (e.g. see OpsCallback).
*/
template<STEPPER_TEMPLATE_ARGS>
class Stepper
{
public:
  /* Typedefs for the Model describing the F in du/dt = F(u, t), the Operations on u, and the type of u.
   * Unwrapped<...> is used to allow the Stepper to either own these components as member variables, or take pointers
   * to externally supplied objects. In the latter case Model_ will be a pointer type, for intance.
   */
  typedef typename Unwrapped<Model_>::type Model;
  typedef typename Unwrapped<Ops_>::type Ops;
  typedef typename Model::State State;

  
private:
  Model_ m_model;
  Ops_ m_ops;
  StepControl *m_ext_ctrl;
  State *m_ext_u; // m_ext_ctrl and this are both owned by the caller, and supplied via doStep.
  bool first_step;

protected:
  
  StepControl& ctrl() { return *m_ext_ctrl; }

  /**
    @brief This function has to do the actual work of advancing the solution. It must be overloaded in subclasses.
    ctrl.dt has been set with the requested step size
    ctrl.dt can be decreased if nessesary
    ctrl.euler_step will be set in the first call to calcSlope/calcSlopesIMEX
    ctrl.t must not be increased
  */
  virtual void doActualStep() = 0;
  virtual void internalInit() {}
  virtual double getStabilityLimitTimeFactor() const { return 1.; }  /// an algorithm specific factor (e.g. RK4 allows slightly larger steps than Euler).
  //TF multiplication throws overflow if ctrl.euler_dt is max doule
  void limitDt(StepControl &ctrl) const {
    if(ctrl.euler_dt > 10000000000000)
    {
      ctrl.dt = ctrl.dt;
    }
    else
    {
      ctrl.dt = std::min(ctrl.dt, 0.99*getStabilityLimitTimeFactor()*ctrl.euler_dt);
    }
  }
  void checkStepWidthAndThrow(StepControl &ctrl) const
  {
    if (ctrl.dt > getStabilityLimitTimeFactor()*ctrl.euler_dt)
      throw StepsizeTooLarge();
  }
  void initCloned(State &a, const State &b) { ops().initFrom(a, b, ConsMode::AS_COPY); }
  void initLike(State &a, const State &b) { ops().initFrom(a, b, ConsMode::CLEAN); }

public:
  Stepper() : m_ext_u(NULL), extrapolate(true), first_step(true), m_model(), m_ops(), m_ext_ctrl()
  {
  }

  virtual ~Stepper() {}

  bool extrapolate;

  /**
   * @brief Initialize the Stepper 
   */
  void set_model(typename Unwrapped<Model_>::argument_type model_) { m_model = model_; }
  void set_ops(typename Unwrapped<Ops_>::argument_type ops_) { m_ops = ops_; }
  /**
   * @brief Query data 
   */
  Model& model() { return Unwrapped<Model_>::get(m_model); }
  Ops& ops() { return Unwrapped<Ops_>::get(m_ops); }
  State& current_state() { return *m_ext_u; }
  
  bool doStep(State &state_, StepControl &ctrl_)
  {
    m_ext_u = &state_;
    m_ext_ctrl = &ctrl_;

    if (first_step)
    {
      internalInit();
      first_step = false;
    }
    
    limitDt(ctrl());

    myAssert(ctrl().dt > 0.);

    try {
      doActualStep();
    }
    catch (StepsizeTooLarge &)
    {
      return false;
    }

    ctrl().update();
    return true;
  };
};


#define STEPPER_DERIVED_SETUP \
  public:\
  typedef Stepper<STEPPER_TEMPLATE_ARGS2> Super; \
  typedef typename Super::State State; \
  using Super::model; \
  using Super::ops; \
  private:\
  using Super::ctrl; \
  using Super::current_state;

  

template<STEPPER_TEMPLATE_ARGS>
class ExplicitEuler : public Stepper<STEPPER_TEMPLATE_ARGS2>
{
  STEPPER_DERIVED_SETUP
public:
  void doActualStep()
  {
    State k1;
    model().calcSlope(current_state(), k1, ctrl());
    Super::limitDt(ctrl());
    ops().addScaled(1.,current_state(), ctrl().dt, k1);
  }
};


template<STEPPER_TEMPLATE_ARGS>
class ImprovedExplicitEuler : public Stepper<STEPPER_TEMPLATE_ARGS2>
{
  STEPPER_DERIVED_SETUP
public:
  void doActualStep()
  {
    State k1, k2, u1;
    State& u = current_state();

    model().calcSlope(u, k1, ctrl());
    
    Super::limitDt(ctrl());

    Super::initCloned(u1, u);
    ops().addScaled(1., u1, ctrl().dt, k1);
    
    model().calcSlope(u1, k2, ctrl());

    Super::checkStepWidthAndThrow(ctrl());
    
    ops().addScaled(1.,u, 0.5*ctrl().dt, k1);
    ops().addScaled(1.,u, 0.5*ctrl().dt, k2);
  }
};



template<STEPPER_TEMPLATE_ARGS>
class Rk4 : public Stepper<STEPPER_TEMPLATE_ARGS2>
{
  STEPPER_DERIVED_SETUP
public: 
  void initState(State &state, const State &u, const State &slope, double dt, int id)
  {
    Super::initCloned(state, u);
    ops().addScaled(1., state, dt, slope);
  }

  void calcSlope(State &slope_v, const State &u, const State &slope_u, double dt, int id)
  {
    State v;
    initState(v, u, slope_u, dt, id);
    model().calcSlope(v, slope_v, ctrl());
  }

  void doActualStep()
  {
    State k[4];
    State& u = current_state();
    
    model().calcSlope(u, k[0], ctrl());
    Super::limitDt(ctrl());
    calcSlope(k[1], u, k[0], 0.5 * ctrl().dt, 1);
    calcSlope(k[2], u, k[1], 0.5 * ctrl().dt, 2);
    calcSlope(k[3], u, k[2], ctrl().dt, 3);
    double weights[4] = { 1., 2., 2., 1. };

    Super::checkStepWidthAndThrow(ctrl());
    
    for (int i=0; i<4; ++i)
      ops().addScaled(1., u, weights[i]*ctrl().dt/6., k[i]);
  }

  virtual double getStabilityLimitTimeFactor() const { return 2.; }
};


template<STEPPER_TEMPLATE_ARGS>
class Rk3 : public Stepper<STEPPER_TEMPLATE_ARGS2>
{
  STEPPER_DERIVED_SETUP
public:
  void correctionStep(State &ul, const State &un, bool first, double eta)
  {
    State s;
    model().calcSlope(ul, s, ctrl());
    if (first)
      Super::limitDt(ctrl());
    ops().addScaled(1., ul, ctrl().dt, s);
    ops().addScaled(1.-eta, ul, eta, un);
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
    State& u = current_state();
    State ul; Super::initCloned(ul, u);
    correctionStep(ul, u, true, 0.);
    correctionStep(ul, u, false, 3./4.);
    correctionStep(ul, u, false, 1./3.);

    Super::checkStepWidthAndThrow(ctrl());
    
    ops().addScaled(0., u, 1., ul);
  }
};



/* This is from [Asher et al. (1997) "Implicit-explicit Runge-Kutta methods for time-dependent partial differential equations"]
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
  virtual void internalInit()
  {
    needs_init = true;
    for (int i=0; i<2; ++i) {
      var_u[i].reset(); var_f[i].reset(); var_g[i].reset();
    }
  }

  
public:
  VSImExBdf2() : Super() {}

  void doActualStep()
  {
    if (ctrl().is_rewind_step) {
      std::cout << "Warning: vsimexbdf2 rewind!" << std::endl;
      internalInit();
    }

    State& u = current_state();
    var_u[1].reset(new State());
    var_g[1].reset(new State());
    var_f[1].reset(new State());

    Super::initCloned(*var_u[1], u);
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
    
    Super::limitDt(ctrl());
    bool is_first = needs_init;

    calcCoefficients();
    State rhs;
    Super::initLike(rhs, *var_u[1]);

    for (int i=0; i<2; ++i)
    {
      myAssert(var_f[i].get() || beta[i] == 0.);
      myAssert(var_g[i].get() || gamma[i] == 0.);
      if (var_u[i].get()) ops().addScaled(1., rhs, -alpha[i], *var_u[i]);
      if (var_f[i].get()) ops().addScaled(1., rhs, step_sizes[1]*beta[i], *var_f[i]);
      if (var_g[i].get()) ops().addScaled(1., rhs, step_sizes[1]*gamma[i], *var_g[i]);
    }

    // free some memory, these correspond to the last step, so they are obsolete now that the rhs has been computed
    var_f[0].reset();
    var_g[0].reset();
    
    State* extrapolation = NULL; 
    if (Super::extrapolate)
    {
      if (is_first)
      {
         // copy current state
        var_u[0].reset(new State());
        Super::initCloned(*var_u[0], *var_u[1]);
      }
      else
      {
        double w = step_sizes[1]/step_sizes[0]; // curr / old step
        ops().addScaled(-w, *var_u[0], 1.+w, *var_u[1]); // last_state = extrapolation = 2x current state - 1x last state
      }
      extrapolation = var_u[0].get();
    }
    else
    {
      // default is extrapolation = current state, model is not allowed to free extrapolation
      extrapolation = var_u[1].get();
      var_u[0].reset();
    }

    Super::checkStepWidthAndThrow(ctrl());

    model().invertImplicitOperator(u, rhs, alpha[2], -step_sizes[1]*gamma[2], ctrl(), *extrapolation);

    var_u[0] = std::move(var_u[1]);
    var_f[0] = std::move(var_f[1]);
    var_g[0] = std::move(var_g[1]);
    step_sizes[0] = step_sizes[1];
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


#if 0
template<STEPPER_TEMPLATE_ARGS>
class RkImEx : public Stepper<STEPPER_TEMPLATE_ARGS2>
{
  STEPPER_DERIVED_SETUP

  int order;
  enum { MAX_ORDER=4 };
  double matrix_a[MAX_ORDER][MAX_ORDER], matrix_b[MAX_ORDER][MAX_ORDER];
  double c_a[MAX_ORDER], c_b[MAX_ORDER], w_a[MAX_ORDER], w_b[MAX_ORDER];

protected:
  void internalInit()
  {
    // pareschi & russo 2003, implicit-explicit runge-kutta schemes and applications to hyperbolic systems with relaxation
    ClearMem(matrix_a, MAX_ORDER*MAX_ORDER);
    ClearMem(matrix_b, MAX_ORDER*MAX_ORDER);
    ClearMem(c_a, MAX_ORDER);
    ClearMem(c_b, MAX_ORDER);
    ClearMem(w_a, MAX_ORDER);
    ClearMem(w_b, MAX_ORDER);
    matrix_a[1,0] = 1.;
    c_a[1] = 1.;
    w_a[0] = w_a[1] = 0.5;
    const double gamma = 0.29289321881345254;
    matrix_b[0,0] = gamma;
    matrix_b[1,0] = 1.-2.*gamma;
    matrix_b[1,1] = gamma;
    c_b[0] = gamma;
    c_b[1] = 1.-gamma;
    w_b[0] = w_b[1] = 0.5;
    order = 2;
  }
  
public:
  void doActualStep()
  {
    // a, and b correspond to the explicit and implicit operators, respectively
    State slope_a[MAX_ORDER], slope_b[MAX_ORDER];
    State& u = current_state();

    // not possible to compute euler_dt during evaluation of the slope since
    // the implicit operator inversion is done first wherefor the time step is needed
    
    for (int i=0; i<order; ++i)
    {
      State u_intermediate, rhs;
      ops().initCloned(rhs, *u);
      for (int j=0; j<i; ++j)
      {
        ops().addScaled(1., rhs, dt*matrix_a[i][j], slope_a[j]);
        ops().addScaled(1., rhs, dt*matrix_b[i][j], slope_b[j]);
      }
      double t_a = ctrl().t + c_a[i] * dt;
      double t_b = ctrl().t + c_b[i] * dt;
      model().invertImplicitOperator(u_intermediate, rhs, 1., -dt*matrix_b[i][i], ctrl(), t_b, u);
      model().calcSlopesIMEX(u_intermediate, slope_a[i], slope_b[i], ctrl());
    }
    for (int i=0; i<order; ++i)
    {
      ops().addScaled(1., *u, dt*w_a[i], slope_a[i]);
      ops().addScaled(1., *u, dt*w_b[i], slope_b[i]);
    }
  }
};
#endif

enum StepperFlags
{
  NO_IMPLICIT = 4,
};

template<STEPPER_TEMPLATE_ARGS, int flags = 0>
struct StepperFactory
{
  typedef Stepper<STEPPER_TEMPLATE_ARGS2> StepperType;
  typedef std::unique_ptr<StepperType> PointerType;

private:
  static PointerType create_(const std::string &name, boost::mpl::false_)
  {
    PointerType p;
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

  static PointerType create_(const std::string &name, boost::mpl::true_)
  {
    PointerType p;
    if (name == "vsimexbdf2")
      p.reset(new VSImExBdf2<STEPPER_TEMPLATE_ARGS2>());
    else
      p = create_(name, boost::mpl::false_());
    return p;
  }

public:
  static PointerType create(const std::string &name)
  {
    PointerType p = create_(name, boost::mpl::bool_<!(flags & NO_IMPLICIT)>());
    if (!p.get()) throw std::runtime_error(boost::str(boost::format("unable to create stepper %s") % name));
    return p;
  }
};


template<STEPPER_TEMPLATE_ARGS, int flags>
inline typename StepperFactory<STEPPER_TEMPLATE_ARGS2, flags>::PointerType create(const std::string &name, typename Unwrapped<Model_>::argument_type model_, typename Unwrapped<Ops_>::argument_type ops_)
{
  typedef StepperFactory<Model_,Ops_,flags> Factory;
  typename Factory::PointerType p = Factory::create(name);
  p->set_model(model_);
  p->set_ops(ops_);
  return p;
}

template<STEPPER_TEMPLATE_ARGS, int flags>
inline typename StepperFactory<STEPPER_TEMPLATE_ARGS2, flags>::PointerType create(const std::string &name)
{
  typedef StepperFactory<Model_,Ops_,flags> Factory;
  typename Factory::PointerType p = Factory::create(name);
  return p;
}






















template<STEPPER_TEMPLATE_ARGS>
class Stepper2
{
public:
  typedef typename Unwrapped<Model_>::type Model;
  typedef typename Unwrapped<Ops_>::type Ops;
  typedef typename Model::State State;


private:
  Model_ m_model;
  Ops_ m_ops;
  bool first_step;

protected:

  virtual void doActualStep(State &u, double t, double &dt) = 0;
  virtual void internalInit() {}
  virtual double getStabilityLimitTimeFactor() const { return 1.; }
  void limitDt(double &dt) const
  {
    dt = std::min(dt, 0.99*getStabilityLimitTimeFactor()*model().getEulerDt());
  }
  void checkStepWidthAndThrow(double dt) const
  {
    if (dt > getStabilityLimitTimeFactor()*model().getEulerDt())
      throw StepsizeTooLarge();
  }
  void initCloned(State &a, const State &b) { ops().initFrom(a, b, ConsMode::AS_COPY); }
  void initLike(State &a, const State &b) { ops().initFrom(a, b, ConsMode::CLEAN); }

public:
  Stepper2() : extrapolate(true), first_step(true), m_model(), m_ops()
  {
  }

  virtual ~Stepper2() {}

  bool extrapolate;

  void set_model(typename Unwrapped<Model_>::argument_type model_) { m_model = model_; }
  void set_ops(typename Unwrapped<Ops_>::argument_type ops_) { m_ops = ops_; }
  Model& model() { return Unwrapped<Model_>::get(m_model); }
  const Model& model() const { return Unwrapped<Model_>::get(m_model); }
  Ops& ops() { return Unwrapped<Ops_>::get(m_ops); }

  bool doStep(State &u, double &t, double dt)
  {
    if (first_step)
    {
      internalInit();
      first_step = false;
    }

    myAssert(dt > 0.);

    try {
      doActualStep(u, t, dt);
    }
    catch (StepsizeTooLarge &)
    {
      return false;
    }

    t += dt;
    return true;
  };
};


#define STEPPER_DERIVED_SETUP2 \
  public:\
  typedef Stepper2<STEPPER_TEMPLATE_ARGS2> Super; \
  typedef typename Super::State State; \
  using Super::model; \
  using Super::ops; \
  private:


template<STEPPER_TEMPLATE_ARGS>
class ExplicitEuler2 : public Stepper2<STEPPER_TEMPLATE_ARGS2>
{
  STEPPER_DERIVED_SETUP2
public:
  void doActualStep(State &u, double t, double &dt)
  {
    Super::limitDt(dt);
    State k1;
    model().calcSlope(u, t, k1);
    Super::limitDt(dt);
    ops().addScaled(1.,u, dt, k1);
  }
};


template<STEPPER_TEMPLATE_ARGS>
class ImprovedExplicitEuler2 : public Stepper2<STEPPER_TEMPLATE_ARGS2>
{
  STEPPER_DERIVED_SETUP2
public:
  void doActualStep(State &u, double t, double &dt)
  {
    Super::limitDt(dt);
    State k1, k2, u1;

    model().calcSlope(u, t, k1);

    Super::limitDt(dt);

    Super::initCloned(u1, u);
    ops().addScaled(1., u1, dt, k1);

    model().calcSlope(u1, t+dt, k2);

    Super::checkStepWidthAndThrow(dt);

    ops().addScaled(1.,u, 0.5*dt, k1);
    ops().addScaled(1.,u, 0.5*dt, k2);
  }
};


template<STEPPER_TEMPLATE_ARGS>
class ImplicitEuler2 : public Stepper2<STEPPER_TEMPLATE_ARGS2>
{
  STEPPER_DERIVED_SETUP2
public:
  void doActualStep(State &u, double t, double &dt)
  {
    State rhs;
    Super::initCloned(rhs, u);

    model().invertImplicitOperator(u, rhs, 1., -dt, t+dt);
  }
};


template<STEPPER_TEMPLATE_ARGS>
class CrankNicholson2 : public Stepper2<STEPPER_TEMPLATE_ARGS2>
{
  STEPPER_DERIVED_SETUP2
public:
  void doActualStep(State &u, double t, double &dt)
  {
    State rhs, k1;
    model().calcSlope(u, t+dt, k1);
    
    Super::initCloned(rhs, u);
    ops().addScaled(1., rhs, 0.5*dt, k1);

    model().invertImplicitOperator(u, rhs, 1., -0.5*dt, t+dt);
  }
};


}

#undef STEPPER_TEMPLATE_ARGS2
#undef STEPPER_TEMPLATE_ARGS
#undef STEPPER_DERIVED_SETUP

#endif
