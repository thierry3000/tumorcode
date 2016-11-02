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
#include <fenv.h>
#include "common.h"
#include "time_stepper_utils_new.h"

using namespace NewSteppers;



typedef double State;


class Model : boost::noncopyable
{
public:
  typedef double State; 
  Model(double lambda) : lambda(lambda), max_dt(1.), go_implicit(false) {}
  
  void calcSlope(const State &state, State &slope, StepControl &ctrl)
  {
    slope = state*lambda;
    ctrl.euler_dt = my::Smallest<double>(0.5/std::abs(lambda));
  }

  void calcSlopesIMEX(const State &state, State &slopef, State &slopeg, StepControl &ctrl)
  {
    if (go_implicit)
    {
      slopeg = state*lambda;
      slopef = 0.;
    }
    else
    {
      //explicit mode
      slopef = state*lambda;
      slopeg = 0.;
    }
    ctrl.euler_dt = my::Smallest<double>(0.5/std::abs(lambda));
  }

  void invertImplicitOperator(State &lhs, const State &rhs, double identity_coeff, double operator_coeff, StepControl &ctrl, State &extrapolation)
  {
    if (go_implicit)
    {
      lhs = rhs/(identity_coeff + operator_coeff*lambda);
    }
    else
    {
      // lhs c * I = rhs
      lhs = rhs / identity_coeff;
    }
  }

  void getDataNames(DynArray<string> &names) const
  {
    names.push_back("x");
    names.push_back("y");
  };

  void storeState(double t, State &state, DynArray<double> &dst)
  {
    dst[0] = state;
  }

  bool go_implicit;
  double lambda;
  double max_dt;
};


enum Method
{
  EULER = 0,
  IMP_EULER,
  RK4,
  RK3Shu,
  VS_IMEX_BDF2,
  VS_IMEX_BDF2_Implicit,
  //DUMKA3,
  N_METHODS
};

static  const char* method_name[] = {
  "euler",
  "imp-euler",
  "rk4",
  "rk3shu",
  "VSImExBdf2",
  "VSImExBdf2_Implicit",
  //"Dumka3 (w. Step Control)"
};



typedef NewSteppers::StepperFactory<Model*, Operations<double> > Factory;
typedef Factory::StepperType MyStepper;

// bool doStep_(MyStepper &stepper, Model::State &state, StepControl &ctrl)
// {
//   ctrl = stepper.doStep(state, ctrl);
//   return true;
// }


void run(h5cpp::Group g, const ptree &params, int id, double lambda)
{
  double y0 = 1.;

  Model model(lambda);
  State state(y0);
  ObserverOde observer(model, params);

  const int method = params.get<int>("method");
  const double dt = params.get<double>("out_intervall");

  observer.g = g = g.require_group(str(format("study%02i") % id));
  g.attrs().set("method", method_name[method]);
  g.attrs().set("dt", dt);

  std::string name;
  switch(method)
  {
    case EULER: name = "euler"; break;
    case IMP_EULER: name = "impeuler"; break;
    case RK4: name = "rk4"; break;
    case RK3Shu: name = "rk3"; break;
    case VS_IMEX_BDF2: name = "vsimexbdf2"; break;
    case VS_IMEX_BDF2_Implicit: 
      name = "vsimexbdf2"; 
      model.go_implicit = true;
      break;
  }


  std::auto_ptr<MyStepper> stepper(Factory::create(name));
  stepper->set_model(&model);

  auto doStep = [&](StepControl& ctrl) -> bool {
    return stepper->doStep(state, ctrl);
  };
  auto doObserve = [&](double t) {
    observer(t, state, model, params);
  };
  
  run(doStep, doObserve, params);
}


void run()
{
  double tend = 30.;
  
  ptree p;
  p.put("fn_out", "stepper_test");
  p.put("tend", tend);
  const int N_INTERVALLS = 8;
  double out_intervall[N_INTERVALLS] = { 0.01, 0.05, 0.1, 0.5, 1., 5. , 10., 15. };
  double lambda = -1./2;
  
  string fn_out = p.get<string>("fn_out")+".h5";
  h5cpp::File f(fn_out.c_str(), "w");

  { // store exact solution
    h5cpp::Group g = f.root().create_group("exact");
    double t = 0., dt = 0.1;
    DynArray<double> ax, ay;
    while (true)
    {
      double y = std::exp(lambda*t);
      ax.push_back(t);
      ay.push_back(y);
      t += dt;
      if (t > tend + 0.5*dt) break;
    }
    h5cpp::create_dataset<>(g, "x", ax);
    h5cpp::create_dataset<>(g, "y", ay);
  }
  f.root().attrs().set("lambda", lambda);

#if 1
  int id = 0;
  for (int i=0; i<N_INTERVALLS; ++i)
  {
    for (int j=0; j<N_METHODS; ++j)
    {
      p.put("method", j);
      p.put("out_intervall", out_intervall[i]);
      run(f.root(), p, id++, lambda);
    }
  }
#else
  // test one method
  p.put("method", Dumka3);
  p.put("out_intervall", 5.);
  run(f.root(), p, 0);
#endif
}


int main(int argc, char **argv)
{
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
  ::run();
}

