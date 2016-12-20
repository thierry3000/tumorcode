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
#ifndef TEST_STEPPER_UTILS_V2_H
#define TEST_STEPPER_UTILS_V2_H

#include "mwlib/time_stepper_new.h"

#include "common.h"
#include "hdfio.h"

using boost::property_tree::ptree;

#if 0
template<class Stepper, class State, class Model, class Observer>
void run(State &state, Model &model, Observer &observer, const ptree &params)
{
  double tend = params.get<double>("tend");
  double out_intervall = params.get<double>("out_intervall");
  double t = 0, dt = 0.;

  cout << format("model run %f steps") % out_intervall << endl;

  Stepper stepper(model, state);

  observer.prepare(model, params);

  int out_num = 0;
  double next_time = 0;
  while (true)
  {
    if (t >= next_time - 0.1 * dt)
    {
      //cout << format("out @t = %f, dt = %f, num = %i") % t % dt % out_num << endl;
      observer.observe(t, out_num, state, model, params);
      // ------
      ++out_num;
      next_time += out_intervall;
    }
    if (stepper.ctrl.stopflag) break;
    if (t >= tend - 0.1 * dt) break;
    dt =  stepper.doStep(next_time - t);
    t += dt;
  }

  observer.postprocess(model, params);
}
#endif

namespace NewSteppers
{


/* Implementes the split step method using the 2nd order accurate strange splitting scheme.
   Given an ODE like
    du/dt = A[u,t] + B[u,t] + C[u,t] + ... 
   with operators A, B, C and so on, the method applies integration steps for A, B, C, etc. individually one after another.
*/
template<class State>
static bool doStrangeSplittingStepsNoBackup(
    State &state, NewSteppers::StepControl &ctrl, int num_operators,
    boost::function<bool (int, State&, NewSteppers::StepControl &)> stepfunc)
{
  typedef NewSteppers::StepControl StepControl;
  ctrl.dt = std::min(0.5*ctrl.dt, 0.90 * ctrl.euler_dt); // individual steps are half of the requested dt
  
  StepControl original_ctrl(ctrl);

  auto exec_step = [&](int idx, State &state, StepControl &ctrl) -> bool
  {
    StepControl c(ctrl);
    bool res = stepfunc(idx, state, c);
    if (!my::compare_approx_less(c.euler_dt, ctrl.euler_dt) ||
        !my::compare_approx_equal(ctrl.dt, c.dt)) res = false;
    // There is a little caveat. If one of the steppers needs a smaller time step than anticipated,
    // the results of previously evaluated steppers is invalid since they used an inconsistently large step size.
    // In case of this we remember the smaller step size and note that we cannot use the results.
    ctrl.euler_dt = c.euler_dt;
    //cout << format("substep %i, t = %f, dt = %f, euler_dt = %f, eval ok = %i, is_rewind = %i") % idx % ctrl.t % c.dt % c.euler_dt % res % ctrl.is_rewind_step << endl;
    return res;
  };

  bool ok = true;
  // strange splitting
  int i = 0;
  for (; i<num_operators; ++i)
  {
    ok &= exec_step(i, state, ctrl);
  }
  ctrl.update();
  myAssert(ok == false || my::compare_approx_equal(ctrl.t, original_ctrl.t + ctrl.dt));
  for (--i; i>=0; --i)
  {
    ok &= exec_step(i, state, ctrl);
  }
  ctrl.update();
  myAssert(ok == false || my::compare_approx_equal(ctrl.t, original_ctrl.t + 2.*ctrl.dt));
  ctrl.dt *= 2.;
  return ok;
}

/* This is a variant of the split step method which can handle mistakenly overestimated step sizes.
   It does so by keeping a copy of the initial state as backup and rewinds if needed, using the smallest
   valid step size reported by the individual steppers. */
template<class State>
static bool doStrangeSplittingSteps(
    State &state, NewSteppers::StepControl &ctrl, int num_operators,
    boost::function<bool (int, State&, NewSteppers::StepControl &)> stepfunc,
    boost::function<void (State &, const State &)> cloneState)
{
  State backup; cloneState(backup, state);
  while(true)
  {
    StepControl original_ctrl(ctrl);
    bool ok = doStrangeSplittingStepsNoBackup(state, ctrl, num_operators, stepfunc);
    if (ok)
      break;
    else
    {
      // basically we keep the euler_dt and try again in the old state
      ctrl.t = original_ctrl.t;
      ctrl.dt = original_ctrl.dt;
      ctrl.is_rewind_step = true; // also set this so the stepper knows that it has to discard state histories
      cout << format("SplitStep Rewind: requested step dt/2=%f, but euler_dt=%f") % ctrl.dt % ctrl.euler_dt << endl;
      cloneState(state, backup);
    }      
  }
  return true;
}



template<class DoStep, class Observer>
static bool run(DoStep &doStep, Observer &observer, const ptree &params)
{
  typedef NewSteppers::StepControl StepControl;
  
  double tend = params.get<double>("tend");
  double out_intervall = params.get<double>("out_intervall");

  StepControl ctrl;
  double next_time = 0;
  bool cont = true;
  while (cont)
  {
    if (ctrl.t >= next_time - 0.1 * ctrl.dt)
    {
      //cout << format("out @t = %f, dt = %f, num = %i") % t % dt % out_num << endl;
      observer(ctrl.t);
      next_time += out_intervall;
    }
    if (ctrl.t >= tend - 0.1 * ctrl.dt) break;
    ctrl.dt = next_time - ctrl.t;
    cont = doStep(ctrl);
  }
  return cont;
}




/*
 * parameters:
 *  fn_out : string
 *  save_images : bool
 *  save_hdf : bool
 *  hdf_clear_file : bool, optional
 *  hdf_group_name_pattern : string, e.g. "out%04i", optional
 *  stepper : string, optional
 */
struct ObserverPde
{
  string fn_out, group_name_pattern;
  bool save_hdf;
  h5cpp::File f;
  my::Time t_real_start;
  bool first, clear_the_file;
  int out_num;

  h5cpp::File openH5File()
  {
    if (first)
    {
      first = false;
      return h5cpp::File(fn_out+".h5", clear_the_file ? "w" : "a");
    }
    else
      return h5cpp::File(fn_out+".h5", "a");
  }

  ObserverPde(const ptree &params)
  {
    fn_out = params.get<string>("fn_out");
    save_hdf = params.get<bool>("save_hdf");
    first = true;
    clear_the_file = params.get<bool>("hdf_clear_file", true);
    group_name_pattern = params.get<string>("hdf_group_name_pattern", "out%04i");
    out_num = 0;
    
    boost::optional<string> stepper = params.get_optional<string>("stepper");
    if (stepper)
    {
      h5cpp::File f = openH5File();
      f.root().attrs().set<string>("stepper", *stepper);
    }
  }

  template<class State, class Model>
  void writeH5(const string &groupname, double t, State &state, Model &model)
  {
    h5cpp::File f = openH5File();
    cout << format("hdf output -> %s") % f.get_file_name() << endl;
    h5cpp::Group g = f.root().create_group(groupname);
    g.attrs().set("time", t);
    g.attrs().set("real_time", (my::Time() - t_real_start).to_s());
    MemUsage memusage = GetMemoryUsage();
    g.attrs().set<uint64>("mem_vsize", memusage.vmem_peak);
    g.attrs().set<uint64>("mem_rss", memusage.rss_peak);
    model.writeH5(f, g, state, t, out_num);
  }

  template<class State, class Model>
  void operator()(double t, State &state, Model &model, const ptree &params)
  {
    if (save_hdf)
    {
      writeH5((format("out%04i") % out_num).str(), t, state, model);
    }
    ++out_num;
  }
};


struct ObserverOde
{
  h5cpp::Group g;
  DynArray<DynArray<double> > data;
  DynArray<string> data_names;
  typedef void result_type;

  template<class Model>
  ObserverOde(Model &model, const ptree &params)
  {
    data.clear();
    data_names.clear();
    model.getDataNames(data_names);
    data.resize(data_names.size());
  }
  
  template<class State, class Model>
  void operator()(double t, State &state, Model &model, const ptree &params)
  {
    DynArray<double> d(data_names.size()-1);

    data[0].push_back(t);

    model.storeState(t, state, d);
    for (int i=0; i<data_names.size()-1; ++i)
    {
      data[i+1].push_back(d[i]);
    }
  }

  void finish()
  {
    for (int i=0; i<data_names.size(); ++i)
    {
      h5cpp::Dataset ds = h5cpp::create_dataset(g, data_names[i], data[i]);
    }
  }
};


}

#endif
