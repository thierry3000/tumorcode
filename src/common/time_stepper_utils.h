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
#ifndef TEST_STEPPER_UTILS_H
#define TEST_STEPPER_UTILS_H

#include "mwlib/time_stepper.h"

#include "common.h"
#include "hdfio.h"

using boost::property_tree::ptree;
using Steppers::StepControl;

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


namespace Steppers
{

inline bool combine_split_step_ctrl_(StepControl &dst, const StepControl &ctrl)
{
  double ddt = dst.dt - ctrl.dt;
  double euler_dt = dst.euler_dt;
  dst = ctrl;
  dst.euler_dt = std::min(dst.euler_dt, euler_dt);
  return std::abs(ddt)/(std::abs(ctrl.dt) + 1.e-18) < 1.e-6;
}
  
}
/**
 * need to match the time step for the drug
 * propergation with the time steps of ??
 */
template<class State>
static StepControl doStrangeSplittingSteps(State &state, const StepControl &ctrl_, int num_operators, boost::function3<StepControl, int, State&, const StepControl &> stepfunc)
{
  State backup; backup.initCloned(state);

  StepControl ctrl = ctrl_;
  // individual steps are half of the requested dt
  ctrl.dt = std::min(0.5*ctrl_.dt, 0.90 * ctrl.euler_dt); 
  while (true)
  {
    // strange splitting
    StepControl combined_ctrl; 
    combined_ctrl.dt = ctrl.dt;
    int i = 0;
    bool ok = true;
    //we are doing this step for several operators
    for (; i<num_operators; ++i)
      ok &= Steppers::combine_split_step_ctrl_(combined_ctrl, stepfunc(i, state, ctrl));
    ctrl.t = combined_ctrl.t; // the time was advanced
    for (--i; i>=0; --i)
      ok &= Steppers::combine_split_step_ctrl_(combined_ctrl, stepfunc(i, state, ctrl));
    // rewind
    // the steppers covered less time than intended and therefore 
    // it has to be assumed that they did steps of differing lengths
    if (!ok) 
    {
      cout << format("SplitStep Rewind: requested step dt/2=%f, but euler_dt=%f") % ctrl.dt % combined_ctrl.euler_dt << endl;
      state.addScaled(1., backup, 0.);
      state.setId(backup.getId());
      // try again
      ctrl.t = ctrl_.t;
      ctrl.dt = std::min(ctrl.dt, 0.98 * combined_ctrl.euler_dt);
      ctrl.is_rewind_step = true;
    }
    else
    {
      ctrl = combined_ctrl;
      break;
    }
  }
  ctrl.dt *= 2.;
  return ctrl;
}


template<class State, class Model, class Observer>
static void run(State &state, Model &model, Observer &observer, Steppers::Stepper<Model> &stepper, const ptree &params)
{
  double tend = params.get<double>("tend");
  double out_intervall = params.get<double>("out_intervall");

  stepper.init(model);

  StepControl ctrl;
  double next_time = 0;
  while (true)
  {
    if (ctrl.t >= next_time - 0.1 * ctrl.dt)
    {
      //cout << format("out @t = %f, dt = %f, num = %i") % t % dt % out_num << endl;
      observer(ctrl.t, state, model, params);
      next_time += out_intervall;
    }
    if (ctrl.stopflag) break;
    if (ctrl.t >= tend - 0.1 * ctrl.dt) break;
    ctrl.dt = next_time - ctrl.t;
    ctrl = stepper.doStep(state, ctrl);
  }
}


template<class State, class Model, class Observer>
static void run_model(State &state, Model &model, Observer &observer, const ptree &params)
{
  double tend = params.get<double>("tend");
  double out_intervall = params.get<double>("out_intervall");

  //cout << format("model run %f steps") % out_intervall << endl;

  Steppers::StepControl ctrl;
  
  double next_time = 0;
  while (true)
  {
    if (ctrl.t >= next_time - 0.1 * ctrl.dt)
    {
      //cout << format("out @t = %f, dt = %f, num = %i") % t % dt % out_num << endl;
      observer(ctrl.t, state, model, params);
      next_time += out_intervall;
    }
    if (ctrl.stopflag) break;
    if (ctrl.t >= tend - 0.1 * ctrl.dt) break;
    ctrl.dt = next_time - ctrl.t;
    ctrl =  model.doStep(state, ctrl);
  }
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
  bool save_hdf, save_image;
  H5::H5File f;
  my::Time t_real_start;
  bool first, clear_the_file;
  int out_num;

  H5::H5File openH5File()
  {
    if (first)
    {
      first = false;
      return H5::H5File(fn_out+".h5", clear_the_file ? H5F_ACC_RDWR : H5F_ACC_TRUNC);
    }
    else
      return H5::H5File(fn_out+".h5", H5F_ACC_TRUNC);
  }

  ObserverPde(const ptree &params)
  {
    fn_out = params.get<string>("fn_out");
    save_image = params.get<bool>("save_images");
    save_hdf = params.get<bool>("save_hdf");
    first = true;
    clear_the_file = params.get<bool>("hdf_clear_file", true);
    group_name_pattern = params.get<string>("hdf_group_name_pattern", "out%04i");
    out_num = 0;
    
    boost::optional<string> stepper = params.get_optional<string>("stepper");
    if (stepper)
    {
      H5::H5File f = openH5File();
      H5::Group root = f.openGroup("/");
      writeAttrToH5(root, string("stepper"), *stepper);
      //f.root().attrs().set<string>("stepper", *stepper);
    }
  }

  template<class State, class Model>
  void writeH5(const string &groupname, double t, State &state, Model &model)
  {
    H5::H5File f = openH5File();
    cout << format("hdf output -> %s") % f.getFileName() << endl;
    H5::Group g = f.createGroup(groupname);
    writeAttrToH5(g, string("time"), t);
    writeAttrToH5(g, string("real_time"), (my::Time() - t_real_start).to_s());
//     g.attrs().set("time", t);
//     g.attrs().set("real_time", (my::Time() - t_real_start).to_s());
    MemUsage memusage = GetMemoryUsage();
    writeAttrToH5(g, string("mem_vsize"), (int)memusage.vmem_peak);
    writeAttrToH5(g, string("mem_rss"), (int)memusage.rss_peak);
//     g.attrs().set<uint64>("mem_vsize", memusage.vmem_peak);
//     g.attrs().set<uint64>("mem_rss", memusage.rss_peak);
    model.writeH5(f, g, state, t, out_num);
  }

  template<class State, class Model>
  void operator()(double t, State &state, Model &model, const ptree &params)
  {
    if (save_hdf)
    {
      writeH5((format("out%04i") % out_num).str(), t, state, model);
    }
    if (save_image)
    {
      string fn = str(format("%s_%04i.png") % fn_out % out_num);
      cout << format("image output -> %s") % fn << endl;
      std::vector<Image> images;
      model.appendToImages(state, images);
      Image bigimg;
      DrawImageGrid(bigimg, images);
      bigimg.Write(fn);
    }
    ++out_num;
  }
};


struct ObserverOde
{
  H5::Group g;
  DynArray<DynArray<double> > data;
  DynArray<string> data_names;

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
      //h5cpp::Dataset ds = h5cpp::create_dataset(g, data_names[i], data[i]);
      H5::DataSet ds = writeDataSetToGroup(g, data_names[i], data[i]);
    }
  }
};


#endif