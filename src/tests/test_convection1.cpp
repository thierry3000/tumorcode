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
#include <boost/array.hpp>
#include <boost/property_tree/info_parser.hpp>

#include "mwlib/ptree_ext.h"

#include "common.h"
#include "time_stepper_utils.h"
#include "continuum-flow.h"




namespace ConvectionModel
{


enum {
  VFILL_CONST = 1,
  VFILL_DIVERGE_CUBIC = 2
};


struct Model
{
  typedef Array3d<float> State;
  enum { stepperFlags = 0 };
  
  void init(const LatticeDataQuad3d &grid_, int ndims_,  const ptree &params)
  {
    grid = grid_;
    ndims = ndims_;
    LatticeIndexRanges ir = LatticeIndexRangesFromCellRange(grid.Box(), ndims);
    for (int i=0; i<ndims; ++i)
    {
      velocity_fields[i].initFromBox(ir.faces[i]);
    }
    use_conservative_term = params.get<bool>("conservative");
    // fill field here
    int mode = params.get<int>("velocity_profile");
    if (mode == VFILL_CONST)
    {
      Float3 v = params.get<Float3>("v");
      for (int i=0; i<ndims; ++i)
      {
        FOR_BBOX3(p, ir.faces[i])
        {
          velocity_fields[i](p) = v[i];
        }
      }
    }
    else if (mode == VFILL_DIVERGE_CUBIC)
    {
      LatticeDataQuad3d ld_face(ir.faces[0], grid.Scale());
      ld_face.SetCellCentering(Vec<bool,3>(false, true, true));
      ld_face.SetWorldPosition(grid.GetWorldBox().min);
      double v = params.get<float>("v");
      FOR_BBOX3(p, ir.faces[0])
      {
        Float3 wp = ld_face.LatticeToWorld(p);
        velocity_fields[0](p) = v * 2.0 * (my::smooth_heaviside<float>(wp[0], ::Size(grid.GetWorldBox())[0]*0.2) - 0.5); //wp[0]*wp[0]*my::sign(wp[0]);
      }
    }

    float max_vel = 0.;
    for (int i=0;i<ndims; ++i)
    {
      FOR_BBOX3(p, ir.faces[i])
      {
        max_vel = std::max(std::abs(velocity_fields[i](p)), max_vel);
      }
    }
    max_dt = 0.5 * grid.Scale() / max_vel / (2.*ndims);
    cout << "init done" << endl;
    cout << "grid = "; grid.print(cout); cout<<endl;
    cout << "dt = " << max_dt << endl;
  }

  void calcSlope(const StateConc &conc, StateConc &slope, StepControl &ctrl)
  {
    if (slope.empty()) {
      slope.initFromBox(ExtendForDim(grid.Box(), ndims, 2));
    }
    slope.fill(0);
    ctrl.dt = std::min<double>(ctrl.dt, max_dt);
    LatticeIndexRanges ir = LatticeIndexRangesFromCellRange(grid.Box(), ndims);
    
    CopyBorder(const_cast<StateConc&>(conc)[ir.cells], ndims, 2);
    //AddFctStepSlope<float>(ir.cells, state.conc, velocity_fields.data(), ndims, ctrl.dt, grid.Scale(), slope.conc, use_conservative_term ? CONSERVATIVE : NON_CONSERVATIVE);
    AddFctSlope2<float>(ir.cells, conc, velocity_fields.data(), ndims, grid.Scale(), slope, use_conservative_term ? CONSERVATIVE : NON_CONSERVATIVE);
  }

  void writeH5(h5cpp::File f, h5cpp::Group g, StateConc &conc, double t, int out_num)
  {
    h5cpp::Group ld_group = RequireLatticeDataGroup(f.root(), "field_ld", grid);
    h5cpp::Dataset ds = WriteScalarField(g, "conc", conc[grid.Box()], grid, ld_group);
  }

  void appendToImages(StateConc &conc, std::vector<Image> &images)
  {
    Image img;
    DrawArray(img, conc[grid.Box()], DrawArrayOpts().outputRange().scaled(0.5, 0.25));
    images.push_back(img);
    StateConc slope;
    StepControl ctrl;
    calcSlope(conc, slope, ctrl);
    DrawArray(img, slope[grid.Box()], DrawArrayOpts().outputRange());
    images.push_back(img);
  }

  bool use_conservative_term;
  double max_dt;
  int ndims;
  LatticeDataQuad3d grid;
  boost::array<Array3d<float>, 3> velocity_fields;
};


void run(const ptree &params)
{
  Model model;
  StateConc state;
  ObserverPde observer;
  
  LatticeDataQuad3d ld;
  int ndims;

  float scale = params.get<float>("scale");
  string init_conc_from = params.get<string>("init_conc");
  if (GetExtension(init_conc_from) == "png")
  {
    string in_fn = init_conc_from;
    Array3d<float> a = ReadArrayFromImage<float>(in_fn,0,1);
    ld.Init(a.size(), scale);
    
    ld.SetOriginPosition(-0.5*ld.GetWorldBox().max);
//     ld.SetWorldPosition(-ld.GetWorldBox().max * 0.5); // set origin = lower left side
    ndims = a.size()[1] > 1 ? 2 : 1;
    model.init(ld, ndims, params);
    state.initFromBox(ExtendForDim(ld.Box(), ndims, 2));
    state[ld.Box()].fill(a);
  }
  if (GetExtension(init_conc_from) == "h5")
  {
    
  }
  else throw std::runtime_error(str(format("dont know what to do with input %s") % init_conc_from));
  
  Steppers::Stepper<Model>::BasePointer stepper = Steppers::Stepper<Model>::create(params.get<string>("stepper"));
  ::run<>(state, model, observer, *stepper, params);
}


}


ptree make_default_params()
{
  ptree p;
  p.put("scale", 10.);
  p.put("size", "set by file or <x,y,z>");
  p.put("init_conc", "file.png or file.h5[:path] or hardcoded type");
  p.put("conservative", true);
  p.put("velocity_profile", "file.h5[:path] or hardcoded type");
  p.put("velocity_vector", Float3(10., 0, 0));
  p.put("fn_out", "out");
  p.put("tend", 10);
  p.put("out_intervall", 1);
  p.put("save_hdf", true);
  p.put("save_images", false);
  p.put("stepper", "rk3");
  p.put("num_threads", 1);
  return p;
}


int main(int argc, char **argv)
{
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
  my::MultiprocessingInitializer mpinit(argc, argv);

  my::ProgramOptions list;
  list.show_help_on();
  string param_fn, fn_out, param_cmdline;
  list.add(param_fn, "-p", "parameter filename (use 'show' to see the default parameters)");
  list.add(param_cmdline, "-pa", "parameter argument string");
  list.add(fn_out, "-o", "output filename");

  if(!list.parse(argv, argc)) return -1;

  // default parameters
  ptree params = make_default_params();

  if (param_fn == "show")
  {
    cout << "parameters:" << endl;
    boost::property_tree::write_info(cout, params);
    return 0;
  }

  if (!param_fn.empty()) {
    std::ifstream f(param_fn.c_str());
    ptree p;
    boost::property_tree::read_info(f, p);
    boost::property_tree::update(params, p);
  }
  if (!param_cmdline.empty()) {
    std::istringstream f(param_cmdline);
    ptree p;
    boost::property_tree::read_info(f, p);
    boost::property_tree::update(params, p);
  }
  my::SetNumThreads(params.get<int>("num_threads", 1));
  ConvectionModel::run(params);
}

