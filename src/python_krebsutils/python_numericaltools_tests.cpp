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
#include "python_helpers.h"

#include "common.h"
#include "continuum-flow.h"
#include "continuum-surfacetension.h"
#include "lattice-data-polymorphic.h"
#include "levelset.h"
#include "shared-objects.h"
#include "convection_diffusion_solver.h"
#include "continuum-stepper-states.h"
#include "time_stepper_utils_new.h"


using NewSteppers::ObserverPde;
using NewSteppers::ObserverOde;
using NewSteppers::StepControl;


template<class T>
void AddAdvection(const BBox3 &bb, const ConstArray3d<T> &val, const ConstArray3d< T >*const u, int ndims, T spacing, T dt, Array3d<T> res, const ptree &params)
{
  string advection_method = params.get<string>("advection_method");
  ConvectionMode conservation = params.get<bool>("conservative",true) ? CONSERVATIVE : NON_CONSERVATIVE;
  if (advection_method == "kt")
  {
    AddKTAdvection<float>(bb, val, u, ndims, spacing, res, conservation);
  }
#if 0
  // thsi scheme is nonsense i think
  else if(advection_method == "fct")
  {
    AddFctStepSlope<float>(bb, val, u, ndims, dt, spacing, res, conservation);
  }
#endif
  else if(advection_method == "upwind1")
  {
    AddENOAdvection(bb, val, u, ndims, spacing, res,  conservation, 1);
  }
  else if(advection_method == "eno2")
  {
    AddENOAdvection(bb, val, u, ndims, spacing, res,  conservation, 2);
  }
  else if(advection_method == "eno3")
  {
    AddENOAdvection(bb, val, u, ndims, spacing, res,  conservation, 3);
  }
  else if(advection_method == "weno5")
  {
    AddWENOAdvection(bb, val, u, ndims, spacing, res,  conservation, 5);
  }
  else if(advection_method == "weno3")
  {
    AddWENOAdvection(bb, val, u, ndims, spacing, res,  conservation, 3);
  }
  else throw std::runtime_error(str(format("advection_method %s unkown!") % advection_method));
}


inline float closest_intersect(float a, float b)
{
  return std::min(a,b);
}


float zalezakDiskFunction(const Float3 &pp, float radius, const Float3 &pos)
{
  // in scale*[0,1]^2 domain
  Float3 p = (pp - pos) / radius;
  
  //const float rad = 15. / 100.;
  //const Float3 pos(50./100.,75./100.,0);
  const float rad = 1.;
  const float width = 0.33;
  const float len = 1.;
  const float len_ext = 0.5;
  const Float3 slot_pos(0, 0 - rad + len/2. - len_ext, 0.);

  float d1 = rad - p.norm();
  float d2 = -my::min(width/2.-std::abs(p[0] - slot_pos[0]), (len+len_ext)/2.-std::abs(p[1]-slot_pos[1]), width/2.-std::abs(p[2]-slot_pos[2]));
  return radius * closest_intersect(d1,d2);
}


Float3 rotationalVelocityField(const Float3 &wp, const Float3 &rot_axis, const Float3 &center, bool shear)
{
  if (shear) {
    return cross<float>(rot_axis, safe_normalized(wp - center));
  } else
    return cross<float>(rot_axis, wp - center);
}



void generateVelocityField(const LatticeDataQuad3d &ld, int ndims, const ptree &params, const py::dict &kwargs, FaceVarArrays &velocities)
{
  LatticeIndexRanges ir = LatticeIndexRanges::FromCellRange(ld.Box(), ndims);

  py::object the_func = kwargs["velocity_profile"];
  for (int axis=0; axis<ndims; ++axis)
  {
    LatticeDataQuad3d ld_face = CellLdToFaceLd(ld, ndims, axis);
    FOR_BBOX3(p, ir.faces[axis])
    {
      Float3 wp = ld_face.LatticeToWorld(p);
      Float3 v = py::extract<Float3>(the_func(wp[0], wp[1], wp[2])); // zuviel gerechnet aber egal ...
      velocities[axis](p) = v[axis];
    }
  }
}


void generateInitialConc(const LatticeDataQuad3d &ld, int ndims, const ptree &params, const py::dict &kwargs, Array3d<float> conc)
{
  py::object the_func = kwargs["conc_profile"];
  FOR_BBOX3(p, ld.Box())
  {
    Float3 wp = ld.LatticeToWorld(p);
    float v = py::extract<float>(the_func(wp[0], wp[1], wp[2]));
    conc(p) = v;
  }
}



namespace ConvectionModelTests
{

enum {
  VFILL_CONST = 1,
  VFILL_DIVERGE_CUBIC = 2
};


typedef Array3d<float> State;
typedef Array3dOps<float> Ops;

struct Model
{
  typedef ConvectionModelTests::State State;
  enum { stepperFlags = 0 };

  void init(const LatticeDataQuad3d &grid_, int ndims_, Ops &ops_, const ptree &params_, FaceVarArrays &velocities)
  {
    grid = grid_;
    ndims = ndims_;
    params = params_;
    ops = &ops_;

    
    LatticeIndexRanges ir = LatticeIndexRanges::FromCellRange(grid.Box(), ndims);
    float max_vel = 0.;
    for (int i=0;i<ndims; ++i)
    {
      velocity_fields[i] = velocities[i];
      FOR_BBOX3(p, ir.faces[i])
      {
        max_vel = std::max(std::abs(velocity_fields[i](p)), max_vel);
      }
    }

    max_dt = 0.5 * grid.Scale() / max_vel / (2.*ndims);
    cout << "dt = " << max_dt << endl;
  }

  void calcSlope(const State &conc, State &slope, StepControl &ctrl)
  {
    ops->init(slope, false);
    ctrl.euler_dt.min(max_dt);
    CopyBorder(const_cast<State&>(conc)[grid.Box()], ndims, 3);
    AddAdvection<float>(grid.Box(), conc, velocity_fields.data(), ndims, grid.Scale(), std::min(max_dt, ctrl.dt), slope, params);
  }

  void writeH5(H5::H5File f, H5::Group g, State &conc, double t, int out_num)
  {
    H5::Group ld_group = RequireLatticeDataGroup(f, "field_ld", grid);
    H5::DataSet ds = WriteScalarField(g, "conc", conc[grid.Box()], grid, ld_group);
  }

  void appendToImages(State &conc, std::vector<Image> &images)
  {
    Image img;
    DrawArray(img, conc[grid.Box()], DrawArrayOpts().outputRange().scaled(0.5, 0.25));
    images.push_back(img);
    State slope;
    StepControl ctrl;
    calcSlope(conc, slope, ctrl);
    DrawArray(img, slope[grid.Box()], DrawArrayOpts().outputRange());
    images.push_back(img);
  }

  ptree params;
  double max_dt;
  int ndims;
  LatticeDataQuad3d grid;
  FaceVarArrays velocity_fields;
  Ops *ops;
};


ptree make_default_params()
{
  ptree p;
  p.put("scale", 10.);
  p.put("size", "set by file or <x,y,z>");
  p.put("conservative", true);
  p.put("conc_profile", "func");
  p.put("velocity_profile", "func");
  p.put("fn_out", "out");
  p.put("tend", 10);
  p.put("out_intervall", 1);
  p.put("save_hdf", true);
  p.put("save_images", false);
  p.put("stepper", "rk3");
  p.put("num_threads", 1);
  return p;
}


void runTest(const py::str &param_info_str, const py::object &py_ld, const py::dict &kwargs)
{
  ptree params = convertInfoStr(param_info_str, ConvectionModelTests::make_default_params());

  LatticeDataQuad3d ld = py::extract<LatticeDataQuad3d>(py_ld);
  Int3 size = ::Size(ld.Box());
  int ndims = size[1]<=1 ? 1 : (size[2]<=1 ? 2 : 3);

  Model model;
  ObserverPde observer(params);

  Ops ops;
  DomainDecomposition mtboxes;
  mtboxes.insert(0,ld.Box()); // single threaded
  ops.init(mtboxes, ld.Box(), ndims, 3);
  
  State state;
  ops.init(state);
  generateInitialConc(ld, ndims, params, kwargs, state);

  FaceVarArrays velocities(ld.Box(), ndims);
  generateVelocityField(ld, ndims, params, kwargs, velocities);
  
  model.init(ld, ndims, ops, params, velocities);

  auto stepper = NewSteppers::create<Model*, Ops*, NewSteppers::NO_IMPLICIT>(
    params.get<string>("stepper"), &model, &ops);

  auto doStep = [&](StepControl& ctrl) -> bool {
    ctrl.euler_dt = my::MakeSmallestMax<double>();
    bool ok = stepper->doStep(state, ctrl);
    return ok && !ctrl.is_rewind_step;
  };
  auto doObserve = [&](double t) {
    observer(t, state, model, params);
  };
  NewSteppers::run(doStep, doObserve, params);
}


}  // ConvectionModelTests



void runTestRedistancing(const py::str &param_info_str, const py::object &py_ld, const py::dict &kwargs)
{
  ptree params = convertInfoStr(param_info_str, ConvectionModelTests::make_default_params());

  LatticeDataQuad3d ld = py::extract<LatticeDataQuad3d>(py_ld);
  Int3 size = ::Size(ld.Box());
  int ndims = size[1]<=1 ? 1 : (size[2]<=1 ? 2 : 3);

  Levelset ls; ls.init(ld, ndims, 0);
  generateInitialConc(ld, ndims, params, kwargs, ls.phi[ld.Box()]);

  reinit_levelset(ls, ld, ndims, 20, params);
}



namespace LevelsetConvection
{

typedef Levelset State;
typedef LevelsetOps Ops;


class Model
{
public:
  typedef Levelset State;
  
  void init(const LatticeDataQuad3d &grid_, int ndims_, const ptree &params_, FaceVarArrays &velocities)
  {
    grid = grid_;
    ndims = ndims_;
    params = params_;
    params.put("conservative", false);
    mtboxes.init(MakeMtBoxGrid(grid.Box(), Int3(64, 32, 32)));
    ops.init(mtboxes, grid.Box(), ndims, 3);
    stepper = NewSteppers::create<Model*,Ops*,NewSteppers::NO_IMPLICIT>(
      params.get<string>("stepper"), this, &ops);

    LatticeIndexRanges ir = LatticeIndexRanges::FromCellRange(grid.Box(), ndims);
    float max_vel = 0.;
    for (int i=0;i<ndims; ++i)
    {
      velocity_fields[i] = velocities[i];
      FOR_BBOX3(p, ir.faces[i])
      {
        max_vel = std::max(std::abs(velocity_fields[i](p)), max_vel);
      }
    }

    max_dt = 0.8 * grid.Scale() / max_vel / (ndims);
    last_levelset_reinit = 0;
    cout << "dt = " << max_dt << endl;
  }

  bool doStep(State &levelset, StepControl ctrl)
  {
    if (params.get<bool>("levelset_reinit", true))
    {
      //if (last_levelset_reinit < ctrl.t - max_dt * 20)
      {
        //double gnorm = calc_levelset_badness(state.levelset, grid, ndims);
        //cout << format("gnorm = %f, maxval = %f") % gnorm % state.levelset.max_val << endl;
        //if (gnorm > 0.0001)
        if (last_levelset_reinit < ctrl.t - max_dt)
        {
          string method = params.get<string>("gradient_approx_method");
          if (method != "none")
          {
            ptree reinit_pt = make_ptree("gradient_approx_method", method);
            ptree pt = reinit_levelset(levelset, grid, ndims, 5, reinit_pt);
            cout << "levelset reinit" << endl;
          }
          last_levelset_reinit = ctrl.t;
        }
      }
    }

    {
      double g = calc_levelset_badness(levelset, grid, ndims);
      gradientnorm.push_back(g);
      time.push_back(ctrl.t);
    }
    
    cout << format("step t=%f") % ctrl.t << endl;
    return stepper->doStep(levelset, ctrl);
  }

  void calcSlope(const State &levelset, State &slope, StepControl &ctrl)
  {
    ctrl.euler_dt.min(max_dt);
    ops.init(slope, true);

    #pragma omp parallel
    {
      BOOST_FOREACH(const BBox3 &bb, mtboxes.getCurrentThreadRange())
      {
        CopyBorder(const_cast<State&>(levelset).phi, bb, grid.Box(), ndims, 3);
        AddAdvection<float>(bb, levelset.phi, velocity_fields.data(), ndims, grid.Scale(), std::min(max_dt, ctrl.dt), slope.phi, params);
      }
    }
    // debug
    #if 0
    std::vector<Image> images;
    images.push_back(DrawArray(facevalues[0], DrawArrayOpts().title("ux").outputRange()));
    images.push_back(DrawArray(facevalues[1], DrawArrayOpts().title("uy").outputRange()));
    images.push_back(DrawArray(state.levelset.phi, DrawArrayOpts().title("u").outputRange()));
    images.push_back(DrawArray(slope.levelset.phi, DrawArrayOpts().title("du").outputRange()));
    Image img;
    DrawImageGrid(img, images);
    img.Write(str(format("eno%010i.png") % int(ctrl.t*1000000)));
    #endif
  }

  void writeH5(H5::H5File f, H5::Group g, State &state, double t, int out_num)
  {
    H5::Group ld_group = RequireLatticeDataGroup(f, "field_ld", grid);
    WriteScalarField(g, "ls", state.phi, grid, ld_group);
    Array3d<float> theta; theta.initFromBox(grid.Box());
    state.ToHeaviside(theta, grid.Box(), grid.Scale());
    WriteScalarField(g, "theta", theta, grid, ld_group);
    float mass = theta.valueStatistics().Sum() * std::pow(grid.Scale(), ndims);
    writeAttrToGroup<float>(g, string("mass"), mass);
    //g.attrs().set("mass", mass);
  }

  void appendToImages(State &state, std::vector<Image> &images) const
  {
    Image img;
    DrawArray(img, state.phi[grid.Box()], DrawArrayOpts().outputRange());
    images.push_back(img);
    Array3d<float> theta; theta.initFromBox(grid.Box());
    state.ToHeaviside(theta, grid.Box(), grid.Scale());
    DrawArray(img, theta[grid.Box()], DrawArrayOpts());
    images.push_back(img);
  }

  double max_dt;
  double last_levelset_reinit;
  int ndims;
  ptree params;
  LatticeDataQuad3d grid;
  boost::array<Array3d<float>, 3> velocity_fields;
  DynArray<double> gradientnorm;
  DynArray<double> time;
  Ops ops;
  DomainDecomposition mtboxes;
  typedef NewSteppers::Stepper<Model*, Ops*> MyStepper;
  std::auto_ptr<MyStepper> stepper;
};


void runLS(const py::str &param_info_str, const py::object &py_ld, const py::dict &kwargs)
{
  ptree params = convertInfoStr(param_info_str, ConvectionModelTests::make_default_params());

  LatticeDataQuad3d ld = py::extract<LatticeDataQuad3d>(py_ld);
  Int3 size = ::Size(ld.Box());
  int ndims = size[1]<=1 ? 1 : (size[2]<=1 ? 2 : 3);
  
  Model model;
  ObserverPde observer(params);

  FaceVarArrays velocities(ld.Box(), ndims);
  generateVelocityField(ld, ndims, params, kwargs, velocities);

  model.init(ld, ndims, params, velocities);
  
  State state;
  model.ops.init(state, true);
  generateInitialConc(ld, ndims, params, kwargs, state.phi);
  state.max_val = state.phi.valueStatistics().MaxAbs();

  auto doStep = [&](StepControl& ctrl) -> bool {
    ctrl.euler_dt = my::MakeSmallestMax<double>();
    return model.doStep(state, ctrl);
  };
  auto doObserve = [&](double t) {
    observer(t, state, model, params);
  };
  NewSteppers::run(doStep, doObserve, params);

  {
    H5::H5File f = observer.openH5File();
    H5::Group g = f.createGroup("gradientmeasure");
    writeDataSetToGroup<DynArray<double>>(g, string("time"), model.time);
    //H5::createDataset(g, "time", model.time);
    //h5cpp::create_dataset(g, "gradnorm", model.gradientnorm);
    writeDataSetToGroup<DynArray<double>>(g, string("gradientmeasure"), model.gradientnorm);
  }
}


} // LevelsetConvection


#if 0
namespace ConvectionDiffusionReaction
{


struct Testmodel
{
  ContinuumGrid grid_storage;
  ContinuumGrid* grid;
  double time;
  ConvectionReactionDiffusionModel<float> crdmodel;
  py::object py_velocity, py_src_linear, py_src_const, py_src_explicit, py_kdiff;
  Array3d<float> conc;

  void run_me(const py::str &param_info_str, const py::object &py_ld, const py::dict &kwargs)
  {
    ptree params = convertInfoStr(param_info_str, ConvectionModelTests::make_default_params());
    //LatticeDataQuad3d ld = dynamic_cast<const polymorphic_latticedata::Derived<LatticeDataQuad3d>*>(extractLatticeData(py_ld).get())->get();
    LatticeDataQuad3d ld = py::extract<LatticeDataQuad3d>(py_ld);
    Int3 size = ::Size(ld.Box());
    int ndims = size[1]<=1 ? 1 : (size[2]<=1 ? 2 : 3);
    grid_storage = ContinuumGrid(ld, ndims);
    grid = &grid_storage;

    time = 0;
    
    conc = MakeArray3dWithBorder<float>(grid->ld.Box(), grid->dim, 3);
    py::object py_init_conc = kwargs["conc_profile"];
    fillCellVariableByPyFunction(grid->ld.Box(), time, conc, py_init_conc, conc);

    py_velocity = kwargs["velocity"];
    py_src_linear = kwargs["src_impl_clinear"];
    py_src_const = kwargs["src_impl_cconst"];
    py_src_explicit = kwargs["src_expl"];
    py_kdiff = kwargs["kdiff"];

    crdmodel.init(*grid, ptree());
    crdmodel.velocity_f        = boost::bind(&Testmodel::velocity_f, this, _1, _2, _3, _4, _5, _6);
    crdmodel.source_explicit_f = boost::bind(&Testmodel::source_explicit_f, this, _1, _2, _3, _4, _5);
    crdmodel.source_implicit_f = boost::bind(&Testmodel::source_implicit_f, this, _1, _2, _3, _4, _5);
    crdmodel.diffusivity_f     = boost::bind(&Testmodel::diffusivity_f, this, _1, _2, _3, _4);

    ObserverPde observer(params);

    my::SetNumThreads(params.get<int>("num_threads"));
    
    run_model<Array3d<float>, Testmodel, ObserverPde>(conc, *this, observer, params);
  }


  Steppers::StepControl doStep(Array3d<float> state, Steppers::StepControl &ctrl)
  {
    return crdmodel.doStep(state, ctrl);
  }
  

  void writeH5(h5cpp::File f, h5cpp::Group g, Array3d<float> &conc, double t, int out_num) const
  {
    h5cpp::Group ld_group = RequireLatticeDataGroup(f.root(), "field_ld", grid->ld);
    h5cpp::Dataset ds = WriteScalarField(g, "conc", conc[grid->ld.Box()], grid->ld, ld_group);
  }

  void appendToImages(Array3d<float> &state, std::vector<Image> &images) const
  {
  }
  
  void fillCellVariableByPyFunction(const BBox3 &bbox, double time, ConstArray3d<float> conc, py::object py_func, Array3d<float> dst)
  {
    FOR_BBOX3(p, bbox)
    {
      Float3 wp = grid->ld.LatticeToWorld(p);
      float v = py::extract<float>(py_func(wp[0], wp[1], wp[2], time, conc(p)));
      dst(p) = v;
    }
  }


  void fillFaceVariableByPyFunction(const BBox3 &bbox, double time, ConstArray3d<float> conc, py::object py_func, FaceVarArrays &dst)
  {
    LatticeIndexRanges ir = LatticeIndexRanges::FromCellRange(bbox, grid->dim);
    for (int axis=0; axis<grid->dim; ++axis)
    {
      LatticeDataQuad3d ld_face = CellLdToFaceLd(grid->ld, grid->dim, axis);
      FOR_BBOX3(p, ir.faces[axis])
      {
        Float3 wp = ld_face.LatticeToWorld(p);
        Float3 v = py::extract<Float3>(py_func(wp[0], wp[1], wp[2], time, conc(p))); // zuviel gerechnet aber egal ...
        dst[axis](p) = v[axis];
      }
    }
  }
  
  
  void velocity_f(int boxindex, const BBox3 &bboxext, const BBox3 &bbox, Array3d<float> conc, FaceVarArrays &velocities_dst, float &max_vel)
  {
    fillFaceVariableByPyFunction(bbox, time, conc, py_velocity, velocities_dst);
    for (int d=0; d<grid->dim; ++d)
      max_vel = std::max(max_vel, velocities_dst[d].maxAbs());
  }

  void source_implicit_f(int boxindex, const BBox3 &bbox, ConstArray3d<float> conc, Array3d<float> clinear, Array3d<float> cconst)
  {
    fillCellVariableByPyFunction(bbox, time, conc, py_src_linear, clinear);
    fillCellVariableByPyFunction(bbox, time, conc, py_src_const, cconst);
  }

  void source_explicit_f(int boxindex, const BBox3 &bbox, ConstArray3d<float> conc, Array3d<float> src, float &srcmax)
  {
    fillCellVariableByPyFunction(bbox, time, conc, py_src_explicit, src);
    // numerical differentiation for return in srcmax
    Array3d<float> diff(bbox);
    Array3d<float> dc(bbox); dc[bbox].fill(conc[bbox]);
    float h = 1.e-6;
    dc += h;
    fillCellVariableByPyFunction(bbox, time, dc, py_src_explicit, diff);
    diff[bbox] -= src[bbox];
    diff *= 1./h;
    srcmax = std::max(srcmax, diff.maxAbs());
  }

  void diffusivity_f(int boxindex, const BBox3 &bbox, ConstArray3d<float> conc, Array3d<float> kdiff)
  {
    fillCellVariableByPyFunction(bbox, time, conc, py_kdiff, kdiff);
  }
};
  

  
void runCDR(const py::str &param_info_str, const py::object &py_ld, const py::dict &kwargs)
{
  Testmodel m;
  m.run_me(param_info_str, py_ld, kwargs);
}
  
}
#endif


namespace SimpleSteppersTest
{

typedef double State;
typedef NewSteppers::StepControl StepControl;

class Model : boost::noncopyable
{
public:
  typedef double State;
  enum  { stepperFlags = Steppers::CAN_SOLVE_IMPLICIT };

  Model(double lambda) : lambda(lambda), go_implicit(false)
  {
    max_dt = std::numeric_limits<double>::max(); //0.5/std::abs(lambda);
  }

  
  void calcSlope(const State &state, State &slope, StepControl &ctrl)
  {
    slope = state*lambda;
    ctrl.euler_dt.min(max_dt);
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
    ctrl.euler_dt.min(max_dt);
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


void run(H5::Group g, const ptree &params, int id, double lambda)
{
  double y0 = 1.;

  Model model(lambda);
  State state(y0);
  ObserverOde observer(model, params);

  const int method = params.get<int>("method");
  const double dt = params.get<double>("out_intervall");

  //observer.g = g = g.require_group(str(format("study%02i") % id));
  observer.g = g = g.openGroup(str(format("study%02i") % id));
  writeAttrToGroup<string>(g, string("method"), string(method_name[method]));
  writeAttrToGroup<double>(g, string("dt"), dt);
//   g.attrs().set("method", method_name[method]);
//   g.attrs().set("dt", dt);

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

  typedef NewSteppers::StepperFactory<Model*, NewSteppers::Operations<double> > Factory;
  typedef Factory::StepperType MyStepper;
  
  std::auto_ptr<MyStepper> stepper(Factory::create(name));
  stepper->set_model(&model);

  auto doStep = [&](StepControl& ctrl) -> bool {
    ctrl.euler_dt = my::MakeSmallestMax<double>();
    bool ok = stepper->doStep(state, ctrl);
    if (!ok) {
      ctrl.update(); // no problem, we want to see the solution diverge
      ctrl.is_rewind_step = false;
    }
    return true;
  };
  auto doObserve = [&](double t) {
    observer(t, state, model, params);
  };

  bool ok = NewSteppers::run(doStep, doObserve, params);

  observer.finish();
}


void runSimpleSteppersTest(const std::string &fn_out)
{
  double tend = 30.;

  ptree p;
  p.put("tend", tend);
  const int N_INTERVALLS = 8;
  double out_intervall[N_INTERVALLS] = { 0.01, 0.05, 0.1, 0.5, 1., 5. , 10., 15. };
  double lambda = -1./2;

  H5::H5File f(fn_out.c_str(), H5F_ACC_RDWR);

  { // store exact solution
    H5::Group g = f.createGroup("exact");
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
    writeDataSetToGroup<DynArray<double>>(g, string("x"), ax);
    writeDataSetToGroup<DynArray<double>>(g, string("y"), ay);
//     h5cpp::create_dataset(g, "x", ax);
//     h5cpp::create_dataset(g, "y", ay);
  }
  writeAttrToGroup<double>(f.openGroup("/"), string("lambda"), lambda);
  //f.root().attrs().set("lambda", lambda);

#if 1
  int id = 0;
  for (int i=0; i<N_INTERVALLS; ++i)
  {
    for (int j=0; j<N_METHODS; ++j)
    {
      p.put("method", j);
      p.put("out_intervall", out_intervall[i]);
      run(f.openGroup("/"), p, id++, lambda);
    }
  }
#else
  // test one method
  p.put("method", Dumka3);
  p.put("out_intervall", 5.);
  run(f.root(), p, 0);
#endif
}


  
}


#if BOOST_VERSION>106300
py::object calcCurvature(const py::object py_ld, np::ndarray py_phi, bool has_ghost_boundary, bool return_stf)
{
  int dim = py_phi.get_nd();
  const LatticeDataQuad3d ld = py::extract<LatticeDataQuad3d>(py_ld);
  const BBox3 extbox = ExtendForDim(ld.Box(), dim, 1);
  Array3d<float> phi = Array3dFromPy<float>(py_phi);
  if (!has_ghost_boundary)
  {
    phi.move(ld.Box().min);
    Array3d<float> phi_ex(extbox);
    phi_ex[ld.Box()].fill(phi[ld.Box()]);
    CopyBorder(phi_ex[ld.Box()], dim, 1);
    phi = phi_ex;
  }
  else
    phi.move(extbox.min);
  
  myAssert(phi.size() == Size(extbox));
  SurfaceTensionForce<float> sft; sft.init(ld.Box(), phi, dim, ld.Scale());
//   np::arrayt<float> py_res = np::zeros(dim, Cast<Py_ssize_t>(Size(ld.Box())).data(), np::getItemtype<float>());
  np::ndarray py_res = np::zeros(py::tuple(Cast<Py_ssize_t>(Size(ld.Box())).data()), np::dtype::get_builtin<float>());
  Array3d<float> res = Array3dFromPy<float>(py_res);
  res.move(ld.Box().min);
  FOR_BBOX3(p, ld.Box())
  {
    res(p) = sft.computeCurvature(p);
  }
  if (!return_stf)
    return py_res;

  const Int3 sz_ = Size(ld.Box());
  const Py_ssize_t dims[4] = { dim, sz_[0], sz_[1], sz_[2] };
  //np::arrayt<float> acc_stf = np::zeros(4, dims, np::getItemtype<float>());
  np::ndarray acc_stf = np::zeros(py::make_tuple(dim, sz_[0], sz_[1], sz_[2]), np::dtype::get_builtin<float>());
  FOR_BBOX3(p, ld.Box())
  {
    Float3 f = sft.computeForce(p);
    const Int3 q = p - ld.Box().min;
    for (int axis=0; axis<dim; ++axis)
    {
//       acc_stf(axis, q[0], q[1], q[2]) = f[axis];
      acc_stf[axis][q[0]][q[1]][q[2]] = f[axis];
    }
  }
  return py::make_tuple(py_res, acc_stf);
}
#else
py::object calcCurvature(const py::object py_ld, np::arrayt<float> py_phi, bool has_ghost_boundary, bool return_stf)
{
  int dim = py_phi.rank();
  const LatticeDataQuad3d ld = py::extract<LatticeDataQuad3d>(py_ld);
  const BBox3 extbox = ExtendForDim(ld.Box(), dim, 1);
  Array3d<float> phi = Array3dFromPy<float>(py_phi);
  if (!has_ghost_boundary)
  {
    phi.move(ld.Box().min);
    Array3d<float> phi_ex(extbox);
    phi_ex[ld.Box()].fill(phi[ld.Box()]);
    CopyBorder(phi_ex[ld.Box()], dim, 1);
    phi = phi_ex;
  }
  else
    phi.move(extbox.min);
  
  myAssert(phi.size() == Size(extbox));
  SurfaceTensionForce<float> sft; sft.init(ld.Box(), phi, dim, ld.Scale());
  np::arrayt<float> py_res = np::zeros(dim, Cast<Py_ssize_t>(Size(ld.Box())).data(), np::getItemtype<float>());
  Array3d<float> res = Array3dFromPy<float>(py_res);
  res.move(ld.Box().min);
  FOR_BBOX3(p, ld.Box())
  {
    res(p) = sft.computeCurvature(p);
  }
  if (!return_stf)
    return py_res.getObject();

  const Int3 sz_ = Size(ld.Box());
  const Py_ssize_t dims[4] = { dim, sz_[0], sz_[1], sz_[2] };
  np::arrayt<float> acc_stf = np::zeros(4, dims, np::getItemtype<float>());
  FOR_BBOX3(p, ld.Box())
  {
    Float3 f = sft.computeForce(p);
    const Int3 q = p - ld.Box().min;
    for (int axis=0; axis<dim; ++axis)
    {
      acc_stf(axis, q[0], q[1], q[2]) = f[axis];
    }
  }
  return py::make_tuple(py_res.getObject(), acc_stf.getObject());
}
#endif

#if BOOST_VERSION>106300
void fill_with_smooth_delta(np::ndarray py_arr, float width)
{
//   np::arrayt<float> arr(py_arr);
//   for (int i=0; i<arr.shape()[0]; ++i)
//   {
//     arr(i) = my::smooth_delta_cos<float>(arr(i), width);
//   }
}
#else
void fill_with_smooth_delta(nm::array py_arr, float width)
{
  np::arrayt<float> arr(py_arr);
  for (int i=0; i<arr.shape()[0]; ++i)
  {
    arr(i) = my::smooth_delta_cos<float>(arr(i), width);
  }
}
#endif




namespace SurfaceTensionTestModeling
{

struct ConstFaceVarFunctor
{
  double operator()(const Int3 &p, const Int3 &nbp, int axis, int side) const { return 1; }
};


struct State
{
  Array3df rho;
  Levelset ls;
};

class Model
{
  double last_levelset_reinit;
  ptree params;
  ContinuumGrid grid;
  double kdiff, kstf, param_dt;
  
  FaceVarArrays velocity_fields;
  FaceVarArrays force_fields;
  double max_dt_diffusion, max_dt_velocity, max_dt_surfacetension, last_euler_dt;
  
public:
  typedef SurfaceTensionTestModeling::State State;
  //typedef ContainerOps<Array3dOps<float>, LevelsetOps> Ops;
  bool fail_flag;
  DomainDecomposition mtboxes;
  enum Border { BORDER = 3 };
  Array3dOps<float> ops_rho;
  LevelsetOps ops_ls;
  typedef NewSteppers::OpsCallback<State> Ops;
  Ops ops;
  
  Model(const LatticeDataQuad3d &ld_, int ndims_, const ptree &params_) :
    grid(ld_, ndims_), params(params_)
  {
    mtboxes.init(MakeMtBoxGrid(grid.Box(), Int3(64, 32, 32)));
    last_levelset_reinit = 0;
    max_dt_diffusion = max_dt_velocity = last_euler_dt = 0;
    velocity_fields.init(grid.ld.Box(), grid.dim);
    force_fields.init(grid.ld.Box(), grid.dim);
    kdiff = params.get("kdiff", 1.);
    kstf  = params.get("kstf", 1.);
    param_dt = params.get("dt", -1.);
    fail_flag = false;

    ops_rho.init(mtboxes, grid.Box(), grid.dim, BORDER);
    ops_ls.init(mtboxes, grid.Box(), grid.dim, BORDER);
    ops.initFrom = [=](State &u, const State &v, ConsMode mode)
    {
      this->ops_rho.initFrom(u.rho, v.rho, mode);
      this->ops_ls.initFrom(u.ls, v.ls, mode);
    };
    ops.addScaled = [=](double fu, State &u, double fv, const State &v)
    {
      this->ops_rho.addScaled(fu, u.rho, fv, v.rho);
      this->ops_ls.addScaled(fu, u.ls, fv, v.ls);
    };
  }

  void initState(State &state, bool clear)
  {
    ops_rho.init(state.rho, clear);
    ops_ls.init(state.ls, clear);
  }

  void preStep(State &state, StepControl &ctrl)
  {
    auto &rho = state.rho;
    auto &ls  = state.ls;
    double maxval  = rho.valueStatistics().MaxAbs();
    if (maxval > 1.e3)
    {
      fail_flag = true;
      return;
    }
    
    if (params.get<bool>("levelset_reinit", true))
    {
      if (last_levelset_reinit < ctrl.t - last_euler_dt*2.)
      {
        reinit_levelset(ls, grid.ld, grid.dim, 5, ptree());
        last_levelset_reinit = ctrl.t;
      }
    }
    cout << format("step t=%f") % ctrl.t << endl;
  }

  void calcSlope(const State &state, State &slope, StepControl &ctrl)
  {
    // grad p = f_sf - 1/k * v
    // => v = k f_sf - k grad p
    ops.initFrom(slope, state, CLEAN);

    #pragma omp parallel
    {
      BOOST_FOREACH(const BBox3 &bb, mtboxes.getCurrentThreadRange())
      {
        CopyBorder(const_cast<Array3d<float>&>(state.rho), bb, grid.Box(), grid.dim, BORDER);
        CopyBorder(const_cast<Array3d<float>&>(state.ls.phi), bb, grid.Box(), grid.dim, BORDER);
      }
    }

    const BBox3 extbox = ExtendForDim(grid.ld.Box(), grid.dim, 1);
    SurfaceTensionForce<float> stforce; stforce.init(extbox, state.ls.phi, grid.dim, grid.ld.Scale());

    // compute velocities
    ptree pt = make_ptree("geometric_coeff_mean", true);
    for (int axis=0; axis<grid.dim; ++axis)
    {
      velocity_fields[axis].fill(0);
      
      AddVelocitiesFromPotential2(grid.ld.Box(), axis, grid.dim, grid.ld.Scale(), ConstFaceVarFunctor(), state.rho, velocity_fields[axis], pt);

      FOR_BBOX3(p, grid.ir.faces[axis])
      {
        float f = -kstf*stforce.computeFaceForce(axis, p);
        velocity_fields[axis](p) += f;
        force_fields[axis](p) = f;
      }
      velocity_fields[axis] *= kdiff;
      // here we could multiply the velocity with a conductivity variable
    }

    max_dt_velocity = ComputeMaximalExplicitAdvectionEulerTimeStep(velocity_fields.data(), grid.dim, grid.ld.Scale());
    max_dt_diffusion = 0.8*0.5/kdiff*my::sqr(grid.ld.Scale())/(grid.dim+1.e-13);
    max_dt_surfacetension = 2./(kdiff*grid.dim*kstf)*my::cubed(grid.ld.Scale());
    last_euler_dt = param_dt>0 ? param_dt : my::min(max_dt_velocity, max_dt_diffusion, max_dt_surfacetension);
    ctrl.euler_dt.min(last_euler_dt);
    cout << "ctrl.euler_dt " << ctrl.euler_dt << endl;
    
    #pragma omp parallel
    {
      BOOST_FOREACH(const BBox3 &bb, mtboxes.getCurrentThreadRange())
      {
        AddKTAdvection(bb, state.ls.phi, velocity_fields.data(), grid.dim, grid.ld.Scale(), slope.ls.phi, NON_CONSERVATIVE);
        AddKTAdvection(bb, state.rho, velocity_fields.data(), grid.dim, grid.ld.Scale(), slope.rho);
      }
    }
    // debug
    #if 0
    std::vector<Image> images;
    images.push_back(DrawArray(velocity_fields[0], DrawArrayOpts().title("ux").outputRange()));
    images.push_back(DrawArray(velocity_fields[1], DrawArrayOpts().title("uy").outputRange()));
    images.push_back(DrawArray(state.levelset.phi, DrawArrayOpts().title("u").outputRange()));
    Image img;
    DrawImageGrid(img, images);
    img.Write(str(format("stf%010i.png") % int(ctrl.t*1000000)));
    #endif
  }

  void writeH5(H5::H5File f, H5::Group g, const State &state, double t, int out_num)
  {
    auto& levelset  = state.ls;
    auto& rho       = state.rho;
    
    //h5cpp::Group ld_group = RequireLatticeDataGroup(f.root(), "field_ld", grid.ld);
    H5::Group ld_group = RequireLatticeDataGroup(f, "field_ld", grid.ld);
    WriteScalarField(g, "ls", levelset.phi, grid.ld, ld_group);
    Array3d<float> theta; theta.initFromBox(grid.ld.Box());
    levelset.ToHeaviside(theta, grid.ld.Box(), grid.ld.Scale());
    writeAttrToGroup<double>(g, string("area"),  theta.valueStatistics().Sum()*std::pow(grid.ld.Scale(), grid.dim));
    //g.attrs().set("area", theta.valueStatistics().Sum()*std::pow(grid.ld.Scale(), grid.dim));
    theta *= rho;
    my::Averaged<double> stats = theta.valueStatistics();
//     g.attrs().set("mass", stats.Sum() * std::pow(grid.ld.Scale(), grid.dim));
//     g.attrs().set("rho_min", stats.Min());
//     g.attrs().set("rho_max", stats.Max());
//     g.attrs().set("max_dt_diffusion", max_dt_diffusion);
//     g.attrs().set("max_dt_velocity", max_dt_velocity);
    
    writeAttrToGroup<double>(g, string("mass"), stats.Sum() * std::pow(grid.ld.Scale(), grid.dim));
    writeAttrToGroup<double>(g, string("rho_min"), stats.Min());
    writeAttrToGroup<double>(g, string("rho_max"), stats.Max());
    writeAttrToGroup<double>(g, string("max_dt_diffusion"),  max_dt_diffusion);
    writeAttrToGroup<double>(g, string("max_dt_velocity"), max_dt_velocity);
    
    WriteScalarField(g, "rho", rho, grid.ld, ld_group);
    for (int axis=0; axis<grid.dim; ++axis)
    {
      LatticeDataQuad3d ldface = CellLdToFaceLd(grid.ld, grid.dim, axis);
      H5::Group ldf_group = RequireLatticeDataGroup(f, str(format("face_ld_%i") % axis), ldface);
      WriteScalarField(g, str(format("vel_%i") % axis), velocity_fields[axis], ldface, ldf_group);
      WriteScalarField(g, str(format("force_%i") % axis), force_fields[axis], ldface, ldf_group);
    }
  }
};


bool runSTF(const py::str &param_info_str, const py::object &py_ld, const py::dict &kwargs)
{
  FpExceptionStateGuard fpexceptguard_;
  ptree params = convertInfoStr(param_info_str, ConvectionModelTests::make_default_params());

  LatticeDataQuad3d ld = py::extract<LatticeDataQuad3d>(py_ld);
  Int3 size = ::Size(ld.Box());
  int ndims = size[1]<=1 ? 1 : (size[2]<=1 ? 2 : 3);

  Model model(ld, ndims, params);
  ObserverPde observer(params);
  
  State state;
  {
    model.initState(state, true);
    generateInitialConc(ld, ndims, params, kwargs, state.ls.phi);
    state.ls.max_val = state.ls.phi.valueStatistics().MaxAbs();
    state.rho.fill(1.);
  }

  auto stepper = NewSteppers::create<Model*,Model::Ops*,NewSteppers::NO_IMPLICIT>(
    params.get<string>("stepper"), &model, &(model.ops));

  auto doStep = [&](StepControl& ctrl) -> bool {
    ctrl.euler_dt = my::MakeSmallestMax<double>();
    model.preStep(state, ctrl);
    bool ok = stepper->doStep(state, ctrl);
    return ok && !model.fail_flag;
  };
  auto doObserve = [&](double t) {
    observer(t, state, model, params);
  };

  Py_BEGIN_ALLOW_THREADS
  NewSteppers::run(doStep, doObserve, params);
  Py_END_ALLOW_THREADS
  writeAttrToGroup<bool>(observer.openH5File().openGroup("/"), string("hard_fail"), model.fail_flag);
  //observer.openH5File().root().attrs().set("hard_fail", model.fail_flag);
  return !model.fail_flag;
}

  
}

#if 0
#include "trilinos_linsys_construction.h"
#include "trilinos_matrixfree.h"

namespace EllipticSolverTest
{

struct ConstFaceVarFunctor
{
  double val;
  ConstFaceVarFunctor(double val) : val(val) {}
  double operator()(const Int3& p, const Int3& q, int axis, int side) const { return val; }
  double operator()(int axis, const Int3& p) const { return val; }
};


// compute \laplace c - lambda c = -r
void testEllipticSolver(const Int3 size, bool write_image, bool apply_operator_instead_of_solve, double &result_time_ms, double &result_mem_bytes_per_dof)
{
  LatticeDataQuad3d ld;
  ld.Init(Int3(size), 1.);
  ld.SetWorldPosition(-ld.GetWorldBox().max*0.5);
  const float rmax = ld.GetWorldBox().max.norm();
  int num_dof = ld.NumSites();

  Epetra_SerialComm epetra_comm;
  Epetra_Map epetra_map(ld.NumSites(), 0, epetra_comm);
  DomainDecomposition mtboxes(MakeMtBoxGrid(ld.Box(), Int3(64, 32, 32)));
  
  typedef TrilinosMatrixFreeEllipticOperator<ConstFaceVarFunctor, Array3d<float> > OperatorType;
    
  Array3d<float> diag_offset(ld.Box(), Cons::DONT);
  Epetra_Vector rhs(epetra_map, false);
  Epetra_Vector lhs(rhs.Map(), false);

  #pragma omp parallel
  {
    BOOST_FOREACH(const BBox3 &bbox, mtboxes.getCurrentThreadRange())
    {
      FOR_BBOX3(p, bbox)
      {
        float r = ld.LatticeToWorld(p).norm();
        //float loc_rhs = my::smooth_heaviside<float>(rmax*0.2 - r, ld.Scale()*2);
        //rhs[ld.LatticeToSite(p)] = my::sqr((1./sqrt(2.))*ld.LatticeToWorld(p).dot(Float3(1.,-1, 0)));
  //      rhs[ld.LatticeToSite(p)] = std::max(std::abs(ld.LatticeToWorld(p)[0]), std::abs(ld.LatticeToWorld(p)[1]));
  //       diag_offset(p) = -10.;

        float loc_rhs = (p == Int3(20, 20, 0)) ? 1 : 0;
        rhs[ld.LatticeToSite(p)] = -loc_rhs;
        diag_offset(p) = -1./my::sqr(100*ld.Scale());
        lhs[ld.LatticeToSite(p)] = 0.;
      }
    }
  }

  rhs.Scale(-1.);
  boost::scoped_ptr<OperatorType> op(
    new OperatorType(
      ld, size[2]==1 ? 2 : 3, mtboxes, -1., ConstFaceVarFunctor(1.), diag_offset, epetra_map));

  boost::property_tree::ptree pt;
  pt.put("output", 1);
  pt.put("max_iter", 100);
  pt.put("throw", false);

  {
    uint64 mem_usage = GetMemoryUsage().rss;
    my::Time t_;
    if (apply_operator_instead_of_solve)
      op->Apply(rhs, lhs);
    else
      SolveEllipticEquation(*op, rhs, lhs, pt);
    result_time_ms = (my::Time() - t_).to_ms();
    result_mem_bytes_per_dof = double(GetMemoryUsage().rss_peak - mem_usage)/num_dof;
  }

  if (write_image)
  {
    Epetra_Vector residual(lhs.Map(), true);
    op->Apply(lhs, residual);

    Array3d<float> lhs_arr(ld.Box());
    Array3d<float> rhs_arr(ld.Box());
    Array3d<float> resid_arr(ld.Box());
    FOR_BBOX3(p, ld.Box())
    {
      lhs_arr(p) = lhs[ld.LatticeToSite(p)];
      rhs_arr(p) = -rhs[ld.LatticeToSite(p)];
      resid_arr(p) = residual[ld.LatticeToSite(p)];
    }

    DynArray<Image> images;
    images.push_back(DrawArray(lhs_arr, DrawArrayOpts().outputRange()));
    images.push_back(DrawArray(rhs_arr, DrawArrayOpts().outputRange()));
    images.push_back(DrawArray(resid_arr, DrawArrayOpts().outputRange()));
    Image bigimg;
    DrawImageGrid(bigimg, images);
    bigimg.Write("elliptic_solver_test_1.png");
  }
}

py::tuple testEllipticSolverPy(const Int3 &size, bool write_image, bool apply_operator_instead_of_solve)
{
  double result_time_ms, result_mem_bytes_per_dof;
  testEllipticSolver(size, write_image, apply_operator_instead_of_solve, result_time_ms, result_mem_bytes_per_dof);
  return py::make_tuple(result_time_ms, result_mem_bytes_per_dof);
}

  
};
#endif


void export_NumericalToolsTests()
{
#if 0
  py::def("convectionDiffusionReaction", ConvectionDiffusionReaction::runCDR);
py::def("testEllipticSolver", EllipticSolverTest::testEllipticSolverPy);
#endif
  py::def("levelsetRedistancing", runTestRedistancing);
  py::def("zalesakDisk", zalezakDiskFunction);
  py::def("rotationalVelocityField", rotationalVelocityField);
  py::def("calcCurvature", calcCurvature);
  py::def("fill_with_smooth_delta", fill_with_smooth_delta);
  py::def("convectionTest", ConvectionModelTests::runTest);
  py::def("convectionLevelset", LevelsetConvection::runLS);
  py::def("surfaceTensionForceTest",SurfaceTensionTestModeling::runSTF);
  py::def("steppersConvergenceTest", SimpleSteppersTest::runSimpleSteppersTest);
}


