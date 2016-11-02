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
#if 0

#include <fenv.h>
#include <boost/array.hpp>
#include <boost/assign/std/vector.hpp>
#include <boost/foreach.hpp>


#include "time_stepper_utils.h"
#include "levelset.h"
#include "continuum-flow.h"
#include "shared-objects.h"

using Steppers::StepControl;



inline float closest_intersect(float a, float b)
{
  return std::min(a,b);
}


void makeZalezakDisk(const LatticeDataQuad3d &ld, Array3d<float> dist)
{
  float rad = 15;
  Float3 pos(50,75,0);
  float width = 5;
  float len = 25;
  float len_ext = 10;
  Float3 slot_pos(0, 0 - rad + len/2. - len_ext, 0.);
  FOR_BBOX3(ip, ld.Box())
  {
    Float3 p = ld.LatticeToWorld(ip)/ld.Scale() - pos;
    float d1 = rad - p.norm();
    float d2 = -my::min(width/2.-std::abs(p[0] - slot_pos[0]), (len+len_ext)/2.-std::abs(p[1]-slot_pos[1]), width/2.-std::abs(p[2]-slot_pos[2]));
    dist(ip) = ld.Scale()*closest_intersect(d1,d2);
  }
}




void run(int dim)
{
  ptree p;
  const string fn = "zalesak";
  p.put("in_fn", fn+".png");
  p.put("fn_out", fn+"_res");

  my::SetNumThreads(2);
  //const ptree &params = p;

  Array3d<float> a;
  LatticeDataQuad3d ld;
  string in_fn = p.get<string>("in_fn");
  if (in_fn == "zalesak.png")
  {
    ld.Init(BBox3(0,0,0,99,99,0), 10.);
    a.initFromBox(ld.Box());
    makeZalezakDisk(ld, a);
    ld.SetWorldPosition(-ld.GetWorldBox().max * 0.5); // set origin = lower left side
    DrawArray(a, DrawArrayOpts().outputRange()).Write(p.get<string>("in_fn"));
  }
  else
  {
    a = ReadArrayFromImage<float>(in_fn,0,1);
    Int3 size = a.size();
    if (dim > 2)
      size[2] = 20;
    ld.Init(size, 10);
    a += -0.5;
    a *= 2.*10*ld.Scale();
    ld.SetWorldPosition(-ld.GetWorldBox().max * 0.5); // set origin = lower left side
  }
  Int3 size = Size(ld.Box());

  Levelset levelset;
  levelset.init(ld, dim, 2);
  for (int z=0; z<size[2]; ++z)
  {
    Array3d<float> p = levelset.phi[Set(ld.Box(), 2, z)];
    p.fill(a);
    if (size[2]>3)
    {
      float f = ld.Scale()*10*2*(std::abs(0.5 - (float(z)/size[2])));
      p -= f;
    }
  }

  ptree pinit = p;  pinit.put("fn_out", fn+"_reinitanim");
  reinit_levelset(levelset, ld, dim, 5, pinit);

  if (size[2]==1)
  {
    Image img;
    levelset.Draw(img, ld, DrawArrayOpts().outputRange());
    img.Write(p.get<string>("fn_out")+"_real.png");
  }
  else
  {
    Image img, bigimg;
    std::vector<int> zs; boost::assign::push_back(zs)(1)(size[2]/2)(size[2]-2);
    std::vector<Image> imgs;
    BOOST_FOREACH(int z, zs)
    {
      levelset.Draw(img, ld, DrawArrayOpts().zslice(z).outputRange().title(str(format("z=%i") % z)));
      imgs.push_back(img);
    }
    DrawImageGrid(bigimg, imgs, 1);
    bigimg.Write(p.get<string>("fn_out")+"_real.png");
  }

  double g = calc_levelset_badness(levelset, ld, dim);
  cout << "badness: " << g << endl;
  
//   p.put("fn_out", fn+"_res2");
//   reinit_levelset(levelset, ld, dim, 5, params);
}



namespace LevelsetConvection
{

struct State
{
  Levelset levelset;

  void addScaled(double s, const State &u, double s_self = 1.)
  {
    levelset.addScaled(s, u.levelset, s_self);
  }

  void setId(int id_) {}

  void init(const LatticeDataQuad3d &grid, int ndims)
  {
    levelset.init(grid, ndims, 3);
  }

  void initCloned(const State &other)
  {
    levelset.initCloned(other.levelset);
  }
};



class Model
{
public:
  typedef LevelsetConvection::State State;
  enum { stepperFlags = Steppers::HAS_PRESTEP };
  
  void init(const LatticeDataQuad3d &grid_, int ndims_,  const ptree &params)
  {
    grid = grid_;
    ndims = ndims_;
    LatticeIndexRanges ir = LatticeIndexRangesFromCellRange(grid.Box(), ndims);
    Float3 rot_axis = ndims==2 ? Float3(0,0,1) : Float3(1,1,1).normalized();
    for (int i=0; i<ndims; ++i)
    {
      velocity_fields[i].initFromBox(ir.faces[i]);
      LatticeDataQuad3d ld_face(ir.faces[i], grid.Scale());
      Vec<bool,3> cc(true); cc[i] = false;
      ld_face.SetCellCentering(cc);
      ld_face.SetWorldPosition(grid.LatticeToWorld(Int3(0))-Float3(0.5*grid.Scale()));
      cout << "ld face" << ld_face << endl;
      FOR_BBOX3(p, ir.faces[0])
      {
        Float3 wp = ld_face.LatticeToWorld(p);
        velocity_fields[i](p) = cross(rot_axis, wp)[i];
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
    max_dt = 0.8 * grid.Scale() / max_vel / (ndims);
    last_levelset_reinit = 0;
    cout << "init done" << endl;
    cout << "grid = "; grid.print(cout); cout<<endl;
    cout << "dt = " << max_dt << endl;
  }

  void preStep(State &state, StepControl &ctrl)
  {
#if 1
    //if (last_levelset_reinit < ctrl.t - max_dt * 20)
    {
      double gnorm = calc_levelset_badness(state.levelset, grid, ndims);
      cout << format("gnorm = %f, maxval = %f") % gnorm % state.levelset.max_val << endl;
      if (gnorm > 0.0001)
      {
        ptree pt = reinit_levelset(state.levelset, grid, ndims, 5);
        cout << "levelset reinit" << endl;
        //last_levelset_reinit = ctrl.t;
      }
    }
    //else
#endif
    cout << format("step t=%f") % ctrl.t << endl;
  }

  void calcSlope(const State &state, State &slope, StepControl &ctrl)
  {
    LatticeIndexRanges ir = LatticeIndexRangesFromCellRange(grid.Box(), ndims);
    ctrl.euler_dt = max_dt;
    DynArray<BBox3> mtboxes = MakeMtBoxGrid(ir.cells);
    
    if (slope.levelset.phi.empty())
      slope.init(grid, ndims);
    slope.levelset.phi.fill(0);

#ifdef USE_HIGH_ORDER_LEVELSET_METHODS
    CopyBorder(const_cast<State&>(state).levelset.phi[ir.cells], ndims, 3);
    #pragma omp parallel
    {
      //FaceVarArrays facevalues;
      #pragma omp for schedule(dynamic, 1)
      for (int i=0; i<mtboxes.size(); ++i)
      {
        const BBox3 bb = mtboxes[i];
        //InitFaceVar(bb, ndims, facevalues);
        //ENOReconstructFaceValues(bb, ndims, grid.Scale(), state.levelset.phi, velocity_fields.data(), facevalues.data());
        //AddAdvection(bb, ndims, grid.Scale(), velocity_fields.data(), facevalues.data(), slope.levelset.phi, NON_CONSERVATIVE);
        AddENOAdvection(bb, state.levelset.phi, velocity_fields.data(), ndims, grid.Scale(), slope.levelset.phi, NON_CONSERVATIVE, 3);
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
#else
    CopyBorder(const_cast<State&>(state).levelset.phi[ir.cells], ndims, 2);
    AddFctSlope2<float>(ir.cells, state.levelset.phi, velocity_fields.data(), ndims, grid.Scale(), slope.levelset.phi, NON_CONSERVATIVE);
#endif
  }

  void writeH5(h5cpp::File f, h5cpp::Group g, State &state, double t, int out_num)
  {
    h5cpp::Group ld_group = RequireLatticeDataGroup(f.root(), "field_ld", grid);
    WriteScalarField(g, "ls", state.levelset.phi, grid, ld_group);
    Array3d<float> theta; theta.initFromBox(grid.Box());
    state.levelset.ToHeaviside(theta, grid.Box(), grid.Scale());
    WriteScalarField(g, "theta", theta, grid, ld_group);
    float mass = theta.valueStatistics().Sum() * std::pow(grid.Scale(), ndims);
    g.attrs().set("mass", mass);
  }

  void appendToImages(State &state, std::vector<Image> &images) const
  {
    Image img;
    DrawArray(img, state.levelset.phi[grid.Box()], DrawArrayOpts().outputRange());
    images.push_back(img);
    Array3d<float> theta; theta.initFromBox(grid.Box());
    state.levelset.ToHeaviside(theta, grid.Box(), grid.Scale());
    DrawArray(img, theta[grid.Box()], DrawArrayOpts());
    images.push_back(img);
  }

  double max_dt;
  double last_levelset_reinit;
  int ndims;
  LatticeDataQuad3d grid;
  boost::array<Array3d<float>, 3> velocity_fields;
};


void run_convective(ptree p)
{
#if 0
  string fn_in = params.get<string>("fn_in");
  Array3d<float> a = ReadArrayFromImage<float>(fn_in,0,1);
  LatticeDataQuad3d ld;
  float scale = params.get<float>("scale");
  ld.Init(a.size(), scale);
  a += -0.5;
  a *= 2.*10*ld.Scale();
  ld.SetWorldPosition(-ld.GetWorldBox().max * 0.5); // set origin = lower left side
  int dim = a.size()[1] > 1 ? 2 : 1;
#else
  const int dim = 2;
  LatticeDataQuad3d ld(BBox3(0,0, dim==2 ? 0 : -50, 99, 99, dim==2 ? 0 : 49), 10.);
  ld.SetCellCentering(Bool3(true));
  Array3d<float> a(ld.Box());
  makeZalezakDisk(ld, a);
  Float3 world_center(-(ld.GetWorldBox().max-ld.GetWorldBox().min)*0.5); world_center[2]=0.;
  ld.SetWorldPosition(world_center); // set origin = lower left side
  {
    h5cpp::File f("zalesak.h5","w");
    h5cpp::Group gld = f.root().create_group("field_ld");
    WriteHdfLd(gld, ld);
    WriteScalarField(f.root(), "distance", a, ld, gld);
  }
  p.put("tend", my::mconst::pi2());
  p.put("out_intervall", my::mconst::pi2()/20.);
  p.put("save_hdf", true);
  p.put("save_images", false);
  p.put("num_threads", 3);
  p.put("fn_out", str(format("zalesak_conv_%i_restart") % dim));
#endif

  cout << "ld: " << ld << endl;
  
  my::SetNumThreads(p.get<int>("num_threads"));
  
  Model model;
  State state;
  ObserverPde observer;

  model.init(ld, dim, p);
  state.init(ld, dim);
  state.levelset.phi[ld.Box()].fill(a);
  state.levelset.max_val = state.levelset.phi.valueStatistics().MaxAbs();
  state.levelset.ToHeaviside(a, ld.Box(), ld.Scale());
  double mass_init = a.valueStatistics().Sum();
  
  reinit_levelset(state.levelset, ld, dim, 5);

#ifdef USE_HIGH_ORDER_LEVELSET_METHODS
  Steppers::Stepper<Model>::BasePointer stepper = Steppers::Stepper<Model>::create(p.get("stepper", "rk3"));
#else
  Steppers::Stepper<Model>::BasePointer stepper = Steppers::Stepper<Model>::create(p.get("stepper", "impeuler"));
#endif
  ::run<>(state, model, observer, *stepper, p);

  state.levelset.ToHeaviside(a, ld.Box(), ld.Scale());
  double mass_end = a.valueStatistics().Sum();

  cout << format("mass error %f -> %f, %f %%") % mass_init % mass_end % ((mass_end-mass_init)*100./mass_init) << endl;
}


void run_convective()
{
  ptree p;
  p.put("scale", 10.);
  p.put("fn_in", "levelset2.png");
  p.put("fn_out", "levelset2_bfecc_conv");
  p.put("tend", 1);
  p.put("out_intervall", 0.01);
  p.put("save_hdf", false);
  p.put("save_images", true);
  run_convective(p);
}


}





int main(int argc, char **argv)
{
  my::MultiprocessingInitializer mpinit(argc, argv);
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
  //run(2);
  LevelsetConvection::run_convective();
}




#if 0
void run_distancemap(const ptree &params)
{
  LatticeDataQuad3d ld;
  int ndims;

  string in_fn = params.get<string>("in_fn");
  Array3d<float> a = ReadArrayFromImage<float>(in_fn,0,1);
  a *= 2.;
  a -= 1.;
  ld.Init(a.size(), 0.2);
  ld.SetWorldPosition(-ld.GetWorldBox().max * 0.5); // set origin = lower left side
  ndims = a.size()[1] > 1 ? 2 : 1;

  FOR_BBOX3(p, ld.Box())
  {
    float val = a(p);
    if (std::abs(val) > 0.5)
      a(p) = std::numeric_limits<float>::max();
    else
      a(p) = std::abs(a(p));
  }

  DistanceFieldComputer c;
  c.Do(ld, a);

  Image img;
  DrawArray(img, a, DrawArrayOpts().outputRange());
  img.Write(params.get<string>("fn_out")+"_distancemap.png");
}
#endif

#endif

int main(int argc, char **argv)
{
  return 0;
}

