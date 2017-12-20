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
#include "levelset.h"
#include "mwlib/time_stepper.h"
#include "hdfio.h"
#include "continuum-flow.h"
#include "mwlib/ptree_ext.h"
#include "shared-objects.h"
#include "time_stepper_utils_new.h"

#include <boost/foreach.hpp>


using boost::property_tree::ptree;
using NewSteppers::StepControl;

inline float smooth_heaviside(float x, float radius)
{
  return my::smooth_heaviside_sin(x, radius*0.5);
}
inline float smooth_sgn(float x, float radius)
{
  return 2.*(smooth_heaviside(x,radius)-0.5);
}
inline float smooth_delta(float x, float radius)
{
  return my::smooth_delta_cos(x, radius*0.5);
}
enum
{
  SITE_INACTIVE=0,
  SITE_ACTIVE=2,
  SITE_BOUNDARY=4
};


void Levelset::Draw(Image &img, const LatticeDataQuad3d &ld, const DrawArrayOpts &opts, const Array3d<char>* active)
{
  BBox3 bb = ld.Box();
  int z;
  if (opts.slice_axis_==2)
    z = opts.slice_;
  else
    z = (bb.min[2]+bb.max[2])/2;
  bb.Set(2, z);
  float scale = 1.;
  if (!opts.normalize_)
    scale = opts.scale_;
  
  //ConstArray3d<float> phi = this->phi[bb];
  Int3 s = Size(bb); //phi.size();
  s[2] = 1;

  my::MinMax<float> mmax;
  img.Init(s[0],s[1]);
  FOR_BBOX3(p, bb)
  {
    float g = phi(p) * scale / ld.Scale() * 0.1;
    mmax.add(phi(p));
    uchar cg = g*255;
    uchar col[3] = { cg, cg, cg };
    if (g < 0) col[2] = 0;
    if (active)
    {
      if ((*active)(p) == SITE_INACTIVE)
        col[1] = col[2] = 0;
      else if((*active)(p) == SITE_BOUNDARY)
        col[0] = col[1] = 0;
    }
    img.SetPixel(p[0], p[1], col[0], col[1], col[2]);
  }
  if (opts.outputRange_)
  {
    img.SetColor(255, 255, 255);
    img.DrawText(2,2, str(format("min=%f\nmax=%f") % mmax.min % mmax.max));
    //img.DrawText(2,2, str(format("max=%f") % this->max_val));
  }
}


struct MakeThetaOp
{
  float w, w_inv;
  MakeThetaOp(float w) : w(w), w_inv(1./w) {}
  void inplace(float &t, float p)
  {
    /* linear no good!
    t = my::cut(p, -w, w) * 0.5*w_inv + 0.5;
    */
    // phi(x) = (1+exp(-x*e/w*2.))^-1 // makes a transition from 0 to 1 over -w/2 to w/2
    //const double e = 2.718;
    //t = 1./(1+std::exp(-p*e*2.0*w_inv));

    // better yet, sin based heaviside; is widely used
    t = smooth_heaviside(p, w);
  }
};


void Levelset::ToHeaviside(Array3d<float> theta, const BBox3 &bbox, double width) const
{
  myAssert(width < this->max_val);
  Array3d<float> dst(theta[bbox]);
  ConstArray3d<float> src(phi[bbox]);
  _ArrayNdInternal::apply(dst, src, MakeThetaOp(width));
}

float Levelset::ToHeaviside(float phi, float width)
{
  return smooth_heaviside(phi, width);
};


struct LevelSetReinit
{
  typedef Levelset State;
  enum { stepperFlags = /*Steppers::HAS_PRESTEP |*/ Steppers::HAS_POSTSTEP };
  
  LatticeDataQuad3d ld;
  int dim;
  Array3d<char> active;
  const DomainDecomposition *mtboxes;
  DomainDecomposition mtboxes_active_list; /// Only perform redistancing in sub-boxes near the zero-level.
  int gradient_approx_method;
  LevelsetOps levelsetops;
#ifdef USE_HIGH_ORDER_LEVELSET_METHODS
  NewSteppers::Rk3<LevelSetReinit*, LevelsetOps*> stepper;
#else
  NewSteppers::ImprovedExplicitEuler<LevelSetReinit*, LevelsetOps*> stepper;
#endif
  
  void init(const LatticeDataQuad3d &ld_, int dim_, const DomainDecomposition &mtboxes_, const State &state, const ptree &params) //, bool do_full_restart_ = false)
  {
    ld = ld_;
    dim = dim_;
    mtboxes = &mtboxes_;
    levelsetops.init(mtboxes_, ld.Box(), dim, 3);
    stepper.set_model(this);
    stepper.set_ops(&levelsetops);
    
    string s = params.get<string>("gradient_approx_method", "weno5");
    if (s == "weno5") {
      gradient_approx_method = 5;
    } else if (s == "weno3") {
      gradient_approx_method = 3;
    } else if (s == "upwind") {
      gradient_approx_method = 1;
    }
    else throw std::runtime_error(str(format("unkown gradient approximation method %s") % s));

    Array3dOps<char>(*mtboxes, ld.Box(), dim, 1).init(active, true);
  }

  

  
#if 0
  Array3d<float> integration_stencil,
                 sussman_constraint_const_factor;
  Array3d<float> init_phi, last_phi;
  bool do_full_restart;

  
  Array3d<float> initIntegrationStencil(int dim)
  {
    const float coeff[3] = { 1./8., 6./8., 1./8. };
    BBox3 b(0,0,0,0,0,0);
    for (int i=0; i<dim; ++i)
      b.max[i] = 2;
    Array3d<float> s(b);
    if (dim > 1)
    {
      Array3d<float> t = initIntegrationStencil(dim-1);
      for (int j=0; j<3; ++j)
      {
        b.min[dim-1] = b.max[dim-1] = j;
        s[b].addScaled(coeff[j], t);
      }
    }
    else
    {
      for (int j=0; j<3; ++j) s(Int3(j,0,0)) = coeff[j];
    }
    return s;
  }
#endif

  /// Determine what boxes are active (considered in the update)
  void fillActiveByValue(ConstArray3d<float> phi, float max_value)
  {
    #pragma omp parallel
    {
      BOOST_FOREACH(const BBox3 &bb, mtboxes->getCurrentThreadRange())
      {
        bool is_box_active = false;
        FOR_BBOX3(p, bb)
        {
          char a = (std::abs(phi(p))<max_value) ? SITE_ACTIVE : SITE_INACTIVE;
          active(p) = a;
          if (a != SITE_INACTIVE)
            is_box_active = true;
        }
        if (is_box_active)
          mtboxes_active_list.insert(my::OmpGetCurrentThread(), bb);
      }
    }
  }

  /// Determine what boxes are active (considered in the update)
  void fillActiveByWidth(ConstArray3d<float> &phi, int active_width_)
  {
    const BBox3 bb = ld.Box();
    active_width_ -= 1;
    
    DynArray<Int3> sites1(1024,ConsTags::RESERVE), sites2(1024,ConsTags::RESERVE);
    FOR_BBOX3(p, bb)
    {
      float v = phi(p);
      bool a = false;
      for (int dir = 0; dir<2*dim && !a; ++dir)
      {
        Int3 nbp = ld.NbLattice(p, dir);
        if ((v>0) != (phi(nbp)>0))
          a = true;
      }
      if (a)
      {
        active(p) = SITE_ACTIVE;
        sites1.push_back(p);
      }
    }

    for (int width_iter=0; width_iter<active_width_+1; ++width_iter)
    {
      for (int i=0; i<sites1.size(); ++i)
      {
        const Int3 p = sites1[i];
        Int3 pp(p);
        for (int axis = 0; axis<dim; ++axis)
        {
          for (int i=-1; i<=1; i+=2)
          {
            pp[axis] = p[axis]+i;
            if (active(pp) == SITE_INACTIVE && ld.IsInsideLattice(pp))
            {
              active(pp) = (width_iter==(active_width_)) ? SITE_BOUNDARY : SITE_ACTIVE;
              sites2.push_back(pp);
            }
          }
          pp[axis] = p[axis];
        }
      }
      sites1.swap(sites2);
      sites2.remove_all();
    }

    tbb::spin_mutex mutex;
    #pragma omp parallel
    {
      DynArray<BBox3> active_boxes(1024, ConsTags::RESERVE);
      BOOST_FOREACH(const BBox3 &bb, mtboxes->getCurrentThreadRange())
      {
        bool is_box_active = false;
        FOR_BBOX3(p, bb)
        {
          if (active(p) != SITE_INACTIVE)
          {
            is_box_active = true;
          }
        }
        if (is_box_active)
        {
          mutex.lock();
          mtboxes_active_list.insert(my::OmpGetCurrentThread(), bb);
          mutex.unlock();
        }
      }
    }
    //cout << "num activated " << mtboxes_active_list.size() << endl;
  }


  /**
  @brief Fill the areas of the distancemap that were excluded from the redistancing with an updated value.
  
  This updated value is taken as the smallest distance encounted on the boundary between the recomputed area
  around the zero-level and the in-active area elsewhere.
  */
  void finalize(State &state, float width)
  {
    //cerr <<  __func__  << endl;
    Array3d<float> phi = state.phi;

    //int num_visited = 0, num_active_mtboxes = 0;
    float min_abs_val = width * 10;
    #pragma omp parallel
    {
      float th_min_abs_val = width * 10; 
      BOOST_FOREACH(const BBox3 &bb, mtboxes_active_list.getCurrentThreadRange())
      {
        FOR_BBOX3(p, bb)
        {
          if (active(p) == SITE_BOUNDARY)
          {
            th_min_abs_val = std::min(std::abs(phi(p)), th_min_abs_val);
          }
        }
//         #pragma omp critical
//         {
//           num_active_mtboxes ++;
//           cerr << "b" << i << " on" << my::OmpGetCurrentThread() << endl;
//         }
      }
      #pragma omp critical
      {
            min_abs_val = std::min(th_min_abs_val, min_abs_val);
      }
      #pragma omp barrier

      BOOST_FOREACH(const BBox3 &bb, mtboxes->getCurrentThreadRange())
      {
        FOR_BBOX3(p, bb)
        {
          float f = phi(p);
          if (!active(p))
            phi(p) = f<0 ? -min_abs_val : min_abs_val;
          else
            phi(p) = my::sign(phi(p)) * std::min(min_abs_val, std::abs(phi(p)));
        }
      }
    }
    //cerr << "num_visited = " << num_visited << ", max val = " << min_abs_val << " amtb = " << num_active_mtboxes << endl;

    state.max_val = min_abs_val;
  }


  /**
  @brief I forgot how this works. Please look in my collection of Levelset papers.
  */
  double globalGradientDeviation(const ConstArray3d<float> &phi)
  {
    double res = 0.;
    double cnt = 0.;

    #pragma omp parallel reduction(+:res) reduction(+:cnt)
    {
      Array3d<float> grad_buffer;

      BOOST_FOREACH(const BBox3 &bb, mtboxes_active_list.getCurrentThreadRange())
      {
        grad_buffer.initFromBox(bb);
        calcGradientNorm(bb, phi, phi, grad_buffer);
        FOR_BBOX3(p, bb)
        {
          if (active(p) != SITE_ACTIVE) continue;
          double g = std::abs(grad_buffer(p) - 1.);
          res += g*g;
          cnt += 1;
        }
//         #pragma omp critical
//         {
//           ++num_grad_boxes
//           cerr << "b" << box_index << " on" << my::OmpGetCurrentThread() << endl;
//         }
      }
    }
//    cout << format("glob grad nb=%i, n=%i") % num_grad_boxes % cnt << endl;
    return std::sqrt(res)/dim / cnt;
  }


  /**
  @brief I computed these coefficients in there with a program called pyweno.
  */
  void calcGradientNorm(const BBox3 &bbox, ConstArray3d<float> phi_sign, ConstArray3d<float> phi, Array3d<float> res) const
  {
    float factor = 1./ld.Scale();

    FOR_BBOX3(p, bbox)
    {
      if(active(p) == SITE_INACTIVE) continue;
      float D_p=0., D_m=0.;
      for (int axis=0; axis<dim; ++axis)
      {
        float d_p = 0, d_m = 0;
        if (gradient_approx_method == 5)
        {
          const int fsi = 1, i = 2;
          float sigma0, sigma1, sigma2, omega0, omega1, omega2, acc, fr0, fr1, fr2, fs0;
          float values_buffer[7];
          float *f = values_buffer+3;
          for (int i=-3; i<=3; ++i) f[i] = phi(add_to_axis(p, axis, i));

          f = values_buffer+1;
          // for bias ((+)) ---
          //--- smoothness ---
          sigma0 = 3.3333333333333333333333333333333333*f[(i+0)*fsi]*f[(i+0)*fsi] - 17.0*f[(i+0)*fsi]*f[(i+1)*fsi] + 14.0*f[(i+0)*fsi]*f[(i+2)*fsi] - 3.6666666666666666666666666666666667*f[(i+0)*fsi]*f[(i+3)*fsi] + 22.0*f[(i+1)*fsi]*f[(i+1)*fsi] - 37.0*f[(i+1)*fsi]*f[(i+2)*fsi] + 10.0*f[(i+1)*fsi]*f[(i+3)*fsi] + 16.0*f[(i+2)*fsi]*f[(i+2)*fsi] - 9.0*f[(i+2)*fsi]*f[(i+3)*fsi] + 1.3333333333333333333333333333333333*f[(i+3)*fsi]*f[(i+3)*fsi];
          sigma1 = 10.0*f[(i+0)*fsi]*f[(i+0)*fsi] - 19.0*f[(i+0)*fsi]*f[(i+1)*fsi] + 6.0*f[(i+0)*fsi]*f[(i+2)*fsi] - 7.0*f[(i+0)*fsi]*f[(i-1)*fsi] + 10.0*f[(i+1)*fsi]*f[(i+1)*fsi] - 7.0*f[(i+1)*fsi]*f[(i+2)*fsi] + 6.0*f[(i+1)*fsi]*f[(i-1)*fsi] + 1.3333333333333333333333333333333333*f[(i+2)*fsi]*f[(i+2)*fsi] - 1.6666666666666666666666666666666667*f[(i+2)*fsi]*f[(i-1)*fsi] + 1.3333333333333333333333333333333333*f[(i-1)*fsi]*f[(i-1)*fsi];
          sigma2 = 22.0*f[(i+0)*fsi]*f[(i+0)*fsi] - 17.0*f[(i+0)*fsi]*f[(i+1)*fsi] - 37.0*f[(i+0)*fsi]*f[(i-1)*fsi] + 10.0*f[(i+0)*fsi]*f[(i-2)*fsi] + 3.3333333333333333333333333333333333*f[(i+1)*fsi]*f[(i+1)*fsi] + 14.0*f[(i+1)*fsi]*f[(i-1)*fsi] - 3.6666666666666666666666666666666667*f[(i+1)*fsi]*f[(i-2)*fsi] + 16.0*f[(i-1)*fsi]*f[(i-1)*fsi] - 9.0*f[(i-1)*fsi]*f[(i-2)*fsi] + 1.3333333333333333333333333333333333*f[(i-2)*fsi]*f[(i-2)*fsi];
          //--- weights ---
          omega0 = 0.1/(1.0e-6*1.0e-6 + 2*1.0e-6*sigma0 + sigma0*sigma0);
          omega1 = 0.6/(1.0e-6*1.0e-6 + 2*1.0e-6*sigma1 + sigma1*sigma1);
          omega2 = 0.3/(1.0e-6*1.0e-6 + 2*1.0e-6*sigma2 + sigma2*sigma2);
          acc = omega0 + omega1 + omega2;
          omega0 = omega0/acc;
          omega1 = omega1/acc;
          omega2 = omega2/acc;
          //--- reconstructions ---
          fr0 = -1.8333333333333333333333333333333333*f[(i+0)*fsi] + 3.0*f[(i+1)*fsi] - 1.5*f[(i+2)*fsi] + 0.33333333333333333333333333333333333*f[(i+3)*fsi];
          fr1 = -0.5*f[(i+0)*fsi] + f[(i+1)*fsi] - 0.16666666666666666666666666666666667*f[(i+2)*fsi] - 0.33333333333333333333333333333333333*f[(i-1)*fsi];
          fr2 = 0.5*f[(i+0)*fsi] + 0.33333333333333333333333333333333333*f[(i+1)*fsi] - f[(i-1)*fsi] + 0.16666666666666666666666666666666667*f[(i-2)*fsi];
          fs0 = fr0*omega0 + fr1*omega1 + fr2*omega2;
          d_p = fs0;

          f = values_buffer;
          // for bias ((-)) ---
          //--- smoothness ---
          sigma0 = 3.3333333333333333333333333333333333*f[(i+0)*fsi]*f[(i+0)*fsi] - 17.0*f[(i+0)*fsi]*f[(i+1)*fsi] + 14.0*f[(i+0)*fsi]*f[(i+2)*fsi] - 3.6666666666666666666666666666666667*f[(i+0)*fsi]*f[(i+3)*fsi] + 22.0*f[(i+1)*fsi]*f[(i+1)*fsi] - 37.0*f[(i+1)*fsi]*f[(i+2)*fsi] + 10.0*f[(i+1)*fsi]*f[(i+3)*fsi] + 16.0*f[(i+2)*fsi]*f[(i+2)*fsi] - 9.0*f[(i+2)*fsi]*f[(i+3)*fsi] + 1.3333333333333333333333333333333333*f[(i+3)*fsi]*f[(i+3)*fsi];
          sigma1 = 10.0*f[(i+0)*fsi]*f[(i+0)*fsi] - 19.0*f[(i+0)*fsi]*f[(i+1)*fsi] + 6.0*f[(i+0)*fsi]*f[(i+2)*fsi] - 7.0*f[(i+0)*fsi]*f[(i-1)*fsi] + 10.0*f[(i+1)*fsi]*f[(i+1)*fsi] - 7.0*f[(i+1)*fsi]*f[(i+2)*fsi] + 6.0*f[(i+1)*fsi]*f[(i-1)*fsi] + 1.3333333333333333333333333333333333*f[(i+2)*fsi]*f[(i+2)*fsi] - 1.6666666666666666666666666666666667*f[(i+2)*fsi]*f[(i-1)*fsi] + 1.3333333333333333333333333333333333*f[(i-1)*fsi]*f[(i-1)*fsi];
          sigma2 = 22.0*f[(i+0)*fsi]*f[(i+0)*fsi] - 17.0*f[(i+0)*fsi]*f[(i+1)*fsi] - 37.0*f[(i+0)*fsi]*f[(i-1)*fsi] + 10.0*f[(i+0)*fsi]*f[(i-2)*fsi] + 3.3333333333333333333333333333333333*f[(i+1)*fsi]*f[(i+1)*fsi] + 14.0*f[(i+1)*fsi]*f[(i-1)*fsi] - 3.6666666666666666666666666666666667*f[(i+1)*fsi]*f[(i-2)*fsi] + 16.0*f[(i-1)*fsi]*f[(i-1)*fsi] - 9.0*f[(i-1)*fsi]*f[(i-2)*fsi] + 1.3333333333333333333333333333333333*f[(i-2)*fsi]*f[(i-2)*fsi];
          //--- weights ---
          omega0 = 0.3/(1.0e-6*1.0e-6 + 2*1.0e-6*sigma0 + sigma0*sigma0);
          omega1 = 0.6/(1.0e-6*1.0e-6 + 2*1.0e-6*sigma1 + sigma1*sigma1);
          omega2 = 0.1/(1.0e-6*1.0e-6 + 2*1.0e-6*sigma2 + sigma2*sigma2);
          acc = omega0 + omega1 + omega2;
          omega0 = omega0/acc;
          omega1 = omega1/acc;
          omega2 = omega2/acc;
          //--- reconstructions ---
          fr0 = -0.33333333333333333333333333333333333*f[(i+0)*fsi] - 0.5*f[(i+1)*fsi] + f[(i+2)*fsi] - 0.16666666666666666666666666666666667*f[(i+3)*fsi];
          fr1 = -f[(i+0)*fsi] + 0.5*f[(i+1)*fsi] + 0.33333333333333333333333333333333333*f[(i+2)*fsi] + 0.16666666666666666666666666666666667*f[(i-1)*fsi];
          fr2 = -3.0*f[(i+0)*fsi] + 1.8333333333333333333333333333333333*f[(i+1)*fsi] + 1.5*f[(i-1)*fsi] - 0.33333333333333333333333333333333333*f[(i-2)*fsi];
          fs0 = fr0*omega0 + fr1*omega1 + fr2*omega2;
          d_m = fs0;
        }
        else if (gradient_approx_method == 3)
        {
          const int fsi = 1, i = 1;
          float sigma0, sigma1, omega0, omega1, acc, fr0, fr1, fs0;
          float values_buffer[5];
          float *f = values_buffer+2;
          for (int i=-2; i<=2; ++i) f[i] = phi(add_to_axis(p, axis, i));

          f = values_buffer+1;
          // for bias ((+)) ---
          //--- smoothness ---
          sigma0 = f[(i+0)*fsi]*f[(i+0)*fsi] - 4.0*f[(i+0)*fsi]*f[(i+1)*fsi] + 2.0*f[(i+0)*fsi]*f[(i+2)*fsi] + 4.0*f[(i+1)*fsi]*f[(i+1)*fsi] - 4.0*f[(i+1)*fsi]*f[(i+2)*fsi] + f[(i+2)*fsi]*f[(i+2)*fsi];
          sigma1 = 4.0*f[(i+0)*fsi]*f[(i+0)*fsi] - 4.0*f[(i+0)*fsi]*f[(i+1)*fsi] - 4.0*f[(i+0)*fsi]*f[(i-1)*fsi] + f[(i+1)*fsi]*f[(i+1)*fsi] + 2.0*f[(i+1)*fsi]*f[(i-1)*fsi] + f[(i-1)*fsi]*f[(i-1)*fsi];
          //--- weights ---
          omega0 = 0.33333333333333333333333333333333333/(1.0e-6*1.0e-6 + 2*1.0e-6*sigma0 + sigma0*sigma0);
          omega1 = 0.66666666666666666666666666666666667/(1.0e-6*1.0e-6 + 2*1.0e-6*sigma1 + sigma1*sigma1);
          acc = omega0 + omega1;
          omega0 = omega0/acc;
          omega1 = omega1/acc;
          //--- reconstructions ---
          fr0 = -1.5*f[(i+0)*fsi] + 2.0*f[(i+1)*fsi] - 0.5*f[(i+2)*fsi];
          fr1 = 0.5*f[(i+1)*fsi] - 0.5*f[(i-1)*fsi];
          fs0 = fr0*omega0 + fr1*omega1;
          d_p = fs0;

          f = values_buffer;
          // for bias ((-)) ---
          //--- smoothness ---
          sigma0 = f[(i+0)*fsi]*f[(i+0)*fsi] - 4.0*f[(i+0)*fsi]*f[(i+1)*fsi] + 2.0*f[(i+0)*fsi]*f[(i+2)*fsi] + 4.0*f[(i+1)*fsi]*f[(i+1)*fsi] - 4.0*f[(i+1)*fsi]*f[(i+2)*fsi] + f[(i+2)*fsi]*f[(i+2)*fsi];
          sigma1 = 4.0*f[(i+0)*fsi]*f[(i+0)*fsi] - 4.0*f[(i+0)*fsi]*f[(i+1)*fsi] - 4.0*f[(i+0)*fsi]*f[(i-1)*fsi] + f[(i+1)*fsi]*f[(i+1)*fsi] + 2.0*f[(i+1)*fsi]*f[(i-1)*fsi] + f[(i-1)*fsi]*f[(i-1)*fsi];
          //--- weights ---
          omega0 = 0.66666666666666666666666666666666667/(1.0e-6*1.0e-6 + 2*1.0e-6*sigma0 + sigma0*sigma0);
          omega1 = 0.33333333333333333333333333333333333/(1.0e-6*1.0e-6 + 2*1.0e-6*sigma1 + sigma1*sigma1);
          acc = omega0 + omega1;
          omega0 = omega0/acc;
          omega1 = omega1/acc;
          //--- reconstructions ---
          fr0 = -0.5*f[(i+0)*fsi] + 0.5*f[(i+2)*fsi];
          fr1 = -2.0*f[(i+0)*fsi] + 1.5*f[(i+1)*fsi] + 0.5*f[(i-1)*fsi];
          fs0 = fr0*omega0 + fr1*omega1;
          d_m = fs0;
        }
        else if (gradient_approx_method == 1)
        {
          float values_buffer[3];
          float *f = values_buffer+1;
          for (int i=-1; i<=1; ++i) f[i] = phi(add_to_axis(p, axis, i));
          for (int k1=-1; k1<=0; ++k1)
          {
            float d1 = (f[k1+1]-f[k1]);
            if (k1 == -1)
              d_m = d1;
            else
              d_p = d1;
          }
        }
        
        D_p += std::max(my::sqr(std::max(d_m, 0.f)), my::sqr(std::min(d_p,0.f)));
        D_m += std::max(my::sqr(std::min(d_m, 0.f)), my::sqr(std::max(d_p,0.f)));
      }
      float D = std::sqrt(phi_sign(p)>0 ? D_p : D_m);
      res(p) = factor*D;
    }
  }

#ifdef USE_CONSTRAINT_LEVELSET_METHOD
  void initSussmannConstraint()
  {
    integration_stencil = initIntegrationStencil(dim);
    sussman_constraint_const_factor.initFromBox(ld.Box());
    #pragma omp parallel for schedule(dynamic, 1)
    for (int box_index=0; box_index<mtboxes.size(); ++box_index)
    {
      if (!mtboxes_active[box_index]) continue;
      calcSussmanConstraintConstFactor(mtboxes[box_index]);
    }
  }

  void calcSussmanConstraintConstFactor(const BBox3 &bbox)
  {
    BBox3 bbox_ext = ExtendForDim(bbox, dim, 1);
    float spacing = ld.Scale();
    ConstArray3d<float> phi0 = init_phi;

    Array3d<float> arr_f2(bbox_ext),
                   arr_lambda(bbox), arr_grad(bbox_ext);

    calcGradientNorm(BBox3(bbox_ext).Clip(ld.Box()), phi0, phi0, arr_grad);

    FOR_BBOX3(p, bbox_ext)
    {
      arr_f2(p) = my::sqr(smooth_delta(phi0(p), spacing))*arr_grad(p);
    }

    // denom in lambda
    convolution(arr_f2[bbox_ext], integration_stencil, arr_lambda[bbox]); // use as temporary

    FOR_BBOX3(p, bbox)
    {
      if(active(p) == SITE_INACTIVE) continue;
      sussman_constraint_const_factor(p) = 1./(arr_lambda(p)+1.e-13) * smooth_delta(phi0(p), spacing) * arr_grad(p);
    }
  }
  
/*
 * this is lambda_ij * H'(phi_n+1)|grad phi_n+1|
 */
  void calcSussmanConstraint(const BBox3 &bbox, float dt, ConstArray3d<float> phi1, Array3d<float> res) const
  {
    ConstArray3d<float> phi0 = init_phi;
    BBox3 bbox_ext = ExtendForDim(bbox, dim, 1);
    float spacing = ld.Scale();
    
    Array3d<float> arr_f1(bbox_ext);

    FOR_BBOX3(p, bbox_ext)
    {
      arr_f1(p) = -smooth_delta(phi0(p), spacing)*(phi1(p)-phi0(p))/dt;
    }

    // nominator in lambda
    convolution(arr_f1[bbox_ext], integration_stencil, res[bbox]);

    FOR_BBOX3(p, bbox)
    {
      if(active(p) == SITE_INACTIVE) continue;
      res(p) *= sussman_constraint_const_factor(p);
    }
  }
#endif

  void calcSlope(const State &state, State &slope, StepControl &ctrl)
  {
    float spacing = ld.Scale();
    
    if (slope.phi.empty()) // use fresh ops with all boxes (not just active)
      LevelsetOps(*mtboxes, ld.Box(), dim, 0).init(slope, true); // zero out the memory (for safety)

    //CopyBorder(const_cast<State&>(state).phi[ld.Box()], dim, 2);

    #pragma omp parallel
    {
      Array3d<float> grad_buffer;     
      BOOST_FOREACH(const BBox3 &bbox, mtboxes_active_list.getCurrentThreadRange())
      {
        CopyBorder(const_cast<Array3d<float>&>(state.phi), bbox, ld.Box(), dim, 3);
        
        grad_buffer.initFromBox(bbox);
        calcGradientNorm(bbox, state.phi, state.phi, grad_buffer);
        FOR_BBOX3(p, bbox)
        {
          if (active(p) == SITE_INACTIVE) continue;

          float v = state.phi(p);
          float sign = v / std::sqrt(v*v + my::sqr(spacing*grad_buffer(p)));
          float locdphi = - sign * (grad_buffer(p) - 1.);

          slope.phi(p) = locdphi;
        }
      }
    }

    // factor in front of 0.8 for accuracy
    ctrl.dt = ctrl.euler_dt.min(std::min(ctrl.dt, 0.25 * 0.8*spacing/(dim)));
  }


  bool doStep(State &state, StepControl &ctrl)
  {
    return stepper.doStep(state, ctrl);
  }
  
  
  void postStep(State &state, const StepControl &ctrl)
  {
#ifdef USE_CONSTRAINT_LEVELSET_METHOD
    if (!do_full_restart)
    {
      if (sussman_constraint_const_factor.empty())
        initSussmannConstraint();
      #pragma omp parallel
      {
        Array3d<float> buffer;
        #pragma omp for schedule(dynamic, 1)
        for (int box_index=0; box_index<mtboxes.size(); ++box_index)
        {
          if (!mtboxes_active[box_index]) continue;
          const BBox3 &bbox = mtboxes[box_index];
          buffer.initFromBox(bbox);
          calcSussmanConstraint(bbox, ctrl.dt, state.phi, buffer);
          state.phi[bbox].addScaled(ctrl.dt, buffer[bbox]);
        }
      }
      CopyBorder(state.phi[ld.Box()], dim, 2);
    }
#endif
  }

  void appendToImages(State &state, std::vector<Image> &images)
  {
    Image img;
    state.Draw(img, ld, DrawArrayOpts().scaled(1., 0).outputRange(), &active);
    images.push_back(img);
    State slope;
    StepControl ctrl;
    ctrl.dt = 0.125;
    calcSlope(state, slope, ctrl);
    DrawArray(img, slope.phi[ld.Box()], DrawArrayOpts().outputRange());
    images.push_back(img);
#ifdef USE_CONSTRAINT_LEVELSET_METHOD
    if (sussman_constraint_const_factor.empty())
      initSussmannConstraint();
    Array3d<float> buffer(ld.Box());
    calcSussmanConstraint(ld.Box(), ctrl.dt, state.phi, buffer);
    DrawArray(img, buffer, DrawArrayOpts().outputRange());
    images.push_back(img);
#endif
  }

  void writeDebugImage(const string &out_fn, State &state)
  {
    std::vector<Image> images;
    appendToImages(state, images);
    Image bigimg;
    DrawImageGrid(bigimg, images);
    bigimg.Write(out_fn+".png");
  }

  void writeDebugH5(const string &out_fn, State &state, int iter)
  {
//     h5cpp::File f(out_fn+".h5", iter==0 ? "w" : "a");
//     h5cpp::Group ld_group = RequireLatticeDataGroup(f.root(),"field_ld", ld);
//     h5cpp::Group g = f.root().create_group(str(format("out%04i") % iter));
//     g.attrs().set("time", iter);
//     WriteScalarField(g, "ls", state.phi[ld.Box()], ld, ld_group);
//     WriteArray3D<float>(g, "ls_full", state.phi[ExtendForDim(ld.Box(), dim, 3)]);
//     if (iter == 0)
//       WriteArray3D<char>(f.root(), "active", active[ld.Box()]);
    H5::H5File f;
    if ( iter==0)
    {
      f = H5::H5File(out_fn+".h5", H5F_ACC_TRUNC);
    }
    else
    {
      f = H5::H5File(out_fn+".h5", H5F_ACC_RDWR);
    }
    //H5::H5File f = H5::H5File(out_fn+".h5", iter==0 ? "w" : "a");
    H5::Group ld_group = RequireLatticeDataGroup(f,"field_ld", ld);
    H5::Group g = f.createGroup(str(format("out%04i") % iter));
    //H5::Group g = f.root().create_group(str(format("out%04i") % iter));
    //H5::Attribute attr = g.createAttribute("time");
    writeAttrToGroup<int>(g, "time", iter);
    
    //g.attrs().set("time", iter);
    
    WriteScalarField(g, "ls", state.phi[ld.Box()], ld, ld_group);
    WriteArray3D<float>(g, "ls_full", state.phi[ExtendForDim(ld.Box(), dim, 3)]);
    if (iter == 0)
      WriteArray3D<char>(f.openGroup("/"), "active", active[ld.Box()]);
  }
};


/*
reinit levelset to a signed distance map
parameters:
fn_out, output filename

*/
ptree reinit_levelset(Levelset &levelset, const LatticeDataQuad3d &ld, int dim, const DomainDecomposition &mtboxes, int width, const ptree &params)
{
  string out_fn = params.get<string>("fn_out", string());
  //string out_fn = "levelset_debug.h5";

  LevelSetReinit::State state;
  LevelsetOps(mtboxes, ld.Box(), dim, 3).init(state, false);

  #pragma omp parallel
  {
    BOOST_FOREACH(const BBox3 &box, mtboxes.getCurrentThreadRange())
    {
      state.phi[box].fill(levelset.phi[box]);
      CopyBorder(state.phi, box, ld.Box(), dim, 3);
    }
  }
  
  LevelSetReinit reinit;
  reinit.init(ld, dim, mtboxes, state, params);
  reinit.fillActiveByWidth(state.phi, width);

  if (!out_fn.empty())
    reinit.writeDebugImage(out_fn+str(format("_%04i") % 0), state);

  NewSteppers::StepControl ctrl;
  int iter = 0;
  while(ctrl.t < (width+2)*ld.Scale())
  {
    ctrl.dt = 100000.;
    ctrl.euler_dt = my::Smallest<double>(std::numeric_limits<double>::max());
    
    reinit.doStep(state, ctrl);

    if (!out_fn.empty())
    {
      cout << format("reinit %i dt=%f, t=%f") % iter % ctrl.dt % (ctrl.t/ld.Scale()) << endl;
      reinit.writeDebugImage(out_fn+str(format("_%04i") % iter), state);
      reinit.writeDebugH5(out_fn, state, iter);
    }
    ++iter;
  }

  reinit.finalize(state, width*ld.Scale());

  if (!out_fn.empty())
  {
    reinit.writeDebugImage(out_fn+str(format("_%04i") % iter), state);
    reinit.writeDebugH5(out_fn, state, iter);
  }

  ptree res;
  bool return_grad_error = params.get<bool>("return_grad_error", false) || !out_fn.empty();
  if (return_grad_error)
  {
    reinit.fillActiveByValue(levelset.phi, state.max_val*0.75);
    double gnorm = reinit.globalGradientDeviation(levelset.phi);
    res.put("grad_error", gnorm);
    cout << format("gradient deviation: %f, max val=%f") % res.get<double>("grad_error") % state.max_val << endl;
  }

  #pragma omp parallel
  {
    BOOST_FOREACH(const BBox3 &box, mtboxes.getCurrentThreadRange())
    {
      levelset.phi[box].fill(state.phi[box]);
    }
  }
  levelset.max_val = state.max_val;
  return res;
}



ptree reinit_levelset(Levelset &levelset, const LatticeDataQuad3d &ld, int dim, int width, const ptree &params)
{
  DomainDecomposition mtboxes(MakeMtBoxGrid(ld.Box(), Int3(16, 16, 16)));
  return reinit_levelset(levelset, ld, dim, mtboxes, width, params);
}



double calc_levelset_badness(const Levelset &levelset, const LatticeDataQuad3d &ld, int dim, const DomainDecomposition &mtboxes)
{
  LevelSetReinit reinit;
  reinit.init(ld, dim, mtboxes, levelset, ptree());
  reinit.fillActiveByValue(levelset.phi, levelset.max_val*0.75);
  return reinit.globalGradientDeviation(levelset.phi);
}

double calc_levelset_badness(const Levelset &levelset, const LatticeDataQuad3d &ld, int dim)
{
  DomainDecomposition mtboxes(MakeMtBoxGrid(ld.Box(), Int3(16, 16, 16)));
  return calc_levelset_badness(levelset, ld, dim, mtboxes);
}


void AddScaled(double f0, Levelset& l0, double f1, const Levelset &l1)
{
  AddScaled(f0, l0.phi, f1, l1.phi);
}

double levelset_ok_badness()
{
  return 0.005;
}
