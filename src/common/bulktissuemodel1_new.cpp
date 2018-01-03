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

#if 1
#include "bulktissuemodel1_new.h"
#include "trilinos_linsys_construction.h"
#include "continuum-surfacetension.h"

#include "mwlib/random.h"
#include "mwlib/log.h"

#include <Sacado_Fad_SFad.hpp>






namespace NewBulkTissueModel
{


enum { NUM_COMPONENTS = 2 };
typedef Sacado::Fad::SFad<float, NUM_COMPONENTS> AdScalar;

typedef Sacado::Fad::SFad<float, 3> AdScalar3;


enum { ID_CELLS = 0, ID_THETA_NECRO = 1 };

Params::Params()
{
  tumor_diameter = 200;
  time_prol = 24;
  time_death = 240;
  time_prol_tumor = 24;
  time_death_tumor = 240;
  time_necrosis = 100.;
  ncells_norm = 0.4;
  ncells_tumor = 0.5;
  ncells_sigma = 0.1;
  ncells_ecm = 0.2;
  cell_mobility = 3600;
  cell_mobility_tumor = 3600;
  cell_mobility_ecmstar = -1; // disabled, do not use ecm dependence
  cell_mobility_ecmstarstar = -1;
  
  cell_compression_modulus = 1.;
 
  timestep_factor = 0.9;
  interface_width = -1.;
  wrong_model = false;
  sigma_model_version = 1;
  source_model_version = 1;

  o2_necro_threshold = -1.;
  o2_prol_threshold = -1.;
  use_necrotic_regions = false;
  
  ecm_noise_std = 0.05;
  ecm_noise_density = -1; // disabled
  random_seed = 123456;

  surface_tension = -1.; // disabled
}

#define PT_ASSIGN(name) boost::property_tree::get(name, #name, pt)

void Params::assign(const ptree &pt)
{
  PT_ASSIGN(tumor_diameter);
  PT_ASSIGN(time_prol);
  PT_ASSIGN(time_death);
  PT_ASSIGN(time_prol_tumor);
  PT_ASSIGN(time_death_tumor);
  PT_ASSIGN(time_necrosis);
  PT_ASSIGN(ncells_norm);
  PT_ASSIGN(ncells_tumor);
  PT_ASSIGN(ncells_sigma);
  PT_ASSIGN(ncells_ecm);
  PT_ASSIGN(cell_mobility);
  PT_ASSIGN(cell_mobility_tumor);
  PT_ASSIGN(cell_mobility_ecmstar);
  PT_ASSIGN(cell_mobility_ecmstarstar);
  //PT_ASSIGN(velocity_smoothing_factor);
  PT_ASSIGN(timestep_factor);
  PT_ASSIGN(interface_width);
  PT_ASSIGN(wrong_model);
  PT_ASSIGN(ecm_noise_std);
  PT_ASSIGN(ecm_noise_density);
  PT_ASSIGN(random_seed);
  PT_ASSIGN(sigma_model_version);
  PT_ASSIGN(source_model_version);
  PT_ASSIGN(o2_necro_threshold);
  PT_ASSIGN(o2_prol_threshold);
  PT_ASSIGN(use_necrotic_regions);
  PT_ASSIGN(surface_tension);
  PT_ASSIGN(cell_compression_modulus);
  //PT_ASSIGN();
}

#define AS_PTREE(name) pt.put(#name, name);

ptree Params::as_ptree() const
{
  ptree pt;
  AS_PTREE(tumor_diameter);
  AS_PTREE(time_prol);
  AS_PTREE(time_death);
  AS_PTREE(time_prol_tumor);
  AS_PTREE(time_death_tumor);
  AS_PTREE(time_necrosis);
  AS_PTREE(ncells_norm);
  AS_PTREE(ncells_tumor);
  AS_PTREE(ncells_sigma);
  AS_PTREE(ncells_ecm);
  AS_PTREE(cell_mobility);
  AS_PTREE(cell_mobility_tumor);
  AS_PTREE(cell_mobility_ecmstar);
  AS_PTREE(cell_mobility_ecmstarstar);
  //AS_PTREE(velocity_smoothing_factor);
  AS_PTREE(timestep_factor);
  pt.put("write_levelset_function", false);
  pt.put("write_face_velocities", false);
  pt.put("write_surface_tension_force", false);
  AS_PTREE(interface_width);
  AS_PTREE(wrong_model);
  AS_PTREE(ecm_noise_std);
  AS_PTREE(ecm_noise_density);
  AS_PTREE(random_seed);
  AS_PTREE(sigma_model_version);
  AS_PTREE(source_model_version);
  AS_PTREE(o2_necro_threshold);
  AS_PTREE(o2_prol_threshold);
  AS_PTREE(use_necrotic_regions);
  AS_PTREE(surface_tension);
  AS_PTREE(cell_compression_modulus);
  
  pt.put("elliptic_solver_params.max_iter", 500);
  pt.put("elliptic_solver_params.max_resid", 1.e-7);
  pt.put("elliptic_solver_params.output", 1);
  pt.put("elliptic_solver_params.use_multigrid", false);

  pt.put("init_shape", "sphere");

  pt.put("levelset_stepper", "rk3");
  pt.put("mixture_stepper", "rk3");
  //AS_PTREE();
  return pt;
}


void Params::init(const LatticeDataQuad3d &ld)
{
  psi_norm = ncells_norm + ncells_ecm;
  psi_norm_tumor = ncells_tumor + ncells_ecm;
  if (interface_width<0.)
    interface_width = ld.Scale();
}


/**
@brief A nearly linear decrease from f(0)=1 to f(sigma) = 0. The kinks are smoothed within a range of 0.1 sigma.
*/
template<class T>
inline T smoothed_clamped_linear_ramp(T x, T sigma)
{
  myAssert(sigma > 0.);
  x /= sigma;
  T f0 = my::smooth_heaviside<T>(x, 0.1);
  T f1 = my::smooth_heaviside<T>(x - 1, 0.1);
  x = (1.-x) * f0 + 1. * (1.-f0);
  x = 0. * f1 + x * (1.-f1);
  return x;
}


const string key_ls = "ls",
             key_cells = "phi_cells",
             key_necro = "phi_necro",
             key_other = "phi_other",
             key_oxy   = "oxy",
             key_checksum = "bulktissue_state_id";


void Model::initStateVars(State &state, bool clear) const
{
  ops_levelset.init(state.ls, clear);
  ops_array3d.init(state.cells, clear);
  ops_array3d.init(state.necro, clear);
}


void Model::cloneState(Model::State& u, const Model::State& v)
{
  ops_levelset.initFrom(u.ls, v.ls, ConsMode::AS_COPY);
  ops_array3d.initFrom(u.cells, v.cells, ConsMode::AS_COPY);
  ops_array3d.initFrom(u.necro, v.necro, ConsMode::AS_COPY);
  u.checksum = v.checksum;
}


void Model::init(const ContinuumGrid &grid_, const DomainDecomposition &mtboxes_,
            Model::State& state,
            DataRequestCallback obtain_vessel_volume_,
            DataRequestCallback obtain_oxygen_conc_,
            const ptree& pt_params_)
{
  grid = &grid_;
  mtboxes = &mtboxes_;
  obtain_oxygen_conc = obtain_oxygen_conc_;
  obtain_vessel_volume = obtain_vessel_volume_;
  
  cout << format("going to use a mt grid subdivision with %i boxes and %i threads") % mtboxes->size() % my::GetNumThreads() << endl;
  
  pt_params = pt_params_;
  params.assign(pt_params);
  params.init(grid->ld);
  
  last_levelset_reinit = -std::numeric_limits<double>::max();

  max_dt_diff = max_dt_src = max_dt_vel = max_dt_stf = my::MakeSmallestMax<double>();
  last_euler_dt = std::numeric_limits<double>::max();
  
  ops_array3d.init(*mtboxes, grid->Box(), grid->dim, 3);
  ops_levelset.init(*mtboxes, grid->Box(), grid->dim, 3);

  notify_auxdata_change();
  
  { // create the initial state
    initStateVars(state, true);
    auto ncells = state.cells;
    auto &ls    = state.ls;
    // necrotic fraction is already filled with zeros
    
    const string init_shape = pt_params.get<string>("init_shape");
    bool init_cylindric = init_shape == "cylinder";
    bool init_rect      = init_shape == "rectangle";
    double init_rad = params.tumor_diameter * 0.5;

    // initialize the cell distribution and the levelset function with
    // a tumor nucleus in the center.
    #pragma omp parallel
    {
      Random rnd;
     
      BOOST_FOREACH(const BBox3 bb, mtboxes->getCurrentThreadRange())
      {
        FOR_BBOX3(p, bb)
        {
          Float3 wp = grid->ld.LatticeToWorld(p);
          if (init_cylindric) wp[2] = 0.;
          if (init_rect) { wp[2] = wp[1] = 0.; }

          float f = my::smooth_heaviside(wp.norm() - init_rad, grid->ld.Scale() * 2.);

          ncells(p) = params.ncells_norm;
      #ifdef DEBUG  // random noise to attempt to drive the system out of the steady state. Used to detect early if the system is unstable.
          ncells(p) += rnd.Get11()*0.01*0.45;
      #endif

          float excess = std::max(0., ncells(p) + other_volume(p) - 0.95);
          ncells(p) = std::max<float>(0., ncells(p) - excess);

          ls.phi(p) = init_rad - wp.norm();
        }

      }
    }
    reinit_levelset(ls, grid->ld, grid->dim, *mtboxes, 10);
  }
  
  // Tell the time steppers how to compute the slope (or derivative) function.
  // I.E. want to compute F(u,t) in du/dt = F(u,t).
  // Also tell the time steppers how to perform mathematical operations on our system state.
  stepper_convection.model().calcSlope = boost::bind(&Model::calcSlopeConvection, this, _1, _2, _3);
  stepper_convection.ops().addScaled = [=](double fu, State &u, double fv, const State &v)
  {
    this->ops_array3d.addScaled(fu, u.cells, fv, v.cells);
    this->ops_array3d.addScaled(fu, u.necro, fv, v.necro);
    this->ops_levelset.addScaled(fu, u.ls, fv, v.ls);
  };
  stepper_convection.ops().initFrom = [=](State &u, const State &v, ConsMode consMode)
  {
    this->ops_array3d.initFrom(u.cells, v.cells, consMode);
    this->ops_array3d.initFrom(u.necro, v.necro, consMode);
    this->ops_levelset.initFrom(u.ls, v.ls, consMode);
  };

  stepper_mixture.model().calcSlope = boost::bind(&Model::calcSlopeMixture, this, _1, _2, _3);
  stepper_mixture.ops().addScaled = [=](double fu, State &u, double fv, const State &v)
  {
    this->ops_array3d.addScaled(fu, u.cells, fv, v.cells);
    this->ops_array3d.addScaled(fu, u.necro, fv, v.necro);
  };
  stepper_mixture.ops().initFrom = [=](State &u, const State &v, ConsMode consMode)
  {
    u.ls = v.ls; // shallow copy
    this->ops_array3d.initFrom(u.cells, v.cells, consMode);
    this->ops_array3d.initFrom(u.necro, v.necro, consMode);
  };
#if 0
  stepper_mixture.model().calcSlopesIMEX         = boost::bind(&Model::calcSlopesIMEXMixture, this, _1, _2, _3, _4);
  stepper_mixture.model().invertImplicitOperator = boost::bind(&Model::invertImplicitOperatorMixture, this, _1, _2, _3, _4, _5, _6);
#endif
}



void Model::notify_auxdata_change()
{
  ops_array3d.init(other_volume, false);
  const bool has_vessels = obtain_vessel_volume != NULL;

  tbb::spin_mutex mutex;
  
  #pragma omp parallel
  {
    Random rnd(params.random_seed * (my::OmpGetCurrentThread()+1)); // might need it
    Array3df buffer(ExtendForDim(mtboxes->getBoxWithLargestVolume(), grid->dim, 1));
    Array3df vessel_volume_field;
    
    BOOST_FOREACH(const BBox3 bb, mtboxes->getCurrentThreadRange())
    {
      if (!has_vessels && params.ecm_noise_density <=0.)
      {
        other_volume[bb].fill(params.ncells_ecm);
      }
      else
      {
        if (params.ecm_noise_density > 0.)
        {
          BBox3 extbox = ExtendForDim(bb, grid->dim, 1);
          buffer.reshape(extbox); // assumes that AddSmoothDelta has 2 sites support
          buffer.fill(0.);
          
          const double num_noise = std::pow(grid->ld.Scale()*params.ecm_noise_density, grid->dim);

          FOR_BBOX3(p, bb)
          {
            double tmp = num_noise;
            while (rnd.Get01() < tmp)
            {
              float val = rnd.Get11f()*0.5*params.ecm_noise_std;
              Float3 pos = grid->ld.LatticeToWorld(p) + 0.49*grid->ld.Scale()*Float3(rnd.Get11f(), rnd.Get11f(), rnd.Get11f());
              AddSmoothDelta(buffer, grid->ir.cells, grid->ld, grid->dim, pos, val);
              tmp -= 1.;
            }
          }

          // here we have to sync because extbox overlaps with other threads and
          // we have to add contributions which go over the local threads box!
          // The lock is for the whole region which is not nessesary but easily programmed.
          mutex.lock();
          other_volume[extbox] += buffer[extbox];
          mutex.unlock();
        }

        if(has_vessels)
          vessel_volume_field = obtain_vessel_volume(bb);
        FOR_BBOX3(p, bb)
        {
          float ecm = other_volume(p);
          ecm += params.ncells_ecm;
          if (has_vessels)
            ecm = std::max(ecm, vessel_volume_field(p));
          ecm = my::cut<float>(ecm, 0., 0.95);
          other_volume(p) = ecm;
        }
      }
      CopyBorder(other_volume, bb, grid->Box(), grid->dim, 1);
    }
  }
}


/**
@brief computes sigma(psi) = solid pressure
* psi = phi = cell volume fraction (tumor or normal)
* psi_other = phi_other = vessels+ecm
* psi0 = threshold above which cells feel pressure, otherwise not pressure.
* 
* attempt was made to respect excluded volume by vessels. Never worked to my satisfaction.
* derivatives are automatically available by automatic differencing
* WARNING: this term MUST be positive for the diffusion coefficient to be positive!!!

Adscalar comes from trilinos. It is a special number type that allows for automatic computation of mathematical
derivatives of expression.
*/
static AdScalar calcSigma(const AdScalar &phi_cells, const AdScalar &phi_other, const Params &params)
{
//   if (params.sigma_model_version == 3)
//   {
//     // AdScalar phi_norm = (1.-theta) * phi_cells;
//     // AdScalar phi_tum = theta * phi_cells;
//     // AdScalar phi_obst = phi_other;
//     //AdScalar ncells_norm = params.ncells_norm - std::max<AdScalar>(0., params.ncells_norm + phi_obst - 1.);
//     //AdScalar phi_total = phi_norm + phi_tum + phi_obst; //ncells_norm + phi_obst;
//     //AdScalar s = (phi_norm * phi_norm / ncells_norm + phi_tum * phi_tum / ncells_norm + phi_obst) / phi_total;
//     //s = my::smooth_heaviside<AdScalar>(s - 0.9, 0.1) * (s-0.9) + 0.9;
//     
//     s *= params.ncells_compression_modulus;
//     return s;
//   }
//   else
  if (params.sigma_model_version == 2)
  {
    const float psi0 = params.ncells_norm + params.ncells_ecm;
    AdScalar psi00 = my::smooth_max<AdScalar>(0., psi0 - phi_other, 0.05);
    return my::smooth_max<AdScalar>(0., params.cell_compression_modulus * (phi_cells - psi00), 0.05);
  }
  else
  {
    // here the baseline pressure is defined without ecm. Therefore the density where proliferation is zero
    // which is params.ncells_norm does not correspond to zero Sigma any more.
    const float psi0 = params.ncells_norm;
    AdScalar psi = phi_cells + phi_other;
  #if 1 // linear
    return std::max<AdScalar>(0., params.cell_compression_modulus * (psi - psi0));
  #else
    if (psi > psi0)
      return params.cell_compression_modulus * (1.-psi0) * (psi - psi0) / (1. - psi + 1.e-3);
    else
      return 0;
  #endif
  }
}


float CalcSolidPressure(float phi_cells, float phi_other, const Params &params)
{
  AdScalar loc_phi_cells(NUM_COMPONENTS, 0, phi_cells);
  AdScalar loc_phi_other(NUM_COMPONENTS, 1, phi_other);
  return calcSigma(loc_phi_cells, loc_phi_other, params).val();
}


/**
@brief computes phi * sigma(psi)
  * derivatives are automatically available by automatic differencing
  * phi = cell volume fraction (tumor or normal)
  * phi_other = vessels+ecm
  */
static AdScalar calcCellVelocityPotential(float theta, const AdScalar &phi_cells, const AdScalar &phi_other, const AdScalar &sigma, const Params &params)
{
  AdScalar result = sigma;
  if (!params.wrong_model)
    result *= phi_cells;  // total force is average over all species within the unit cell. Therefore phi_cells (volume fraction) appears as weighting factor (?) See BulkTissue papers.
  return result;
}


/**
@brief Compute mobility coefficient, which is the proportionality constant between velocity and force.
  Allows for nonlinear behavior where a certain force threshold must be overcome until cells start to move (normally not used).
*/
static AdScalar calcCellMobility(float theta, const AdScalar &phi_cells, const AdScalar &phi_other, const Params &params)
{
  AdScalar kstar = (params.cell_mobility * (1 - theta) + params.cell_mobility_tumor * theta);
  
  if (params.cell_mobility_ecmstar <= 0.)
  {
    myAssert(kstar.val() >= 0.);
    return kstar;
  }
  else
  {
    float a = params.cell_mobility_ecmstar;
    float b = params.cell_mobility_ecmstarstar;
    if (b < 0) b = 1.;
    
    AdScalar f = (phi_other-a)/(b - a);
    kstar *= 1. - my::clamp<AdScalar>(f, 0., 1.);

    if (!params.wrong_model)
    {
      kstar /= (phi_cells + 1.e-13);
    }
    
    myAssert(kstar.val() >= 0.);
    return kstar;
  }
}

/**
@brief Compute cell volume fraction inserted into the system per time and unit volume.

  There are several variants, none of which i remember how exactly they work. An important principle is
  growth inhibition by limited availability of space (or due to solid pressure). Also the growth rate is
  modulated by oxygen (o2) levels. It can become negative if oxygen is too low.   
  I also tried to incorporate exclusion of space due to blood vessels, with limited success.
*/
static void calcSource(float theta, const AdScalar3 &phi_cells, const AdScalar3 &phi_other, const AdScalar3 &phi_necro, float o2, const Params &params, AdScalar3 *sources)
{
  // the sources
  if (params.source_model_version == 1)
  {
    float tc_death_rate = 1./params.time_death_tumor;
    float norm_death_rate = 1./params.time_death;
    float tc_prol_rate = 1./params.time_prol_tumor;
    float norm_prol_rate = 1./params.time_prol;
    float tc_density_thres = params.psi_norm_tumor;
    float norm_density_thres = params.psi_norm;
    float sigma_density = params.ncells_sigma;
    AdScalar3 psi = phi_cells + phi_other;
    AdScalar3 q_tum = tc_prol_rate*smoothed_clamped_linear_ramp<AdScalar>(psi - tc_density_thres, sigma_density)*o2 - tc_death_rate;
    AdScalar3 q_norm = norm_prol_rate*smoothed_clamped_linear_ramp<AdScalar3>(psi - norm_density_thres, sigma_density)*o2 - norm_death_rate;
    AdScalar3 s  = phi_cells * (theta*q_tum + (1-theta)*q_norm);
    sources[ID_CELLS] = s;
  }

  else if(params.source_model_version == 2)
  {
    std::cerr << "bulktissue-tumor-model: params.source_model_version == 2 not implemented!" << endl;
    exit(0);
//     enum { TUMOR = 0, NORMAL = 1 };
//     float death_rates[2] = { 1./params.time_death_tumor, 1./params.time_death };
//     float prol_rates[2] = { 1./params.time_prol_tumor, 1./params.time_prol };
//     float prol_density_thres[2] = { params.psi_norm_tumor, params.psi_norm };
//     float sigma_density = params.ncells_sigma;
//     float necrosis_rates[2] = { 1./params.time_necrosis, 1./params.time_necrosis };
//     AdScalar3 psi = phi_cells + phi_other;
//     AdScalar3 q_species[2], q_species_necro[2];
// 
//     float f_o2_prol  = my::smooth_heaviside<float>(o2 - params.o2_prol_threshold, params.o2_prol_threshold*0.1);
//     float f_o2_necro  = my::smooth_heaviside<float>(o2 - params.o2_necro_threshold, params.o2_necro_threshold*0.1);
//     
//     for (int j=0; j<2; ++j)
//     {
//       q_species[j] = prol_rates[j]*smoothed_clamped_linear_ramp<AdScalar3>(psi - prol_density_thres[j], sigma_density)*f_o2_prol - death_rates[j]*f_o2_necro;  // - necrosis_rates[j]*(1.-f_o2_necro) 
//       q_species[j] *= (1. - phi_frac_necro); // only the viable fraction of the cells can change
//       q_species_necro[j] = (1. - phi_frac_necro)*necrosis_rates[j]*(1.-f_o2_necro); // this is the fraction of necrotic tissue generated
//     }
//     sources[ID_CELLS]  = phi_cells * (theta*q_species[0] + (1-theta)*q_species[1]);
//     sources[ID_THETA_NECRO] = phi_cells * theta*q_species_necro[0] + (1.-theta)*q_species_necro[1]; // now we have the volume fraction change
    //sources[ID_THETA_NECRO] = (sources[ID_THETA_NECRO]*phi_cells - (phi_frac_necro*phi_cells)*sources[ID_CELLS])/(1.e-13 + my::sqr<AdScalar3>(phi_cells)); // this is the rate of change of the fraction of the volume fraction
  }
  
  else if (params.source_model_version == 3)
  {
    float death_rates[2] = { 1.f/params.time_death_tumor, 1.f/params.time_death };
    float prol_rates[2] = { 1.f/params.time_prol_tumor, 1.f/params.time_prol };
    float prol_density_thres[2] = { params.psi_norm_tumor, params.psi_norm };
    float sigma_density = params.ncells_sigma;
    float necrosis_rates[2] = { 1.f/params.time_necrosis, 1.f/params.time_necrosis };
    AdScalar3 psi = phi_cells + phi_other;
    AdScalar3 q_species[2] = { 0, 0 }, q_species_necro[2] = { 0, 0 };

    float f_o2_prol  = params.o2_prol_threshold>0.  ?  my::smooth_heaviside<float>(o2 - params.o2_prol_threshold, params.o2_prol_threshold*0.1) : 1.;
    float f_o2_necro = params.o2_necro_threshold>0. ?  my::smooth_heaviside<float>(o2 - params.o2_necro_threshold, params.o2_necro_threshold*0.1) : 1.;

    for (int j=0; j<2; ++j)
    {
      AdScalar3 phi = (phi_cells-phi_necro)*(j == 0 ? theta : (1-theta)); // j=0 = tumor
      //AdScalar3 necro = (phi_frac_necro)*phi_cells*(j == 0 ? theta : (1-theta)); // j=0 = tumor
      AdScalar3 dphi_necro = necrosis_rates[j]*(1. - f_o2_necro);
      AdScalar3 dphi_apoptosis = death_rates[j];
      AdScalar3 dphi_prol = my::clamp<AdScalar3>(prol_rates[j]*(prol_density_thres[j]-psi)/sigma_density + death_rates[j], 0., prol_rates[j]*f_o2_prol);
      q_species[j] = phi * (- dphi_apoptosis + dphi_prol);
      if (params.use_necrotic_regions)
        q_species_necro[j] = phi * dphi_necro; // we transfer stuff to the necrotic phase. q_species[j] is unchanged since it includes necrotic tissue
      else
        q_species[j] -= phi * dphi_necro; // we remove tissue, since here necrotic tissue becomes water
    }
    sources[ID_CELLS] = q_species[0] + q_species[1];
    sources[ID_THETA_NECRO] = q_species_necro[0] + q_species_necro[1]; // now we have the volume fraction change
  }
}

/// wrapper to get rid of AdScalar.
void calcSource(float theta, float phi_cells, float phi_other, float phi_necro, float o2, const Params& params, float *outrates)
{
  AdScalar3 q[2];
  calcSource(theta, AdScalar3(phi_cells), AdScalar3(phi_other), AdScalar3(phi_necro), o2, params, q);
  outrates[ID_CELLS] = q[ID_CELLS].val();
  outrates[ID_THETA_NECRO] = q[ID_THETA_NECRO].val();
}



struct ThreadLocal;

struct CellDiffusionCoeffFunctor
{
  const ThreadLocal &master;
  const int component;
  CellDiffusionCoeffFunctor(const int component, const ThreadLocal &master_) : component(component), master(master_) {}
  double operator()(const Int3 &p, const Int3 &nbp, int axis, int side) const;
};

struct FaceVarAveragingFunctor
{
  ConstArray3d<float> coefficient;
  FaceVarAveragingFunctor(ConstArray3d<float> coefficient_) : coefficient(coefficient_) {}

  float operator()(int axis, const Int3 &p) const
  {
    // p is also the grid cell index to the right of the face
    // q is to the left
    Int3 q(p); --q[axis];
    return 0.5*(coefficient(q)+coefficient(p));
  }
};


/**
@brief Compute various terms of the model equations for a sub-box of the simulation grid.

 * Comprise all variables and things needed to compute things for a sub-box within the entire the simulation grid.
 * Called from a thread, to enable parallel processing across the simulation domain (See DomainDecomposition).
 */
struct ThreadLocal
{
  BBox3 bbox, bbext;
  LatticeIndexRanges ir;
  int dim;
  const LatticeDataQuad3d &ld;
  Bool3 is_at_right_boundary;
  const Params &params;
  // following vars used to determine time step length  
  double max_vel,   // maximal velocity
	 max_diff,  // estimation of maximal diffuion coefficient (?)
	 max_source;// maximal cell distribution source rate

  // external data
  ConstArray3d<float> phi_cells, phi_other, phi_necro, o2field;
  // intermediate data caches
  Array3d<float> theta, pressure;
  Array3d<float> diff_coeff[NUM_COMPONENTS];
  Array3d<float> vel_pot, vel_coeff;
  SurfaceTensionForce<float> stforce;

  ThreadLocal(const LatticeDataQuad3d &ld_, int dim_, const Params &params_) :
    params(params_), dim(dim_), ld(ld_),
    max_vel(0.), max_diff(0.), max_source(0.)
    {}

  void reinit(const BBox3 &bbox_, const State &state, const ConstArray3d<float> &o2field_, const ConstArray3d<float> &other_volume_,  int flags)
  {
    bool must_reinit = !(bbox == bbox_);
    bbox = bbox_;
    phi_cells = state.cells;
    phi_other = other_volume_;
    phi_necro = state.necro;
    o2field = o2field_;
    const Levelset &ls = state.ls;

    ir = LatticeIndexRanges::FromCellRange(bbox, dim);
    bbext = ExtendForDim(bbox, dim, 1);
    is_at_right_boundary = cwiseGe(bbox.max, ld.Box().max);

    theta.initFromBox(bbext);
    ls.ToHeaviside(theta, bbext, params.interface_width);
    
    flags &= ~(OutState::THETA | OutState::PHI_CELLS); // optimization: return if we just need the heaviside of the levelset function
    if (flags == 0) return;

    if (must_reinit)
    {
      if (flags & OutState::PRESSURE)
        pressure.initFromBox(bbext);
      if (flags & OutState::VELOCITIES)
      {
        diff_coeff[0].initFromBox(bbext);
        diff_coeff[1].initFromBox(bbext);
        vel_pot.initFromBox(bbext);
        vel_coeff.initFromBox(bbext);
      }
    }

    if (params.surface_tension > 0) {
      // bbext is used here because computing the force on grid faces needs an extra layer of boundary data
      stforce.init(bbext, ls.phi, dim, ld.Scale());
    }

    FOR_BBOX3(p, bbext)
    {
      AdScalar loc_phi_cells(NUM_COMPONENTS, 0, phi_cells(p));
      AdScalar loc_phi_other(NUM_COMPONENTS, 1, phi_other(p));
      float loc_theta = theta(p);
      AdScalar loc_mobility = calcCellMobility(loc_theta, loc_phi_cells, loc_phi_other, params);
      AdScalar loc_sigma = calcSigma(loc_phi_cells, loc_phi_other, params); // solid pressure 
      AdScalar loc_vpot = calcCellVelocityPotential(loc_theta, loc_phi_cells, loc_phi_other, loc_sigma, params); // velocity potential depending on pressure
      // "diffusion coefficient"
      if (flags & OutState::VELOCITIES)
      {
        float c1 = diff_coeff[0](p) = (loc_mobility * loc_vpot.dx(0)).val();
        float c2 = diff_coeff[1](p) = (loc_mobility * loc_vpot.dx(1)).val();
        max_diff = my::max<double>(loc_phi_cells.val()*(c1 + c2), max_diff);
        //myAssert(c1 >= 0 && c2 >= 0);
        myAssert(c1 + c2 >= 0);
        vel_pot(p) = loc_vpot.val();
        vel_coeff(p) = loc_mobility.val();
      }
      if (flags & OutState::PRESSURE)
        pressure(p) = loc_sigma.val();
      float loc_phi_total = (loc_phi_cells.val() + loc_phi_other.val());
    }
  }

  void AddSourceContribution(Array3d<float> res_ncells, Array3d<float> res_necro)
  {
    FOR_BBOX3(p, bbox)
    {
      AdScalar3 loc_phi_cells(3, 0, phi_cells(p));
      AdScalar3 loc_phi_other(3, 1, phi_other(p));
      AdScalar3 loc_theta_necro(3, 2, params.use_necrotic_regions ? phi_necro(p) : 0.);
      float loc_theta = theta(p);
      AdScalar3 q[2];
      calcSource(loc_theta, loc_phi_cells, loc_phi_other, loc_theta_necro, o2field(p), params, q);
      max_source = std::max<double>(max_source,
                                    std::abs(q[0].dx(0)) + std::abs(q[0].dx(2)) + std::abs(q[1].dx(0)) + std::abs(q[1].dx(2)));
      res_ncells(p) += q[ID_CELLS].val();
      if (params.use_necrotic_regions)
        res_necro(p) += q[ID_THETA_NECRO].val();
    }
  }

#if 0
  void AddDiffusionContribution(int component, double prefactor, Array3d<float> res)
  {
    throw std::runtime_error("implicit method for bulktissue not implemented because surface tension cannot be done");
    double spacing = ld.Scale();
    ConstArray3d<float> phi_alpha = component==0 ? phi_cells : phi_other;
    ptree pt = make_ptree("prefactor", prefactor); //("max_kdiff", max_diff); //("stencil","isotropic");
    AddLaplacian2<float, float>(bbox, phi_alpha, CellDiffusionCoeffFunctor(component, *this), dim, spacing, res, pt);
    //max_diff = std::max(max_diff, pt.get<double>("max_kdiff"));
  }

  void AddVelocityContribution(int component, double prefactor, FaceVarArrays &velocities)
  {
    throw std::runtime_error("implicit method for bulktissue not implemented because surface tension cannot be done");
    double spacing = ld.Scale();
    ConstArray3d<float> phi_alpha = component==0 ? phi_cells : phi_other;

    ptree pt = make_ptree("geometric_coeff_mean", false)("max_vel", max_vel); //("stencil","isotropic");
    for (int axis=0; axis<dim; ++axis)
    {
      AddVelocitiesFromPotential(bbox, axis, dim, spacing, diff_coeff[component], phi_other, velocities[axis], pt);
    }
    //note: max_vel is obtained from contributions from phi_cells and phi_other individually!!!
    max_vel = std::max(max_vel, pt.get<double>("max_vel"));
  }

  void AddDiffusionContribution(double prefactor, FiniteVolumeMatrixBuilder &mb)
  {
    throw std::runtime_error("implicit method for bulktissue not implemented because surface tension cannot be done");
    mb.AddDiffusion2(bbox, dim, CellDiffusionCoeffFunctor(0, *this), prefactor, false); //true);
  }
#endif
  
  void SetSurfaceTensionForce(FaceVarArrays &stf)
  {
    if (params.surface_tension > 0.)
    {
      for (int axis=0; axis<dim; ++axis)
      {
        FOR_BBOX3(p, ir.faces[axis])
        {
          stf[axis](p) = -params.surface_tension * stforce.computeFaceForce(axis, p);
        }
      }
    }
  }

  void AddVelocities(FaceVarArrays &velocities, bool include_all_borders = true)
  {
    double spacing = ld.Scale();
    for (int axis=0; axis<dim; ++axis)
    {
      BBox3 bf = ir.faces[axis];
      if (!include_all_borders)
        bf.max[axis]--;

      FaceVarAveragingFunctor coeff_func(vel_coeff);
      FOR_BBOX3(p, bf)
      {
        Int3 q(p); q[axis]--;
        float v = -(vel_pot(p)-vel_pot(q))/spacing;
        if (params.surface_tension > 0.)
          v -= params.surface_tension * stforce.computeFaceForce(axis, p);
        v *= coeff_func(axis, p);
        max_vel = std::max<double>(max_vel, std::abs(v));
        velocities[axis](p) += v;
      }
    }
  }

};


double CellDiffusionCoeffFunctor::operator()(const Int3& p, const Int3& q, int axis, int side) const
{
  //Int3 q(p); --q[axis];
  float k0 = master.diff_coeff[component](p);
  float k1 = master.diff_coeff[component](q);
  float phi0 = master.phi_cells(p);
  float phi1 = master.phi_cells(q);
  float c = 0.5*(k1+k0); //CalcInterfaceCondCoeff(k0, k1);
  return c*(phi0+phi1)*0.5;
}



/**
  *@brief Calculate the advection term.
  *
  * a new implementation with different coefficients, using simply
  * div KE grad (phi sigma(Psi))
  */
void Model::calcSlopeConvection(const State& state_, State& slope, NewSteppers::StepControl& ctrl)
{
  initStateVars(slope, true);
  State &state = const_cast<State&>(state_); // for the external scope the state will look unchanged. However we need to modify 'ghost' cells around the boundary so differential operators can be computed.

  #pragma omp parallel
  {
    BOOST_FOREACH(const BBox3 bbox, mtboxes->getCurrentThreadRange())
    {
      CopyBorder(state.ls.phi, bbox, grid->Box(), grid->dim, 3);
      CopyBorder(state.cells, bbox, grid->Box(), grid->dim, 3);
      CopyBorder(state.necro, bbox, grid->Box(), grid->dim, 3);
    }

    #pragma omp barrier

    ThreadLocal terms(grid->ld, grid->dim, params); // every thread has its own terms variable, reinitializing does NOT clear the max_xxx variables.
    FaceVarArrays velocities;

    BOOST_FOREACH(const BBox3 bbox, mtboxes->getCurrentThreadRange())
    {
      velocities.init(bbox, grid->dim);
      terms.reinit(bbox, state, Array3df(), other_volume, OutState::VELOCITIES);
      terms.AddVelocities(velocities);

      AddKTAdvection(bbox, state.ls.phi, velocities.data(), grid->dim, grid->ld.Scale(), slope.ls.phi, NON_CONSERVATIVE);
      AddKTAdvection(bbox, state.cells, velocities.data(), grid->dim, grid->ld.Scale(), slope.cells, CONSERVATIVE);
      if (params.use_necrotic_regions)
        AddKTAdvection(bbox, state.necro, velocities.data(), grid->dim, grid->ld.Scale(), slope.necro, CONSERVATIVE);
    }
    #pragma omp critical
    {
//       max_dt_vel = std::min(max_dt_vel, 0.8 * 0.5*grid->ld.Scale()/(grid->dim*terms.max_vel+1.e-13));
//       max_dt_diff = std::min(max_dt_diff, 0.8 * 0.5*my::sqr(grid->ld.Scale())/(grid->dim*terms.max_diff+1.e-13));
//       max_dt_stf = (params.surface_tension > 0.) ?
//                       std::min(max_dt_stf, 2./(terms.max_diff*grid->dim*params.surface_tension+1.e-13)*my::cubed(grid->ld.Scale())) :
//                       std::numeric_limits<double>::max();

      max_dt_vel.min(EulerDt::Convection(terms.max_vel, grid->dim, grid->ld.Scale()));
      max_dt_diff.min(EulerDt::Diffusion(terms.max_diff, grid->dim, grid->ld.Scale()));
      if (params.surface_tension > 0.)
        max_dt_stf.min(EulerDt::ByOrder(0.25*terms.max_diff*params.surface_tension, grid->dim, grid->ld.Scale(), 3.));
    }
  }
  ctrl.euler_dt.min(params.timestep_factor * my::min<double>(max_dt_vel, max_dt_src, max_dt_diff, max_dt_stf));
}


void Model::calcSlopeMixture(const State& state_, State& slope, NewSteppers::StepControl& ctrl)
{
  /* a new implementation with different coefficients, using simply
   * div KE grad (\phi \sigma(\Psi))
   */
  initStateVars(slope, true);
  State &state = const_cast<State&>(state_); // for the external scope the state will look unchanged

  #pragma omp parallel
  {
    ThreadLocal terms(grid->ld, grid->dim, params); // every thread has its own terms variable, reinitializing does NOT clear the max_xxx variables.

    BOOST_FOREACH(const BBox3 bbox, mtboxes->getCurrentThreadRange())
    {
      Array3df o2field = obtain_oxygen_conc(bbox);
      terms.reinit(bbox, state, o2field, other_volume, OutState::SOURCES);
      terms.AddSourceContribution(slope.cells, slope.necro);
    }
    
    #pragma omp critical
    { // combines the data from all threads
      max_dt_src.min(EulerDt::ByEigenvalue(terms.max_source));
    }
  }
  ctrl.euler_dt.min(params.timestep_factor * my::min<double>(max_dt_vel, max_dt_src, max_dt_diff, max_dt_stf));
}


#if 0
void Model::calcSlopesIMEXMixture(const State &state, State &slopef, State &slopeg, StepControl &ctrl)
{
  throw std::runtime_error("implicit method for bulktissue not implemented because surface tension cannot be done");
  if (slopef.phi_field[ID_CELLS].empty())
    slopef.init(state.num_phi_fields, state.num_levelsets, grid->ld, grid->dim, 0);
  if (slopeg.phi_field[ID_CELLS].empty())
    slopeg.init(state.num_phi_fields, state.num_levelsets, grid->ld, grid->dim, 0);

  /*
   * let U = phi_cells * Sigma(Psi(phi_cells, phi_other))
   * first compute coefficients in
   *  div (phi_cells * K * dU/dphi_cells * grad phi_cells) + div (phi_cells * K * dU/dphi_other * grad phi_other)
   * then add this whole term to the slope of the implicit operator g
   */

  for (int i=0; i<state.num_phi_fields; ++i)
    CopyBorder(const_cast<State&>(state).phi_field[i][grid->ir.cells], grid->dim, 2);

  #pragma omp parallel
  {
    ThreadLocal terms(grid->ld, grid->dim, params); // every thread has its own terms variable, reinitializing does NOT clear the max_xxx variables.
    FaceVarArrays velocities;
    #pragma omp for nowait schedule(dynamic, 1)
    for (int i=0; i<grid->mtboxes.size(); ++i)
    {
      const BBox3 &bbox = grid->mtboxes[i];
      terms.reinit(bbox, state, other_volume, o2field, OutState::VELOCITIES | OutState::SOURCES);
      velocities.init(bbox, grid->dim);

      // contribution from ecm variations
      terms.AddVelocityContribution(1, 1., velocities);
      AddWENOAdvection(bbox, state.phi_field[ID_CELLS], velocities.data(), grid->dim, grid->ld.Scale(), slopef.phi_field[ID_CELLS], CONSERVATIVE, 3);

      // contribution from cell variations
      terms.AddDiffusionContribution(0, 1., slopeg.phi_field[ID_CELLS]);

      // the source term
      terms.AddSourceContribution(slopef.phi_field[ID_CELLS], slopef.phi_field[ID_THETA_NECRO]);

      if (params.use_necrotic_regions) // advection of necrotic regions
      {
        for (int d=0; d<grid->dim; ++d) velocities[d].fill(0);
        terms.AddVelocities(velocities);
        AddWENOAdvection(bbox, state.phi_field[ID_THETA_NECRO], velocities.data(), grid->dim, grid->ld.Scale(), slopef.phi_field[ID_THETA_NECRO], CONSERVATIVE, 3);
      }
    }
    #pragma omp critical
    { // combines the data from all threads
      //max_dt_vel = std::min(max_dt_vel, 0.8 * 0.5*ld.Scale()/(dim*terms.max_vel+1.e-13));
      max_dt_diff = std::min(max_dt_diff, 0.8 * 0.5*my::sqr(grid->ld.Scale())/(grid->dim*terms.max_diff+1.e-13));
      max_dt_src = std::min(max_dt_src, 1./(1.e-13+terms.max_source));
      fail_flag |= terms.fail_flag;
    }
  }

  ctrl.euler_dt = params.timestep_factor * my::min<double>(max_dt_vel, max_dt_src, max_dt_vel);
}


void Model::invertImplicitOperatorMixture(State &lhs, const State &rhs, double identity_coeff, double operator_coeff, StepControl &ctrl, State &extrapolation)
{
  throw std::runtime_error("implicit method for bulktissue not implemented because surface tension cannot be done");
  /*
   * the implicit equation is now given by
   *  identity_coeff * I + operator_coeff * [ D_cells lhs.ncells + D_other other_volume ] = rhs
   * if existent, the contribution by other_volume is brought to the rhs and the elliptic equation for lhs.ncells is solved
   */

  {
    for (int i=0; i<lhs.num_phi_fields; ++i)
      CopyBorder(const_cast<State&>(extrapolation).phi_field[i][grid->ir.cells], grid->dim, 1);

    FiniteVolumeMatrixBuilder mb;
    mb.Init7Point(grid->ld, grid->dim);
    
    #pragma omp parallel
    {
      ThreadLocal terms(grid->ld, grid->dim, params); // every thread has its own terms variable, reinitializing does NOT clear the max_xxx variables.
      #pragma omp for nowait schedule(dynamic, 1)
      for (int i=0; i<grid->mtboxes.size(); ++i)
      {
        const BBox3 &bbox = grid->mtboxes[i];
        terms.reinit(bbox, extrapolation, other_volume, o2field, OutState::VELOCITIES);
        terms.AddDiffusionContribution(operator_coeff, mb);

        FOR_BBOX3(p, bbox)
          mb.AddLocally(p, identity_coeff, rhs.phi_field[ID_CELLS](p));
      }
    }
    extrapolation.clear();

    Epetra_Vector lhsv(mb.rhs->Map());

    boost::property_tree::ptree pt = pt_params.get_child("elliptic_solver_params");
    SolveEllipticEquation(*mb.m, *mb.rhs, lhsv, pt);

    FOR_BBOX3(p, grid->ir.cells)
    {
      lhs.phi_field[ID_CELLS](p) = lhsv[grid->ld.LatticeToSite(p)];
    }
  }

  lhs.levelset[0].addScaled(1./identity_coeff, rhs.levelset[0], 0.); // identity
  if (params.use_necrotic_regions)
    lhs.phi_field[ID_THETA_NECRO].addScaled(1./identity_coeff, rhs.phi_field[ID_THETA_NECRO], 0.);

  for (int i=0; i<lhs.num_phi_fields; ++i)
    CopyBorder(lhs.phi_field[i][grid->ir.cells], grid->dim, 3);
  CopyBorder(lhs.levelset[0].phi[grid->ir.cells], grid->dim, 3);
  
  //cout << format("dt = %f, by vel: %f, by diff = %f, by src = %f") % ctrl.euler_dt % max_dt_vel % max_dt_diff % max_dt_src << endl;
}
#endif

/// google split-step method
bool Model::doSplitStep(int num, State& state, NewSteppers::StepControl &ctrl)
{
  state.checksum += 1; // because the state will change
  if (num == 0)
  {
    my::LogScope logscope(my::log(), "levelset:");
    return stepper_convection.doStep(state, ctrl);
  }
  else if (num == 1)
  {
    my::LogScope logscope(my::log(), "mixture:");
    return stepper_mixture.doStep(state, ctrl);
  }
  else throw std::runtime_error("BulkTissueModel: unkown splitstep id");
}


void Model::doPreStep(State& state, NewSteppers::StepControl& ctrl)
{
  //printf("fuckreinit\n");
  if (last_levelset_reinit < ctrl.t - max_dt_vel)
  {
    ptree pt;
    reinit_levelset(state.ls, grid->ld, grid->dim, 5);
    last_levelset_reinit = ctrl.t;
    cout << format("levelset reinit, min dist: %f, requested: %f") % state.ls.max_val % (5.*grid->ld.Scale()) << endl;
  }
  max_dt_vel = max_dt_src = max_dt_diff = max_dt_stf = my::MakeSmallestMax<double>();
}


/*
  analytic formula for propagation speed:
KE = 3600
delta = 0.1/24
gamma = 1./24
psi_a = 0.8
psi_n = 0.6
import math
v = math.sqrt(2.*KE*delta*(1.-delta/gamma)*(psia-psin))
see preziosi et al., mathematical modelling of the loss of tissue compression responsiveness, p 220
// results in
2.32379000772
and for
>>> psin = 0.7
1.64316767252
*/

/** @brief Writes stuff to file, try to document the parameters here
 * 
 * @param conc total volume fraction of cells
 * @param ptc phase function tumor,
 * 		1 inside tumor
 * 		0 outside tumor
 * 		in between at boundndary
 * @param press tissue pressure
 * @param source ???
 * @param obstacle ???
 * @param stf surface tension force
 * @param vel 3D velocity of cells
 * @param ls level set function
 * @param necro fraction of necrotic cells?
 */
void Model::writeH5(H5::Group g, State &state, double t, H5::Group &ld_group) const
{
  const OutState &ost = getOutState(state, OutState::ALL);
  writeAttrToH5(g, string("TYPE"), string("BulkTissueFormat1"));
  //g.attrs().set("TYPE", "BulkTissueFormat1");
  const BBox3 cellrange = grid->ld.Box();
  WriteScalarField(g, "conc", ost.phi_cells[cellrange], grid->ld, ld_group);
  WriteScalarField(g, "ptc", ost.theta[cellrange], grid->ld, ld_group);
  WriteScalarField(g, "press", ost.pressure[cellrange], grid->ld, ld_group);
  WriteScalarField(g, "sources", ost.sources[cellrange], grid->ld, ld_group);
  if (pt_params.get<bool>("write_face_velocities"))
  {
    for (int i=0; i<grid->dim; ++i)
    {
      H5::DataSet ds = WriteArray3D(g, str(format("vel_%i") % i), ost.velocities[i][grid->ir.faces[i]]);
    }
  }
  if (pt_params.get<bool>("write_surface_tension_force", false))
  {
    WriteAveragedFaceVariableField(g, "stf", grid->dim, ost.surface_tension_force.data(), grid->ld, ld_group);
  }
  if (pt_params.get<bool>("write_levelset_function"))
    WriteScalarField(g, "ls", state.ls.phi[grid->ir.cells], grid->ld, ld_group);
  
  WriteScalarField(g, "obstacle", other_volume[grid->ir.cells], grid->ld, ld_group);
  WriteAveragedFaceVariableField(g, "vel", grid->dim, ost.velocities.data(), grid->ld, ld_group);
  
  if (params.use_necrotic_regions)
  {
    WriteScalarField(g, "necro", ost.phi_necro[grid->ir.cells], grid->ld, ld_group);
  }
  writeAttrToH5(g, "max_dt_diff", (double)max_dt_diff );
  writeAttrToH5(g, "max_dt_src", (double)max_dt_src );
  writeAttrToH5(g, "max_dt_vel", (double)max_dt_vel );
  writeAttrToH5(g, "max_dt_stf", (double)max_dt_stf );
}

const OutState& Model::getOutState(const State &state_, int flags) const
{
  State &state = const_cast<State&>(state_); // for the external scope the state will look unchanged
  
  int state_id = state.checksum;
  if (outstate.state_id == state_id) {
    if ((outstate.flags & flags) == flags) return outstate;
  }

  auto &ls = state.ls;
  auto cells = state.cells;
  auto necro = state.necro;
  outstate = OutState();

  //cout << format("recomputing output with flags %i") % flags << endl;

  if (flags & OutState::PHI_CELLS)
    outstate.phi_cells.initFromBox(grid->ir.cells);
  if (flags & OutState::THETA) {
    outstate.theta.initFromBox(grid->ir.cells);
    if (params.use_necrotic_regions)
      outstate.phi_necro.initFromBox(grid->ir.cells);
  }
  if (flags & OutState::PRESSURE)
    outstate.pressure.initFromBox(grid->ir.cells);
  if (flags & OutState::SOURCES)
    outstate.sources.initFromBox(grid->ir.cells);
  if (flags & OutState::VELOCITIES)
    outstate.velocities.init(grid->ir, grid->Dim());
  if (flags & OutState::SURFACE_TENSION_FORCE)
    outstate.surface_tension_force.init(grid->ir, grid->Dim());
  
  #pragma omp parallel
  {
    BOOST_FOREACH(const BBox3 bbox, mtboxes->getCurrentThreadRange())
    {
      CopyBorder(ls.phi, bbox, grid->Box(), grid->dim, 3);
      CopyBorder(cells, bbox, grid->Box(), grid->dim, 3);
      CopyBorder(necro, bbox, grid->Box(), grid->dim, 3);
    }

    #pragma omp barrier
    
    ThreadLocal terms(grid->ld, grid->dim, params); // every thread has its own terms variable
    Array3d<float> tmp;
    BOOST_FOREACH(const BBox3 bbox, mtboxes->getCurrentThreadRange())
    {
      Array3df o2field = obtain_oxygen_conc(bbox);
      myAssert(!o2field.empty());
      terms.reinit(bbox, state, o2field, other_volume, flags);
      // copy the stuff from the temporary storage to the global array
      if (flags & OutState::PHI_CELLS)
        outstate.phi_cells[bbox].fill(terms.phi_cells[bbox]);
      if (flags & OutState::THETA)
      {
        outstate.theta[bbox].fill(terms.theta[bbox]);
        if (params.use_necrotic_regions)
          outstate.phi_necro[bbox].fill(terms.phi_necro[bbox]);
      }
      if (flags & OutState::PRESSURE)
        outstate.pressure[bbox].fill(terms.pressure[bbox]);
      if (flags & OutState::SOURCES)
      {
        tmp.initFromBox(bbox);
        terms.AddSourceContribution(outstate.sources, tmp);
        outstate.sources[bbox] -= tmp;
      }
      if (flags & OutState::VELOCITIES)
      {
        terms.AddVelocities(outstate.velocities, false);
      }
      if (flags & OutState::SURFACE_TENSION_FORCE)
      {
        terms.SetSurfaceTensionForce(outstate.surface_tension_force);
      }
    }
  }
      
  outstate.state_id = state_id;
  outstate.flags = flags;
  return outstate;
}


/* returns normal, tumor and necrotic volume fractions
 */
boost::tuple<Array3df,Array3df, Array3df> Model::getTissuePhases(const BBox3 &bbox, const State &state) const
{
  const auto cells = state.cells;
  const auto necro = state.necro;
  const auto& ls   = state.ls;
  Array3df ph_tum(bbox, Cons::DONT), ph_normal(bbox, Cons::DONT), ph_necro(bbox, Cons::DONT);
  FOR_BBOX3(p, bbox)
  {
    float theta = Levelset::ToHeaviside(ls.phi(p), params.interface_width);
    ph_tum(p) = theta * (cells(p) - necro(p));
    ph_normal(p) = (1. - theta) * (cells(p) - necro(p));
    ph_necro(p) = necro(p);
  }
  return boost::make_tuple(
    ph_normal, ph_tum, ph_necro);
}


Array3df Model::getObstructionPhase(const BBox3 &bbox, const State &state) const
{
  return other_volume;
}


}




#if 0 // Trilinos::NOX test, does not work at all because of a bug in trilinos and also the system is not symmetric so NonlinearCG is bound to fail

#include <NOX.H>
#include <NOX_Epetra.H>
#include <Epetra_SerialComm.h>

#include <Epetra_Map.h>
#include <Epetra_Vector.h>

class SimpleProblemInterface : public NOX::Epetra::Interface::Required
{
public:
  Model &model;
  State &state;
  const State &rhs;
  double identity_coefficient, operator_coefficient;
  SimpleProblemInterface(Model &model_, State &state_, const State &rhs_,
                         double identity_coefficient_, double operator_coefficient_)
  : model(model_), state(state_), rhs(rhs_), identity_coefficient(identity_coefficient_), operator_coefficient(operator_coefficient_) {}

  void stateToVector(Epetra_Vector &x, const State &s)
  {
    int k = 0;
    FOR_BBOX3(p, model.ld.Box())
    {
      x[k++] = s.ncells(p);
    }
  }

  void vectorToState(State &s, const Epetra_Vector &x)
  {
    int k=0;
    FOR_BBOX3(p, model.ld.Box())
    {
      s.ncells(p) = x[k++];
    }
  }

  bool computeF(const Epetra_Vector & x, Epetra_Vector & f,
                NOX::Epetra::Interface::Required::FillType F )
  {
    State x_state; x_state.initCloned(state);

    vectorToState(x_state, x);

    model.evaluateImplicitOperatorMixture(x_state, rhs, identity_coefficient, operator_coefficient);

    stateToVector(f, x_state);

    return true;
  };
};


void Model::evaluateImplicitOperatorMixture(State &lhs, const State &rhs, double identity_coeff, double operator_coeff)
{
  Array3d<float> slope_g_ncells(ld.Box());

  CopyBorder(lhs.ncells[ir.cells], dim, 2);

  #pragma omp parallel
  {
    TempDataIMEX terms(ld, dim, params); // every thread has its own terms variable, reinitializing does NOT clear the max_xxx variables.
    FaceVarArrays velocities;
    #pragma omp for nowait schedule(dynamic, 1)
    for (int i=0; i<mtboxes.size(); ++i)
    {
      const BBox3 &bbox = mtboxes[i];
      terms.reinit(bbox, lhs, other_volume, o2field, OutState::VELOCITIES | OutState::SOURCES);
      InitFaceVar(bbox, dim, velocities);
      terms.AddVelocities(velocities);
      AddWENOAdvection(bbox, lhs.ncells, velocities.data(), dim, ld.Scale(), slope_g_ncells, CONSERVATIVE, 5);

      FOR_BBOX3(p, bbox)
      {
        //slope_g_ncells(p) = lhs.ncells(p);
        slope_g_ncells(p) = slope_g_ncells(p)*operator_coeff + identity_coeff*lhs.ncells(p) - rhs.ncells(p);
      }
    }
  }

  lhs.ncells[ir.cells].fill(slope_g_ncells[ir.cells]);
}


void Model::invertImplicitOperatorMixture(State &lhs, const State &rhs, double identity_coeff, double operator_coeff, StepControl &ctrl, State &extrapolation)
{
  /*
   * the implicit equation is now given by
   *  identity_coeff * I + operator_coeff * [ D_cells lhs.ncells + D_other other_volume ] = rhs
   * if existent, the contribution by other_volume is brought to the rhs and the elliptic equation for lhs.ncells is solved
   */
  {
  Teuchos::RCP<SimpleProblemInterface> problem(
    new SimpleProblemInterface(*this, lhs, rhs, identity_coeff, operator_coeff));

  Epetra_SerialComm epetra_comm;
  Epetra_Map epetra_map(product_of_components(::Size(ld.Box())), 0, epetra_comm);

  Epetra_Vector x_init(epetra_map);
  problem->stateToVector(x_init, lhs);

  Teuchos::RCP<Teuchos::ParameterList> nlParamsPtr =
    Teuchos::rcp(new Teuchos::ParameterList);
  Teuchos::ParameterList& nlParams = *(nlParamsPtr.get());
  nlParams.set("Nonlinear Solver", "Line Search Based");
  Teuchos::ParameterList& printParams = nlParams.sublist("Printing");
  printParams.set("MyPID", epetra_comm.MyPID());
  printParams.set("Output Precision", 3);
  printParams.set("Output Processor", 0);
  printParams.set("Output Information",
            NOX::Utils::OuterIteration +
            NOX::Utils::OuterIterationStatusTest +
            NOX::Utils::InnerIteration +
            NOX::Utils::Parameters +
            NOX::Utils::Details +
            NOX::Utils::Warning);
  Teuchos::ParameterList& searchParams = nlParams.sublist("Line Search");
  searchParams.set("Method", "Full Step");
  Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
  dirParams.set("Method", "NonlinearCG");

  NOX::Epetra::Vector noxInitGuess(x_init, NOX::DeepCopy);

  Teuchos::RCP<NOX::Epetra::Interface::Required> iReq = problem;

  Teuchos::RCP<NOX::Epetra::Group> grpPtr =
    Teuchos::rcp(new NOX::Epetra::Group(printParams,
                    iReq,
                    noxInitGuess));

  // Set up the status tests
  Teuchos::RCP<NOX::StatusTest::NormF> testNormF =
    Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-4));
  Teuchos::RCP<NOX::StatusTest::MaxIters> testMaxIters =
    Teuchos::rcp(new NOX::StatusTest::MaxIters(20));
  // this will be the convergence test to be used
  Teuchos::RCP<NOX::StatusTest::Combo> combo =
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR,
                        testNormF, testMaxIters));

  Teuchos::RCP<NOX::Solver::Generic> solver =
    NOX::Solver::buildSolver(grpPtr, combo, nlParamsPtr);

  // Solve the nonlinesar system
  NOX::StatusTest::StatusType status = solver->solve();

  if( NOX::StatusTest::Converged  == status )
    cout << "\n" << "-- NOX solver converged --" << "\n";
  else
    cout << "\n" << "-- NOX solver did not converge --" << "\n";

  // Print the answer
  cout << "\n" << "-- Parameter List From Solver --" << "\n";
  solver->getList().print(cout);

  // Get the Epetra_Vector with the final solution from the solver
  const NOX::Epetra::Group & finalGroup =
    dynamic_cast<const NOX::Epetra::Group&>(solver->getSolutionGroup());
  const Epetra_Vector & finalSolution =
      (dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getX())).getEpetraVector();

  problem->vectorToState(lhs, finalSolution);
  }
#if 0
  {
    CopyBorder(const_cast<State&>(extrapolation).ncells[ir.cells], dim, 1);

    FiniteVolumeMatrixBuilder mb;
    mb.Init7Point(ld, dim);

    #pragma omp parallel
    {
      TempDataIMEX terms(ld, dim, params); // every thread has its own terms variable, reinitializing does NOT clear the max_xxx variables.
      #pragma omp for nowait schedule(dynamic, 1)
      for (int i=0; i<mtboxes.size(); ++i)
      {
        const BBox3 &bbox = mtboxes[i];
        terms.reinit(bbox, extrapolation, other_volume, o2field, OutState::VELOCITIES);
        terms.AddDiffusionContribution(operator_coeff, mb);
        //Array3d<float> tmp(terms.bbox);
        //terms.AddDiffusionContribution(1, 1., tmp);
        FOR_BBOX3(p, bbox)
          mb.AddLocally(p, identity_coeff, rhs.ncells(p));
        //  mb.AddLocally(p, identity_coeff, rhs.ncells(p) - operator_coeff * tmp(p));

      }
    }
    extrapolation.clear();

    Epetra_Vector lhsv(*mb.rhs);

    boost::property_tree::ptree pt;
    pt.put("max_iter", 500);
    pt.put("max_resid", 1.e-7);
    pt.put("output", 1);
    boost::property_tree::update(pt, pt_params.get_child("elliptic_solver_params", ptree()));
    SolveEllipticEquation(*mb.m, *mb.rhs, lhsv, pt);

    FOR_BBOX3(p, ir.cells)
    {
      lhs.ncells(p) = lhsv[ld.LatticeToSite(p)];
    }
  }
#endif

  lhs.levelset().addScaled(1./identity_coeff, rhs.levelset(), 0.); // identity

  CopyBorder(lhs.ncells[ir.cells], dim, 1);
  CopyBorder(lhs.levelset().phi[ir.cells], dim, 1);

  //cout << format("dt = %f, by vel: %f, by diff = %f, by src = %f") % ctrl.euler_dt % max_dt_vel % max_dt_diff % max_dt_src << endl;
}


#endif

#endif


