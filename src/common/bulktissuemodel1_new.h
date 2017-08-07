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
#pragma once // include this file only once per compilation unit (see https://en.wikipedia.org/wiki/Pragma_once)

#include "common.h"
#include "hdfio.h"
#include "shared-objects.h"

#include "mwlib/time_stepper_new.h"
#include "mwlib/any_tree.h"

#include "levelset.h"
#include "continuum-flow.h"
#include "continuum-stepper-states.h"

#include <boost/array.hpp>

using boost::property_tree::ptree;
//using NewSteppers::StepControl;

/**
@brief Following the papers from Preziosi et al.. 
See also the Interstitial fluid flow paper 2013. 
The model represents the fractional volume of tissue 
cells within unit volumes (phi_cells). This cell polulation
is divided into tumor cells and normal cells by 
a sharp interface described by the zero-level of a levelset function. 
In addition a population of necrotic cells (dead cell debris) is explicitely represented
(phi_necro). The volume fraction of viable cells is therefore phi_cells-phi_necro.
The entire cell population moves with a common velocity field.

Note: The model implies essentially a diffusive motion of biological cells.
Here is why: Advection terms goes like  div * (v * u), 
where v is the velocity and u the cell density distribution.
But v is essentially proportional to grad u, so in most cases v = const * grad u,
and thus du/dt = ... + div * ((const u) * grad u), 
with the variable diffusion coefficient (u * const). This must
be taken into account for the selection of time steps and 
algorithms for the numerical solution.
*/
namespace NewBulkTissueModel
{

class Model;


struct TempData;

/**
@brief Interfaces with the world outside of the continuum model.
  units:
    time - h
    space - micron
*/
struct OutState
{
  enum Flags {
    PHI_CELLS = 1<<0, 
    THETA = 1<<1, 
    PRESSURE = 1<<2, 
    SOURCES = 1<<3, 
    VELOCITIES = 1<<4, 
    SURFACE_TENSION_FORCE = 1<<5, // DEBUG_INFO = 1<<5,
    ALL = ~0
  };
  Array3d<float> phi_cells, /// total volume fraction of all cell populations
		 theta,     /// level set function
		 pressure,  /// solid pressure, gradient of it generates velocity of cell movement. Function of phi_cells and theta (?). 
		 sources,   /// source strength of volume fraction distribution (rate of proliferation in volume added per time and unit volume)
		 phi_necro; 
  // See Interstitial fluid flow paper 2013, and papers from BulkTissue!
  FaceVarArrays velocities, /// velocity of cell movement, defined on all faces of the cells of the numerical grid (cubic lattice). See "staggered grid".
	        surface_tension_force; /// due to cell-cell adhesion. Derived from approximation of local curvature of the interface defined by levelset function.
  int state_id,
      flags;
  OutState() : state_id(-1), flags(0) {}
};

struct Params
{
  Params();
  void assign(const ptree &pt);
  ptree as_ptree() const;
  void init(const LatticeDataQuad3d &ld);
  
  float  time_prol, /// mean division times (?)
         time_death, /// mean survival time before apoptosis (?)
         time_prol_tumor,
         time_death_tumor,
         time_necrosis, 
         tumor_diameter,  /// stop condition (?)
         cell_compression_modulus, // proportionality between pressure and cell volume fraction
         sigma_move, /// cell population only moves if force is greater than this threshold (?), normally set to 0
         ncells_norm, /// ?
         cell_mobility, /// proportionality constant between velocity and force
         cell_mobility_tumor,
         cell_mobility_ecmstar, /// Parameterization of the cell-velocity-vs-force relation
         cell_mobility_ecmstarstar,
         ncells_tumor,
         ncells_ecm,
         ncells_sigma,
         timestep_factor, /// additional scaling of time step (use < 1.)
         interface_width, /// compute distance map (= levelset function) within interface_with units of the zero level.
         o2_necro_threshold, /// threshold below which cells die
         o2_prol_threshold,  /// threshold above which cells can proliferate
         surface_tension;   /// coefficient, proportionality between curvature of the interface and force
  bool wrong_model, use_necrotic_regions;
  int sigma_model_version, source_model_version;
  float psi_norm, psi_norm_tumor; /// ??
  float ecm_noise_std; /// std deviation
  float ecm_noise_density; /// percentage of sites which get a noise offset, 0 = no noise, between 0 and 1 -> salt and pepper like
  int random_seed;
};

/**
@brief compute the fractional volume of cells generated per time per unit volume
*/
void calcSource(float theta, float phi_cells, float phi_other, float phi_necro, float o2, const Params& params, float *outrates);

// typedef boost::tuple<Array3dOps<float>, Ops2, Ops3> OpsTuple;
// 
// struct OpsStruct
// {
//   void do()
//   {
//     ops1.do(state.a1);
//     ops2.do(state.a2);
//   }
// };

/**
@brief Internal state representation of the model
*/
struct State
{
  Levelset ls;            /// describes interface between tumor cells and normal cells (see papers on Level Set Method)
  Array3df cells, necro;  /// cells = volume fraction of all cell types (normal, tumor, necrotic), necro = volume fraction of necrotic cells
  int checksum;           /// increased when the state changes
  State() : checksum(0) {}
};


#if 0
class PolyStateBase
{
};

class PolyStateList : public PolyStateBase
{
  std::vector<PolyStateBase*> children;
  
  PolyStateBase& get(const std::string &key) { return children[key]; }
};

template<class T>
class PolyStateArray3d
{
};


class PolyStateBaseOps
{
  void initCloned(PolyStateBase* dst, PolyStateBase *src); // copies data
  void copyShared(PolyStateBase* dst, PolyStateBase *src); // init as a shared view
};

class PolyStateListOps
{
};

class PolyStateArray3dOps
{
};

void initCloned(PolyStateBase* dst, PolyStateBase *src)
{
  // full_ops->copyShared(dst, src) // copy shared pointers from src to dst (dst becomes a shallow copy of src)
  // alternatively just make shallow copies of stuff that is needed but not changed
  // full_ops[mystuff]->initCloned(dst[mystuff], src[mystuff]) // the parts which are changed by the stepper need deep copies
};

void addScaled(double fa, PolyStateBase* dst, double fb, PolyStateBase* src)
{
  //full_ops[mystuff]->addScaled(fa, dst[mystuff], fb, src[mystuff])
}

void constructOps()
{
}
#endif


/**
@brief Used to obtain spatial distributions from the outside world, e.g. of oxygen.
*/
typedef boost::function<Array3df (const BBox3 &)> DataRequestCallback;

/**
@brief does all the work.
*/
class Model : public boost::noncopyable
{
public:
  typedef NewBulkTissueModel::State State;
public:
  void init(const ContinuumGrid &grid_, const DomainDecomposition &mtboxes_,
            Model::State& state,
            DataRequestCallback obtain_vessel_volume_,
            DataRequestCallback obtain_oxygen_conc_,
            const ptree& pt_params_);
  void notify_auxdata_change();
  
  const OutState& getOutState(const State &state, int flags = OutState::ALL) const;
  boost::tuple<Array3df,Array3df, Array3df> getTissuePhases(const BBox3 &bbox, const State &state) const;
  Array3df getObstructionPhase(const BBox3 &bbox, const State &state) const;
  void cloneState(State &dst, const State &src);
  
  void initStateVars(State& state, bool clear) const;
  void calcSlopeConvection(const State &state, State &slope, NewSteppers::StepControl &ctrl);
  void calcSlopeMixture(const State &state, State &slope, NewSteppers::StepControl &ctrl);
#if 0
  void calcSlopesIMEXMixture(const State &state, State &slopef, State &slopeg, StepControl &ctrl);
  void invertImplicitOperatorMixture(State &lhs, const State &rhs, double identity_coeff, double operator_coeff, StepControl &ctrl, State &extrapolation);
  void evaluateImplicitOperatorMixture(State &lhs, const State &rhs, double identity_coeff, double operator_coeff);
#endif
  void writeH5(h5cpp::Group g, State &state, double t, h5cpp::Group &ld_group) const; 
  //void appendToImages(State &state, std::vector<Image> &images) const;

  void doPreStep(State &state, NewSteppers::StepControl &ctrl);
  bool doSplitStep(int num, State& state, NewSteppers::StepControl& ctrl);
  
  // read only plx
  Params params;
  ptree pt_params;
  /**
  @brief Current maximal admissible time steps due to individual terms in the model equations: 
    
    Advection velocity, diffusive behavior implied by the dependence of the velocity on the cell volume fraction, Source Rate, Surface tension force.
    The smallest of these is the time step actually taken.
  */
  my::Smallest<double> max_dt_vel, max_dt_diff, max_dt_src, max_dt_stf;

private:
  const ContinuumGrid *grid;  /// information about the computational lattice
  const DomainDecomposition *mtboxes;  /// describes the division of the computational lattice into smaller boxes to facilitate multi threading (one box per thread).
  
  Array3dOps<float> ops_array3d;  /// encapsulate mathematical operations needed by time steppers (see NewSteppers)
  LevelsetOps ops_levelset;

  /**
  @brief The split-step method is applied to integrate the advection terms separately from the rest of the model equations. These are the individual integrators
  */
  NewSteppers::ImprovedExplicitEuler<NewSteppers::Callback<State>, NewSteppers::OpsCallback<State> > stepper_mixture;
  NewSteppers::Rk3<NewSteppers::Callback<State>, NewSteppers::OpsCallback<State> > stepper_convection;
  DataRequestCallback obtain_oxygen_conc, obtain_vessel_volume;
  Array3df other_volume; /// of blood vessels (not used; experimental.
  
  mutable OutState outstate;

  double last_levelset_reinit; /// last time the redistancing operation was performed on the levelset function
  double last_euler_dt; /// last time step taken
};


/**
@brief A public wrapper around calcSigma. Calculates solid pressure.
  @param phi_cells: cell volume fraction of all cell components
  @param phi_other: volume fraction of other solid components. NewBulkTissueModel puts the ECM here. Optionally, also the vessel volume fraction.
  @param params: Parameters
*/
float CalcSolidPressure(float phi_cells, float phi_other, const Params &params);


} // namespace NewBulkTissueModel
