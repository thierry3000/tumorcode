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
#include "calc-ifdrug.h"
#include "time_stepper_utils.h"
#include "mwlib/log.h"

#ifdef USE_IFDRUGSIM

namespace IfDrug
{



Params::Params()
{
  inject_mode             = DF_INJECT_MODE_EXP;
  max_uptake              = 0.0;
  chalf_uptake            = 0.5;
  capillary_permeability_normal  = 0.00017;
  capillary_permeability_tumor  = 0.017;
  inject_t                = 1.0f; // hours
  kdiff                   = 16;
  inject_max              = 1.0;
  linear_uptake_coeff     =  2./(3600.0); // seems reaonable, based on diffusion distance and mean residence time
  uptake_mode             = DF_UPTAKE_LINEAR;
  stepper                 = "vsimexbdf2";
  stepper_compartments    = "vsimexbdf2";
  comprates_k12           = 0.;
  comprates_k21           = 0.;
  convective_transport_enabled = true;
}


ptree Params::as_ptree() const
{
  boost::property_tree::ptree pt;
  #define DOPT(name) pt.put(#name, name)
  #define PTSWITCH(name, caseid) case caseid: pt.put<string>(#name, #caseid); break;
  switch (inject_mode)
  {
    PTSWITCH(inject_mode, DF_INJECT_MODE_EXP)
    PTSWITCH(inject_mode, DF_INJECT_MODE_JUMP)
    PTSWITCH(inject_mode, DF_INJECT_MODE_JUMP_UP)
  }
  DOPT(max_uptake);
  DOPT(chalf_uptake);
  DOPT(capillary_permeability_normal);
  DOPT(capillary_permeability_tumor);
  DOPT(inject_t);
  DOPT(kdiff);
  DOPT(inject_max);
  DOPT(linear_uptake_coeff);
  switch (uptake_mode)
  {
    PTSWITCH(uptake_mode, DF_UPTAKE_LINEAR)
    PTSWITCH(uptake_mode, DF_UPTAKE_MICHAELIS_MENTEN)
  }
  DOPT(stepper);
  DOPT(stepper_compartments);
  DOPT(comprates_k12);
  DOPT(comprates_k21);
  DOPT(convective_transport_enabled);
  #undef DOPT
  #undef PTSWITCH
  return pt;
}


void Params::assign(const ptree& pt)
{
  #define DOPT(name) boost::property_tree::get(name, #name, pt)
  #define PTSWITCH(name, caseid) (pt.get<string>(#name) == #caseid) name = caseid;
  if PTSWITCH(inject_mode, DF_INJECT_MODE_EXP)
  else if PTSWITCH(inject_mode, DF_INJECT_MODE_JUMP)
  else if PTSWITCH(inject_mode, DF_INJECT_MODE_JUMP_UP)
  DOPT(max_uptake);
  DOPT(chalf_uptake);
  DOPT(capillary_permeability_normal);
  DOPT(capillary_permeability_tumor);
  DOPT(inject_t);
  DOPT(kdiff);
  DOPT(inject_max);
  DOPT(linear_uptake_coeff);
  if PTSWITCH(uptake_mode, DF_UPTAKE_LINEAR)
  else if PTSWITCH(uptake_mode, DF_UPTAKE_MICHAELIS_MENTEN)
  DOPT(stepper);
  DOPT(comprates_k12);
  DOPT(comprates_k21);
  DOPT(stepper_compartments);
  DOPT(convective_transport_enabled);
  #undef DOPT
  #undef PTSWITCH
}


void Params::WriteH5(H5::Group g) const
{
//#define SETP(id) g.attrs().set(#id,id)
//#define SETP(id) writeAttrToH5<string>(g, #id, id)
//#define SETP2(id,val) g.attrs().set(#id,val)
  switch(inject_mode)
  {
  case Params::DF_INJECT_MODE_EXP:
    writeAttrToH5(g, string("DF_INJECT_MODE"),string("DF_INJECT_MODE_EXP")); 
    //SETP2(DF_INJECT_MODE,"DF_INJECT_MODE_EXP");
    break;
  case Params::DF_INJECT_MODE_JUMP:
    //SETP2(DF_INJECT_MODE,"DF_INJECT_MODE_JUMP");
    writeAttrToH5(g, string("DF_INJECT_MODE"),string("DF_INJECT_MODE_JUMP"));
    break;
  }
  writeAttrToH5(g, string("inject_t"),inject_t);
  writeAttrToH5(g, string("kdiff"),kdiff);
  writeAttrToH5(g, string("comprates_k12"),comprates_k12);
  writeAttrToH5(g, string("comprates_k21"),comprates_k21);
  writeAttrToH5(g, string("capillary_permeability_normal"),capillary_permeability_normal);
  writeAttrToH5(g, string("capillary_permeability_tumor"),capillary_permeability_tumor);

  writeAttrToH5(g, string("stepper"),stepper);
  writeAttrToH5(g, string("stepper_compartments"),stepper_compartments);
  
  writeAttrToH5(g, string("convective_transport_enabled"),convective_transport_enabled);
 
  

}


static const PdeReal maturation_at_r5 = GetInitialThickness(5.0f);
static const PdeReal maturation_at_r20= GetInitialThickness(20.0f);

void VesselDrugExavasationModel::Init(const VesselList3d& vl_, const ContinuumGrid& grid_, const VesselsInBoxes& vesselsinboxes_, const Params& params_)
{
  VesselExavasationModel::Init(vl_, grid_, vesselsinboxes_);
  params = &params_;
}

/** @brief 
 * sets values for ceffa, ceffb --> concentration effective at lattice point a/b
 * currently vala, valb is not used
 * why did he implement it? -> copy past from FluidExavasationModel, there one needs this
 * calculates continously intermitant value for resistance between
 * tumor and normal resistance dependent on vessels maturation
 * 
 * for vala, valb see getSourceImplicit, inflow_conc
 */
bool VesselDrugExavasationModel::GetVesselInfo(const Vessel* v, double& coeffa, double& coeffb, double& vala, double& valb) const
{
  coeffa = coeffb = vala = valb = 0;
  if(!v->IsCirculated()) return false;
  double r = 0;
  if (params->capillary_permeability_tumor > 1.e-13 && params->capillary_permeability_normal > 1.e-13)
  {
    double r0 = 1./params->capillary_permeability_tumor;
    double r1 = 1./params->capillary_permeability_normal;
    r = std::max<double>(r0, r1 * v->maturation / maturation_at_r5);
    r = 1./r;
    
//     r = 1./params->capillary_permeability;
//     r = std::max(0.1, v->maturation / maturation_at_r5) * r;
//     r = 1./r;
  }
  else r = 0.;
  coeffa = coeffb = r;
  myAssert(coeffa == coeffb && coeffa >= 0.);
  return true;
}

/**
 * begin the calcualtion
 * \params ContinuumGrid, BoundingBoxed for threading,
 * \params IfFlowState, tissueCompositionCallback_
 * \params vesselModel_, const IfDrug::Params& params_
 */
void Calculator::Init(const ContinuumGrid& grid_, const DynArray< BBox3 >& mtboxes_, const IfFlowState& iffstate_, IfDrug::Calculator::Callback1 tissueCompositionCallback_, IfDrug::Calculator::VesselModel& vesselModel_, const IfDrug::Params& params_)
{
  grid = &grid_;
  mtboxes = &mtboxes_;
  iffstate = &iffstate_;
  params = &params_;
  time = 0;
  //cell_fraction_f = &cell_fraction_f_;
  //water_fraction_f = &water_fraction_f_;
  tissueCompositionCallback = tissueCompositionCallback_;
  
//   extracellular_model.init(*grid, make_ptree("stepper", params->stepper));
//   extracellular_model.diffusivity_f = boost::bind(&Calculator::getDiffusivity,this, _1, _2, _3, _4);
//   extracellular_model.source_implicit_f = boost::bind(&Calculator::getSourceImplicit, this, (int)FLAG_ALL, _1, _2, _3, _4, _5);
//   extracellular_model.velocity_f = boost::bind(&Calculator::getVelocity, this, _1, _2, _3, _4, _5, _6);

  stepper_interaction                                 = Steppers::createCallbackBased<State>(params->stepper_compartments);
  stepper_interaction->model().calcSlope              = boost::bind(&Calculator::calcSlopeEc, this, _1, _2, _3);
  stepper_interaction->model().calcSlopesIMEX         = boost::bind(&Calculator::calcSlopesEc, this, _1, _2, _3, _4);
  stepper_interaction->model().invertImplicitOperator = boost::bind(&Calculator::invertImplicitOperatorEc, this, _1, _2, _3, _4, _5, _6);
  euler_dt_interact = 0.; max_rate_interact = 0.;

  stepper_interstitial                                 = Steppers::createCallbackBased<State1>(params->stepper);
  stepper_interstitial->model().calcSlope              = boost::bind(&Calculator::calcSlopeIf, this, _1, _2, _3);
  stepper_interstitial->model().calcSlopesIMEX         = boost::bind(&Calculator::calcSlopesIf, this, _1, _2, _3, _4);
#if 1
  stepper_interstitial->model().invertImplicitOperator = boost::bind(&Calculator::invertImplicitOperatorIf, this, _1, _2, _3, _4, _5, _6);
#endif
  max_kdiff = 0.; max_src_expl = 0.; max_src_impl = 0.; max_vel = 0.; euler_dt_kdiff = 0.; euler_dt_srcexpl = 0.; euler_dt_srcimpl = 0.; euler_dt_vel = 0.;

  {
    vessel_source_coefficient_buffer.initFromBox(grid->ir.cells);
    Array3d<double> dummy(grid->ir.cells);
    //VesselDrugExavasationModel m; m.Init(vl, *grid, vesselsinboxes_, *params);

    #pragma omp parallel for schedule(dynamic, 1)
    for (int boxindex = 0; boxindex<mtboxes->size(); ++boxindex)
    {
      VirtualGridFunctions<double,2>::List arr = {{ vessel_source_coefficient_buffer, dummy }};
      vesselModel_.GetValues(boxindex, (*mtboxes)[boxindex], arr);
    }
  }
}


#if 0
void Calculator::InitVessels(const VesselList3d &vl, const VesselsInBoxes& vesselsinboxes_)
{
  vessel_source_coefficient_buffer.initFromBox(grid->ld.Box());
  Array3d<double> dummy(grid->ld.Box());
  
  VesselDrugExavasationModel m; m.Init(vl, *grid, vesselsinboxes_, *params);

  #pragma omp parallel for schedule(dynamic, 1)
  for (int boxindex = 0; boxindex<grid->mtboxes.size(); ++boxindex)
  {
    VirtualGridFunctions<double,2>::List arr = {{ vessel_source_coefficient_buffer, dummy }};
    m.GetValues(boxindex, grid->mtboxes[boxindex], arr);
  }
}
#endif


enum {
  MAX_KDIFF = 0,
  MAX_VEL,
  MAX_SRC_IMPL,
  MAX_SRC_EXPL,
  MAX_EXCHANGE,
  MAX_CNT
};
/**
 * this is what each thread needs to know
 */
struct ThreadLocal
{
  typedef Calculator::T T; ///float
  const Calculator *calculator;
  int boxindex;
  const DynArray<BBox3> *mtboxes;
  const ContinuumGrid *grid;
  ConstArray3d<T> conc;
  BBox3 bbox, bboxext;
  const Params *params;
  const IfFlowState *iffstate;
  Array3d<float> phi_cells, phi_water, phi_ecm; ///fraction of each at point in space
  Array3d<double> vessel_src_coeff; ///source due to vessels

  struct DiffusivityFunctor
  {
    Array3d<float> kdiff;
    double operator()(const Int3& p, const Int3& q, int axis, int side) const
    {
      return 0.5*(kdiff(p)+kdiff(q));
    }
  };
  ///Constructor for a single thread
  ThreadLocal(const Calculator* calculator_, ConstArray3d<T> conc_) : calculator(calculator_), conc(conc_), boxindex(-1)
  {
    params   = calculator->params;
    iffstate = calculator->iffstate;
    grid     = calculator->grid;
    mtboxes  = calculator->mtboxes;
  }

  void init(int boxindex_)
  {
    boxindex = boxindex_;
    bbox     = (*mtboxes)[boxindex];
    bboxext  = ExtendForDim(bbox, grid->dim, 1);
    phi_cells.initFromBox(bboxext);
    phi_water.initFromBox(bboxext);
    phi_ecm.initFromBox(bboxext);
    calculator->tissueCompositionCallback(boxindex, bboxext, phi_cells, phi_water, phi_ecm);
    vessel_src_coeff = calculator->vessel_source_coefficient_buffer; 
  }

  /**
    * Q_s1 equation 20
    * write to clinear and cconst array
    */
  void getSourceImplicit(int flags, Array3d<T> clinear, Array3d<T> cconst) const
  {
    T inflow_conc = calculator->GetInflowVesselConc();
    FOR_BBOX3(p, bbox)
    {
      //T c = conc(p);
      //myAssert(!(c<-0.1 || c>params->inject_max*1.1));
      T vessel_source = iffstate->vessel_inflow_field(p);
      T vessel_sink = iffstate->vessel_outflow_field(p);
      T lymph_sink = iffstate->lymph_outflow_field(p);
      //T ncells = phi_cells(p);
      float phi = phi_water(p) + phi_ecm(p);

      T loc_clinear=0., loc_cconst=0.;

      // superseeded by 2 component model
  #if 0
      if (flags & FLAG_UPTAKE)
      {
        // uptake
        float uptake;
        if (params->uptake_mode == Params::DF_UPTAKE_MICHAELIS_MENTEN) {
          uptake = (-params->max_uptake) / (c + params->chalf_uptake); // -lambda * c^n+1 / (c^n + c05)
        }
        else {
          uptake = (-params->linear_uptake_coeff);
        }
        loc_clinear += ncells * uptake;
      }
  #endif
      // contribution diffusion wall
      if (flags & Calculator::FLAG_SRC_DIFF_VESS)
      {
        // vessel wall diffusive sources
        float vc = -vessel_src_coeff(p); // vessel_src_coeff is (-q) from  q (c_v - c)
        loc_cconst += vc * inflow_conc * phi;///diffusivity of wall times inflow_conc times amount of liquid -->const background
        loc_clinear -= vc;
      }
      // contribution flow vesesls as source
      if (flags & Calculator::FLAG_SRC_FLOW_VESS)
        loc_cconst += vessel_source * inflow_conc;
      // contribution flow vesesls as sink
      if (flags & Calculator::FLAG_SINK_FLOW_VESS)
        loc_clinear += vessel_sink; // / (phi + 1.e-13);
      // contribution lympatics
      if (flags & Calculator::FLAG_SINK_FLOW_LYMPH)
        loc_clinear += lymph_sink; // / (phi + 1.e-13);

      clinear(p) = loc_clinear;
      cconst(p) = loc_cconst;
      myAssert(std::isfinite(loc_clinear) && loc_clinear <= 0.);
      myAssert(std::isfinite(loc_cconst) && loc_cconst >= 0.);
    }
  }//getSourceImplicit


  void calcSlopesIf(Calculator::State1 &slope_expl, Calculator::State1 &slope_impl, double *thread_max_val) const
  {
    { //sources (implicit)
      Array3d<T> src_val(bbox), src_clin(bbox);
      getSourceImplicit(Calculator::FLAG_ALL, src_clin, src_val);
      FOR_BBOX3(p, bbox)
      {
        slope_impl(p) += src_clin(p) * conc(p) + src_val(p);
        thread_max_val[MAX_SRC_IMPL] = std::max<double>(thread_max_val[MAX_SRC_IMPL], std::abs(src_clin(p)));
      }
    }

    { 
      /**
       * diffusion (explicit)
       * slope is know here, so we can calculate second derivative (diffusion)
       */
      Array3d<float> kdiff(bboxext);
      Array3d<T> buffer(bboxext);
      FOR_BBOX3(p, bboxext)
      {
	///amount of [volume]
	/// l, Interstistial fraction l is here called a eq. 20
        float a = (phi_ecm(p)+phi_water(p)+1.e-13); // this is not going to end well if it becomes close to 0
        buffer(p) = conc(p) / a;  // this makes the operator non-symmetric! ///amount of sustance
        kdiff(p) = params->kdiff * a;/// D_s in equation 20
      }
      ptree pt = make_ptree("geometric_coeff_mean", false);
      AddLaplacian<T>(bbox, buffer, kdiff, grid->dim, grid->ld.Scale(), slope_impl, pt);
      thread_max_val[MAX_KDIFF] = 1.5*params->kdiff; // hope params->kdiff works ...
    }
    /**
     * LHS of equation 20
     */
    if (params->convective_transport_enabled)
    {
      const LatticeIndexRanges ir = LatticeIndexRanges::FromCellRange(bbox);
      FaceVarArrays velocities(ir, grid->Dim());
      for (int d=0; d<grid->dim; ++d)
      {
        FOR_BBOX3(p, ir.faces[d])
        {
          const Int3 q = add_to_axis(((Eigen::Matrix<int,3,1>)p), d, -1);
          float v = iffstate->velocities[d](p);
          float alpha_l = phi_water(q);
          float alpha_m = phi_ecm(q);
          v *= alpha_l / (alpha_l+alpha_m + 1.e-13);
          velocities[d](p) = v;
          thread_max_val[MAX_VEL] = std::max<double>(std::abs(v), thread_max_val[MAX_VEL]);
        }
      }
      //AddWENOAdvection<T>(bbox, conc, velocities.data(), grid->dim, grid->ld.Scale(), slope_expl, CONSERVATIVE, 5);
      AddKTAdvection<T>(bbox, conc, velocities.data(), grid->dim, grid->ld.Scale(), slope_expl, CONSERVATIVE);
    }
  }

  void calculateSources(int flags, Array3d<T> res) const
  {
    Array3d<T> clinear(bbox), cconst(bbox);
    getSourceImplicit(flags, clinear, cconst);
    FOR_BBOX3(p, bbox)
      res(p) = clinear(p)*conc(p) + cconst(p);
  }

  void AddDiffusion(double operator_coeff, FiniteVolumeMatrixBuilder& mb)
  {
    Array3d<float> kdiff(bboxext);
    Array3d<float> buffer(bboxext);
    FOR_BBOX3(p, bboxext)
    {
      //l in paper
      float a = (phi_ecm(p)+phi_water(p)+1.e-13); // this is not going to end well if it becomes close to 0
      buffer(p) = 1. / a;  // this makes the operator non-symmetric!
      kdiff(p) = params->kdiff * a;
    }
      
    DiffusivityFunctor f_kdiff = { kdiff };

    mb.AddDiffusion2(bbox, grid->dim, f_kdiff, operator_coeff, false, buffer);
  }
};





void Calculator::calcSlopeIf(const State1& state, State1& slope, StepControl& ctrl)
{
  State1 tmp;
  calcSlopesIf(state, slope, tmp, ctrl);
  slope.addScaled(1., tmp);

  ctrl.euler_dt = my::min(ctrl.euler_dt, euler_dt_kdiff, euler_dt_srcimpl);
}


void Calculator::calcSlopesIf(const State1& state, State1& slope_expl, State1& slope_impl, StepControl& ctrl)
{
  slope_expl.initFromBox(grid->ir.cells);
  slope_impl.initFromBox(grid->ir.cells);

  const Array3d<T> conc = state;
  CopyBorder(const_cast<Array3d<T>&>(conc)[grid->ir.cells], grid->dim, 3);
  
  #pragma omp parallel
  {
    ThreadLocal local(this, state);
    Eigen::Matrix<double, MAX_CNT, 1> thread_max_val = Eigen::Matrix<double, MAX_CNT, 1>::Zero();
    
    #pragma omp for nowait schedule(dynamic, 1)
    for (int boxindex=0; boxindex<mtboxes->size(); ++boxindex)
    {
      local.init(boxindex);
      local.calcSlopesIf(slope_expl, slope_impl, thread_max_val.data());
    }

    #pragma omp critical
    {
      max_src_expl = std::max(max_src_expl, thread_max_val[MAX_SRC_EXPL]);
      max_src_impl = std::max(max_src_impl, thread_max_val[MAX_SRC_IMPL]);
      max_kdiff = std::max(max_kdiff, thread_max_val[MAX_KDIFF]);
      max_vel = std::max(max_vel, thread_max_val[MAX_VEL]);
    }
  }

  euler_dt_srcexpl = 0.5/(1.e-13+max_src_expl);
  euler_dt_vel  = 0.8 * 0.5*grid->ld.Scale()/(grid->dim*max_vel+1.e-13);
  euler_dt_kdiff = 0.8 * 0.5*my::sqr(grid->ld.Scale())/(grid->dim*max_kdiff+1.e-13);
  euler_dt_srcimpl = 0.5/(1.e-13+max_src_impl);

  ctrl.euler_dt = my::min(euler_dt_srcexpl, euler_dt_vel, 60.);//euler_dt_kdiff //euler_dt_srcimpl
  
}


#if 1
void Calculator::insertLinearSystemCoefficientsIf(int boxindex, const BBox3& bbox, const State1& rhs, const State1& extrapolation, double identity_coeff, double operator_coeff, FiniteVolumeMatrixBuilder& mb)
{
  ThreadLocal local(this, extrapolation);
  local.init(boxindex);
  
  local.AddDiffusion(operator_coeff, mb);

  Array3d<T> src_clin(bbox), src_val(bbox);
  local.getSourceImplicit(Calculator::FLAG_ALL, src_clin, src_val);
  
  FOR_BBOX3(p, bbox)
    mb.AddLocally(p, operator_coeff*src_clin(p) + identity_coeff, -operator_coeff*src_val(p) + rhs(p));
}


void Calculator::invertImplicitOperatorIf(State1& lhs, const State1& rhs, double identity_coeff, double operator_coeff, StepControl& ctrl, State1& extrapolation)
{
  // solve  (identity_coeff * I + operator_coeff * A) lhs = rhs
  ptree pt_params = make_ptree("preconditioner","multigrid")("output", 1)("max_iter", 100)("solver","bicgstab")("max_resid", 1.e-7);
  StationaryDiffusionSolve(grid->ld, *mtboxes, grid->dim,
                           lhs,
                           boost::bind(&Calculator::insertLinearSystemCoefficientsIf, this, _1, _2, rhs, extrapolation, identity_coeff, operator_coeff, _3),
                           pt_params);
}
#endif

#if 0
void Calculator::invertImplicitOperatorIf(State1& lhs, const State1& rhs, double identity_coeff, double operator_coeff, StepControl& ctrl, State1& extrapolation)
{
  // solve  (identity_coeff * I + operator_coeff * A) lhs = rhs
  
  #pragma omp parallel
  {
    ThreadLocal local(this, extrapolation);
    Array3d<T> src_linear, src_const;

    #pragma omp for schedule(dynamic, 1)
    for (int boxindex=0; boxindex<grid->mtboxes.size(); ++boxindex)
    {
      local.init(boxindex);
      src_linear.initFromBox(local.bbox);
      src_const.initFromBox(local.bbox);
      local.getSourceImplicit(FLAG_ALL, src_linear, src_const);

      FOR_BBOX3(p, local.bbox)
      {
        // solution to local affine linear equation
        const T loc_src_const = src_const(p),
                loc_src_linear = src_linear(p),
                loc_rhs        = rhs(p);
        const T loc_lhs = (loc_rhs - operator_coeff * loc_src_const) /
                            (identity_coeff + operator_coeff * loc_src_linear);
        lhs(p) = loc_lhs;
        myAssert(loc_lhs > -1 && loc_lhs < 1.e13);
      }
    }
  }
}
#endif




//-------   interstitial-cell interaction code comes here ----------
/** @brief diffusion reaction
 * most simple case:
 * we have 2 fields that can interact, 
 * e.g. CONC_EXTRACELLULAR and CONC_CELL
 * in this example the matrix has dimension 2x2 and entries
 * which descripe the interaction.
 * of course we need such a matrix for every lattice point in space, so quite a number
 */
typedef Eigen::Matrix<float, NUM_FIELDS, NUM_FIELDS> MatEc;
typedef Eigen::Matrix<float, NUM_FIELDS, 1> VecEc;


void Calculator::calcSlopeEc(const Calculator::State& state, Calculator::State& slope, StepControl& ctrl)
{
  State tmp;
  calcSlopesEc(state, slope, tmp, ctrl);
  slope.addScaled(1., tmp);

  ctrl.euler_dt = std::min(ctrl.euler_dt, euler_dt_interact);
}
/**
 * Diffusion reaction stuff
 */
void Calculator::calcSlopesEcHelper(const State &state, boost::function<void (int, const BBox3 &, const Array3d<MatEc> )> func)
{
  float k12 = params->comprates_k12;
  float k21 = params->comprates_k21;

  // lulz ... this uses Gershgorin's Theorem about eigenvalues.
  max_rate_interact = my::max(std::abs(k12 + my::sign(k12)*std::abs(k21)),
                              std::abs(k21 + my::sign(k21)*std::abs(k12)));
  euler_dt_interact = 0.5/(max_rate_interact + 1.e-13);
  
  #pragma omp parallel
  {
    ThreadLocal local(this, state.field[ID_CONC_EXTRACELLULAR]);
    Array3d<MatEc> clinear;

    #pragma omp for schedule(dynamic, 1)
    for (int boxindex=0; boxindex<mtboxes->size(); ++boxindex)
    {
      local.init(boxindex);
      clinear.initFromBox(local.bbox);
      
      //----- reaction part -----
      FOR_BBOX3(p, local.bbox)// for every point on the sublattice in the current thread
      {
	//get the diffusion reaction matrix
        MatEc& m = clinear(p);
	//get fraction from different fields
        float phi_cells = local.phi_cells(p) + 1.e-13;
        float phi_lm    = local.phi_ecm(p) + local.phi_water(p) + 1.e-13;
        // the division is there because of the conc is the volume average
        // but we need the component average --> makes sense
	// especially: this make the reaction dependent on the local concentration of substances!
        m(ID_CONC_EXTRACELLULAR, ID_CONC_EXTRACELLULAR) = -k12 / phi_lm;
        m(ID_CONC_EXTRACELLULAR, ID_CONC_CELL) = k21 / phi_cells;
        m(ID_CONC_CELL, ID_CONC_EXTRACELLULAR) = k12 / phi_lm;
        m(ID_CONC_CELL, ID_CONC_CELL) = -k21 / phi_cells;
	/** every entry of m is multiplied the the fraction present
	 * the product is due to the contact surface area
	 */
        m *= phi_cells*phi_lm;
      }
      //---- diffusion part ----
      // clinear was filled within the reaction part
      func(boxindex, local.bbox, clinear);
    }
  }
}

/** together with
 */
static void explicitSlopeEcCallback(int boxindex, const BBox3 &bbox, const Array3d<MatEc> clinear, const State &state, State &slope)
{
  FOR_BBOX3(p, bbox)
  {
    // fucking stupid but ...
    VecEc ds, s;
    for (int i=0; i<NUM_FIELDS; ++i)
      s[i] = state.field[i](p);
    ds = clinear(p) * s; // just a matrix multiplication s' = A * s,  A contains the rates
    for (int i=0; i<NUM_FIELDS; ++i)
      slope.field[i](p) = ds[i];
  }
}


static void implicitSlopeEcCallback(int boxindex, const BBox3 &bbox, const Array3d<MatEc> clinear, State& lhs, const State& rhs, double identity_coeff, double operator_coeff)
{
  FOR_BBOX3(p, bbox)
  {
    // fucking stupid but ...
    VecEc ds, s;
    for (int i=0; i<NUM_FIELDS; ++i)
      s[i] = rhs.field[i](p);
    MatEc m = operator_coeff * clinear(p);
    m += identity_coeff * MatEc::Identity();
    m = m.inverse().eval();
    ds = m * s;
    for (int i=0; i<NUM_FIELDS; ++i)
      lhs.field[i](p) = ds[i];
  }
}

/**
 * some how the slopes are calculated with an implicit and an explicit
 * method
 */
void Calculator::calcSlopesEc(const State& state, State& slope_expl, State& slope_impl, StepControl& ctrl)
{
  slope_expl.init(NUM_FIELDS, grid->ld, grid->dim, 0);
  slope_impl.init(NUM_FIELDS, grid->ld, grid->dim, 0);

  calcSlopesEcHelper(state, boost::bind(&explicitSlopeEcCallback, _1, _2, _3, boost::cref(state), boost::ref(slope_impl)));

  ctrl.euler_dt = 0.5 * std::numeric_limits<double>::max();
}


void Calculator::invertImplicitOperatorEc(State& lhs, const State& rhs, double identity_coeff, double operator_coeff, StepControl& ctrl, State& extrapolation)
{
  calcSlopesEcHelper(extrapolation, boost::bind(&implicitSlopeEcCallback, _1, _2, _3, boost::ref(lhs), boost::cref(rhs), identity_coeff, operator_coeff));
}

/** @brief this is handeled to the time stepper
 */
StepControl Calculator::doSplitStep(int num, State& state, const StepControl& ctrl)
{
  if (num == 0)// propergate IFF
  {
    return stepper_interstitial->doStep(state.field[0], ctrl);
  }
  else// propergate drug
  {
    myAssert(num == 1);
    return stepper_interaction->doStep(state, ctrl);
  }
}

/**
 * @brief push model forward in time
 */
Steppers::StepControl Calculator::doStep(State& state, const Steppers::StepControl& ctrl)
{
  my::Time t_;
#if 0
  time = ctrl.t;
  Steppers::StepControl res = extracellular_model.doStep(state, ctrl);
  time = res.t;
#ifdef DEBUG
  FOR_BBOX3(p, grid->ld.Box())
  {
    float c = state(p);
    myAssert(c > -0.1 && c < params->inject_max * 1.1);
  }
#endif
  return res;
#endif
  my::LogScope logscope(my::log(), "drug: ");
  
  time = ctrl.t;
  Steppers::StepControl res = doStrangeSplittingSteps<State>(state, ctrl, 2, boost::bind(&Calculator::doSplitStep, this, _1, _2, _3));
  time = res.t;

#ifdef DEBUG
//   FOR_BBOX3(p, grid->ld.Box())
//   {
//     for (int i=0; i<state.num_fields; ++i)
//     {
//       float c = state.field[i](p);
//       myAssert(c > -0.01);
//     }
//   }
  //cout << "conc_e = " << state.field[0].valueStatistics() << endl;
  //cout << "conc_c = " << state.field[1].valueStatistics() << endl;
#endif
  cout << format("t=%f, euler_dt=%f | dt_kdiff=%f, dt_src_impl=%f, dt_vel=%f, dt_k12=%f")
    % ctrl.t
    % res.euler_dt
    % euler_dt_kdiff
    //% euler_dt_srcexpl
    % euler_dt_srcimpl
    % euler_dt_vel
    % euler_dt_interact << endl;
  cout << format("step time: %s") % (my::Time()-t_) << endl;
  return res;
}
/** choose the injetion method here */
Calculator::T Calculator::GetInflowVesselConc() const
{
  T conc_in = 0;
  switch(params->inject_mode)
  {
    case Params::DF_INJECT_MODE_EXP:
      conc_in = std::exp(-time*std::log(2.0f)/params->inject_t);
      break;
    case Params::DF_INJECT_MODE_JUMP:
      if(time<params->inject_t) conc_in = 1.;
      break;
    case Params::DF_INJECT_MODE_JUMP_UP:
      if (time>params->inject_t) conc_in = 1.;
      break;
  }
  return conc_in * params->inject_max;
}


#if 0
void Calculator::CalculateSources(int boxindex, const BBox3& bbox, const State &state, int flags, Array3d< T > res) const
{
  ConstArray3d<T> conc = state.field[ID_CONC_EXTRACELLULAR];
  Array3d<T> clinear(bbox), cconst(bbox);
  getSourceImplicit(flags, boxindex, bbox, conc, clinear, cconst);
  FOR_BBOX3(p, bbox)
    res(p) = clinear(p)*conc(p) + cconst(p);
}
#endif


void Calculator::writeH5(H5::H5File f, H5::Group g, const State& state, double t, H5::Group ld_group) const
{
  const Calculator& drugcalc = *this;
  
  
//   h5::Attributes a = g.attrs();
//   a.set("DRUG_INFLOW_CONC", drugcalc.GetInflowVesselConc());

  //writeAttrToH5<T>(g, string("DRUG_INFLOW_CONC"),drugcalc.GetInflowVesselConc() );
  const LatticeDataQuad3d &ld = grid->ld;

  Array3d<T> conc_total(ld.Box());
  //Array3d<T> uptake_field(ld.Box());
  Array3d<T> src_diff_vess(ld.Box());
  Array3d<T> src_flow_vess(ld.Box());
  Array3d<T> src_lymph(ld.Box());
  Array3d<T> sink_flow_vess(ld.Box());

  #pragma omp parallel for schedule(dynamic, 1)
  for (int boxindex=0; boxindex<mtboxes->size(); ++boxindex)
  {
    ThreadLocal local(this, state.field[ID_CONC_EXTRACELLULAR]);
    local.init(boxindex);
    
    //local.calculateSources(IfDrug::Calculator::FLAG_UPTAKE, uptake_field);
    local.calculateSources(IfDrug::Calculator::FLAG_SRC_DIFF_VESS, src_diff_vess);
    local.calculateSources(IfDrug::Calculator::FLAG_SRC_FLOW_VESS, src_flow_vess);
    local.calculateSources(IfDrug::Calculator::FLAG_SINK_FLOW_LYMPH, src_lymph);
    local.calculateSources(IfDrug::Calculator::FLAG_SINK_FLOW_VESS, sink_flow_vess);

    FOR_BBOX3(p, local.bbox)
    {
      conc_total(p) = state.field[ID_CONC_EXTRACELLULAR](p) + state.field[ID_CONC_CELL](p);
    }
  }

  //h5cpp::Datatype disktype = h5cpp::get_disktype<float>();
  H5::DataType disktype = H5::PredType::NATIVE_FLOAT;
//   WriteScalarField(g, "conc", conc_total, ld, ld_group, disktype);
//   WriteScalarField(g, "conc_ex", state.field[ID_CONC_EXTRACELLULAR], ld, ld_group, disktype);
//   WriteScalarField(g, "conc_cell", state.field[ID_CONC_CELL], ld, ld_group, disktype);
//   //WriteScalarField(g, "uptake", uptake_field, ld, ld_group, disktype);
//   WriteScalarField(g, "src_diff_vess", src_diff_vess, ld, ld_group, disktype);
//   WriteScalarField(g, "src_flow_vess", src_flow_vess, ld, ld_group, disktype);
//   WriteScalarField(g, "src_lymph", src_lymph, ld, ld_group, disktype);
  
  
  WriteScalarField(g, "conc", conc_total, ld, ld_group);
  WriteScalarField(g, "conc_ex", state.field[ID_CONC_EXTRACELLULAR], ld, ld_group);
  WriteScalarField(g, "conc_cell", state.field[ID_CONC_CELL], ld, ld_group);
  //WriteScalarField(g, "uptake", uptake_field, ld, ld_group, disktype);
  WriteScalarField(g, "src_diff_vess", src_diff_vess, ld, ld_group);
  WriteScalarField(g, "src_flow_vess", src_flow_vess, ld, ld_group);
  WriteScalarField(g, "src_lymph", src_lymph, ld, ld_group);
}



}

#endif
