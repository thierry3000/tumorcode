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

#include "calc-ifflow.h"

//#ifdef USE_IFFLOW

#include "mwlib/lattice-data.h"
#include "mwlib/field.h"
#include "shared-objects.h"
#include "mwlib/histogram.h"
#include "continuum-flow.h"
#include "mwlib/field.h"

#include "mwlib/log.h"

#include "trilinos_linsys_construction.h"
#include "convection_diffusion_solver.h"
#include "calc-ifdrug.h"



#define H_TO_SEC 3600.
#define SEC_TO_H 0.00027777777777777778

// atmospheric pressure = 100 kPa

/* parameters from

  levick 1987
  FL_COND_TISSUE = 4.5 e-12 cm^4/dyn/s = 4.5 micron^2/kPa/s

  Jain 2007:
  FL_COND_TUMOR = FL_COND_TISSUE = 2.5 E-7 cm^2  /s /mmHg = 187 micron^2/kPa/s

  Zhao, Salmon, Sarntinoranont (2007) Effect of heterogeneous 
  vasculature on interstitial transport within a solid tumor

    Baxter and Jain 1989
  FL_COND_TISSUE = 6.39 micron^2/kPa_s  |   8.53 E-9  cm^2 / mmHg s
  FL_COND_TUMOR  = 30.9 micron^2/kPa_s  |   4.13 E-8 cm^2 / mmHg s
  PERM_CAPI_TISSUE = 2.7E-3 micron/kPa_s | 0.36 E-7 cm /mmHg s
  PERM_CAPI_TUMOR  = 2.1E-2             | 2.8 E-7 cm / mmHg s
  POROSITY = 0.28

  P_VESS = 2.08 kPa  | 15.6 mmHg
  PI_VESS = 2.67 kPa |  20 mmHg
  PI_INT  = 1.33 kPa |  10 mmHg
  SIGMA_T_TISSUE = 0.91
  SIGMA_T_TUMOR = 0.82 

  Walter F., PhD. Boron. Medical Physiology: A Cellular And Molecular Approaoch. Elsevier/Saunders
  P_VESS_ARTERY     =  4.67 kPa   |  35 mmHg
  P_VESS_VEIN       =  2. kPa     |  15 mmHg
  P_INT             =  -0.27 kPa  | -2 mmHg (interstitial pressure!)
  SIGMA_T_X_PI_VESS        = 3.73 kPa    |  28 mmHg
  SIGMA_T_X_PI_INT_ARTERY  = 0.01 kPa    |  0.1 mmHg
  SIGMA_T_X_PI_INT_VEIN    = 0.4 kPa     |  3 mmHg

  
    Foy and Blake 2001
  DF_KDIFF = 6.0E3 micron^2/s (Gd-DTPA) 
    Chary and Jain 1989
  DF_KDIFF = 1.3E-1 micron^2/s (albumin)

  using the values O2 diffusivity from sinek_2009 and scaling by molecular weight ratio sqrt(18/532)
  DF_KDIFF(dox) = 180 micron^2/s
  they actually use 16 micron^2/s because it is that low due to ecm binding

  fun fact: diffusion coefficient should scale with 1/sqrt(m), where m is the mass of the particle (in gases, but we)
  
  Graff BA, Bjrnaes I, Rofstad EK (2001) Microvascular permeability of human melanoma
  xenografts to macromolecules: relationships to tumor volumetric growth rate, tumor angiogenesis, and VEGF expression.  
  DF_PERM_CAPI = 1E-6 cm/s = 1E-2 micron/s
  
  DF_INJECT_T = 3h (albumin), 1.6h (Gd-DTPA) // half life
*/  

/*
K = k * c / ( c + c' )
where, k is the maximum rate when c goes to infinity, and c' is the concentration where K = 1/2 k.
The slope at c=0 is k/c'.
*/

IffParams::IffParams()
{
  ecm_friction = 0.05;
  cell_friction = 0.5;
  iff_capillary_permeability_tissue = 2.7e-3 / .2;
  iff_capillary_permeability_tumor  = 2.1e-2 / .2;

  iff_cond_tissue = 6.4 / 0.2;
  iff_cond_tumor = 30.9 / 0.2;
  iff_cond_necro = 30.9 / 0.2;
  
  lymph_surf_per_vol = 6.0 * 2.0*(my::mconst::pi())/(1000.0); //6.0 * 2.0*(PI*10.*100.)/(100.0*100.0*100.0);
  lymph_surf_per_vol_tumor = 0.;
  lymph_press     = 0.0f;
  lymph_permeability      = 2.0/1000.0f; // similar to normal capillary permeability
  // @ permeability 2./1000 we have for a 10^2 micron patch with deltaP = 3 kPa a
  // volume outflow of 0.6 micron^3 / s. Way too much, no?

  debugOut = false;
  coupled_flow = false;

  capillary_oncotic_pressure = 2.7;
  interstitial_oncotic_pressure = 1.33;
  osmotic_reflection_coeff = 0.9;
  osmotic_reflection_coeff_tumor = 0.;
  better_permeability_model = true;
}


void IffParams::assign(const ptree &pt)
{
  #define DOPT(name) boost::property_tree::get(name, #name, pt)
  DOPT(ecm_friction);
  DOPT(cell_friction);
  DOPT(iff_cond_tissue);
  DOPT(iff_cond_tumor);
  DOPT(iff_cond_necro);
  DOPT(iff_capillary_permeability_tissue);
  DOPT(iff_capillary_permeability_tumor);
  DOPT(lymph_surf_per_vol);
  DOPT(lymph_surf_per_vol_tumor);
  DOPT(lymph_press);
  DOPT(lymph_permeability);
  DOPT(debugOut);
  DOPT(capillary_oncotic_pressure);
  DOPT(interstitial_oncotic_pressure);
  DOPT(osmotic_reflection_coeff);
  DOPT(better_permeability_model);
  DOPT(osmotic_reflection_coeff_tumor);
  //DOPT(coupled_flow);
  #undef DOPT
}

ptree IffParams::as_ptree() const
{
  boost::property_tree::ptree pt;
  #define DOPT(name) pt.put(#name, name)
  DOPT(ecm_friction);
  DOPT(cell_friction);
  DOPT(iff_cond_tissue);
  DOPT(iff_cond_tumor);
  DOPT(iff_cond_necro);
  DOPT(iff_capillary_permeability_tissue);
  DOPT(iff_capillary_permeability_tumor);
  DOPT(lymph_surf_per_vol);
  DOPT(lymph_surf_per_vol_tumor);
  DOPT(lymph_press);
  DOPT(lymph_permeability);
  DOPT(debugOut);
  DOPT(capillary_oncotic_pressure);
  DOPT(interstitial_oncotic_pressure);
  DOPT(osmotic_reflection_coeff);
  DOPT(better_permeability_model);
  DOPT(osmotic_reflection_coeff_tumor);
  DOPT(coupled_flow);
#undef DOPT
  return pt;
}





void VesselExavasationModel::Init(const VesselList3d &vl_, const ContinuumGrid& grid_, const VesselsInBoxes& vesselsinboxes_)
{
  grid = &grid_;
  vl = &vl_;
  vess_ld = &vl->Ld();
  vesselsinboxes = &vesselsinboxes_;
}



void VesselExavasationModel::GetValues(int boxindex, const BBox3 &bbox, VirtualGridFunctions<double, 2>::List &values) const
{
  values[0][bbox].fill(0);
  values[1][bbox].fill(0);

  CylinderNetworkSampler sampler; sampler.Init(grid->ld.Scale(), ptree());
  for (int j=0; j<(*vesselsinboxes)[boxindex].size(); ++j)
  {
    const Vessel* v = (*vesselsinboxes)[boxindex][j];

    //if(!GetVesselInfo(v,vinfo)) continue;
    double coeffa=0, coeffb=0, vala=0, valb=0;
    if (!GetVesselInfo(v, coeffa, coeffb, vala, valb)) continue;

    Set(sampler, v, *vess_ld);
    for(int s=sampler.GenerateSurfaceSamples()-1; s>=0; --s)
    {
      CylinderNetworkSampler::Sample smpl = sampler.GetSample(s);

      //double press = my::lerp<double, double>(smpl.fpos[2], v->NodeA()->press, v->NodeB()->press); // linear interpolation
      //double coeff = sampler.weight_per_volume * vinfo.wall_conductivity;

      double val = my::lerp<double, double>(smpl.fpos[2], vala, valb); // linear interpolation
      double coeff = my::lerp<double, double>(smpl.fpos[2], coeffa, coeffb); // linear interpolation
      coeff *= sampler.weight_per_volume;

#if 0
      Int3 lpos = grid->ld.WorldToLattice(smpl.wpos);
      if (bbox.Overlaps(lpos))
      {
        values[IFF_LINEAR_COEFF](lpos) -= coeff;
        values[IFF_CONST_COEFF](lpos) += coeff*val;
      }
#else
      AddSmoothDelta<double>(values[IFF_LINEAR_COEFF], bbox, grid->ld, grid->dim, smpl.wpos, -coeff);
      AddSmoothDelta<double>(values[IFF_CONST_COEFF], bbox, grid->ld, grid->dim, smpl.wpos, coeff*val);
#endif
    }
  }
}






static const PdeReal maturation_at_r5 = GetInitialThickness(5.0f);
static const PdeReal maturation_at_r20= GetInitialThickness(20.0f);

void VesselFluidExavasationModel::Init(const VesselList3d& vl_, const ContinuumGrid& grid_, const VesselsInBoxes& vesselsinboxes_, const IffParams& params_)
{
  VesselExavasationModel::Init(vl_, grid_, vesselsinboxes_);
  params = &params_;
}


bool VesselFluidExavasationModel::GetVesselInfo(const Vessel* v, double& coeffa, double& coeffb, double& vala, double& valb) const
{
  if(!v->IsCirculated()) { return false; }
  // mfactor is 0 for m = 0, goes linearly to 1 at m = maturation_at_r5
  //FlReal mfactor = 1.+((v->maturation-maturation_at_r5)/maturation_at_r5);
  //r = std::max(0., params.iff_capillary_permeability_tumor + (params.iff_capillary_permeability_tissue-params.iff_capillary_permeability_tumor)*mfactor);
  double r = 0;
  double r0 = 1./params->iff_capillary_permeability_tumor;
  double r1 = 1./params->iff_capillary_permeability_tissue;
  double sigma = params->osmotic_reflection_coeff;
  double sigma_tum = params->osmotic_reflection_coeff_tumor;
  if (params->better_permeability_model)
  {
    r = std::max<double>(r0, r1 * v->maturation / maturation_at_r5);
    sigma = sigma_tum + my::cut(v->maturation / maturation_at_r5, 0., 1.) * (sigma - sigma_tum);
  }
  else
  {
    r = r0 + (r1 - r0)*(v->maturation)/(maturation_at_r5);
  }
  coeffa = coeffb = 1./r;
  // Starling law
  double osmotic_stuff = sigma * (params->capillary_oncotic_pressure - params->interstitial_oncotic_pressure);
  vala = v->NodeA()->press - osmotic_stuff;
  valb = v->NodeB()->press - osmotic_stuff;
  myAssert(coeffa >= 0.);
  return true;
}



/*------------------------------------------------------
// calculate flow special with wall permeability
------------------------------------------------------*/

FlReal CalcPrecisePressure2( FlReal r, FlReal l, FlReal viscos, FlReal perm, float apos, FlReal p1, FlReal p2, FlReal extp1, FlReal extp2);
void CalcSpecialFlowCoeff2( FlReal r, FlReal l, FlReal visc, FlReal perm, FlReal fpos, FlReal &c_1, FlReal &c_2, FlReal &c_pe1, FlReal &c_pe2 );
FlReal CalcFlowRate2( FlReal r, FlReal l, FlReal visc, FlReal perm, FlReal fpos, FlReal p1, FlReal p2, FlReal pext1, FlReal pext2 );


/*------------------------------------------------------
------------------------------------------------------*/


void IfFlowState::Init(const BBox3 &bbox, int dim, int flags)
{
  if (flags & PRESSURE)
    pfield = MakeArray3dWithBorder<double>(bbox, dim, 1);

  if (flags & SOURCE_DETAIL) {
    vessel_inflow_field.initFromBox(bbox);
    vessel_outflow_field.initFromBox(bbox);
    lymph_outflow_field.initFromBox(bbox);
    flags |= SOURCE;
  }
  if (flags & SOURCE)
    source_field.initFromBox(bbox);

  if (flags & VELOCITIES)
  {
    LatticeIndexRanges ir = LatticeIndexRanges::FromCellRange(bbox, dim);
    for (int i=0; i<dim; ++i)
    {
      velocities[i].initFromBox(ir.faces[i]);
    }
  }
  this->flags = flags;
  this->max_velocity = 0.;
}



/*------------------------------------------------------
------------------------------------------------------*/

enum {
  IFF_THETA_TUMOR = 0,
  IFF_THETA_NECRO = 1,
  IFF_PHI_WATER = 2,
};


struct IffVelocityConductivityFaceVarFunctor
{
  // v = c * grad p
  ConstArray3d<double> conductivity, phi_water;
  
  IffVelocityConductivityFaceVarFunctor(ConstArray3d<double> &conductivity_, ConstArray3d<double> &phi_water_)
    : conductivity(conductivity_), phi_water(phi_water_) {}

  /** @brief: Compute conductivity coefficient c for connection between neighboring lattice sites
   *  p and q. Also return water volume fraction in l0 and l1 (?).
   **/
  void get_stuff(const Int3 &p, const Int3 &q, int axis, int side, double &c, float &l0, float &l1) const
  {
    //Int3 q(p); --q[axis];
    double k0, k1;
    // p and q reference neighboring lattice sites.
    k0 = conductivity(p);
    k1 = conductivity(q);
    l0 = phi_water(p);
    l1 = phi_water(q);
    c = CalcInterfaceCondCoeff(k0, k1);
  }

  double operator()(const Int3 &p, const Int3 &nbp, int axis, int side) const
  {
    double c; float l0, l1;
    get_stuff(p, nbp, axis, side, c, l0, l1);
    return c;
  }
};

struct IffDiffusionConductivityFaceVarFunctor
{
  // div c l grad p = ...
  IffVelocityConductivityFaceVarFunctor f;
  
  IffDiffusionConductivityFaceVarFunctor(ConstArray3d<double> &conductivity_, ConstArray3d<double> &phi_water_) : f(conductivity_, phi_water_) {}
  
  //double operator()(int axis, const Int3 &p) const
  double operator()(const Int3 &p, const Int3 &nbp, int axis, int side) const
  {
    //Int3 fp(p); fp[axis]+=side;
    double c; float l0, l1;
    f.get_stuff(p, nbp, axis, side, c, l0, l1);
    return c*0.5*(l0+l1);
  }
};


/*------------------------------------------------------
------------------------------------------------------*/


void UncoupledIfModel::Init(const ContinuumGrid &grid_,
                            const DynArray<BBox3> &mtboxes_,
                            const VirtualGridFunctions<double,2> &exavasation_model_,
                            const VirtualGridFunctions<double,2> &lymph_model_,
                            const VirtualGridFunctions<double,2> &tissue_model_,
                            const IffParams &params_)
{
  grid = &grid_;
  mtboxes = &mtboxes_;
  exavasation_model = &exavasation_model_;
  lymph_model = &lymph_model_;
  tissue_model = &tissue_model_;
  params = &params_;

  pfield = MakeArray3dWithBorder<double>(grid->ld.Box(), grid->dim, 1);
}



void UncoupledIfModel::insertFieldCoefficients(int box_index, const BBox3& bbox, FiniteVolumeMatrixBuilder& mb)
{
  Array3d<double> conductivities(ExtendForDim(bbox, grid->dim, 1)),
                 phi_water(ExtendForDim(bbox, grid->dim, 1));
  VirtualGridFunctions<double,2>::List arr1 = {{ conductivities, phi_water }};
  
  tissue_model->GetValues(box_index, conductivities.getBox(), arr1);

  mb.AddDiffusion2(bbox, grid->dim, IffDiffusionConductivityFaceVarFunctor(conductivities, phi_water), 1.); // diffusion term

  // vessel sources
  VirtualGridFunctions<double,2>::List coeff = {{ Array3d<double>(bbox), Array3d<double>(bbox) }};
  exavasation_model->GetValues(box_index, bbox, coeff);

  FOR_BBOX3(p, bbox)
    mb.AddLocally(p, phi_water(p)*coeff[IFF_LINEAR_COEFF](p), -phi_water(p)*coeff[IFF_CONST_COEFF](p));

  // lymphatic sources
  lymph_model->GetValues(box_index, bbox, coeff);

  FOR_BBOX3(p, bbox)
    mb.AddLocally(p, coeff[IFF_LINEAR_COEFF](p), -coeff[IFF_CONST_COEFF](p));
}




void UncoupledIfModel::calculate()
{
  my::log().push("pfield: ");
  //ptree pt = make_ptree("output", 1)("system_scaling_factor", -1.)("max_iter",200);
  ptree pt = make_ptree	("preconditioner","multigrid")
			("verbosity", "normal")
			("use_smoothed_aggregation", false)
			("max_iter", 500)
			("conv","rhs")
			("max_resid",1.e-8)
			("system_scaling_factor", -1.);
  /*
   * calling diffusion solver with member variables from class UncoupledIfModel
   * 
   * _1: box_index
   * _2: bbox
   * _3: finite volume matrix builder
   */
  StationaryDiffusionSolve(grid->ld, *mtboxes, grid->dim, pfield,
                           boost::bind(&UncoupledIfModel::insertFieldCoefficients, this, _1, _2, _3),
                           pt);

  CopyBorder(pfield, grid->dim, 1);
  my::log().pop();
}






void UncoupledIfModel::fillOutState(IfFlowState& state, int box_index, tbb::spin_mutex& mutex)
{
  const BBox3 bbox = (*mtboxes)[box_index];
  
  if (state.flags & IfFlowState::PRESSURE)
  {
    state.pfield[bbox].fill(pfield[bbox]);
  }
  
  if (state.flags & (IfFlowState::SOURCE | IfFlowState::SOURCE_DETAIL))
  {
    Array3d<double> vessel_inflow_field, vessel_outflow_field, lymph_outflow_field;
    vessel_inflow_field.initFromBox(bbox);
    vessel_outflow_field.initFromBox(bbox);
    lymph_outflow_field.initFromBox(bbox);

    Array3d<double> conductivities(ExtendForDim(bbox, grid->dim, 1)),
                  phi_water(ExtendForDim(bbox, grid->dim, 1));
    VirtualGridFunctions<double,2>::List arr1 = {{ conductivities, phi_water }};
    tissue_model->GetValues(box_index, conductivities.getBox(), arr1);
    
    // vessel sources
    VirtualGridFunctions<double, 2>::List coeff = {{ Array3d<double>(bbox), Array3d<double>(bbox) }};
    exavasation_model->GetValues(box_index, bbox, coeff);

    FOR_BBOX3(p, bbox)
    {
      double q = phi_water(p)*(coeff[IFF_LINEAR_COEFF](p)*pfield(p) + coeff[IFF_CONST_COEFF](p));
      if (q > 0.)
        vessel_inflow_field(p) = q;
      else
        vessel_outflow_field(p) = q;
    }

    lymph_model->GetValues(box_index, bbox, coeff);
    FOR_BBOX3(p, bbox)
    {
      double c = coeff[IFF_LINEAR_COEFF](p), value = coeff[IFF_CONST_COEFF](p);
      double q = value + c*pfield(p);
      myAssert(q <= 0.);
      lymph_outflow_field(p) = q < 0 ? q : 0.;
    }

    state.source_field[bbox] += vessel_inflow_field[bbox];
    state.source_field[bbox] += vessel_outflow_field[bbox];
    state.source_field[bbox] += lymph_outflow_field[bbox];
    if (state.flags & IfFlowState::SOURCE_DETAIL)
    {
      state.vessel_inflow_field[bbox].fill(vessel_inflow_field[bbox]);
      state.vessel_outflow_field[bbox].fill(vessel_outflow_field[bbox]);
      state.lymph_outflow_field[bbox].fill(lymph_outflow_field[bbox]);
    }
  }

  if (state.flags & IfFlowState::VELOCITIES)
  {
    Array3d<double> conductivities(ExtendForDim(bbox, grid->dim, 1)),
                    phi_water(ExtendForDim(bbox, grid->dim, 1));
    VirtualGridFunctions<double, 2>::List arr1 = {{ conductivities, phi_water }};
  
    tissue_model->GetValues(box_index, conductivities.getBox(), arr1);
    
    ptree pt = make_ptree("max_vel", 0.);
    for (int axis=0; axis<grid->dim; ++axis)
    {
      pt.put("exclude_right_boundary", bbox.max[axis]<grid->ld.Box().max[axis]);
      AddVelocitiesFromPotential2(bbox,
                                  axis,
                                  grid->dim,
                                  grid->ld.Scale(),
                                  IffVelocityConductivityFaceVarFunctor(conductivities, phi_water),
                                  pfield,
                                  state.velocities[axis], pt);
    }
    mutex.lock();
    state.max_velocity = std::max(state.max_velocity, pt.get<double>("max_vel"));
    mutex.unlock();
  }
}




void UncoupledIfModel::getOutState(IfFlowState &state, int flags)
{
  my::Time t_;
  state.Init(grid->ld.Box(), grid->dim, flags);

  tbb::spin_mutex mutex;
  #pragma omp parallel for schedule(dynamic, 1)
  for (int box_index=0; box_index<mtboxes->size(); ++box_index)
  {
    fillOutState(state, box_index, mutex);
  }
  cout << "getOutState: " << (my::Time()-t_).to_ms() << " ms" << endl;
  //ifstate->velocity_cfl_dt = 0.5*field_ld->Scale()/(dim*ifstate->max_velocity+1.e-13);
  if (params->debugOut)
  {
    std::vector<std::vector<Image> > images(4, std::vector<Image>(std::max<int>(3,grid->dim)));
    Image img;
    BBox3 bb = BBox3(grid->ir.cells).Set(2, (grid->ir.cells.min[2] + grid->ir.cells.max[2])/2);

    DrawArray(img, state.source_field[bb], DrawArrayOpts().title("src").outputRange());
    images[0][0] = img;

    DrawArray(img, pfield[bb], DrawArrayOpts().title("press").outputRange());
    images[0][1] = img;

    for (int i=0; i<grid->dim; ++i)
    {
      BBox3 bb = BBox3(grid->ir.faces[i]).Set(2, (grid->ir.faces[i].min[2] + grid->ir.faces[i].max[2])/2);

      DrawArray(img, state.velocities[i][bb], DrawArrayOpts().title(str(format("vel_%i") % i)).outputRange());
      images[1][i] = img;
    }

    DrawArray(img, state.vessel_inflow_field[bb], DrawArrayOpts().title("v_in").outputRange());
    images[2][0] = img;

    DrawArray(img, state.vessel_outflow_field[bb], DrawArrayOpts().title("v_out").outputRange());
    images[2][1] = img;

    DrawArray(img, state.lymph_outflow_field[bb], DrawArrayOpts().title("l_out").outputRange());
    images[2][2] = img;

//     images[3][0] = DrawArray(ifbase->phi_cells[bb], DrawArrayOpts().title("cells").outputRange());
//     images[3][1] = DrawArray(ifbase->phi_water[bb], DrawArrayOpts().title("water").outputRange());
//     images[3][2] = DrawArray(ifbase->theta_tumor[bb], DrawArrayOpts().title("tumor").outputRange());
    
    Image imgbig;
    DrawImageGrid(imgbig, images);
    imgbig.WritePng("iff-pressure.png");
  }
}



//#endif

