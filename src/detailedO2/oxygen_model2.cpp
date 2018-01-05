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
#include "common/shared-objects.h"
#include "common/continuum-flow.h"
#include "common/time_stepper_utils_new.h"
#include "common/trilinos_linsys_construction.h"
#include "common/vessels3d.h"

#include "mwlib/math_ext.h"

#include <boost/math/tools/roots.hpp>
#include <boost/math/tools/tuple.hpp>
#include <boost/foreach.hpp>

#define OUTPUT_PO2INTEGRATION(x)
#define OUTPUT_PO2MEASURECOMP(x)

#include "oxygen_model2.h"
#include <exception>


#define APPROXIMATE_FEM_TRANSVASCULAR_EXCHANGE_TERMS 1
#define USE_TRANSVASCULAR_EXCHANGE_SMOOTH_DISTRIBUTION_KERNEL 0

namespace DetailedPO2
{
// #define PT_ASSIGN(name) boost::property_tree::get(name, #name, pt)
// #define AS_PTREE(name) pt.put(#name, name);
void Parameters::assign(const ptree &pt)
{
  #define DOPT(name) boost::property_tree::get(name, #name, pt)
  
  DOPT(po2init_r0);
  DOPT(po2init_dr);
  DOPT(po2init_cutoff);
  DOPT(plasma_solubility); 
  ////DOPT(c0);
  //DOPT(S_n);
  DOPT(sat_curve_exponent);
  //DOPT(S_p50);
  DOPT(sat_curve_p50);
  //DOPT(D_tissue);
  //DOPT(kd);
  DOPT(D_plasma);
  DOPT(tissue_solubility);
  //#dM = 0.05,
  DOPT(rd_norm);
  DOPT(rd_tum);
  DOPT(rd_necro);
  DOPT(max_iter);
  DOPT(num_threads);
  DOPT(convergence_tolerance);
  DOPT(axial_integration_step_factor);
  DOPT(debug_zero_o2field);
  //DOPT(grid_lattice_const);
  
  DOPT(michaelis_menten_uptake);
  DOPT(massTransferCoefficientModelNumber);
  DOPT(conductivity_coeff1);
  DOPT(conductivity_coeff2);
  DOPT(conductivity_coeff3);
  DOPT(detailedO2name);
  DOPT(loglevel);
  DetailedPO2::Parameters::UpdateInternalValues();
  #undef DOPT
}
ptree Parameters::as_ptree() const
{
  boost::property_tree::ptree pt;
  
  #define DOPT(name) pt.put(#name, name)
  DOPT(po2init_r0);
  DOPT(po2init_dr);
  DOPT(po2init_cutoff);
  DOPT(plasma_solubility); 
  ////DOPT(c0);
  //DOPT(S_n);
  DOPT(sat_curve_exponent);
  //DOPT(S_p50);
  DOPT(sat_curve_p50);
  DOPT(po2_mmcons_k[TISSUE]);
  DOPT(po2_mmcons_k[TCS]);
  DOPT(po2_mmcons_k[DEAD]);
  //DOPT(D_tissue);
  //DOPT(kd);
  DOPT(D_plasma);
  DOPT(tissue_solubility);
  //#dM = 0.05,
  DOPT(rd_norm);
  DOPT(rd_tum);
  DOPT(rd_necro);
  DOPT(max_iter);
  DOPT(num_threads);
  DOPT(convergence_tolerance);
  DOPT(axial_integration_step_factor);
  DOPT(debug_zero_o2field);
  //DOPT(grid_lattice_const);
  
  DOPT(michaelis_menten_uptake);
  DOPT(massTransferCoefficientModelNumber);
  DOPT(conductivity_coeff1);
  DOPT(conductivity_coeff2);
  DOPT(conductivity_coeff3);
  DOPT(detailedO2name);
  DOPT(loglevel);
  #undef DOPT
  return pt;
}
  
void WriteOutput(H5::Group basegroup,
                 const VesselList3d &vl,
                 const Parameters &params,
                 const boost::optional<const VesselPO2Storage&> vesselpo2,
                 const boost::optional<DynArray<const Vessel*>&> sorted_vessels,
                 const boost::optional<ContinuumGrid&> grid,
                 const boost::optional<Array3df> po2field,
                 const boost::optional<const FiniteVolumeMatrixBuilder&> mbopt);
  

inline double NANd() { return std::numeric_limits<double>::quiet_NaN(); }
inline float NANf() { return std::numeric_limits<float>::quiet_NaN(); }
inline bool isFinite(double x) { return std::isfinite(x); }

/*-------------------------------------------------------------*/
// paramters and basic equation
#if 0
double Parameters::Permeability(double r, double w) const
{
  double t1 = std::log((r+15.)/r)*r; // Mist!!! Wandstärkewert hardcoded. 
  double t2 = po2_kdiff/t1*tissue_solubility; // mlO2 / ml / mmHg * um^2/s
  // -> [j] = um^3 O2 / mmHg / s -> 1.e-12 mlO2 / mmHg / s
  return t2;
}
#endif




/**
 * @brief Saturation of dependent on pressure
 * 
 * This functions calculates the amount of oxygen docked to the RBCs dependent 
 * on the partial pressure. The model assumes a hill equation which is an
 * improvement to the the michaeli menten kinetics.
 **/
double Parameters::Saturation( double p )  const //calculates saturation dependent on the partial pressure,
{
  if (p<=0.) 
    return 0.;
  const double n = sat_curve_exponent;
  const double mu = p/sat_curve_p50;
#ifdef DEBUG
//   std::cout << "p: " << p << std::endl;
//   std::cout << "mu: " << mu << std::endl;
//   std::cout << "sat_curve_exponent: " << sat_curve_exponent << std::endl;
//   std::cout << "sat_curve_p50: " << sat_curve_p50 << std::endl;
#endif
  const double mun = std::pow(mu, n);
  double S = mun/(mun+1.);
#ifdef DEBUG
  //std::cout << "S: " << S << std::endl;
#endif
  return S;
}

std::pair<double,double> Parameters::DiffSaturation( double p )  const // derivative of saturation
{
  if (p<=0.) return std::make_pair(0.,0.);
  const double n = sat_curve_exponent;
  const double mu = p/sat_curve_p50;
  const double mun = std::pow(mu, n);
  const double mun_1 = std::pow(mu, n-1.);
  double S = mun/(mun+1.);
  double dS = n*mun_1/my::sqr(mun+1.)*(1./sat_curve_p50);
  return std::make_pair(S,dS);
}


boost::tuple<double, double, double> Parameters::DiffSaturation2(double p) const  // 2nd derivative
{
  if (p<=0.) return boost::make_tuple(0., 0., 0.);
  const double n = sat_curve_exponent;
  const double mu = p/sat_curve_p50;
  const double mun = std::pow(mu, n);
  const double mun_1 = std::pow(mu, n-1.);
  double S = mun/(mun+1.);
  double dS = n*mun_1/my::sqr(mun+1.);
  //double dS2 = 6.*mu*(1.-2.*mu3)/my::cubed(mu3+1.)/my::sqr(PHalf());
  double dS2 = ((n*mun_1*(n-1.)*(mun+1.)-2.*n*n*mun_1*mun)/mu/my::sqr(mun+1.))*(1./sat_curve_p50);
  return boost::make_tuple(S,dS,dS2);
}


Parameters::Parameters()
{
  sat_curve_exponent = 2.;
  sat_curve_p50 = 38.;
  max_iter = 100;
  axial_integration_step_factor = 0.5;
  debug_zero_o2field = false;
  po2init_r0 = 0.;
  po2init_dr = 0.;
  po2init_cutoff = 0.;
  transvascular_ring_size = 0.5;
  tissue_boundary_condition_flags = 0; // neumann bc
  tissue_boundary_value = 0.;
  po2_mmcons_m0[TISSUE] = 4.5/6.e4;    // corresponds to M_0 in oxygen paper
  po2_mmcons_k[TISSUE]  = 4.; // mmHg  // corresponds to P_{M50} in oxygen paper
  po2_mmcons_m0[TCS] = 4.5/6.e4*2.;
  po2_mmcons_k[TCS]  = 2.; // mmHg
  po2_mmcons_m0[DEAD] = 0.;
  po2_mmcons_k[DEAD]  = 2.; // mmHg
  michaelis_menten_uptake = false;
  plasma_solubility = 3.1e-5;
  tissue_solubility = 2.8e-5;
  haemoglobin_binding_capacity = 0.5; /*mlO2/cm^3*/
  extra_tissue_source_linear = 0.;
  extra_tissue_source_const = 0.;
  loglevel = 1;
  num_threads = 4;
  convergence_tolerance = 1.e-3;
  approximateInsignificantTransvascularFlux = false;
  
  rd_norm = 200.;
  rd_tum = 100.;
  rd_necro = 500.;
  D_plasma = 2000.;
  massTransferCoefficientModelNumber = 0;
  conductivity_coeff1 = 0; // must set this as parameter, no default
  conductivity_coeff2 = 0;
  conductivity_coeff3 = 0;
  detailedO2name = "none";
  
  UpdateInternalValues();
}

void Parameters::UpdateInternalValues()
{
  const double n = sat_curve_exponent;
  double t1 = 2.*n*n-2.;
  double t2 = std::sqrt(3.*n*n*n*n-3.*n*n);
  double t3 = n*n + 2. + 3.*n;
  ds2_zeros_[0] = std::pow((t1+t2)/t3, 1./n)*sat_curve_p50;
  ds2_zeros_[1] = std::pow((t1-t2)/t3, 1./n)*sat_curve_p50;
  double max_r[2] = {
    boost::get<2>(DiffSaturation2(ds2_zeros_[0])),
    boost::get<2>(DiffSaturation2(ds2_zeros_[1]))
  };
  ds2_max_ = std::max(std::abs(max_r[0]), std::abs(max_r[1]));

  conc_neglect_s = -1.;
  double p = sat_curve_p50;
  while (true)
  {
    double f = haemoglobin_binding_capacity*Saturation(p)/(plasma_solubility*p);
    if (f < 1.e-6)
      break;
    p *= 0.5;
  }
  conc_neglect_s  = BloodPO2ToConc(p, 1.);
  SetTissueParamsByDiffusionRadius(D_plasma, tissue_solubility, rd_norm, rd_tum, rd_necro);
}


double Parameters::DiffSaturationMaxRateOfChange(double p) const
{
  return ds2_max_;
}

/* *
 * @brief Eqation (1) in the pdf
 * 
 */
double Parameters::BloodPO2ToConc(double p, double h)  const
{
  return h*haemoglobin_binding_capacity*Saturation(p) +  plasma_solubility*p;  // mlO2/cm^3
}


std::pair<double, double> Parameters::BloodPO2ToHematocritAndPlasmaConc(double p, double h) const
{
  return std::make_pair(haemoglobin_binding_capacity*Saturation(p), plasma_solubility*p);
}

/**
 * @brief Inversion of Eqation (1)
 * 
 * This function numerically inverts eqation (1).
 * Since hill curve is almost constant and zeros for low pressure, we archieve a better
 * stability when neglegting these regime 
 */
double Parameters::ConcToBloodPO2(double conc, double h)  const
{ 
  const double a = h*haemoglobin_binding_capacity, b = plasma_solubility;

  if (conc < conc_neglect_s || (h <= 0.))
    return conc/b;
  
  auto f_and_df = [=](double p) -> boost::math::tuple<double, double>
  {
    double S, dS; boost::tie(S,dS) = DiffSaturation(p);
    double f0 = a*S+b*p - conc;
    double f1 = a*dS+b;
    return boost::math::make_tuple(f0, f1);
  };

  auto f = [=](double p) -> double {
    const double S = Saturation(p);
    return a*S+b*p - conc;
  };

  const double rbound = conc/b*1.001;

  auto tol_bisect = [=](double l, double r) { return (r-l) < 0.001*(r+l); };
  //auto tol_newton = [=](double l, double r) { return (r-l) < 1.e-6*PHalf(); };

#ifdef DEBUG
  {
    double f0 = f(0.);
    double f1 = f(rbound);
    myAssert(f0 < f1);
    myAssert(f0 <= 0. && f1 >= 0.);
  }
#endif
  
  double l, r;
  boost::tie(l,r) = boost::math::tools::bisect(f, 0., rbound, tol_bisect);
  //double res = 0.5*(l+r);
  //boost::uintmax_t max_iter = 3;
  double res = boost::math::tools::newton_raphson_iterate(f_and_df, 0.5*(l+r), l, r, 10); //, max_iter);
  return res;
}

/* *
 * @brief Needed for ratial tumor growth
 * 
 * This sets
 */

void Parameters::SetTissueParamsByDiffusionRadius(double kdiff_, double alpha_, double rdiff_norm_, double rdiff_tum_, double rdiff_necro_)
{
  po2_kdiff = kdiff_; // um^2/s
  po2_cons_coeff[0] = po2_kdiff/my::sqr(rdiff_norm_)*tissue_solubility; // added TissueSolutbility
  po2_cons_coeff[1] = po2_kdiff/my::sqr(rdiff_tum_)*tissue_solubility;
  po2_cons_coeff[2] = po2_kdiff/my::sqr(rdiff_necro_)*tissue_solubility;
}


std::pair<double, double> Parameters::ComputeUptake(double po2, float *tissue_phases, int phases_count) const
{
  double dm_total = 0., m_total=0.;
  if (!michaelis_menten_uptake)
  {
    for (int i=0; i<phases_count; ++i)
    {
      dm_total += po2_cons_coeff[i]*tissue_phases[i];
    }
    m_total = dm_total*po2;
  }
  else
  {
    po2 = std::max(0., po2); // errors can make po2 negative, but we still want to use the consumption rates for at least po2=0
    for (int i=0; i<phases_count; ++i)
    {
      double t1 = 1./(po2+po2_mmcons_k[i]);
      double m  = po2_mmcons_m0[i]*po2*t1;
      double dm = po2_mmcons_m0[i]*po2_mmcons_k[i]*t1*t1; // Fixed formula for first derivative
      dm_total += dm*tissue_phases[i];
      m_total += m*tissue_phases[i];
    }
  }
  return std::make_pair(m_total, dm_total);
}


#if USE_TRANSVASCULAR_EXCHANGE_SMOOTH_DISTRIBUTION_KERNEL

namespace SmallKernelConvolutionInternal
{
  BBox3 ComputeConvolutionBBox(const BBox3 &bbox, const LatticeDataQuad3d &ld, int dim, const Float3 &pos, float filter_radius)
  {
    Int3 ip; Float3 q;
    boost::tie(ip, q) = ld.WorldToFractionalCoordinate(pos);
    BBox3 convbox(Cons::DONT);
    for (int i = 0; i < dim; ++i)
    {
      convbox.min[i] = std::max(bbox.min[i], my::iceil<float>(ip[i]+q[i] - filter_radius));
      convbox.max[i] = std::min(bbox.max[i], my::ifloor<float>(ip[i]+q[i] + filter_radius));
    }
    return convbox;
  }
  
  template<class Func, int dim>
  struct Eval
  {
    static void eval(const BBox3 &bbox, Int3 &ip, Func func)
    {
      for (int i=bbox.min[dim]; i<=bbox.max[dim]; ++i) // iterate over neighbor sites of ip along the current dimension
      {
        ip[dim] = i;  // set the current point
        Eval<Func, dim-1>::eval(bbox, ip, func); // one lower dimension
      }
    }
  };
  
  
  template<class Func>
  struct Eval<Func, -1>
  {
    static void eval(const BBox3 &bbox, const Int3 &ip, Func func)
    {
      func(ip);
    }
  }; 

}


template<class Func>
static void SmallKernelConvolution(const BBox3 &convbox, int dim, Func func)
{
  Int3 p;
  if (dim == 3)
    SmallKernelConvolutionInternal::Eval<Func, 2>::eval(convbox, p, func);
  else if (dim == 2) {
    p[2] = 0;
    SmallKernelConvolutionInternal::Eval<Func, 1>::eval(convbox, p, func);
  } else if (dim == 1) {
    p[2] = p[1] = 0;
    SmallKernelConvolutionInternal::Eval<Func, 0>::eval(convbox, p, func);
  }
}

struct SmoothCosKernel : public boost::noncopyable
{
  float kernel[3][128]; // separable function by dimension, careful with the memory here. Make sure that the kernel width is less than 128 lattice sites!
  BBox3 bbox;
  const float filter_radius;

  SmoothCosKernel(const Float3 &wp, const LatticeDataQuad3d &ld) : filter_radius(ld.Scale()*2.0)
  {
    bbox = SmallKernelConvolutionInternal::ComputeConvolutionBBox(ld.Box(), ld, 3, wp, filter_radius/ld.Scale());
    const Float3 boxOffset = ld.LatticeToWorld(bbox.min);
    for (int i = 0; i < 3; ++i)
    {
      int k = 0;
      for (int q = bbox.min[i]; q <= bbox.max[i]; ++q, ++k)
      {
        float r = std::abs(boxOffset[i] + k*ld.Scale() - wp[i]);  // this is the distance of the current lattice point from the center point of the kernel
        kernel[i][k] = my::smooth_delta_cos<float, float>(r, filter_radius*0.5f);  // for some reason i designed this that the actual interval where that function is nonzero is two times the radius argument
      }
    }
  }
  float operator()(const Int3 &p) const
  {
    float res = kernel[0][p[0] - bbox.min[0]];
    res *= kernel[1][p[1] - bbox.min[1]];
    res *= kernel[2][p[2] - bbox.min[2]];
    return res;
  }
  
  std::ostream& print(std::ostream &os, const LatticeDataQuad3d& ld) const
  {
    const auto & kernel = *this;
    os << "bbox = " << kernel.bbox << endl;
    os << "wbox = " << ld.LatticeToWorld(kernel.bbox.min) << " - " << ld.LatticeToWorld(kernel.bbox.max) << endl;
    //os << "filter_radius = " << kernel.filter_radius << endl;
    float coeffsum = 0.;
    for (int i = 0; i < 3; ++i)
    {
      os << format("kernel[dim = %i] = ") % i;
      int k = 0;
      for (int q = kernel.bbox.min[i]; q <= kernel.bbox.max[i]; ++q, ++k)
      {
        os << kernel.kernel[i][k] << ", ";
        coeffsum += kernel.kernel[i][k];
      }
      os << endl;
    }
    coeffsum *= ld.Scale() / 3.0; // should be = 1
    os << "coefficient sum = " << coeffsum;
    return os;
  }
};
#endif  //USE_TRANSVASCULAR_EXCHANGE_SMOOTH_DISTRIBUTION_KERNEL


double InterpolateField(const Parameters &params, const LatticeDataQuad3d &ld, const Float3 &wp, const Array3df &po2field)
{
#if USE_TRANSVASCULAR_EXCHANGE_SMOOTH_DISTRIBUTION_KERNEL
  SmoothCosKernel kernel(wp, ld);
  float result = 0.f;
  auto func = [&](const Int3 &p)
  {
    result += kernel(p)*po2field(p);
  };
  SmallKernelConvolution<decltype(func)>(kernel.bbox, 3, func);
  result *= my::cubed(ld.Scale());
  return result;
#else
  return FieldInterpolate::Value<float>(po2field, ld, FieldInterpolate::Extrapolate(), wp); 
#endif
}

void AddSourceContributionsTo(const Parameters &params, const LatticeDataQuad3d &ld, const Float3 &wp, FiniteVolumeMatrixBuilder &matrix_builder, double coeffLinear, double coeffConst)
{
#if USE_TRANSVASCULAR_EXCHANGE_SMOOTH_DISTRIBUTION_KERNEL
  SmoothCosKernel kernel(wp, ld);
  const double w = 1./params.tissue_solubility;
  auto func = [&](const Int3 &p)
  {
    float f = kernel(p)*w;
    matrix_builder.AddLocally(p, -f*coeffLinear, f*coeffConst);
  };
  SmallKernelConvolution<decltype(func)>(kernel.bbox, 3, func);
#else
#if APPROXIMATE_FEM_TRANSVASCULAR_EXCHANGE_TERMS
  const double w = 1./params.tissue_solubility/my::cubed(ld.Scale());
  Int3 ip; Float3 q;
  boost::tie(ip, q) = ld.WorldToFractionalCoordinate(wp);
  FOR_BBOX3(iq, BBox3(ip[0],ip[1],ip[2],ip[0]+1,ip[1]+1,ip[2]+1))
  {
    if (!ld.IsInsideLattice(iq)) continue;
    double f = ((iq[0]==ip[0]) ? (1.-q[0]) : q[0]) *
                ((iq[1]==ip[1]) ? (1.-q[1]) : q[1]) *
                ((iq[2]==ip[2]) ? (1.-q[2]) : q[2]);
    // source = linearCoeff * po2field + constCoeff = S/V*K * (po2vessel - po2field)
    matrix_builder.AddLocally(iq, -w*f*coeffLinear, w*f*coeffConst);
  }
  //printf("adding source @ %f, %f, %f\n", wp[0], wp[1], wp[2]);
#else
  const double w = 1./params.tissue_solubility/my::cubed(ld.Scale());
  Int3 ip; Float3 q;
  boost::tie(ip, q) = ld.WorldToFractionalCoordinate(wp);
  FOR_BBOX3(iq, BBox3(ip[0],ip[1],ip[2],ip[0]+1,ip[1]+1,ip[2]+1))
  {
    if (!ld.IsInsideLattice(iq)) continue;
    double f = ((iq[0]==ip[0]) ? (1.-q[0]) : q[0]) *
                ((iq[1]==ip[1]) ? (1.-q[1]) : q[1]) *
                ((iq[2]==ip[2]) ? (1.-q[2]) : q[2]);
    FOR_BBOX3(ir, BBox3(ip[0], ip[1], ip[2], ip[0]+1, ip[1]+1, ip[2]+1))
    {
      if (!ld.IsInsideLattice(ir)) continue;
      double g = ((ir[0]==ip[0]) ? (1.-q[0]) : q[0]) *
                  ((ir[1]==ip[1]) ? (1.-q[1]) : q[1]) *
                  ((ir[2]==ip[2]) ? (1.-q[2]) : q[2]);
      matrix_builder.Add(ld.LatticeToSite(iq), ld.LatticeToSite(ir), -f*g*w*coeffLinear);
    }
    matrix_builder.AddLocally(iq, 0, w*f*coeffConst);
  }    
#endif
#endif
}


// this should compute 2pi r MTC = pi Nu D alpha
//---------------------------------------------
// Anm Welter 18.11.2015: Vorher stand hier 2pi r MTC = 2 pi nu D alpha,
// mit zwei pi. Die Zwei war falsch. Weil die Nusselt zahl ja ueber den
// Durchmesser definitert ist. Fuer die Ergebnisse im O2 Paper 2016, ist
// dies jedoch zweitranging, da massTransferCoefficientModelNumber = 0 fuer
// alles ausser den "Single-Vessel-Validation" Plots verwendet wurde.

// this should compute 2pi r MTC = 2pi nu D alpha
// with  nu  = p2 (1 - exp(-r/p1)
// MTC = D alpha /2 r * nu
double ComputeCircumferentialMassTransferCoeff(const Parameters &params, double r)
{
  const double p1 = params.conductivity_coeff1;
  const double p2 = params.conductivity_coeff2;
  const double p3 = params.conductivity_coeff3;
  
  if (params.massTransferCoefficientModelNumber == 1)
  {
    const double nusseltNumber = p2*(1.0 - std::exp(-r/p1)) + p3 * r;
    const double kd = params.plasma_solubility*params.D_plasma;
    const double intravascularConductivity = my::mconst::pi()*nusseltNumber*kd;
    //printf("p1=%f, p2=%f, p3=%f, nu=%f, c=%f\n", p1, p2, p3, nusseltNumber, intravascularConductivity);
    return intravascularConductivity;
    
    // with extravascular conductivity
    //const double t1 = my::mconst::pi2()/std::log((r+ld.Scale())/r);
    //const double t1 = my::mconst::pi2()*r / (ld.Scale()*0.5);
    //const double extravascularConductivity = (params.po2_kdiff*params.tissue_solubility)*t1; // mlO2 / ml / mmHg * um^2/s
    //return 1.0/(1.0/intravascularConductivity + 1.0/extravascularConductivity);
  }
  else
  {
    const double p = p1 + std::exp(-r/p2)*p3;
    return my::mconst::pi2()*r*p;
  }
}

/** @brief used for the measurement
 * 
 * \param po2tissue -> local buffer variable,initalized with NANd
 * \param wp position in real world coordinates
 * \param p position in 3d space
 * \param dp gradient
 * \param tissuePo2Field this will be filled
 */
struct ComputeWallFluxesDiffusive : boost::noncopyable
{
  double po2tissue;
  Float3 wp;

  const Float3 p, dp;
  const Parameters &params;
  const LatticeDataQuad3d &ld;
  const double r;
  const TissuePhases &phases;
  const Array3df &tissuePo2Field;

  ComputeWallFluxesDiffusive(const Parameters &params_, const LatticeDataQuad3d &ld_, double r_, const Array3df &tissuePo2Field_, const TissuePhases &phases_, const Float3 &p_, const Float3 &dp_)
    : params(params_), ld(ld_), r(r_), phases(phases_), p(p_), dp(dp_), tissuePo2Field(tissuePo2Field_)
  {
    po2tissue = NANd();
  }

  void StartNewPosition(double x)
  {
      wp = p + x * dp;
      po2tissue = DetailedPO2::InterpolateField(params, ld, wp, tissuePo2Field);
  }

  std::pair<double, double> ComputeFluxes(double po2) const
  {
    const double prefactor = ComputeCircumferentialMassTransferCoeff(params, r);
    const double jtotal = prefactor * (po2 - po2tissue);
    return std::make_pair(jtotal, jtotal);
  }

  void AddSourceContributionsTo(FiniteVolumeMatrixBuilder &matrix_builder, double po2, double lengthWeight) const
  {
    // note: solubility cancels out
    const double prefactor = ComputeCircumferentialMassTransferCoeff(params, r);
    DetailedPO2::AddSourceContributionsTo(params, ld, wp, matrix_builder, -lengthWeight*prefactor, lengthWeight*prefactor*po2);
  }
};



#if 0
struct ComputeWallFluxesWithConsumption : boost::noncopyable
{
  //double m0, dm; // cached consumption rate
  double po2tissue;
  Float3 wp;
  Float3 phases_local;
  
  const Float3 p, dp;
  const Parameters &params;
  const LatticeDataQuad3d &ld;
  const double r;
  const TissuePhases &phases;
  const Array3df &tissuePo2Field;

  ComputeWallFluxesWithConsumption(const Parameters &params_, const LatticeDataQuad3d &ld_, double r_, const Array3df &tissuePo2Field_, const TissuePhases &phases_, const Float3 &p_, const Float3 &dp_)
    : params(params_), ld(ld_), r(r_), phases(phases_), p(p_), dp(dp_), tissuePo2Field(tissuePo2Field_)
  {
    //last_x = std::numeric_limits<double>::quiet_NaN();
  }

  void StartNewPosition(double x)
  {
      wp = p + x * dp;
      // use constant uptake approximation computed for the current tissue po2 level
      po2tissue = DetailedPO2::InterpolateField(params, ld, wp, tissuePo2Field);
      for (int i=0; i<phases.count; ++i)
        phases_local[i] = FieldInterpolate::Value(phases.phase_arrays[i], ld, FieldInterpolate::OutOfDomainExtrapolateT(), wp);
      //tie(m0, dm) = params.ComputeUptake(po2tissue, phases_local.data(), phases.count);
  }

  boost::tuple<double, double, double, double, double, double>
  ComputeFluxCoefficients(double po2)
  {
    double m0, dm;
    tie(m0, dm) = params.ComputeUptake(0.5*(po2tissue+po2), phases_local.data(), phases.count);
    const double xi0 = 1;
    const double xi1 = (r+ld.Scale()*params.transvascular_ring_size)/r;
    const double L   = r;
    const double f0  = po2;
    const double f1  = po2tissue;
    const double b   = -(L*L*m0)/(params.po2_kdiff*params.tissue_solubility);
    const double c   =  (params.po2_kdiff*params.tissue_solubility)/L*r*my::mconst::pi2();
    const double t1 = std::log(xi1);
    const double t2 = 1./t1;
    const double t3 = 1./(xi1 * t1);
    const double k0_inner = -0.25*c*b*t2*(-2.*t1         + xi1*xi1 - 1.);
    const double k0_outer = -0.25*c*b*t2*(-2.*t1*xi1*xi1 + xi1*xi1 - 1.);
    const double kv_inner = c*t2;
    const double kv_outer = c*t3;
    const double kt_inner = -c*t2;
    const double kt_outer = -c*t3;
    return boost::make_tuple(k0_inner, kv_inner, kt_inner, k0_outer, kv_outer, kt_outer); // first 3 coefficients for inner flux, then 3 coefficients for outer flux
  }
  
  std::pair<double, double> ComputeFluxes(double po2)
  {
    double k0inner, kvinner, ktinner, k0outer, kvouter, ktouter;
    tie(k0inner, kvinner, ktinner, k0outer, kvouter, ktouter) = ComputeFluxCoefficients(po2);
    const double jinner = k0inner + kvinner*po2 + ktinner*po2tissue;
    const double jouter = k0outer + kvouter*po2 + ktouter*po2tissue;
    return std::make_pair(jinner, jouter);
  }

  void AddSourceContributionsTo(Array3df srcLinearCoeff, Array3df srcConstantCoeff, double po2, double lengthWeight)
  {
    double k0inner, kvinner, ktinner, k0outer, kvouter, ktouter;
    tie(k0inner, kvinner, ktinner, k0outer, kvouter, ktouter) = ComputeFluxCoefficients(po2);
    const double prefactor = lengthWeight / (params.tissue_solubility);
    k0inner *= prefactor;
    kvinner *= prefactor;
    ktinner *= prefactor;
    k0outer *= prefactor;
    kvouter *= prefactor;
    ktouter *= prefactor;
    DetailedPO2::AddSourceContributionsTo(params, ld, wp, srcLinearCoeff, srcConstantCoeff, ktouter, kvouter*po2 + k0outer);
  }
};
#endif

typedef ComputeWallFluxesDiffusive ComputeRadialFluxes;


static double ComputePO2RateOfChange(const Parameters &params, double po2, double flow_rate, double h, double j_tv)
{
  double S, dS; boost::tie(S,dS) = params.DiffSaturation(po2);
  double t1 = dS*h*params.haemoglobin_binding_capacity+params.plasma_solubility;
  double t2 = 1./t1;
  double t3 = 1./flow_rate;
  double slope = -j_tv*t2*t3;
  return slope;
}


class VascularPO2PropagationModel : boost::noncopyable
{
public:
  // make one numerical integration step of length dx along the vascular tube axis, using the given starting po2 (po2), also update the po2 and x value.
  virtual void doStep(double &po2, double &x, double dx) = 0;
  virtual void reset(double flowRate, double h, ComputeRadialFluxes &computeRadialFluxes) = 0; // the shit with the fluxes is not so nice
  static std::unique_ptr<VascularPO2PropagationModel> allocate(const Parameters &params);
};

/**
 * @brief The integration of equation (5) happens here.
 * -----------------------------------------------
 *  Experimentation with different methods result:
 * -- Cranc-Nicholson: Oszillations, very costly, non-physical values during inveversion of the slope operator
 * -- ImplicitEuler: 1st Order Accurate. Need complicated function inversion, but works with strong coupling / slow flow rates / fast rates of change.
 * -- ExplicitEuler: 1st Order Accurate. Inefficient for strong coupling but hack is possible, where po2 is set to external po2 if strong coupling is detected
 */

class VascularPO2PropagationImplicitEuler  : public VascularPO2PropagationModel
{
public:
  VascularPO2PropagationImplicitEuler(const Parameters &params_) :
    params(params_), flow_rate(NANd()), h(NANd()), computeFlux(nullptr)
  {
    stepper.set_model(this);
  }

  virtual void reset(double flowRate_, double h_, ComputeRadialFluxes &computeFlux_)
  {
    computeFlux = &computeFlux_;
    flow_rate = flowRate_;
    h = h_;
  }

  virtual void doStep(double &po2, double &x, double dx)
  {
    stepper.doStep(po2, x, dx);
  }
  
  
public: // should be private but ImplicitEuler stepper needs some of this stuff
  typedef VascularPO2PropagationImplicitEuler SelfType;
  typedef double State; // needed for NewSteppers::ImplicitEuler
  typedef NewSteppers::ImplicitEuler2<SelfType*,  NewSteppers::Operations<double> > ImplicitEulerStepper;

private:  
  const Parameters &params;
  ComputeRadialFluxes *computeFlux;
  double flow_rate, h, wall_flux_coeff;
  ImplicitEulerStepper stepper;

public:
  void invertImplicitOperator(double &x, const double &rhs, double identity_factor, double operator_factor, double t)
  {
    /* g := -2 pi r K*(p - p_tissue) / [ q*(H c_0 dS/dp(p) + alpha ) ]
     * f := identity_factor * p + operator_factor * g
     * we want to invert f w.r.t. p
     * x = new p
     * rhs = previous p
     * t = axial position for x
     * 
     * we are looking for the value of po2
     */

    computeFlux->StartNewPosition(t);
    const double extpo2 = computeFlux->po2tissue;
    if (flow_rate <= 0.) {  x = extpo2;  return;  }
    
    /*
     * this function is zero where po2 inverts the equation
     */
    auto objective_function = [=](double po2) -> double
    {
      double S, dS; boost::tie(S,dS) = params.DiffSaturation(po2);
      double t1 = dS*h*params.haemoglobin_binding_capacity+params.plasma_solubility;
      t1 *= flow_rate;
      double jinner, jouter;
      tie(jinner, jouter) = computeFlux->ComputeFluxes(po2);
      double t2 = -jinner;
      double res = identity_factor * po2 + operator_factor * t2/t1;
      return res - rhs;
    };

    double lbound = 0.;
    double rbound = std::max(rhs, extpo2);
    //termination condition 
    auto tol_bisect = [=](double l, double r) -> bool { return (r-l) <= 1.e-9*(r+l); };
    
    if (tol_bisect(lbound, rbound)) // protect against (approximately) equal l-and rbounds.
    {
      x = 0.5*(lbound+rbound);
      return;
    }
    
  #ifdef DEBUG
    {
      double f0 = objective_function(lbound);
      double f1 = objective_function(rbound);
      myAssert(f0 < f1);
      myAssert(f0 <= 0. && f1 >= 0.);
    }
  #endif
  
    double l, r;
    boost::tie(l,r) = boost::math::tools::bisect(objective_function, lbound, rbound, tol_bisect);
    x = 0.5*(l+r);
  }
  
  /* *
   * Eqation (6) in pdf
   */
  void calcSlope(const double &po2, double t, double &slope)
  {
    computeFlux->StartNewPosition(t);
    double jinner, jouter;
    tie(jinner, jouter) = computeFlux->ComputeFluxes(po2);
    slope = ComputePO2RateOfChange(params, po2, flow_rate, h, jinner);
  }
};


// just give vessels a constant O2 level, no integration, no change no nothing
class VascularPO2PropagationNone : public VascularPO2PropagationModel
{
public:
  VascularPO2PropagationNone(const Parameters &params_)
  {
  }

  virtual void reset(double flowRate_, double h_, ComputeRadialFluxes &computeFlux_)
  {
  }

  virtual void doStep(double &po2, double &x, double dx)
  {
    x += dx;
  }
};


std::unique_ptr<VascularPO2PropagationModel> VascularPO2PropagationModel::allocate(const Parameters &params)
{
  if (params.approximateInsignificantTransvascularFlux)
    return std::unique_ptr<VascularPO2PropagationNone>(new VascularPO2PropagationNone(params));
  else
    return std::unique_ptr<VascularPO2PropagationImplicitEuler>(new VascularPO2PropagationImplicitEuler(params));
}


/**
 * @brief section 2 in pdf
 * 
 * Uses mass conservation and equality a of in flow and outflow
 * qblood: 	Blood Flow
 * qrbc:	Flow of red blood cells
 * mo2flux:	mean o2 flux
 */
void ComputeVesselO2Conc(const VesselNode* node, const Parameters &params, VesselPO2Storage &vesselpo2, const VesselList3d &vl)
{
  /*better idea the system stays is equilibrium, so that the po2 is the same in all outflow vessels*/
  OUTPUT_PO2INTEGRATION(fprintf(stderr,"calcpo2 node %i: ", node->Index());)
  double qblood = 0., qrbc = 0., mo2flux = 0.;
  int num_inflow_nodes = 0, num_outflow_nodes = 0;
  
  for (int i=0; i<node->Count(); ++i)//loop over incoming edges see eq. (12)
  {
    NodeNeighbor<const VesselNode*, const Vessel*> nb = node->GetNode(i);
    if (!nb.edge->IsCirculated()) continue;
    if (nb.node->press <= node->press) continue;
    ++num_inflow_nodes;
    int side_index = (nb.edge->NodeA() == node) ? 0 : 1;
    double po2_vess = vesselpo2[nb.edge->Index()][side_index];
    myAssert(isFinite(po2_vess));
    const double q = nb.edge->q;
    const double h = nb.edge->hematocrit;
    const double conc = params.BloodPO2ToConc(po2_vess, h);
    qblood += q;
    qrbc += h * q;
    mo2flux += conc*q;//particle or mass flux
    OUTPUT_PO2INTEGRATION(fprintf(stderr,"e%i,", nb.edge->Index());)
  }
  OUTPUT_PO2INTEGRATION(fprintf(stderr,"\n");)

  double po2 = 0.;
  if (qblood > 0.)//as long as there is blood flow
  {
    double heff = qrbc / qblood;//effective hemoglobin flow
    double ceff = mo2flux / qblood;//effective concentration
    po2 = params.ConcToBloodPO2(ceff, heff);//is P^{\tilde} on pdf
    myAssert(po2 <= 1000.);
  }
  
  for (int i=0; i<node->Count(); ++i)//loop over outgoing edges
  {
    NodeNeighbor<const VesselNode*, const Vessel*> nb = node->GetNode(i);
    if (!nb.edge->IsCirculated()) continue;
    if (nb.node->press >= node->press) continue;
    ++num_outflow_nodes;
    int side_index = (nb.edge->NodeA() == node) ? 0 : 1;
    vesselpo2[nb.edge->Index()][side_index] = po2;//write the calculated pressure to the corresponding lattice sites
  }
  //according to mw, exception handling is computational expensive-> therefore only do it in debug
  //maybe none circulated destroy this assert ???
  //myAssert(num_inflow_nodes>0 or num_outflow_nodes>0);
}


/**
 * @brief Transforms information from the VesselList into physical world
 */
boost::tuple<Float3, Float3, float> GetSegmentLineParameters(const VesselList3d &vl, const Vessel* v, const VesselNode *startnode,bool world)
{
  Float3 p_[2];
  if(!world)// old stuff
  {
    p_[0] = vl.Ld().LatticeToWorld(v->LPosA());
    p_[1] = vl.Ld().LatticeToWorld(v->LPosB());
  }
  else// use the new worldpos array
  {
    p_[0] = v->NodeA()->worldpos;
    p_[1] = v->NodeB()->worldpos;
  }
  Float3 p,dp;
  if (startnode == v->NodeA())
  {
    p = p_[0];
    dp = p_[1]-p_[0];
  }
  else
  {
    p = p_[1];
    dp = p_[0]-p_[1];
  }
  float len = dp.norm();
  dp /= len;
  return boost::make_tuple(p, dp, len);
}


template<class Callback>
static void NumericallyIntegrateVesselPO2(const Parameters &params, 
					  const TissuePhases &phases, 
					  const ContinuumGrid &grid, 
					  const Array3df extpo2, 
					  VascularPO2PropagationModel *vascularPropagationModel,
					  const VesselList3d &vl, 
					  const Vessel* v, // the vessel in question
					  const VesselNode *upstream_node, // pointer to upstream node (attached to one of the ends)
					  double po2start, // po2 at vessel inlet
					  Callback &callback,  // called for each integration point, like callback(i, Nsteps+1, x, weight, po2, computeFlux);
					  double &po2end, // return value or po2 at vessel end
					  bool world
					  )
{
    Float3 p, dp; float len;
    boost::tie(p, dp, len) = GetSegmentLineParameters(vl, v, upstream_node,world);
    ComputeRadialFluxes computeFlux(params, grid.ld, v->r, extpo2, phases, p, dp);
    vascularPropagationModel->reset(v->q, v->hematocrit, computeFlux);
    const int Nsteps = std::max<int>(1, 0.5+len/(params.axial_integration_step_factor*grid.ld.Scale()));
    const double dx = len/Nsteps;
    double x = 0.;
    int i = 0;
    double po2 = po2start;
    
    OUTPUT_PO2MEASURECOMP(printf("sim %i, start = %lf\n",v->Index(),po2start);)

    while(true)
    {
      double weight = (i==0 || i==Nsteps) ? dx*0.5 : dx;
      computeFlux.StartNewPosition(x);
      callback(i, Nsteps+1, x, weight, po2, computeFlux);
      if (i == Nsteps) break;
      vascularPropagationModel->doStep(po2, x, dx);
      OUTPUT_PO2MEASURECOMP(printf("x = %lf, po2 = %lf, extpo2 = %f\n", x, po2, lookup_po2field(x));)
      i += 1;
    }
    po2end = po2;
}



/**
 * @brief integrates the po2 along a vessel
 * 
 * The hiarachial ordering of the vessel segments is provided, in the sorted_vessels array.
 * As Starting point of the integration, an initial value dependent on the  
 */
void IntegrateVesselPO2(const Parameters &params, 
			VesselPO2Storage &vesselpo2,
			const VesselList3d &vl, 
			DynArray<const Vessel*> &sorted_vessels,
			DynArray<const VesselNode*> &arterial_roots,
			const ContinuumGrid &grid,
			const Array3df extpo2,
			const TissuePhases &phases,
			FiniteVolumeMatrixBuilder &matrix_builder,
			bool world
 		      )
{
  if (params.loglevel > 0)
    cout << format("Start IntegrateVesselPO2 !") << endl;
  my::Time t_;
  //DynArray<bool> nodal_o2ready(vl.GetNCount(), false);
  
  /** tumor can occupy nodes with higher indexes,
   * not yet present. should we make this a map?.
   */
  boost::unordered_map<int, bool> nodal_o2ready;
  for(int i =0;i<vl.GetNCount();++i)
  {
    const VesselNode *nd = vl.GetNode(i);
    nodal_o2ready[nd->Index()] = false;
  }
  BOOST_FOREACH(const VesselNode* nd, arterial_roots)
  {//loop over all arterial roots, we follow the blood stream starting here
    nodal_o2ready[nd->Index()] = true;
    for (int i=0; i<nd->Count(); ++i)
    {
      const Vessel* v = nd->GetEdge(i);
      int side_index = (v->NodeA() == nd) ? 0 : 1;
      // setting initial value depending on vessel radius
      vesselpo2[v->Index()][side_index] = params.PInit(v->r);
    }
  }
  
  std::unique_ptr<VascularPO2PropagationModel> vascularPO2PropagationModel = VascularPO2PropagationModel::allocate(params);
  my::Averaged<double> stats;
  BOOST_FOREACH(const Vessel* v, sorted_vessels)
  {
    if (!v->IsCirculated()) //setting not perufded vessels to zero.
    {
      vesselpo2[v->Index()][0] = 0.;
      vesselpo2[v->Index()][1] = 0.;
      continue;
    }
    const VesselNode* upstream_node = GetUpstreamNode(v);
    if (!nodal_o2ready[upstream_node->Index()])
    {
      ComputeVesselO2Conc(upstream_node, params, vesselpo2, vl);
      nodal_o2ready[upstream_node->Index()] = true;
    }
    const int side_idx = (upstream_node == v->NodeA()) ? 0 : 1;
    const double po2start = vesselpo2[v->Index()][side_idx];//read out starting point
    double po2end = NANd();
    auto sourceGenerationCallback = [&](int i, int numPoints, double x, double weight, double po2, ComputeRadialFluxes &computeFlux) -> void
    {
      computeFlux.AddSourceContributionsTo(matrix_builder, po2, weight);
    };

    OUTPUT_PO2MEASURECOMP(printf("sim %i, start = %lf\n",v->Index(),po2start);)
    
    NumericallyIntegrateVesselPO2(params, phases, grid, extpo2, vascularPO2PropagationModel.get(), vl, v, upstream_node, po2start, sourceGenerationCallback, po2end, world);
    
    vesselpo2[v->Index()][side_idx  ] = po2start;
    vesselpo2[v->Index()][side_idx^1] = po2end;

    stats.Add(po2start);
    stats.Add(po2end);
    
    OUTPUT_PO2INTEGRATION(fprintf(stderr, "v%i (%i,%i) -> (%f, %f)\n", v->Index(), v->NodeA()->Index(), v->NodeB()->Index(), po2start, po2);)
  }//for each vessel in sorted_vessels list
  
  if (params.loglevel > 0)
    cout << format("IntegrateVesselPO2: %s in %f ms") % stats % (my::Time()-t_).to_ms() << endl;
}



void DetailedP02Sim::PrepareNetworkInfo(const VesselList3d &vl, DynArray<const Vessel*> &sorted_vessels, DynArray<const VesselNode*> &roots)
{
  DynArray<int> order;
  TopoSortVessels(vl, order);
  CheckToposort(vl, order);

  sorted_vessels.resize(vl.GetECount());
  for (int i=0; i<order.size(); ++i)
  {
    myAssert(order[i]>=0);
    sorted_vessels[order[i]] = vl.GetEdge(i);
  }

  roots.clear();
  /* in case some roots are not circulated,
   * we do not want to keep them in the list.
   */
  //uint count_roots = 0;
  for (int i=0; i<vl.GetECount(); ++i)
  {
    const Vessel* v = vl.GetEdge(i);
    if (!v->IsCirculated()) continue;
    
    if (v->NodeA()->IsBoundary())
    {
      //count_roots++;
      roots.push_back(v->NodeA());
      OUTPUT_PO2INTEGRATION(fprintf(stderr,"vessel %i/%i, order %i -> init po2\n",v->Index(), v->NodeA()->Index(), order[v->Index()]);)
    }
    if (v->NodeB()->IsBoundary())
    {
      //count_roots++;
      roots.push_back(v->NodeB());
      OUTPUT_PO2INTEGRATION(fprintf(stderr,"vessel %i/%i, order %i -> init po2\n",v->Index(), v->NodeB()->Index(), order[v->Index()]);)
    }
  }
  //roots.resize(count_roots);
}



template<class T>
struct ConvergenceCriteriumAccumulator2Norm
{
  T val;
  int n;
  ConvergenceCriteriumAccumulator2Norm() : n(0), val() {}
  void Add(const T &x) { val += my::sqr(x); ++n; }
//   T operator()() const { return n>0 ? (std::sqrt(val)/n) : std::numeric_limits<T>::quiet_NaN(); }
  // T.F.
  // Intel compiler is not respecting the quiet_NaN and throw, I do not know why --> intel bug
  T operator()() const { return n>0 ? (std::sqrt(val)/n) : std::numeric_limits<T>::max(); }
};

template<class T>
struct ConvergenceCriteriumAccumulatorMaxNorm
{
  T val;
  int n;
  ConvergenceCriteriumAccumulatorMaxNorm() : val(), n(0) {}
  void Add(const T &x) { val = std::max(val, std::abs(x)); ++n; }
//   T operator()() const { return n>0 ? val : std::numeric_limits<T>::quiet_NaN(); }
  // T.F.
  // Intel compiler is not respecting the quiet_NaN and throw, I do not know why --> intel bug
  T operator()() const { return n>0 ? val : std::numeric_limits<T>::max(); }
};



typedef boost::optional<TissuePhases> OptTissuePhases;

/**
 * @brief Diffusion of oxygen from vessels to tissue
 * 
 * These code propagates the oxygen from the vessel network into the tissue.
 * See section 3 in pdf
 * this is the work horse called by every iteration in the loop
 * 
 * data from previous runs is provided
 */
void ComputePo2Field(const Parameters &params, 
		     const ContinuumGrid &grid, 
		     DomainDecomposition &mtboxes, 
		     const TissuePhases &phases, 
		     Array3df po2field, 
		     FiniteVolumeMatrixBuilder &mb,  
		     bool keep_preconditioner)
{
  if (params.debug_zero_o2field) return;
  
  //I think this initializes the matrix with the values given before
  my::Time t_;
  #pragma omp parallel
  {
    const double consumption_prefactor = 1./params.tissue_solubility;
    BOOST_FOREACH(const DomainDecomposition::ThreadBox bbox, mtboxes.getCurrentThreadRange())
    {
      //diffusion part of differential equation
      mb.AddDiffusion<> (bbox, ConstValueFunctor<float>(1.), -params.po2_kdiff);

      FOR_BBOX3(p, bbox)
      {
        double m, dm, po2 = po2field(p);
        Float3 phases_loc = phases(p);
        //oxygen uptake of tissue according michalis menten model or simpler
        boost::tie(m, dm) = params.ComputeUptake(po2, phases_loc.data(), phases.count);
        double cons_coeff = -dm;
        double cons_const = -(m-dm*std::max(po2, 0.)); // fixed: negative po2 would lead to even more consumption, resulting in more negative po2. So limit evaluation to po2>0.
                
        double lin_coeff = consumption_prefactor*cons_coeff + params.extra_tissue_source_linear;
        double src_const_loc = consumption_prefactor*cons_const + params.extra_tissue_source_const;
        double rhs = -src_const_loc;

        mb.AddLocally(p, -lin_coeff, -rhs);
      }

      if (params.tissue_boundary_condition_flags !=  FiniteVolumeMatrixBuilder::NEUMANN)
        mb.SetDirichletBoundaryConditions(bbox, params.tissue_boundary_condition_flags, params.tissue_boundary_value);
    }
  }
  
  RCP<Epetra_Vector> lhs= rcp(new Epetra_Vector(mb.rhs->Map()));
  #pragma omp parallel
  {
    BOOST_FOREACH(const BBox3 bbox, mtboxes.getCurrentThreadRange())
    {
      FOR_BBOX3(p, bbox)
      {
        (*lhs)[grid.ld.LatticeToSite(p)] = po2field(p);
      }
    }
  }
  //trilionos solver settings
  ptree solver_params = make_ptree  ("preconditioner","multigrid")
                                    ("verbosity", (params.loglevel>0 ? (params.loglevel>1 ? "full" : "normal") : "silent"))("use_smoothed_aggregation", false)
                                    ("max_iter", 500)
                                    ("conv","rhs")
                                    ("max_resid",1.e-8)
                                    ("keep_preconditioner", keep_preconditioner);

  try {
    EllipticEquationSolver &&solver = EllipticEquationSolver{mb.m, mb.rhs, solver_params};
    solver.solve(lhs);
  }
  catch (const ConvergenceFailureException &e)
  {
      if (e.reason == ConvergenceFailureException::MAX_ITERATIONS)
      {
        solver_params.put("keep_preconditioner", false);
        EllipticEquationSolver &&solver = EllipticEquationSolver{mb.m, mb.rhs, solver_params};
        solver.solve(lhs);
      }
      else throw e;
  }
  /** @brief write back the results from matrix solve */
  #pragma omp parallel
  {
    BOOST_FOREACH(const BBox3 bbox, mtboxes.getCurrentThreadRange())
    {
      FOR_BBOX3(p, bbox)
      {
        po2field(p) = (*lhs)[grid.ld.LatticeToSite(p)];
      }
    }
  }

  if (params.loglevel > 0)
    cout << format("Po2Field: %s in  %f ms") % po2field.valueStatistics() % (my::Time()-t_).to_ms() << endl;
}


void DetailedP02Sim::init(Parameters &params_, 
                          BloodFlowParameters &bfparams_,
                          VesselList3d &vl, 
                          double grid_lattice_const, 
                          double safety_layer_size, 
                          boost::optional<Int3> grid_lattice_size,
                          boost::optional<H5::Group> tumorgroup,
                          boost::optional<Array3df> previous_po2field,
                          boost::optional<DetailedPO2::VesselPO2Storage> previous_po2vessels
)
{
  bfparams = bfparams_;
  
  CalcFlow(vl, bfparams);
  params = params_;
  world = !vl.HasLattice();
  
  //multithreading
  my::SetNumThreads(params.num_threads);
  {
    //this worked only for lattices
    //int dim = (::Size(vl->Ld().Box())[2]<=1) ? 2 : 3;
    int dim=0;
    if (world)
    {
      dim = 3;
    }
    else
    {
      dim = (::Size(vl.Ld().Box())[2]<=1) ? 2 : 3;
    }
    
    if (grid_lattice_size)
    {
      ld.Init(*grid_lattice_size, grid_lattice_const);
      ld.SetCellCentering(Bool3(1, 1, dim>2));
      wbox = vl.Ld().GetWorldBox();
      worldCenter = 0.5*(wbox.max + wbox.min);
      gridCenter = 0.5*(ld.GetWorldBox().max + ld.GetWorldBox().min);
      ld.SetOriginPosition(ld.GetOriginPosition() + (worldCenter - gridCenter));
    }
    else
    {
      //added safety space to reduce boundary errors
      if (world)
      {
        SetupFieldLattice(vl.GetWorldBoxFromVesselsOnly(), dim, grid_lattice_const, safety_layer_size, ld);
      }
      else
      {
        SetupFieldLattice(vl.Ld().GetWorldBox(), dim, grid_lattice_const, safety_layer_size, ld);
      }
    }
    grid.init(ld, dim);
    mtboxes.init(MakeMtBoxGrid(grid.Box(), Int3(32, 32, 32)));
    if (params.loglevel > 0)
    {
      cout << "continuum grid:" << endl;
      grid.ld.print(cout);
      cout << endl;
      cout << "multithreading:" << endl;
      cout << my::GetNumThreads()<<endl;
    }
    if (params.loglevel > 0)
    {
      cout << "vessel lattice" << endl;
      if (world)
      {
        cout<<vl.GetWorldBoxFromVesselsOnly()<<endl;
      }
      else
      {
        vl.Ld().print(cout);
      }
      cout << endl;
    }
    isVesselListGood(vl);
    if(tumorgroup)
    {
      SetupTissuePhases(phases, grid, mtboxes, tumorgroup);//filling
    }
    
    //set up field with same discrete points as in the given grid, leave data memory uninitialized
    //po2field is declared by mother function by not initialized, that happening here
    po2field = Array3df(grid.Box(), Cons::DONT);
    po2vessels.resize(vl.GetECount(), Float2(NANf()));
    // I try this hardcoded before do it with variables: read in po2field from previous iteration 
    if(previous_po2field and previous_po2vessels)
    {
      cout<<"works"<<endl;
      po2field.fill(*previous_po2field);
      if(vl.GetECount() <= previous_po2vessels->size())
      {
        for(int i = 0;i<vl.GetECount();++i)
        {
          po2vessels[i] = previous_po2vessels->operator[](i);
        }
      }
      else
      {
	printf("Warning: more vessels suggested than in previous run present!\n");
      }
      //previous_po2vessels
//       h5cpp::File f2("/localdisk/thierry/output_milotti/with_o2_sim/fakeTumMTS-default-typeI-sample00-milotti_detailed.h5", "r");
//       h5cpp::Group po2Group = f2.root().open_group("out0000/po2");
//       h5cpp::Dataset ds;
//       
//       ds = po2Group.open_dataset("po2field");
// 
//       // check reading the dataspace dimensions first
//       h5cpp::Dataspace sp = ds.get_dataspace();
//       hsize_t dims[H5S_MAX_RANK];
//       int rank = sp.get_dims(dims);
//       cout << "dataset rank = " << rank << ", dims = ";
//       for (int i = 0; i < rank; ++i)
//       {
//         cout << dims[i] << " ";
//       }
//       cout << endl;
// 
//       { // check reading to memory block pointed to by a c-style pointer
//         vector<float> data(sp.get_npoints());
//         ds.read<float>(&data[0]);
// 
//         cout << "dataset = ";
//         for (auto it = data.begin(); it != data.end(); ++it)
//         {
//           cout << *it << ", ";
//         }
//         cout << endl;
//       }
//       po2field.fill(ds);
//     { // check reading to vector
//       ds = root.open_dataset("dataset_from_list");
//       vector<int> data; h5::read_dataset<int>(ds, data);
//     }
      
      
    }
    else
    {
      po2field.fill(params.debug_zero_o2field ? 0.f : params.po2init_cutoff);
    }
    
    
    
    
    /** executes topological ordering
     * potential remodelling due to tumor should be considered
     */
    sorted_vessels.clear();
    roots.clear();
    arterial_roots.clear();
    PrepareNetworkInfo(vl, sorted_vessels, roots);
    BOOST_FOREACH(const VesselNode* nd, roots)
    {
      if(nd->Count()>0)// there needs to be a connected vessel to check for type
      {
	if(nd->GetEdge(0)->IsArtery())
	{
	  arterial_roots.push_back(nd);
	}
      }
    }   //note vl is lost after this call!!!!
  }
}

/**
 * @brief Head function called by python interface
 */
int DetailedP02Sim::run(VesselList3d &vl)
{
  //sets up the linear trilionos matrix system, builder is implemented as struct
  //could use for example different stencils
  FiniteVolumeMatrixBuilder tissue_diff_matrix_builder;

#if APPROXIMATE_FEM_TRANSVASCULAR_EXCHANGE_TERMS
  tissue_diff_matrix_builder.Init7Point(grid.ld, grid.dim);
#else
  tissue_diff_matrix_builder.Init27Point(grid.ld, grid.dim);
#endif
  //buffer within consecutive runs
  Array3df last_po2field(grid.Box());
  //maybe 0 or prams.po2init_cutoff
  last_po2field.fill(po2field);
  // an other buffer
  DynArray<Float2> last_vessel_po2(vl.GetECount(), Float2(0.));//begining and end 0.
  //metadata.add_child("iterations", ptree());
  /* when doing iterative calls, we need to override this */
  metadata.put_child("iterations", ptree());
  
  ConvergenceCriteriumAccumulator2Norm<double> delta_field2, delta_vess2;
  ConvergenceCriteriumAccumulatorMaxNorm<double> delta_fieldM, delta_vessM;
  /*
   * ****** MAIN LOOP **********
   */
  for (int iteration_num = 0;; ++iteration_num)
  {
    if (!params.debug_fn.empty() && ((iteration_num % 1) == 0) && iteration_num>0)
    {
      //h5cpp::File f(params.debug_fn, iteration_num==0 ? "w" : "a");
      H5::H5File f(params.debug_fn,iteration_num==0 ? H5F_ACC_RDWR : H5F_ACC_TRUNC);
      WriteOutput(f.createGroup(str(format("out%04i-a") % iteration_num)),
                  vl, 
		  params,
                  //vesselpo2,
                  po2vessels,
                  sorted_vessels,
                  grid, 
		  po2field,
                  tissue_diff_matrix_builder);
      f.close();
    }
    const double tolerance = params.convergence_tolerance;
    // break if convergent
    const double res_delta_fieldM = delta_fieldM();
    const double res_delta_vessM = delta_vessM();
    if (iteration_num > params.max_iter || (res_delta_fieldM<tolerance && res_delta_vessM<tolerance))
      break;
    
    //initialize output with zero
    tissue_diff_matrix_builder.ZeroOut();
    /*
     * 1) propagate the oxygen along the blood stream
     */
    //IntegrateVesselPO2(params, po2vessels, vl, sorted_vessels, roots, grid.ld, po2field, phases, tissue_diff_matrix_builder,world);
    IntegrateVesselPO2(params, po2vessels, vl, sorted_vessels, arterial_roots, grid.ld, po2field, phases, tissue_diff_matrix_builder,world);
    if (!params.debug_fn.empty() && ((iteration_num % 1) == 0) && iteration_num>0)
    {
      //h5cpp::File f(params.debug_fn, iteration_num==0 ? "w" : "a");
      H5::H5File f(params.debug_fn,iteration_num==0 ? H5F_ACC_RDWR : H5F_ACC_TRUNC);
      WriteOutput(f.createGroup(str(format("out%04i-b") % iteration_num)),
                  vl, 
		  params,
                  po2vessels,
                  sorted_vessels,
                  grid, 
		  po2field,
                  tissue_diff_matrix_builder);
      f.close();
    }
    //bool keep_preconditioner = (iteration_num>2 && tissue_diff_solver.iteration_count<25);
    bool keep_preconditioner = true;
    /*
     * 2) propagate the oxygen from the blood stream to the tissue
     */
    ComputePo2Field(params, grid, mtboxes, phases, po2field, tissue_diff_matrix_builder, keep_preconditioner);
    
    /*
     * From here on the results are handled
     */
    {//interupt
      // dampening: new values are a linear combination of previous and currently computed values. 
      // f is the fractional share of the previous value.
      const double f = 0.3; 
      delta_fieldM = delta_vessM = ConvergenceCriteriumAccumulatorMaxNorm<double>();
      delta_field2 = delta_vess2 = ConvergenceCriteriumAccumulator2Norm<double>();
      //save the changes to prior run
      for (int i=0; i<po2vessels.size(); ++i)
      {
        delta_vessM.Add(po2vessels[i][0]-last_vessel_po2[i][0]);
        delta_vessM.Add(po2vessels[i][1]-last_vessel_po2[i][1]);
        delta_vess2.Add(po2vessels[i][0]-last_vessel_po2[i][0]);
        delta_vess2.Add(po2vessels[i][1]-last_vessel_po2[i][1]);
        po2vessels[i] = po2vessels[i]*(1.-f)+f*last_vessel_po2[i];
        last_vessel_po2[i] = po2vessels[i];
      }
      //loop over all points in the ContinuumGrid
      FOR_BBOX3(p, grid.Box())
      {
        float &current = po2field(p);
        float &last    = last_po2field(p);
        // michaelis menten solution with large zero-order term can undershoot the po2 field in negative values
        // which the vessel po2 inegration routine cannot stand. Therefor the po2 field is limited here from below.
        // Realistic solutions should have positive values without cutoff ofc.
        current = std::max(0.f, current);
        delta_fieldM.Add(last-current);
        delta_field2.Add(last-current);
        current = current*(1.-f) + last*f;
        last = current;
      }
      if (params.loglevel > 0)
        cout << format("iteration %i: dvM=%f, dfM=%f, dv2=%f, df2=%f") %  iteration_num % delta_vessM() % delta_fieldM() % delta_vess2() % delta_field2() << endl;
      else
        cout << ".";
      {
        ptree node;
        node.put("iteration", iteration_num);
        node.put("delta_vessM", delta_vessM());
        node.put("delta_fieldM", delta_fieldM());
        node.put("delta_vess2", delta_vess2());
        node.put("delta_field2", delta_field2());
        node.put("dampening", f);
        metadata.get_child("iterations").add_child(str(format("iteration%04i") % iteration_num), node);
      }
    }//end interupt
    last_po2field.initCopy(po2field);
    //last_vessel_po2 = vesselpo2;
    last_vessel_po2 = po2vessels;
    if (my::checkAbort())
      return 1;
  }//end loop nointerations
  if (params.loglevel == 0)
    cout << endl;
  else
    cout << "computing final results" << endl;
  
  tissue_diff_matrix_builder.ZeroOut();
  IntegrateVesselPO2(params, po2vessels, vl, sorted_vessels, arterial_roots, grid.ld, po2field, phases, tissue_diff_matrix_builder, world);
  //IntegrateVesselPO2(params, po2vessels, vl, sorted_vessels, roots, grid.ld, po2field, phases, tissue_diff_matrix_builder, world);
  
  ComputePo2Field(params, grid, mtboxes, phases, po2field, tissue_diff_matrix_builder, true);
  return 0;
}



#if 1

Measurement::Measurement(std::auto_ptr< VesselList3d > vl_ptr_, const DetailedPO2::Parameters& params_, const ContinuumGrid& grid_, const TissuePhases &phases_, Array3df& po2field_, DetailedPO2::VesselPO2Storage& vesselpo2_) :
  vl(vl_ptr_), po2field(po2field_), vesselpo2(vesselpo2_), params(params_), phases(phases_), grid(grid_)
{
  o2mass_in_root = o2mass_out_root = o2mass_out_vessels = o2mass_consumed_transvascular = 0.;
  vascularPO2PropagationModel = std::auto_ptr<VascularPO2PropagationModel>(VascularPO2PropagationModel::allocate(params).release());
  
  for (int i=0; i<vl->GetECount(); ++i)
  {
    const Vessel* v = vl->GetEdge(i);
    const VesselNode *nda = v->NodeA(), *ndb = v->NodeB();
    if (!((nda->flags & BOUNDARY) || (ndb->flags & BOUNDARY))) continue;
    const VesselNode* upstream_node = GetUpstreamNode(v);
    if (nda->flags & BOUNDARY)
    {
      double flx = v->q * params.BloodPO2ToConc(vesselpo2[i][0], v->hematocrit);
      if (nda == upstream_node) o2mass_in_root += flx;
      else o2mass_out_root += flx;
    }
    if (ndb->flags & BOUNDARY)
    {
      double flx = v->q * params.BloodPO2ToConc(vesselpo2[i][1], v->hematocrit);
      if (ndb == upstream_node) o2mass_in_root += flx;
      else o2mass_out_root += flx;
    }
  }
}
/** @brief
 * this calculated the quantities
 * \param x obviously x is the length of vessel in units of segments
 * \param po2 oxygen partial pressure along vessel
 * \param extpo2 external tissue po2, pressure extern of vessel
 * \param jtv transvascular flux
 * \param dS_dx change of saturation along vessel?
 * see sampleVessels in __init__.py of detailedo2
 */
void Measurement::computeVesselSolution(int idx, DynArray< VesselPO2SolutionRecord >& sol, bool world)
{
    sol.remove_all();
    const Vessel* v = vl->GetEdge(idx);
    
    if (!v->IsCirculated()) return;

    sol.reserve(128);
    
    const VesselNode* upstream_node = GetUpstreamNode(v);
    const int side_idx = (upstream_node == v->NodeA()) ? 0 : 1;
    const double po2start = vesselpo2[v->Index()][side_idx];
    const double po2end   = vesselpo2[v->Index()][side_idx^1];
    const bool reverse = v->NodeA() != upstream_node;
    
    double wl;
    if( world )
      wl = (v->NodeA()->worldpos-v->NodeB()->worldpos).norm();
    else
      wl = v->WorldLength(vl->Ld());
    myAssert(wl>0);
    const double len = wl;
    
    auto sourceGenerationCallback = [&](int i, int numPoints, double x, double weight, double po2, ComputeRadialFluxes &computeFlux) -> void
    {
      double jinner, jouter;
      tie(jinner, jouter) = computeFlux.ComputeFluxes(po2);
      VesselPO2SolutionRecord r;
      r[0] = reverse ? len-x : x;
      r[1] = po2;
      r[2] = computeFlux.po2tissue;
      r[3] = jinner / (my::mconst::pi2() * v->r); // FIX: division by circumference because this should be j_tv which is defined as flux per surface area. With this fix i get quite correct transvascular flux & gradient measurements!
      r[4] = std::abs(params.DiffSaturation(po2).second*ComputePO2RateOfChange(params, po2, v->q, v->hematocrit, jinner));
      sol.push_back(r);
      this->o2mass_out_vessels += jinner*weight;
      this->o2mass_consumed_transvascular += (jinner - jouter)*weight;
    };

    double newpo2end = 0;
    NumericallyIntegrateVesselPO2<decltype(sourceGenerationCallback)>(params, phases, grid, po2field, vascularPO2PropagationModel.get(), *vl, v, upstream_node, po2start, sourceGenerationCallback, newpo2end,world);
    if (std::abs(newpo2end - po2end) > 0.1)
    {
      cout << format("warning po2 difference at ends of %i for %f (loaded) and %f (measure) is way too large") % v->Index() % po2end % newpo2end << endl;
    }

    
    OUTPUT_PO2MEASURECOMP(printf("measure %i, start = %lf, end = %lf\n",v->Index(),po2start,po2end);)
    
    if (reverse)
      std::reverse(sol.begin(), sol.end());
}
#endif



/*--------------------------------------------------------------------------
 * debug output
---------------------------------------------------------------------------- */

void WriteOutput(H5::Group basegroup,
                 const VesselList3d &vl,
                 const Parameters &params,
                 const boost::optional<const VesselPO2Storage&> vesselpo2,
                 const boost::optional<DynArray<const Vessel*>&> sorted_vessels,
                 const boost::optional<ContinuumGrid&> grid,
                 const boost::optional<Array3df> po2field,
                 const boost::optional<const FiniteVolumeMatrixBuilder&> mbopt)
{
//     h5cpp::File f(fn,"w");
//     h5cpp::Group basegroup = grpname.empty() ? f.root() : f.root().create_group(grpname);
    H5::Group g = basegroup.createGroup("vessels");
    WriteVesselList3d(vl, g, make_ptree("w_all",false)("w_pressure",true));

    if (sorted_vessels)
    {
      DynArray<int> toposort_indices(vl.GetECount());
      for (int i=0; i<vl.GetECount(); ++i)
      {
        const Vessel* v = (*sorted_vessels)[i];
        toposort_indices[v->Index()] = i;
      }
      //h5cpp::create_dataset(g.open_group("edges"), "topoorder", toposort_indices);
      //h5cpp::create_dataset(g.open_group("edges"), "topoorder", toposort_indices);
      writeDataSetToGroup(g.openGroup("edges"), string("topoorder"), toposort_indices);
    }
    if (vesselpo2)
    {
      DynArray<float> avg_po2(vl.GetECount());
      for (int i=0; i<vl.GetECount(); ++i)
      {
        float po2 = ((*vesselpo2)[i][0]+(*vesselpo2)[i][1])*0.5;
        avg_po2[i] = isFinite(po2) ? po2 : -1.f;
      }
      //h5cpp::create_dataset(g.open_group("edges"), "avgpo2", avg_po2);
      writeDataSetToGroup(g.openGroup("edges"), string("avgpo2"), avg_po2);
    }

    if (po2field)
    {
      //h5cpp::Group ld_group = RequireLatticeDataGroup(basegroup, "field_ld", grid->ld);
      H5::Group ld_group = RequireLatticeDataGroup(basegroup, grid->ld);
      g = basegroup.createGroup("fields");
      if (po2field)
        WriteScalarField(g, "po2field", *po2field, grid->ld, ld_group);
    }

    if (mbopt)
    {
      H5::Group ld_group = RequireLatticeDataGroup(basegroup, grid->ld);
//       g = basegroup.require_group("fields");
      g = basegroup.openGroup("fields");
      
      const FiniteVolumeMatrixBuilder& mb = *mbopt;
      const Epetra_CrsMatrix& mat = *mb.m;
      const Epetra_Vector& rhs = *mb.rhs;
      const BBox3 bbox = grid->Box();
      Array3d<double> diagArr(bbox),
                      rhsArr(bbox),
                      rowsumArr(bbox);
      int numEntries;
      double values[128];
      int    indices[128];
      FOR_BBOX3(p, bbox)
      {
        int row = grid->ld.LatticeToSite(p);
        double rowsum = 0;
        mat.ExtractGlobalRowCopy(row, 27, numEntries, values, indices);
        for (int i=0; i<numEntries; ++i)
        {
          if (indices[i] == row) diagArr(p) = values[i];
          rowsum += values[i];
        }
        rhsArr(p) = rhs[row];
        rowsumArr(p) = rowsum;
      }

      WriteScalarField(g, "diag", diagArr, grid->ld, ld_group);
      WriteScalarField(g, "rowsum", rowsumArr, grid->ld, ld_group);
      WriteScalarField(g, "rhs", rhsArr, grid->ld, ld_group);
    }
}




/*--------------------------------------------------------------------------
 * test functions go here
---------------------------------------------------------------------------- */

void TestSaturationCurve()
{
  H5::H5File f("detailedpo2_saturationtest.h5", H5F_ACC_RDWR);
  H5::Group root = f.openGroup("/");
  int n = 100000;
  double h = 0.45;
  Parameters params;
  DynArray<double> xy(n*2);

  {
  my::Time t_;
  for( int i=0; i<n; ++i )
  {
    double p = (params.sat_curve_p50*5*i)/n;
    double conc = params.BloodPO2ToConc(p, h);
    xy[i*2+0] = p;
    xy[i*2+1] = conc;
  }
  double t_ms = (my::Time()-t_).to_ms();

//   h5cpp::Dataset ds = h5cpp::create_dataset<double>(root, "p_to_conc", h5cpp::Dataspace::simple_dims(n, 2), get_ptr(xy));
//   h5cpp::Dataset ds = h5cpp::create_dataset<double>(root, "p_to_conc", h5cpp::Dataspace::simple_dims(n, 2), get_ptr(xy));
  H5::DataSet ds = writeDataSetToGroup(root, string("p_to_conc"), xy );
  writeAttrToH5(ds, string("time"), t_ms);
  
//   ds.attrs().set("time", t_ms);
  }

  {
  my::Time t_ = my::Time();
  double maxconc = params.BloodPO2ToConc(params.sat_curve_p50*5., h);
  for( int i=0; i<n; ++i )
  {
    double conc = (maxconc*i)/n;
    double p = params.ConcToBloodPO2(conc, h);
    xy[i*2+0] = p;
    xy[i*2+1] = conc;
  }
  double t_ms = (my::Time()-t_).to_ms();

  //h5cpp::Dataset ds = h5cpp::create_dataset(f.root(), "conc_to_p", h5cpp::Dataspace::simple_dims(n, 2), get_ptr(xy));
  //ds.attrs().set("time", t_ms);
  H5::DataSet ds = writeDataSetToGroup(root, string("conc_to_p"), xy );
  writeAttrToH5(ds, string("time"), t_ms);
  }

  {
  for( int i=0; i<n; ++i )
  {
    double p = (params.sat_curve_p50*5*i)/n;
    double conc = params.BloodPO2ToConc(p, h);
    double pback = params.ConcToBloodPO2(conc, h);
    xy[i*2+0] = p;
    xy[i*2+1] = pback;
  }
  //h5cpp::Dataset ds = h5cpp::create_dataset(f.root(), "diff", h5cpp::Dataspace::simple_dims(n, 2), get_ptr(xy));
  H5::DataSet ds = writeDataSetToGroup(root, string("diff"), xy );
  }

  {
  xy.resize(n*3);
  for( int i=0; i<n; ++i )
  {
    double p = (params.sat_curve_p50*5*i)/n;
    double s, ds; tie(s,ds) = params.DiffSaturation(p);
    xy[i*3+0] = p;
    xy[i*3+1] = s;
    xy[i*3+2] = ds;
  }
  //h5cpp::Dataset ds = h5cpp::create_dataset(f.root(), "ds", h5cpp::Dataspace::simple_dims(n, 3), get_ptr(xy));
  H5::DataSet ds = writeDataSetToGroup(root, string("ds"), xy );
  }
  
  f.close();
}


void TestSingleVesselPO2Integration()
{
#if 0
  double len = 100.;
  double h = 0.45;
  double q_scale = 1.e-18;
  double flow_rate = q_scale*(3.14*4.*4.*1000.); // 1mm/s velocity, through 4um radius capillary
  // permeability: 0.18 nlO2 / (cm^2 sec mmHg)  (Kiani 2003)
  // = 0.18 * 1e-9 * 1e3 * 1e12 um^3 O2 / (1e8 um^2 sec mmHg)
  // = 0.18 * 1e-2 um^3 O2 / (um^2 sec mmHg)
  double permeability = 0.18e-2;
  Parameters params;

  auto extpo2callback = [=](double x)
  {
    return 0.;
  };
  ModelPO2Integration<decltype(extpo2callback)> model(params, flow_rate, h, permeability, extpo2callback, 30.);

  std::vector<double> xyz;
  const int N = 10;
  double dx = len/(N-1);
  double x = 0., po2 = params.PInit(100.), euler_dt = 0.;
  for (int i=0; i<N; ++i)
  {
    xyz.push_back(x);
    xyz.push_back(po2);
    xyz.push_back(euler_dt);
    model.doStep(po2, x, dx);
    euler_dt = model.getEulerDt();
  }

  h5cpp::File f("integratepo2.h5","w");
  h5cpp::Group root = f.root();
  h5cpp::create_dataset(root, "data", h5cpp::Dataspace::simple_dims(N, 3), get_ptr(xyz));
#endif
}
Array3df DetailedP02Sim::getPo2field()
{
  return po2field;
}
DetailedPO2::VesselPO2Storage DetailedP02Sim::getVesselPO2Storrage()
{
  return po2vessels;
}

}//namespace DetailedPO2
