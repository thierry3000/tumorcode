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
#ifndef IFDRUG_SIM
#define IFDRUG_SIM

#include "calc-ifflow.h"
#include "convection_diffusion_solver.h"
#include <hdfio.h>

#define USE_IFDRUGSIM
#ifdef USE_IFDRUGSIM
namespace IfDrug
{


struct Params
{
  Params();
  void WriteH5(H5::Group g) const;
  void assign(const ptree &pt);
  ptree as_ptree() const;

  FlReal inject_t,
         inject_max,
         kdiff,
         capillary_permeability_normal,
         capillary_permeability_tumor,
         max_uptake,
         chalf_uptake,
         linear_uptake_coeff,
         comprates_k12,
         comprates_k21;
  int inject_mode, uptake_mode;
  string stepper;
  string stepper_compartments;
  bool convective_transport_enabled;

  enum DfInjectModes
  {
    DF_INJECT_MODE_EXP = 0,
    DF_INJECT_MODE_JUMP = 1,
    DF_INJECT_MODE_JUMP_UP = 2,
  };
  enum DfUptakeModes
  {
    DF_UPTAKE_LINEAR = 0,
    DF_UPTAKE_MICHAELIS_MENTEN = 1
  };
};



class VesselDrugExavasationModel : public VesselExavasationModel
{
    const Params *params;
  public:
    void Init(const VesselList3d &vl_, const ContinuumGrid &grid_, const VesselsInBoxes &vesselsinboxes_, const Params &param_);
    bool GetVesselInfo(const Vessel* v, double &coeffa, double &coeffb, double &vala, double &valb) const;
};


#if 0
double getWallPermeability(double r, double m, const Params *params);

class ContinuumVesselExavasationModel : public VirtualGridFunctions<double, 2>
{
    boost::function<double (double, double)> get_permeability;
    boost::shared_ptr<VirtualGridFunctions<float, 4> > get_vessel; // radii, maturation, S/V, intravascular value
    
  public:
    void GetValues(int boxindex, const BBox3 &bbox, VirtualGridFunctions<double, 2>::List &values) const
    {
      Array3d<float> radii(bbox), maturation(bbox), s_per_v, iv_value(bbox);
      VirtualGridFunctions<float, 4>::List arr = {{ radii, maturation, s_per_v, iv_value }};
      get_vessel(boxindex, bbox, arr);
            
      FOR_BBOX3(p, bbox)
      {
        double w = get_permeability(radii(p), maturation(p));
        values[IFF_LINEAR_COEFF](p) = -w * s_per_v(p);
        values[IFF_CONST_COEFF](p) = w * s_per_v(p) * iv_value(p);
      }
    }
};
#endif

/**
 * @brief drug stuff
 * Here comes the information about the drug.
 * At the moment we consider two fields, the concentration of drugs in the 
 * extracellular space and inside the cell.
 * 
 * ID_CONC_EXTRACELLULAR
 * ID_CONC_CELL
 */
enum {
  ID_CONC_EXTRACELLULAR = 0,
  ID_CONC_CELL = 1,
  NUM_FIELDS
};
struct State : boost::noncopyable
{
  int id;
  Array3d<float> field[2];
  int num_fields;

  State() : id(-1), num_fields(0) {}

  void addScaled(double s, const State &u, double s_self = 1.)
  {
    for (int i=0; i<num_fields; ++i)
      AddScaled(s_self, field[i], s, u.field[i]);
  }

  void setId(int id_) { id = id_; }
  int getId() const { return id; }

  void init(int num_fields_, const LatticeDataQuad3d &grid, int ndims, int border)
  {
    num_fields = num_fields_;
    for (int i=0; i<num_fields; ++i)
      field[i] = MakeArray3dWithBorder<float>(grid.Box(), ndims, border);
    id = 0;
  }

  void initCloned(const State &other)
  {
    num_fields = other.num_fields;
    id = other.id;
    for (int i=0; i<num_fields; ++i)
      field[i].initDeepCopy<>(other.field[i]);
  }

  void initLike(const State &other)
  {
    num_fields = other.num_fields;
    id = 0;
    for (int i=0; i<num_fields; ++i)
      field[i].initFromBox(other.field[i].getBox());
  }

  void clear()
  {
    num_fields = 0;
    for (int i=0; i<num_fields; ++i)
      field[i].clear();
  }
};

/**
 * @brief calculte the drug transport.
 * executed by the assigned callback.
 */
class Calculator
{
  public:
    typedef float T;
    typedef IfDrug::State State; ///State is for the drug
    typedef Array3d<T> State1; ///State1 is for the IFF
    typedef boost::function<void (int, const BBox3 &, Array3d<float> &, Array3d<float> &, Array3d<float> &)> Callback1;
    /// arguments are boxindex, box, phi_cells, phi_water, phi_ecm
    typedef VirtualGridFunctions<double, 2> VesselModel;

  public:
    Calculator() : grid(NULL), iffstate(NULL), time(-1), params(NULL) {}
    //void Init(const ContinuumGrid &grid_, const IfFlowState &iffstate_, const VirtualGridFunctions<float,1> &cell_fraction_f_, const VirtualGridFunctions<float,1> &water_fraction_f_, const Params &params_);
    void Init(const ContinuumGrid &grid_, const DynArray<BBox3> &mtboxes_, const IfFlowState &iffstate_, Callback1 tissueCompositionCallback_, VesselModel &vesselModel_, const Params &params_);
    //void InitVessels(const VesselList3d& vl, const VesselsInBoxes& vesselsinboxes_);
    
    Steppers::StepControl doStep(State &state, const Steppers::StepControl &ctrl);

    void writeH5(H5::H5File f, H5::Group g, const State &state, double t, H5::Group ld_group) const;

    enum {
      FLAG_UPTAKE = 1<<0,
      FLAG_SRC_DIFF_VESS = 1<<1,
      FLAG_SRC_FLOW_VESS = 1<<2,
      FLAG_SINK_FLOW_VESS = 1<<3,
      FLAG_SINK_FLOW_LYMPH = 1<<4,
      FLAG_ALL = ~0
    };

    T GetInflowVesselConc() const;  /// time dependent vessel drug conc.
    
  private:
    friend struct ThreadLocal;
    const Params *params;
    Array3d<double> vessel_source_coefficient_buffer;
    const ContinuumGrid *grid;
    const DynArray<BBox3> *mtboxes;
    const IfFlowState *iffstate;
    Callback1 tissueCompositionCallback;
    
    double time, euler_dt_interact, max_rate_interact;
    double max_kdiff, max_src_expl, max_src_impl, max_vel, euler_dt_kdiff, euler_dt_srcexpl, euler_dt_srcimpl, euler_dt_vel;
    
    /// interaction with intra-cellular compartment; callback functions for time steppers
    void calcSlopeEc(const State &state,
                  State &slope,
                  Steppers::StepControl &ctrl);
    void calcSlopesEc(const State &state,
                        State &slope_expl,
                        State &slope_impl,
                        Steppers::StepControl &ctrl);
    /// for the IMEX method we need two slopes
    void invertImplicitOperatorEc(State &lhs,
                                const State &rhs,
                                double identity_coeff,
                                double operator_coeff,
                                Steppers::StepControl &ctrl,
                                State &extrapolation);
    // called by calcSlopesEC and invertImplicitOperatorEC. Computes rate coefficient matrix and hands the actual work of computing the slope to the callback func.
    void calcSlopesEcHelper(const State &state, boost::function<void (int, const BBox3 &, const Array3d<Eigen::Matrix<float, NUM_FIELDS, NUM_FIELDS> > )> func);


    // interstitial compartment; callback functions for time steppers
    void calcSlopeIf(const State1 &state, State1 &slope, Steppers::StepControl &ctrl);
    void calcSlopesIf(const State1 &state, State1 &slope_expl, State1 &slope_impl, Steppers::StepControl &ctrl);
#if 1
    void insertLinearSystemCoefficientsIf(int boxindex,
                                          const BBox3 &bbox,
                                          const State1 &rhs,
                                          const State1& extrapolation,
                                          double identity_coeff,
                                          double operator_coeff,
                                          FiniteVolumeMatrixBuilder &mb);
#endif
#if 1
    void invertImplicitOperatorIf(State1 &lhs,
                                const State1 &rhs,
                                double identity_coeff,
                                double operator_coeff,
                                Steppers::StepControl &ctrl,
                                State1 &extrapolation);
#endif
    /**
     * these Callbacks are present for the Steppers
     * preStep = preStepDefault;
     * postStep = postStepDefault;
     * calcSlope = calcSlopeDefault;
     * calcSlopesIMEX = calcSlopesIMEXDefault;
     * invertImplicitOperator = invertImplicitOperatorDefault;
     */ 
    typedef Steppers::Stepper<Steppers::Callback<State>, Steppers::STORAGE_COPY>  StepperType2;
    std::auto_ptr<StepperType2> stepper_interaction;//interaction with drug!
    typedef Steppers::Stepper<Steppers::Callback<State1>, Steppers::STORAGE_COPY> StepperType1;
    std::auto_ptr<StepperType1> stepper_interstitial;// interstitial fluid flow

    // for operator splitting
    // num=0 computes interstitial compartment,
    // num=1 computes interstitial-cellular interaction
    Steppers::StepControl doSplitStep(int num, State& state, const Steppers::StepControl& ctrl);
};

}
#endif


#endif
