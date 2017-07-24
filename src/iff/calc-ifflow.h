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
#ifndef CALC_IFFLOW_H
#define CALC_IFFLOW_H

//#define USE_IFFLOW
//#ifdef USE_IFFLOW

#include "mwlib/field.h"
#include "continuum-flow.h"
#include "calcflow.h"
#include "vessels3d.h"
#include "hdf_wrapper.h"
#include "mwlib/time_stepper.h"
#include "convection_diffusion_solver.h"
#include "shared-objects.h"

//#include <boost/array.hpp>

//#define USE_IFDRUGSIM

/**
  units:
    time - sec
    space - micron
    pressure - kPa
*/


#if 0
struct IFNodeInfo
{
  FlReal press;
  bool bBoundary;
  Float3 wpos;
};

struct IFVesselSeg
{
  Float3 wpos;
  float weight_per_volume;
  float axial_pos;
};
#endif


struct IffParams
{
  IffParams();
  void WriteH5(h5cpp::Group g) const;
  void assign(const ptree &pt);
  ptree as_ptree() const;
  
  // interstitial fluid flow parameters
  bool debugOut, coupled_flow, better_permeability_model;
  double iff_cond_tissue,
         iff_cond_tumor,
         iff_cond_necro,
         ecm_friction,
         cell_friction,
         iff_capillary_permeability_tissue,
         iff_capillary_permeability_tumor,
         lymph_surf_per_vol,
         lymph_surf_per_vol_tumor,
         lymph_press,
         lymph_permeability,
         capillary_oncotic_pressure,
         interstitial_oncotic_pressure,
         osmotic_reflection_coeff,
         osmotic_reflection_coeff_tumor;
};

/**
 * info about: 	pfield (IFP field)
 * 		vessel_inflow_field
 * 		vesesl_outflow_field
 * 		lymph_outflow_field
 * 		source_field
 */

struct IfFlowState : boost::noncopyable
{
  enum {
    PRESSURE = 1,
    VELOCITIES,
    SOURCE,
    SOURCE_DETAIL,
    ALL = ~0
  };

  // continuum fields
  Array3d<double> pfield;       // fluid pressure, border = 1
  Array3d<double> vessel_inflow_field, vessel_outflow_field, lymph_outflow_field;
  Array3d<double> source_field;
  
  boost::array<Array3d<float>, 3> velocities;
  double max_velocity;
  int flags;

  IfFlowState() : flags(0), max_velocity(0) {}
  void Init(const BBox3& bbox, int dim, int flags);
};


/*------------------------------------------------------
------------------------------------------------------*/


class VesselExavasationModel : boost::noncopyable, public VirtualGridFunctions<double, 2>
{
  /** @brief
   * This should be a concrete implementation of a abstract base class.
   * A similar class could be used for oxygen exavasation
   */
  const ContinuumGrid *grid;

  typedef VesselList3d::LatticeData VLD;
  const VesselList3d *vl;
  const VLD *vess_ld;
  const VesselsInBoxes *vesselsinboxes;
  //const IffParams *params;
  
public:
  VesselExavasationModel() : grid(NULL), vl(NULL), vess_ld(NULL) {}
  virtual ~VesselExavasationModel() {}

  void Init(const VesselList3d &vl_, const ContinuumGrid &grid_, const VesselsInBoxes &vesselsinboxes_);

  /** @brief
   * Get info about wall permeability, etc, The source term is coeff * (value - outside_value)
   * This term is evaluated for some sampling points, where *a and *b are linearly interpolated along the
   * centeral axis of the vessel
   */
  virtual bool GetVesselInfo(const Vessel* v, double &coeffa, double &coeffb, double &vala, double &valb) const = 0;
  
  /** @brief
   * computes the factors in  source_i = sum_j c_ij (p_j - pe_i),
   * where clinear_i = - sum_j c_ij, and cconst_i = sum_j c_ij p_j
   * The j-summation is over samples on the vessel tree, and the i-index stands for a grid cell
   * clinear = values[0]
   * cconst = values[1]
   */
  virtual void GetValues(int boxindex, const BBox3 &bbox, VirtualGridFunctions<double, 2>::List &values) const;
};


enum {
  IFF_LINEAR_COEFF = 0,
  IFF_CONST_COEFF = 1
};


class VesselFluidExavasationModel : public VesselExavasationModel
{
  const IffParams *params;
public:
  void Init(const VesselList3d &vl_, const ContinuumGrid &grid_, const VesselsInBoxes &vesselsinboxes_, const IffParams &params_);
  bool GetVesselInfo(const Vessel* v, double &coeffa, double &coeffb, double &vala, double &valb) const;
};




/*------------------------------------------------------
------------------------------------------------------*/


class SimpleLymphModel : boost::noncopyable, public VirtualGridFunctions<double, 2>
{
public:
  Array3d<float> theta_tumor;
  IffParams *params;

  virtual void GetValues(int boxindex, const BBox3 &bbox, VirtualGridFunctions<double, 2>::List &values) const
  {
    FOR_BBOX3(p, bbox)
    {
      double ftum   = theta_tumor(p);
      double ftiss  = (1. - ftum);
      double lymph_surf_per_vol = ftiss * params->lymph_surf_per_vol + ftum * params->lymph_surf_per_vol_tumor;
      double lcoeff = -lymph_surf_per_vol*params->lymph_permeability;
      double ccoeff = (-lcoeff) * params->lymph_press;
      values[IFF_LINEAR_COEFF](p) = lcoeff;
      values[IFF_CONST_COEFF](p) = ccoeff;
    }
  }
};


/*------------------------------------------------------
------------------------------------------------------*/


enum {
  IFF_TISS_CONDUCTIVITY = 0,  // total conductivity
  IFF_TISS_PHI_WATER = 1 // water volume fraction
};


/** @note
 * this class is defined header only.
 */
class SimpleTissueModel : boost::noncopyable, public VirtualGridFunctions<double, 2>
{
public:
  Array3d<float> theta_tumor, theta_necro, phi_water;
  IffParams *params;

  // might add the water volume fraction as a result, and rename this class to TissueModel or so ...
  virtual void GetValues(int boxindex, const BBox3 &bbox, VirtualGridFunctions<double, 2>::List &values) const
  {
    FOR_BBOX3(p, bbox)
    {
      double ftum = theta_tumor(p);
      double fnecro = theta_necro(p) * ftum;
      double ftiss = 1.-ftum;
      double coeff = ftum * params->iff_cond_tumor + ftiss * params->iff_cond_tissue + fnecro * params->iff_cond_necro;
      values[IFF_TISS_CONDUCTIVITY](p) = coeff;
      values[IFF_TISS_PHI_WATER](p) = phi_water(p);
    }
  }
};



/**
 * @brief This class calculated iff in tissue
 * @note Future project make these coupled like the o2 calculation
 * 
 * This class models the transport of IF from the 
 * vessels to the surrounding tissue. The exchange with 
 * the plasma in the bloodvessels is neglected. This 
 * means the loss in the vessels is not modeled!
 */

class UncoupledIfModel
{
private:
  const ContinuumGrid *grid;
  const DynArray<BBox3> *mtboxes;
  const VirtualGridFunctions<double, 2> *exavasation_model;
  const VirtualGridFunctions<double, 2> *lymph_model;
  const VirtualGridFunctions<double, 2> *tissue_model;
  const IffParams *params;

  Array3d<double> pfield;

  void insertFieldCoefficients(int box_index, const BBox3 &bbox, FiniteVolumeMatrixBuilder &mb);

public:
  /// Create an UncoupledIfModel
  UncoupledIfModel() : grid(NULL), mtboxes(NULL), exavasation_model(NULL), lymph_model(NULL), tissue_model(NULL) {}
  void Init(const ContinuumGrid &grid_,
            const DynArray<BBox3> &mtboxes_,
            const VirtualGridFunctions<double, 2> &exavasation_model_,
            const VirtualGridFunctions<double, 2> &lymph_model_,
            const VirtualGridFunctions<double, 2> &tissue_model_,
            const IffParams &params_);
  
  /**
   * initialize and compute the pressure field
   */
  void calculate();
  /**
   * "postprocessing":
   * initialize and fill the IfFlowState structure with data determined by flags
   * mtbox, if >= 0, determines which area to fill, corresponding to the mt boxes in ifbase
   */
  void getOutState(IfFlowState &state, int flags);
  /**
   * "postprocessing" :
   * fill the state structure, but only the region indicated by mtbox,
   * mutex is to protect shared data of state since multiple regions may be filled simultaneously
   */
  void fillOutState(IfFlowState& state, int mtbox, tbb::spin_mutex& mutex);
};


/*------------------------------------------------------
------------------------------------------------------*/

/*------------------------------------------------------
------------------------------------------------------*/

#if 0 ///this is work in progress TF
class CoupledIfModel 
{
private:
  const ContinuumGrid *grid;
  const DynArray<BBox3> *mtboxes;
  const VirtualGridFunctions<double, 2> *exavasation_model;
  const VirtualGridFunctions<double, 2> *lymph_model;
  const VirtualGridFunctions<double, 2> *tissue_model;
  const IffParams *params;

  Array3d<double> pfield;

  void insertFieldCoefficients(int box_index, const BBox3 &bbox, FiniteVolumeMatrixBuilder &mb);

public:
  UncoupledIfModel() : grid(NULL), mtboxes(NULL), exavasation_model(NULL), lymph_model(NULL), tissue_model(NULL) {}
  void Init(const ContinuumGrid &grid_,
            const DynArray<BBox3> &mtboxes_,
            const VirtualGridFunctions<double, 2> &exavasation_model_,
            const VirtualGridFunctions<double, 2> &lymph_model_,
            const VirtualGridFunctions<double, 2> &tissue_model_,
            const IffParams &params_);
  
  /*
   * initialize and compute the pressure field
   */
  void calculate();
  /*
   * "postprocessing":
   * initialize and fill the IfFlowState structure with data determined by flags
   * mtbox, if >= 0, determines which area to fill, corresponding to the mt boxes in ifbase
   */
  void getOutState(IfFlowState &state, int flags);
  /*
   * "postprocessing" :
   * fill the state structure, but only the region indicated by mtbox,
   * mutex is to protect shared data of state since multiple regions may be filled simultaneously
   */
  void fillOutState(IfFlowState& state, int mtbox, tbb::spin_mutex& mutex);
};

#endif

#endif

//#endif

