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
#ifndef _OXYGEN_MODEL2_H_
#define _OXYGEN_MODEL2_H_
//#define mwOMP
#include "../python_krebsutils/python_helpers.h"

#include "common/shared-objects.h"
#include "common/continuum-grid.h"
#include "common/vessels3d.h"
#include "common/time_stepper_utils_new.h"

#include "common/trilinos_linsys_construction.h"

#include "common/shared-objects.h"
#include "common/continuum-flow.h"

#include "mwlib/math_ext.h"

#include <boost/math/tools/roots.hpp>
#include <boost/math/tools/tuple.hpp>
#include <boost/foreach.hpp>
#include <exception>
#include <omp.h>

namespace DetailedPO2
{
  
// enum TissueIDs
// {
//   NORMAL = 0,
//   TUMOR  = 1,
//   NECRO  = 2
// };
// enum TumorTypes
// {
//   FAKE = 0,
//   BULKTISSUE = 1
// };
// TumorTypes determineTumorType(h5cpp::Group &tumorgroup);
/**
  * theses are options which I do not yet support by the python input_file_name
  */
struct OutputOptions
{
  bool writeTopoorder = false;
  bool writeAvgPo2 = false;
  bool writeMatrixBuilderInfo = false;
};
/* todo:
 * turn parameters into variables and make it so that parameters can be given from the calling routine
 * */
class Parameters
{

public:
  int getIterationCount(){return iteration_count;};
  void increaseIterationCount(){iteration_count++;};
  Parameters();
  
  void assign(const ptree &pt);
  ptree as_ptree() const;
  
  void UpdateInternalValues();

  // Note: mlO2/ml = mlO2/cm^3 = um^3 O2 / um^3!

  //inline double PHalf()  const { return sat_curve_p50; }
  //inline double MaxConc()  const { return haemoglobin_binding_capacity;  }
  inline double PInit(double r) const { return std::min(po2init_cutoff, po2init_r0 + po2init_dr*r); /*mmHg*/ }

  double Saturation( double p ) const;
  std::pair<double,double> DiffSaturation( double p ) const; // saturation, and derivative of saturation (2nd entry)
  boost::tuple<double, double, double> DiffSaturation2(double p) const; // 2nd derivative
  double DiffSaturationMaxRateOfChange(double p) const;
  // conc is in mlO2/cm^3, and po2 in mmHg
  double BloodPO2ToConc(double p, double h) const;
  // returns concentrations ch = c0*S(p), cp = alpha*p
  std::pair<double, double> BloodPO2ToHematocritAndPlasmaConc(double p, double h) const;
  double ConcToBloodPO2(double conc, double h)  const;

 /**
 * @brief Computes oxygen consumption in units of [ml O2/ ml Tisse /s]
 *
 * Arguments are the tissue partial pressure and the tissue composition.
 * If the boolean michaelis_menten_uptake is set it uses the Michaelis Menten Kurve with
 * corresponding codeded parameters. Otherwise, it uses a linear model, i.e. consumption is proportional to po2.
 */
  //std::pair<double, double> ComputeUptake(double po2, float *tissue_phases, int phases_count) const; // (value, first derivative)
  //void SetTissueParamsByDiffusionRadius(double kdiff_, double alpha_, double rdiff_norm_, double rdiff_tum_, double rdiff_necro_);
  //void writeParametersToHDF(H5::Group &parameter_out_group);
  
  int max_iter;
  double ds2_zeros_[2];
  double ds2_max_;
  double conc_neglect_s;
  //Parameters& operator=(const Parameters&);
  int iteration_count;
  double haemoglobin_binding_capacity;
  double solubility_plasma, solubility_tissue; /* mlO2/cm^3/mmHg = (mlO2 / ml tissue) / mmHg */
  double D_plasma;
  double rd_norm;
  double rd_tum;
  double rd_necro;
  double po2init_r0, po2init_dr, po2init_cutoff;
  //double po2_kdiff;
  double sat_curve_exponent; // the n in S = p^n / (p^n + p50^n)
  double sat_curve_p50;
  double axial_integration_step_factor; // measured in fraction of numerical grid cell size.
  double transvascular_ring_size; // measured in fraction of numerical grid cell size.
  string debug_fn;
  bool debug_zero_o2field; // zero out the o2 field
  bool michaelis_menten_uptake;
  bool useCellBasedUptake;
  double po2_mmcons_k_norm;
  double po2_mmcons_k_tum;
  double po2_mmcons_k_necro;
  double po2_mmcons_m0_norm;
  double po2_mmcons_m0_tum;
  double po2_mmcons_m0_necro;
  //double po2_mmcons_m0[3], po2_mmcons_k[3];
  int tissue_boundary_condition_flags;
  double tissue_boundary_value;
  double extra_tissue_source_linear; // this is for comparison with moschandreou (2011), where a "background" of supplying capillaries must be around
  double extra_tissue_source_const; // dito
  double conductivity_coeff1, conductivity_coeff2, conductivity_coeff3; // for transvascular transport
  bool approximateInsignificantTransvascularFlux;
  int massTransferCoefficientModelNumber;
  //general simulation parameters
  int loglevel;
  int num_threads;
  string detailedO2name;
  double convergence_tolerance; // field and vessel po2 differences from iteration to iteration must both be lower than this for the iteration to stop
  string tissue_po2_boundary_condition;
  string input_file_name;
  string input_group_path;
  string output_file_name;
  string output_group_path;
  string tumor_file_name;
  string tumor_group_path;
  string vessel_group_path;
  OutputOptions currentOutputOptions;
};


// partial pressure in blood plasma at the endpoints of each segment
typedef DynArray<Float2> VesselPO2Storage;




//void ComputePO2(const Parameters &params, VesselList3d& vl, ContinuumGrid &grid, DomainDecomposition &mtboxes, Array3df &po2field, VesselPO2Storage &storage, const TissuePhases &phases, ptree &metadata, bool world);
double ComputeCircumferentialMassTransferCoeff(const Parameters &params, double r);

typedef Eigen::Matrix<float, 5, 1> VesselPO2SolutionRecord; //x, po2, ext_po2, conc_flux, dS/dx;

class VascularPO2PropagationModel;

struct Measurement : boost::noncopyable
{
  const Parameters params;
  std::auto_ptr<VesselList3d> vl;
  Array3df po2field;
  ContinuumGrid grid;
  VesselPO2Storage vesselpo2;
  const TissuePhases &phases;
  std::auto_ptr<VascularPO2PropagationModel> vascularPO2PropagationModel;
  
  double o2mass_in_root, o2mass_out_root, o2mass_out_vessels, o2mass_consumed_transvascular;
  
  Measurement(std::auto_ptr<VesselList3d> vl_ptr_, const Parameters &params, const ContinuumGrid& grid_, const TissuePhases &phases_, Array3df &po2field_, VesselPO2Storage &vesselpo2_);
  void computeVesselSolution(int idx, DynArray<VesselPO2SolutionRecord> &sol, bool world);
};


void TestSaturationCurve();
void TestSingleVesselPO2Integration();

struct DetailedPO2Sim : public boost::noncopyable
{
  bool world;
  std::shared_ptr<VesselList3d> vl;
  Parameters params;
  BloodFlowParameters bfparams;
  
  ContinuumGrid grid;
  DomainDecomposition mtboxes;
  LatticeDataQuad3d ld;
  FloatBBox3 wbox;
  Float3 worldCenter;
  Float3 gridCenter;
  
  double po2_cons_coeff[3];
  void SetTissueParamsByDiffusionRadius();
  std::pair<double, double> ComputeUptake(double po2, float *tissue_phases, int phases_count) const; // (value, first derivative)
  void ComputePo2Field( 
		     const ContinuumGrid &grid, 
		     DomainDecomposition &mtboxes, 
		     const TissuePhases &phases, 
		     Array3df po2field,
         boost::optional<Array3d<float>> cell_based_o2_uptake,
		     FiniteVolumeMatrixBuilder &mb,  
		     bool keep_preconditioner);
  //sets up the linear trilionos matrix system, builder is implemented as struct
  //could use for example different stencils
  FiniteVolumeMatrixBuilder tissue_diff_matrix_builder;
  
  boost::optional<Array3df> cell_based_o2_uptake;
  Array3df po2field;
  DetailedPO2::VesselPO2Storage po2vessels;
  ptree metadata;
  // Get to know where tumor is and where normal tissue is.
  // I.e. get volume fractions of each cell type.
  // after this call the 3D field phases is filled with
  // 3 vallues giving the portion of corresponding tissue type
  TissuePhases phases;//Declaration
  void init(
            BloodFlowParameters &bfparams, 
            //VesselList3d &vl, 
            double grid_lattice_const, 
            double safety_layer_size, 
            boost::optional<Int3> grid_lattice_size, 
            boost::optional<H5::Group> tumorgroup,
            boost::optional<Array3df> previous_po2field, boost::optional<DetailedPO2::VesselPO2Storage> previous_po2vessels, boost::optional<Array3d<float>> cell_based_o2_uptake);
  int run();
  void PrepareNetworkInfo(const VesselList3d &vl, DynArray<const Vessel*> &sorted_vessels, DynArray<const VesselNode*> &roots);
  
  DynArray<const Vessel*> sorted_vessels;
  DynArray<const VesselNode*> roots;
  DynArray<const VesselNode*> arterial_roots;
  
  Array3df getPo2field();
  DetailedPO2::VesselPO2Storage getVesselPO2Storrage();
//   void WriteOutput(H5::Group &basegroup,
//                  const VesselList3d &vl,
//                  const Parameters &params,
//                  const boost::optional<const VesselPO2Storage&> vesselpo2,
//                  const boost::optional<DynArray<const Vessel*>&> sorted_vessels,
//                  const boost::optional<ContinuumGrid&> grid,
//                  const boost::optional<Array3df> po2field,
//                  const boost::optional<const FiniteVolumeMatrixBuilder&> mbopt);
  //void WriteOutput_new(H5::H5File &outfile);
  void WriteDataOutput(H5::Group &out_grp);
  void WriteParametersToHDF(H5::Group &out_params);  
  /**
 * @brief integrates the po2 along a vessel
 * 
 * The hiarachial ordering of the vessel segments is provided, in the sorted_vessels array.
 * As Starting point of the integration, an initial value dependent on the inlet hemoglobin and
 * arterial root pressure
 * @note this function is not very parallel
 */
void IntegrateVesselPO2(const Parameters &params, 
			VesselPO2Storage &vesselpo2,
			DynArray<const Vessel*> &sorted_vessels,
			DynArray<const VesselNode*> &arterial_roots,
			const ContinuumGrid &grid,
			const Array3df extpo2,
			const TissuePhases &phases,
			FiniteVolumeMatrixBuilder &matrix_builder,
			bool world
 		      );
};
// template<class T>
// static T checkedExtractFromDict(const py::dict &d, const char* name);
void InitParameters(DetailedPO2::Parameters &params, const py::dict &py_parameters);
};//end namespace

#endif
