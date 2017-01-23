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


#include "../python_krebsutils/python-helpers.h"
#include "numpy.hpp"
#include "oxygen_model2.h"
#include "calcflow.h"

#ifdef EPETRA_MPI
  #include "mpi.h"
#endif

#define BOOST_RESULT_OF_USE_DECLTYPE 1

#include <boost/foreach.hpp>
#include <boost/python/errors.hpp>
#include <boost/iterator/transform_iterator.hpp>

namespace py = boost::python;
namespace np = boost::python::numpy;
namespace nm = boost::python::numeric;
namespace h5 = h5cpp;

/**
 * @brief Sets tissue phases on the lattice sites
 * 
 * Checks wether a tumor calculation is in the h5 file or not.
 * In case it is, the tissue phases are interpolated to lattice points.
 * Otherwise it starts with a homogeneous distribution of normal tissue.
 */

namespace DetailedPO2
{

void SetupTissuePhases(DetailedPO2::TissuePhases &phases, const ContinuumGrid &grid, DomainDecomposition &mtboxes, py::object py_tumorgroup)
{
  // there is a tumor group provided by the python side
  if (!py_tumorgroup.is_none())
  {
    h5cpp::Group   tumorgroup      = PythonToCppGroup(py_tumorgroup);
    //if (tumorgroup.attrs().exists("TYPE") && tumorgroup.attrs().get<string>("TYPE")=="faketumor")
    if( determineTumorType(tumorgroup) == TumorTypes::FAKE)
    {
      double tumor_radius = tumorgroup.attrs().get<double>("TUMOR_RADIUS");
      phases = DetailedPO2::TissuePhases(2, grid.Box());
      #pragma omp parallel
      {
        BOOST_FOREACH(const BBox3 &bbox, mtboxes.getCurrentThreadRange())
        {
          FOR_BBOX3(p, bbox)
          {
            Float3 wp =  grid.ld.LatticeToWorld(p);
            double r = wp.norm();
            double f = my::smooth_heaviside_sin(-r+tumor_radius, grid.Spacing()); // is = 1 within the tumor
            phases.phase_arrays[DetailedPO2::NORMAL](p) = 1.f-f;
            phases.phase_arrays[DetailedPO2::TUMOR](p) = f;
          }
        }
      }
    }
    else if (determineTumorType(tumorgroup) == TumorTypes::BULKTISSUE)
    {
      h5cpp::Dataset cell_vol_fraction_ds = tumorgroup.open_dataset("conc");
      h5cpp::Dataset tumor_fraction_ds = tumorgroup.open_dataset("ptc");
      h5cpp::Dataset necro_fraction_ds = tumorgroup.open_dataset("necro");
      h5cpp::Group   ldgroup        = tumor_fraction_ds.get_file().root().open_group(tumor_fraction_ds.attrs().get<string>("LATTICE_PATH"));

      Array3df cell_vol_fraction;
      Array3df tumor_fraction;
      Array3df necro_fraction;
      LatticeDataQuad3d fieldld;
      ReadHdfLd(ldgroup, fieldld);
      ReadArray3D(cell_vol_fraction_ds, cell_vol_fraction);
      ReadArray3D(tumor_fraction_ds, tumor_fraction);
      ReadArray3D(necro_fraction_ds, necro_fraction);

      phases = DetailedPO2::TissuePhases(3, grid.Box());

      #pragma omp parallel
      {
        BOOST_FOREACH(const BBox3 &bbox, mtboxes.getCurrentThreadRange())
        {
          FOR_BBOX3(p, bbox)
          {
	    //interpolate fields to point p on the lattice
	    float this_cell_vol_fraction = FieldInterpolate::ValueAveraged(cell_vol_fraction, fieldld, FieldInterpolate::Const(0.f), grid.ld.LatticeToWorld(p));
            float this_tumor_fraction = FieldInterpolate::ValueAveraged(tumor_fraction, fieldld, FieldInterpolate::Const(0.f), grid.ld.LatticeToWorld(p));
            float this_necro_fraction = FieldInterpolate::ValueAveraged(necro_fraction, fieldld, FieldInterpolate::Const(0.f), grid.ld.LatticeToWorld(p));
	    
	    phases.phase_arrays[DetailedPO2::NORMAL](p) = (1.f-this_tumor_fraction)*this_cell_vol_fraction;
            phases.phase_arrays[DetailedPO2::TUMOR](p) = this_tumor_fraction*this_cell_vol_fraction;
	    phases.phase_arrays[DetailedPO2::NECRO](p) = this_necro_fraction*this_cell_vol_fraction;
          }
        }
      }
    }
  }
  // if there is no tumorgroup provided by the python side
  // we fill everything with normal tissue
  else 
  {
    phases = DetailedPO2::TissuePhases(1, grid.Box());
    phases.phase_arrays[DetailedPO2::NORMAL].fill(1.);
  }
}

template<class T>
static T checkedExtractFromDict(const py::dict &d, const char* name)
{
  try
  {
    return py::extract<T>(d.get(name));
  }
  catch (const py::error_already_set &e) 
  {
    std::cerr << format("unable to extract parameter '%s': ") % name;
    throw e; // don't every try to handle this!
  }
}

void InitParameters(DetailedPO2::Parameters &params, py::dict py_parameters)
{
#define GET_PO2_PARAM_FROM_DICT(TYPE, NAME) checkedExtractFromDict<TYPE>(py_parameters, NAME);
#define GET_PO2_PARAM_IF_NONNONE(TARGET, TYPE, NAME) { py::object o(py_parameters.get(NAME)); if (!o.is_none()) TARGET=py::extract<TYPE>(o); }
  params.max_iter = GET_PO2_PARAM_FROM_DICT(int, "max_iter");
  params.sat_curve_exponent = GET_PO2_PARAM_FROM_DICT(double, "S_n");
  params.sat_curve_p50 = GET_PO2_PARAM_FROM_DICT(double, "S_p50");
  GET_PO2_PARAM_IF_NONNONE(params.axial_integration_step_factor, double, "axial_integration_step_factor");
  GET_PO2_PARAM_IF_NONNONE(params.debug_zero_o2field, bool, "debug_zero_o2field");
  GET_PO2_PARAM_IF_NONNONE(params.debug_fn, string, "debug_fn");
  GET_PO2_PARAM_IF_NONNONE(params.transvascular_ring_size, double, "transvascular_ring_size");
  double kd, rd_norm = 100., rd_tum = 100., rd_necro = 100.;
  kd = GET_PO2_PARAM_FROM_DICT(double, "D_tissue");
  params.tissue_solubility = GET_PO2_PARAM_FROM_DICT(double, "solubility_tissue");
  params.plasma_solubility = GET_PO2_PARAM_FROM_DICT(double, "solubility_plasma");
//   try 
//   {
//     kd = GET_PO2_PARAM_FROM_DICT(double, "D_tissue");
//   }
//   catch (py::error_already_set &e) // load legacy name
//   {
//     std::cerr << "fallback to kD_tissue" << endl;
//     kd = GET_PO2_PARAM_FROM_DICT(double, "kD_tissue");
//   }
//   try
//   {
//     params.tissue_solubility = GET_PO2_PARAM_FROM_DICT(double, "solubility_tissue");
//   }
//   catch (py::error_already_set &e) // load legacy name
//   {
//     std::cerr << "fallback to alpha_t" << endl;
//     params.tissue_solubility = GET_PO2_PARAM_FROM_DICT(double, "alpha_t");
//   }
//   try 
//   {
//     params.plasma_solubility = GET_PO2_PARAM_FROM_DICT(double, "solubility_plasma");
//   }
//   catch (py::error_already_set &e) // load legacy name
//   {
//     std::cerr << "fallback to alpha_p" << endl;
//     params.plasma_solubility = GET_PO2_PARAM_FROM_DICT(double, "alpha_p");
//   }
  
  GET_PO2_PARAM_IF_NONNONE(rd_norm, double, "rd_norm");
  GET_PO2_PARAM_IF_NONNONE(rd_tum, double, "rd_tum");
  GET_PO2_PARAM_IF_NONNONE(rd_necro, double, "rd_necro");
  //rd_tum = GET_PO2_PARAM_FROM_DICT(double, "rd_tum");
  //rd_necro = GET_PO2_PARAM_FROM_DICT(double, "rd_necro");
  params.SetTissueParamsByDiffusionRadius(kd, params.tissue_solubility, rd_norm, rd_tum, rd_necro);
  params.po2init_r0 = GET_PO2_PARAM_FROM_DICT(double, "po2init_r0");
  params.po2init_dr = GET_PO2_PARAM_FROM_DICT(double, "po2init_dr");
  params.po2init_cutoff = GET_PO2_PARAM_FROM_DICT(double, "po2init_cutoff");
  GET_PO2_PARAM_IF_NONNONE(params.michaelis_menten_uptake, bool, "michaelis_menten_uptake");
  GET_PO2_PARAM_IF_NONNONE(params.po2_mmcons_k[DetailedPO2::NORMAL], double, "mmcons_k_norm");
  GET_PO2_PARAM_IF_NONNONE(params.po2_mmcons_k[DetailedPO2::TUMOR], double, "mmcons_k_tum");
  GET_PO2_PARAM_IF_NONNONE(params.po2_mmcons_k[DetailedPO2::NECRO], double, "mmcons_k_necro");
  GET_PO2_PARAM_IF_NONNONE(params.po2_mmcons_m0[DetailedPO2::NORMAL], double, "mmcons_m0_norm");
  GET_PO2_PARAM_IF_NONNONE(params.po2_mmcons_m0[DetailedPO2::TUMOR], double, "mmcons_m0_tum");
  GET_PO2_PARAM_IF_NONNONE(params.po2_mmcons_m0[DetailedPO2::NECRO], double, "mmcons_m0_necro");
  GET_PO2_PARAM_IF_NONNONE(params.haemoglobin_binding_capacity, double, "haemoglobin_binding_capacity");
  string bc_type("neumann");
  GET_PO2_PARAM_IF_NONNONE(bc_type, string, "tissue_po2_boundary_condition");
  if (bc_type == "dirichlet_x") 
    params.tissue_boundary_condition_flags = 1; //FiniteVolumeMatrixBuilder;
  else if (bc_type == "dirichlet_yz") 
    params.tissue_boundary_condition_flags = 2;
  else if (bc_type == "dirichlet") 
    params.tissue_boundary_condition_flags = 3;
  else if (bc_type == "neumann") 
    params.tissue_boundary_condition_flags = 0;
  else 
    throw std::invalid_argument("tissue_po2_boundary_condition must be 'dirichlet','dirichlet_x', 'dirichlet_yz' or 'neumann'");
  GET_PO2_PARAM_IF_NONNONE(params.tissue_boundary_value, double, "tissue_boundary_value");
  GET_PO2_PARAM_IF_NONNONE(params.extra_tissue_source_linear, double, "extra_tissue_source_linear");
  GET_PO2_PARAM_IF_NONNONE(params.extra_tissue_source_const, double, "extra_tissue_source_const");
  // required parameters for transvascular transport
  params.conductivity_coeff1 = GET_PO2_PARAM_FROM_DICT(double, "conductivity_coeff1"); 
  params.conductivity_coeff2 = GET_PO2_PARAM_FROM_DICT(double, "conductivity_coeff2");
  params.conductivity_coeff3 = GET_PO2_PARAM_FROM_DICT(double, "conductivity_coeff3");
  GET_PO2_PARAM_IF_NONNONE(params.convergence_tolerance, double, "convergence_tolerance");
  GET_PO2_PARAM_IF_NONNONE(params.loglevel, int, "loglevel");
  GET_PO2_PARAM_IF_NONNONE(params.approximateInsignificantTransvascularFlux, bool, "approximateInsignificantTransvascularFlux");
  GET_PO2_PARAM_IF_NONNONE(params.D_plasma, double, "D_plasma");
  GET_PO2_PARAM_IF_NONNONE(params.massTransferCoefficientModelNumber, int, "massTransferCoefficientModelNumber");
#undef GET_PO2_PARAM_FROM_DICT
#undef GET_PO2_PARAM_IF_NONNONE
  params.UpdateInternalValues();
}


template<class T>
inline boost::optional<T> getOptional(const char* name, py::dict &d)
{
  py::object o = d.get(name);
  if (!o.is_none())
    return boost::optional<T>(py::extract<T>(o));
  else
    return boost::optional<T>();
}


/** 
 * @brief Callback for the main functions
 * 
 * Calls the important stuff.
 */
static void PyComputePO2(py::object py_vesselgroup, py::object py_tumorgroup, py::dict py_parameters, py::object py_bfparams, py::object py_h5outputGroup)
{
  bool world = false;
  DetailedPO2::Parameters params;
  InitParameters(params, py_parameters);
  
  //h5cpp::Group group = PythonToCppGroup(py_group);
  //h5cpp::Group vesselgroup = group.open_group(path_vessels);
  h5cpp::Group vesselgroup = PythonToCppGroup(py_vesselgroup);
  
  // THIIIIRYYYYY, filter muss = false sein sonst stimmt in der Ausgabe in der Hdf5 Datei die Anzahl der Vessels nicht mehr mit den daten im recomputed_flow Verzeichnis ueberein!
  std::auto_ptr<VesselList3d> vl = ReadVesselList3d(vesselgroup, make_ptree("filter",false));
  
  
  ContinuumGrid grid;
  DomainDecomposition mtboxes;
  double grid_lattice_const              = py::extract<double>(py_parameters.get("grid_lattice_const", 30.));
  double safety_layer_size               = py::extract<double>(py_parameters.get("safety_layer_size", grid_lattice_const*3.));
  boost::optional<Int3> grid_lattice_size = getOptional<Int3>("grid_lattice_size", py_parameters);
  {
<<<<<<< HEAD
    int dim = (::Size(vl->Ld().Box())[2]<=1) ? 2 : 3;
=======
    //this worked only for lattices
    //int dim = (::Size(vl->Ld().Box())[2]<=1) ? 2 : 3;
    int dim=0;
    if (world)
    {
      dim = 3;
    }
    else
    {
      dim = (::Size(vl->Ld().Box())[2]<=1) ? 2 : 3;
    }
>>>>>>> b91767e... bugfix: define dimension in global scope
    LatticeDataQuad3d ld;
    if (grid_lattice_size)
    {
      ld.Init(*grid_lattice_size, grid_lattice_const);
      ld.SetCellCentering(Bool3(1, 1, dim>2));
      FloatBBox3 wbox = vl->Ld().GetWorldBox();
      Float3 worldCenter = 0.5*(wbox.max + wbox.min);
      Float3 gridCenter = 0.5*(ld.GetWorldBox().max + ld.GetWorldBox().min);
      ld.SetOriginPosition(ld.GetOriginPosition() + (worldCenter - gridCenter));
    }
    else
    {
      //added safety space to reduce boundary errors
      SetupFieldLattice(vl->Ld().GetWorldBox(), dim, grid_lattice_const, safety_layer_size, ld); 
    }
    grid.init(ld, dim);
    mtboxes.init(MakeMtBoxGrid(grid.Box(), Int3(32, 32, 32)));
    if (params.loglevel > 0)
    {
      cout << "continuum grid:" << endl;
      grid.ld.print(cout);
      cout << endl;
    }
  }

  if (params.loglevel > 0)
  {
    cout << "vessel lattice" << endl;
    vl->Ld().print(cout);
    cout << endl;
  }

  //if (!py_bfparams.is_none())
  if (py_bfparams)
  {
    BloodFlowParameters bfparams = py::extract<BloodFlowParameters>(py_bfparams);
    CalcFlow(*vl, bfparams);
  }
  isVesselListGood(*vl);

  Array3df po2field;
  DetailedPO2::VesselPO2Storage po2vessels;
  ptree metadata;

  //cout << format("in c++: %.20f %.20f %.20f\n") % params.conductivity_coeff1 % params.conductivity_coeff2 % params.conductivity_coeff_gamma;
  
  // Get to know where tumor is and where normal tissue is.
  // I.e. get volume fractions of each cell type.
  // after this call the 3D field phases is filled with
  // 3 vallues giving the portion of corresponding tissue type
  DetailedPO2::TissuePhases phases;//Declaration
  SetupTissuePhases(phases, grid, mtboxes, py_tumorgroup);//filling
  
  { //ReleaseGIL unlock(); // allow the python interpreter to do things while this is running
  /**
   * MOST IMPORTANT CALL
   * grid is conitnuum grid, mtboxes is decomposition into threads
   */
  DetailedPO2::ComputePO2(params, *vl, grid, mtboxes, po2field, po2vessels, phases, metadata, world);
  }
  if (PyErr_Occurred() != NULL) return; // don't save stuff
  
  //This writes the results back to python
#if 0
  Int2 sz(2,po2vessels.size());
  nm::array  py_po2vessels = np::copy<float,2>(sz.data(), (float*)po2vessels.data(), calc_strides::first_dim_varies_fastest(sz).data());
  nm::array  py_po2field   = np::copy<float,3>(po2field.size().data(), po2field.getPtr(), po2field.strides().data());
  py::object py_ld(grid.ld);
  
  return py::make_tuple(py_po2vessels, py_ld, py_po2field);
#endif
  {
    h5cpp::Group outputGroup = PythonToCppGroup(py_h5outputGroup);
    h5cpp::Group ldgroup = outputGroup.create_group("field_ld");
    WriteHdfLd(ldgroup, grid.ld);
    WriteScalarField<float>(outputGroup, "po2field", po2field, grid.ld, ldgroup);
    h5cpp::create_dataset<float>(outputGroup, "po2vessels", h5cpp::Dataspace::simple_dims(po2vessels.size(), 2), (float*)po2vessels[0].data(), h5cpp::CREATE_DS_COMPRESSED); // FIX ME: transpose the array!
    WriteHdfPtree(outputGroup, metadata, HDF_WRITE_PTREE_AS_DATASETS);
  }
}


static py::object PyComputeSaturation(nm::array py_po2, py::dict py_parameters)
{
  DetailedPO2::Parameters params;
  InitParameters(params, py_parameters);
  
  np::arrayt<float> po2(py_po2);

  if (!(po2.rank() == 1 && po2.isCContiguous())) throw std::invalid_argument("rank 1 and contiguous expected");
  np::arrayt<float> result(np::empty(1, po2.shape(), np::getItemtype<float>()));

  for (int i=0; i<po2.shape()[0]; ++i)
  {
    result(i) = params.Saturation(po2(i));
  }
  return result.getObject();
}


static py::object PyComputeMassTransferCoefficient(nm::array py_radius, py::dict py_parameters)
{
  DetailedPO2::Parameters params;
  InitParameters(params, py_parameters);
  
  np::arrayt<float> radius(py_radius);
  if (!(radius.rank() == 1 && radius.isCContiguous())) throw std::invalid_argument("rank 1 and contiguous expected");
  
  np::arrayt<double> mtc(np::empty(1, radius.shape(), np::getItemtype<double>()));
  for (int i=0; i<radius.shape()[0]; ++i)
  {
    double r = radius(i);
    if (r <= 0)
      mtc(i) = 0;
    else
    {
        mtc(i) = (1.0/(my::mconst::pi2()*r)) * ComputeCircumferentialMassTransferCoeff(params, r);
    }
  }
  return mtc.getObject();
}


static py::object PyComputeConcentration(nm::array py_po2, nm::array py_hematocrit, py::dict py_parameters)
{
  DetailedPO2::Parameters params;
  InitParameters(params, py_parameters);
  
  np::arrayt<float> po2(py_po2);
  np::arrayt<float> hema(py_hematocrit);

  if (!(po2.rank() == 1 && po2.isCContiguous())) throw std::invalid_argument("rank 1 and contiguous expected");
  if (!(hema.rank() == 1 && hema.isCContiguous())) throw std::invalid_argument("rank 1 and contiguous expected");

  np::arrayt<float> result(np::empty(1, po2.shape(), np::getItemtype<float>()));

  for (int i=0; i<po2.shape()[0]; ++i)
  {
    result(i) = params.BloodPO2ToConc(po2(i), hema(i));
  }
  return result.getObject();
}


static py::object PyComputePO2FromConc(nm::array py_conc, nm::array py_hematocrit, py::dict py_parameters)
{
  DetailedPO2::Parameters params;
  InitParameters(params, py_parameters);
  
  np::arrayt<float> conc(py_conc);
  np::arrayt<float> hema(py_hematocrit);
  
  if (!(conc.rank() == 1 && conc.isCContiguous())) throw std::invalid_argument("rank 1 and contiguous expected");
  if (!(hema.rank() == 1 && hema.isCContiguous())) throw std::invalid_argument("rank 1 and contiguous expected");
  
  np::arrayt<float> result(np::empty(1, conc.shape(), np::getItemtype<float>()));
  
  for (int i=0; i<conc.shape()[0]; ++i)
  {
    result(i) = params.ConcToBloodPO2(conc(i), hema(i));
  }
  return result.getObject();
}


static py::object PyComputeUptake(nm::array py_po2field, const LatticeDataQuad3d &field_ld, py::object py_tumorgroup,  py::dict py_parameters)
{
  DetailedPO2::Parameters params;
  InitParameters(params, py_parameters);
  ContinuumGrid grid;
  grid.init(field_ld, field_ld.Size()[2]<=1 ? 2 : 3);
  DomainDecomposition mtboxes;
  mtboxes.init(MakeMtBoxGrid(grid.Box(), Int3(32, 32, 32)));
  DetailedPO2::TissuePhases phases;
  SetupTissuePhases(phases, grid, mtboxes, py_tumorgroup);
  np::arrayt<float> po2field(py_po2field);
  np::arrayt<float> consumption(np::empty(3, po2field.shape(), np::getItemtype<float>()));
  #pragma omp parallel
  {
    const double extra_factor = params.tissue_solubility;
    BOOST_FOREACH(const BBox3 &bbox, mtboxes.getCurrentThreadRange())
    {
      FOR_BBOX3(p, bbox)
      {
        double po2 = po2field(p[0],p[1],p[2]);
        Float3 phases_loc = phases(p);
        double m, dm;
        boost::tie(m, dm) = params.ComputeUptake(po2, phases_loc.data(), phases.count);
        double mextra = -extra_factor*(params.extra_tissue_source_linear * po2 + params.extra_tissue_source_const);
        m += mextra;
        consumption(p[0],p[1],p[2]) = m; // now already reported in [mlO2 / mlTissue]
      }
    }
  }

  return consumption.getObject();
}


template<int rows>
static Eigen::Matrix<float, rows, 1> LinearInterpolation(float xeval, const DynArray<Eigen::Matrix<float, rows, 1> > &sol)
{
  typedef Eigen::Matrix<float, rows, 1> Record;
  float x0 = sol[0][0];
  float x1 = sol[sol.size()-1][0];
  myAssert(xeval>=x0 && xeval<=x1);

  struct getx { 
    typedef float result_type;
    float operator()(const Record &r) const { return r[0]; } 
  };
  
  auto it_begin = boost::make_transform_iterator(sol.begin(), getx());
  auto it_end   = boost::make_transform_iterator(sol.end(), getx());
  auto iter = std::upper_bound<>(it_begin, it_end, xeval);

  if (iter == it_end)
    return sol[sol.size()-1];
  if (iter == it_begin)
    return sol[0];

  int idx = iter - it_begin;
  const Record &r0 = sol[idx-1];
  const Record &r1 = sol[idx];
  float f = (xeval-r0[0])/(r1[0]-r0[0]);
  myAssert(f>=-0.01 && f<=1.001);
  return (1.-f)*r0 + f*r1;
}


// may be a measurement class can come back later when it makes more sense to store persistent data between analysis steps
py::object PySampleVessels(py::object py_vesselgroup, py::object py_tumorgroup, py::dict py_parameters, nm::array py_vesselpo2, nm::array py_po2field, const LatticeDataQuad3d &field_ld, float sample_len)
{
  bool world = false;
  h5cpp::Group vesselgroup = PythonToCppGroup(py_vesselgroup);
  std::auto_ptr<VesselList3d> vl_ = ReadVesselList3d(vesselgroup, make_ptree("filter",false));
  const VesselList3d &vl = *vl_;
  
  DetailedPO2::VesselPO2Storage vesselpo2(vl.GetECount());
  Array3df po2field(field_ld.Box());

  ContinuumGrid grid(field_ld);
  DomainDecomposition mtboxes(MakeMtBoxGrid(grid.Box(), Int3(32, 32, 32)));
  
  DetailedPO2::TissuePhases phases;//Declaration
  SetupTissuePhases(phases, grid, mtboxes, py_tumorgroup);//filling
  
  np::copy<float,2>((float*)get_ptr(vesselpo2), Int2(2,vesselpo2.size()).data(), calc_strides::first_dim_varies_fastest(Int2(2, vesselpo2.size())).data(), py_vesselpo2);
  np::copy<float,3>(po2field.getPtr(), po2field.size().data(), po2field.strides().data(), py_po2field);
  
  DetailedPO2::Parameters params;
  InitParameters(params, py_parameters);

  DetailedPO2::Measurement measure(vl_, params, grid, phases, po2field, vesselpo2);
  // vl_ is now NULL!
  
  /* VesselPO2SolutionRecord; //x, po2, ext_po2, conc_flux, dS/dx; 
      I dont know what is x at the moment since it is a single float
  */
  typedef DetailedPO2::VesselPO2SolutionRecord Rec;

  int cnt = vl.GetECount();

  DynArray<float> tmp(1024, ConsTags::RESERVE);//buffer to write back to python
  DynArray<Rec> sol; // vessel solution buffer

  int num_total_samples = 0;

  CylinderNetworkSampler sampler;
  sampler.Init(sample_len, ptree());

  for(int i=0; i<cnt && !my::checkAbort(); ++i)
  {
    const Vessel* v = vl.GetEdge(i);
    Float3 p0 = vl.Ld().LatticeToWorld(v->LPosA()),
          p1 = vl.Ld().LatticeToWorld(v->LPosB());
    //cout << format("visit edge %s %s") % p0 % p1 << endl;

    sampler.Set(p0, p1, 0.); // 3rd arg is the radius
    int num_samples = sampler.GenerateLineSamples();

    measure.computeVesselSolution(i, sol, world); // this is done for every vessel
//     cout << "sol for vessel " << i << endl;
//     for (int j=0; j<sol.size(); ++j)
//       cout << sol[j] << endl;

    for(int j=0; j<num_samples; ++j)
    {
      CylinderNetworkSampler::Sample ls = sampler.GetSample(j);

      auto smpl = (sol.size()>1) ? LinearInterpolation(ls.fpos[2]*(sol[sol.size()-1][0]-sol[0][0]), sol) : DetailedPO2::VesselPO2SolutionRecord::Zero();

      for (int i=1; i<Rec::RowsAtCompileTime; ++i)//note not from 0
        tmp.push_back(smpl[i]);
      tmp.push_back(sampler.weight); // weight is the length per sample point
    }
    num_total_samples += num_samples;
  }

  enum { ncomps = Rec::RowsAtCompileTime };
  Int2 sz(ncomps-1+1, num_total_samples);
  np::arrayt<float> samples = np::copy<float,2>(sz.data(), get_ptr(tmp), calc_strides::first_dim_varies_fastest(sz).data());

  py::tuple fluxes = py::make_tuple(measure.o2mass_in_root,
                                    measure.o2mass_out_root,
                                    measure.o2mass_out_vessels,
                                    measure.o2mass_consumed_transvascular);

  return py::make_tuple(samples, fluxes);
}

// may be a measurement class can come back later when it makes more sense to store persistent data between analysis steps
py::object PySampleVesselsWorld(py::object py_vesselgroup, py::object py_tumorgroup, py::dict py_parameters, np::arrayt<float> py_vesselpo2, np::arrayt<float> py_po2field, const LatticeDataQuad3d &field_ld, float sample_len)
{
  bool world = true;
  h5cpp::Group vesselgroup = PythonToCppGroup(py_vesselgroup);
  std::auto_ptr<VesselList3d> vl_ = ReadVesselList3d(vesselgroup, make_ptree("filter",false));
  const VesselList3d &vl = *vl_;
  
  DetailedPO2::VesselPO2Storage vesselpo2(vl.GetECount());
  Array3df po2field(field_ld.Box());

  ContinuumGrid grid(field_ld);
  DomainDecomposition mtboxes(MakeMtBoxGrid(grid.Box(), Int3(32, 32, 32)));
  
  DetailedPO2::TissuePhases phases;//Declaration
  SetupTissuePhases(phases, grid, mtboxes, py_tumorgroup);//filling
  
  np::copy<float,2>((float*)get_ptr(vesselpo2), Int2(2,vesselpo2.size()).data(), calc_strides::first_dim_varies_fastest(Int2(2, vesselpo2.size())).data(), py::object(py_vesselpo2));
  np::copy<float,3>(po2field.getPtr(), po2field.size().data(), po2field.strides().data(), py::object(py_po2field));
  
  DetailedPO2::Parameters params;
  InitParameters(params, py_parameters);

  DetailedPO2::Measurement measure(vl_, params, grid, phases, po2field, vesselpo2);
  // vl_ is now NULL!
  
  // here we go (sampling)
  typedef DetailedPO2::VesselPO2SolutionRecord Rec;

  int cnt = vl.GetECount();

  DynArray<float> tmp(1024, ConsTags::RESERVE);
  DynArray<Rec> sol; // vessel solution buffer

  int num_total_samples = 0;

  CylinderNetworkSampler sampler;
  sampler.Init(sample_len, ptree());

  for(int i=0; i<cnt && !my::checkAbort(); ++i)
  {
    const Vessel* v = vl.GetEdge(i);
    Float3 p0 = v->NodeA()->worldpos;
    Float3 p1 = v->NodeB()->worldpos;
   
    //cout << format("visit edge %s %s") % p0 % p1 << endl;

    sampler.Set(p0, p1, 0.); // 3rd arg is the radius
    int num_samples = sampler.GenerateLineSamples();

    measure.computeVesselSolution(i, sol, world);
//     cout << "sol for vessel " << i << endl;
//     for (int j=0; j<sol.size(); ++j)
//       cout << sol[j] << endl;

    for(int j=0; j<num_samples; ++j)
    {
      CylinderNetworkSampler::Sample ls = sampler.GetSample(j);

      auto smpl = (sol.size()>1) ? LinearInterpolation(ls.fpos[2]*(sol[sol.size()-1][0]-sol[0][0]), sol) : DetailedPO2::VesselPO2SolutionRecord::Zero();

      for (int i=1; i<Rec::RowsAtCompileTime; ++i)
        tmp.push_back(smpl[i]);
      tmp.push_back(sampler.weight); // weight is the length per sample point
    }
    num_total_samples += num_samples;
  }

  enum { ncomps = Rec::RowsAtCompileTime };
  Int2 sz(ncomps-1+1, num_total_samples);
  py::object samples = np::copy<float,2>(sz.data(), get_ptr(tmp), calc_strides::first_dim_varies_fastest(sz).data());

  py::tuple fluxes = py::make_tuple(measure.o2mass_in_root,
                                    measure.o2mass_out_root,
                                    measure.o2mass_out_vessels,
                                    measure.o2mass_consumed_transvascular);

  return py::make_tuple(samples, fluxes);
}



DetailedPO2::Parameters* AllocateParametersFromDict(const py::dict &d)
{
  std::auto_ptr<DetailedPO2::Parameters> p(new DetailedPO2::Parameters());
  InitParameters(*p, d);
  return p.release();
}


void TestLinearInterpolation()
{
  DynArray<Eigen::Matrix<float, 3, 1> > data;
  Random rnd(235657);
  float x = 0.;
  for (int i=0;; ++i)
  {
    Eigen::Matrix<float, 3, 1> d; d << x, std::sin(x), std::cos(x);
    data.push_back(d);

    if (i >= 10) break;

    float dx = rnd.Get01f();
    x += dx;
  }


  const float xmax = x;

  DynArray< Eigen::Matrix<float, 4, 1> > samples;

  for (int smpl_cnt=0; smpl_cnt<50; ++smpl_cnt)
  {
    x = rnd.Get01f()*xmax;
    const Eigen::Matrix<float, 3, 1> r = LinearInterpolation<>(x, data);

    Eigen::Matrix<float, 4, 1> s; s << x, r[0], r[1], r[2];
    samples.push_back(s);
  }

  h5cpp::File f("interpolationtest.h5");
  h5cpp::create_dataset(f.root(), "data", h5cpp::Dataspace::simple_dims(data.size(),3), (float*)get_ptr(data));
  h5cpp::create_dataset(f.root(), "samples", h5cpp::Dataspace::simple_dims(samples.size(), 4), (float*)get_ptr(samples));
}




void export_oxygen_computation()
{
  py::class_<Parameters>("DetailedO2Parameters", py::no_init)
    .def("Saturation", &Parameters::Saturation)
    .def("BloodPO2ToConc", &Parameters::BloodPO2ToConc)
    .def("ConcToBloodPO2", &Parameters::ConcToBloodPO2)
    .def("PInit", &Parameters::PInit);
  py::def("AllocateDetailedO2ParametersFromDict", &AllocateParametersFromDict, py::return_value_policy<py::manage_new_object>());
  py::def("computePO2", PyComputePO2);
  py::def("computeSaturation_", PyComputeSaturation);
  py::def("computeConcentration_", PyComputeConcentration);
  py::def("computeMassTransferCoefficient_", PyComputeMassTransferCoefficient);
  py::def("computePO2FromConc_", PyComputePO2FromConc);
  py::def("computeO2Uptake", PyComputeUptake);
  py::def("sampleVessels", PySampleVessels);
  py::def("sampleVesselsWorld", PySampleVesselsWorld);
  // just some tests
  py::def("testCalcOxy", TestSaturationCurve);
  py::def("testCalcOxy2", TestSingleVesselPO2Integration);
}


} // namespace


#ifdef DEBUG
BOOST_PYTHON_MODULE(libdetailedo2_d)
#else
BOOST_PYTHON_MODULE(libdetailedo2_)
#endif
{
  PyEval_InitThreads();
  DetailedPO2::export_oxygen_computation();
}