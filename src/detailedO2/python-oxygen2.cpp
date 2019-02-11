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


#include "H5Cpp.h"
//#include "numpy.hpp"
#include "oxygen_model2.h"
#include "calcflow.h"

#define BOOST_RESULT_OF_USE_DECLTYPE 1

#include <boost/foreach.hpp>
#include <boost/python/errors.hpp>
#include <boost/iterator/transform_iterator.hpp>

/**
 * @brief Sets tissue phases on the lattice sites
 * 
 * Checks wether a tumor calculation is in the h5 file or not.
 * In case it is, the tissue phases are interpolated to lattice points.
 * Otherwise it starts with a homogeneous distribution of normal tissue.
 */

namespace DetailedPO2
{

void InitParameters(DetailedPO2::Parameters &params, const py::dict &py_parameters)
{
#ifndef NDEBUG
  cout << "begin init params" << endl;
#endif

  checkedExtractFromDict(py_parameters, "po2init_r0", params.po2init_r0 );
  checkedExtractFromDict(py_parameters, "po2init_dr", params.po2init_dr );
  checkedExtractFromDict(py_parameters, "po2init_cutoff", params.po2init_cutoff );
  //checkedExtractFromDict(py_parameters, "num_threads", params.num_threads );
  checkedExtractFromDict(py_parameters, "solubility_tissue", params.solubility_tissue );
  checkedExtractFromDict(py_parameters, "solubility_plasma", params.solubility_plasma );
  checkedExtractFromDict(py_parameters, "max_iter", params.max_iter);
  checkedExtractFromDict(py_parameters, "sat_curve_exponent", params.sat_curve_exponent);
  checkedExtractFromDict(py_parameters, "sat_curve_p50", params.sat_curve_p50);
  checkedExtractFromDict(py_parameters, "conductivity_coeff1", params.conductivity_coeff1);
  checkedExtractFromDict(py_parameters, "conductivity_coeff2", params.conductivity_coeff2);
  checkedExtractFromDict(py_parameters, "conductivity_coeff3", params.conductivity_coeff3);
  checkedExtractFromDict(py_parameters, "axial_integration_step_factor", params.axial_integration_step_factor);
  checkedExtractFromDict(py_parameters, "debug_zero_o2field", params.debug_zero_o2field);
  checkedExtractFromDict(py_parameters, "debug_fn", params.debug_fn);
  checkedExtractFromDict(py_parameters, "transvascular_ring_size", params.transvascular_ring_size);
  checkedExtractFromDict(py_parameters, "rd_norm", params.rd_norm);
  checkedExtractFromDict(py_parameters, "rd_tum", params.rd_tum);
  checkedExtractFromDict(py_parameters, "rd_necro", params.rd_necro);
  checkedExtractFromDict(py_parameters, "michaelis_menten_uptake",params.michaelis_menten_uptake);
  
  checkedExtractFromDict(py_parameters, "haemoglobin_binding_capacity",params.haemoglobin_binding_capacity);
  checkedExtractFromDict(py_parameters, "tissue_boundary_value",params.tissue_boundary_value);
  checkedExtractFromDict(py_parameters, "extra_tissue_source_linear",params.extra_tissue_source_linear);
  checkedExtractFromDict(py_parameters, "extra_tissue_source_const",params.extra_tissue_source_const);
  checkedExtractFromDict(py_parameters, "convergence_tolerance",params.convergence_tolerance);
  checkedExtractFromDict(py_parameters, "loglevel",params.loglevel);
  checkedExtractFromDict(py_parameters, "approximateInsignificantTransvascularFlux",params.approximateInsignificantTransvascularFlux);
  checkedExtractFromDict(py_parameters, "D_plasma",params.D_plasma);
  checkedExtractFromDict(py_parameters, "massTransferCoefficientModelNumber",params.massTransferCoefficientModelNumber);
  checkedExtractFromDict(py_parameters, "tissue_po2_boundary_condition",params.tissue_po2_boundary_condition);
  
  checkedExtractFromDict(py_parameters, "grid_lattice_const",params.grid_lattice_const);
  checkedExtractFromDict(py_parameters, "safety_layer_size",params.safety_layer_size);
 
  
  // required parameters for transvascular transport
  if (params.tissue_po2_boundary_condition == "dirichlet_x") 
    params.tissue_boundary_condition_flags = 1; //FiniteVolumeMatrixBuilder;
  else if (params.tissue_po2_boundary_condition == "dirichlet_yz") 
    params.tissue_boundary_condition_flags = 2;
  else if (params.tissue_po2_boundary_condition == "dirichlet") 
    params.tissue_boundary_condition_flags = 3;
  else if (params.tissue_po2_boundary_condition == "neumann") 
    params.tissue_boundary_condition_flags = 0;
  else
  {
    std::cout << "bc_type: " << params.tissue_po2_boundary_condition << std::endl;
    throw std::invalid_argument("tissue_po2_boundary_condition must be 'dirichlet','dirichlet_x', 'dirichlet_yz' or 'neumann'");
  }
  
  
  //calculate stuff parameters dependent on the given parameters
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
//static void PyComputePO2(py::object py_vesselgroup, py::object py_tumorgroup, py::dict py_parameters, py::object py_bfparams, py::object py_h5outputGroup)
//static void PyComputePO2(string fn, string vesselgroup_path, string tumorgroup_path, py::dict py_parameters, py::object py_bfparams, string h5_out_path)
static void PyComputePO2(py::dict &py_parameters, py::object &py_bfparams)
{
  DetailedPO2::DetailedPO2Sim s;
  s.params = py::extract<DetailedPO2::Parameters>(py_parameters);
  //calculate stuff parameters dependent on the given parameters
  s.params.UpdateInternalValues();
  
  
  
  //Parameters params;
  //InitParameters(s.params, py_parameters);
  cout << "parameters initialized" << std::endl;
  
  //h5cpp::Group vesselgroup = PythonToCppGroup(py_vesselgroup);
  // this was fn
  s.params.input_file_name = py::extract<string>(py_parameters.get("input_file_name", "none"));
  s.params.input_group_path = py::extract<string>(py_parameters.get("input_group_path", "none")); 
  s.params.vessel_group_path = py::extract<string>(py_parameters.get("vessel_group_path", "none"));
  const string output_file_name = py::extract<string>(py_parameters.get("output_file_name", "none"));
  const string tumor_file_name = py::extract<string>(py_parameters.get("tumor_file_name", "none"));
  const string tumor_group_path = py::extract<string>(py_parameters.get("tumor_group_path", "none"));
  const string out_grp_path = py::extract<string>(py_parameters.get("output_group_path", "none"));
  string vesselgroup_path = py::extract<string>(py_parameters.get("vessel_group_path", "none"));
  //vesselgroup_path = "/" + vesselgroup_path;
  //input_group_path = "/" + input_group_path;
  //h5cpp::File *o2File = new h5cpp::File(fn,"a");
  //std::shared_ptr<VesselList3d> vl;
  //s.params.input_group_path = input_group_path;
  H5::H5File o2File;
  H5::H5File vesselInputFile;
  try
  {
    vesselInputFile = H5::H5File( s.params.input_file_name , H5F_ACC_RDONLY);
    H5::Group vesselgroup = vesselInputFile.openGroup(string("/") + s.params.input_group_path);
    s.vl = ReadVesselList3d(vesselgroup, make_ptree("filter",false));
    o2File = H5::H5File( output_file_name, H5F_ACC_TRUNC );
  }  
  catch( H5::FileIException &error )
  {
    error.printErrorStack();
  }
    
  // THIIIIRYYYYY, filter muss = false sein sonst stimmt in der Ausgabe in der Hdf5 Datei die Anzahl der Vessels nicht mehr mit den daten im recomputed_flow Verzeichnis ueberein!
  
  //double grid_lattice_const               = py::extract<double>(py_parameters.get("grid_lattice_const", 30.));
  //double safety_layer_size                = py::extract<double>(py_parameters.get("safety_layer_size", grid_lattice_const*3.));
  //boost::optional<Int3> grid_lattice_size = getOptional<Int3>("grid_lattice_size", py_parameters);
  

  //if (!py_bfparams.is_none())
  BloodFlowParameters bfparams;
  try
  {
    bfparams = py::extract<BloodFlowParameters>(py_bfparams);
  }
  catch(H5::Exception &e)
  {
    e.printErrorStack();
    cerr << "could not extract blood flow parameters from python... using default values from constructor" << endl;
  }
  
  CalcFlow(*s.vl, bfparams);
  
  try //create recomputed_flow in hdf
  {
    H5::Group vess_recomp = H5::Group(o2File.createGroup("recomputed_flow").createGroup("vessels")); // groupname should end by vesselgroup
    ptree getEverytingPossible = make_ptree("w_adaption", false);
    WriteVesselList3d(*s.vl, vess_recomp, getEverytingPossible);
    vess_recomp.close();
  }
  catch(H5::Exception &e)
  {
    e.printErrorStack();
    cerr << "could not create vess_recomp " << endl;
  }
  
  
  
  
  //cout << format("in c++: %.20f %.20f %.20f\n") % params.conductivity_coeff1 % params.conductivity_coeff2 % params.conductivity_coeff_gamma;
  boost::optional<H5::H5File> h5_tumor_file;
  boost::optional<H5::Group> tumorgroup;
  boost::optional<Array3df> previous_po2field;
  boost::optional<DetailedPO2::VesselPO2Storage> previous_po2vessels;
  boost::optional<Array3d<float>> cell_based_o2_uptake;
  
  //boost::optional<std::string> tumorgroup_path;
  //h5cpp::Group   *tumorgroup = new h5cpp::Group();
  
//   if (!py_tumorgroup.is_none())
//   {
//     tumorgroup = PythonToCppGroup(py_tumorgroup);
//   }
  if (tumor_group_path != "none")
  {
    h5_tumor_file = boost::optional<H5::H5File>(H5::H5File(tumor_file_name, H5F_ACC_RDONLY));
    tumorgroup = boost::optional<H5::Group>(h5_tumor_file->openGroup(tumor_group_path));
  }
  
//  s.init(params, bfparams,*vl,grid_lattice_const, safety_layer_size, grid_lattice_size, tumorgroup, previous_po2field, previous_po2vessels);
//   else
//   {
//     tumorgroup = nullptr;
//   }
  
  
  
  //{ //ReleaseGIL unlock(); // allow the python interpreter to do things while this is running
  /**
   * MOST IMPORTANT CALL
   * grid is conitnuum grid, mtboxes is decomposition into threads
   */
  try
  {
#ifdef EPETRA_MPI
    std::cout << "EPETRA_MPI flag is set!\n" << std::endl;
    int mpi_is_initialized = 0;
    int prov;
    MPI_Initialized(&mpi_is_initialized);
    if (!mpi_is_initialized)
      //MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE,&prov);
      MPI_Init_thread(0, NULL, 1,&prov);
#endif
    //DetailedP02Sim s;
    s.init(bfparams, tumorgroup, previous_po2field, previous_po2vessels,cell_based_o2_uptake);
    s.run();
  }
  catch(std::exception &ex)
  {
    std::cout << ex.what();
  }
  catch(H5::Exception &e)
  {
    e.printErrorStack();
  }
    
    /* OUTPUT */
    
  H5::Group po2_out_group = o2File.createGroup(string("/po2"));
  H5::Group h5_out_params = o2File.createGroup(string("/parameters"));
  s.WriteDataOutput(po2_out_group);
  s.WriteParametersToHDF(h5_out_params);
  //s.WriteOutput_new(po2_out_group, h5_out_params);
    
//     //H5::Group outputgroup = o2File.createGroup(out_grp_path);
//     H5::Group po2_out_group = o2File.createGroup(string("/po2"));
//     H5::Group outputgroup = o2File.createGroup(string("/") + string("po2/") + input_group_path);
//     
//     writeAttrToH5(outputgroup, string("SOURCE_VESSELS_FILE"), input_file_name);
//     writeAttrToH5(outputgroup, string("SOURCE_VESSELS_PATH"), input_group_path);
//     writeAttrToH5(outputgroup, string("SOURCE_TISSUE_FILE"), string("none"));
//     writeAttrToH5(outputgroup, string("SOURCE_TISSUE_PATH"), string("none"));
//     H5::Group h5_o2_lattice = outputgroup.createGroup("field_ld");
//     s.grid.ld.WriteHdfLd(h5_o2_lattice);
//     WriteScalarField(outputgroup, string("po2field"), s.po2field, s.grid.ld, h5_o2_lattice);
//     writeDataSetToGroup(outputgroup, string("po2vessels"), s.po2vessels);
//     //WriteHdfPtree(outputgroup, s.metadata, HDF_WRITE_PTREE_AS_DATASETS);
//     WriteHdfPtree(outputgroup, s.metadata, HDF_WRITE_PTREE_AS_ATTRIBUTE);
//     H5::Group h5_params = outputgroup.createGroup("parameters");
//     H5::Group h5_meta_data = h5_params.createGroup("metadata");
//     WriteHdfPtree(h5_meta_data, s.metadata, HDF_WRITE_PTREE_AS_ATTRIBUTE);
//     H5::Group h5_o2_params = h5_params.createGroup("o2");
//     WriteHdfPtree(h5_o2_params, s.params.as_ptree());
//     //s.params.writeParametersToHDF(h5_o2_params);
//     //WriteHdfPtree(h5_o2_params, s.params.as_ptree(), HDF_WRITE_PTREE_AS_ATTRIBUTE);
//     H5::Group h5_bf_params = h5_params.createGroup("calcflow");
//     WriteHdfPtree(h5_bf_params, s.bfparams.as_ptree());
//   }
//   catch(std::exception &ex)
//   {
//     std::cout << ex.what();
//   }
//   catch(H5::Exception e)
//   {
//     e.printErrorStack();
//   }
  
  //cout << "c++ part of detailedO2 finished" << std::endl;
  //DetailedPO2::ComputePO2(params, *vl, grid, mtboxes, po2field, po2vessels, phases, metadata, world);
//   }
//   if (PyErr_Occurred() != NULL) return; // don't save stuff
  
  //This writes the results back to python
#if 0
  Int2 sz(2,po2vessels.size());
  nm::array  py_po2vessels = np::copy<float,2>(sz.data(), (float*)po2vessels.data(), calc_strides::first_dim_varies_fastest(sz).data());
  nm::array  py_po2field   = np::copy<float,3>(po2field.size().data(), po2field.getPtr(), po2field.strides().data());
  py::object py_ld(grid.ld);
  
  return py::make_tuple(py_po2vessels, py_ld, py_po2field);
#endif
  
//    readInFile->close();
//    h5cpp::File *outputFile = new h5cpp::File(fn,"a");
//    try{
//     H5::Group outputgroup = o2File.createGroup();
//     H5::Group h5_o2_lattice = outputgroup.createGroup("field_ld");
//     s.grid.ld.WriteHdfLd(h5_o2_lattice);
//     WriteScalarField(outputgroup, "po2field", s.po2field, s.grid.ld, ldgroup);
//     WriteHdfPtree(outputgroup, s.metadata, HDF_WRITE_PTREE_AS_DATASETS);
//     //WriteHdfLd(h5_o2_lattice, s.
//    }
//    catch(H5::Exception e)
//    {
//      e.printErrorStack();
//    }

//     h5cpp::Group root = o2File->root();
//     h5cpp::Group gout = root.create_group("data");
//     h5cpp::Group g_o2;
//     h5cpp::Group po2outputGroup = gout.create_group("po2");
//     h5cpp::Group ldgroup = po2outputGroup.create_group("field_ld");
//     cout<<"start writiong!"<<endl;
    //WriteHdfLd(ldgroup, s.grid.ld);
    //o2File->flush();
    //outputFile.flush();
    //outputFile.close();
    //h5cpp::Group outputGroup = outputFile.root().create_group("data");
    //h5cpp::Group outputGroup = h5cpp::Group(outputFile.root().create_group("data")); // groupname should end by vesselgroup
// //     h5cpp::Group outputGroup = PythonToCppGroup(py_h5outputGroup);
    //h5cpp::Group ldgroup = po2outputGroup.create_group("field_ld");
    //WriteHdfLd(ldgroup, s.grid.ld);
    //WriteScalarField<float>(outputGroup, "po2field", s.po2field, s.grid.ld, ldgroup);
//      h5cpp::create_dataset<float>(outputGroup, "po2vessels", h5cpp::Dataspace::simple_dims(s.po2vessels.size(), 2), (float*)s.po2vessels[0].data(), h5cpp::CREATE_DS_COMPRESSED); // FIX ME: transpose the array!
//      WriteHdfPtree(outputGroup, s.metadata, HDF_WRITE_PTREE_AS_DATASETS);
    //outputFile.close();
    std::cout << boost::format("wrote o2 to file %s at %s") % output_file_name % out_grp_path;
}

#if BOOST_VERSION>106300
np::ndarray PyComputeSaturation(const np::ndarray &py_po2, const py::dict &py_parameters)
{
  cout << "PyComputeSaturation called" << endl;
  cout << "BOOST_VERSION>106300 found" << endl;
  /**
   * this needs to behave differently wether py::dict everyting is 
   * in strings or with correct data type, which is new 
   */
  //std::string D_plasma_cpp_string = py::extract<std::string>(py_parameters["vessel_group_path"]);
  //double D_plasma_cpp_double = py::extract<double>(py_parameters["D_plasma"]);
  DetailedPO2::Parameters params = py::extract<DetailedPO2::Parameters>(py_parameters);
  const int return_length = py_po2.get_shape()[0];
  cout << "shape of py_po2: " << endl
  << py_po2.get_shape()[0] << endl 
  <<py_po2.get_shape()[1] << endl;
//   if(!(py_po2.get_nd() == 1))
//     throw std::invalid_argument("rank 1 and contiguous expected");
#ifndef NDEBUG
  cout << "make empty in PyComputeSaturation: " << return_length << endl;
#endif
  cout << "return_length: " << return_length << endl;
  //py::tuple shape = py::make_tuple(3, 1);
  //np::dtype dtype = np::dtype::get_builtin<double>();
  //np::ndarray a = np::zeros(shape, dtype);
  np::ndarray result = np::zeros(py::make_tuple(return_length,1), np::dtype::get_builtin<float>());
  //np::ndarray result = np::zeros(py::make_tuple(py_po2.shape(0)), np::dtype::get_builtin<float>());
  cout << "bla" << endl;cout.flush();
  for (int i=0; i<return_length; i++)
  {
    //data[i] = params.Saturation(py::extract<float>(py_po2[i]));
    auto bla = params.Saturation(py::extract<float>(py_po2[i]));
    result[i] = bla;
  }
  cout << "returning from PyComputeSaturation" << endl;
  return result;
}
#else
static py::object PyComputeSaturation(nm::array &py_po2, py::dict &py_parameters)
{
  cout << "PyComputeSaturation called" << endl;
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
#endif


#if BOOST_VERSION>106300
#else
static py::object PyComputeMassTransferCoefficient(nm::array py_radius, py::dict py_parameters)
{
  cout << "PyComputeMassTransferCoefficient called" << endl;
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
#endif


#if BOOST_VERSION>106300
static py::object PyComputeConcentration(np::ndarray py_po2, np::ndarray py_hematocrit, py::dict py_parameters)
{
  cout << "PyComputeConcentration called" << endl;
  DetailedPO2::Parameters params;
  InitParameters(params, py_parameters);
  
//   np::arrayt<float> po2(py_po2);
//   np::arrayt<float> hema(py_hematocrit);

//   if (!(po2.rank() == 1 && po2.isCContiguous())) throw std::invalid_argument("rank 1 and contiguous expected");
//   if (!(hema.rank() == 1 && hema.isCContiguous())) throw std::invalid_argument("rank 1 and contiguous expected");
  if (!(py_po2.get_nd() == 1 )) throw std::invalid_argument("1 d array expected");
  if (!(py_hematocrit.get_nd() == 1)) throw std::invalid_argument("1 d array expected");

//   np::arrayt<float> result(np::empty(1, po2.shape(), np::getItemtype<float>()));
  np::ndarray result = np::empty(py::make_tuple(py_po2.shape(0)), np::dtype::get_builtin<float>());

  for (int i=0; i<py_po2.shape(0); ++i)
  {
    result[i] = params.BloodPO2ToConc(py::extract<float>(py_po2[i]), py::extract<float>(py_hematocrit[i]));
  }
  return result;
}
#else
static py::object PyComputeConcentration(nm::array py_po2, nm::array py_hematocrit, py::dict py_parameters)
{
  cout << "PyComputeConcentration called" << endl;
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
#endif


#if BOOST_VERSION>106300
#else
static py::object PyComputePO2FromConc(nm::array py_conc, nm::array py_hematocrit, py::dict py_parameters)
{
  cout << "PyComputePO2FromConc called" << endl;
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
#endif


#if 0// not H5Cpp ready
#if BOOST_VERSION>106300
#else
static py::object PyComputeUptake(nm::array py_po2field, const LatticeDataQuad3d &field_ld, py::object py_tumorgroup,  py::dict py_parameters)
{
  DetailedPO2::Parameters params;
  InitParameters(params, py_parameters);
  ContinuumGrid grid;
  grid.init(field_ld, field_ld.Size()[2]<=1 ? 2 : 3);
  DomainDecomposition mtboxes;
  mtboxes.init(MakeMtBoxGrid(grid.Box(), Int3(32, 32, 32)));
  TissuePhases phases;
  /* this needs to to properly done later*/
  h5cpp::Group tumorgroup;
  if (!py_tumorgroup.is_none())
  {
    tumorgroup = PythonToCppGroup(py_tumorgroup);
  }
  
  SetupTissuePhases(phases, grid, mtboxes, tumorgroup);
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
#endif
#endif

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

#if 0// not H5Cpp ready
#if BOOST_VERSION>106300
#else
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
  
  TissuePhases phases;//Declaration
  h5cpp::Group tumorgroup;
  if (!py_tumorgroup.is_none())
  {
    tumorgroup = PythonToCppGroup(py_tumorgroup);
  }
  SetupTissuePhases(phases, grid, mtboxes, tumorgroup);//filling
  
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
#endif
#endif

#if 0//not H5Cpp ready
#if BOOST_VERSION>106300
#else
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
  
  TissuePhases phases;//Declaration
  
  h5cpp::Group tumorgroup;
  if (!py_tumorgroup.is_none())
  {
    tumorgroup = PythonToCppGroup(py_tumorgroup);
  }
  SetupTissuePhases(phases, grid, mtboxes, tumorgroup);//filling
  
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
#endif
#endif

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

//   h5cpp::File f("interpolationtest.h5");
//   h5cpp::create_dataset(f.root(), "data", h5cpp::Dataspace::simple_dims(data.size(),3), (float*)get_ptr(data));
//   h5cpp::create_dataset(f.root(), "samples", h5cpp::Dataspace::simple_dims(samples.size(), 4), (float*)get_ptr(samples));
}


struct DetailedO2ParamsFromPy
{
    static void Register()
    {
      boost::python::converter::registry::push_back(
        &convertible,
        &construct,
        boost::python::type_id<Parameters>());
    }

    static void* convertible(PyObject* obj_ptr)
    {
      try
      {
        py::dict o(py::handle<>(py::borrowed(obj_ptr)));
      }
      catch (std::exception &e)
      {
        std::cerr << "Conversion to DetailedPO2::Parameters error: " << e.what();
        return NULL;
      }
      return obj_ptr;
    }

    static void construct(
      PyObject* obj_ptr,
      boost::python::converter::rvalue_from_python_stage1_data* data)
    {
      py::dict o(py::handle<>(py::borrowed(obj_ptr)));
      Parameters o2params; // use default values if no equivalent is found
      /* lets check whether D_plasma is stored as double or as
       * string 
       */
      bool isEverythingInDictAString;
      try
      {
        std::string buffer = py::extract<std::string>(o["D_plasma"]);
        isEverythingInDictAString = true;
      }
      catch(std::exception &e)
      {
        isEverythingInDictAString = false;
        e.what();
      }
      if (isEverythingInDictAString)
      {
        cout << "Your O2 parameters are stored as strings." << endl;
        #define DOPT(name) checkedExtractFromDict_all_strings(o, #name, o2params.name)
        DOPT(po2init_r0);
        DOPT(po2init_dr);
        DOPT(po2init_cutoff);
        DOPT(solubility_plasma); 
        DOPT(sat_curve_exponent);
        DOPT(sat_curve_p50);
        DOPT(po2_mmcons_k_norm);
        DOPT(po2_mmcons_k_tum);
        DOPT(po2_mmcons_k_necro);
        DOPT(po2_mmcons_m0_norm);
        DOPT(po2_mmcons_m0_tum);
        DOPT(po2_mmcons_m0_necro);
        
        DOPT(D_plasma);
        DOPT(solubility_tissue);
        DOPT(rd_norm);
        DOPT(rd_tum);
        DOPT(rd_necro);
        DOPT(max_iter);
        DOPT(num_threads);
        DOPT(convergence_tolerance);
        DOPT(axial_integration_step_factor);
        DOPT(debug_zero_o2field);
        
        DOPT(michaelis_menten_uptake);
        DOPT(useCellBasedUptake);
        DOPT(massTransferCoefficientModelNumber);
        DOPT(conductivity_coeff1);
        DOPT(conductivity_coeff2);
        DOPT(conductivity_coeff3);
        DOPT(detailedO2name);
        DOPT(loglevel);
        DOPT(tissue_po2_boundary_condition);
        DOPT(extra_tissue_source_const);
        DOPT(extra_tissue_source_linear);
        
        DOPT(input_file_name);
        DOPT(input_group_path);
        DOPT(output_file_name);
        DOPT(output_group_path);
        DOPT(tumor_file_name);
        DOPT(tumor_group_path);
        DOPT(vessel_group_path);
        
        DOPT(debug_fn);
        DOPT(transvascular_ring_size);
        DOPT(haemoglobin_binding_capacity);
        DOPT(tissue_boundary_value);
        DOPT(approximateInsignificantTransvascularFlux);
        
        DOPT(grid_lattice_const);
        DOPT(safety_layer_size);
        
        #undef DOPT
      }
      if (not isEverythingInDictAString)
      {
        //for debug PURPOSE only
        //checkedExtractFromDict_all_strings(o, "po2init_r0", o2params.po2init_r0);
        cout << "Your O2 parameters are with appropriate data types." << endl;
        #define DOPT(name) checkedExtractFromDict(o, #name, o2params.name)
        DOPT(po2init_r0);
        DOPT(po2init_dr);
        DOPT(po2init_cutoff);
        DOPT(solubility_plasma); 
        DOPT(sat_curve_exponent);
        DOPT(sat_curve_p50);
        DOPT(po2_mmcons_k_norm);
        DOPT(po2_mmcons_k_tum);
        DOPT(po2_mmcons_k_necro);
        DOPT(po2_mmcons_m0_norm);
        DOPT(po2_mmcons_m0_tum);
        DOPT(po2_mmcons_m0_necro);
        
        DOPT(D_plasma);
        DOPT(solubility_tissue);
        DOPT(rd_norm);
        DOPT(rd_tum);
        DOPT(rd_necro);
        DOPT(max_iter);
        DOPT(num_threads);
        DOPT(convergence_tolerance);
        DOPT(axial_integration_step_factor);
        DOPT(debug_zero_o2field);
        
        DOPT(michaelis_menten_uptake);
        DOPT(useCellBasedUptake);
        DOPT(massTransferCoefficientModelNumber);
        DOPT(conductivity_coeff1);
        DOPT(conductivity_coeff2);
        DOPT(conductivity_coeff3);
        DOPT(detailedO2name);
        DOPT(loglevel);
        DOPT(tissue_po2_boundary_condition);
        DOPT(extra_tissue_source_const);
        DOPT(extra_tissue_source_linear);
        
        DOPT(input_file_name);
        DOPT(input_group_path);
        DOPT(output_file_name);
        DOPT(output_group_path);
        DOPT(tumor_file_name);
        DOPT(tumor_group_path);
        DOPT(vessel_group_path);
        
        DOPT(debug_fn);
        DOPT(transvascular_ring_size);
        DOPT(haemoglobin_binding_capacity);
        DOPT(tissue_boundary_value);
        DOPT(approximateInsignificantTransvascularFlux);
        
        DOPT(grid_lattice_const);
        DOPT(safety_layer_size);
        
        #undef DOPT
      }

      void* storage = ((boost::python::converter::rvalue_from_python_storage<Parameters>*)data)->storage.bytes;
      new (storage) Parameters(o2params);
      data->convertible = storage;
    }
};


struct DetailedO2ParamsToPy
{

  static PyObject* convert(const Parameters& p)
  {
    #define DOPT(name) d[#name] = p.name
    py::dict d;
    DOPT(po2init_r0);
    DOPT(po2init_dr);
    DOPT(po2init_cutoff);
    DOPT(solubility_plasma); 
    DOPT(sat_curve_exponent);
    DOPT(sat_curve_p50);
    DOPT(po2_mmcons_k_norm);
    DOPT(po2_mmcons_k_tum);
    DOPT(po2_mmcons_k_necro);
    DOPT(po2_mmcons_m0_norm);
    DOPT(po2_mmcons_m0_tum);
    DOPT(po2_mmcons_m0_necro);
    
    DOPT(D_plasma);
    DOPT(solubility_tissue);
    DOPT(rd_norm);
    DOPT(rd_tum);
    DOPT(rd_necro);
    DOPT(max_iter);
    DOPT(num_threads);
    DOPT(convergence_tolerance);
    DOPT(axial_integration_step_factor);
    DOPT(debug_zero_o2field);
    
    DOPT(michaelis_menten_uptake);
    DOPT(useCellBasedUptake);
    DOPT(massTransferCoefficientModelNumber);
    DOPT(conductivity_coeff1);
    DOPT(conductivity_coeff2);
    DOPT(conductivity_coeff3);
    DOPT(detailedO2name);
    DOPT(loglevel);
    DOPT(tissue_po2_boundary_condition);
    DOPT(extra_tissue_source_const);
    DOPT(extra_tissue_source_linear);
    
    DOPT(input_file_name);
    DOPT(input_group_path);
    DOPT(output_file_name);
    DOPT(output_group_path);
    DOPT(tumor_file_name);
    DOPT(tumor_group_path);
    DOPT(vessel_group_path);
    
    DOPT(debug_fn);
    DOPT(transvascular_ring_size);
    DOPT(haemoglobin_binding_capacity);
    DOPT(tissue_boundary_value);
    DOPT(approximateInsignificantTransvascularFlux);
    
    DOPT(grid_lattice_const);
    DOPT(safety_layer_size);
    
    #undef DOPT
    return py::incref(d.ptr());
  }

  static void Register()
  {
    py::to_python_converter<Parameters, DetailedO2ParamsToPy>();
  }
};

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
  //those are not yet boost::numpy ready
//  py::def("computeMassTransferCoefficient_", PyComputeMassTransferCoefficient);
//   py::def("computePO2FromConc_", PyComputePO2FromConc);
//   py::def("computeO2Uptake", PyComputeUptake);
//   py::def("sampleVessels", PySampleVessels);
//   py::def("sampleVesselsWorld", PySampleVesselsWorld);
  // just some tests
  py::def("testCalcOxy", TestSaturationCurve);
  py::def("testCalcOxy2", TestSingleVesselPO2Integration);
}

void export_parameter_converters()
{
  DetailedO2ParamsFromPy::Register();
  DetailedO2ParamsToPy::Register();
}
} // namespace


#ifndef NDEBUG
BOOST_PYTHON_MODULE(libdetailedo2_d)
#else
BOOST_PYTHON_MODULE(libdetailedo2_)
#endif
{
  PyEval_InitThreads();
  my::checkAbort = PyCheckAbort; // since this is the python module, this is set to use the python signal check function
  DetailedPO2::export_oxygen_computation();
  DetailedPO2::export_parameter_converters();
  //bla
}
