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

#include "python_helpers.h"

#include "pylatticedata.h"

#include "calcflow.h"
#include "shared-objects.h"
#include "vessels3d.h"


py::list calc_vessel_hydrodynamics(const string fn, const string vesselgroup_path ,bool return_flags, const py::object &py_bfparams, bool simple, bool storeCalculationInHDF)
{
  FpExceptionStateGuard exception_state_guard(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
  
  const BloodFlowParameters bfparams = py::extract<BloodFlowParameters>(py_bfparams);

  H5::H5File readInFile;
  H5::Group g_vess;
  try{
    readInFile = H5::H5File(fn, H5F_ACC_RDONLY);
    g_vess = readInFile.openGroup(vesselgroup_path); // groupname should end by vesselgroup
  }
  catch(H5::Exception &e)
  {
    e.printErrorStack();
  }
  //h5cpp::Group g_vess = h5cpp::Group(readInFile->root().open_group(vesselgroup_path)); // groupname should end by vesselgroup
  
  std::shared_ptr<VesselList3d> vl = ReadVesselList3d(g_vess, make_ptree("filter", false));

  Py_ssize_t num_nodes = vl->GetNCount();
  Py_ssize_t num_edges = vl->GetECount();

  //for some reasons this was needed before calcflow,
  //however for secomb2hdf.py this is bad.
  //ComputeCirculatedComponents(vl.get());
//   if (simple)
//     CalcFlowSimple(*vl, bfparams, true);
//   else
  try{
#ifdef EPETRA_MPI
    std::cout << "EPETRA_MPI flag is set!\n" << std::endl;
    int mpi_is_initialized = 0;
    int prov;
    MPI_Initialized(&mpi_is_initialized);
    if (!mpi_is_initialized)
      //MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE,&prov);
      MPI_Init_thread(0, NULL, 1,&prov);
#endif
  CalcFlow(*vl, bfparams);
  }
  catch(std::exception &ex)
  {
    std::cout << ex.what();
  }
  
  if( storeCalculationInHDF )
  {
    try
    {
      g_vess.openGroup("recomputed");
    }
    catch(H5::Exception &e)
    {
      H5::Group recomp = g_vess.createGroup("recomputed");
      ptree getEverytingPossible = make_ptree("w_adaption", false);
      WriteVesselList3d(*vl, recomp, getEverytingPossible);
    }
//     if( not g_vess.exists("recomputed") )
//     {
//       h5::Group grp_temp;
//       grp_temp = g_vess.create_group("recomputed");
//       ptree getEverytingPossible = make_ptree("w_adaption", false);
//       WriteVesselList3d(*vl, grp_temp, getEverytingPossible);
//     }
  }

#if BOOST_VERSION>106300
// boost::numpy
  py::tuple shape = py::make_tuple(num_edges);
  np::dtype dtype = np::dtype::get_builtin<double>();
  np::ndarray pya_flow = np::empty(shape, dtype);
  np::ndarray pya_force = np::empty(shape, dtype);
  np::ndarray pya_press = np::empty(shape, dtype);

  for (int i=0; i<num_edges; ++i)
  {
    const Vessel* v = vl->GetEdge(i);
    pya_flow[i] = v->q;
    pya_force[i] = v->f;
  }
  for (int i=0; i<num_nodes; ++i)
  {
    pya_press[i] = vl->GetNode(i)->press;
  }
  py::list l;
  l.append(pya_press);
  l.append(pya_flow);
  l.append(pya_force);
  if (!simple)
  {
    np::ndarray pya_hema = np::empty(shape, dtype);
    for (int i=0; i<num_edges; ++i)
    {
      const Vessel* v = vl->GetEdge(i);
      pya_hema[i] = v->hematocrit;
    }
    l.append(pya_hema);
  }
  if (return_flags)
  {
    np::ndarray pya_flags = np::empty(shape, np::dtype::get_builtin<int>());
    for (int i=0; i<num_edges; ++i)
    {
      const Vessel* v = vl->GetEdge(i);
      pya_flags[i] = (int) v->flags;
    }
    l.append(pya_flags);
  }
  
#else
  //welter numpycpp
  np::arrayt<double> pya_flow(np::empty(1, &num_edges, np::getItemtype<double>()));
  np::arrayt<double> pya_force(np::empty(1, &num_edges, np::getItemtype<double>()));
  np::arrayt<double> pya_press(np::empty(1, &num_nodes, np::getItemtype<double>()));

  for (int i=0; i<num_edges; ++i)
  {
    const Vessel* v = vl->GetEdge(i);
    pya_flow(i) = v->q;
    pya_force(i) = v->f;
  }
  for (int i=0; i<num_nodes; ++i)
  {
    pya_press(i) = vl->GetNode(i)->press;
  }
  py::list l;
  l.append(pya_press.getObject());
  l.append(pya_flow.getObject());
  l.append(pya_force.getObject());
  if (!simple)
  {
    np::arrayt<double> pya_hema(np::empty(1, &num_edges, np::getItemtype<double>()));
    for (int i=0; i<num_edges; ++i)
    {
      const Vessel* v = vl->GetEdge(i);
      pya_hema(i) = v->hematocrit;
    }
    l.append(pya_hema.getObject());
  }
  if (return_flags)
  {
    np::arrayt<uint> pya_flags(np::empty(1, &num_edges, np::getItemtype<uint>()));
    for (int i=0; i<num_edges; ++i)
    {
      const Vessel* v = vl->GetEdge(i);
      pya_flags(i) = v->flags;
    }
    l.append(pya_flags.getObject());
  }
#endif
  return l;
};


struct BloodFlowParamsFromPy
{
    static void Register()
    {
      boost::python::converter::registry::push_back(
        &convertible,
        &construct,
        boost::python::type_id<BloodFlowParameters>());
    }

    static void* convertible(PyObject* obj_ptr)
    {
      try
      {
        py::dict o(py::handle<>(py::borrowed(obj_ptr)));
      }
      catch (std::exception &e)
      {
        std::cerr << "Conversion to BloodFlowParameters error: " << e.what();
        return NULL;
      }
      return obj_ptr;
    }

    static void construct(
      PyObject* obj_ptr,
      boost::python::converter::rvalue_from_python_stage1_data* data)
    {
      py::dict o(py::handle<>(py::borrowed(obj_ptr)));
      BloodFlowParameters bfparams; // use default values if no equivalent is found in the dict
      bfparams.rheology = strTo<Rheology>(py::extract<string>(o.get("rheology", toStr(bfparams.rheology))));
      bfparams.viscosityPlasma = py::extract<double>(o.get("viscosityPlasma", (double)bfparams.viscosityPlasma));
      bfparams.inletHematocrit = py::extract<double>(o.get("inletHematocrit", bfparams.inletHematocrit));
      bfparams.includePhaseSeparationEffect = py::extract<bool>(o.get("includePhaseSeparationEffect", bfparams.includePhaseSeparationEffect));

      void* storage = ((boost::python::converter::rvalue_from_python_storage<BloodFlowParameters>*)data)->storage.bytes;
      new (storage) BloodFlowParameters(bfparams);
      data->convertible = storage;
    }
};


struct BloodFlowParamsToPy
{

  static PyObject* convert(const BloodFlowParameters& p)
  {
    py::dict d;
    d["rheology"] = toStr(p.rheology);
    d["viscosityPlasma"] = p.viscosityPlasma;
    d["inletHematocrit"] = p.inletHematocrit;
    d["includePhaseSeparationEffect"] = p.includePhaseSeparationEffect;
    return py::incref(d.ptr());
  }

  static void Register()
  {
    py::to_python_converter<BloodFlowParameters, BloodFlowParamsToPy>();
  }
};


static FlReal PyCalcRelViscosity( FlReal r, FlReal h, string rheologyStr)
{
  Rheology rheology = strTo<Rheology>(rheologyStr);
  return CalcRelViscosity(r, h, rheology);
};


#if BOOST_VERSION>106300
static py::object PyCalcViscosities(np::ndarray pyRad, np::ndarray pyHema, const BloodFlowParameters &bloodFlowParameters)
{
  CheckArray1(FlReal, pyRad, 0); // do not check dimension 
  //np::arrayt<FlReal> rad(pyRad);
  const int ecnt = pyRad.get_shape()[0];
  CheckArray1(FlReal, pyHema, ecnt); // same shape as pyHema
  //np::arrayt<FlReal> hema(pyHema);
  // output
  np::ndarray visc = np::empty(py::tuple(pyRad.get_shape()), np::dtype::get_builtin<FlReal>());
  
  // WARNING: this need cleaning up because it is copy pasted from CalcViscosities
  #pragma omp parallel for
  for(int i=0; i<ecnt; ++i)
  {
    float x = bloodFlowParameters.viscosityPlasma;
    x *= CalcRelViscosity(py::extract<FlReal>(pyRad[i]), py::extract<FlReal>(pyHema[i]), bloodFlowParameters.rheology);
    myAssert(x > 0. && std::isfinite(x));
    visc[i] = x;
  }
  return visc;
}
#else
static py::object PyCalcViscosities(nm::array pyRad, nm::array pyHema, const BloodFlowParameters &bloodFlowParameters)
{
  CheckArray1(FlReal, pyRad, 0); // do not check dimension 
  np::arrayt<FlReal> rad(pyRad);
  const int ecnt = rad.shape()[0];
  CheckArray1(FlReal, pyHema, ecnt); // same shape as pyHema
  np::arrayt<FlReal> hema(pyHema);
  // output
  np::arrayt<FlReal> visc(np::empty(1, rad.shape(), np::getItemtype<FlReal>()));
  
  // WARNING: this need cleaning up because it is copy pasted from CalcViscosities
  #pragma omp parallel for
  for(int i=0; i<ecnt; ++i)
  {
    float x = bloodFlowParameters.viscosityPlasma;
    x *= CalcRelViscosity(rad[i], hema[i], bloodFlowParameters.rheology);
    myAssert(x > 0. && std::isfinite(x));
    visc[i] = x;
  }
  return visc.getObject();
}
#endif


#if BOOST_VERSION>106300
static py::object PyCalcConductivities(np::ndarray pyRad, np::ndarray pyLen, np::ndarray pyVisc)
{
  CheckArray1(FlReal, pyRad, 0); // do not check dimension 
  //np::arrayt<FlReal> rad(pyRad);
  const int ecnt = pyRad.get_shape()[0];
  CheckArray1(FlReal, pyLen, ecnt);
  CheckArray1(FlReal, pyVisc, ecnt);
  //np::arrayt<FlReal> len(pyLen);
  //np::arrayt<FlReal> visc(pyVisc);
  // output
  //np::arrayt<FlReal> cond(np::empty(1, rad.shape(), np::getItemtype<FlReal>()));
  np::ndarray cond = np::empty(py::tuple(pyRad.get_shape()), np::dtype::get_builtin<FlReal>());
  // WARNING: this need cleaning up because it is copy pasted from CalcConductivities
  #pragma omp parallel for
  for(int i=0; i<ecnt; ++i)
  {
    FlReal coeff = CalcFlowCoeff(py::extract<FlReal>(pyVisc[i]),py::extract<FlReal>(pyRad[i]),py::extract<FlReal>(pyLen[i]));
    myAssert(coeff > 0. && std::isfinite(coeff));
    cond[i] = coeff;
  }
  return cond;
}
#else
static py::object PyCalcConductivities(nm::array pyRad, nm::array pyLen, nm::array pyVisc)
{
  CheckArray1(FlReal, pyRad, 0); // do not check dimension 
  np::arrayt<FlReal> rad(pyRad);
  const int ecnt = rad.shape()[0];
  CheckArray1(FlReal, pyLen, ecnt);
  CheckArray1(FlReal, pyVisc, ecnt);
  np::arrayt<FlReal> len(pyLen);
  np::arrayt<FlReal> visc(pyVisc);
  // output
  np::arrayt<FlReal> cond(np::empty(1, rad.shape(), np::getItemtype<FlReal>()));
  
  // WARNING: this need cleaning up because it is copy pasted from CalcConductivities
  #pragma omp parallel for
  for(int i=0; i<ecnt; ++i)
  {
    FlReal coeff = CalcFlowCoeff(visc[i],rad[i],len[i]);
    myAssert(coeff > 0. && std::isfinite(coeff));
    cond[i] = coeff;
  }
  return cond.getObject();
}
#endif


void export_calcflow()
{
  BloodFlowParamsFromPy::Register();
  BloodFlowParamsToPy::Register();
  py::def("calc_vessel_hydrodynamics", calc_vessel_hydrodynamics);
  py::def("CalcFahraeusEffect", CalcFahraeusEffect);
  py::def("CalcRelativeViscosity", PyCalcRelViscosity);
  py::def("PressureRadiusRelation", PressureRadiusRelation);
  py::def("CalcViscosities", PyCalcViscosities);
  py::def("CalcConductivities", PyCalcConductivities);
}

