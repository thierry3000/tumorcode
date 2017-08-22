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
#include "../../src/python_krebsutils/python_helpers.h"
#include "numpy.hpp"
#include "adaption_model2.h"

#ifdef USE_PAGMO
#include "adaption_as_pagmo_problem.h"
#endif

#include "../common/calcflow.h"
#include <algorithm>

#define BOOST_RESULT_OF_USE_DECLTYPE 1

#include <boost/foreach.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>

namespace py = boost::python;
namespace nm = boost::python::numeric;

/**
 * @brief Calculates radii due to metabolic and topological demand.
 * 
 * Paper by Secomb et al.
 */

namespace Adaption
{
//bool IsTrue(bool pbool){return pbool;} 

void InitBFParameters(BloodFlowParameters *params, const py::dict *py_parameters)
{
#define GET_ADAPTION_PARAM_FROM_DICT(TYPE, NAME) py::extract<TYPE>(py_parameters->get(NAME))
#define GET_ADAPTION_PARAM_IF_NONNONE(TARGET, TYPE, NAME) { py::object o(py_parameters->get(NAME)); if (!o.is_none()) TARGET=py::extract<TYPE>(o); }
  
#ifdef DEBUG
  printf("entered InitBFParameters\n");
#endif
  //this could possibly done for every parameter
  try
  {
    double viscosityPlasma = GET_ADAPTION_PARAM_FROM_DICT(double,"viscosityPlasma");
    if(viscosityPlasma<0)
      throw std::runtime_error("got bad viscosity from python2");
    params->viscosityPlasma = viscosityPlasma;
    double inletHematocrit = GET_ADAPTION_PARAM_FROM_DICT(double,"inletHematocrit");
    if(inletHematocrit <=0 or inletHematocrit >=1)
      throw std::runtime_error("got bad inletHematocrit from python2");
    bool includePhaseSeparationEffect = GET_ADAPTION_PARAM_FROM_DICT(bool,"includePhaseSeparationEffect");
#ifdef DEBUG
      cout << "got as bool " << includePhaseSeparationEffect << endl;
#endif
    if(includePhaseSeparationEffect == true or includePhaseSeparationEffect == false)
      params->includePhaseSeparationEffect = includePhaseSeparationEffect;
    string rheology = GET_ADAPTION_PARAM_FROM_DICT(string,"rheology");
#ifdef DEBUG
    cout << "got rheology string" << rheology << endl;
#endif
    params->rheology = strTo<Rheology>(rheology);
    //bfparams.rheology = strTo<Rheology>(py::extract<string>(o.get("rheology", toStr(bfparams.rheology))));
    //cout << std::numeric_limits<Rheology>::max() << endl;
    //cout << "end" << toStr(params->rheology) << endl;
    bool foundRheology=false;
    for(int i=0;i<RheologyEnd; i++)
    {
      if(toStr(Rheology(i)) == rheology)
      {
        foundRheology = true;
      }
    }
    if(not foundRheology)
      throw std::runtime_error("got bad rheology from python2");
  }
  catch(std::runtime_error &e)
  {
    e.what();
  }

#undef GET_ADAPTION_PARAM_FROM_DICT
#undef GET_ADAPTION_PARAM_IF_NONNONE
#ifdef DEBUG
  printf("leave InitBFParameters\n");
#endif
}

void InitParameters(Adaption::Parameters *params, const py::dict *py_parameters)
{
#define GET_ADAPTION_PARAM_FROM_DICT(TYPE, NAME) py::extract<TYPE>(py_parameters->get(NAME))
#define GET_ADAPTION_PARAM_IF_NONNONE(TARGET, TYPE, NAME) { py::object o(py_parameters->get(NAME)); if (!o.is_none()) TARGET=py::extract<TYPE>(o); }
  
#ifdef DEBUG
  printf("entered InitParameters\n");
#endif
  //this could possibly done for every parameter
  try
  {
    double k_c = GET_ADAPTION_PARAM_FROM_DICT(double,"k_c");
    if(k_c<0)
      throw std::runtime_error("got bad k_c from python2");
    params->k_c = k_c;
  }
  catch(std::runtime_error &e)
  {
    e.what();
  }
  params->k_m = GET_ADAPTION_PARAM_FROM_DICT(double,"k_m");
  params->k_s = GET_ADAPTION_PARAM_FROM_DICT(double,"k_s");
  params->Q_refdot = GET_ADAPTION_PARAM_FROM_DICT(double,"Q_refdot");
  params->S_0 = GET_ADAPTION_PARAM_FROM_DICT(double,"S_0");
  params->max_nun_iterations = GET_ADAPTION_PARAM_FROM_DICT(double,"max_nun_iterations");
  params->qdev = GET_ADAPTION_PARAM_FROM_DICT(double,"qdev");
  params->starting_radii = GET_ADAPTION_PARAM_FROM_DICT(double,"starting_radii");
  params->delta_t = GET_ADAPTION_PARAM_FROM_DICT(double,"delta_t");
  params->cond_length = GET_ADAPTION_PARAM_FROM_DICT(double,"cond_length");
  GET_ADAPTION_PARAM_IF_NONNONE(params->avgRootNodeConductivity, double, "avgRootNodeConductivity");
  GET_ADAPTION_PARAM_IF_NONNONE(params->radMin_for_kill, double, "radMin_for_kill");
  GET_ADAPTION_PARAM_IF_NONNONE(params->boundary_Condition_handling, uint, "boundary_Condition_handling");
  GET_ADAPTION_PARAM_IF_NONNONE(params->a_pressure, double, "a_pressure");
  GET_ADAPTION_PARAM_IF_NONNONE(params->a_flow, double, "a_flow");
  //std::cout<<params->write2File<<std::endl;
  GET_ADAPTION_PARAM_IF_NONNONE(params->write2File,bool,"write2File")
  GET_ADAPTION_PARAM_IF_NONNONE(params->outputFileName,string,"outputFileName")
  //std::cout<<params->write2File<<std::endl;
  //std::cout<<params->tum_manitulate_s1<<std::endl;
  params->tum_manitulate_s1 = GET_ADAPTION_PARAM_FROM_DICT(bool,"tum_manitulate_s1");
  //std::cout<<params->tum_manitulate_s1<<std::endl;
  params->tum_manitulate_s2 = GET_ADAPTION_PARAM_FROM_DICT(bool,"tum_manitulate_s2");
  params->tum_manitulate_s3 = GET_ADAPTION_PARAM_FROM_DICT(bool,"tum_manitulate_s3");
  params->tum_manitulate_s4 = GET_ADAPTION_PARAM_FROM_DICT(bool,"tum_manitulate_s4");
  params->tum_manitulate_s5 = GET_ADAPTION_PARAM_FROM_DICT(bool,"tum_manitulate_s5");
  
  params->vesselFileName = GET_ADAPTION_PARAM_FROM_DICT(string,"vesselFileName");
  params->vesselGroupName = GET_ADAPTION_PARAM_FROM_DICT(string,"vesselGroupName");
  
  GET_ADAPTION_PARAM_IF_NONNONE(params->pop, int, "pop");
  GET_ADAPTION_PARAM_IF_NONNONE(params->individuals, int, "individuals");
  GET_ADAPTION_PARAM_IF_NONNONE(params->opt_iter, int, "opt_iter");
  
  //GET_ADAPTION_PARAM_IF_NONNONE(params->pop, int, "pop");
  //GET_ADAPTION_PARAM_IF_NONNONE(params->individuals, int, "individuals");
#undef GET_ADAPTION_PARAM_FROM_DICT
#undef GET_ADAPTION_PARAM_IF_NONNONE
#ifdef DEBUG
  printf("leave InitParameters\n");
#endif
}

/*
 * we will consider the circulated vessels only,
 * therefore we write them, in the h5 file before we starting
 */
// static void PyPrepareForAdaptation(std::auto_ptr<VesselList3d> vl_,h5cpp::Group vesselgroup_, BloodFlowParameters *bfparams_, Adaption::Parameters *adap_params)
// {
//   vl_ = ReadVesselList3d(vesselgroup_, make_ptree("filter", true));
//   //to create the bc array, we need proper flow and pressure values!
//   //radii could then be overwriten
//   CalcFlow(vl_, *bfparams_);
//   h5cpp::Group grp_temp;
// #if 0
//   if( adap_params.write2File )
//   {
//     if(not out_.exists("recomputed"))
//     {
//       grp_temp = out_.create_group("recomputed");
//       ptree getEverytingPossible = make_ptree("w_adaption", true);
//       WriteVesselList3d(*vl_, grp_temp, getEverytingPossible);
//     }
//     else
//     {
//       // dont need this
//       //grp_temp = out_.open_group("recomputed");
//     }
//   }
// #endif
// }
//static py::object PyComputeAdaption(py::object py_vesselgroup, py::dict py_parameters, py::dict py_bfparams, py::object py_h5outputGroup)
static py::object PyComputeAdaption(const py::dict py_parameters, const py::dict py_bfparams, bool doOutput)
{
#ifndef TOTAL_SILENCE
  cout<<" PyComputeAdaption is called "<<endl;
#endif
#if 1
  //BloodFlowParameters bfparams = py::extract<BloodFlowParameters>(py_bfparams);
  BloodFlowParameters bfparams;
  InitBFParameters(&bfparams, &py_bfparams);
  Adaption::Parameters params;
  InitParameters(&params, &py_parameters);

#ifndef TOTAL_SILENCE
  cout<<" Parameters initialized "<<endl;
#endif
  

  std::tuple<uint,FlReal,FlReal, FlReal> return_state;
  
  return_state = runAdaption_Loop(params, bfparams, doOutput);
  
  
  return py::make_tuple(std::get<0>(return_state), std::get<1>(return_state), std::get<2>(return_state),std::get<3>(return_state));

#endif
} 


/* we decide to parse a string name here and NOT the h5py handle
 * this caused some major trouble in the past
 * and strings are easy to serialize!
 */
#ifdef USE_PAGMO
static py::object PydoAdaptionOptimization(py::str py_vesselgroup_str, py::dict py_parameters, py::object py_bfparams)
{
  #if 1
  cout<<" PydoAdaptionOptimization is called "<<endl;
#endif
  BloodFlowParameters bfparams_buffer = py::extract<BloodFlowParameters>(py_bfparams);
  BloodFlowParameters *bfparams = &bfparams_buffer;
  
  Adaption::Parameters *params = new Adaption::Parameters();
  InitParameters(params, py_parameters);
  

  // so kann man es besser machen:
//  std::string vesselListClass = "GRAPH";
//   try{
//     vesselListClass = vesselgroup.attrs().get<string>("CLASS");
//   }
//   catch(h5cpp::Exception &e)  // will programm abbruch fuer alle fehler ausser wenn CLASS attribute fehlt. Pruefe ob H5 exception. Andere exceptions machen programm abbruch.
//   {
//     cerr << "PyComputeAdaption: fall back to default vessel list reader (lattice based) due to error: ";
//     cerr << e.what();
//   }
  
//   if(vesselListClass == "REALWORLD")
//   {
//    // world = true;
//   }
//   else if (vesselListClass == "GRAPH")
//   {
//     // nothing to do
//   }
//   else
//   { // output the name of the class, too
//     cerr<<"PyComputeAdaption: Unknows CLASS: "<< vesselListClass << "; fall back to lattice based vessellist" << endl;
//   }
  
  //gain some topological info for adaptive parameters
//   int no_of_roots = 0;
//   for(int i=0;i<vl->GetNCount();++i)
//   {
//     VesselNode *nd = vl->GetNode(i);
//     if(nd->IsBoundary())no_of_roots++;
//   }
//   params->no_of_roots = no_of_roots;
  
  
  //starting with equal radii with parameters starting value is greate than 0.
  //otherwise we use the radii from the input file vary them a little bit
//   if(params->starting_radii>0.)
//   {
//     for(int i=0;i<vl->GetECount();++i)
//     {
//       vl->GetEdge(i)->r = params->starting_radii;
//     }
//   }
//   else
//   {
// #if 0
//     std::default_random_engine generator;
//     double initial_variation = 0.0;
//     std::uniform_real_distribution<double> distribution(-initial_variation,initial_variation);
//     for(int i=0;i<vl->GetECount();++i)
//     {
//       Vessel* v=vl->GetEdge(i);
// #if 1//leafes boundary radius as they are
//       auto it = vl->GetBCMap().find(v->NodeA());
//       auto it2 = vl->GetBCMap().find(v->NodeB());
//       if(it !=vl->GetBCMap().end() or it2 !=vl->GetBCMap().end())
//       {
// #ifdef DEBUG
// 	cout<<"skipping vessel #"<<v->Index()<< " since it has boundary nodes!!!"<<endl;
// #endif
// 	continue;
//       }
// #endif
//       v->r = v->r+v->r*distribution(generator);
//     }
// #endif
//   }
  
  std::cout<<"running real stuff.. starting pagmo"<<std::endl;
  // Initialise the MPI environment.
  // maybe it is a good advice to choose the same number of mpi island as in the mpirun -n command

  // Create a problem and an algorithm.
  std::string myVesselFileString = boost::python::extract<std::string>(py_vesselgroup_str);
  params->vesselFileName = myVesselFileString;
  pagmo::problem::adaption_problem prob(*params,*bfparams);
  
  //output what is read from the parameters
  std::printf("pop: %i, individuals: %i", params->pop, params->individuals);
  pagmo::population pop(prob,params->pop);
  pagmo::algorithm::pso algo(params->individuals);
//  pagmo::population pop(prob,3);
//  pagmo::algorithm::pso algo(5);
  algo.set_screen_output(true);
  algo.human_readable();
#if 0
  /* we do the evolution without islands
   * this does not crashe on the calcflow stuff!
   */
  for(int i=0;i<params->opt_iter;i++)
  {
    algo.evolve(pop);
    std::cout<<"Loop .. the champion is ..." <<std::endl;
    std::cout<< pop.champion().x <<std::endl;
  }
#else
/* in contrast this crashed when solving the linear system
 * ML_memory_free warning : header/tail mismatch
 * AZ_manage_memory
 * mostly at
 * solver_impl->Iterate(params.get<int>("max_iter", 50), params.get<double>("max_resid", 1.e-9));
 */
//   pagmo::island isl = pagmo::island(algo,pop);
//   for(int i=0;i<5;++i)
//   {
//     isl.evolve(1);
//     std::cout<<"Loop .. the champion is ..." <<std::endl;
//     std::cout<< pop.champion().x <<std::endl;
//   }
  
  // Create an archipelago of 10 MPI islands.
  std::printf("creating archipel\n");
#ifdef PAGMO_ENABLE_MPI
  std::printf("we assume this is mpi ready");
  pagmo::mpi_environment env;
  
  printf("you chose %i MPI processes.\n", env.get_size());
  printf("this is rank: %i\n", env.get_rank());
  printf("is this mt: %i\n", env.is_multithread());
  
  pagmo::archipelago a; // = pagmo::archipelago(algo,prob);
  for (int i = 0; i < env.get_size(); ++i) {
    a.push_back(pagmo::mpi_island(algo,prob,params->pop));
  }
#else
  pagmo::archipelago a = pagmo::archipelago(algo,prob,3,10);
#endif
  
  a.set_topology(pagmo::topology::ring());
  // Evolve the archipelago 10 times.
  a.human_readable();
  std::printf("start to evolve\n");
  a.evolve(3);
  std::printf("join islands\n");
  //a.join();
#endif
//   for(auto isl: a.get_islands())
//   {
//     std::cout<< isl->get_population().champion().x <<std::endl;
//   }
  std::cout<<"Done .. the champion is ..." <<std::endl;
  // OUTPUT OF SUCESSFULLY NETWORK WOULD BE FINE
  //std::cout<< pop.champion().x <<std::endl;
  //we return the results of the optimization to python and proceed there
  std::cout<<"return this back to python"<<std::endl;
  //return py::make_tuple(pop.champion().x[0],pop.champion().x[1],pop.champion().x[2]);
  return py::make_tuple(1.,2.,3.);
}
#endif


// Adaption::Parameters* AllocateParametersFromDict(const py::dict &d)
// {
//   std::auto_ptr<Adaption::Parameters> p(new Adaption::Parameters());
//   InitParameters(p.get(), d);
//   return p.release();
// }

void export_adaption_computation()
{ 
//   py::def("AllocateAdaptionParametersFromDict", &AllocateParametersFromDict, py::return_value_policy<py::manage_new_object>());
  py::def("computeAdaption", PyComputeAdaption);
  py::def("testAdaption", TestAdaption);
#ifdef USE_PAGMO
  py::def("doAdaptionOptimization",PydoAdaptionOptimization);
#endif
}
  
}//namespace

#ifdef DEBUG
BOOST_PYTHON_MODULE(libadaption_d)
#else
BOOST_PYTHON_MODULE(libadaption_)
#endif
{
  PyEval_InitThreads(); // need for release of the GIL (http://stackoverflow.com/questions/8009613/boost-python-not-supporting-parallelism)
  my::checkAbort = PyCheckAbort; // since this is the python module, this is set to use the python signal check function
  Adaption::export_adaption_computation();
}
