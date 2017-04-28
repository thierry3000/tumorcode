#include "adaption_as_pagmo_problem.h"
#ifdef USE_PAGMO
int doAdaptionOptimization_without_py()
{
  #if 1
  cout<<" PydoAdaptionOptimization is called "<<endl;
#endif
  //BloodFlowParameters bfparams_buffer = py::extract<BloodFlowParameters>(py_bfparams);
  //BloodFlowParameters *bfparams = &bfparams_buffer;
  BloodFlowParameters *bfparams = new BloodFlowParameters();
  
  Adaption::Parameters *params = new Adaption::Parameters();
  //InitParameters(params, py_parameters);
  

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
  int mc_steps=0;
  int dim=0;
#ifdef PAGMO_ENABLE_MPI
  pagmo::mpi_environment env;
  mc_steps = 10000000;
  dim = 400;
#else
  mc_steps = 10;
  dim = 4;
#endif
  // Create a problem and an algorithm.
  //std::string myVesselFileString = boost::python::extract<std::string>(py_vesselgroup_str);
  std::string myVesselFileString = "/localdisk/thierry/vessel_trees_better/my_chosen/PSO_data_vessels-large_2d-typeE-17x1L600-sample05_adption_p_human_guess.h5";
  pagmo::problem::adaption_problem prob(myVesselFileString,*params,*bfparams);
  std::printf("pop: %i, individuals: %i", params->pop, params->individuals);
//  pagmo::population pop(prob,params->pop);
//  pagmo::algorithm::pso algo(params->individuals);
  pagmo::population pop(prob,3);
  pagmo::algorithm::pso algo(5);
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
  pagmo::archipelago a; // = pagmo::archipelago(algo,prob);
  for (int i = 0; i < 3; ++i) {
    a.push_back(pagmo::mpi_island(algo,prob,1));
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
  //std::cout<< pop.champion().x <<std::endl;
  //we return the results of the optimization to python and proceed there
  std::cout<<"return this back to python"<<std::endl;
  //return py::make_tuple(pop.champion().x[0],pop.champion().x[1],pop.champion().x[2]);
  //return py::make_tuple(1.,2.,3.);
  return 42;
}
#endif

int main()
{
  doAdaptionOptimization_without_py();
}
