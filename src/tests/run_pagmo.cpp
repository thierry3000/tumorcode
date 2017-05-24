#include "../adaption/adaption_as_pagmo_problem.h"
#ifdef PAGMO_ENABLE_MPI
  //std::printf("PAGMO_ENABLE_MPI enabled");
  pagmo::mpi_environment env;
#endif
#ifdef USE_PAGMO
void doAdaptionOptimization_without_py()
{

  //BloodFlowParameters bfparams_buffer = py::extract<BloodFlowParameters>(py_bfparams);
  //BloodFlowParameters *bfparams = &bfparams_buffer;
  BloodFlowParameters *bfparams = new BloodFlowParameters();
  
  Adaption::Parameters *params = new Adaption::Parameters();
  params->k_c = 1.42;
  params->k_m = 1.42;
  params->k_s = 1.42;
  params->Q_refdot = 42;
  params->max_nun_iterations=150;
  params->qdev=0.1;
  params->starting_radii=5;
  params->delta_t=0.1;
  params->cond_length=1500;
  params->vesselGroupName = "adaption/vessels_after_adaption";
  params->vesselFileName = "/localdisk/thierry/vessel_trees_better/my_chosen/PSO_data_vessels-large_2d-typeE-17x1L600-sample05_adption_p_human_guess.h5";

//   no_of_roots = 42;
//   max_flow = 42;
//   min_flow = 42;
//   avgRootNodeConductivity = 0;
//   S_0=42;
//   cond_length=4242;
//   tum_manitulate_s1=tum_manitulate_s2=tum_manitulate_s3=tum_manitulate_s4=tum_manitulate_s5=false;
//   /*
//    * important value!
//    * if adaption without tumor this goes to this default value
//    */
//   radMin_for_kill = 2.5;
//   write2File = true;
//   boundary_Condition_handling = KEEP;
//   a_pressure = 1.8;
//   a_flow = 200000.;
//   pop = 5;
//   individuals = 42;
//   opt_iter = 1;
  
  //InitParameters(params, py_parameters);
  
  std::cout<<"running real stuff.. starting pagmo"<<std::endl;
  // Initialise the MPI environment.
  
  // Create a problem and an algorithm.
  //std::string myVesselFileString = boost::python::extract<std::string>(py_vesselgroup_str);
  //std::string myVesselFileString = "/localdisk/thierry/vessel_trees_better/my_chosen/PSO_data_vessels-large_2d-typeE-17x1L600-sample05_adption_p_human_guess.h5";
  pagmo::problem::adaption_problem prob(*params,*bfparams);
  //std::printf("pop in params: %i, individuals in params: %i", params->pop, params->individuals);
  //std::printf("not used in run_pagmo");
  pagmo::population pop(prob,2);
//  pagmo::algorithm::pso algo(params->individuals);
  //pagmo::population pop(prob,1);
  pagmo::algorithm::pso algo(2);
  algo.set_screen_output(true);
  algo.human_readable();
  // this follows
  //https://github.com/esa/pagmo/blob/master/examples/evolve_spheres.cpp
  // Create an archipelago of 10 MPI islands.
  std::printf("creating archipel\n");
#ifdef PAGMO_ENABLE_MPI
//   std::printf("PAGMO_ENABLE_MPI enabled");
//   pagmo::mpi_environment env;
  pagmo::archipelago a; // = pagmo::archipelago(algo,prob);
  for (int i = 0; i < env.get_size(); ++i) {
    //a.push_back(pagmo::mpi_island(algo,prob,12));
    a.push_back(pagmo::mpi_island(algo,pop));
  }
#else
  pagmo::archipelago a = pagmo::archipelago(algo,prob,1,1);
#endif
  
  a.set_topology(pagmo::topology::fully_connected());
  // Evolve the archipelago 10 times.
  a.human_readable();
  std::printf("start to evolve\n");
  //a.evolve(30);
  
  for(int i=0;i<2;i++)
  {
    a.evolve(1);
    for(auto isl: a.get_islands())
    {
      std::cout<<"current pop: "<<std::endl;
      std::cout<< isl->get_population().champion().x <<std::endl;
    }
  }
  
  std::printf("join islands\n");
  a.join();
  
  int idx = archi_best_idx(a);
  std::cout << "and the winner is ......" << "\n" << a.get_island(idx)->get_population().champion().x << std::endl;
  //we return the results of the optimization to python and proceed there
}
#endif

int main()
{
#ifdef USE_PAGMO
  doAdaptionOptimization_without_py();
#endif
}
