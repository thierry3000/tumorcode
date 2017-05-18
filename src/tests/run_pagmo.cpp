#include "../adaption/adaption_as_pagmo_problem.h"
#ifdef USE_PAGMO
void doAdaptionOptimization_without_py()
{
  cout<<" PydoAdaptionOptimization is called "<<endl;

  //BloodFlowParameters bfparams_buffer = py::extract<BloodFlowParameters>(py_bfparams);
  //BloodFlowParameters *bfparams = &bfparams_buffer;
  BloodFlowParameters *bfparams = new BloodFlowParameters();
  
  Adaption::Parameters *params = new Adaption::Parameters();
  //InitParameters(params, py_parameters);
  
  std::cout<<"running real stuff.. starting pagmo"<<std::endl;
  // Initialise the MPI environment.
  
  // Create a problem and an algorithm.
  //std::string myVesselFileString = boost::python::extract<std::string>(py_vesselgroup_str);
  std::string myVesselFileString = "/localdisk/thierry/vessel_trees_better/my_chosen/PSO_data_vessels-large_2d-typeE-17x1L600-sample05_adption_p_human_guess.h5";
  pagmo::problem::adaption_problem prob(myVesselFileString,*params,*bfparams);
  std::printf("pop in params: %i, individuals in params: %i", params->pop, params->individuals);
  std::printf("not used in run_pagmo");
//  pagmo::population pop(prob,params->pop);
//  pagmo::algorithm::pso algo(params->individuals);
  //pagmo::population pop(prob,1);
  pagmo::algorithm::pso algo(1);
  algo.set_screen_output(true);
  algo.human_readable();
  
  // Create an archipelago of 10 MPI islands.
  std::printf("creating archipel\n");
#ifdef PAGMO_ENABLE_MPI
  std::printf("PAGMO_ENABLE_MPI enabled");
  pagmo::mpi_environment env;
  pagmo::archipelago a; // = pagmo::archipelago(algo,prob);
  for (int i = 0; i < env.get_size(); ++i) {
    a.push_back(pagmo::mpi_island(algo,prob,12));
  }
#else
  pagmo::archipelago a = pagmo::archipelago(algo,prob,1,1);
#endif
  
  a.set_topology(pagmo::topology::ring());
  // Evolve the archipelago 10 times.
  a.human_readable();
  std::printf("start to evolve\n");
  a.evolve(10);
  std::printf("join islands\n");
  a.join();

  for(auto isl: a.get_islands())
  {
    std::cout<< isl->get_population().champion().x <<std::endl;
  }
  
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
