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
#include <iostream>
//#include "../python_krebsutils/python-helpers.h"
#include "numpy.hpp"
//#include "../adaption/adaption_model2.h"
//#include "../common/calcflow.h"
#include <algorithm>

#include "hdf5.h"
#define FILE "/localdisk/thierry/vessel_trees_better/my_chosen/PSO_data_vessels-large_2d-typeE-17x1L600-sample05_adption_p_human_guess.h5"
//#define PAGMO_ENABLE_MPI --> in CMakeLists
#include <pagmo/pagmo.h>

using namespace pagmo;

void test_2()
{
  std::cout<<"running pagmo test 2"<<std::endl;
//   hid_t       file_id, group_id;  /* identifiers */
//   herr_t      status;
//   
//   /* Create a new file using default properties. */
//   //file_id = H5Fcreate(FILE, H5F_ACC_RDONLY, H5P_DEFAULT, H5P_DEFAULT);
//   file_id = H5Fopen(FILE, H5F_ACC_RDONLY, H5P_DEFAULT);
//   group_id = H5Gopen1(file_id,"parameters");
//   //H5Gget_info(group_id);
//   status = H5Gclose(group_id);
//   status = H5Fclose(file_id);
}

void test_1()
{
  std::cout<<"running pagmo test 1"<<std::endl;
  // Initialise the MPI environment.
  int mc_steps=0;
  int dim=0;
#ifdef PAGMO_ENABLE_MPI
  
  mpi_environment env;
  
  printf("you chose %i MPI processes.\n", env.get_size());
  printf("this is rank: %i\n", env.get_rank());
  printf("is this mt: %i\n", env.is_multithread());
  
  mc_steps = 100000;
  dim = 400;
#else
  mc_steps = 100000000;
  dim = 40;
#endif
  // Create a problem and an algorithm.
  problem::dejong prob(dim);
  algorithm::monte_carlo algo(mc_steps);
  algo.set_screen_output(true);
  // Create an archipelago of 48 MPI islands.
  archipelago a;
  int population = 3;
  a.set_topology(topology::ring());
  for (int i = 0; i < env.get_size(); ++i) {
#ifdef PAGMO_ENABLE_MPI
    std::cout<<"I will solve this problem!!!"<<std::endl;
	  a.push_back(mpi_island(algo,prob,population));
#else
	  a.push_back(island(algo,prob,1));
#endif
  }
  // Evolve the archipelago 10 times.
  a.evolve(5);
  a.join();
  
}



int main(int argc, char **argv)
{
  //throw std::runtime_error("implement parameter handling with boost program options");
  std::cout<<"begin test_pagmo"<<std::endl;
  
  test_1();
  test_2();
#if 0
  my::MultiprocessingInitializer mpinit(argc, argv);
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
  Int3 size(100, 100, 10);
  int repetitions = 1;
  int num_threads = 1;
  my::ProgramOptions opts;
  opts.add(size, "--size");
  opts.add(repetitions, "--repeat");
  opts.add(num_threads, "--num-threads");
  opts.parse(argv, argc);
  my::SetNumThreads(num_threads);
  cout << "size = " << size << ", " << num_threads << " threads" << endl;
  for (int i=0; i<repetitions; ++i)
  {
//     cout << "WITH SPARSE MATRIX" << endl;
//     doit(size);
    cout << "MATRIX FREE" << endl;
    doit(size);
    cout << "*** memory cons: " << (double(GetMemoryUsage().rss_peak)/(1024*1024)) << " mb" << endl;
  }
#endif
  return 0;
}
