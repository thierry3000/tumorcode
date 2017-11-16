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
#include "trilinos_linsys_construction.h"
#include "trilinos_matrixfree.h"
#include <fenv.h>

#include "mwlib/random.h"

#include "shared-objects.h"
#include <boost/program_options.hpp>

namespace po=boost::program_options;

struct ConstFaceVarFunctor
{
  double val;
  ConstFaceVarFunctor(double val) : val(val) {}
  double operator()(const Int3& p, const Int3& q, int axis, int side) const { return val; }
  double operator()(int axis, const Int3& p) const { return val; }
};


// compute \laplace c - lambda c = -r
void doit(const Int3 size)
{
  LatticeDataQuad3d ld;
  ld.Init(Int3(size), 1.);
  ld.SetOriginPosition(-ld.GetWorldBox().max*0.5);

  //const float dt = 100;
  const float rmax = ld.GetWorldBox().max.norm();

  Random rnd;
  Array3d<float> randomness(ld.Box());
  FOR_BBOX3(p, ld.Box())
  {
    randomness(p) = rnd.Get01()*0.1 + 1.;
  }
  int num_dof = ld.NumSites();
  
  uint64 mem_usage1 = GetMemoryUsage().rss;
  my::Time t_;
  FiniteVolumeMatrixBuilder mb;
  mb.Init7Point(ld, 3);

  DynArray<BBox3> mtboxes = MakeMtBoxGrid(ld.Box());
  #pragma omp parallel
  {
    //Random rnd(1234 + 23467 * omp_get_thread_num());
    #pragma omp for nowait
    for (int i=0; i<mtboxes.size(); ++i)
    {
      mb.AddDiffusion2(mtboxes[i], 3, ConstFaceVarFunctor(1), 1., false, randomness);
      FOR_BBOX3(p, mtboxes[i])
      {
        float r = ld.LatticeToWorld(p).norm();
        float rhs = my::smooth_heaviside<float>(rmax*0.2 - r, ld.Scale()*2);
        //mb.AddLocally(p, prefactor*(-1.)/my::sqr(200.), prefactor*rhs);
        //float q = 0.05 * rnd.Get11();
        //float lc = 0., cc = 0.;
        //if (q < 0.)
//           lc = q;
//         else
//           cc = q;
//        mb.AddLocally(p, 1 - dt * lc, rhs + dt * cc);
        mb.AddLocally(p, -1./my::sqr(10*ld.Scale()), -rhs);
      }
    }
  }
  cout << "cons time: " << (my::Time() - t_).to_ms() << endl;
  cout << "matrix memory cons: " << double(GetMemoryUsage().rss - mem_usage1)/num_dof << " b/dof" << endl;
  mem_usage1 = GetMemoryUsage().rss;
  
  //double scaling = 1./mb.m->NormInf();
  double scaling = -1.;
  
  mb.m->Scale(scaling);
  mb.rhs->Scale(scaling);
  Teuchos::RCP<Epetra_Vector> lhs = Teuchos::rcp(new Epetra_Vector(mb.rhs->Map()));
  
  boost::property_tree::ptree pt;
  pt.put("output", 1);
  pt.put("use_multigrid", false);
  pt.put("output_matrix", false);
  pt.put("use_smoothed_aggregation", true);
  pt.put("solver", "cg");
  pt.put("max_iter", 100);
  pt.put("throw", false);

  SolveEllipticEquation(mb.m, mb.rhs, lhs, pt);

  //EllipticEquationSolver solver;
  //solver.init(*mb.m, *mb.rhs, pt);
  //solver.solve(lhs);

  cout << "solve memory cons: " << double(GetMemoryUsage().rss_peak - mem_usage1)/num_dof << " b/dof" << endl;
  
  Array3d<float> res(ld.Box());
  FOR_BBOX3(p, ld.Box())
  {
    res(p) = (*lhs)[ld.LatticeToSite(p)];
  }

  Image img;
  DrawArray(img, res, DrawArrayOpts().outputRange());
  img.Write("elliptic_solver_test_1.png");

#if 0
  mb.m->PutScalar(0);
  #pragma omp parallel
  {
    //Random rnd(1234 + 23467 * omp_get_thread_num());
    #pragma omp for nowait
    for (int i=0; i<mtboxes.size(); ++i)
    {
      mb.AddDiffusion2(mtboxes[i], 3, ConstFaceVarFunctor(1), 1., false, randomness);
      FOR_BBOX3(p, mtboxes[i])
      {
        float r = ld.LatticeToWorld(p).norm();
        float rhs = my::smooth_heaviside<float>(rmax*0.2 - r, ld.Scale()*2);
        mb.AddLocally(p, -1./my::sqr(30*ld.Scale()),0.);
      }
    }
  }
  mb.m->Scale(scaling);
  //solver.init(*mb.m, *mb.rhs, pt);
  //solver.solve(lhs);
  SolveEllipticEquation(*mb.m, *mb.rhs, lhs, pt);
  
  FOR_BBOX3(p, ld.Box())
  {
    res(p) = lhs[ld.LatticeToSite(p)];
  }
  DrawArray(img, res, DrawArrayOpts().outputRange());
  img.Write("elliptic_solver_test_2.png");
#endif
}



/*

// compute \laplace c - lambda c = -r
void doit2(const Int3 size)
{
  LatticeDataQuad3d ld;
  ld.Init(Int3(size), 1.);
  ld.SetWorldPosition(-ld.GetWorldBox().max*0.5);

  //const float dt = 100;
  const float rmax = ld.GetWorldBox().max.norm();

  Random rnd;
  Array3d<float> randomness(ld.Box());
  FOR_BBOX3(p, ld.Box())
  {
    randomness(p) = rnd.Get01()*0.1 + 1.;
  }
  int num_dof = ld.NumSites();

  uint64 mem_usage1 = GetMemoryUsage().rss;
  my::Time t_;


  Epetra_SerialComm epetra_comm;
  Epetra_Map epetra_map(ld.NumSites(), 0, epetra_comm);

  typedef TrilinosMatrixFreeEllipticOperator<ConstFaceVarFunctor, Array3d<float> > OperatorType;

  Array3d<float> diag_offset(ld.Box());
  Epetra_Vector rhs(epetra_map);
  {
    FOR_BBOX3(p, ld.Box())
    {
      float r = ld.LatticeToWorld(p).norm();
      //float loc_rhs = my::smooth_heaviside<float>(rmax*0.2 - r, ld.Scale()*2);
      //rhs[ld.LatticeToSite(p)] = my::sqr((1./sqrt(2.))*ld.LatticeToWorld(p).dot(Float3(1.,-1, 0)));
//      rhs[ld.LatticeToSite(p)] = std::max(std::abs(ld.LatticeToWorld(p)[0]), std::abs(ld.LatticeToWorld(p)[1]));
//       diag_offset(p) = -10.;

      float loc_rhs = (p == Int3(20, 20, 0)) ? 1 : 0;
      rhs[ld.LatticeToSite(p)] = -loc_rhs;
      diag_offset(p) = -1./my::sqr(100*ld.Scale());
    }
  }

  rhs.Scale(-1.);
  boost::scoped_ptr<OperatorType> op(
    new OperatorType(
      ld, size[2]==1 ? 2 : 3, -1., ConstFaceVarFunctor(1.), diag_offset, epetra_map));

  mem_usage1 = GetMemoryUsage().rss;
  
  Epetra_Vector lhs(rhs.Map(), true);

  boost::property_tree::ptree pt;
  pt.put("output", 1);
  pt.put("use_multigrid", false);
  pt.put("max_iter", 100);
  pt.put("throw", false);

  { my::Time t_;
  //SolveEllipticEquation(*op, rhs, lhs, pt);
  op->Apply(rhs, lhs);
  cout << "time: " << (my::Time() - t_) << endl;}

  cout << "solve memory cons: " << double(GetMemoryUsage().rss_peak - mem_usage1)/num_dof << " b/dof" << endl;

  Epetra_Vector residual(lhs.Map(), true);
  op->Apply(lhs, residual);
  
  Array3d<float> lhs_arr(ld.Box());
  Array3d<float> rhs_arr(ld.Box());
  Array3d<float> resid_arr(ld.Box());
  FOR_BBOX3(p, ld.Box())
  {
    lhs_arr(p) = lhs[ld.LatticeToSite(p)];
    rhs_arr(p) = -rhs[ld.LatticeToSite(p)];
    resid_arr(p) = residual[ld.LatticeToSite(p)];
  }

  DynArray<Image> images;
  images.push_back(DrawArray(lhs_arr, DrawArrayOpts().outputRange()));
  images.push_back(DrawArray(rhs_arr, DrawArrayOpts().outputRange()));
  images.push_back(DrawArray(resid_arr, DrawArrayOpts().outputRange()));
  Image bigimg;
  DrawImageGrid(bigimg, images);
  bigimg.Write("elliptic_solver_test_1.png");
}*/




int main(int argc, char **argv)
{
  // Declare the supported options.
  po::options_description desc("Allowed options");
  desc.add_options()
      ("help,h", "produce help message")
      ("size,s", po::value<Int3>(), "set size")
      ("repetitions,r" ,po::value<int>(), "set repeat")
      ("num_threads,n" ,po::value<int>(), "set threads")
  ;
  //throw std::runtime_error("implement parameter handling with boost program options");
#if 1
  my::MultiprocessingInitializer mpinit(argc, argv);
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
  Int3 size(100, 100, 10);
  int repetitions = 1;
  int num_threads = 1;
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);
  
  if (vm.count("help")) 
  {
    cout << desc << "\n";
    return 1;
  }
  if (vm.count("size")) 
  {
    cout << "size was set to " << vm["size"].as<Int3>() << ".\n";
  } else 
  {
    cout << "size was not set using 100,100,10.\n";
  }
  if (vm.count("repetitions")) 
  {
    cout << "repetitions was set to " << vm["repetitions"].as<int>() << ".\n";
  } else 
  {
    cout << "repetitions was not set using 1.\n";
  }
  if (vm.count("num_threads")) 
  {
    cout << "num_threads was set to " << vm["num_threads"].as<int>() << ".\n";
  } else 
  {
    cout << "num_threads was not set using 1.\n";
  }
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
}
