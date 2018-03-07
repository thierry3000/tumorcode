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
#include "shared-objects.h"
#include "continuum-utils.h"
#include "continuum-grid.h"
#include "common/trilinos_linsys_construction.h"
#include "common/vessels3d.h"
#include "mwlib/math_ext.h"

/*
 * solve: laplace U + lambda (U0 - U) = 0
 * where U = local volume average of pysource
 * and lambda = local volume fraction clipped at 1 times source_factor
 */
#if BOOST_VERSION>106300
np::ndarray CalcIntervascularInterpolationField(
                                      np::ndarray pypos, 
                                      np::ndarray pyedges, 
                                      np::ndarray pyradius,
                                      np::ndarray pysource,
                                      const py::object &py_ldfield,
                                      double source_factor,
                                      int samples_per_cell)
{
  LatticeDataQuad3d ld = py::extract<LatticeDataQuad3d>(py_ldfield);
  ContinuumGrid grid(ld);
  DomainDecomposition mtboxes;
  mtboxes.init(MakeMtBoxGrid(grid.Box(), Int3(32, 32, 32)));
  
//   np::arrayt<float> pos(pypos);
//   np::arrayt<int> edges(pyedges);
//   np::arrayt<float> radius(pyradius);
//   np::arrayt<float> source(pysource);
  int cnt = pyedges.get_shape()[0];
  
  Array3d<float> source_gridded(ld.Box(), 0.);
  Array3d<float> weight_gridded(ld.Box(), 0.);
  
  int num_total_samples = 0;
  CylinderNetworkSampler sampler;
  sampler.Init(ld.Scale(), make_ptree("samples_per_cell", samples_per_cell));
  
  for(int i=0; i<cnt; ++i)
  {
    Float3 p0, p1;
    float local_source0, local_source1;
    int a = py::extract<int>(pyedges[i][0]);
    int b = py::extract<int>(pyedges[i][1]);
    for (int j=0; j<3; ++j)
    {
//       p0[j] = pos(edges(i,0), j);
//       p1[j] = pos(edges(i,1), j);
      p0[j] = py::extract<float>(pypos[a][j]);
      p1[j] = py::extract<float>(pypos[b][j]);
    }
//     local_source0 = source(i, 0);
//     local_source1 = source(i, 1);
    local_source0 = py::extract<float>(pysource[i][0]);
    local_source1 = py::extract<float>(pysource[i][1]);
    
    sampler.Set(p0, p1, py::extract<float>(pyradius[i])); // 3rd arg is the radius
    int num_samples = sampler.GenerateVolumeSamples();
    for (int k=0; k<num_samples; ++k)
    {
      float local_source = local_source0 + sampler.GetSample(k).fpos[2] * (local_source1 - local_source0); // linear interpolation along axis
      AddSmoothDelta(source_gridded, ld.Box(), ld, 3, sampler.GetSample(k).wpos, local_source*sampler.weight_per_volume);
      AddSmoothDelta(weight_gridded, ld.Box(), ld, 3, sampler.GetSample(k).wpos, sampler.weight_per_volume);
    }
  }
  
#if 1
  FiniteVolumeMatrixBuilder mb;
  mb.Init7Point(grid.ld, grid.dim);
  #pragma omp parallel
  {
    BOOST_FOREACH(const DomainDecomposition::ThreadBox bbox, mtboxes.getCurrentThreadRange())
    {
      mb.AddDiffusion<> (bbox, ConstValueFunctor<float>(1.), -1.0);
      FOR_BBOX3(p, bbox)
      {
        float w = weight_gridded(p);
        float s = source_gridded(p);
        if (w > 0)
        {
          s /= w;
          w = source_factor * my::cut(w, 0.f, 1.f);
        }
        mb.AddLocally(p, w, w*s);
      }
    }
  }
  
  weight_gridded.clear();
  source_gridded.clear();
  
  //Epetra_Vector lhs(mb.rhs->Map());
  Teuchos::RCP<Epetra_Vector> lhs = Teuchos::rcp(new Epetra_Vector(mb.rhs->Map()));
  ptree solver_params = make_ptree("preconditioner","multigrid")("verbosity", "normal")("use_smoothed_aggregation", false)("max_iter", 500)("conv","rhs")("max_resid",1.e-8);
  
  {
     EllipticEquationSolver solver(mb.m, mb.rhs, solver_params);
     //solver.init(mb.m, mb.rhs, solver_params);
//     solver.solve(lhs);
//    EllipticEquationSolver solver(mb.m, mb.rhs, solver_params);
    //solver.init(mb.m, mb.rhs, solver_params);
    solver.solve(lhs);
  }

//   np::arrayt<float> res(np::zeros(3, Cast<Py_ssize_t>(::Size(ld.Box())).data(), np::getItemtype<float>()));
  np::ndarray res = np::zeros(py::tuple(Cast<Py_ssize_t>(::Size(ld.Box())).data()), np::dtype::get_builtin<float>());

  #pragma omp parallel
  {
    BOOST_FOREACH(const DomainDecomposition::ThreadBox bbox, mtboxes.getCurrentThreadRange())
    {
      FOR_BBOX3(p, bbox)
      {
        res[p[0]][p[1]][p[2]] = (*lhs)[grid.ld.LatticeToSite(p)];
        //res[p[0]][p[1]][p[2]] = (*lhs)[0];
      }
    }
  }
#else 
  /* this is some debugging code to just write down the source contribution. You should get a picture
   * which looks like a povray rendering of the blood pressure for instance if you used the blood pressure.
   */
  np::arrayt<float> res(np::zeros(3, ::Size(ld.Box()).data(), np::getItemtype<float>()));
  #pragma omp parallel
  {
    BOOST_FOREACH(const DomainDecomposition::ThreadBox bbox, mtboxes.getCurrentThreadRange())
    {
      FOR_BBOX3(p, bbox)
      {
        if (weight_gridded(p) > 0)
          res(p[0],p[1],p[2]) = source_gridded(p) / weight_gridded(p);
        else
          res(p[0],p[1],p[2]) = 0;
      }
    }
  }
#endif
  return res;  
}
#else
np::arraytbase CalcIntervascularInterpolationField(
                                      nm::array pypos, 
                                      nm::array pyedges, 
                                      nm::array pyradius,
                                      nm::array pysource,
                                      const py::object &py_ldfield,
                                      double source_factor,
                                      int samples_per_cell)
{
  LatticeDataQuad3d ld = py::extract<LatticeDataQuad3d>(py_ldfield);
  ContinuumGrid grid(ld);
  DomainDecomposition mtboxes;
  mtboxes.init(MakeMtBoxGrid(grid.Box(), Int3(32, 32, 32)));
  
  np::arrayt<float> pos(pypos);
  np::arrayt<int> edges(pyedges);
  np::arrayt<float> radius(pyradius);
  np::arrayt<float> source(pysource);
  int cnt = edges.shape()[0];
  
  Array3d<float> source_gridded(ld.Box(), 0.);
  Array3d<float> weight_gridded(ld.Box(), 0.);
  
  int num_total_samples = 0;
  CylinderNetworkSampler sampler;
  sampler.Init(ld.Scale(), make_ptree("samples_per_cell", samples_per_cell));
  
  for(int i=0; i<cnt; ++i)
  {
    Float3 p0, p1;
    float local_source0, local_source1;
    for (int j=0; j<3; ++j)
    {
      p0[j] = pos(edges(i,0), j);
      p1[j] = pos(edges(i,1), j);
    }
    local_source0 = source(i, 0);
    local_source1 = source(i, 1);
    
    sampler.Set(p0, p1, radius(i)); // 3rd arg is the radius
    int num_samples = sampler.GenerateVolumeSamples();
    for (int k=0; k<num_samples; ++k)
    {
      float local_source = local_source0 + sampler.GetSample(k).fpos[2] * (local_source1 - local_source0); // linear interpolation along axis
      AddSmoothDelta(source_gridded, ld.Box(), ld, 3, sampler.GetSample(k).wpos, local_source*sampler.weight_per_volume);
      AddSmoothDelta(weight_gridded, ld.Box(), ld, 3, sampler.GetSample(k).wpos, sampler.weight_per_volume);
    }
  }
  
#if 1
  FiniteVolumeMatrixBuilder mb;
  mb.Init7Point(grid.ld, grid.dim);
  #pragma omp parallel
  {
    BOOST_FOREACH(const DomainDecomposition::ThreadBox bbox, mtboxes.getCurrentThreadRange())
    {
      mb.AddDiffusion<> (bbox, ConstValueFunctor<float>(1.), -1.0);
      FOR_BBOX3(p, bbox)
      {
        float w = weight_gridded(p);
        float s = source_gridded(p);
        if (w > 0)
        {
          s /= w;
          w = source_factor * my::cut(w, 0.f, 1.f);
        }
        mb.AddLocally(p, w, w*s);
      }
    }
  }
  
  weight_gridded.clear();
  source_gridded.clear();
  
  //Epetra_Vector lhs(mb.rhs->Map());
  Teuchos::RCP<Epetra_Vector> lhs = Teuchos::rcp(new Epetra_Vector(mb.rhs->Map()));
  ptree solver_params = make_ptree("preconditioner","multigrid")("verbosity", "normal")("use_smoothed_aggregation", false)("max_iter", 500)("conv","rhs")("max_resid",1.e-8);
  
  {
    EllipticEquationSolver solver(mb.m, mb.rhs, solver_params);
    //solver.init(mb.m, mb.rhs, solver_params);
//     solver.solve(lhs);
//    EllipticEquationSolver solver(mb.m, mb.rhs, solver_params);
    //solver.init(mb.m, mb.rhs, solver_params);
    solver.solve(lhs);
  }

  np::arrayt<float> res(np::zeros(3, Cast<Py_ssize_t>(::Size(ld.Box())).data(), np::getItemtype<float>()));

  #pragma omp parallel
  {
    BOOST_FOREACH(const DomainDecomposition::ThreadBox bbox, mtboxes.getCurrentThreadRange())
    {
      FOR_BBOX3(p, bbox)
      {
        res(p[0],p[1],p[2]) = (*lhs)[grid.ld.LatticeToSite(p)];
      }
    }
  }
#else 
  /* this is some debugging code to just write down the source contribution. You should get a picture
   * which looks like a povray rendering of the blood pressure for instance if you used the blood pressure.
   */
  np::arrayt<float> res(np::zeros(3, ::Size(ld.Box()).data(), np::getItemtype<float>()));
  #pragma omp parallel
  {
    BOOST_FOREACH(const DomainDecomposition::ThreadBox bbox, mtboxes.getCurrentThreadRange())
    {
      FOR_BBOX3(p, bbox)
      {
        if (weight_gridded(p) > 0)
          res(p[0],p[1],p[2]) = source_gridded(p) / weight_gridded(p);
        else
          res(p[0],p[1],p[2]) = 0;
      }
    }
  }
#endif
  return res;  
}
#endif

void export_compute_interpolation_field()
{
  py::def("CalcIntervascularInterpolationField_", CalcIntervascularInterpolationField);
}
