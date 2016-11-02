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
#include "python-helpers.h"
#include "trilinos_matrixfree.h"
#include "trilinos_linsys_construction.h"
#include "common.h"


namespace EllipticSolverTest
{

struct ConstFaceVarFunctor
{
  double val;
  ConstFaceVarFunctor(double val) : val(val) {}
  double operator()(const Int3& p, const Int3& q, int axis, int side) const { return val; }
  double operator()(int axis, const Int3& p) const { return val; }
};


static void InitGrid(const Int3 &size, ContinuumGrid &grid, DomainDecomposition &mtboxes)
{
  LatticeDataQuad3d ld;
  ld.Init(Int3(size), 1.);
  ld.SetOriginPosition(-ld.GetWorldBox().max*0.5);
  grid.init(ld, size[2]>1 ? 3 : size[1]>1 ? 2 : 1);
  mtboxes.init(MakeMtBoxGrid(ld.Box(), Int3(64, 32, 32)));
}




py::object testOperatorApplication(const Int3 &size, bool write_image)
{
  ContinuumGrid grid;
  DomainDecomposition mtboxes;
  InitGrid(size, grid, mtboxes);
  mtboxes.print(cout);
  Array3df diag_offset(grid.Box(), Cons::DONT);
  Array3dd  lhs(grid.Box(), Cons::DONT),
            result(grid.Box(), Cons::DONT);

  
  const float rmax = grid.ld.GetWorldBox().max.norm();
  #pragma omp parallel
  {
    BOOST_FOREACH(const BBox3 &bbox, mtboxes.getCurrentThreadRange())
    {
      FOR_BBOX3(p, bbox)
      {
        float r = grid.ld.LatticeToWorld(p).norm();
        //float loc_rhs = my::smooth_heaviside<float>(rmax*0.2 - r, ld.Scale()*2);
        //rhs[ld.LatticeToSite(p)] = my::sqr((1./sqrt(2.))*ld.LatticeToWorld(p).dot(Float3(1.,-1, 0)));
  //      rhs[ld.LatticeToSite(p)] = std::max(std::abs(ld.LatticeToWorld(p)[0]), std::abs(ld.LatticeToWorld(p)[1]));
  //       diag_offset(p) = -10.;

        float loc_lhs = (p == Int3(20, 20, 0)) ? 1 : 0;
        lhs(p) = loc_lhs;
        diag_offset(p) = -1./my::sqr(100*grid.ld.Scale());
        result(p) = 0.;
      }
    }
  }
  
  typedef TrilinosMatrixFreeEllipticOperator<ConstFaceVarFunctor, Array3df> EpetraOperatorType;

#ifdef EPETRA_MPI
  Epetra_MpiComm epetra_comm(MPI_COMM_SELF);
#endif
#ifndef EPETRA_MPI
  Epetra_SerialComm epetra_comm;
#endif
  Epetra_Map epetra_map((int)grid.ld.NumSites(), 0, epetra_comm);
  
  boost::scoped_ptr<EpetraOperatorType> op(new EpetraOperatorType(
      grid.ld, grid.dim, mtboxes, -1., ConstFaceVarFunctor(1.), diag_offset, epetra_map));

  Epetra_Vector epetra_result(Epetra_DataAccess::View, epetra_map, result.getPtr());
  Epetra_Vector epetra_lhs(Epetra_DataAccess::View, epetra_map, lhs.getPtr());
  
  uint64 mem_usage = GetMemoryUsage().rss;
  my::Time t_;
  
  op->Apply(epetra_lhs, epetra_result);

  double result_time_ms = (my::Time() - t_).to_ms();
  double result_mem_bytes_per_dof = double(GetMemoryUsage().rss_peak - mem_usage)/grid.ld.NumSites();

  if (write_image)
  {
    DynArray<Image> images;
    images.push_back(DrawArray(lhs, DrawArrayOpts().outputRange()));
    images.push_back(DrawArray(result, DrawArrayOpts().outputRange()));
    Image bigimg;
    DrawImageGrid(bigimg, images);
    bigimg.Write("elliptic_operator_application.png");
  }

  return py::make_tuple(result_time_ms, result_mem_bytes_per_dof);
}



// compute \laplace c - lambda c = -r
py::object testOperatorSolve(const Int3 &size, bool write_image)
{
  ContinuumGrid grid;
  DomainDecomposition mtboxes;
  InitGrid(size, grid, mtboxes);

  Array3df diag_offset(grid.Box(), Cons::DONT);
  Array3dd lhs(grid.Box(), Cons::DONT),
            rhs(grid.Box(), Cons::DONT);


  const float rmax = grid.ld.GetWorldBox().max.norm();
  #pragma omp parallel
  {
    BOOST_FOREACH(const BBox3 &bbox, mtboxes.getCurrentThreadRange())
    {
      FOR_BBOX3(p, bbox)
      {
        float r = grid.ld.LatticeToWorld(p).norm();
        //float loc_rhs = my::smooth_heaviside<float>(rmax*0.2 - r, ld.Scale()*2);
        //rhs[ld.LatticeToSite(p)] = my::sqr((1./sqrt(2.))*ld.LatticeToWorld(p).dot(Float3(1.,-1, 0)));
  //      rhs[ld.LatticeToSite(p)] = std::max(std::abs(ld.LatticeToWorld(p)[0]), std::abs(ld.LatticeToWorld(p)[1]));
  //       diag_offset(p) = -10.;

        float loc_rhs = (p == Int3(20, 20, 0)) ? 1 : 0;
        rhs(p) = loc_rhs; // a minus should be here but the whole equation is multiplied by -1 so both signs cancel
        diag_offset(p) = -1./my::sqr(100*grid.ld.Scale());
        lhs(p) = 0.;
      }
    }
  }

  typedef TrilinosMatrixFreeEllipticOperator<ConstFaceVarFunctor, Array3df> EpetraOperatorType;

//trilinos_linsys_construction.h in already include at this stage
#ifdef EPETRA_MPI
  Epetra_MpiComm epetra_comm(MPI_COMM_SELF);
#endif
#ifndef EPETRA_MPI
  Epetra_SerialComm epetra_comm;
#endif
  Epetra_Map epetra_map((int)grid.ld.NumSites(), 0, epetra_comm);

  boost::scoped_ptr<EpetraOperatorType> op(new EpetraOperatorType(
      grid.ld, grid.dim, mtboxes, -1., ConstFaceVarFunctor(1.), diag_offset, epetra_map));

  Epetra_Vector epetra_rhs(Epetra_DataAccess::View, epetra_map, rhs.getPtr());
  Epetra_Vector epetra_lhs(Epetra_DataAccess::View, epetra_map, lhs.getPtr());

  boost::property_tree::ptree params;
  params.put("output", 1);
  params.put("max_iter", 100);
  params.put("throw", false);

  uint64 mem_usage = GetMemoryUsage().rss;
  my::Time t_;
  
  SolveEllipticEquation(*op, epetra_rhs, epetra_lhs, params);

  double result_time_ms = (my::Time() - t_).to_ms();
  double result_mem_bytes_per_dof = double(GetMemoryUsage().rss_peak - mem_usage)/grid.ld.NumSites();

  if (write_image)
  {
    Array3dd residual(grid.Box(), Cons::DONT);
    Epetra_Vector epetra_residual(Epetra_DataAccess::View, epetra_map, residual.getPtr());
    op->Apply(epetra_lhs, epetra_residual);

    DynArray<Image> images;
    images.push_back(DrawArray(lhs, DrawArrayOpts().outputRange()));
    images.push_back(DrawArray(rhs, DrawArrayOpts().outputRange()));
    images.push_back(DrawArray(residual, DrawArrayOpts().outputRange()));
    Image bigimg;
    DrawImageGrid(bigimg, images);
    bigimg.Write("elliptic_solver_test.png");
  }

  return py::make_tuple(result_time_ms, result_mem_bytes_per_dof);
}



};


void export_elliptic_solver_test()
{
  py::def("EST_testOperatorApplication", EllipticSolverTest::testOperatorApplication);
  py::def("EST_testOperatorSolve", EllipticSolverTest::testOperatorSolve);
}

