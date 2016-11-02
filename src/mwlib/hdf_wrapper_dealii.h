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
#ifndef HDF_WRAPPER_DEALII_H
#define HDF_WRAPPER_DEALII_H

#include <exception>
#include <stdexcept>

#include <grid/grid_out.h>
#include <grid/grid_in.h>
#include <grid/persistent_tria.h>

#include "hdf_wrapper_stl.h"



template<int dim>
static void write_triangulation(h5cpp::Group gparent, const std::string &name, const dealii::Triangulation<dim> &coarse_triangulation, const dealii::PersistentTriangulation<dim> &triangulation)
{
  using namespace dealii;
  using namespace my;

  Assert(coarse_triangulation.n_levels() == 1, ExcInvalidState());
  
  h5::Group g = gparent.create_group("tria");
  h5::Attributes a = g.attrs();
  a.set<std::string>("type", "deal.ii refined mesh");
  a.set<int>("dim", dim);
  a.set<std::string>("format", "msh");

  {
    GridOut grid_out;
    grid_out.set_flags(GridOutFlags::Msh(true, true));
    std::ostringstream f;
    grid_out.write_msh(coarse_triangulation, f);
    h5::create_dataset_range (g, "mesh", f.str());
    //cout << "wrote mesh" << endl << f.str() << endl;
  }

  {
    std::ostringstream f;
    triangulation.write_flags(f);
    h5::create_dataset_range(g, "refineflags", f.str());
    //cout << "wrote flags" << endl << f.str() << endl;
  }
}

template<int dim>
static void read_triangulation(h5cpp::Group gparent, const std::string &name, dealii::Triangulation<dim> &coarse_triangulation, dealii::PersistentTriangulation<dim> &triangulation)
{
  using namespace dealii;
  using namespace my;
  
  h5::Group g = gparent.open_group("tria");
  h5::Attributes a = g.attrs();
  if (a.get<std::string>("type") != "deal.ii refined mesh")
    throw std::runtime_error("wrong tria type string");
  if (a.get<int>("dim") != Triangulation<dim>::dimension)
    throw std::runtime_error("wrong space dimension");
  if (a.get<std::string>("format") != "msh")
    throw std::runtime_error("wrong mesh format");

  std::string data;
  {
    h5::Dataset meshds = g.open_dataset("mesh");
    h5::read_dataset_stl_resizing(meshds, data);
    //cout << "read mesh" << endl << data << endl;
    std::istringstream f(data);
    data.clear();

    GridIn<dim> grid_in;
    grid_in.attach_triangulation(coarse_triangulation);
    grid_in.read_msh(f);
  }
  {
    h5::Dataset flagds = g.open_dataset("refineflags");
    h5::read_dataset_stl_resizing(flagds, data);
    //cout << "read flags" << endl << data << endl;
    std::istringstream f(data);
    data.clear();

    triangulation.read_flags(f);
    triangulation.restore();
  }
}


#endif
