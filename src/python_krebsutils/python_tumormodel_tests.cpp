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
//#include "numpy.hpp"
#include "bulktissuemodel1_new.h"

#include <boost/property_tree/info_parser.hpp>

/* WARNING: M Welter switched the header from bulktissuemodel1 to bulktissuemodel1_new without testing if it still works.
 * The reason is the removal of the legacy file. */

np::arraytbase calcBulkTissueSourceTerm(const py::list &arg_list, const py::str param_info_str)
{
  /*
   * arg_list is a list of state variables for the source term
   * param_info_str is a boost .info file as string with the same parameters
   * as in the full simulation
   */
  string param_string = py::extract<string>(param_info_str);
  std::istringstream param_stream(param_string);
  ptree pt_tmp, pt_params = NewBulkTissueModel::Params().as_ptree();
  boost::property_tree::read_info(param_stream, pt_tmp);
  boost::property_tree::update(pt_params, pt_tmp);
  NewBulkTissueModel::Params params; params.assign(pt_params);
  params.init(LatticeDataQuad3d(BBox3(0,0,0,1,1,1), 30.));
  
  Py_ssize_t n = py::len(arg_list);
  np::arrayt<float> res = np::zeros(1, &n, np::getItemtype<float>());
  
  for (int i=0; i<n; ++i)
  {
    const py::tuple a = py::extract<py::tuple>(arg_list[i]);
    float theta = py::extract<float>(a[0]);
    float phi_cells = py::extract<float>(a[1]);
    float phi_other = py::extract<float>(a[2]);
    float phi_necro = py::extract<float>(a[3]);
    float o2 = py::extract<float>(a[4]);
    float rates[2]; // 0 = cells+necro, 1 = necro
    NewBulkTissueModel::calcSource(theta, phi_cells, phi_other, phi_necro, o2, params, rates);
    res[i] = rates[0];
  }
  return res;
}


void export_model_helpers()
{
  py::def("calcBulkTissueSourceTerm", calcBulkTissueSourceTerm);
}
