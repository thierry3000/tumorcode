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
#include "iff_drug_app.h"

#include "../python_krebsutils/python-helpers.h"
#include "numpy.hpp"
#include "../common/calcflow.h"
#include <algorithm>

#define BOOST_RESULT_OF_USE_DECLTYPE 1

namespace py = boost::python;
namespace nm = boost::python::numeric;
namespace h5 = h5cpp;

/**
 * @brief Calculates Interstistial pressure, flow and drug stuff
 * 
 * 1.Welter, M. & Rieger, H. 
 * Interstitial Fluid Flow and Drug Delivery in Vascularized Tumors: A Computational Model. 
 * PLoS ONE 8, e70395 (2013).
 */
namespace Iff
{
/** @brief solving python2 interface problems with bools */
bool IsTrue(bool pbool){return pbool;}
/** @brief main
 * taking care of the simulation
 * \param param_info_str prams from the python dict piped into string from boostpython BOOST_PYTHON_MODULE
 * \param outfilename filename for the simulation results
 * \param drug_measurement_function a python class measuring certain extra stuff during runtime (quite adavanced ;-)
 */
void run_iffsim(const py::str &param_info_str, const string &outfilename, py::object &drug_measurement_function)
{
  FpExceptionStateGuard exception_state_guard(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
  ptree pt_params;
  {
    string param_string = py::extract<string>(param_info_str);
    std::istringstream param_stream(param_string);
    boost::property_tree::read_info(param_stream, pt_params);
  }
  IffDrugApp3d app;
  app.Main(pt_params, outfilename, drug_measurement_function);
}



void export_iffsim()
{
  py::def("run_iffsim", run_iffsim);
}

}//namespace
#ifdef DEBUG
BOOST_PYTHON_MODULE(libiff_d)
#else
BOOST_PYTHON_MODULE(libiff_)
#endif
{
  Iff::export_iffsim();
}
