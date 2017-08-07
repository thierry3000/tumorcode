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
#include "vesselgen/vesselgen.h"
#include "calcflow.h"
#include <boost/property_tree/info_parser.hpp>

#include <fenv.h>
#include "mwlib/hdf_wrapper_ptree.h"
#include "mwlib/log.h"

/** @brief Vesselgenerator stuff
 */
void run_vesselgen(const py::str &param_info_str)
{
  ptree pt_params = convertInfoStr(param_info_str, ptree());
  VesselGenerator::run(pt_params);
}
void export_vesselgen()
{
  py::def("run_vesselgen", run_vesselgen);
  //FIX ME: vesselgen_generate_* got a new parameter which needs to be exported properly
  py::def("vesselgen_generate_grid", VesselGenerator::vesselgen_generate_grid);
  py::def("vesselgen_generate_single", VesselGenerator::vesselgen_generate_single);
  py::def("vesselgen_generate_symmetric", VesselGenerator::vesselgen_generate_symmetric);
}
