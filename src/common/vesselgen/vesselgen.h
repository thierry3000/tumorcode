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
#ifndef _VESSELGEN_H_
#define _VESSELGEN_H_

#include "../common.h"
#include <H5Cpp.h>

struct BloodFlowParameters;

namespace VesselGenerator
{

void run(const boost::property_tree::ptree &parameters);
ptree getDefaultParameters();
ptree getVariableArgumentCountParameters();

void vesselgen_generate_grid(H5::Group outgroup, const Int3 &size, float scale, const string &ld_type, float capillary_radius, float press_min, float press_max, const BloodFlowParameters &bfparams);
void vesselgen_generate_single(H5::Group outgroup, const Int3 &size, const int direction_mode, float scale, const string &ld_type, float capillary_radius, float press_min, float press_max, int segment_size, const BloodFlowParameters &bfparams);
void vesselgen_generate_grid_no_flow(H5::Group outgroup, const Int3 &size, float scale, const string &ld_type, float capillary_radius, float press_min, float press_max, const BloodFlowParameters &bfparams);
void vesselgen_generate_symmetric(const H5::Group outgroup, const int &exponent_of_two, float scale, const BloodFlowParameters &bfparams, const bool &only2D);
}

#endif