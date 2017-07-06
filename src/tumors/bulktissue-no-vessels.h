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
#pragma once // include this file only once per compilation unit (see https://en.wikipedia.org/wiki/Pragma_once)

#include <fenv.h>

#include "common.h"
#include "time_stepper_utils_new.h"
#include "continuum-stepper-states.h"
#include "mwlib/log.h"
#include "trilinos_linsys_construction.h"
#include "bulktissuemodel1_new.h"
#include "vesselmodel1.h"//needed for parameter handling!!!!
#include "shared-objects.h"
#include "oxygen_model.h"
#include <boost/property_tree/info_parser.hpp>
#ifdef USE_ADAPTION
  #include "adaption/adaption_model2.h"
#endif

namespace BulkTissueWithoutVessels
{

struct State
{
  NewBulkTissueModel::State tumor_state;
  Array3df o2sources, o2field;
};

/*------------------------------
 *    parameter handling
 * -----------------------------*/

enum TissueId {
  TCS = 0,
  DEAD = 1,
  TISSUE = 2,
};

/*
parameters:
lattice_size, Int3
lattice_scale, float
tend, float
out_intervall, float
save_hdf, bool,
save_images, bool,
fn_out, string
*/
//this goes now along with tum-only-vessels
struct SimulationParameters
{
  Int3 lattice_size;
  Int3 mtboxsize;
  double lattice_scale;
  double out_intervall;
  int num_threads;
  double tend;
  bool save_hdf;
  bool safe_images;
  string fn_out,fn_vessel;
  string paramset_name;
  SimulationParameters();
  void assign(const ptree &pt);
  ptree as_ptree() const;
  static void update_ptree(ptree &dst, const ptree &src);
};

void run(const ptree &params);

}//end namespace BulkTissueWithoutVessels
