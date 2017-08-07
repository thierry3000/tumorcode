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
#ifndef OXYGEN_MODEL_H
#define OXYGEN_MODEL_H

//#include "mwlib/field.h"
//#include "mwlib/math_ext.h"
#include "vessels3d.h"
//#include "mwlib/ptree_ext.h"
#include "common.h"

namespace O2Model
{

/*
 * add source distribution of vessels
 * params:
 *  capillary_wall_permeability
 *  hematocrit_init: vessel hematocrit is divided by this to obtain the vessel o2 level
 * 
 * vessel permeability scales down as vessel radius increases
 *
 * it needs the parameters
 *   capillary_wall_permeability
 *   hematocrit_init
 * the blood o2 level is computed as v->hematocrit/hematocrit_init
 */
struct SimpleO2Params
{
  SimpleO2Params();
  void assign(const ptree &pt);
  ptree as_ptree() const;
  static void update_ptree(ptree &dst, const ptree &src);
  
  double
         o2_range[3],
         o2_cons_coeff[3],
         o2_rel_tumor_source_density,
         o2_level_normal,
         hematocrit_init,
         reference_intercapillary_distance,
         capillary_wall_permeability;
  int test_obstacle;
  bool use_o2_source_decay;
};

// void AddSourceDistribution(Array3d<float> clinear_field, Array3d<float> rhs_field,
//                                  const LatticeDataQuad3d &field_ld,
//                                  int dim,
//                                  const VesselList3d &vl,
//                                  const boost::property_tree::ptree &params);
void AddSourceDistribution(const BBox3 &bbox,
                           const LatticeDataQuad3d &field_ld,
                           int dim,
                           Array3d<float> l_coeff,
                           Array3d<float> rhs,
                           const VesselList3d::LatticeData &vessld,
                           const DynArray<const Vessel*> &vessels,
                           const ptree &params);

void assignRangeParam(double &range, double &cons_coeff, const string &parameter_postfix, const ptree &pt);
void assignWallPermeabilityCoefficient(double &capillary_wall_permeability, double &o2_level_normal, double consumption_coeff, double reference_intercapillary_distance, double capillary_radius, const ptree &pt);

/*
  c = consumption coefficient = 1./Rd^2
  q = vessel wall permeability
  vess_rad = radius
  vess_spacing = intercapillary distance
*/
double CalcHomogeneousMaxOxy(double c, double q, double vess_rad, double vess_spacing);
double CalcHomogeneousCoeffOxy(double c, double o2_level_normal, double vess_rad, double vess_spacing);

}//end namespace
#endif
