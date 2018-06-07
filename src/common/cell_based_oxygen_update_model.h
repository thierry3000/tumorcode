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
#ifndef CELL_BASED_OXYGEN_UPTAKE_MODEL_H
#define CELL_BASED_OXYGEN_UPTAKE_MODEL_H

#include <boost/foreach.hpp>

#include "mwlib/field.h"
#include "mwlib/math_ext.h"
#include "mwlib/ptree_ext.h"
#include "continuum-utils.h"

#include <omp.h>

/** analog to growthfactor_model.h
 * 
 * data is only updated, when the change to the previous
 * step is more than 1%.
 * 
 */
struct CellBasedO2Uptake
{
  Array3df fieldLastSources, o2_lut;
  Array3dOps<float> ops;
  const DomainDecomposition *mtboxes;
  const ContinuumGrid *grid;
  
  /*
   * initialize with lattice and parameters
   * params:
   *  rGf = the diffusion radius
   * will compute superposition of cone-shaped decay functions around sources
   * fieldGf is public, the other stuff is private
   */
  void init(const ContinuumGrid &grid, const DomainDecomposition &mtboxes_, const boost::property_tree::ptree &params);
  

  void initField(Array3df &o2uptakefield);
  //void initField();

  void update(Array3df o2uptakefield, ConstArray3d<float> sources);
  
  private:
  struct _conv_addfunc
  {
    float s;
    void inplace( float &a, const float v )
    {
      a += s*v;
    }
  };

};


#endif  //CELL_BASED_OXYGEN_UPTAKE_MODEL_H
