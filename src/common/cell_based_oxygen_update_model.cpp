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
#include "cell_based_oxygen_update_model.h"

/** the source should be at this stage,
 * here we add how this source propergates some GF
 */

void CellBasedO2Uptake::init(const ContinuumGrid& grid_, const DomainDecomposition& mtboxes_, const boost::property_tree::ptree& params)
{
  mtboxes = &mtboxes_;
  grid = &grid_;
  ops.init(*mtboxes, grid_.Box(), grid_.dim, 0);
  ops.init(fieldLastSources, true);
  
  float scale = grid->ld.Scale();
  float rO2Consumtion = params.get<float>("rO2Consumtion");

  int lut_size_max = 2*my::iceil(rO2Consumtion/(scale))+1;
  Int3 lut_size(lut_size_max);
  if (grid->Dim()<=2) lut_size[2]=1; // 2 dimensions
  if (grid->Dim()<=1) lut_size[1]=1; // 1
  
  o2_lut.init(lut_size);
  float sum = 0;
  for( int z=0; z<lut_size[2]; ++z )
  {
    for( int y=0; y<lut_size[1]; ++y )
    {
      for( int x=0; x<lut_size[0]; ++x )
      {
        const Float3 wp = Float3(x-lut_size[0]/2,y-lut_size[1]/2,z-lut_size[2]/2)*scale;
        float val = std::max<float>(0.f,1.0f-norm(wp)/rO2Consumtion);
        sum += val;
        o2_lut(Int3(x,y,z)) = val;
      }
    }
  }
  sum = 1.0f/sum;
  myAssert(o2_lut.isContiguous());
  o2_lut *= sum;
}

//void CellBasedO2Uptake::initField(Array3df& o2uptakefield)
void CellBasedO2Uptake::initField(Array3df &o2uptakefield)
{
  ops.init(o2uptakefield, true);
}

void CellBasedO2Uptake::update(Array3df o2uptakefield, ConstArray3d<float> sources)
{
  #pragma omp parallel
    {
      BOOST_FOREACH(const BBox3 &bbox, mtboxes->getCurrentThreadRange())
      {
        FOR_BBOX3(p, bbox)
        {
          float last = fieldLastSources(p);
          float src = sources(p);
          float diff = src-last;
          if(std::abs(diff)>1.0e-2)
          {
            {
              _conv_addfunc f; f.s = diff;
              array_apply_offset3d<float, float,_conv_addfunc>(
                o2uptakefield,
                o2_lut,
                f,
                p-(o2_lut.size()/2),
                0
              );
            }
            fieldLastSources(p) = src;
          }
        }
      }
    }
}
