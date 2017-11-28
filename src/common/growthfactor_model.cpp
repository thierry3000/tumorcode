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
#include "growthfactor_model.h"

void GfModel_Cell::AddSourceDistribution_from_cell( std::vector<double> x,
                                      std::vector<double> y,
                                      std::vector<double> z,
                            const BBox3 &bbox,
                            const LatticeDataQuad3d &field_ld, 
                            int dim, 
                            Array3d<float> l_coeff, 
                            Array3d<float> rhs, 
                            const ptree &params)
{
//   //const double capillary_wall_permeability = params.get<double>("capillary_wall_permeability");
//   const double coeffBloodGlucose   = 1.0f;
// 
//   const double maturation_at_r5 = GetInitialThickness(5.0f);
//   //const double maturation_at_r20= GetInitialThickness(20.0f);
//   //const double hematocrit_init = params.get<double>("hematocrit_init");
// /** 
//  * In 2d, the samples are projected on the xy plane.
//  * It is assumed that a 2d slice is 200 micron in height. Hence the scaling factor. 
//  * 
//  * 3D case this is just 1.
//  */
//   const double dim_coeff = dim==2 ? (field_ld.Scale()/params.get<double>("fake_height_2d",200.)) : 1.;
//   
//   CylinderNetworkSampler sampler; sampler.Init(field_ld.Scale(), params);
//   //each thread could execute its own loop independently
//   for (int i=0; i<vessels.size(); ++i)
//   {
//     const Vessel* v = vessels[i];
//     if( !v->IsCirculated() ) continue;
//     
//     //float o2level = coeffBloodOxy * v->hematocrit / hematocrit_init;
//     float glucoseLevel = coeffBloodGlucose;
// //     myAssert(v->hematocrit > 0.);
//     //float perm = capillary_wall_permeability * maturation_at_r5 / std::max<float>(maturation_at_r20, v->maturation);  // maturation has the value of vessel wall thickness, have to limit it so the coefficient does not become infinite
    float cell_perm = 1e-10;
    float loc_l_coef = cell_perm*2;
    float loc_rhs = cell_perm*4;
    loc_l_coef = 1;
    loc_rhs = 5;
  //     //float perm = capillary_wall_permeability * 1./std::max<float>(1., v->maturation / maturation_at_r5);
// //     float perm = capillary_wall_permeability;
//     float perm = 1.0;
// //     myAssert(v->maturation > 0.);
//     float loc_l_coef = -perm*dim_coeff;
//     float loc_rhs = -perm*glucoseLevel*dim_coeff;
// 
//     myAssert(std::isfinite(loc_l_coef) && std::isfinite(loc_rhs));
// 
//     sampler.Set(vessld.LatticeToWorld(v->LPosA()), vessld.LatticeToWorld(v->LPosB()), v->r);
//     int cnt = sampler.GenerateSurfaceSamples();
//     for (int j=0; j<cnt; ++j)
//     {
//       const CylinderNetworkSampler::Sample &s = sampler.GetSample(j);
// 
//       AddSmoothDelta(l_coeff, bbox, field_ld, dim, s.wpos, sampler.weight_per_volume*loc_l_coef);
//       AddSmoothDelta(rhs, bbox, field_ld, dim, s.wpos, sampler.weight_per_volume*loc_rhs);
//     }
//   }
  printf("lets design the GF source distribution\n");
  for(int i=0; i<x.size();++i)
  {
    float offset = 10*15;
    offset = 0.0;
    Float3 pos(x[i]+offset,y[i]+offset,z[i]+offset);
    AddSmoothDelta(l_coeff, bbox, field_ld, dim, pos, loc_l_coef);
    AddSmoothDelta(rhs, bbox, field_ld, dim, pos, loc_rhs);
  }
  
}
