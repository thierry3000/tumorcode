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
#ifndef GROWTHFACTOR_MODEL_H
#define GROWTHFACTOR_MODEL_H

#include <boost/foreach.hpp>

#include "mwlib/field.h"
#include "mwlib/math_ext.h"
#include "mwlib/ptree_ext.h"
#include "continuum-utils.h"

#include <omp.h>

struct GfModel
{
  Array3df fieldLastSources, gf_lut;
  Array3dOps<float> ops;
  const DomainDecomposition *mtboxes;
  
  /*
   * initialize with lattice and parameters
   * params:
   *  rGf = the diffusion radius
   * will compute superposition of cone-shaped decay functions around sources
   * fieldGf is public, the other stuff is private
   */
  void init(const ContinuumGrid &grid, const DomainDecomposition &mtboxes_, const boost::property_tree::ptree &params)
  {
    mtboxes = &mtboxes_;
    ops.init(*mtboxes, grid.Box(), grid.dim, 0);
    ops.init(fieldLastSources, true);

    float scale = grid.ld.Scale();
    float rGf = params.get<float>("rGf");

    int lut_size_max = 2*my::iceil(rGf/(scale))+1;
    Int3 lut_size(lut_size_max);
    if (grid.dim<=2) lut_size[2]=1; // 2 dimensions
    if (grid.dim<=1) lut_size[1]=1; // 1
    
    gf_lut.init(lut_size);
    float sum = 0;
    for( int z=0; z<lut_size[2]; ++z )
    {
      for( int y=0; y<lut_size[1]; ++y )
      {
        for( int x=0; x<lut_size[0]; ++x )
        {
          const Float3 wp = Float3(x-lut_size[0]/2,y-lut_size[1]/2,z-lut_size[2]/2)*scale;
          float val = std::max<float>(0.f,1.0f-norm(wp)/rGf);
          sum += val;
          gf_lut(Int3(x,y,z)) = val;
        }
      }
    }
    sum = 1.0f/sum;
    myAssert(gf_lut.isContiguous());
    gf_lut *= sum;
  }

  void initField(Array3df &gffield)
  {
    ops.init(gffield, true);
  }

  void update(Array3df gffield, ConstArray3d<float> sources)
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
                gffield,
                gf_lut,
                f,
                p-(gf_lut.size()/2),
                0
              );
            }
            fieldLastSources(p) = src;
          }
        }
      }
    }
  }


private:
  struct _conv_addfunc
  {
    float s;
    void inplace( float &a, const float v )
    {
      a += s*v;
    }
  };

//   void add(const Int3& pos, float src_value)
//   {
//     _conv_addfunc f; f.s = src_value;
//     array_apply_offset3d<float, float,_conv_addfunc>(
//       fieldGf,
//       gf_lut,
//       f,
//       pos-(gf_lut.size()/2),
//       0
//     );
//   }
};

struct GfModel_Cell
{
  Array3df fieldLastSources, gf_lut;
  Array3dOps<float> ops;
  const DomainDecomposition *mtboxes;
  
  /*
   * initialize with lattice and parameters
   * params:
   *  rGf = the diffusion radius
   * will compute superposition of cone-shaped decay functions around sources
   * fieldGf is public, the other stuff is private
   */
  void init(const ContinuumGrid &grid, const DomainDecomposition &mtboxes_, const boost::property_tree::ptree &params)
  {
    mtboxes = &mtboxes_;
    ops.init(*mtboxes, grid.Box(), grid.dim, 0);
    ops.init(fieldLastSources, true);

    float scale = grid.ld.Scale();
    float rGf = params.get<float>("rGf");

    int lut_size_max = 2*my::iceil(rGf/(scale))+1;
    Int3 lut_size(lut_size_max);
    if (grid.dim<=2) lut_size[2]=1; // 2 dimensions
    if (grid.dim<=1) lut_size[1]=1; // 1
    
    gf_lut.init(lut_size);
    float sum = 0;
    for( int z=0; z<lut_size[2]; ++z )
    {
      for( int y=0; y<lut_size[1]; ++y )
      {
        for( int x=0; x<lut_size[0]; ++x )
        {
          const Float3 wp = Float3(x-lut_size[0]/2,y-lut_size[1]/2,z-lut_size[2]/2)*scale;
          float val = std::max<float>(0.f,1.0f-norm(wp)/rGf);
          sum += val;
          gf_lut(Int3(x,y,z)) = val;
        }
      }
    }
    sum = 1.0f/sum;
    myAssert(gf_lut.isContiguous());
    gf_lut *= sum;
  }

  void initField(Array3df &gffield)
  {
    ops.init(gffield, true);
  }

  void update(Array3df gffield, ConstArray3d<float> sources)
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
              //static void array_apply_offset3d(Array3d<T> &a, const ConstArray3d<S> &b, F &func, const Int3 &boffset, uchar flags=0)
              _conv_addfunc f; f.s = diff;
              array_apply_offset3d<float, float,_conv_addfunc>(
                gffield,// Array3d<T>
                gf_lut, // ConstArray3d<T>
                f,// F func
                p-(gf_lut.size()/2),// Int3 boffset
                0// flags
              );
            }
            fieldLastSources(p) = src;
          }
        }
      }
    }
  }
void AddSourceDistribution_from_cell( std::vector<double> x,
                                      std::vector<double> y,
                                      std::vector<double> z,
                            const BBox3 &bbox,
                            const LatticeDataQuad3d &field_ld, 
                            int dim, 
                            Array3d<float> l_coeff, 
                            const ptree &params);

private:
  struct _conv_addfunc
  {
    float s;
    void inplace( float &a, const float v )
    {
      a += s*v;
    }
  };

//   void add(const Int3& pos, float src_value)
//   {
//     _conv_addfunc f; f.s = src_value;
//     array_apply_offset3d<float, float,_conv_addfunc>(
//       fieldGf,
//       gf_lut,
//       f,
//       pos-(gf_lut.size()/2),
//       0
//     );
//   }

};


#endif
