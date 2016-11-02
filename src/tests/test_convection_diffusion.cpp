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

#include "convection_diffusion_solver.h"
#include "mwlib/ptree_ext.h"

#include <boost/format.hpp>
#include "mwlib/timer.h"

using boost::str;
using boost::format;


int main(int argc, char **argv)
{
  return 0;
}

#if 0
  const int dim = 2;
  LatticeDataQuad3d ld;
  ld.Init(Int3(100, 100, 1), 1.);
  ld.SetWorldPosition(-0.5*ld.GetWorldBox().max);
  LatticeIndexRanges ir = LatticeIndexRangesFromCellRange(ld.Box(), dim);

  double p_velocity = 10.;
  double p_kdiff = 2;
  double p_tend = 5;
  double p_source_strength = 50.;
  
  boost::array<Array3d<float>, 3> velocity_fields;
  for(int axis=0; axis<dim; ++axis)
  {
    velocity_fields[axis].initFromBox(ir.faces[axis]);
  }
  velocity_fields[0].fill(p_velocity);

  Array3d<float> conc(ExtendForDim(ir.cells, dim, 2));
  Array3d<float> lc_field(ir.cells), cc_field(ir.cells);
  Array3d<float> kdiff(ExtendForDim(ir.cells, dim, 1), p_kdiff);
  conc.setBox(ir.cells);
  FOR_BBOX3(p, ir.cells)
  {
    Float3 wp = ld.LatticeToWorld(p);
    float r = wp.norm();
    float val = my::smooth_heaviside<float>(5 - r, ld.Scale());
    //conc(p) = val;
    //cc_field(p) = val;
    cc_field(p) = val * p_source_strength;
    lc_field(p) = -val * p_source_strength;
  }

  //double max_dt = 0.75 * ld.Scale() / (dim*std::max(p_velocity, 1.));
  //double max_dt = 0.1;

  typedef Array3d<float> AT;
  typedef ConvectionDiffusionStepper<AT, AT, AT, AT, float> Stepper;
  Stepper stepper; stepper.init(ld, dim);

  my::Time rt_;
  
  int outnum = 0, iteration_num=0;
  double t = 0, next_out_t = 0.;
  while (true)
  {
    if (t >= next_out_t) {
      const string fn = str(format("convdiff-imexbdf-%02i.png") % outnum);
      cout << format("t = %f -> %s") % t % fn << endl;
      Image img;
      DrawArray(img, conc[ir.cells], DrawArrayOpts().outputRange());
      img.Write(fn);
      ++outnum;
      next_out_t += 1.;
    }

    if (t > p_tend) break;

//     FOR_BBOX3(p, ir.cells)
//     {
//       lc_field(p) = p_sinkcoeff * (conc(p) + p_sinkc05);
//     }

    boost::property_tree::ptree params = boost::property_tree::make_ptree("output", 1);
    double dt = stepper.doStep( 1.,
                                conc,
                                velocity_fields,
                                kdiff,
                                lc_field,
                                cc_field,
                                params);

    t += dt;
    ++iteration_num;
  }

  cout << "runtime = " << (my::Time() - rt_) << endl;

}
#endif