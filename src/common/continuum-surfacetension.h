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
#ifndef CONTINUUM_SFT_H
#define CONTINUUM_SFT_H

#include "common.h"
#include "continuum-utils.h"
#include "continuum-grid.h"

/* The surface tension force is computed for cell centers, i.e.
 * Fsf = delta(phi)kappa(x)*n,
 * which is approximated as smooth_delta(phi)*[div (grad phi / |grad phi|)] grad phi.
 * 
 * Phi gives the signed distance to the interface (|grad phi| = 1)
 * bbox is a sub-region for which the results will be available, because
 * internal buffers for intermediate results are used.
 *
 * The difference scheme computes the normals (grad phi / |grad phi|) on grid vertices,
 * and their divergence again on cell centers.
*/
template<class T>
class SurfaceTensionForce
{
  BBox3 bcells, bverts;
  int ndim;
  float spacing;
  Array3d<Float3> normVec; // vertex centered
  Array3d<T> phi; // distance function

  DynArray<int> stidx_diff_cv, stidx_diff_vc, stidx_diff_cc[3];
  DynArray<float> stcoeff_diff_cvvc[3], stcoeff_diff_cc[3];
  void ComputeNormal();
  
public:
  void init(const BBox3 &bbox, const Array3d<T> &phi, int ndim, T spacing);
  void clear();

  float computeCurvature(const Int3 &p) const;
  Float3 computeForce(const Int3 &p) const;
  // on cell faces; p is the face index.
  float computeFaceForce(int axis, const Int3 &p) const;
};

#endif
  