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

#include "continuum-surfacetension.h"

template<class T>
void SurfaceTensionForce<T>::init(const BBox3 &bbox_, const Array3d<T> &phi_, int ndim_, T spacing_)
{
  phi = phi_;
  ndim = ndim_;
  spacing = spacing_;
  bcells = bbox_;
  bverts = CellRangeToVertexRange(bcells, ndim);

  Int3 st_cv_offset(0), st;
  for (int i=0; i<ndim; ++i) st_cv_offset[i]=-1;

  normVec.initFromBox(bverts);

  Array3d<T> stencil;
  for (int i=0; i<ndim; ++i)
  {
    stencil = Stencils::Diff1<float>(i, ndim, spacing, Stencils::W2, Stencils::FWD);
    if (i == 0)
    { // since the size of this stencil is 2^dim we can get way with one inde array for all orientations
      stidx_diff_cv = ComputeConvolutionOffsets(BBox3(stencil.getBox()).Move(st_cv_offset), phi.strides());
      stidx_diff_vc = ComputeConvolutionOffsets(stencil.getBox(), normVec.strides());
    }
    stcoeff_diff_cvvc[i] = ComputeConvolutionValues<float>(stencil);
    stencil = Stencils::Diff1<float>(i, ndim, spacing, Stencils::W1, Stencils::CENTERED);
    stidx_diff_cc[i] = ComputeConvolutionOffsets(BBox3(stencil.getBox()).Move(i,-1), phi.strides());
    stcoeff_diff_cc[i] = ComputeConvolutionValues<float>(stencil);
  }

  ComputeNormal();

#if 0
  Image bigimg;
  std::vector<Image> images;
  for (int axis=0; axis<3; ++axis)
    for (int i=0; i<3; ++i)
      images.push_back(DrawArray(Array3d<float>::getComponent(normVec,axis)[bverts], DrawArrayOpts().slice(i, 0).scaled(0.5, 0.5).title(str(format("v_%i (x_%i=0)") % axis % i))));
  DrawImageGrid(bigimg, images);
  bigimg.Write("sft_normals.png");
#endif
}


template<class T>
void SurfaceTensionForce<T>::clear()
{
  Destruct(*this);
  new (this) SurfaceTensionForce<T>;
}


template<class T>
void SurfaceTensionForce<T>::ComputeNormal()
{
  FOR_BBOX3(p, bverts)
  {
    const int p_offset = phi.offset(p);
    Float3 n(0.);
    for (int axis=0; axis<ndim; ++axis)
    {
      // convolute with the difference operator
      for (int i=0; i<stidx_diff_cv.size(); ++i)
        n[axis] += stcoeff_diff_cvvc[axis][i] * phi.dereference(p_offset+stidx_diff_cv[i]);
    }
    float l = n.norm();
    if (l > 0) n *= (1./l);
    normVec(p) = n;
  }
}



template<class T>
float SurfaceTensionForce<T>::computeCurvature(const Int3 &p) const
{
  //FOR_BBOX3(p, bcells)
  {
    const int p_offset = normVec.offset(p);
    float loc_curv = 0.; // div normal
    for (int axis=0; axis<ndim; ++axis)
    {
      float d_a = 0.; // (d normal / d axis)
      for (int i=0; i<stidx_diff_vc.size(); ++i)
        d_a += stcoeff_diff_cvvc[axis][i] * normVec.dereference(p_offset+stidx_diff_vc[i])[axis];
      loc_curv += d_a;
    }
    //curv_(p) = loc_curv;
    return loc_curv;
  }
}


template<class T>
Float3 SurfaceTensionForce<T>::computeForce(const Int3 &p) const
{
  Float3 grad_phi(0.);
  const int p_offset = phi.offset(p);
  for (int axis=0; axis<ndim; ++axis)
  {
    for (int i=0; i<stidx_diff_cc[axis].size(); ++i)
      grad_phi[axis] += stcoeff_diff_cc[axis][i] * phi.dereference(p_offset+stidx_diff_cc[axis][i]);
  }
  float curv = computeCurvature(p);
  float delta = my::smooth_delta_cos<float>(phi.dereference(p_offset), spacing);
  return delta * curv * grad_phi;
}


template<class T>
float SurfaceTensionForce<T>::computeFaceForce(int axis, const Int3 &p) const
{
  Int3 q(p); q[axis]--;
  float phi_qp[2] = { phi(q), phi(p) };
  float curv[2] = {
    computeCurvature(q), // better use a small cache of some sort
    computeCurvature(p)
  };
  float grad = (phi_qp[1]-phi_qp[0])/spacing;
  float avg  = 0.5*(phi_qp[1]+phi_qp[0]);
  float delta = my::smooth_delta_cos<float>(avg, spacing);
  return delta * (curv[0]+curv[1])*0.5 * grad;
}



#define INSTANTIATE(T)\
  template struct SurfaceTensionForce<T>;

INSTANTIATE(float)

  