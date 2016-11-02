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
#ifndef FD_STENCILS_H
#define FD_STENCILS_H

#include "continuum-grid.h"
#include "mwlib/field.h"

inline PdeReal CalcInterfaceCondCoeff(double k0, double k1) { return 2.0*(k0*k1)/(k0+k1+1.e-30); }

/**
@brief Finite difference operators for regular grids
*/
namespace Stencils
{

// Width determines the size and coefficients perpendicular to the differentiation axis.
// It is 1 site for W1, 2^dim for W2, 3^dim for W3.
enum Width
{
  W1,
  W2,
  W3,
  W3_ISOTROPIC,
  W3_SPARSE,
};

// Style determines the length of the filter along the differentiation axis.
// It is 2 cells for FWD (standard forward difference), and 3 for CENTERED.
enum Style
{
  FWD,
  CENTERED
};


template<class T>
const Array3d<T> Diff1(int ax, int ndim, const T spacing, Width width, Style style);

template<class T>
const Array3d<T> Diff2(int ax, int ndim, const T spacing, Width width);

template<class T>
const Array3d<T> Laplace(int ndim, const T spacing, Width width);

template<class T>
const Array3d<T> Box(int ax, int ndim, Width width, int l);


// from function values f_i, and filter stencil values a_k, b_k, k=0..2
// build new filter c that corresponds to
// c x f = a x b x f = sum_i=0..2 sum_j=0..2 a_i * b_j * f_i+j
template<class T>
void Concatenate(const ConstArray3d<T> a, const ConstArray3d<T> b, Array3d<T> out);

template<class T> 
void Concatenate(const ConstArray3d<T> a, const ConstArray3d<T> values, const ConstArray3d<T> b, Array3d<T> out);

template<class T>
inline Array3d<T> Concatenate(const ConstArray3d<T> a, const ConstArray3d<T> b)
{
  Array3d<T> c(a.size()+b.size()-Int3(1));
  Concatenate(a, b, c);
  return c;
}

template<class T>
inline Array3d<T> Align(const ConstArray3d<T> &a, int axis)
{
  ConstArray3d<T> t(a);
  t.swapAxes(0, axis);
  return Array3d<T>(t, Cons::COPY);
}


} // namespace


// build 3x3x3 filters
void BuildLaplaceFilterQuad3d( Array3d<float> &r );
void BuildDiff1FilterQuad3d( Array3d<float> &b, int ax );
void BuildDiff2FilterQuad3d( Array3d<float> &b, int ax );

// convolution, make sure there is a apropriate region around the site p
template<class T,class S>
S Convolve3x3x3( const ConstArray3d<T> &filter, const ConstArray3d<S> &a, const Int3 &p );

/*------------------------------------------------------
------------------------------------------------------*/

template<class T>
void AddSmoothDelta(Array3d<T> &arr, const BBox3 &bbox, const LatticeDataQuad3d &ld, int dim, const Float3 &pos, T weight);

/*------------------------------------------------------
------------------------------------------------------*/

template<class T>
class Array3dOps
{
  const DomainDecomposition *mtboxes;
  BBox3 bbox, bbext;
  int border, dim;

public:
  typedef Array3d<T> state_type;
  Array3dOps() : mtboxes(NULL), border(0) {}
  Array3dOps(const DomainDecomposition &mtboxes_, const BBox3 &bbox_, int dim_, int border_) : mtboxes(&mtboxes_), bbox(bbox_), dim(dim_), border(border_)
  {
    bbext = ExtendForDim(bbox, dim, border);
  }
  void init(const DomainDecomposition &mtboxes_,  const BBox3 &bbox_, int dim_, int border_)
  {
    Destruct(*this);
    new (this) Array3dOps<T>(mtboxes_, bbox_, dim_, border_);
  }
  // perform allocation and multithreaded (ccNUMA aware) clearing
  void init(state_type &u, bool clear = true) const;
  void addScaled(double fa, state_type &u, double fb, const state_type &v) const;
  void initFrom(state_type &u, const state_type &other, ConsMode mode) const;
};


#endif