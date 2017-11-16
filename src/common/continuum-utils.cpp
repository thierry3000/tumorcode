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

#include "continuum-utils.h"
#include "levelset.h"
#include <boost/foreach.hpp>

namespace Stencils
{


template<class T>
const Array3d<T> ImplDiff(int ax, int ndim, int order, T spacing, Width width, Style style)
{
  typedef Vec<T,3> VecT;
  int w, l = (style==FWD) ? 2 : 3;
  VecT coeff;
  switch(width)
  {
  case W3_ISOTROPIC:
    if(order == 1 && style==CENTERED)
    {
      coeff = VecT(1./6., 4./6., 1./6.);
    }
    else
    {
      coeff = VecT(1./12.,10./12.,1./12.);
    }
    w = 3;
    break;
  case W3_SPARSE:
    coeff = VecT(0., 1., 0.);
    w = 3;
    break;
  case W3:
    coeff = VecT(1./3., 1./3., 1./3.);
    w = 3;
    break;
  case W2:
    coeff = VecT(1./2., 1./2., 0.);
    w = 2;
    break;
  default:
    coeff = VecT(1., 0., 0.);
    w = 1;
  }
  Array3d<T> d;
  d.init(Int3(l, ndim>1 ? w : 1, ndim>2 ? w : 1), 0);
  if(order == 1)
  {
    T f = T(1.)/(l-1);
    d(l-1, 0, 0) = f;
    d(0, 0, 0) = -f;
  }
  else
  {
    myAssert(l == 3);
    d(0, 0, 0) = d(2, 0, 0) = 1.;
    d(1, 0, 0) = -2.;
  }
  if(ndim > 1)
  {  
    Array3d<T> t(d.size());
    for(int i=0; i<w; ++i) for(int j=0; j<l; ++j)
    {
      t(j,i,0) = coeff[i] * d(j,0,0);
    }
    d = t;
  }
  if(ndim > 2)
  {
    Array3d<T> t(d.size());
    for(int i=0; i<w; ++i) for(int j=0; j<w; ++j) for(int k=0; k<l; ++k)
    {
      t(k,j,i) = coeff[i] * d(k,j,0);
    }
    d = t;
  }
  if(spacing != 1.)
  {
    d *= (1./(order==1 ? spacing : spacing*spacing));
  }
  if(ax != 0)
  {
    d.swapAxes(0, ax);
    return Array3d<T>(d, Cons::COPY);
  }
  else return d;
}


template<class T>
const Array3d<T> Diff1(int ax, int ndim, const T spacing, Width width, Style style)
{
  return ImplDiff(ax, ndim, 1, spacing, width, style);
}

template<class T>
const Array3d<T> Diff2(int ax, int ndim, const T spacing, Width width)
{
  return ImplDiff(ax, ndim, 2, spacing, width, Stencils::CENTERED);
}

template<class T>
const Array3d<T> Laplace(int ndim, const T spacing, Width width)
{
  myAssert(width==Stencils::W3 || width==Stencils::W3_ISOTROPIC || width==Stencils::W3_SPARSE);
  Array3d<T> a = Diff2(0, ndim, spacing, width);
  //printf("--\n");
  if(ndim>1)
    a += Diff2(1, ndim, spacing, width);
  //printf("--\n");
  if(ndim>2)
    a += Diff2(2, ndim, spacing, width);
  //printf("--\n");
  return a;
}

template<class T>
const Array3d<T> Box(int ax, int ndim, Width width, int l)
{
  Array3d<T> d;
  int w = 1;
  switch(width)
  {
  case W3_ISOTROPIC:
  case W3_SPARSE:
  case W3:
    w = 3;
    break;
  case W2:
    w = 2;
    break;
  default:
    w = 1;
  }
  d.init(Int3(l, ndim>1 ? w : 1, ndim>2 ? w : 1), 0);
  d.fill(1./d.size().prod());
  if(ax != 0)
  {
    d.swapAxes(0, ax);
    return Array3d<T>(d, Cons::COPY);
  }
  else return d;
}


template<class T>
void Concatenate(const ConstArray3d<T> a, const ConstArray3d<T> b, Array3d<T> out)
{
  // from function values f_i, and filter stencil values a_k, b_k, k=0..2
  // build new filter c that corresponds to  
  // c x f = a x b x f = sum_i=0..2 sum_j=0..2 a_i * b_j * f_i+j
  out.fill(0);
  myAssert(out.size() == a.size() + b.size()-Int3(1));
  FOR_REG3V2(i, a.size())
  {
    FOR_REG3V2(j, b.size())
    {
      out(i+j) += a(i) * b(j);
    }
  }
}

template<class T> 
void Concatenate(const ConstArray3d<T> a, const ConstArray3d<T> values, const ConstArray3d<T> b, Array3d<T> out)
{
  // from function values f_i, and filter stencil values a_k, b_k, k=0..2
  // build new filter c that corresponds to  
  // c x f = a x values * b x f = sum_i=0..2 sum_j=0..2 a_i * values_i * b_j * f_i+j
  // for example to build the operator \nabla \cdot (values \nabla)
  myAssert(values.size() == a.size());
  myAssert(out.size() == a.size() + b.size()-Int3(1));
  //Array3d<T> c(a.size() + b.size()-Int3(1));
  out.fill(0);
  FOR_REG3V2(i, a.size())
  {
    FOR_REG3V2(j, b.size())
    {
      out(i+j) += a(i) * values(i) * b(j);
    }
  }
}

#define INSTANTIATE(T)\
    template const Array3d<T> Diff1(int ax, int ndim, const T spacing, Width width, Style style);\
    template const Array3d<T> Diff2(int ax, int ndim, const T spacing, Width width);\
    template const Array3d<T> Laplace(int ndim, const T spacing, Width width);\
    template void Concatenate(const ConstArray3d<T> a, const ConstArray3d<T> b, Array3d<T> out);\
    template void Concatenate(const ConstArray3d<T> a, const ConstArray3d<T> values, const ConstArray3d<T> b, Array3d<T> out); \
    template const Array3d<T> Box(int ax, int ndim, Width width, int l);

INSTANTIATE(float)
INSTANTIATE(double)

} // namespace stencils





#define RESTRICT
/*------------------------------------------------------
   Quad 3D convolution filters
------------------------------------------------------*/

template<class T,class S>
S Convolve3x3x3( const ConstArray3d<T> &filter, const ConstArray3d<S> &a, const Int3 &p )
{
  myAssert(filter.size()==Int3(3));
  const T* RESTRICT pf = filter.getPtr();
  const Int3 m = a.strides();
  float r = 0;
  Int3 q;
  for(q[2]=p[2]-1; q[2]<=p[2]+1; ++q[2])
  {
    for(q[1]=p[1]-1; q[1]<=p[1]+1; ++q[1])
    {
      const S* RESTRICT pa = a.getPtr() + a.offset(p[0]-1,q[1],q[2]);
      r += (*pf)*(*pa);
      ++pf; pa+=m[0];
      r += (*pf)*(*pa);
      ++pf; pa+=m[0];
      r += (*pf)*(*pa);
      ++pf;
    }
  }
  return r;
}

template float Convolve3x3x3<float,float>( const ConstArray3d<float> &filter, const ConstArray3d<float> &a, const Int3 &p );
template double Convolve3x3x3<float,double>( const ConstArray3d<float> &filter, const ConstArray3d<double> &a, const Int3 &p );
template uchar Convolve3x3x3<float,uchar>( const ConstArray3d<float> &filter, const ConstArray3d<uchar> &a, const Int3 &p );
template int Convolve3x3x3<float,int>( const ConstArray3d<float> &filter, const ConstArray3d<int> &a, const Int3 &p );


void PrintFilter( const Array3d<float> &r)
{
  for(int k=0; k<3; ++k)
  {
    printf("k=%i\n",k);
    for(int j=0; j<3; ++j)
    {
      printf("%05f %05f %05f\n",r(0,j,k),r(1,j,k),r(2,j,k));
    }
    printf("\n");
  }
}


static void _BuildFilterCopyZAndOrient( Array3d<float> &a, Array3d<float> &b, float factor, int ax )
{
  for(int j=0; j<3; ++j) for(int i=0; i<3; ++i)
  {
    a(i,j,2) = a(i,j,0) = (factor)*a(i,j,1);
    a(i,j,1) = (1.0-2.0*factor)*a(i,j,1);
  }
  b.init(Int3(3));
  if(ax==2)
    FOR_REG3(int,i,0,2,j,0,2,k,0,2) { b(k,j,i) = a(i,j,k); }
  else if(ax==1)
    FOR_REG3(int,i,0,2,j,0,2,k,0,2) { b(j,i,k) = a(i,j,k); }
  else
    FOR_REG3(int,i,0,2,j,0,2,k,0,2) { b(i,j,k) = a(i,j,k); }
}


void BuildDiff2FilterQuad3d( Array3d<float> &b, int ax )
{
  // second derivative in ax direction times the grid step
  Array3d<float> a(Int3(3));
  a(2,0,1) = a(0,2,1) = a(2,2,1) = a(0,0,1) = (1.0f/12.0f);
  a(1,1,1) = -(10.0f/6.0f);
  a(0,1,1) = a(2,1,1) = (10.0f/12.0f);
  a(1,0,1) = a(1,2,1) = -(1.0/6.0f);
  _BuildFilterCopyZAndOrient(a,b,(1.0f/12.0f),ax);
}


void BuildDiff1FilterQuad3d( Array3d<float> &b, int ax )
{
  // first derivative in ax direction times the grid step
  Array3d<float> a(Int3(3));
  a(1,1,1) = 0;
  a(0,1,1) = -(4.0f/12.0f);
  a(2,1,1) =  (4.0f/12.0f);
  a(0,2,1) = a(0,0,1) = -(1.0f/12.0f);
  a(2,2,1) = a(2,0,1) =  (1.0f/12.0f);
  _BuildFilterCopyZAndOrient(a,b,(1.0f/6.0f),ax);
}


void BuildLaplaceFilterQuad3d( Array3d<float> &r )
{
  // laplace operator,  times the grid step squared
  r.init(Int3(3));
  Array3d<float> dd;
  for(int i=0; i<3; ++i)
  {
    BuildDiff2FilterQuad3d(dd,i);
    //printf("d%c filter:\n",("xyz"[i]));
    //PrintFilter(dd);
    FOR_REG3(int,i,0,2,j,0,2,k,0,2) { r(i,j,k) += dd(i,j,k); }
  }
  //printf("laplace filter:\n");
  //PrintFilter(r);
}


//-----------------------------------
//-----------------------------------


template<class T, int dim>
struct AddSmoothDeltaInternal
{
  static void eval(Array3d<T> &arr, const BBox3 &bbox, const Int3 &ip, const Float3 &fp, T value)
  {
    int a = std::max(bbox.min[dim], ip[dim] + (fp[dim]<0.5f ? -2 : -1));
    int b = std::min(bbox.max[dim], ip[dim] + (fp[dim]>0.5f ?  2 :  1));
    Int3 my_ip(ip);
    for (int i=a; i<=b; ++i) // iterate over neighbor sites of ip along the current dimension
    {
      my_ip[dim] = i;  // set the current point
      float r = i-ip[dim]+0.5f-fp[dim]; // and the distance from the delta function center
      T my_value = my::smooth_delta_cos<T,T>(r, 1.) * value;  // delta is product of 1d delta functions
      AddSmoothDeltaInternal<T, dim-1>::eval(arr, bbox, my_ip, fp, my_value); // one lower dimension
    }
  }
};


template<class T>
struct AddSmoothDeltaInternal<T, -1>
{
  static void eval(Array3d<T> &arr, const BBox3 &bbox, const Int3 &ip, const Float3 &fp, T value)
  {
    arr(ip) += value;
  }
};


template<class T>
void AddSmoothDelta(Array3d<T> &arr, const BBox3 &bbox,  const LatticeDataQuad3d &ld, int dim, const Float3 &pos, T weight)
{
  Int3 ip; Float3 q;
  boost::tie(ip, q) = ld.WorldToLatticeCell(pos);
  // delta func has support of (-2, 2)
  // have to consider [ip-2, ip+2]
  if (dim == 3)
    AddSmoothDeltaInternal<T,2>::eval(arr, bbox, ip, q, weight);
  else if (dim == 2) {
    ip[2] = 0;
    AddSmoothDeltaInternal<T,1>::eval(arr, bbox, ip, q, weight);
  } else if (dim == 1) {
    ip[2] = ip[1] = 0;
    AddSmoothDeltaInternal<T,0>::eval(arr, bbox, ip, q, weight);
  }
}


#define INSTANTIATE_ADDSMOOTHDELTA(T)\
  template void AddSmoothDelta<T>(Array3d<T> &arr, const BBox3 &bbox, const LatticeDataQuad3d &ld, int dim, const Float3 &pos, T weight);


INSTANTIATE_ADDSMOOTHDELTA(float)
INSTANTIATE_ADDSMOOTHDELTA(double)


//-----------------------------------
//-----------------------------------




template<class T>
void Array3dOps<T>::init(Array3d<T>&arr, bool clear) const
{
  //const BBox3 bbext = ExtendForDim(bbox, dim, border);
  arr.initFromBox(bbext, Cons::DONT);
  arr.setBox(bbox);

  if (!clear) return;

  #pragma omp parallel
  {
    BOOST_FOREACH(BBox3 b, mtboxes->getCurrentThreadRange())
    {
      for (int axis=0; axis<dim; ++axis) // extend teh box to include the boundary
      {
        if (b.min[axis]==bbox.min[axis]) b.min[axis] -= border;
        if (b.max[axis]==bbox.max[axis]) b.max[axis] += border;
      }
      arr[b].fill(0.);
    }
  }
}


template<class T>
void Array3dOps<T>::addScaled(double fa, Array3d<T> &u, double fb, const Array3d<T> &v) const
{
  #pragma omp parallel
  {
    BOOST_FOREACH(const BBox3 &b, mtboxes->getCurrentThreadRange())
    {
      u[b].addScaled(fb, v[b], fa);
    }
  }
}


template<class T>
void Array3dOps<T>::initFrom(Array3d<T> &arr, const Array3d<T> &other, ConsMode mode) const
{
  init(arr, false);
  #pragma omp parallel
  {
    BOOST_FOREACH(const BBox3 &b, mtboxes->getCurrentThreadRange())
    {
      if (mode == ConsMode::AS_COPY)
        arr[b].fill(other[b]);
      else
        arr[b].fill(0.);
    }
  }
}


#define INSTANTIATE_ARRAY3DOPS(T) \
  template class Array3dOps<T>;

INSTANTIATE_ARRAY3DOPS(float)
INSTANTIATE_ARRAY3DOPS(double)
INSTANTIATE_ARRAY3DOPS(char)
INSTANTIATE_ARRAY3DOPS(int)
