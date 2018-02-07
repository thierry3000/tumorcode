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

#include "continuum-grid.h"

LatticeIndexRanges LatticeIndexRanges::FromCellRange(const BBox3 &b, int ndim)
{
  LatticeIndexRanges ir;
  ir.cells = b;
  ir.verts = b;
  for (int i=0; i<ndim; ++i)
  {
    ir.verts.max[i] += 1;
    ir.faces[i] = BBox3(b).Move(1, i, 1);
  }
  return ir;
}



ContinuumGrid::ContinuumGrid(const LatticeDataQuad3d &ld_, int dim_) :
  ld(ld_), dim(dim_)
{
  if (dim < 0)
  {
    dim = ld.Size()[2]<=1 ? (ld.Size()[1]<=1 ? (ld.Size()[0]<=1 ? 0 : 1) : 2) : 3;
  }
  ir = LatticeIndexRanges::FromCellRange(ld.Box());
}

void ContinuumGrid::init(const LatticeDataQuad3d& ld_, int dim_)
{
  Destruct(*this);
  new (this) ContinuumGrid(ld_, dim_);
}

BBox3 CellRangeToFaceRange(const BBox3 &b, int ndim, int axis)
{
  myAssert(axis < ndim);
  return BBox3(b).Move(1, axis, 1);
}

BBox3 CellRangeToVertexRange(const BBox3 &b, int ndim)
{
  BBox3 c(b);
  for (int i=0; i<ndim; ++i) c.max[i]++;
  return c;
}

Bool3 CellCenteringForDim(int ndim)
{
  Bool3 b(false);
  for (int i=0; i<ndim; ++i) b[i] = true;
  return b;
}



LatticeDataQuad3d CellLdToFaceLd(const LatticeDataQuad3d &ld, int ndim, int axis)
{
  myAssert(ld.GetCellCentering() == CellCenteringForDim(ndim));
  BBox3 b = CellRangeToFaceRange(ld.Box(), ndim, axis);
  LatticeDataQuad3d ld_face(b, ld.Scale());
  Bool3 is_cell_centered(true); is_cell_centered[axis]=false;
  ld_face.SetCellCentering(is_cell_centered);
  ld_face.SetOriginPosition(ld.GetOriginPosition());
  return ld_face;
}

LatticeDataQuad3d CellLdToVertexLd(const LatticeDataQuad3d &ld, int ndim)
{
  myAssert(ld.GetCellCentering() == CellCenteringForDim(ndim));
  BBox3 b = CellRangeToVertexRange(ld.Box(), ndim);
  LatticeDataQuad3d ld_vert(b, ld.Scale());
  ld_vert.SetCellCentering(Bool3(true));
  ld_vert.SetOriginPosition(ld.GetOriginPosition());
  return ld_vert;
}

//-----------------------------------
//-----------------------------------

void FaceVarArrays::init(const BBox3& bb, int dim)
{
  LatticeIndexRanges ir = LatticeIndexRanges::FromCellRange(bb, dim);
  init(ir, dim);
}

void FaceVarArrays::init(const LatticeIndexRanges& ir, int dim)
{
  FaceVarArrays::Super &arr = *this;
  for (int i=0; i<dim; ++i)
  {
    arr[i].initFromBox(ir.faces[i]);
  }
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
  //myAssert(ld.GetCellCentering() ==Bool3(true));
  Int3 ip; Float3 q;
  boost::tie(ip, q) = ld.WorldToFractionalCoordinate(pos);
  // delta func has support of (-2, 2)
  // have to consider [ip-2, ip+2]
  if (dim == 3)
    AddSmoothDeltaInternal<T,2>::eval(arr, bbox, ip, q, weight);
  else if (dim == 2)
    AddSmoothDeltaInternal<T,1>::eval(arr, bbox, ip, q, weight);
  else if (dim == 1)
    AddSmoothDeltaInternal<T,0>::eval(arr, bbox, ip, q, weight);
}


#define INSTANTIATE_ADDSMOOTHDELTA(T)\
  template void AddSmoothDelta<T>(Array3d<T> &arr, const BBox3 &bbox, const LatticeDataQuad3d &ld, int dim, const Float3 &pos, T weight);


INSTANTIATE_ADDSMOOTHDELTA(float)
INSTANTIATE_ADDSMOOTHDELTA(double)



//-----------------------------------
//-----------------------------------
/** @brief does the index splitting
 * bb: total box elements
 */
void FillBoxGridArray(const BBox3 &bb, const Int3 &cnt, const Int3 &bs, DynArray<BBox3> &res)
{
  Int3 ijk;
  for (ijk[2]=0; ijk[2]<cnt[2]; ++ijk[2])
    for (ijk[1]=0; ijk[1]<cnt[1]; ++ijk[1])
      for (ijk[0]=0; ijk[0]<cnt[0]; ++ijk[0])
      {
        BBox3 r;
        for (int a=0; a<3; ++a)
        {
          r.min[a] = bb.min[a] + ijk[a]*bs[a];
          if (ijk[a] == cnt[a]-1)
            r.max[a] = bb.max[a];
          else
            r.max[a] = r.min[a] + bs[a] - 1;
        }
        res.push_back(r);
      }
}


/*
 * makes a list of blocks which form a grid and their union is bb
 * bs -> hint for block size
 * e.g.: bb={42,39,32}   bs={8,8,8}
 */
DynArray<BBox3> MakeMtBoxGrid(const BBox3 &bb, const Int3 &bs)
{
  Int3 cnt = Size(bb).cwiseQuotient(bs).cwiseMax(Int3(1)); // e.g. cnt={5,5,4}
  DynArray<BBox3> res; res.reserve(cnt.prod());
  FillBoxGridArray(bb, cnt, bs, res);
  return res;
}

DynArray<BBox3> MakeMtBoxGrid(const BBox3 &bb)
{
  return MakeMtBoxGrid(bb, Int3(IS_DEBUG ? 8 : 64));//if debug Int3(8) = {8,8,8}
}


DynArray<BBox3> MakeMtBoxGridLarge(const BBox3 &mainbb, int max_size)
{
  //HACK2018
  int num_threads = 1;
  //int num_threads = my::GetNumThreads();
  int bs = std::min(max_size, std::max(16, maxCoeff(Size(mainbb)) / (2 * num_threads)));
  Int3 cnt = (Size(mainbb)/bs).cwiseMax(Int3(1));
  int total_cnt = cnt.prod();

  DynArray<BBox3> boxes;
  boxes.reserve(total_cnt);
  FillBoxGridArray(mainbb, cnt, Int3(bs), boxes);
  return boxes;
}


DomainDecomposition::DomainDecomposition(const DynArray< BBox3 >& boxes_) : largestBox(BBox3(), -1) //, boost::optional< DynArray< double >& > weights)
{
  //if (weights.is_initialized()) throw std::runtime_error("DomainDecomposition weights not implemented");
  //HACK2018
  //int n = my::OmpGetMaxThreadCount();
  int n = 1;
  for (int i=0; i<boxes_.size(); ++i)
  {
    int thread_num = i % n;
    insert(thread_num, boxes_[i]);
  }
}

DomainDecomposition::DomainDecomposition() : largestBox(BBox3(), -1)
{
  // nothing to do
}


void DomainDecomposition::init(const DynArray< BBox3 >& boxes_)
{
  Destruct(*this);
  new (this) DomainDecomposition(boxes_);
}


DomainDecomposition::range_type DomainDecomposition::getCurrentThreadRange() const
{
  //HACK2018
  //const DynArray<ThreadBox> &ar = by_thread[my::OmpGetCurrentThread()];
  const DynArray<ThreadBox> &ar = by_thread[1];
  return range_type(ar.begin(), ar.end());
}

void DomainDecomposition::insert(int thread_num, const BBox3& bbox_)
{
  ThreadBox bbox(bbox_, boxes.size());
  boxes.push_back(bbox);
  by_thread[thread_num].push_back(bbox);
  if (Volume(bbox) > Volume(largestBox))
    largestBox = bbox;
}

const DomainDecomposition::ThreadBox DomainDecomposition::getBoxWithLargestVolume() const
{
  return largestBox;
}


std::ostream& DomainDecomposition::print(std::ostream &os) const
{
  //HACK2018
  //os << format("DomainDecomposition: %i threads\n") % my::OmpGetMaxThreadCount();
  os << format("DomainDecomposition: %i threads\n") % 1;
  for (int threadidx=0; threadidx<1; ++threadidx)
  //for (int threadidx=0; threadidx<my::OmpGetMaxThreadCount(); ++threadidx)
  {
    os << format("  thread %i\n  ") % threadidx;
    for (int i=0; i<by_thread[threadidx].size(); ++i)
    {
      const ThreadBox &bb = by_thread[threadidx][i];
      os << format(" #%i %s (gid %i),") % i % ((BBox3)(bb)) % bb.global_index;
    }
    os << "\n";
  }
  return os;
}


//-----------------------------------
//-----------------------------------

template<class T>
void AddScaled(double f0, Array3d<T> arr0, double f1, ConstArray3d<T> arr1)
{
  const BBox3 bbox = arr0.getBox();
  assert(arr1.getBox() == bbox); // both arrays must have the same box
  assert(!arr0.hasSameMem(arr1));

  DynArray<BBox3> mpboxes = MakeMtBoxGridLarge(bbox);

  #pragma omp parallel for schedule(dynamic, 1)
  for (int i=0; i<mpboxes.size(); ++i)
  {
    const BBox3 b = mpboxes[i];
    arr0[b].addScaled(f1, arr1[b], f0);
  }
}

#define INSTANTIATE_ADDSCALED(T) \
  template void AddScaled<T>(double, Array3d<T>, double, ConstArray3d<T>);

INSTANTIATE_ADDSCALED(double)
INSTANTIATE_ADDSCALED(float)





