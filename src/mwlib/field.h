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
#ifndef _FIELD_H
#define _FIELD_H

#include <boost/iterator/iterator_facade.hpp>
#include <boost/range/iterator_range.hpp>
//#include <boost/range/algorithm_ext/for_each.hpp>
//#include <boost/range/algorithm/generate.hpp>

#include "dynamicarray.h"
#include "helpers-vec.h"
#include "operators.h"
#include "helpers-sys.h"
#include "helpers-mem.h"
#include "lattice-data.h"
#include "refcounted_array.h"
#include "drawing.h"


template<class U> class ConstArray3d;
template<class U> class Array3d;

#if 0
#define FOR_AREG3(i,mo,x0,x1,y0,y1,z0,z1)\
  for(int _moz=(mo)[2],\
          _moy=(mo)[1],\
          _mox=(mo)[0],\
           _z=(z0)*_moz,\
          _z1=(z1)*_moz,\
          _y0=(y0)*_moy,\
          _y1=(y1)*_moy,\
          _x0=(x0)*_mox,\
          _x1=(x1)*_mox;\
          _z<=_z1; _z+=_moz)\
    for(int _j=_z+_y0,_j1=_z+_y1; _j<=_j1; _j+=_moy)\
      for(int i=_j+_x0,_i1=_j+_x1; i<=_i1;  i+=_mox)


//for looping over array regions, mo is an int[3] of memory coefficients
#define FOR_ARRAY3(i, a)\
  FOR_AREG3(i,a.strides(),0,a.size().x()-1,0,a.size().y()-1,0,a.size().z()-1)

#define FOR_ARRAY3_EX(i,t,a,aexpr)\
  for(bool __t=true; __t;)\
    for(t a = aexpr; __t; __t=false)\
      FOR_ARRAY3(i, a)
#endif


namespace _ArrayNdInternal
{
  template<int axis>
  struct applicator_dim
  {
    template<class F, class P0>
    inline static void apply(P0 p0, const int *strides0, const int *size, F &f)
    {
      const P0 end = p0 + size[axis] * strides0[axis];
      for(; p0 < end; p0 += strides0[axis])
        applicator_dim<axis-1>::apply(p0, strides0, size, f);
    }
    template<class F, class P0, class P1>
    inline static void apply(P0 p0, const int *strides0, P1 p1, const int *strides1, const int *size, F &f)
    {
      const P0 end = p0 + size[axis] * strides0[axis];
      for(; p0 < end; p0 += strides0[axis], p1 += strides1[axis])
        applicator_dim<axis-1>::apply(p0, strides0, p1, strides1, size, f);
    }
  };

  template<>
  struct applicator_dim<-1>
  {
    template<class F, class P0>
    inline static void apply(P0 p0, const int *strides0, const int *size, F &f)
    {
      f.inplace(*p0);
    }
    template<class F, class P0, class P1>
    inline static void apply(P0 p0, const int *strides0, P1 p1, const int *strides1, const int *size, F &f)
    {
      f.inplace(*p0, *p1);
    }
  };

  template<class F, class A>
  F apply(A &a, F f)
  {
    const Int3 size = a.size();
    typedef typename A::pointer_type P0;
    const Int3 strides0 = a.strides();
    P0 basep0 = a.getPtr() + a.offset(a.getBox().min);
#ifdef DEBUG
    if (allG(size, 0))
    {
      a.debugCheck(basep0);
      a.debugCheck(basep0 + a.offset(size - Int3(1)));
    }
#endif
    applicator_dim<2>::apply<F,P0>(basep0, strides0.data(), size.data(), f);
    return f;
  }

  template<class F, class A0, class A1>
  F apply(A0 &a0, A1 &a1, F f)
  {
    const Int3 size = a0.size();
    myAssert(size == a1.size());
    typedef typename A0::pointer_type P0;
    const Int3 strides0 = a0.strides();
    P0 basep0 = a0.getPtr() + a0.offset(a0.getBox().min);
    typedef typename A1::pointer_type P1;
    const Int3 strides1 = a1.strides();
    P1 basep1 = a1.getPtr() + a1.offset(a1.getBox().min);
#ifdef DEBUG
    if (allG(size,0))
    {
      a0.debugCheck(basep0);
      a0.debugCheck(basep0 + a0.offset(size - Int3(1)));
      a1.debugCheck(basep1);
      a1.debugCheck(basep1 + a1.offset(size - Int3(1)));
    }
#endif
    applicator_dim<2>::apply<F,P0,P1>(basep0, strides0.data(), basep1, strides1.data(), size.data(), f);
    return f;
  }
  
  template<class F, class A0, class A1, class A2>
  F apply(A0 &a0, A1 &a1, A2 &a2, F f)
  {
    const Int3 size = a0.size();
    myAssert(size == a1.size());
    myAssert(size == a2.size());
    typedef typename A0::pointer_type P0;
    const Int3 strides0 = a0.strides();
    P0 basep0 = a0.getPtr() + a0.offset(a0.getBox().min);
    typedef typename A1::pointer_type P1;
    const Int3 strides1 = a1.strides();
    P1 basep1 = a1.getPtr() + a1.offset(a1.getBox().min);
    typedef typename A2::pointer_type P2;
    const Int3 strides2 = a2.strides();
    P2 basep2 = a2.getPtr() + a2.offset(a2.getBox().min);
#ifdef DEBUG
    if (allG(size,0))
    {
      a0.debugCheck(basep0);
      a0.debugCheck(basep0 + a0.offset(size - Int3(1)));
      a1.debugCheck(basep1);
      a1.debugCheck(basep1 + a1.offset(size - Int3(1)));
      a2.debugCheck(basep2);
      a2.debugCheck(basep2 + a2.offset(size - Int3(1)));
    }
#endif
    applicator_dim<2>::apply<F,P0>(basep0, strides0.data(), basep1, strides1.data(), basep2, strides2.data(), size.data(), f);
    return f;
  }
}



namespace _ArrayNdInternal
{

template<class T>
class Array3dIterator_ : 
  public boost::iterator_facade<Array3dIterator_<T>, T, boost::forward_traversal_tag>
{  
  struct enabler {};
public:
  Array3dIterator_() : ptr(NULL), base(NULL) {}
  
  Array3dIterator_(const Int3 &p, T* ptr_, const Int3 &m, const Int3 &l
#ifdef ARRAYND_DEBUG
    , const T* blockBegin_, const T* blockEnd_
#endif
    ) : m(m), l(l), p(p), base(ptr), ptr(ptr_)
  {
#ifdef ARRAYND_DEBUG
    blockBegin = blockBegin_;
    blockEnd = blockEnd_;
#endif
    ptr += p.dot(m);
  }

  template<class OtherValue>
  Array3dIterator_(Array3dIterator_<OtherValue> const& other) :
    m(other.m), l(other.l), p(other.p), base(other.base), ptr(other.ptr)
#ifdef ARRAYND_DEBUG
      ,blockBegin(other.blockBegin), blockEnd(other.blockEnd)
#endif
    {}
#ifdef DEBUG
    const T* blockBegin, *blockEnd;
#endif

private:
  friend class boost::iterator_core_access;
  template<class> friend class Array3dIterator_;

  void increment()
  {
    ++p.x();
    ptr += m.x();
    if(p.x()>=l.x()) 
    {
      p.x() = 0;
      ++p.y();
      if(p.y()>=l.y())
      {
        p.y() = 0;
        ++p.z();
        if(p.z()>=l.z())
        {
          ptr = base + l.dot(m);
          return;
        }
      }
      ptr = base + p.y()*m.y() + p.z()*m.z();
    }
  }
  
  T& dereference() const
  {
    ARRAYND_ASSERT(ptr>=blockBegin && ptr<blockEnd);
    return *ptr;
  }

  template<class OtherValue>
  bool equal(Array3dIterator_<OtherValue> const& other) const
  {
    return other.ptr == ptr;
  }
  
  Int3 m,l,p;
  T* ptr, *base;
};

template<class T, class U>
struct CopyCons
{
  void operator()(T &a, const U &b) const { new (&a) T(b); }
};

template<class T, class U>
struct StatsOp
{
  my::Averaged<T> &avg;
  StatsOp(my::Averaged<T> &avg) : avg(avg) {}
  void inplace(const U &value) const { avg.Add(value); }
};

template<class T>
struct MaxAbsOp
{
  T &res;
  MaxAbsOp(T &res_) : res(res_) {}
  void inplace(const T &value) const { res = std::max(std::abs(value), res); }
};

}


template<class T>
class ConstArray3d
{
  template<class U> friend class ConstArray3d;
  template<class U> friend class Array3d;

public:
  typedef ConstArray3d<T> This;
  typedef LatticeIndexing<3, const T*, int> Lattice;
  typedef T const value_type;
  typedef T const* pointer_type;
  typedef T const* const_pointer_type;
  //typedef int index_t;
  typedef int offset_t;
  enum { MAX_DIM = 3 };

protected:
  ConstArray3d(const Int3 &_l, const Int3 &_m, offset_t moffset, const ArrayNdRef &newMem) : memRef(newMem)
  {
    lattice.Init(_l, getBasePtr() + moffset, _m);
  }

  ConstArray3d(const Int3 &l_)
  {
    if (allEq(l_,0)) return;
    Int3 l = l_.cwiseMax(1);
    memRef.alloc(l.prod(), sizeof(T),this->destructFunc);
    lattice.Init(l, getBasePtr());
  }

  T* getNonConstBasePtr()
  {
    return static_cast<T*>(memRef.storage());
  }

  T* getNonConstPtr()
  {
    return const_cast<T*>(lattice.SiteAtOrigin());
  }

public:
  ConstArray3d()
  {
  }

  ConstArray3d(const ConstArray3d &a) : memRef(a.memRef),lattice(a.lattice)
  {
  }

  struct MoveHelperConst
  {
    MoveHelperConst(const ConstArray3d<T> &x_) : x(x_) {}
    MoveHelperConst(const MoveHelperConst &mh) : x(mh.x) { mh.clear(); }
    void clear() const { const_cast<ConstArray3d<T>&>(x).clear(); }
    ConstArray3d<T> x;
  };

  MoveHelperConst release()
  {
    MoveHelperConst mh(*this);
    this->clear();
    return mh;
  }
  
  ConstArray3d(const MoveHelperConst &mh) : memRef(mh.x.memRef), lattice(mh.x.lattice)
  {
    mh.clear();
  }

  const ConstArray3d<T>& operator=(const MoveHelperConst &mh)
  {
    *this = mh.x;
    mh.clear();
    return *this;
  }

  static Int3 defaultStrides() { return Int3::Constant(-1); }

  ConstArray3d(const Int3 &_l, const Int3 &_m, const T* mem, std::size_t buffer_size, bool owned)
  {
    init(_l, _m, mem, buffer_size, owned);
  }

  ConstArray3d(const BBox3 &bbox_, const Int3 &strides_, const T* addr_of_bbox_min, std::size_t buffer_size, bool owned)
  {
    init(bbox_, strides_, addr_of_bbox_min, buffer_size, owned);
  }

  void init(const Int3 &_l, const Int3 &_m, const T* mem, std::size_t buffer_size, bool owned)
  {
    myAssert(_l.prod() <= buffer_size);
    memRef.take((void*)(mem), buffer_size, sizeof(T),owned ? destructFunc : NULL);
    lattice.Init(_l, getBasePtr(), _m);
  }

  void init(const BBox3 &bbox_, const Int3 &strides_, const T* addr_of_bbox_min, std::size_t buffer_size, bool owned)
  {
    myAssert(Volume(bbox_) <= buffer_size);
    memRef.take((void*)(addr_of_bbox_min), buffer_size, sizeof(T), owned ? destructFunc : NULL);
    lattice.Init(bbox_, addr_of_bbox_min, strides_, Lattice::AT_BBOX_MIN);
  }

  ~ConstArray3d()
  {
  }

  const ConstArray3d<T>& operator=(const ConstArray3d<T>&a)
  {
    if(&a == this) return *this;
    memRef = a.memRef;
    lattice = a.lattice;
    return *this;
  }

  void move(const Int3 &ioffset)
  {
    this->debugCheck(offset(getBox().min));
    this->debugCheck(offset(getBox().max));

    lattice.MoveSiteRange(ioffset);
    lattice.SetBox(Move(lattice.Box(), ioffset));

    this->debugCheck(offset(getBox().min));
    this->debugCheck(offset(getBox().max));
  }

  // Init with a new index range (like initFromBox) but keep the memory allocation.
  // YOU MUST make sure that the memory allocation is large enough!
  void reshape(const BBox3 &bbox)
  {
    myAssert(memRef.allocation_size() <= Volume(bbox));
    lattice.Init(bbox, static_cast<T*>(memRef.storage()), defaultStrides(), Lattice::AT_BBOX_MIN);
  }
  
  void setSlice(const BBox3 &_bb)
  {
    this->debugCheck(offset(_bb.min));
    this->debugCheck(offset(_bb.max));

#ifdef DEBUG
    BBox3 bbold = lattice.Box();
#endif
    
    // _bb.min becomes the new origin, index range starts at 0
    const BBox3 bo(Int3(0), _bb.max - _bb.min);
    lattice.SetBox(bo);
    lattice.MoveSiteRange(-_bb.min);
    
#ifdef DEBUG
    this->debugCheck(offset(bbold.min-_bb.min));
    this->debugCheck(offset(bbold.max-_bb.min));
    this->debugCheck(offset(Int3(0)));
    this->debugCheck(offset(_bb.max - _bb.min));
#endif
  }

  void setBox(const BBox3 &_bb)
  {
    // set the region box, index range does not change
    this->debugCheck(offset(_bb.min));
    this->debugCheck(offset(_bb.max));

    lattice.SetBox(_bb);

    this->debugCheck(offset(_bb.min));
    this->debugCheck(offset(_bb.max));
  }

  ConstArray3d<T> constSlice(const BBox3 &_bbox) const
  {
    ConstArray3d<T> a(*this);
    a.setSlice(_bbox);
    return a;
  }

  ConstArray3d<T> operator[](const BBox3 &_bbox) const
  {
    return constSlice(_bbox);
  }

  void setFlat()
  {
    myAssert(isContiguous());

    Lattice old(lattice);
    lattice.Init(Int3(Size(lattice.Box()),1,1));
    lattice.SetBase(Int3(0), old.LatticeToSite(old.Box().min));

    this->debugCheck(offset(Int3(0)));
    this->debugCheck(offset(Int3(size()[0]-1, 0, 0)));
  }

  Int3 getCoords(offset_t id) const
  {
    myAssert(isContiguous());
    return lattice.SiteToLattice(lattice.Offset()+id);
  }


  void RemoveSingularDims()
  {
    Int3 dims = Size(getBox());
    int j = 0;
    // for example in some iteration
    // 00x00xx
    // j i
    // becomes
    // x0000xx
    //  j i
    for (int i=0; i<3; ++i)
    {
      if (dims[i] <= 1) continue;
      if (i != j) swapAxes(i, j);
      j++;
    }
  }

  // compute offset into the mem array from coordinates
  int offset(int x, int y, int z) const
  {
    const Int3 &m = lattice.Strides();
    return x*m[0]+y*m[1]+z*m[2];
  }

  template<class Index>
  int offset(const Index &p) const
  {
    return offset(p[0],p[1],p[2]);
  }

  T operator()(int x, int y, int z) const
  {
    return getPtr()[this->debugCheck(offset(x,y,z))];
  }

  template<class Index>
  T operator()(const Index &p) const
  {
    return getPtr()[this->debugCheck(offset<Index>(p))];
  }

  T dereference(int i) const
  {
    return getPtr()[this->debugCheck(i)];
  }

  const T* getBasePtr() const
  {
    return const_cast<This*>(this)->getNonConstBasePtr(); 
  }

  const T* getPtr() const
  {
    return const_cast<This*>(this)->getNonConstPtr();
  }

  bool isContiguous() const
  {
    return lattice.Strides() == Lattice::CalcStrides(lattice.Size());
  }

  bool empty() const {
    return size() == Vec<int,3>(0);
  }

  bool hasSameLattice(const ConstArray3d<T> &a) const
  {
    return lattice == a.lattice;
  }

  bool hasSameContent(const ConstArray3d<T> &a) const
  {
    if (!(lattice == a.lattice)) return false;

    struct Compare
    {
      bool ok;
      Compare() : ok(true) {}
      void inplace(const T &a, const T &b) { ok &= (a==b); }
    };

    Compare cmp = _ArrayNdInternal::apply(a, *this, Compare());
    return cmp.ok;
  }

  template<class U>
  bool hasSameMem(const ConstArray3d<U> &a) const
  {
    if (memRef.header()==NULL && a.memRef.header()==NULL) return false;
    return memRef == a.memRef;
  }

  offset_t count() const
  {
    return lattice.NumSites();
  }

  uint memSize() const
  {
    return sizeof(*this)+(memRef.header() ? memRef.header()->cnt*sizeof(T) : 0);
  }

  int refCount() const
  {
    return memRef.header() ? memRef.header()->ref_cnt : 0;
  }

  void clear()
  {
    memRef = ArrayNdRef();
    lattice = Lattice();
  }

  Int3 strides() const
  {
    return lattice.Strides();
  }

  const Int3& size() const
  {
    return lattice.Size();
  }

  Int2 size2() const
  {
    const Int3 &s = size();
    return Vec<int,2>(s[0],s[1]);
  }
  
  BBox3 getBox() const
  {
    return lattice.Box();
  }
  
  BBox2 getBox2() const
  {
    const BBox3 &b = lattice.Box();
    return BBox2(b.min[0],b.min[1],b.max[0],b.max[1]);
  }

  void swap(ConstArray3d<T> &a)
  {
    a.memRef.swap(memRef);
    std::swap(lattice, a.lattice);
  }

#ifdef ARRAYND_DEBUG
  const T* blockBegin() const
  {
    return (T*)memRef.storage();
  }
  const T* blockEnd() const
  {
    if (memRef.storage())
      return (T*)((char*)memRef.storage() + memRef.header()->cnt*memRef.header()->itemsize);
    else
      return NULL;
  }
  const T* debugCheck(const T* p) const
  {
    AssertAlways(p>=blockBegin() && p<blockEnd());
    return p;
  }
  offset_t debugCheck(const offset_t i) const
  {
    AssertAlways(getPtr()+i>=blockBegin() && getPtr()+i<blockEnd());
    return i;
  }
#else
  offset_t debugCheck(const offset_t i) const { return i; }
  const T* debugCheck(const T* p) const { return p; }
#endif

  void swapAxes(int a, int b)
  {
    lattice.SwapAxes(a, b);
  }

  typedef _ArrayNdInternal::Array3dIterator_<T const> const_iterator;
  const const_iterator getIter(const Int3 &pos) const
  {
    return const_iterator( pos
                          ,lattice.LatticeToSite(lattice.Box().min)
                          ,strides()
                          ,size()
#ifdef ARRAYND_DEBUG
                          ,blockBegin()
                          ,blockEnd()
#endif
                          );
  }
  const const_iterator begin() const {  return getIter(Int3(0));  }
  const const_iterator end() const  {  return getIter(size());  }
  typedef boost::iterator_range<const_iterator> ConstRange;
  const ConstRange range() const { return ConstRange(begin(),end()); }

  my::Averaged<double> valueStatistics() const
  {
    my::Averaged<double> avg;
    _ArrayNdInternal::apply(*this, _ArrayNdInternal::StatsOp<double,T>(avg));
    return avg;
  }

  T maxAbs() const
  {
    T res = 0.;
    _ArrayNdInternal::apply(*this, _ArrayNdInternal::MaxAbsOp<T>(res));
    return res;
  }

  std::size_t estimateMemoryUsage() const
  {
    std::size_t sz = sizeof(Array3d<T>);
    if (!memRef.header()) return sz;
    std::size_t cnt = memRef.header()->cnt;
    if (cnt == 0) return sz;
    return sz + cnt * estimateMemoryUsage(getBasePtr()[0]);
  }

protected:
  static void destructFunc(void *mem, size_t cnt)
  {
    DestructVec<T>((T*)mem,cnt);
  }

  ArrayNdRef memRef;
  Lattice lattice;
};


template <class T>
inline std::size_t estimateMemoryUsage(const ConstArray3d<T> &a)
{
  return a.estimateMemoryUsage();
}





template<class T>
class Array3d : public ConstArray3d<T>
{
  template<class U> friend class ConstArray3d;
  template<class U> friend class Array3d;
public:
  typedef ConstArray3d<T> Super;
  typedef T value_type;
  typedef T * pointer_type;
  typedef typename Super::offset_t offset_t;
  //typedef typename Super::index_t index_t;    
  using Super::count;
protected:
  using Super::memRef;
  using Super::lattice;
  //using Super::mptr;
  //using Super::l;
  //using Super::m; 
#ifdef ARRAYND_DEBUG  
  using Super::blockBegin;
  using Super::blockEnd;
#endif  
  
  Array3d(const Int3 &l,const Int3 &m,offset_t moffset,ArrayNdRef &newMem) : Super(l,m,moffset,newMem)
  {
  } 
//   template<class RangeType>
//   void constructFromRange(const RangeType &range_)
//   {
//     boost::for_each(std::make_pair(getPtr(),getPtr()+count()), range_, _ArrayNdInternal::CopyCons<T,typename RangeType::value_type>());
//   }  
public:

  //basic constructors
  Array3d() : Super()
  {
  }

  Array3d(const Array3d &a) : Super(a)
  {
  }



  Array3d(const Int3 &_l, const Int3 &_m, const T* mem, std::size_t buffer_size, bool owned) : Super(_l, _m, mem, buffer_size, owned)
  {
#ifdef DEBUG
    printf("Am I using this?");
#endif
  }

  Array3d(const BBox3 &bbox_, const Int3 &strides_, const T* addr_of_bbox_min, std::size_t buffer_size, bool owned) : Super(bbox_, strides_, addr_of_bbox_min, buffer_size, owned)
  {
  }

  struct MoveHelper : public ConstArray3d<T>::MoveHelperConst
  {
    MoveHelper(const Array3d<T> &x_) : ConstArray3d<T>::MoveHelperConst(x_) {}
    MoveHelper(const typename ConstArray3d<T>::MoveHelperConst &mh) : ConstArray3d<T>::MoveHelperConst(mh) {}
  };

  MoveHelper release()
  {
    MoveHelper mh(*this);
    this->clear();
    return mh;
  }

  Array3d(const MoveHelper &mh) : Super(mh)
  {
  }

  const Array3d& operator=(const MoveHelper &mh)
  {
    Super::operator=(mh);
    return *this;
  }

  template<class U>
  Array3d(const ConstArray3d<U> &a, Cons::COPY_TYPE) : Super(a.size())
  {
    this->fill(a[a.getBox()]);
    this->move(a.getBox().min);
  }

  template<class U>
  Array3d(const ConstArray3d<U> &a, Cons::DEEP_COPY_TYPE) : Super()
  {
    if (!a.memRef.header()) return;
    memRef.alloc(a.memRef.header()->cnt,sizeof(T),ConstArray3d<T>::destructFunc);
    lattice.Init(a.lattice.Box(), Super::getBasePtr() + (a.getPtr() - a.getBasePtr()), a.lattice.Strides());
    ConstructVecCopy<T,U>(getBasePtr(), a.memRef.header()->cnt, a.getBasePtr());
  }

  const Array3d& operator=(const Array3d &a)
  {
    Super::operator=(a);
    return *this;
  }


  Array3d(const Int3 &s, Cons::DONT_TYPE) : Super(s)
  {
  }

  explicit Array3d(const Int3 &s, const T& filler=T()) : Super(s)
  {
    ConstructVec<T>(getPtr(), count(), filler);
  }
  
  void init(const Int3 &s, const T &filler=T())
  {
#if 1
    bool mem_reused = false;
    if (memRef.header() && memRef.header()->ref_cnt == 1)
    {
      int n = s.prod();
      if (memRef.header()->cnt == n)
      {
        lattice.Init(s, Super::getBasePtr());
        fill(filler);
        mem_reused  = true;
      }
    }
    if (!mem_reused)
#endif
    {
      Destruct(*this);
      new (this) Array3d(s, filler);
    }
  }

  void init(const Int3 &s, Cons::DONT_TYPE)
  {
    Destruct(*this);
    new (this) Array3d(s, Cons::DONT);
  }
  
  void initFromBox(const BBox3 &bb, const T &filler=T())
  {
    Destruct(*this);
    new (this) Array3d(bb, filler);
  }

  void initFromBox(const BBox3 &bb, Cons::DONT_TYPE)
  {
    Destruct(*this);
    new (this) Array3d(bb, Cons::DONT);
  }
  
  explicit Array3d(const BBox3 &bb, const T &filler=T()) : Super(Size(bb))
  {
    ConstructVec<T>(getPtr(), count(), filler);
    if (this->memRef.allocation_size()>0)
      this->move(bb.min);
  }

  Array3d(const BBox3 &bb, Cons::DONT_TYPE) : Super(Size(bb))
  {
    this->move(bb.min);
  }
  
  // init funcs
//   template<class S>
//   void initFromShape(const S &s, const T& filler=T())
//   {
//     Int3 ss(1);
//     myAssert(s.size() <= 3);
//     for (int i=0; i<std::min(3,s.size()); ++i) ss[i] = s[i];
//     init(ss, filler);
//   }
  
//   template<class RangeType>
//   void initFromRange(const Int3 &s, const RangeType &range_)
//   {
//     Destruct(*this);
//     new (this) Array3d(s, Cons::DONT);
//     boost::for_each(std::make_pair(this->mptr,this->mptr+this->count()), range_, _ArrayNdInternal::CopyCons<T,typename RangeType::value_type>());
//   }
// 
//   template<class RangeType>
//   void initFromBoxRange(const BBox3 &b, const RangeType &range_)
//   {
//     initFromRange(::Size(b), range_);
//     Super::moveIndexRange(b.min);
//   }

  template<class U>
  void initCopy(const ConstArray3d<U> &a)
  {
    Destruct(*this);
    new (this) Array3d(a,Cons::COPY);
  }

  template<class U>
  void initDeepCopy(const ConstArray3d<U> &a)
  {
    Destruct(*this);
    new (this) Array3d(a,Cons::DEEP_COPY);
  }



  // static functions, constructing arrays based on init funcs
//   template<class S>
//   static Array3d<T> fromShape(const S &s, const T &filler=T())
//   {
//     Array3d<T> a;
//     a.initFromShape(s, filler);
//     return a;
//   }

//   static Array3d<T> fromBox(const BBox3 &box, const T &filler=T())
//   {
//     Array3d<T> a;
//     a.initFromBox(box, filler);
//     return a;
//   }

//   template<class RangeType>
//   static Array3d<T> fromRange(const Int3 &s, const RangeType &range_)
//   {
//     Array3d<T> a;
//     a.initFromRange(s, range_);
//     return a;
//   }
// 
//   template<class RangeType>
//   static Array3d<T> fromBoxRange(const BBox3 &b, const RangeType &range_)
//   {
//     Array3d<T> a;
//     a.initFromBoxRange(b, range_);
//     return a;
//   }

  //template<class U, class Map>
  //static const Array3d<T> map(const ConstArray3d<U> &a, const Map &map)
  //{
  //  struct F {
  //    const Map& m;
  //    F(const Map& m) : m(m) {}
  //    void inplace(T &a, const U &b) const { new (&a) T(m(b)); }
  //  };
  //  Array3d<T> res(a.size(), Cons::DONT);
  //  _ArrayNdInternal::apply<F, Array3d<T>, const ConstArray3d<U> >(res, a, F(map));
  //  return res;
  //}
  
  static Array3d<T> getComponent(Array3d< Vec<T,3> > &a, int i)
  {
    offset_t off = 3 * (a.getPtr()-(Vec<T,3>*)a.memRef.storage()) + i;
    Array3d<T> res(a.lattice.Size(), 3 * a.strides(), off, a.memRef);
    res.setBox(a.getBox());
    return res;
  }

  Array3d<T> slice(const BBox3 &_bbox)
  {
    Array3d<T> a(*this);
    a.setSlice(_bbox);
    return a;
  } 

  Array3d<T> operator[](const BBox3 &_bbox)
  {
    Array3d<T> a(*this);
    a.setSlice(_bbox);
    return a;
  }

//   DynArray<T> dynArray() const
//   { 
//     return DynArray<T>(memRef, getPtr(), count());
//   }

  // inherited operators, brought to this namespace
  using Super::operator[];
  //using Super::operator();
  using Super::getPtr;
  using Super::offset;
  using Super::dereference;
  
  // [] access by uint indices
//   T& operator[](int id)
//   {
//     return Super::getNonConstPtr()[debugCheck(offset(id))];
//   }  

  T& dereference(offset_t i)
  {
    return Super::getNonConstPtr()[this->debugCheck(i)];
  }
  
  T& operator()(int x, int y, int z)
  {
    return Super::getNonConstPtr()[this->debugCheck(offset(x, y, z))];
  }

  template<class Index>
  T& operator()(const Index &p)
  {
    return Super::getNonConstPtr()[this->debugCheck(offset(p))];
  }
  
  template<class Index>
  T operator()(const Index &p) const
  {
    return Super::operator()(p);
  }

  T operator()(int x, int y, int z) const
  {
    return Super::operator()(x, y, z);
  }

  T* getPtr()
  {
    return Super::getNonConstPtr();
  }

  T* getBasePtr()
  {
    return Super::getNonConstBasePtr();
  }

  void fill(const T &fill) 
  {
    _ArrayNdInternal::apply(*this,Ops::Bind2nd(Ops::Assign<T>(),fill));
  }

  template<class U>
  void fill( const ConstArray3d<U> &a )
  {
    _ArrayNdInternal::apply(*this,a,Ops::Assign<T>());
  }

  Array3d& operator*=(const T &q)
  {
    _ArrayNdInternal::apply(*this,Ops::Bind2nd(Ops::Mult<T>(),q));
    return *this;
  }

  template<class U>
  Array3d& operator*=(const ConstArray3d<U>& a)
  {
    _ArrayNdInternal::apply(*this,a,Ops::Mult<T>());
    return *this;
  }

  Array3d& operator+=(const T &q)
  {
    _ArrayNdInternal::apply(*this,Ops::Bind2nd(Ops::Add<T>(),q));
    return *this;
  }

  template<class U>
  Array3d& operator+=(const ConstArray3d<U>& a)
  {
    _ArrayNdInternal::apply(*this,a,Ops::Add<T>());
    return *this;
  }

  Array3d& operator-=(const T &q)
  {
    _ArrayNdInternal::apply(*this,Ops::Bind2nd(Ops::Subtract<T>(),q));
    return *this;
  }

  template<class U>
  Array3d& operator-=(const ConstArray3d<U>& a)
  {
    _ArrayNdInternal::apply(*this,a,Ops::Subtract<T>());
    return *this;
  }  

  template<class U>
  Array3d& addScaled(T s, const ConstArray3d<U> &a, T s_self = T(1))
  {
    if (s_self == T(1))
      _ArrayNdInternal::apply(*this, a, Ops::MultAdd<T>(s));
    else if (s_self == T(0) && s == T(1))
      this->fill(a);
    else
      _ArrayNdInternal::apply(*this, a, Ops::MultAddScaleSelf<T>(s, s_self));
    return *this;
  }

  typedef _ArrayNdInternal::Array3dIterator_<T> iterator;
  const iterator getIter(const Int3 &pos) const
  {
          return iterator( pos
                          ,Super::lattice.LatticeToSite(lattice.Box().min)
                          ,Super::strides()
                          ,Super::size()
#ifdef ARRAYND_DEBUG
                          ,Super::blockBegin()
                          ,Super::blockEnd()
#endif
                          );
  }
  const iterator begin() {  return getIter(Int3(0));  }
  const iterator end()  {  return getIter(Super::size());  }
  typedef boost::iterator_range<iterator> Range;
  Range range() { return Range(begin(),end()); }
};


typedef Array3d<float> Array3df;
typedef Array3d<double> Array3dd;



#define DECL_ARRAY_CWISE_MATH_BIN(op, op2)\
template<class T>\
inline Array3d<T> operator op (const ConstArray3d<T> &a, const ConstArray3d<T> &b)\
{\
  Array3d<T> q(a, Cons::COPY);\
  q op2 b;\
  return q;\
}

DECL_ARRAY_CWISE_MATH_BIN(+,+=)
DECL_ARRAY_CWISE_MATH_BIN(-,-=)
DECL_ARRAY_CWISE_MATH_BIN(*,*=)
DECL_ARRAY_CWISE_MATH_BIN(/,/=)


template<class T, class U>
Array3d<T> operator*(const ConstArray3d<T> &a, const U &u)
{
  Array3d<T> tmp(a, Cons::COPY);
  tmp *= u;
  return tmp;
}

template<class T, class U>
Array3d<T> operator*(const U &u, const ConstArray3d<T> &a)
{
  return a*u;
}

template<class T, class U>
Array3d<T> operator+(const ConstArray3d<T> &a, const U &u)
{
  Array3d<T> tmp = Array3d<T>::copy(a);
  tmp += u;
  return tmp;
}

template<class T, class U>
Array3d<T> operator+(const U &u, const ConstArray3d<T> &a)
{
  return a+u;
}


class Image;

template<class T>
Array3d<T> ReadArrayFromImage(const string &fn, int channel, int zsize);


template<class T>
void DrawArrayInt(Image &img, const ConstArray3d<T> &arr, bool normalized, double val_scale, double val_offset, bool bText, const string &title);

struct DrawArrayOpts
{
  DrawArrayOpts() : slice_(0), slice_axis_(-1), outputRange_(false), scale_(1.), offset_(0.), normalize_(true) {}
  DrawArrayOpts& scaled(double scale, double offset) { normalize_=false; scale_=scale; offset_=offset; return *this; } // scaled value := val * scale + offset
  DrawArrayOpts& title(const string &title) { title_=title; return *this; }
  DrawArrayOpts& outputRange() { outputRange_=true; return *this; }
  DrawArrayOpts& zslice(int i) { slice_axis_=2; slice_=i; return *this; }
  DrawArrayOpts& slice(int axis, int coord) { slice_axis_ = axis; slice_ = coord; return *this; }

  int slice_axis_, slice_;
  bool outputRange_, normalize_;
  string title_;
  double scale_, offset_;
};

// new draw functions
template<class T>
inline void DrawArray(Image &img, ConstArray3d<T> arr, const DrawArrayOpts &opts = DrawArrayOpts())
{
  if (opts.slice_axis_>=0)
  {
    arr = arr[BBox3(arr.getBox()).Set(opts.slice_axis_, opts.slice_)];
    arr.RemoveSingularDims();
  }
    //arr = arr[BBox3(arr.getBox()).Set(2, opts.zslice_)];
  DrawArrayInt(img, arr, opts.normalize_, opts.scale_, opts.offset_, opts.outputRange_, opts.title_);
}

template<class T>
inline Image DrawArray(ConstArray3d<T> arr, const DrawArrayOpts &opts = DrawArrayOpts())
{
  Image img;
  DrawArray(img, arr, opts);
  return img;
}


// old draw functions
template<class T>
inline void DrawArrayNormalized(Image &img, const ConstArray3d<T> &arr, bool bText)
{
  DrawArrayInt(img, arr, true, 1., 0., bText, string());
}

template<class T>
inline void DrawArray(Image &img, const ConstArray3d<T> &arr, double val_scale, double val_offset, bool bText)
{
  DrawArrayInt(img, arr, false, val_scale, val_offset, bText, string());
}

template<class T>
inline void WriteArrayAsImage(const string &s, const ConstArray3d<T> &arr, double val_scale, double val_offset, bool bText)
{
  Image img;
  DrawArray<T>(img,arr, val_scale, val_offset, bText);
  img.Write(s);
}

template<class T>
inline void WriteArrayAsImageNormalized(const string &s, const ConstArray3d<T> &arr, bool bText = false)
{
  Image img;
  DrawArrayNormalized<T>(img,arr,bText);
  img.Write(s);
}



template<class T> const T* get_ptr( const ConstArray3d<T> &a ) { return a.getPtr(); }
template<class T> T* get_ptr( Array3d<T> &a ) { return a.getPtr(); }



enum {
  CONVOLVE_NOCLIP = 1
};


template<class A> 
struct get_ptr_type
{
    typedef typename Select<IsConst<A>::result,typename A::const_pointer_type,typename A::pointer_type>::Result Result;
};


template<class T,class S,class F>
static void array_apply_offset3d(Array3d<T> &a, const ConstArray3d<S> &b, F &func, const Int3 &boffset, uchar flags=0)
{
  BBox3 bbb = b.getBox();
  BBox3 abb = bbb; abb.Move(boffset);
  if(!(flags&CONVOLVE_NOCLIP))
  {
    BBox3 xbb = a.getBox();
    abb.Intersection(xbb);
    bbb = abb;
    bbb.Move(-boffset);
  }

  Array3d<T> aa(a); aa.setSlice(abb);
  ConstArray3d<S> bb(b); bb.setSlice(bbb);
  _ArrayNdInternal::apply(aa, bb, func);
}


DynArray<int> ComputeConvolutionOffsets(const BBox3 &stencil_box, const Int3 &strides, ConstArray3d<bool> mask = Array3d<bool>());

template<class T>
DynArray<T> ComputeConvolutionValues(ConstArray3d<T> &stencil, ConstArray3d<bool> mask = Array3d<bool>());



template<class T, class S, class OpMult = Ops::Mult<S>, class OpCombine = Ops::Assign<S> >
class Convolutor
{
  DynArray<int> filter_offsets;
  DynArray<T> filter_coeff;
  BBox3 filter_box;
  ConstArray3d<S> src;
public:
  OpMult multOp;
  OpCombine combineOp;
  Convolutor() {}
  Convolutor(ConstArray3d<T> filter, ConstArray3d<S> a, bool centered_stencil = false) : src(a)
  {
    myAssert(filter.isContiguous());
    filter_box  = filter.getBox();
    if (centered_stencil)
    {
      Int3 c = (filter.size()-Int3(1))/2;
      filter_box.Move(-c);
    }
    filter_offsets = ComputeConvolutionOffsets(filter_box, src.strides());
    filter_coeff   = ComputeConvolutionValues(filter);
  }

  void Init(ConstArray3d<T> filter_, ConstArray3d<S> a_, bool centered_stencil = false)
  {
    this->~Convolutor<T,S,OpMult,OpCombine>();
    new (this) Convolutor<T,S,OpMult,OpCombine>(filter_, a_, centered_stencil);
  }

  S point_convolve(const Int3 &p) const
  {
    S ret = S();
    const int p_offset = src.offset(p);
    for (int i=0; i<filter_offsets.size(); ++i)
    {
      ret += multOp(filter_coeff[i], src.dereference(p_offset + filter_offsets[i]));
    }
    return ret;
  }

  template<class U>
  void convolve(Array3d<U> dst) const
  {
    myAssert(src.size() - Size(filter_box) + Int3::Constant(1) == dst.size());
    FOR_BBOX3(p, dst.getBox())
    {
      combineOp.inplace(dst(p), point_convolve(p));
    }
  }
};

// template<class S, class T, class U>
// static void convolution(ConstArray3d<S> a, ConstArray3d<T> filter, const Int3 filterOffset, Array3d<U> dst)
// {
//   if (filterOffset != Int3(0))
//     a = a[BBox3(a.getBox()).Move(filterOffset)];
//   Convolutor<T, S>(filter, a).convolve(dst);
// }

template<class S, class T, class U>
static void convolution(ConstArray3d<S> a, ConstArray3d<T> filter, Array3d<U> dst)
{
  Convolutor<T, S, Ops::Mult<S>, Ops::Assign<U> >(filter, a).convolve(dst);
}

template<class S, class T, class U>
static void add_convolution(ConstArray3d<S> a, ConstArray3d<T> filter, Array3d<U> dst)
{
  Convolutor<T, S, Ops::Mult<S>, Ops::Add<U> >(filter, a).convolve(dst);
}

/*
template<class S, class T, class U, class OpMult, class OpCombine>
static void convolution(ConstArray3d<S> a, ConstArray3d<T> filter, const Int3 filterOffset, Array3d<U> dst, const OpMult multOp = Ops::Mult<S>(), const OpCombine combineOp = Ops::Assign<U>())
{
  myAssert(filter.isContiguous());
  const Int3 ma = a.strides();
  const Int3 md = dst.strides();
  const Int3 s = dst.size();
  myAssert(a.size() - filter.size() + Int3::Constant(1) == dst.size());

  ConstArray3d<T> flatFilter = filter; flatFilter.setFlat();
  DynArray<int> foff(flatFilter.size()[0]);
  int i=0;
  FOR_REG3V2(p, filter.size())
  {
    foff[i++] = a.offset(p);
  }

  for(int z=0; z<s.z(); ++z)
  {
    for(int y=0; y<s.y(); ++y)
    {
      const S *pa = a.getPtr() + a.offset(filterOffset.x(), y+filterOffset.y(), z+filterOffset.z());
      const S *paEnd = pa + ma.x() * s.x();
      U* pd = dst.getPtr() + dst.offset(0, y, z);      

      for(; pa < paEnd; pa += ma.x(), pd += md.x())
      {        
        dst.debugCheck(pd);
        S q = S();
        for(i=0; i<foff.size(); ++i)
        {
          a.debugCheck(pa+foff[i]);
          q += multOp(flatFilter[i], pa[foff[i]]);
        }
        combineOp.inplace(*pd, q);
      }
    }
  }
}


template<class S, class StencilGen, class U, class OpMult, class OpCombine>
static void convolution_pos_dependent(ConstArray3d<S> a, const StencilGen &stencilgen, const Int3 filterOffset, Array3d<U> dst, const OpMult multOp = Ops::Mult<S>(), const OpCombine assignOp = Ops::Assign<U>())
{
  const Int3 ma = Int3::Map(a.strides(), 3);
  const Int3 md = Int3::Map(dst.strides(), 3);
  const Int3 s = a.size();
  const Int3 stencilSize = stencilgen.size();

  // allocated memory
  Array3d<typename StencilGen::value_type> stencil(stencilSize);
  const DynArray<typename StencilGen::value_type> flatStencil = stencil.constDynArray();
  // offsets in main array
  DynArray<int> foff(flatStencil.size());
  int i=0;
  FOR_REG3V2(p, stencilSize)
  {
    foff[i++] = a.offset(p);
  }

  Int3 p;
  for(p.z()=0; p.z()<s.z(); ++p.z())
  {
    for(p.y()=0; p.y()<s.y(); ++p.y())
    {
      const S *pa = a.getPtr(filterOffset.x(), p.y()+filterOffset.y(), p.z()+filterOffset.z());
      const S *paEnd = pa + ma.x() * s.x();
      U* pd = dst.getPtr(0, p.y(), p.z());

      for(p.x()=0; pa < paEnd; pa += ma.x(), pd += md.x(), ++p.x())
      {
        dst.debugCheck(pd);
        stencilgen(p, stencil); // update stencil for current position
        S q = S();
        for(i=0; i<foff.size(); ++i)
        {
          a.debugCheck(pa+foff[i]);
          q += multOp(flatStencil[i], pa[foff[i]]);
        }
        assignOp.inplace(*pd, q);
      }
    }
  }
}
*/


template<class T>
Array3d<T> resizedArray3d(ConstArray3d<T> arr, const Int3 &nz, float filterwidth=1.0);


template<class T>
void PrintArray(const ConstArray3d<T> a, std::ostream &os = std::cout);


namespace Steppers { 
  
template<class Q> struct Operations;
  
template<class T>
struct Operations<Array3d<T> >
{
  typedef Array3d<T> State;
  static void addScaled(double fa, State &a, double fb, const State &b) { a.addScaled(fb, b, fa); }
  static void initCloned(State &a, const State &b) { a.initDeepCopy(b); }
  static void initLike(State &a, const State &b) { a.initFromBox(b.getBox()); }
  static void setId(State &a, int id) {}
  static int getId(const State &a) { return -1; }
  static double errorMeasure(double atol, double rtol, const State &y0, const State &y1, const State &err)
  {
      double z = 0.;
      FOR_BBOX3(p, y0.getBox())
      {
        z += my::sqr(1./(rtol*std::max(std::abs(y0(p)), std::abs(y1(p)))+atol) * err(p));
      }
      z = std::sqrt(z / Volume(y0.getBox()));
      return z;
  }
};

}



#endif
