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
#ifndef GLOBALS_H_
#define GLOBALS_H_

#include "helpers-vec.h"
#include "helpers-mem.h"


template<int d, class SiteType = int, class StrideType = SiteType>
class LatticeIndexing // Uniform lattice -> integer type index or pointer
{
public:
  typedef BBoxd<int, d> IBox;
  typedef Vec<int, d> IVec;
  typedef Vec<StrideType, d> SVec;
  enum { dim = d };

private:
  IVec l;
  SVec strides;
  IBox box;
  SiteType offset;
  // see http://eigen.tuxfamily.org/dox-devel/TopicTemplateKeyword.html
  StrideType dotStrides(const IVec &x) const
  {
    return strides.dot(x.template cast<StrideType>());
  }

public:
  const IVec& Size() const { return l; }
  const IBox& Box() const { return box; }
  StrideType NumSites() const { return strides[dim-1]*l[dim-1]; }
  const SVec& Strides() const { return strides; }
  SiteType SiteAtOrigin() const { return offset; }
  
  SiteType SiteAtBoxMin() const { return offset + dotStrides(box.min); }
  SiteType SiteAtBoxMax() const { return offset + dotStrides(box.max); }

  static SVec CalcStrides(const IVec &l)
  {
    SVec r;
    r[0] = 1;
    for (int i=1; i<dim; ++i)
      r[i] = l[i-1] * r[i-1];
    return r;
  }

  enum StartSiteIndexMode  {
    AT_ORIGIN,
    AT_BBOX_MIN
  };
  
  LatticeIndexing() : offset(0) {}

  void Init( const IBox &bb_, const SiteType site_idx = 0, const SVec &strides_ = SVec(-1), StartSiteIndexMode start_site_mode = AT_ORIGIN)
  {
    box = bb_;
    l = ::Size(bb_);
    strides = strides_;
    if (strides[0] == -1)
      strides = CalcStrides(l);
    offset = site_idx;
    if (start_site_mode == AT_BBOX_MIN)
      // site_idx equals offset + bb.min.dot(strides), hence offset = site_idx - bb.min.dot(strides)
      offset -= dotStrides(box.min);
  }

  void Init(const IVec &l_, const SiteType site_idx_at_origin = 0, const SVec &strides_ = IVec(-1))
  {
    offset = site_idx_at_origin;
    l = l_;
    box.min = IVec(0);
    box.max = l - IVec(1);
    strides = strides_;
    if (strides[0] == -1)
      strides = CalcStrides(l);
  }

  SiteType LatticeToSite( const IVec &p ) const
  {
    // see http://eigen.tuxfamily.org/dox-devel/TopicTemplateKeyword.html
    return dotStrides(p) + offset;
  }

  const IVec SiteToLattice( SiteType site_ ) const
  {
    myAssert(strides == CalcStrides(l));
    IVec r;
    StrideType site = site_ - SiteAtBoxMin();
    for (int i=dim-1; i>=1; --i)
    {
      r[i] = site/strides[i] + box.min[i];
      site %= strides[i];
    }
    r[0] = site + box.min[0];
    return r;
  }

   bool IsInsideLattice( const Int3 &p ) const { return box.Overlaps(p); }
   bool IsInsideLattice( SiteType site ) const { return IsInsideLattice(SiteToLattice(site)); }

   void SwapAxes(int i, int j)
   {
     myAssert(i>=0 && i<dim && j>=0 && j<dim);
     std::swap(strides[i], strides[j]);
     std::swap(l[i], l[j]);
     for (int s=0; s<=1; ++s)
      std::swap(box[s][i], box[s][j]);
   }

  void SetBox(const BBox3 &b)
  {
    box.min = b.min;
    box.max = b.max;
    l = ::Size(b);
  }

  void MoveSiteRange(const Int3 p)
  {
    offset -= dotStrides(p);
  }
};



struct AxisDirLen
{
  int dir;
  uint len;
  bool isValid() const { return dir>=0; }
  AxisDirLen() : dir(-1), len(0) {}
  AxisDirLen(int dir, uint len) : dir(dir), len(len) {}
};




template<int dim_>
struct LatticeWorldTransform
{
  enum { dim = dim_ };
protected:
  typedef Vec<int, dim> IVec;
  typedef Vec<float, dim> FVec;
  FVec wo, centering;
  float scale, scale_inv;

public:
  LatticeWorldTransform() : scale(0.) {}
  LatticeWorldTransform(float scale) { Init(scale); }
  void Init(float scale);
  
  // Origin refers to index vector zero. Hence OriginPositin is the position of index zero in world coordinates.
  void SetOriginPosition(const FVec &pos);
  FVec GetOriginPosition() const { return wo; }
  
  void SetCellCentering(const Vec<bool, dim> &cc);
  Vec<bool, dim_> GetCellCentering() const;

   float Scale() const { return scale; }
   void Scale(float s) { scale = s; scale_inv = 1./scale; }
  float ScaleInv() const { return scale_inv; }
  
  FVec LatticeToWorld( const IVec &p ) const
  {
    // centering == 0.5 means cell centered coordinate, 0 is vertex centered
    FVec r;
    for (int i=0; i<dim; ++i)
    {
      r[i] = ((int(p[i]) + centering[i]) * scale) + wo[i];
    }
    return r;
  }

  IVec WorldToLattice(const FVec &wp) const
  {
    // Cubic regions around lattice sites(cells) are mapped to the corresponding integer coordinate.
    return WorldToLatticeCell(wp).first;
  }

  /* this gives the index vector of the closes cell/vertex and
     a fractional coordinate of the location within the cell/vertex
     so that 0 is the left border, and 1 is the right border of the
     cell/dual cell around the vertex  */
  std::pair<IVec,FVec> WorldToLatticeCell(const FVec &wp) const;

  /* this gives the index vector of the closest cell/vertex to the left of the
   * world position wp (denoted x0) and furthermore a fractional coordinate
   * representing the distance to the neighbor cell center/vertex to the right
   * (denoted x1), where 0 is exactly on x0 and 1 is exactly on x1 */
  std::pair<IVec,FVec> WorldToFractionalCoordinate(const FVec &wp) const;
  
  BBoxd<float, dim_> GetWorldBox(const BBoxd<int, dim> &lbb) const;
};





struct LatticeDataQuad3d : public LatticeIndexing<3, int64, int64>, public LatticeWorldTransform<3>
{
  typedef LatticeIndexing<3,int64,int64> LI;
  typedef LatticeWorldTransform<3> LWT;
  typedef Int3 LatticeIndexType;
  typedef int64 SiteType;
  //Float3 wo; // world coordinate offset, added in lattice -> world
  //Float3 centering;

  enum LatDirQuad3d  // direction form pos1 to pos2 for use with nn[]
  {
    DIR_MX  = 0,  // left
    DIR_PX  = 1, // right
    DIR_MY  = 2,
    DIR_PY  = 3,
    DIR_MZ  = 4,
    DIR_PZ  = 5,
    DIR_CNT = 6,
    DIR_CNT2D = 4
  };    
  
  LatticeDataQuad3d();
  LatticeDataQuad3d(const LatticeDataQuad3d &ld) : LI(), LWT() { CopyMem(&ld, this, 1); }
  LatticeDataQuad3d(const Int3 &l, float scale=1.0f ) { Init(l, scale); }
  LatticeDataQuad3d(const BBox3 &bb, float scale=1.0f ) { Init(bb, scale); }

  void Init( const Int3 &l, float scale=1.0f);
  void Init( const BBox3 &bb, float scale=1.0);

  int NbCount() const { return DIR_CNT; }
  SiteType NbSite( SiteType site, int dir ) const 	{ return site + nb[dir]; }
  SiteType NbSite(  const Int3 &p, int dir ) const { return LatticeToSite(p)+nb[dir]; }
  Int3 NbLattice( const Int3 &p, int dir ) const { return p + vnb[dir]; }

  FloatBBox3 GetWorldBox() const { return LWT::GetWorldBox(LI::Box()); }

  static void DirToAxisSide(int dir, int &axis, int &side);
  static int ReversedDir( int dir ) { return (dir&1) ? (dir-1) : (dir+1); }

  AxisDirLen GetAxisDirLen(const Int3 &p1, const Int3 &p2 ) const;
  Int3 GetLatticeIndexOnRefinedGrid(const Int3 &pos, int refinement_subdivision) const;

  void print(std::ostream &os) const;

protected:
  // i dont bother with thread safety. If two LatticeDatas happen to
  // be constructed at the same time, they will just initialize this
  // stuff twice with the same data, which is not harmful
  static volatile bool nbs_init;
  int nb[DIR_CNT];
  static Int3 vnb[DIR_CNT];
  //float scale, scale_inv; // lattice spacing and 1/spacing
};



struct LatticeDataFCC : public LatticeIndexing<3, int64, int64>
{
  typedef LatticeIndexing<3,int64> Base;
  typedef Int3 LatticeIndexType;
  typedef int64 SiteType;
  enum { DIR_CNT = 12 };
  
  Float3 wo; // world coordinate offset, added in lattice -> world

  LatticeDataFCC() { ClearMem( this, 1 ); }
  LatticeDataFCC(const LatticeDataFCC &ld) { CopyMem(&ld, this, 1); }
  LatticeDataFCC(const BBox3 &bb, float scale=1.0f ) { Init(bb, scale); }
  void Init( const BBox3 &bb, float scale=1.0);

  void SetOriginPosition(const Float3 &pos) { wo = pos; }
  Float3 GetOriginPosition() const { return wo; }
  
  FloatBBox3 GetWorldBox() const;
  float Scale() const { return scale; }
  void Scale(float s) { scale = s; scale_inv = 1.f/scale; }

  int NbCount() const { return DIR_CNT; }
  SiteType NbSite(SiteType site, int dir) const;
  SiteType NbSite(const Int3 &p, int dir) const;
  Int3 NbLattice( const Int3 &p, int dir ) const;
  Float3 LatticeToWorld( const Int3 &p ) const;

  AxisDirLen GetAxisDirLen(const Int3 &p1, const Int3 &p2 ) const;
  Int3 GetLatticeIndexOnRefinedGrid(const Int3 &pos, int refinement_subdivision) const;

  
  void print(std::ostream &os) const;

protected:
  static volatile bool nbs_init;
  // cases y&1=0,1 * mod(z,3)=0,1,2
  int nb[2][3][DIR_CNT];
  static Int3 vnb[2][3][DIR_CNT];
  float scale, scale_inv; // lattice spacing and 1/spacing
};



enum {
  GET_WORLD_DIRECTION_AXES_NORMALIZE = 1,
};

template<class LD>
inline void GetWorldDirectionAxes(const LD &ld, Float3 *wdelta, int flags=0)
{
  const int N = ld.NbCount();
  const Float3 w0 = ld.LatticeToWorld(Int3(0));
  for (int i=0; i<N; ++i)
  {
    wdelta[i] = ld.LatticeToWorld(ld.NbLattice(Int3(0), i)) - w0;
    if (flags & GET_WORLD_DIRECTION_AXES_NORMALIZE)
      wdelta[i].normalize();
  }
}

template<class LD>
inline void GetReverseDir(const LD &ld, int *dir)
{
  Float3 wdelta[128];
  GetWorldDirectionAxes(ld, wdelta, GET_WORLD_DIRECTION_AXES_NORMALIZE);
  for (int i=0; i<ld.NbCount(); ++i)
  {
    for (int j=i+1; j<ld.NbCount(); ++j)
    {
      float d = wdelta[i].dot(wdelta[j]) + 1.;
      if (d*d > 1.e-3) continue;
      dir[i] = j;
      dir[j] = i;
    }
  }
}


bool operator==( const LatticeDataQuad3d& a, const LatticeDataQuad3d& b );
bool operator==( const LatticeDataFCC& a, const LatticeDataFCC& b );
inline std::ostream& operator<<(std::ostream &os, const LatticeDataQuad3d &ld) { ld.print(os); return os; }
inline std::ostream& operator<<(std::ostream &os, const LatticeDataFCC &ld) { ld.print(os); return os; }


#endif
