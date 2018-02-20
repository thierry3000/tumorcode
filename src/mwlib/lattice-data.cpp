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

#include "lattice-data.h"
#include "common/hdfio.h"

#include <float.h>

#include "math_ext.h"
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/export.hpp>

#define LatticeWorldTransform_TEMPL_ARG1 template<int dim_>
#define LatticeWorldTransform_CLASS LatticeWorldTransform<dim_>

LatticeWorldTransform_TEMPL_ARG1
void LatticeWorldTransform_CLASS::Init(float scale_)
{
  Scale(scale_);
}


LatticeWorldTransform_TEMPL_ARG1
void LatticeWorldTransform_CLASS::SetOriginPosition(const FVec &pos)
{
  wo = pos;
}


LatticeWorldTransform_TEMPL_ARG1
void LatticeWorldTransform_CLASS::SetCellCentering(const Vec<bool, dim> &cc)
{
  for (int i=0; i<dim; ++i)
    centering[i] = cc[i] ? 0.5-std::numeric_limits<float>::epsilon() : 0.;
}


LatticeWorldTransform_TEMPL_ARG1
Vec<bool, dim_>
  LatticeWorldTransform_CLASS::GetCellCentering() const
{
  return cwiseG(centering, 0.);
}



LatticeWorldTransform_TEMPL_ARG1
std::pair<typename LatticeWorldTransform_CLASS::IVec, typename LatticeWorldTransform_CLASS::FVec>
  LatticeWorldTransform_CLASS::WorldToLatticeCell(const FVec &wp) const
{
  FVec q;
  IVec ip;
  for (int i=0; i<dim; ++i)
  {
    q[i] = (wp[i] - wo[i]) * scale_inv - centering[i] + 0.5f;
    ip[i] = my::ifloor(q[i]);
    q[i] -= float(ip[i]);
  }
  return std::make_pair(ip, q);
}

LatticeWorldTransform_TEMPL_ARG1
std::pair<typename LatticeWorldTransform_CLASS::IVec, typename LatticeWorldTransform_CLASS::FVec>
LatticeWorldTransform_CLASS::WorldToFractionalCoordinate(const FVec &wp) const
{
  FVec q;
  IVec ip;
  for (int i=0; i<dim; ++i)
  {
    q[i] = (wp[i] - wo[i]) * scale_inv - centering[i]; // the difference to WorldToLatticeCell is that the 0.5 is missing here. It makes a big different in the floor operation.
    ip[i] = my::ifloor(q[i]);
    q[i] -= float(ip[i]);
  }
  return std::make_pair(ip, q);
}


LatticeWorldTransform_TEMPL_ARG1
BBoxd<float, dim_>
  LatticeWorldTransform_CLASS::GetWorldBox(const BBoxd<int, LatticeWorldTransform_CLASS::dim> &lbb) const
{
  BBoxd<float, dim> r;
  r.min = LatticeToWorld(lbb.min);
  r.max = LatticeToWorld(lbb.max);
  for (int i=0; i<dim; ++i)
  {
    r.min[i] -= centering[i] * scale;
    r.max[i] += centering[i] * scale;
  }
  return r;
}


template struct LatticeWorldTransform<1>;
template struct LatticeWorldTransform<2>;
template struct LatticeWorldTransform<3>;


volatile bool LatticeDataQuad3d::nbs_init = false;
Int3 LatticeDataQuad3d::vnb[LatticeDataQuad3d::DIR_CNT];

LatticeDataQuad3d::LatticeDataQuad3d()
{
  ClearMem( this, 1 );
}


void LatticeDataQuad3d::Init( const Int3 &l_, float scale_)
{
  Init(BBox3(Int3(0), l_-Int3(1)), scale_);
}


void LatticeDataQuad3d::Init(const BBox3 &bb, float scale_)
{
  type = string("QUAD3D");
  myAssert(!bb.IsEmpty());
  myAssert(scale_ > 0.);

  LI::Init(bb, 0);
  LI::MoveSiteRange(bb.min);
  LWT::Init(scale_);
  SetOriginPosition(Float3(0));

  if (!nbs_init)
  {
    int k =0;
    for (int ax=0; ax<3; ++ax)
    {
      for (int side=-1; side<=1; side+=2)
      {
        Int3 q(0); q[ax]=side;
        vnb[k++] = q;
      }
    }
    nbs_init = true;
  }

  for (int i=0; i<DIR_CNT; ++i)
  {
    nb[i] = LatticeToSite(NbLattice(Int3(0,0,0),i));
  }
}


void LatticeDataQuad3d::DirToAxisSide(int dir, int &axis, int &side)
{
  axis = dir/2;
  side = dir%2;
}


void LatticeDataQuad3d::print(std::ostream &os) const
{
  os << "[LatticeDataQuad3d " << std::endl;
  os << "  box:    " << Box()  << std::endl;
  os << "  size:   " << Size() << " x " << scale << std::endl;
  os << "  stride: " << Strides() << std::endl;
  os << "  offset: " << SiteAtBoxMin() << std::endl;
  os << "  cell-centering: " << Bool3(cwiseG(centering, 0.f)) << std::endl;
  os << "  world-offset of <0,0,0>: " << wo << std::endl;
  os << "  world-box: " << GetWorldBox() << std::endl;
  os << "]";
}


AxisDirLen LatticeDataQuad3d::GetAxisDirLen(const Int3 &p1, const Int3 &p2 ) const
{
  typedef LatticeDataQuad3d Ld;
  //myAssert( ld.IsInsideLattice( p1 ) && ld.IsInsideLattice( p2 ) && p1!=p2 );
  const bool beq[3] = { p1.x()==p2.x(), p1.y()==p2.y(), p1.z()==p2.z() };
  const int  dirs[3] = { Ld::DIR_MX, Ld::DIR_MY, Ld::DIR_MZ };
  AxisDirLen r;
  for( std::size_t i=0; i<3; ++i )
  {
    if( !beq[i] )
    {
      myAssert( beq[(i+1)%3]==true && beq[(i+2)%3]==true );
      r.dir = dirs[i];
      if( p1[i] < p2[i] ) ++r.dir; // go to minus dir by default, switch to plus if needed
      r.len = std::abs(p1[i]-p2[i]);
      return r;
    }
  }
  //myAssert((!"Not aligned to bonds"));
  return r;
}


Int3 LatticeDataQuad3d::GetLatticeIndexOnRefinedGrid(const Int3& pos, int refinement_subdivision) const
{
  return pos * (1 + refinement_subdivision);
}

volatile bool LatticeDataFCC::nbs_init = false;
Int3 LatticeDataFCC::vnb[2][3][LatticeDataFCC::DIR_CNT];

void LatticeDataFCC::Init(const BBox3 &bb, float scale_)
{
  type = string("FCC");
  myAssert(!bb.IsEmpty());
  myAssert(scale_ > 0.);

  Base::Init(bb, 0);
  Base::MoveSiteRange(bb.min);

  scale = scale_;
  scale_inv = 1./scale;

  SetOriginPosition(Float3(0));

  if (!nbs_init)
  {
    // auto generated with python script which iterates over all possible nearest neighbors and extracts the relative indices
    int k=0;
    vnb[0][0][k++] = Int3( 0,-1,-1);
    vnb[0][0][k++] = Int3( 1,-1,-1);
    vnb[0][0][k++] = Int3( 0, 0,-1);
    vnb[0][0][k++] = Int3( 0,-1, 0);
    vnb[0][0][k++] = Int3( 1,-1, 0);
    vnb[0][0][k++] = Int3(-1, 0, 0);
    vnb[0][0][k++] = Int3( 1, 0, 0);
    vnb[0][0][k++] = Int3( 0, 1, 0);
    vnb[0][0][k++] = Int3( 1, 1, 0);
    vnb[0][0][k++] = Int3( 0,-1, 1);
    vnb[0][0][k++] = Int3(-1, 0, 1);
    vnb[0][0][k++] = Int3( 0, 0, 1);
    k=0;
    vnb[0][1][k++] = Int3( 0, 0,-1);
    vnb[0][1][k++] = Int3( 1, 0,-1);
    vnb[0][1][k++] = Int3( 1, 1,-1);
    vnb[0][1][k++] = Int3( 0,-1, 0);
    vnb[0][1][k++] = Int3( 1,-1, 0);
    vnb[0][1][k++] = Int3(-1, 0, 0);
    vnb[0][1][k++] = Int3( 1, 0, 0);
    vnb[0][1][k++] = Int3( 0, 1, 0);
    vnb[0][1][k++] = Int3( 1, 1, 0);
    vnb[0][1][k++] = Int3( 1,-1, 1);
    vnb[0][1][k++] = Int3( 0, 0, 1);
    vnb[0][1][k++] = Int3( 1, 0, 1);
    k=0;
    vnb[0][2][k++] = Int3(-1, 0,-1);
    vnb[0][2][k++] = Int3( 0, 0,-1);
    vnb[0][2][k++] = Int3( 0, 1,-1);
    vnb[0][2][k++] = Int3( 0,-1, 0);
    vnb[0][2][k++] = Int3( 1,-1, 0);
    vnb[0][2][k++] = Int3(-1, 0, 0);
    vnb[0][2][k++] = Int3( 1, 0, 0);
    vnb[0][2][k++] = Int3( 0, 1, 0);
    vnb[0][2][k++] = Int3( 1, 1, 0);
    vnb[0][2][k++] = Int3( 0, 0, 1);
    vnb[0][2][k++] = Int3( 0, 1, 1);
    vnb[0][2][k++] = Int3( 1, 1, 1);
    k=0;
    vnb[1][0][k++] = Int3(-1,-1,-1);
    vnb[1][0][k++] = Int3( 0,-1,-1);
    vnb[1][0][k++] = Int3( 0, 0,-1);
    vnb[1][0][k++] = Int3(-1,-1, 0);
    vnb[1][0][k++] = Int3( 0,-1, 0);
    vnb[1][0][k++] = Int3(-1, 0, 0);
    vnb[1][0][k++] = Int3( 1, 0, 0);
    vnb[1][0][k++] = Int3(-1, 1, 0);
    vnb[1][0][k++] = Int3( 0, 1, 0);
    vnb[1][0][k++] = Int3(-1,-1, 1);
    vnb[1][0][k++] = Int3(-1, 0, 1);
    vnb[1][0][k++] = Int3( 0, 0, 1);
    k=0;
    vnb[1][1][k++] = Int3( 0, 0,-1);
    vnb[1][1][k++] = Int3( 1, 0,-1);
    vnb[1][1][k++] = Int3( 0, 1,-1);
    vnb[1][1][k++] = Int3(-1,-1, 0);
    vnb[1][1][k++] = Int3( 0,-1, 0);
    vnb[1][1][k++] = Int3(-1, 0, 0);
    vnb[1][1][k++] = Int3( 1, 0, 0);
    vnb[1][1][k++] = Int3(-1, 1, 0);
    vnb[1][1][k++] = Int3( 0, 1, 0);
    vnb[1][1][k++] = Int3( 0,-1, 1);
    vnb[1][1][k++] = Int3( 0, 0, 1);
    vnb[1][1][k++] = Int3( 1, 0, 1);
    k=0;
    vnb[1][2][k++] = Int3(-1, 0,-1);
    vnb[1][2][k++] = Int3( 0, 0,-1);
    vnb[1][2][k++] = Int3(-1, 1,-1);
    vnb[1][2][k++] = Int3(-1,-1, 0);
    vnb[1][2][k++] = Int3( 0,-1, 0);
    vnb[1][2][k++] = Int3(-1, 0, 0);
    vnb[1][2][k++] = Int3( 1, 0, 0);
    vnb[1][2][k++] = Int3(-1, 1, 0);
    vnb[1][2][k++] = Int3( 0, 1, 0);
    vnb[1][2][k++] = Int3( 0, 0, 1);
    vnb[1][2][k++] = Int3(-1, 1, 1);
    vnb[1][2][k++] = Int3( 0, 1, 1);
    
    nbs_init = true;
  }

  for (int y=0; y<2; ++y)
    for (int z=0; z<3; ++z)
    {
      for (int k=0; k<DIR_CNT; ++k)
        nb[y][z][k] = LatticeToSite(vnb[y][z][k] + Int3(0,y,z)) - LatticeToSite(Int3(0,y,z));
    }
}


LatticeDataFCC::SiteType LatticeDataFCC::NbSite(const Int3 &p, int dir) const
{
  return LatticeToSite(p) + nb[p[1]&1][my::mod(p[2],3)][dir];
  
}


LatticeDataFCC::SiteType LatticeDataFCC::NbSite(SiteType site, int dir) const
{
  return LatticeToSite(NbLattice(SiteToLattice(site), dir));
  //return NbSite(SiteToLattice(site), dir);
}


Int3 LatticeDataFCC::NbLattice( const Int3 &p, int dir ) const
{
  return p + vnb[p[1]&1][my::mod(p[2],3)][dir];
  
}

static Int3 fcc2cartesian(const Int3 &p)
{
  /* transform hcp lattice indices to indices for non-uniform cartesian lattice
   * with 1/6 the lattice spacing as the hcp lattice
   */
  Int3 r(p);
  r *= 6;
  if (p[1]&1) { r[0]-=3; }
  if (my::mod(p[2],3)==1) { r[0]+=3; }
  r[1] += (my::mod(p[2],3))*2;
  return r;
}


static Int3 cartesian2fcc(const Int3 &q)
{
  Int3 f, r;
  const int divs[3] = { 6,12,18 };
  for (int i=0; i<3; ++i) {
    f[i] = my::mod(q[i], divs[i]);
    r[i] = q[i]/divs[i] - (q[i]<0 && f[i]!=0 ? 1 : 0);
  }
  r[1] *= 2;
  r[2] *= 3;
  //std::cout << q << " " << f << " " << r;
  if (f == Int3(3,6,0))
    r += Int3(1,1,0);
  else if (f == Int3(3,2,6))
    r += Int3(0,0,1);
  else if (f == Int3(0,8,6))
    r += Int3(0,1,1);
  else if (f == Int3(0,4,12))
    r += Int3(0,0,2);
  else if (f == Int3(3,10,12))
    r += Int3(1,1,2);
  else {
    assert(f == Int3(0));
  }
  //std::cout << " " <<  r << std::endl;
  return r;
}



Float3 LatticeDataFCC::LatticeToWorld( const Int3 &p ) const
{
  Int3 q = fcc2cartesian(p);
  Float3 w;
  float f = scale/6.;
  w[0] = q[0]*f;
  w[1] = q[1]*f*0.8660254037844386;  // sqrt(3)/2
  w[2] = q[2]*f*0.8164965809277259; // sqrt(6)/3
  w += wo;
  return w;
}


FloatBBox3 LatticeDataFCC::GetWorldBox() const
{
  FloatBBox3 r;
  const BBox3 b = Box();
  BBox3 q;
  for (int ax=0; ax<3; ++ax)
  {
    q.min[ax] = 0;
    q.max[ax] = std::min(b.max[ax]-b.min[ax], 2);
  }
  FOR_BBOX3(p, q)
  {
    r.Add(LatticeToWorld(b.min + p));
    r.Add(LatticeToWorld(b.max - p));
  }
  return r;
}


void LatticeDataFCC::print(std::ostream &os) const
{
  os << "[LatticeDataHCP " << std::endl;
  os << "  box:    " << Box()  << std::endl;
  os << "  size:   " << Size() << " x " << scale << std::endl;
  os << "  stride: " << Strides() << std::endl;
  os << "  offset: " << SiteAtBoxMin() << std::endl;
  os << "  world-offset of <0,0,0>: " << wo << std::endl;
  os << "  world-box: " << GetWorldBox() << std::endl;
  os << "]";
}



AxisDirLen LatticeDataFCC::GetAxisDirLen(const Int3 &p1, const Int3 &p2 ) const
{
  /*
   * given two lattice indices, find the direction index dir and
   * the number of times one has to march in this direction
   * in order to get from p1 to p2.
   * Uses integer vector math to probe all possible lattice axes.
   */
  Int3 q1 = fcc2cartesian(p1),
       q2 = fcc2cartesian(p2);
  //std::cout << p2 << std::endl;
  Int3 v = q2-q1;
  for (int d=0; d<DIR_CNT; ++d)
  {
    Int3 w = fcc2cartesian(NbLattice(p1, d)) - q1;
    int a = w.dot(v);
    int lv2 = v.squaredNorm();
    int lw2 = w.squaredNorm();
    if (lv2*lw2 == a*a && a>0) // vectors are parallel
    {
      int l = w[0]!=0 ? v[0]/w[0] : v[1]/w[1];
      return AxisDirLen(d, l);
    }
  }
  return AxisDirLen();
}


Int3 LatticeDataFCC::GetLatticeIndexOnRefinedGrid(const Int3 &pos, int refinement_subdivision) const
{
  Int3 q = fcc2cartesian(pos);
  q *= 1 + refinement_subdivision; // 1 subdiv means the edge length is halved, n subdivs means edge length = l_coarse / (1+subdivs)
  return cartesian2fcc(q);
}




bool operator==( const LatticeDataQuad3d& a, const LatticeDataQuad3d& b )
{
  return a.Box()==b.Box() && a.Scale()==b.Scale() && a.LatticeToSite(a.Box().min) == b.LatticeToSite(b.Box().min);
}

bool operator==( const LatticeDataFCC& a, const LatticeDataFCC& b )
{
  return a.Box()==b.Box() && a.Scale()==b.Scale() && a.LatticeToSite(a.Box().min) == b.LatticeToSite(b.Box().min);
}

void LatticeDataFCC::WriteHdfLd(H5::Group& g) const
{
  //string latticeType = getType();
  writeAttrToH5(g, "TYPE", getType());
  writeAttrToH5(g,"SIZEX", Size()[0]);
  writeAttrToH5(g,"SIZEY", Size()[1]);
  writeAttrToH5(g,"SIZEZ", Size()[2]);
  writeAttrToH5(g, "SIZE", Size());

  BBox3 bb = Box();
  // box is stored in a 6 component vector, xxyyzz, must match python code
  Int6 bv;
  //Int6 bv;
  for (int i=0; i<3; ++i)
  {
    bv[i*2  ] = bb.min[i];
    bv[i*2+1] = bb.max[i];
  }
  writeAttrToH5(g, string("BOX"), bv);
  writeAttrToH5(g, string("WORLD_OFFSET"),GetOriginPosition() );
  writeAttrToH5(g, string("SCALE"), Scale());
}

void LatticeDataQuad3d::WriteHdfLd(H5::Group& g) const
{
  Bool3 centering = GetCellCentering();
  writeAttrToH5(g, "CENTERING", centering);
  
  writeAttrToH5(g, "TYPE", getType());
  writeAttrToH5(g,"SIZEX", Size()[0]);
  writeAttrToH5(g,"SIZEY", Size()[1]);
  writeAttrToH5(g,"SIZEZ", Size()[2]);
  writeAttrToH5(g, "SIZE", Size());

  BBox3 bb = Box();
  // box is stored in a 6 component vector, xxyyzz, must match python code
  Int6 bv;
  //Int6 bv;
  for (int i=0; i<3; ++i)
  {
    bv[i*2  ] = bb.min[i];
    bv[i*2+1] = bb.max[i];
  }
  writeAttrToH5(g, string("BOX"), bv);
  writeAttrToH5(g, string("WORLD_OFFSET"),GetOriginPosition() );
  writeAttrToH5(g, string("SCALE"), Scale());
}
