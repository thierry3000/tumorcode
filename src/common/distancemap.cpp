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
#include "distancemap.h"
#include <limits>

DistanceFieldComputer::DistanceFieldComputer() :
  FLAG_OUTER(-1),
  FLAG_INNER(-2),
  DIST_MAX(std::numeric_limits<float>::max()*0.5f),
  heap(1024, ConsTags::RESERVE)
{
}


void DistanceFieldComputer::InsertActive( const Int3 &p )
{
  flags(p)  = heap.size();
  int indx = heap.insert(p);
  myAssert(indx>=0);
  flags(p) = indx;    
}


void DistanceFieldComputer::RemoveActive( const Int3 &p )
{
  int indx = flags(p);
  myAssert(indx>=0);
  heap.remove(indx);
  flags(p)=FLAG_INNER;
}


float DistanceFieldComputer::SolveDistance( float q, float x, float y, float z, bool bx, bool by, bool bz ) const
{
  myAssert(bx || by || bz);
  Float3 coeff(0.f);
  if(bx) {
    myAssert(x<DIST_MAX*0.5f);
    coeff[0] += 1.0f;
    coeff[1] += -2.0f*x;
    coeff[2] += x*x;
  }
  if(by) {
    myAssert(y<DIST_MAX*0.5f);
    coeff[0] += 1.0f;
    coeff[1] += -2.0f*y;
    coeff[2] += y*y;      
  }
  if(bz) {
    myAssert(z<DIST_MAX*0.5f);
    coeff[0] += 1.0f;
    coeff[1] += -2.0f*z;
    coeff[2] += z*z;
  }
  coeff[2] -= q;
  float disc = coeff[1]*coeff[1] - 4.0f*coeff[0]*coeff[2];
  if(disc<0) { return std::numeric_limits<float>::max(); }
  disc = std::sqrt(disc);
  float s0 = 0.5f*(-coeff[1] + disc)/coeff[0];
  float s1 = 0.5f*(-coeff[1] - disc)/coeff[0];  
  float res = std::max(s0,s1);
  //myAssert((!bx || res>x) && (!by || res>y));
  return res;
}
  
  
void DistanceFieldComputer::CheckOkDir(const Int3 &pp, int dir, bool &ok, float &val) const
{
  const Int3 p = ld.NbLattice(pp,dir);
  if(ld.IsInsideLattice(p)) {
    float tmp = f(p);
    val= tmp;
    ok = (std::abs(val)<DIST_MAX*0.5f) ? true : false;
  }
  else 
  {
    val = DIST_MAX;
    ok  = false;
  }  
}


float DistanceFieldComputer::ComputeCenterGradient( const Int3 &p )
{
  float q  = 1.0f * ld.Scale()*ld.Scale();
  Float3 val[2], minval;
  Vec<bool,3> ok[2], bothok, okaxis;
  bool anyok = false;
  int n_pos=0, n_neg=0;
  //float _x[2],_y[2];
  //bool _xok[2],_yok[2];

  for (int axis=0; axis<dim; ++axis)
  {
    for (int side=0; side<2; ++side)
    {
      int dir = axis*2+side;
      CheckOkDir(p, dir, ok[side][axis], val[side][axis]);
//       if (val[side][axis] >= 0.)
//         ++n_pos;
//       else
//       {
//         ++n_neg;
//         val[side][axis] = -val[side][axis];
//       }
    }
    bothok[axis] = ok[0][axis] && ok[1][axis];
    minval[axis] = bothok[axis] ? std::min(val[0][axis],val[1][axis]) : (ok[0][axis] ? val[0][axis] : (ok[1][axis] ? val[1][axis] : std::numeric_limits<float>::max()));
    okaxis[axis] = minval[axis] < DIST_MAX*0.5;
    anyok |= okaxis[axis];
//     anyok &= n_pos == 0 || n_neg == 0;
  }

  if (anyok)
    return SolveDistance(q, minval[0], minval[1], minval[2], okaxis[0], okaxis[1], okaxis[2]);
  else {
    return DIST_MAX;
  }
}


void DistanceFieldComputer::Do(const LatticeDataQuad3d &ld_, Array3d<float> &_f )
{
  f = _f;
  ld = ld_;
  const BBox3 bb = ld.Box();
  dim = ::Size(bb)[2] == 1 ? (::Size(bb)[1] == 1 ? 1 : 2) : 3;
  
  flags.initFromBox(bb, FLAG_OUTER);

  heap.lessfunc.f = f;
  heap.swapfunc.flags = flags;

  DynArray<Int3> inner(1024,ConsTags::RESERVE);
  FOR_BBOX3(p,bb)
  {
    if(std::abs(f(p)) < std::numeric_limits<float>::max()*0.5) {
      flags(p) = FLAG_INNER;
      inner.push_back(p);
    } else {
      f(p) = DIST_MAX;
    }
  }

  for(int i=0; i<inner.size(); ++i)
  {      
    const Int3 p = inner[i];
    //bool insert = false;
    for(int k=0; k<LD::DIR_CNT; ++k)
    {
      const Int3 nb = ld.NbLattice(p,k);
      if(!ld.IsInsideLattice(nb)) continue;
      if(flags(nb)==FLAG_OUTER)
      {
        //insert = true;
        float g = ComputeCenterGradient(nb);
        if (g < DIST_MAX*0.5) // hack! because i dont find the bug where it cannot compute the value from the gradient
          f(nb) = g;
        else
          f(nb) = 0;
        //printf("added %i %i %i value = %f\n",nb.x(),nb.y(), nb.z(),f(nb));
        InsertActive(nb);
      }
    }
    //if (!insert) continue;
    //InsertActive(p);
  }

  clear_and_free_memory(inner);
  
  while(heap.size()>0)
  {
    Int3 p = heap.extremum();
    RemoveActive(p);

   // printf("added %i %i %i value = %f\n",p.x(),p.y(), p.z(),f(p));
    //myAssert(f(p) < DIST_MAX*0.5);

    for(int k=0; k<LD::DIR_CNT; ++k)
    {
      const Int3 nb = ld.NbLattice(p,k);
      if(!ld.IsInsideLattice(nb)) continue;
      if(flags(nb)==FLAG_INNER) continue;
      float g = ComputeCenterGradient(nb);
      if (g < DIST_MAX*0.5)
        f(nb) = g;
      //f(nb) = std::min(f(nb), g);
      //cout << "visited " << nb << " new value is " << f(nb) << endl;
      if(flags(nb)>=0) RemoveActive(nb);
      InsertActive(nb);
    }    
  }
  //GeDebugOut("compute distance field: %i ms",GeGetTimer()-_t);
}

