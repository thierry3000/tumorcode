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
#include "remodeler.h"
#include "../calcflow.h"

#include <boost/scoped_ptr.hpp>


namespace VesselGenerator
{
  
void vesselgen_generate_grid(const H5::Group outgroup, const Int3 &size, float scale, const string &ld_type, float capillary_radius, float press_min, float press_max, const BloodFlowParameters &bfparams)
{
  typedef VesselList3d::LatticeData LatticeData;
  std::auto_ptr<LatticeData> ldp = LatticeData::Make(ld_type.c_str(), BBox3().Add(Int3(0)).Add(size-Int3(1)), scale);

  std::auto_ptr<VesselList3d> vl( new VesselList3d(ldp) );
  vl->Init(*ldp);

  TreeRootList tree_roots;

  const Float3 offset = ldp->GetWorldBox().min;
  const Float3 diagonal = ldp->GetWorldBox().max-ldp->GetWorldBox().min;
  const float one_over_diag_sqrnorm = 1./diagonal.squaredNorm();
  
  const BBox3 bbox = ldp->Box();
  FOR_BBOX3(p, bbox)
  {
    VesselNode* nd = vl->InsertNode(p);
    nd->flags |= CAPILLARY | CIRCULATED | CONNECTED;
    if ((p.array() == bbox.min.array()).any() || (p.array() == bbox.max.array()).any()) // boundary
    {
      nd->flags |= BOUNDARY;
      tree_roots.insert(std::make_pair(p, TreeRoot(p, -1, nd->flags, 0)));
      Float3 wp = ldp->LatticeToWorld(p);
      nd->press = diagonal.dot(wp-offset)*one_over_diag_sqrnorm*(press_max-press_min) + press_min;
    }
  }

  FOR_BBOX3(p, bbox)
  {
    for (int dir = 0; dir<ldp->NbCount(); ++dir)
    {
      const Int3 q = ldp->NbLattice(p, dir);
      if (!ldp->IsInsideLattice(q)) continue;
      if (vl->FindVessel(p, q)) continue;
      Vessel* v = vl->InsertVessel(p, q);
      v->flags |= CAPILLARY | CIRCULATED | CONNECTED;
      v->r = v->reference_r = capillary_radius;
      v->hematocrit = bfparams.inletHematocrit;
    }
  }

  CalcFlow(*vl, bfparams);
  
  {
    DoOutput(outgroup, *vl, tree_roots);
  }
}

void vesselgen_generate_grid_no_flow(const H5::Group outgroup, const Int3 &size, float scale, const string &ld_type, float capillary_radius, float press_min, float press_max, const BloodFlowParameters &bfparams)
{
  typedef VesselList3d::LatticeData LatticeData;
  std::auto_ptr<LatticeData> ldp = LatticeData::Make(ld_type.c_str(), BBox3().Add(Int3(0)).Add(size-Int3(1)), scale);

  std::auto_ptr<VesselList3d> vl( new VesselList3d(ldp) );
  vl->Init(*ldp);

  TreeRootList tree_roots;

  const Float3 offset = ldp->GetWorldBox().min;
  const Float3 diagonal = ldp->GetWorldBox().max-ldp->GetWorldBox().min;
  const float one_over_diag_sqrnorm = 1./diagonal.squaredNorm();
  
  const BBox3 bbox = ldp->Box();
  FOR_BBOX3(p, bbox)
  {
    VesselNode* nd = vl->InsertNode(p);
    nd->flags |= CAPILLARY | CIRCULATED | CONNECTED;
    if ((p.array() == bbox.min.array()).any() || (p.array() == bbox.max.array()).any()) // boundary
    {
      nd->flags |= BOUNDARY;
      tree_roots.insert(std::make_pair(p, TreeRoot(p, -1, nd->flags, 0)));
      Float3 wp = ldp->LatticeToWorld(p);
      nd->press = diagonal.dot(wp-offset)*one_over_diag_sqrnorm*(press_max-press_min) + press_min;
    }
  }

  FOR_BBOX3(p, bbox)
  {
    for (int dir = 0; dir<ldp->NbCount(); ++dir)
    {
      const Int3 q = ldp->NbLattice(p, dir);
      if (!ldp->IsInsideLattice(q)) continue;
      if (vl->FindVessel(p, q)) continue;
      Vessel* v = vl->InsertVessel(p, q);
      v->flags |= CAPILLARY | CIRCULATED | CONNECTED;
      v->r = v->reference_r = capillary_radius;
      v->hematocrit = bfparams.inletHematocrit;
    }
  }

  //CalcFlow(*vl,NULL, bfparams);
  
  {
    DoOutput(outgroup, *vl, tree_roots);
  }
}


enum {
  DIRECTION_STRAIGHT = 0,
  DIRECTION_DIAG_XY = 1,
  DIRECTION_DIAG_XYZ = 2
};


void vesselgen_generate_single(H5::Group outgroup, const Int3 &size, const int direction_mode, float scale, const string &ld_type, float capillary_radius, float press_min, float press_max, int segment_size, const BloodFlowParameters &bfparams)
{
  typedef VesselList3d::LatticeData LatticeData;
  std::auto_ptr<LatticeData> ldp = LatticeData::Make(ld_type.c_str(), BBox3().Add(Int3(0)).Add(size-Int3(1)), scale);

  boost::scoped_ptr<VesselList3d> vl( new VesselList3d(ldp) );
  vl->Init(*ldp);

  auto InsertVessel = [&](VesselNode* nd1, VesselNode* nd2) -> Vessel*
  {
      Vessel* v = vl->InsertVessel(nd1, nd2);
      v->flags |= CAPILLARY | CIRCULATED | CONNECTED;
      v->r = v->reference_r = capillary_radius;
      v->hematocrit = bfparams.inletHematocrit;
      return v;
  };

  auto InsertNode = [&](const Int3 &p) -> VesselNode*
  {
      VesselNode* nd = vl->InsertNode(p);
      nd->flags |= CAPILLARY | CIRCULATED | CONNECTED;
      return nd;
  };
  
  VesselNode* root[2] = { NULL, NULL };
  int dir = 0;
  
//   if (dynamic_cast<polymorphic_latticedata::Derived<LatticeDataQuad3d>*>(ldp.get()) != NULL)
//   {
//     int step = std::max(1, segment_size);
//     Int3 p(0, size[1]/2, size[2]/2);
// 
//     VesselNode* nd_last = NULL;
//     while (true)
//     {
//       VesselNode* nd =InsertNode(p);
//       if (nd_last)
//         InsertVessel(nd_last, nd);
// 
//       if (p[0] == size[0]-1) break;
// 
//       p[0] += step;
//       p[0] = std::min(p[0], size[0]-1);
//       nd_last = nd;
//     }
// 
//     p[0] = 0;
//     root[0] = vl->FindNode(p);
//     p[0] = size[0]-1;
//     root[1] = vl->FindNode(p);
//   }
//   else
//   {
//     const LatticeDataFCC& ld = dynamic_cast<const polymorphic_latticedata::Derived<LatticeDataFCC>*>(ldp.get())->get();
//     int reversedDirs[12];
//     GetReverseDir(ld, reversedDirs);
//     
//     Int3 icenter = size/2;
//     VesselNode* node_center = InsertNode(icenter);
// 
//     dir = 6;
//     switch (direction_mode) // got the numbers from latticedatatest.py
//     {
//       case DIRECTION_DIAG_XY: dir = 8; break;
//       case DIRECTION_DIAG_XYZ: dir = 11; break;
//     }
// 
//     {
//       Int3 p         = icenter;
//       VesselNode* nd = node_center;
//       int length_counter = 0; // how many steps since the last node
//       while (true)
//       {
//         Int3 nextp = ld.NbLattice(p, dir);
//         ++length_counter;
//         
//         if (!ld.IsInsideLattice(nextp)) break;
// 
//         if (length_counter >= segment_size)
//         {
//           VesselNode* nextnd = vl->InsertNode(nextp);
//           InsertVessel(nd, nextnd);
// 
//           length_counter = 0;
//           nd = nextnd;
//         }
// 
//         p  = nextp;
//       }
//       root[0] = nd;
//     }
// 
//     dir = reversedDirs[dir];
// 
//     { // copy & past from the brackets above
//       Int3 p         = icenter;
//       VesselNode* nd = node_center;
//       int length_counter = 0; // how many steps since the last node
//       while (true)
//       {
//         Int3 nextp = ld.NbLattice(p, dir);
//         ++length_counter;
// 
//         if (!ld.IsInsideLattice(nextp)) break;
// 
//         if (length_counter >= segment_size)
//         {
//           VesselNode* nextnd = vl->InsertNode(nextp);
//           InsertVessel(nd, nextnd);
// 
//           length_counter = 0;
//           nd = nextnd;
//         }
//         p  = nextp;
//       }
//       root[1] = nd; // exception: not copy pasted
//     } // end copy pasted
//   }

  if (root[0]->lpos[0] > root[1]->lpos[1])
    std::swap(root[0], root[1]);
  
  root[0]->press = press_max;
  root[1]->press = press_min;

  //cout << "--- vesselgen_grid ---" << endl;
  //ldp->print(cout); cout << endl;
  //cout << "root0 pos " << ldp->LatticeToWorld(root[0]->lpos) << endl;
  //cout << "root1 pos " << ldp->LatticeToWorld(root[1]->lpos) << endl;
  
  TreeRootList tree_roots;
  for (int i=0; i<2; ++i)
  {
    root[i]->flags |= BOUNDARY;
    tree_roots.insert(std::make_pair(root[i]->lpos, TreeRoot(root[i]->lpos, -1, root[i]->flags, 0)));
  }

  CalcFlow(*vl, bfparams);

  {
    DoOutput(outgroup, *vl, tree_roots);
    H5::Group g = outgroup.createGroup("parameters");
//     h5cpp::Group g = outgroup.require_group("parameters");
//     h5cpp::Attributes a = g.attrs();
//     a.set("dir", dir);
    writeAttrToGroup<int>(g, "dir", dir);
  }
}
void vesselgen_generate_symmetric(const H5::Group outgroup, const int &exponent_of_two, float scale, const BloodFlowParameters &bfparams, const bool &only2D)
{
  typedef VesselList3d::LatticeData LatticeData;
  const int points_per_site = pow(2,exponent_of_two)+1;
  Int3 size(points_per_site,points_per_site,points_per_site);
  //not working, do not know why yet???
  //if( only2D )
  if( true )
  {
    size[2] = 1;
  }
  std::auto_ptr<LatticeData> ldp = LatticeData::Make("quad", BBox3().Add(Int3(0)).Add(size-Int3(1)), scale);

  boost::scoped_ptr<VesselList3d> vl( new VesselList3d(ldp) );
  vl->Init(*ldp);

  TreeRootList tree_roots;

  //const Float3 offset = ldp->GetWorldBox().min;
  const Float3 diagonal = ldp->GetWorldBox().max-ldp->GetWorldBox().min;
  const float one_over_diag_sqrnorm = 1./diagonal.squaredNorm();
  
  const BBox3 bbox = ldp->Box();
  FOR_BBOX3(p, bbox)
  {
    VesselNode* nd = vl->InsertNode(p);
    nd->flags |= CAPILLARY | CIRCULATED | CONNECTED;
    if( p[0] == 0 and p[1] == 0 and p[2]==0 )
    {
      //printf("wow");
      nd->flags |= BOUNDARY;
      tree_roots.insert(std::make_pair(p, TreeRoot(p, -1, nd->flags, 0)));
      //Float3 wp = ldp->LatticeToWorld(p);
      nd->press = 3.;
    }
    if( p[0] == points_per_site-1 and p[1] == points_per_site-1 and p[2]==points_per_site-1 )
    {
      //printf("wow");
      nd->flags |= BOUNDARY;
      tree_roots.insert(std::make_pair(p, TreeRoot(p, -1, nd->flags, 0)));
      //Float3 wp = ldp->LatticeToWorld(p);
      nd->press = 5.;
    }
    
//     if ((p.array() == bbox.min.array()).any() || (p.array() == bbox.max.array()).any()) // boundary
//     {
//       nd->flags |= BOUNDARY;
//       tree_roots.insert(std::make_pair(p, TreeRoot(p, -1, nd->flags, 0)));
//       Float3 wp = ldp->LatticeToWorld(p);
//       double press_max=6.;
//       double press_min=3.;
//       nd->press = diagonal.dot(wp-offset)*one_over_diag_sqrnorm*(press_max-press_min) + press_min;
//     }
  }
  fflush(stdout);
  FOR_BBOX3(p, bbox)
  {
    for (int dir = 0; dir<ldp->NbCount(); ++dir)
    {
      const Int3 q = ldp->NbLattice(p, dir);
      if (!ldp->IsInsideLattice(q)) continue;
      if (vl->FindVessel(p, q)) continue;
      Vessel* v = vl->InsertVessel(p, q);
      v->flags |= CAPILLARY | CIRCULATED | CONNECTED;
      v->r = v->reference_r = 3.;
      v->hematocrit = bfparams.inletHematocrit;
    }
  }

  CalcFlow(*vl, bfparams);
  
  {
    DoOutput(outgroup, *vl, tree_roots);
  }
}

}//end Vesselgenerator
