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


#include "shared-objects.h"
#include "vessels3d.h"
#include "mwlib/math_ext.h"
#include "boost/foreach.hpp"

#include <memory>
#include <csignal>


//#define DEBUG_CHECKS_ON

// hahaha where dafuq is this from???!!!!
float GetInitialThickness( float rad )
{
  const float d = 2.0f*rad;
  return 0.65*d - 0.2f*log10f( d )*d;
}

/* thickness table from Marieb Human Anatomy & Physiology
 * diameter [um],  thickness [um]
 * 15000           1000              elastic artery
 * 6000            1000              muscular artery
 * 37              6                 arteriole
 * 9               0.5               capillary
 * 20              1                 venule
 * 5000            500               vein
 */


float Vessel::WorldLength( const polymorphic_latticedata::LatticeData &ld ) const
{
  return (len-1)*ld.Scale();
}

// VesselList3d::VesselList3d()
//   : m_ld(nullptr) //,bclist(nullptr)
// {
//   this->bclist.reset(new BCList);
// }
/* a forced constructor */
// VesselList3d::VesselList3d(std::unique_ptr< VesselList3d::LD > &this_ld)
// {
//   this->bclist.clear();
//   this->m_ld = std::move(this_ld);
//   this->Flush();
//   this->lookup_site.Init(*m_ld);
//   this->lookup_bond.Init(*m_ld);
//   this->g.Reserve( 1024 );
// }
VesselList3d::VesselList3d::VesselList3d()
{
}

// std::auto_ptr< polymorphic_latticedata::LatticeData > VesselList3d::getLD() const
// {
//   return m_ld;
// }

void VesselList3d::ClearLattice()
{
  delete m_ld.get();
}


std::unique_ptr< VesselList3d > VesselList3d::Clone()
{
  #ifdef DEBUG
  cout << "this lattice data will be passed to: " << endl;
  this->m_ld->print(cout);
  #endif
  std::unique_ptr<VesselList3d> my_return_list(new VesselList3d());
  //my_return_list.reset(new VesselList3d);
//   my_return_list->Init(this->Ld());
  my_return_list->Init(m_ld);
  //my_return_list->bclist.reset(new BCList);
  #ifdef DEBUG
  cout << "recognized: " << endl;
  my_return_list->m_ld->print(cout);
#endif
  //my_return_list->Init(Ld());
//   #ifdef DEBUG
//   cout << "recognized: " << endl;
//   my_return_list->Ld().print(cout);
// #endif
  //my_return_list->bclist.reset(new BCList);
  for( int i=0;i<GetNCount();i++)
  {
    const VesselNode *p_currentNode = GetNode(i);
    VesselNode *newNode = my_return_list->InsertNode(p_currentNode);//to new vessellist
    newNode->flags = p_currentNode->flags;
    for(auto entry: GetBCMap())
    {
      if(p_currentNode == entry.first)
      {
	my_return_list->SetBC(newNode, entry.second);
      }
    }

  }
  for( int i=0;i<GetECount();i++)
  {
    const Vessel *p_currentVessel = GetEdge(i);
    my_return_list->InsertVessel(p_currentVessel);
  }
  
  return my_return_list;
}

// void VesselList3d::Init( const LD &_ld )
// {
//     this->bclist.clear();
//     this->Flush();
//     this->m_ld = _ld.Clone();
//     //this->m_ld.reset(&_ld);
//     this->lookup_site.Init(*m_ld);
//     this->lookup_bond.Init(*m_ld);
//     this->g.Reserve( 1024 );
// }
void VesselList3d::Init( std::unique_ptr<LD> &_ld )
{
    this->bclist.clear();
    this->Flush();
    this->m_ld = std::move(_ld);
    //this->m_ld.reset(&_ld);
    this->lookup_site.Init(*m_ld);
    this->lookup_bond.Init(*m_ld);
    this->g.Reserve( 1024 );
}

void VesselList3d::Flush()
{
  this->lookup_site.Clear();
  this->lookup_bond.Clear();
  this->g.Flush();
}

void VesselList3d::SetBC(const VesselNode* node, FlowBC &fbc)
{
  // see if the node is in the list of managed nodes.
  //Don't want dangling nodes here.
  myAssert(GetNode(node->Index()) == node); 
  myAssert(node->IsBoundary());
  
  /** I dont know why yet, but emplace is not doing the job */
  //this->bclist.emplace(std::make_pair(node,fbc));
  
  bclist[node] = fbc;
//   for( auto it = GetBCMap().begin(); it!= GetBCMap().end(); ++it)
//   {
//     if( it->first->Index() == node->Index() )
//     {
//       
//     }
// //     cout << "here?" << endl;
// //     cout << "index: " << it->first->Index() << "bcvalue: " << it->second.val <<endl;
//   }
}


void VesselList3d::ClearBC(VesselNode* node)
{
  myAssert(GetNode(node->Index()) == node); // see if the node is in the list of managed nodes. Don't want dangling nodes here.
  myAssert(node->IsBoundary());
  this->bclist.erase(node);
}


VesselNode* VesselList3d::InsertNode( const Int3 &a ) 
{
  myAssert(Ld().IsInsideLattice(a));
#ifdef DEBUG_CHECKS_ON
  {
    Vessel* vtest;
    VesselNode* vctest;
    FindAny(a, vtest, vctest);
    myAssert(vtest == NULL && vctest == NULL);
  }
#endif
    VesselNode* n = g.InsertNode();
    n->lpos = a;
    lookup_site.InsertAtSite(Ld().LatticeToSite(a),n);
    return n;
}

VesselNode* VesselList3d::InsertNode( const Float3 &a ) 
{
  //myAssert(Ld().IsInsideLattice(a));
#ifdef DEBUG_CHECKS_ON
  {
    Vessel* vtest;
    VesselNode* vctest;
    FindAny(a, vtest, vctest);
    myAssert(vtest == NULL && vctest == NULL);
  }
#endif
    VesselNode* n = g.InsertNode();
    n->worldpos = a;
    return n;
}
VesselNode* VesselList3d::InsertNode(const VesselNode *p_n)
{
  VesselNode *newNode 		= this->InsertNode(p_n->lpos);
  newNode->press 		= p_n->press;
  newNode->residual 		= p_n->residual;
  newNode->worldpos 		= p_n->worldpos;
  newNode->bctyp 		= p_n->bctyp;
  newNode->value_of_bc 		= p_n->value_of_bc;
  newNode->flags 		= p_n->flags;
  newNode->has_been_visited 	= p_n->has_been_visited;
  //if pn in boundarymap
  // add also to this->bclist;
  return newNode; 
}
void VesselList3d::InsertVessel(const Vessel* p_v)
{
  Vessel *newVessel 		= this->InsertVessel(p_v->NodeA()->lpos,p_v->NodeB()->lpos);
  newVessel->r 			= p_v->r;
  newVessel->flags 		= p_v->flags;
  newVessel->q 			= p_v->q;
  newVessel->f 			= p_v->f;
  newVessel->hematocrit 	= p_v->hematocrit;
  newVessel->conductivitySignal = p_v->conductivitySignal;
  newVessel->metabolicSignal 	= p_v->metabolicSignal;
  newVessel->S_total 		= p_v->S_total;
  newVessel->maturation 	= p_v->maturation;
  newVessel->reference_r 	= p_v->reference_r;
  newVessel->timeSprout 	= p_v->timeSprout;
  newVessel->timeInTumor 	= p_v->timeInTumor;
  
  newVessel->len 		= p_v->len;
  newVessel->dir 		= p_v->dir;
}


Vessel* VesselList3d::InsertVessel( VesselNode* na, VesselNode* nb )
{
  myAssert( na!=NULL && nb!=NULL );
  Vessel* v;
  if(HasLattice())
  {
    v = g.InsertEdge( na, nb );
#ifdef DEBUG
    //printf("na->index: %i, nb->index: %i\n", na->Index(), nb->Index());
#endif
    const AxisDirLen x = Ld().GetAxisDirLen(na->lpos, nb->lpos);
    myAssert(x.dir>=0);
    myAssert(x.len>0 );
    v->len = x.len+1;
    v->dir = x.dir;
#ifdef DEBUG_CHECKS_ON
    {
      Vessel* vtest;
      VesselNode* vctest;
      Int3 plast, p = na->lpos;
      for (int i=0; i<v->len-1; ++i)
      {
        if (i>0 && i<v->len-1)
        {
          Vessel* vtest;
          VesselNode* vctest;
          FindAny(p, vtest, vctest);
          myAssert(vtest == NULL && vctest == NULL);
        }
        if (i>0)
        {
          vtest = FindVessel(p, plast);
          myAssert(vtest == NULL);
        }
        plast = p;
        p = Ld().NbLattice(p, v->dir);
      }
    }
#endif
    lookup_bond.Insert(Ld().LatticeToSite(na->lpos),x.dir,x.len,v);
  }
  else //this was InsertVesselWorld
  {
    v = g.InsertEdge( na, nb );
#ifdef DEBUG
    printf("na->index: %i, nb->index: %i\n", na->Index(), nb->Index());
#endif
  }
  return v;    
}



Vessel* VesselList3d::InsertVessel( const Int3 &a, const Int3 &b )
{
    myAssert(Ld().IsInsideLattice(a) && Ld().IsInsideLattice(b));
    VesselNode* na = lookup_site.FindAtSite(Ld().LatticeToSite(a));
    VesselNode* nb = lookup_site.FindAtSite(Ld().LatticeToSite(b));
    if( !na ) 
      na = InsertNode(a);
    if( !nb ) 
      nb = InsertNode(b);
    return InsertVessel( na, nb );
}

void VesselList3d::DeleteVessel( Vessel* v, bool bDeleteNodes ) 
{
  if(HasLattice())
  {
    SiteType sa = Ld().LatticeToSite(v->LPosA());
    SiteType sb = Ld().LatticeToSite(v->LPosB());
    const AxisDirLen x = Ld().GetAxisDirLen(v->LPosA(), v->LPosB());
    lookup_bond.Remove(sa, x.dir, x.len);
    
    v->NodeA()->RemoveEdge(v);
    v->NodeB()->RemoveEdge(v);
    if( bDeleteNodes && v->NodeA()->Count()<=0 )
    {
      DeleteUnusedNode(v->NodeA(), sa); // note: NOT calling ListGraph::DeleteUnusedNode! Instead this version also takes care of node removal from BCMap. See VesselList3d::DeleteUnusedNode
    }
    if( bDeleteNodes && v->NodeB()->Count()<=0 )
    {
      DeleteUnusedNode(v->NodeB(), sb);
    }
    g.DeleteEdgeButDoNotTouchNodePointers(v);
  }
  else//this was DeleteVesselWorld
  {
    g.DeleteEdge(v, bDeleteNodes);
  }
}


Vessel* VesselList3d::FindVessel( const Int3 &a, const Int3 &b ) 
{ 
    myAssert(Ld().IsInsideLattice(a) && Ld().IsInsideLattice(b));
    const AxisDirLen x = Ld().GetAxisDirLen(a, b);
    return lookup_bond.Find(Ld().LatticeToSite(a), x.dir, x.len);
}

VesselNode* VesselList3d::FindNode( const Int3 &a ) 
{ 
    myAssert(Ld().IsInsideLattice(a));
    return lookup_site.FindAtSite(Ld().LatticeToSite(a));
}

bool VesselList3d::FindAny( const Int3 &a, Vessel* &v, VesselNode* &vc )
{
    SiteType site = Ld().LatticeToSite(a);
    vc = lookup_site.FindAtSite(site);
    v  = NULL;
    if (vc) return true;
    for( int i=0; i<Ld().NbCount(); ++i )
    {
      v = lookup_bond.Find(site, i, 1);
      if( v ) return true;
    }
    return false;
}

void VesselList3d::DeleteUnusedNode(VesselNode* vc, int site)
{
  myAssert(vc->Count()==0);
  lookup_site.RemoveAtSite(site);
  // remove node from list of boundary conditions
  if (vc->IsBoundary())
  {
    auto it = bclist.find(vc);
    //myAssert (it != bclist.end());  // this is bullshit. Of course, the vessel does not have to be in the BC list because if it is not, pressure BCs are simply assumed.
    if (it != bclist.end())
      bclist.erase(it);
  }
#ifdef DEBUG
  myAssert(bclist.find(vc) == bclist.end()); // not a boundary node -> not in the list of boundary nodes
#endif
  g.DeleteUnusedNode(vc);
}

void VesselList3d::DeleteUnusedNode(VesselNode* vc)
{
  if(HasLattice())
  {
    DeleteUnusedNode(vc, Ld().LatticeToSite(vc->lpos));
  }
  else//this was DeleteUnusedNodeWorld
  {
    myAssert(vc->Count()==0);
    g.DeleteUnusedNode(vc);
    vc = NULL;
  }
}

FloatBBox3 VesselList3d::GetWorldBoxFromVesselsOnly()
{
  FloatBBox3 r;
  //const BBox3 b = Box();
  //BBox3 q;
  for (int ax=0; ax<3; ++ax)
  {
    r.min[ax] = 1000000000000.0;
    r.max[ax] = -100000000000.0;
    for(int i=0;i<this->GetNCount();++i)
    {
      VesselNode *nd = this->GetNode(i);
      if(nd->worldpos[ax]<r.min[ax])
      {
	r.min[ax] = nd->worldpos[ax];
      }
      if(nd->worldpos[ax]>r.max[ax])
      {
	r.max[ax] = nd->worldpos[ax];
      }
    }
  }
//   FOR_BBOX3(p, q)
//   {
//     r.Add(LatticeToWorld(b.min + p));
//     r.Add(LatticeToWorld(b.max - p));
//   }
  return r;
}



void VesselList3d::SplitVessel( Vessel* v, int pos, Vessel* &vnew, VesselNode* &vcnew )
{
  myAssert(pos>0 && pos<v->len-1);
  Int3 lpos = v->LPosA();
  for( int i=0; i<pos; ++i )
  {
    lpos=Ld().NbLattice(lpos,v->dir);
  }
  myAssert(Ld().IsInsideLattice(lpos));

  SiteType sitea = Ld().LatticeToSite(v->LPosA());
  //int sitepos = Ld().LatticeToSite(pos);
  
  lookup_bond.Remove(sitea, v->dir, v->len-1);
  
  // store pointers
  VesselNode* vc1 = v->NodeA();
  VesselNode* vc2 = v->NodeB();
  // insert node, also into lookup table
  vcnew = InsertNode(lpos);
  // correct links
  vc1->RemoveEdge(v);
  vc2->RemoveEdge(v);
  v->AttachNodes(vc1,vcnew);
  v->len = pos+1;

  lookup_bond.Insert(sitea, v->dir, v->len-1, v);
  
//   const AxisDirLen x = GetAxisDirLen(ld, vc1->lpos,vcnew->lpos);
//   v->len = x.len+1;
//   v->dir = x.dir;
  // insert vessel also into lookup table, hopefully overwriting the old vessel entry
  vnew  = InsertVessel(vcnew,vc2);
  *static_cast<VData*>(vnew) = *static_cast<VData*>(v);
}


#define CHECK(exp) (void)( (exp) || (_Assert(#exp, __FILE__, __LINE__), 0) )

void VesselList3d::IntegrityCheck(int check_lookup_)
{
  bool check_lookup = check_lookup_==1; // -1 = default = false, 0 = false, 1 = true
  
  //lookup_bond.PrintContents(cout);
  //lookup_site.PrintContents(cout);
  int num_sites = Volume(Ld().Box());
  DynArray<uchar> markers(num_sites);
  for( int i=0; i<GetECount(); ++i )
  {
    Vessel* v = GetEdge(i);
    if (check_lookup) // check the bonds in the lookup table
    {
      CHECK(v->len >= 2);
      Int3 p = v->LPosA();
      for( int j=0; j<v->len; ++j, p=Ld().NbLattice(p,v->dir) )
      {
        Vessel* vtest; VesselNode* vctest;
        bool bfound = FindAny( p, vtest, vctest );
        CHECK( bfound );
        CHECK( j!=0 || (vctest && vctest==v->NodeA()) );
        CHECK( j!=v->len-1 || (vctest && vctest==v->NodeB()) );
        CHECK( j==0 || j==v->len-1 || vtest==v );
        markers[Ld().LatticeToSite(p)] = 1;
      }
    }
    // check edge lists
    CHECK(v->NodeA()->FindEdgeIndex(v)>=0);
    CHECK(v->NodeB()->FindEdgeIndex(v)>=0);
  }
  for (int i=0; i<GetNCount(); ++i)
  {
    VesselNode* vc = GetNode(i);
    if (check_lookup) // mark occupied sites
    {
      Int3 p = vc->lpos;
      CHECK(Ld().IsInsideLattice(p));
      markers[Ld().LatticeToSite(p)] = 1;
    }
    for (int j=0; j<vc->Count(); ++j) 
    {
      CHECK(vc->GetEdge(j)->NodeA()==vc || vc->GetEdge(j)->NodeB()==vc); // check edge list
      for (int k=j+1; k<vc->Count(); ++k) // double edge
      {
        CHECK(vc->GetEdge(k)->GetOther(vc) != vc->GetEdge(j)->GetOther(vc));
      }
    }
  }

  if (check_lookup) // check if occupied sites are really in the lookup table
  {
    for( int i=0; i<num_sites; ++i )
    {
      if( markers[i] ) continue;
      Int3 p = Ld().SiteToLattice(i);
      Vessel* vtest; VesselNode* vctest;
      bool bfound = FindAny( p, vtest, vctest );
      CHECK( bfound == false );
    }
  }
}

#undef CHECK

void VesselList3d::FillLookup()
{
  lookup_site.Init(Ld());
  lookup_bond.Init(Ld());
  for( int i=0; i<GetNCount(); ++i )
  {
    lookup_site.InsertAtSite(Ld().LatticeToSite(GetNode(i)->lpos), GetNode(i) );
  }
  for( int i=0; i<GetECount(); ++i )
  {
    Vessel* v = GetEdge(i);
    const AxisDirLen x = Ld().GetAxisDirLen(v->NodeA()->lpos, v->NodeB()->lpos);
    lookup_bond.Insert(Ld().LatticeToSite(v->NodeA()->lpos), x.dir, x.len, v);
  }
}

// std::size_t VesselList3d::estimateMemoryUsage() const
// {
//   std::size_t sz1 = lookup_site.estimateMemoryUsage(),
//               sz2 = lookup_bond.estimateMemoryUsage(),
//               sz3 = g.estimateMemoryUsage(),
//               sz4 = sizeof(*this);
//   return sz1+sz2+sz3+sz4;
// }


/*------------------------------------------------------
------------------------------------------------------*/

void HemodynamicBounds::Add(const VesselList3d *vl, bool bClear)
{
  if( bClear )
  {
    press.clear();
    flow.clear();
    velocity.clear();
    force.clear();
  }

  {
    int ecnt = vl->GetECount();
    for( size_t i=0; i<ecnt; ++i )
    {
      const Vessel* v = vl->GetEdge(i);
      force.add( v->f );
      flow.add( v->q );
      velocity.add( v->q/(v->r*v->r*my::mconst::fpi()) );
      press.add( v->NodeA()->press );
      press.add( v->NodeB()->press );
    } 
  }
}

//-----------------------------------
//-----------------------------------

// namespace polymorphic_latticedata {
//   std::auto_ptr<LatticeData> GetSubdivided(const LatticeData &ld, int multi, int additional_boundary)
//   {
//     BBox3 newbox;
//     newbox.min = ld.GetLatticeIndexOnRefinedGrid(ld.Box().min, multi-1);
//     newbox.min = ld.GetLatticeIndexOnRefinedGrid(ld.Box().min, multi-1);
//     newbox.Extend(additional_boundary);
//     
//     std::auto_ptr<LatticeData> newldp(ld.Clone());
//     newldp->Init(newbox,  ld.Scale() / multi);
//     return newldp;
//   };
//   
// }


void GetSubdivided( std::unique_ptr<VesselList3d> &vl, int multi, float newscale, int safety_boundary)
{
  typedef VesselList3d::LatticeData LatticeData;
  const LatticeData &ld = vl.get()->Ld();
  
#ifndef NDEBUG
  std::cout << " in GetSubdivided" << std::endl;
#endif
  //if(multi == 1) return vl;

  int ecnt = vl->GetECount();
  int ncnt = vl->GetNCount();
  
  BBox3 newbox;
  DynArray<Int3> newpos(ncnt);
  for (int i=0; i<ncnt; ++i)
  {
    newpos[i] = ld.GetLatticeIndexOnRefinedGrid(vl->GetNode(i)->lpos, multi-1);
    newbox.Add(newpos[i]);
  }
  newbox.Extend(safety_boundary);

  std::unique_ptr<LatticeData> newldp = polymorphic_latticedata::Make_ld("FCC", vl->Ld().Box(), vl->Ld().Scale());
  //std::unique_ptr<LatticeData> newldp(ld.Clone());
  newldp.get()->Init(newbox,  newscale);

  std::unique_ptr<VesselList3d> vlnew( new VesselList3d() );
  //vlnew->Init(*newldp);
  vlnew->Init(newldp);

  // insert new nodes into empty new vessellist
  for (int i=0; i<ncnt; ++i)
  {
    VesselNode* vcold = vl->GetNode(i);
    VesselNode* vc = vlnew->InsertNode(newpos[i]);
    *(VNodeData*)vc = *(VNodeData*)vcold;
    myAssert(vc->Index() == vcold->Index());
  }

  for (int i=0; i<ecnt; ++i)
  {
    const Vessel* v = vl->GetEdge(i);
    int a = v->NodeA()->Index(),
        b = v->NodeB()->Index();
    Vessel* newv = vlnew->InsertVessel(vlnew->GetNode(a), vlnew->GetNode(b));
    *(VData*)newv = *(VData*)v;
  }
#ifndef NDEBUG
  std::cout << "start adding vl->GetBCMap() : " << vl->GetBCMap().size() << std::endl;
#endif
  for (auto it = vl->GetBCMap().begin(); it != vl->GetBCMap().end(); ++it)
  {
    FlowBC aCondition = FlowBC(it->second.typeOfInstance, it->second.val);
    vlnew->SetBC(vlnew->GetNode(it->first->Index()), aCondition);
    //std::cout << "added condition: " << aCondition.val << " at " << it->first->Index() << std::endl;
  }
  //vl->ClearLattice();
  vl = std::move(vlnew);
#ifndef NDEBUG
  std::cout << "moved pointer" << std::endl;
#endif
  //return vlnew;
}


void GetSubdivided(std::unique_ptr<VesselList3d> &vl, float scale)
{
  const int multi = std::max(int( my::round( vl->Ld().Scale()/scale ) ),1);
  if( multi == 1)
  {
    //return vl;
  }
  else
  {
    GetSubdivided(vl, multi, scale);
  }
}


void SplitSegmentsToOneLatticeBond(VesselList3d &vl)
{
  using LatticeData = VesselList3d::LatticeData; // the c++11 way!
  const LatticeData &ld = vl.Ld();
  int ecnt = vl.GetECount();
  std::vector<Vessel*> to_process; to_process.reserve(ecnt);
  for (int i=0; i<ecnt; ++i)
  {
    Vessel* v = vl.GetEdge(i);
    if (v->len <= 2) continue;
    to_process.push_back(v);
  }
  
  for (Vessel* v: to_process)
  {
    VData data_backup = *v;
    auto *nodea = v->NodeA();
    auto *nodeb = v->NodeB();
    int num_bonds = v->len - 1;
    int dir = v->dir;
    vl.DeleteVessel(v, false); // keep attached nodes
    
    VesselNode* node = nodea;
    Int3 pos1 = nodea->lpos;
    for (int i=0; i<num_bonds; ++i)
    {
      Int3 pos2 = ld.NbLattice(pos1, dir);
      myAssert(vl.FindNode(pos2) == nullptr || vl.FindNode(pos2) == nodeb);
      Vessel* vnew = vl.InsertVessel(pos1, pos2);
      node = vnew->GetOther(node);
      if (node != nodeb)
      {
	*static_cast<VNodeData*>(node) = *nodea;
	node->flags = node->flags & nodeb->flags;
	node->flags.DelBits(BOUNDARY);
      }
      *static_cast<VData*>(vnew) = data_backup;
      pos1 = pos2;
    }
  }
}


//-----------------------------------
//-----------------------------------
static inline float EffectiveRadius( float r1, float r2, float l1, float l2 )
{
  return r1*r2*std::pow((l1+l2)/(my::sqr(my::sqr(r1))*l2+my::sqr(my::sqr(r2))*l1),0.25f);
}

// try to merge two small segments
uint Optimize( VesselList3d *vl )
{
  int reversedir[128];
  GetReverseDir(vl->Ld(), reversedir);
  //FUNC_TIMING_START
  uint num_nodes2 = 0;
  uint num_nodesgr2 = 0;
  uint num_merged = 0;
  static const int LEN_THRES = 5;
  static const float RAD_THRES = 10;
  static const int TIME_SPROUT_THRES = 0;
  static const int TIME_IN_TUMOR_THRES = 1;
  for( int i=0; i<vl->GetNCount(); ++i )
  {
    VesselNode* nd = vl->GetNode(i);
    uint num_perf_vess = 0;
    for( int j=0; j<nd->Count(); ++j ) if(nd->GetEdge(j)->IsCirculated()) ++num_perf_vess;
    if( num_perf_vess==2 ) ++num_nodes2;
    else if( num_perf_vess>2 ) ++num_nodesgr2;
    if( nd->Count()!=2 ) continue;
    Vessel* v1 = nd->GetEdge(0);
    Vessel* v2 = nd->GetEdge(1);
    if( v1->dir!=v2->dir && v1->dir!=reversedir[v2->dir]) continue;
    if( v1->len>LEN_THRES || v2->len>LEN_THRES ) continue;
    if( v1->flags!=v2->flags ) continue;
    if( std::abs(v2->r-v1->r)>RAD_THRES ) continue;
    if( std::abs(v2->timeInTumor-v1->timeInTumor)>TIME_IN_TUMOR_THRES ) continue;
    if( std::abs(v2->timeSprout-v1->timeSprout)>TIME_SPROUT_THRES ) continue;

    VData c;
    int l1 = v1->len-1;
    int l2 = v2->len-1;
    c.r           = EffectiveRadius(v1->r,v2->r,v1->WorldLength(vl->Ld()),v2->WorldLength(vl->Ld()));
    c.q           = 0.5f*(v1->q+v2->q);
    c.f           = 0.5f*(v1->f+v2->f);
    c.flags       = v1->flags;
    c.initial_f   = 0.5f*(v1->initial_f+v2->initial_f);
    c.maturation  = 0.5f*(v1->maturation+v2->maturation);
    c.hematocrit  = 0.5f*(v1->hematocrit+v2->hematocrit);
    //c.viscos      = 0.5f*(v1->viscos+v2->viscos);
    c.timeInTumor = std::max(v1->timeInTumor,v2->timeInTumor);
    c.timeSprout  = std::min(v1->timeSprout,v2->timeSprout);
    c.reference_r  = ((v1->reference_r/v1->r)+(v2->reference_r/v2->r))*0.5f*c.r;
    VesselNode* nd1 = v1->GetOther(nd);
    VesselNode* nd2 = v2->GetOther(nd);
    vl->DeleteVessel(v1,false);
    vl->DeleteVessel(v2,false);
    vl->DeleteUnusedNode(nd);
    Vessel* vnew = vl->InsertVessel( nd1, nd2 );
    *static_cast<VData*>(vnew) = c;
    ++num_merged;
  }
  
  vl->OptimizeMemLayout();

  //printf("num_nodes2 = %i, num_nodesgr2 = %i\n",num_nodes2,num_nodesgr2);
  //printf("num_merged = %u\n",num_merged);
  return num_merged;
  //FUNC_TIMING_END
}

//-----------------------------------
//-----------------------------------

void BondLookup::PrintContents(std::ostream &os) const
{
  typedef Map::const_iterator It;
  os << "Bond Lookup Table = [" << std::endl;
  for (It it = map.begin(); it != map.end(); ++it)
  {
    K k = it->first;
    os << boost::format("%s,%s -> %s") % ld().SiteToLattice(k.first) % ld().SiteToLattice(k.second) % it->second << std::endl;
  }
  os << "]" << std::endl;
}

void BondLookup::Insert(S site, int dir, int len, V x)
{
  S nbsite;
  for( int i=0; i<len; ++i )
  {
    nbsite=ld().NbSite(site,dir);
    map[make_key(site, nbsite)] = x;
    site = nbsite;
  }
}


void SiteLookup::PrintContents(std::ostream &os) const
{
  typedef Map::const_iterator It;
  os << "Site Lookup Table = [" << std::endl;
  for (It it = map.begin(); it != map.end(); ++it)
  {
    os << boost::format("%s -> %s") % ld().SiteToLattice(it->first) % it->second << std::endl;
  }
  os << "]" << std::endl;
}

//-----------------------------------
//-----------------------------------

// sort vessels based on pressure, so that blood flows through vessels in the list in ascending order 
// void TopoSortVesselsRecursive(VesselList3d &vl, VesselNode* n0, DynArray<VesselNode*> &sorted, DynArray<uchar> &visited )
// {
//   visited[n0->Index()] = 1;
//   for( int k=0; k<n0->Count(); ++k )
//   {
//     NodeNeighbor<VesselNode*,Vessel*> nb( n0->GetNode( k ) );
//     if( visited[nb.node->Index()] || !nb.edge->IsCirculated() || nb.node->press<=n0->press ) continue; // goto perfused upstream segments only
//     TopoSortVesselsRecursive(vl, nb.node, sorted, visited );
//     //sorted.push_back( nb.v );
//   }
//   sorted.push_back( n0 );
// }
// 
// 
// void TopoSortVessels( VesselList3d &vl, DynArray<VesselNode*> &sorted )
// {
//   uint time = GetTimer();
//   DynArray<uchar> visited( vl.GetNCount(), 0 );
//   sorted.reserve( vl.GetECount() );
//   for( int i=0; i<vl.GetNCount(); ++i )
//   {
//     if( visited[i] ) continue;
//     TopoSortVesselsRecursive( vl, vl.GetNode(i), sorted, visited);
//   }
//   printf("topo sort time %u ms\n",GetTimer()-time);
// }



void TopoSortVesselsRecursive(const VesselList3d &vl,const VesselNode* n0, const Vessel* v0, int &current_idx, DynArray<int> &order, DynArray<uchar> &visited )
{
  visited[n0->Index()] = 1;
  int n = 0;
  for( int k=0; k<n0->Count(); ++k )
  {
    NodeNeighbor<const VesselNode*, const Vessel*> nb(n0->GetNode(k));
    if (!nb.edge->IsCirculated()) continue;
    if (nb.edge == v0) continue;
    if (nb.node->press <= n0->press) continue;

    if (!visited[nb.node->Index()])
      TopoSortVesselsRecursive( vl, nb.node, nb.edge, current_idx, order, visited );
    myAssert(order[nb.edge->Index()]==-1);
    order[nb.edge->Index()] = current_idx++;
    ++n;
  }
}


/* Order so that order[vessel->index] increases in down stream direction.
 */
void TopoSortVessels(const VesselList3d& vl, DynArray< int >& order )
{
  uint time = GetTimer();

  DynArray<uchar> visited( vl.GetNCount(), 0 );
  order.resize(vl.GetECount(), -1);
  int current_idx = 0;
  
  for( int i=0; i<vl.GetNCount(); ++i )
  {
    if( visited[i] ) continue;
    TopoSortVesselsRecursive( vl, vl.GetNode(i), nullptr, current_idx, order, visited );
  }

  // uncirculated vessels or vessels with equal pressure at the end nodes are processed last
  for (int i=0; i<vl.GetECount(); ++i)
  {
    int idx = vl.GetEdge(i)->Index();
    if (order[idx] < 0)
      order[idx] = current_idx++;
  }
  
  //printf("topo sort time %u ms\n",GetTimer()-time);
}


void CheckToposort(const VesselList3d &vl, const DynArray<int> &order)
{
  for (int i=0; i<vl.GetNCount(); ++i)
  {
    const VesselNode* nd = vl.GetNode(i);
    for (int k=0; k<nd->Count(); ++k)
    {
      NodeNeighbor<const VesselNode*,const Vessel*> nbk(nd->GetNode(k));
      if (!nbk.edge->IsCirculated()) continue;
      bool inflowk = nbk.node->press>nd->press;
      for (int l=0; l<nd->Count(); ++l)
      {
        //bool err = false;
        NodeNeighbor<const VesselNode*,const Vessel*> nbl(nd->GetNode(l));
        bool outflowl = nbl.node->press<nd->press;
        if (!nbl.edge->IsCirculated()) continue;
        
        if (inflowk && outflowl)
        {
          int vidx_in, vidx_out;
          vidx_out = nbl.edge->Index();
          vidx_in  = nbk.edge->Index();
          if (order[vidx_in] >= order[vidx_out])
          {
            fprintf(stderr, "error: order[%i]=%i < order[%i]=%i but flow goes from %i to %i\n", vidx_out, order[vidx_out], vidx_in, order[vidx_in], vidx_in, vidx_out);
            fprintf(stderr, "%i flags = %u, %i flags = %u\n", vidx_in, (uint)vl.GetEdge(vidx_in)->flags, vidx_out, (uint)vl.GetEdge(vidx_out)->flags);
            fprintf(stderr, "%i-pressures = %e, %e\n", vidx_in, vl.GetEdge(vidx_in)->NodeA()->press-nd->press, vl.GetEdge(vidx_in)->NodeB()->press-nd->press);
            fprintf(stderr, "%i-pressures = %e, %e\n", vidx_out, vl.GetEdge(vidx_out)->NodeA()->press-nd->press, vl.GetEdge(vidx_out)->NodeB()->press-nd->press);
            fprintf(stderr, "%i-flow = %e\n", vidx_in, vl.GetEdge(vidx_in)->q);
            fprintf(stderr, "%i-flow = %e\n", vidx_out, vl.GetEdge(vidx_out)->q);
            AssertAlways(false);
          }
        }
        //AssertAlways(!err);
      }
    }
  }
}
