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
#ifndef VESSELS3D_H_
#define VESSELS3D_H_

#include <unordered_map>
#include <boost/unordered_map.hpp>
#include "lattice-data-polymorphic.h"
#include "mwlib/listgraph.h"
//#include "shared-objects.h"
#include "mwlib/mempool.h"
#include "mwlib/math_ext.h"
#include "calcflow.h"

struct HemodynamicBounds;
class VesselList3d;
// namespace polymorphic_latticedata {
//   class LatticeData;
// }


#define DEFINE_PROPERTY( name, varname, Type ) \
struct Prop##name {  \
  typedef Type value_type; \
  inline Type operator()( const ThisType* x ) const { return x->varname; }  \
  inline void operator()( ThisType* x, const Type q ) const { x->varname=q; } \
};


enum Flags
{
  BOUNDARY = 1<<0, // 1
  CIRCULATED = 1<<1, // 2
  CONNECTED  = 1<<2, // 4
  ARTERY = 1<<3, // 8
  VEIN = 1<<4, // 16
  CAPILLARY = 1<<5, // 32
  WITHIN_TUMOR = 1<<6, // 64
  HIDDEN = 1<<7 // 128
};

enum
{
  //VN_LPOS = 1, // Int3
  VN_WPOS,     // Float3
  VN_PRESS,    // float

  //V_INT_LEN,   // int
  V_WORLD_LEN, // float
  V_MVR,    // float
  V_FLOW,      // float
  V_FORCE,     // float
  V_HEMATO,    // float
  V_FLAGS,     // uint
  V_VISCOS,    // float
  V_TIME_IN_TUMOR, // float
  V_TIME_SPROUT, // float
  V_MATURATION, // float
  V_DEFORMED_R, // float
  V_INITIAL_F,  // float
  V_REFERENCE_R // float
};
enum {
  MAX_NUM_ADJACENT_EDGES = 6 // need to reserve more than 3 for the biconnected components algorithm
};

class Vessel;
class VesselNode;

struct VNodeData
{
  typedef VNodeData ThisType;
  VNodeData() : press(0),  worldpos(0,0,0) {}
  double press, residual;
  Float3 worldpos;
  uint bctyp;
  double value_of_bc;
  my::Bitfield<uchar> flags;
  bool has_been_visited;
  DEFINE_PROPERTY( Press, press, double )
  bool IsBoundary() const { return flags.GetBits(BOUNDARY); }
  bool IsCirculated() const { return flags.GetBits(CIRCULATED);}
  void SetWorldPos( Float3 a) { this->worldpos = a;}
  template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
	ar & press;
	ar & residual;
	ar & worldpos;
	ar & bctyp;
	ar & value_of_bc;
	ar & flags;
	ar & has_been_visited;
    } 
};

struct VData
{
  typedef VData ThisType;
  VData() : r(5.0f),q(0),f(0),/*viscos(0),*/hematocrit(0.45),maturation(0),initial_f(0),timeSprout(0),timeInTumor(0),reference_r(0) {}
  float r;
  my::Bitfield<uchar> flags;
  double q,f,/*viscos,*/hematocrit;
  double conductivitySignal, metabolicSignal, S_total;
  float maturation,initial_f,reference_r;
  int timeSprout,timeInTumor;
  bool IsCirculated() const { return flags.GetBits(CIRCULATED); }
  bool IsArtery() const { return flags.GetBits(ARTERY); }
  bool IsVein() const { return flags.GetBits(VEIN); }
  bool IsCapillary() const { return flags.GetBits(CAPILLARY); }
  DEFINE_PROPERTY( Flow, q, double )
  DEFINE_PROPERTY( Force, f, double )
  DEFINE_PROPERTY( Hematocrit, hematocrit, double )
  DEFINE_PROPERTY( Radius, r, float )
  DEFINE_PROPERTY( TimeSprout,timeSprout,int )
  DEFINE_PROPERTY( TimeInTumor,timeInTumor,int )
  DEFINE_PROPERTY( Maturation,maturation,float )
  DEFINE_PROPERTY( InitialF,initial_f,float )
  DEFINE_PROPERTY( ReferenceR,reference_r,float )
  DEFINE_PROPERTY( ConductivitySignal, conductivitySignal, double )
  DEFINE_PROPERTY( MetabolicSignal, metabolicSignal, double )
  
};



class VesselNode: public VNodeData, public ListTreeNode<VesselNode,Vessel,MAX_NUM_ADJACENT_EDGES>, public FixedSizeAllocated<VesselNode>
{
public:
  Int3 lpos;  
  VesselNode() : VNodeData() {}
private:
  friend class boost::serialization::access;
  template <class Archive>
    void serialize(Archive &ar, const unsigned int)
    {
	ar & boost::serialization::base_object<VNodeData>(*this);
	ar & lpos;
    }
};

class Vessel: public VData, public ListTreeEdge<VesselNode,Vessel>, public FixedSizeAllocated<Vessel>
{
public:   
  word len;
  word dir;
public:
  Vessel() : len(-1),dir(-1),VData() {}
  float WorldLength( const polymorphic_latticedata::LatticeData &ld ) const;
  const Int3 LPosA() const { return NodeA()->lpos; }
  const Int3 LPosB() const { return NodeB()->lpos; }
  float getWorldLength(){return sqrt(norm(NodeA()->worldpos - NodeB()->worldpos));}
  template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
      ar & len;
      ar & dir;
    }
};


struct HemodynamicBounds
{
  my::MinMax<double> flow;
  my::MinMax<double> force;
  my::MinMax<double> velocity;
  my::MinMax<double> press;
  my::MinMax<double> condSignal;
  HemodynamicBounds() {}
  HemodynamicBounds( const VesselList3d *g ) { Add(g,false); }
  void Add( const VesselList3d *g, bool bClear=true );
};

void ComputeCirculatedComponents( VesselList3d *list );


float GetInitialThickness( float rad );



#include "lookup3dquad.h"

#if 0
// idea how to store application dependent information
struct UserDataArrays
{
  typedef boost::variant<
    DynArray<int>,
    DynArray<float>,
    DynArray<double>,
    DynArray<uchar>
  > ArrayType;
  std::hash_map<string, ArrayType> data;

  template<class T>
  void add_data(const string &name);

  template<class T>
  DynArray<T>& get_data(const string &name);

  void push_back();
  void pop_back();
  void swap(int i, int j);
  void copy(int i, int j);
};
#endif




class VesselList3d
{
  typedef polymorphic_latticedata::LatticeData LD;
public:
  VesselList3d(boost::shared_ptr<LD> this_ld);
  typedef LD LatticeData; // nicer name for the outside
  typedef LD::SiteType SiteType;
  typedef ListGraph<VesselNode,Vessel> Graphtype;
  typedef boost::unordered_map <const VesselNode*, FlowBC> BCList; // boundary condition list

  void Flush();
  void Reset();
  void Init( const LD &ld );
  inline const LD& Ld() const { return *m_ld; }
  inline bool HasLattice() const { return m_ld.get() != nullptr; }
  Graphtype& Graph() { return g; } 
  const Graphtype& Graph() const { return g; }
  void SetDomainOrigin(const Float3 &pos) { m_ld->SetOriginPosition(pos); }

  inline Vessel*           GetEdge( size_t i )       { return g.edge(i); }
  inline const Vessel*     GetEdge( size_t i ) const { return g.edge(i); }
  inline VesselNode*       GetNode( size_t i )       { return g.node(i); }
  inline const VesselNode* GetNode( size_t i ) const { return g.node(i); }
  int GetECount() const { return g.num_edges(); }
  int GetNCount() const { return g.num_nodes(); }

  VesselNode* 	InsertNode( const Int3 &a );
  VesselNode* 	InsertNode( const Float3 &a );
  VesselNode*	InsertNode( const VesselNode *p_n);
  Vessel*     	InsertVessel( const Int3 &a, const Int3 &b );
  //this was created with world stuff, maybe I need to redo that
  //most MW stuff is not bad
  Vessel*     	InsertVessel( VesselNode* a, VesselNode* b );
  Vessel*     	InsertVessel( const Float3 &a, const Float3 &b );
  void 		InsertVessel( const Vessel *p_v);
  
  void        DeleteVessel( Vessel* v, bool bDeleteNodes = true );
  void        DeleteVesselWorld( Vessel* v, bool bDeleteNodes = true );
  FloatBBox3  GetWorldBoxFromVesselsOnly();
  void        DeleteUnusedNode( VesselNode* vc );
  void        DeleteUnusedNodeWorld( VesselNode* vc);
  void        SplitVessel( Vessel* v, int pos, Vessel* &vnew, VesselNode* &vcnew );
  Vessel*     FindVessel( const Int3 &a, const Int3 &b );
  VesselNode* FindNode( const Int3 &a );
  bool        FindAny( const Int3 &a, Vessel* &v, VesselNode* &vc );
  Vessel*     InsertTempEdge( VesselNode* a, VesselNode* b ) { return g.InsertEdge(a,b); }
  void        DeleteTempEdge( Vessel* v ) { g.DeleteEdge(v,false); }
  void        OptimizeMemLayout() { g.Optimize(); }
  //inline const LD& Ld() const { return *m_ld; }
  inline const BCList& GetBCMap() const {return bclist;}
  //const BCList& GetBCMap() const { return bclist; }
  void          SetBC(const VesselNode* node, FlowBC bc);
  void          ClearBC(VesselNode* node);

  std::size_t estimateMemoryUsage() const;
  void        IntegrityCheck(int check_lookup = -1);
  
  std::auto_ptr<VesselList3d> Clone();
  boost::shared_ptr<LD> getLD() const;
  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive &ar, const unsigned int version);
    SiteLookup 			lookup_site;
    BondLookup 			lookup_bond;
    boost::shared_ptr<LD>	m_ld;
    ListGraph<VesselNode,Vessel>g;
    void FillLookup();
    BCList bclist; // boundary conditions
    void DeleteUnusedNode(VesselNode* vc, int site);
};


inline std::size_t estimateMemoryUsage(const VesselList3d &vl) { return vl.estimateMemoryUsage(); }


/*------------------------------------------------------
------------------------------------------------------*/

uint Optimize( VesselList3d *vl );
std::auto_ptr<VesselList3d> GetSubdivided(std::auto_ptr<VesselList3d> vl, int multi, float newscale, int safety_boundary = 1);
std::auto_ptr<VesselList3d> GetSubdivided(std::auto_ptr<VesselList3d> vl, float scale);
/* Make it so that one vessels covers one and only one lattice bonds. 
 * Removes any vessel which is longer than one bond and replaces it with
 * several shorter vessels. VData is copied from the original to the smaller ones.
 * New nodes are created in the process. Data is copied from one of the two orignal
 * nodes. From which is not specified and can vary. The flags member is an exception
 * as flags are and'ed between the two originals. The boundary flags is always removed. 
 * This is to preserve arterio-venous flags.
 */
void SplitSegmentsToOneLatticeBond(VesselList3d &vl);
void TopoSortVessels( VesselList3d &vl, DynArray<VesselNode*> &sorted );
void TopoSortVessels(const VesselList3d &vl, DynArray<int> &order );

inline const VesselNode* GetUpstreamNode(const Vessel* v)
{
  const VesselNode *a = v->NodeA(),
                   *b = v->NodeB();
  return (a->press > b->press) ? a : b;
}
inline const VesselNode* GetDownstreamNode(const Vessel* v)
{
  const VesselNode *a = v->NodeA(),
		   *b = v->NodeB();
  return (a->press >b->press) ? b : a;
}
inline VesselNode* GetDownstreamNode( Vessel* v)
{
  VesselNode *a = v->NodeA(),
		   *b = v->NodeB();
  return (a->press >b->press) ? b : a;
}

void CheckToposort(const VesselList3d &vl, const DynArray<int> &order);

template <class Archive>
void VesselList3d::serialize(Archive& ar, const unsigned int version)
{
    ar & lookup_site;
    ar & lookup_bond;
    ar & m_ld;
    ar & g;
    ar & bclist;
}

/* note:
 * it is important to place the overwrites of 
 * save_construct_data and load_construct_data
 * here in the cpp. If placing direcly in the 
 * header the "friendship" between boost::serialization
 * and the VesselList3d is not well established 
 *
 * hm this is anyway not the case,
 */
namespace boost{namespace serialization{
template <class Archive>
inline void save_construct_data(
  Archive &ar, const VesselList3d *t, const unsigned int file_version
	    )
{
  // save data required to construct instance
  boost::shared_ptr<VesselList3d::LatticeData> an_other_ld = t->getLD();
  ar & an_other_ld;
}
template <class Archive>
inline void load_construct_data(
  Archive &ar, VesselList3d *t, const unsigned int file_version
	    )
{
  // retrieve data from archive required to construct new instance
  boost::shared_ptr<VesselList3d::LatticeData> m_ld;
  ar & m_ld;
  // invoke inplace constructor to initialize instance of my_class
  ::new(t)VesselList3d(m_ld);
  //t->Init(*m_ld);
}
}}//namespace boost{namespace serialization{

#endif //#define VESSELS3D_H_
