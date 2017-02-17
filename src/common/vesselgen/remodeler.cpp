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
#include "mwlib/timer.h"
#include "mwlib/vector_ext.h"
#include "../continuum-utils.h"

#include <boost/assign/std/vector.hpp>
#include <boost/foreach.hpp>
#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>
#include <boost/foreach.hpp>


using namespace boost::property_tree;




float getF(const VesselNode* vc)
{
  float f = 0;
  for (int j=0; j<vc->Count(); ++j)
  {
    f += vc->GetEdge(j)->f;
  }
  f /= (vc->Count() + 1.e-13);
  return f;
}


float getR(const VesselNode* vc)
{
  float f = 0;
  for (int j=0; j<vc->Count(); ++j)
  {
    f += vc->GetEdge(j)->r;
  }
  f /= (vc->Count() + 1.e-13);
  return f;
}


int nonCapillaryCount(const VesselNode* vc)
{
  int num_neighbors = 0;
  for (int j=0; j<vc->Count(); ++j)
  {
    if (vc->GetEdge(j)->IsCapillary()) continue;
    ++num_neighbors;
  }
  return num_neighbors;
}


bool isTip(const VesselNode* vc)
{
  return nonCapillaryCount(vc) == 1;
}


bool isCirculated(const VesselNode* vc)
{
  for (int j=0; j<vc->Count(); ++j)
  {
    if (vc->GetEdge(j)->IsCirculated()) return true;
  }
  return false;
}


bool hasTooMuchNeighborsForAttachment(const VesselNode* vc)
{
  return vc->Count() >= (vc->IsBoundary() ? 1 : 3);
}



enum {
  VESSGEN_OPT_MODE_V1,
  VESSGEN_OPT_MODE_V2,
  VESSGEN_OPT_MODE_V3
};


void Grower::Init(const ptree &settings)
{
  Int3 size = settings.get<Int3>("lattice_size");
  double scale = settings.get<double>("lattice_spacing");
  uint seed = settings.get<uint>("seed");
  string ld_type = settings.get<string>("lattice_type");
  if (ld_type == "fcc")
    lattice_type_id = LATTICE_TYPE_FCC;
  else if (ld_type == "quad")
    lattice_type_id = LATTICE_TYPE_QUAD;
  else
    throw std::invalid_argument(str(format("invalid lattice type string %s") % ld_type));
  gf_range = settings.get<double>("o2.diffusion_range");
  max_sprout_radius_artery = settings.get<double>("max_sprout_radius_artery");
  max_sprout_radius_vein = settings.get<double>("max_sprout_radius_vein");
  radius_vein = settings.get<double>("radius_vein");
  radius_capi = settings.get<double>("radius_capi");
  radius_artery = settings.get<double>("radius_artery");
  murray_alpha_vein = settings.get<double>("murray_alpha_vein");
  murray_alpha_artery = settings.get<double>("murray_alpha_artery");
  debug_output_every_configuration = settings.get<bool>("full_debug_output",false);
  
  iteration_number = 0;
  hierarchy_level = 0;
  iteration_number_on_level = 0;
  max_num_iter = settings.get<int>("max_num_iter");
  max_hierarchy_level = settings.get<int>("num_hierarchical_iterations");
  generate_more_capillaries = settings.get<bool>("generate_more_capillaries", false);
  capillariesUntilLevel = settings.get<int>("capillariesUntilLevel");
  const auto bloodFlowParametersPtree = settings.get_child_optional("calcflow");
  if (bloodFlowParametersPtree)
    bloodFlowParameters.assign(*bloodFlowParametersPtree);
  
  if (seed == 0)
    seed = 7731 * my::Time().to_ms();
  rand.Init(seed);

  std::auto_ptr<LatticeData> ldp = LatticeData::Make(ld_type.c_str(), BBox3().Add(Int3(0)).Add(size-Int3(1)), scale);

  vl.reset( new VesselList3d() );
  vl->Init(*ldp);

  cout << "vesselgen init ..." << endl;
  //cout << "  size " << size << " scale " << scale << endl;
  ldp->print(cout); cout << endl;
  cout << "  seed " << seed << endl;

#if GFFIELD_ENABLE
  InitGfDistrib();
#endif

  GetWorldDirectionAxes(vl->Ld(), wdiraxes, GET_WORLD_DIRECTION_AXES_NORMALIZE);
  for (int dir=0; dir<get_ld().NbCount(); ++dir)
  {
    dir_is_forbidden[dir] = (abs(wdiraxes[dir][2]) > 1.e-12) && (dim()==2);
  }
  PrecomputeElements();
}


#if GFFIELD_ENABLE
void Grower::InitGfDistrib()
{
  const LatticeData &ld = get_ld();
  
  const float cellsize = 50.;
  Float3 domain_size = Size(ld.GetWorldBox());
  Float3 fncells = domain_size / cellsize; //(0.5 * ld.Scale());
  BBox3 bb;
  for (int i=0; i<3; ++i)
  {
    bb.min[i] = -1;
    bb.max[i] = int(fncells[i]+0.99);
  }
  if (dim() == 2)
  {
    bb.min[2] = bb.max[2] = 0;
  }

  LatticeDataQuad3d field_ld;
  field_ld.Init(bb, cellsize);
  field_ld.SetCellCentering(Bool3(true, true, dim()==3));
  field_ld.SetOriginPosition(ld.GetWorldBox().min);
  grid.init(field_ld, dim());
  mtboxes.insert(0, field_ld.Box());
  cout << "field lattice:  " << field_ld << endl;
  ptree pt = make_ptree("rGf", gf_range);
  gfmodel.init(grid, mtboxes, pt);
  gfmodel.initField(fieldgf);
}





Array3d<float> Grower::ComputeGfSources() const
{
  Array3d<float> gfsource2(get_field_ld().Box());
  CylinderNetworkSampler sampler; sampler.Init(get_field_ld().Scale(), make_ptree("samples_per_cell", 3));
  for(int i=0; i<vl->GetECount(); ++i)
  {
    const Vessel* v = vl->GetEdge(i);
    if(!v->IsCirculated()) continue;
    sampler.Set(get_ld().LatticeToWorld(v->LPosA()), get_ld().LatticeToWorld(v->LPosB()), v->r);
    int cnt = sampler.GenerateLineSamples();
    for(int k=0; k<cnt; ++k)
    {
      AddSmoothDelta<float>(gfsource2, get_field_ld().Box(), get_field_ld(), dim(), sampler.GetSample(k).wpos, 1.);
    }
  }
  return gfsource2;
}


void Grower::UpdateGf()
{
  Array3d<float> gfsource = ComputeGfSources();
  gfmodel.update(fieldgf, gfsource);
}


inline Float3 Grower::GetGfGradient(const Float3& pos) const
{
  return FieldInterpolate::Gradient(fieldgf, get_field_ld(), FieldInterpolate::Extrapolate(), pos);
}

inline float Grower::GetGf( const Float3 &pos ) const
{
  return FieldInterpolate::ValueAveraged(fieldgf, get_field_ld(), FieldInterpolate::Extrapolate(), pos);
}

#endif

// called by Init
void Grower::PrecomputeElements()
{
  const LatticeData &ld = get_ld();
  
  PrecElem elem;
  for (int dir0=0; dir0<ld.NbCount(); ++dir0)
  {
    if (dir_is_forbidden[dir0]) continue; // in 2d skip configurations which are out of the z=0 plane

    elem.dirs[0] = dir0;
    elem.t = ELEM_I;
    precomputed_elements[std::pair<int,int>(ELEM_I,dir0)].push_back(elem);
    //elem.remove_all();

    // ELEM_Y and _T
    for (int dir10=0; dir10<ld.NbCount(); ++dir10)
    {
      for (int dir11=0; dir11<ld.NbCount(); ++dir11)
      {
        if (dir11 <= dir10) continue; // remove permuations where dir10 and dir11 are swapped
        Float3 v0 = wdiraxes[dir0],
               v10 = wdiraxes[dir10],
               v11 = wdiraxes[dir11];
        float det = v11.cross(v10).dot(v0);
        if (std::abs(det) > 1.e-3) continue;  // non-planar configuration
        if (v0.dot(v10) < -1.e-3) continue; // backward dir
        if (v0.dot(v11) < -1.e-3) continue; // backward dir
        if (std::abs(v0.dot(v10)-v0.dot(v11)) > 1.e-3) continue; // v11 and v10 must have same angle to v0

        if (dir_is_forbidden[dir0]) continue; // in 2d skip configurations which are out of the z=0 plane
        if (dir_is_forbidden[dir0]) continue; // in 2d skip configurations which are out of the z=0 plane
             
        elem.dirs[1] = dir10;
        elem.dirs[2] = dir11;
        elem.t = ELEM_Y;
        precomputed_elements[std::pair<int,int>(ELEM_Y,dir0)].push_back(elem);

        // ELEM_T
        if (lattice_type_id == LATTICE_TYPE_QUAD) // special element for quad lattice
        {
          myAssert(false); // this piece of code seems to be wrong obviously
          elem.t = ELEM_T;
          precomputed_elements[std::pair<int,int>(ELEM_T,dir0)].push_back(elem);
        }
      }
    }
  }

#if 0
  {
    typedef std::pair<int,int> K;
    typedef boost::unordered_map<K, DynArray<PrecElem> > Map;
    for (Map::const_iterator it = precomputed_elements.begin(); it != precomputed_elements.end(); ++it)
    {
      int el_type, dir;
      boost::tie(el_type, dir) = it->first;
      cout << format("el_type=%i, dir=%i -> %i variants =") % el_type % dir % it->second.size() << endl;
      const DynArray<PrecElem> &elems = it->second;
      for (int i=0; i<elems.size(); ++i)
      {
        cout << format("el%02i: ") % i;
        for (int j=0; j<3; ++j)
        {
          cout << format("%i,") % elems[i].dirs[j];
        }
        cout << endl;
      }
    }
  }
#endif
}

/*
  |
  
  -+-
   |

 |   |
 +-+-+
   |
*/

int Grower::GetNumTrialElements(int dir, ElemType type) const
{
  return const_cast<Grower*>(this)->precomputed_elements[std::pair<int,int>(type,dir)].size();
}


void Grower::MakeElement( const Int3 &pos, int dir, int index, ElemType type, Element &edges ) const
{
  const PrecElem &pe = const_cast<Grower*>(this)->precomputed_elements[std::pair<int,int>(type,dir)][index];
  // this lookup will crash in 2d if one tries to access a 3d configuration
  
  const LatticeData &ld = get_ld();
  edges.remove_all();
  edges.reserve(5);
  Int3 p0 = pos;
  Int3 p1 = ld.NbLattice( p0, pe.dirs[0] );
  edges.push_back( ElemEdge( p0, p1, pe.dirs[0], EL_START ) );
  if( type == ELEM_I )
  {
    edges[0].flags.AddBits(EL_END);
    return;
  }
  int dr = pe.dirs[1];
  Int3 p1r = ld.NbLattice( p1, dr );
  int dl = pe.dirs[2];
  Int3 p1l = ld.NbLattice( p1, dl );
  edges.push_back( ElemEdge( p1, p1r, dr, type==ELEM_Y ? EL_END : 0) );
  edges.push_back( ElemEdge( p1, p1l, dl, type==ELEM_Y ? EL_END : 0) );
  if( type == ELEM_T )
  {
    edges.push_back( ElemEdge( p1r, ld.NbLattice(p1r,dir), dir, EL_END ) );
    edges.push_back( ElemEdge( p1l, ld.NbLattice(p1l,dir), dir, EL_END ) );
  }
}


bool Grower::TestFitElement( const Element &elem, uchar avflags)
{
  const LatticeData &ld = get_ld();
  for( int i=0; i<elem.size(); ++i )
  {
    const ElemEdge& e = elem[i];
    if( !ld.IsInsideLattice(e.a) ) return false;
    if( !ld.IsInsideLattice(e.b) ) return false;
    VesselNode* vca = vl->FindNode(e.a);
    if( !e.flags.GetBits(EL_START) && vca ) return false;
    VesselNode* vc = NULL; Vessel *v = NULL;
    vl->FindAny(e.b, v, vc); // b is the endpoint of the segment, as seen from the point where the element is attached
    if (v || vc) return false;
    // now check if there is one site of space from here to some vessels of the same type.
    // Check the endpoints only, should be good.
#if 0
    const Int3 start_point = elem[0].a;
    if (dim() < 3 && !(i == 0 &&  vca && vca->Count()>1))
    { // TODO: perhaps do this only in 2d
      for (int dir=0; dir<ld.NbCount(); ++dir)
      {
        if (dir_is_forbidden[dir]) continue; // skip out-of-plane directions in 2d
        Int3 neighbor_site = ld.NbLattice(e.b, dir);
        if (i == 0 && neighbor_site == start_point) continue; // the end point of the first segment has the start point of the element as neighbor. Thats okay.
        if (!ld.IsInsideLattice(neighbor_site)) continue; // is okay, we don't put anything there
        VesselNode* neighbor = vl->FindNode(neighbor_site);
        // check if there is a neighbor which has the same type of arterial or venous. If so, then we cannot put stuff here, to leave some room for growth
        if (neighbor && (neighbor->flags & avflags)!=0)
          return false;
      }
    }
#endif
  }
  return true;
}

void Grower::AddElement( const Element &elem, uchar flags, OpenEnds *ends )
{
  flags = flags & (ARTERY | VEIN | CAPILLARY);
  for( int i=0; i<elem.size(); ++i )
  {
    const ElemEdge& e = elem[i];
    Vessel* vv = vl->InsertVessel( e.a, e.b );
    vv->flags.AddBits(flags);
    vv->NodeA()->flags.AddBits(flags);
    vv->NodeB()->flags.AddBits(flags);
    if (flags & ARTERY)
      vv->r = radius_artery;
    else if (flags & VEIN)
      vv->r = radius_vein;
    else
      vv->r  = radius_capi;
    if( e.flags.GetBits(EL_END) && ends )
    {
      ends->add(e.b,e.dir,flags);
    }
  }
}

const OpenEnd Grower::PickEnd(OpenEnds &ends)
{
  int idx = rand.Get(ends.size());
  return ends.pop(idx);
}


int Grower::GetDirection(const VesselNode* vc)
{
  if (vc->IsBoundary())
  {
    auto it = tree_roots.find(vc->lpos);
    myAssert(it != tree_roots.end()); // if vc->Boundary says it is a root node, there should be an entry in the list of root sites
    return it->second.dir;
  }
  else
  {
    if (vc->Count() == 0)
      return -1;
    else if (vc->Count() == 1)
      return vc->GetEdge(0)->dir;
    else if (vc->Count() == 2)
    {
      // return the direction of the up-stream vessel
      if (vc->GetNode(0).node->press > vc->GetNode(1).node->press)
        return vc->GetEdge(0)->dir;
      else
        return vc->GetEdge(1)->dir;
    }
    else
      return -1;
  }
}

/*---------------------------------------------------------
  --------------------------------------------------------*/

void Grower::CalcFlow()
{
  ::CalcFlow(*vl, bloodFlowParameters);
  hemodynBounds.Add(vl.get());
}

/*---------------------------------------------------------
  --------------------------------------------------------*/

static int GetTreeDownVessels( VesselNode* vc, Vessel* parent, Vessel* *down )
{
  int cnt=0;
  for( int i=0; i<vc->Count(); ++i ) {
    Vessel* v = vc->GetEdge(i);
    if( v->flags.GetBits(CAPILLARY) || v==parent ) continue;
    //if (!v->IsCirculated()) continue;
    down[cnt++] = v;
  }
  return cnt;
}


void Grower::CalcRadi()
{
  int backup_capillariesUntilLevel = capillariesUntilLevel;
  capillariesUntilLevel = my::round(my::lerp(float(hierarchy_level)/max_hierarchy_level, 0, backup_capillariesUntilLevel));
  
  Vessel* rootvessels[32];
  for( int i=0; i<vl->GetNCount(); ++i )
  {
    VesselNode* root = vl->GetNode(i);
    if( !root->flags.GetBits( BOUNDARY ) ) continue; // start from boundary nodes and work the way down/up-stream
    myAssert(root->Count() <= 1);
    int n = GetTreeDownVessels(root, NULL, rootvessels);
    double radius = 0.;
    for( int k=0; k<n; ++k )
    {     
      CalcRadi(rootvessels[k], root);
      radius = std::max<double>(rootvessels[k]->r, radius);
      //press = std::max<double>(press, PressureRadiusRelation( rootvessels[k]->r, rootvessels[k]->IsArtery() ));
      //rootvessels[k]->flags |= BOUNDARY;
    }
    //root->press = press;
    root->press = PressureRadiusRelation(radius, root->flags&ARTERY);
  }

  capillariesUntilLevel = backup_capillariesUntilLevel;
}


struct TreeStackData {
  TreeStackData( Vessel* v, VesselNode* vcparent, int parent_index) : v(v),vcparent(vcparent), level(0), parent_index(parent_index) {}
  Vessel* v;
  VesselNode* vcparent;
  int level;
  int parent_index;
};

void Grower::CalcRadi( Vessel* vroot, VesselNode* nroot  )
{
  const double alphaVein = murray_alpha_vein;
  const double alpha     = murray_alpha_artery;
	const double rArtery   = radius_artery;
	const double rVein     = radius_vein;
  
  DynArray<TreeStackData> stack;
  DynArray<uchar> visited;
  stack.reserve( vl->GetECount() );
  visited.resize( vl->GetECount(), 0 );
  stack.push_back( TreeStackData(vroot,nroot, -1) );
  while( !stack.empty() )
  {
    int index            = stack.size()-1;
    TreeStackData &sd    = stack[index];
    Vessel*            v = sd.v;
    VesselNode* vcparent = sd.vcparent;
    VesselNode* vcbranch = v->GetOther(vcparent);
    myAssert(!v->IsCapillary());
    Vessel* d[32]; d[0]=NULL; d[1]=NULL; // init only the first two elements, the rest of the reserved memory is there for safety
    int numChildren = GetTreeDownVessels( vcbranch, v, d );
    myAssert(numChildren <= 2);
    if( visited[v->Index()]==0 ) // unvisited
    {
      if( d[0]==NULL && d[1]==NULL )
      {
        if (capillariesUntilLevel > 0) // this is level 0, so ...
          v->r = radius_capi;
        else
          v->r = (v->flags&ARTERY) ? rArtery : rVein;
        stack.pop_back();
      }
      else // if this is not a leaf, add its children onto the stack
      {
        if( d[0] ) {
          stack.push_back( TreeStackData(d[0],vcbranch, index) );
        }
        if( d[1] ) {
          stack.push_back( TreeStackData(d[1],vcbranch, index) );
        }
      }
      visited[v->Index()] = 1;
    }
    else // this vessel has been visited before
    {
      // therefore we can compute its radius, knowing the radi of the children
      float r[2]={0,0};
      if( d[0] && d[1] )
      {
        r[0] = d[0]->r;
        r[1] = d[1]->r;
      }
      else if( d[1] )
      {
        r[0] = d[1]->r;
      }
      else if( d[0] )
      {
        r[0] = d[0]->r;
      }
      myAssert(d[0] == NULL || visited[d[0]->Index()]); // should have already visited children
      myAssert(d[1] == NULL || visited[d[1]->Index()]);
      
      const double a = (v->flags&VEIN) ? alphaVein : alpha;
      if (sd.level == capillariesUntilLevel)
      {
        v->r = (v->flags&ARTERY) ? rArtery : rVein;
      }
      else if (sd.level < capillariesUntilLevel)
      {
        v->r = radius_capi;
      }
      /*
      //modiefied murray law (see Safaeian thesis) // added by thierry
      // should check if this leads to measurable differences in rBF, rBV
      else if (sd.level >2)
      {
        v->r = powf(powf(r[0],3.6) + powf(r[1],3.6), 1.0f/3.6 );
      }
      */
      else
      {
        v->r = powf(powf(r[0],a) + powf(r[1],a), 1.0f/a );
      }

      if (sd.parent_index >= 0)
      {
        TreeStackData &parent_sd = stack[sd.parent_index];
        parent_sd.level = std::max(sd.level+1, parent_sd.level);
      }
      myAssert(sd.level>=capillariesUntilLevel || (v->r <= rArtery || v->r <= rVein || v->r <= radius_capi));
      
      v->timeSprout = sd.level; // for debugging
      stack.pop_back();
    }
  }
}

/*---------------------------------------------------------
  --------------------------------------------------------*/

#if 1
struct CapillaryGenerator
{
  VesselList3d &vl;
  const VesselList3d::LatticeData &ld;
  boost::unordered_set<int> visited;
  DynArray<DynArray<Int3> > paths;
  uint source_type, opposing_type;
  float max_radius_artery, max_radius_vein;

  CapillaryGenerator(VesselList3d &vl, float max_radius_artery, float max_radius_vein) :
    vl(vl), ld(vl.Ld()), max_radius_artery(max_radius_artery), max_radius_vein(max_radius_vein)
  {
  }

  bool canConnect(const VesselNode* nd)
  {
    if (hasTooMuchNeighborsForAttachment(nd)) return false;
    if (nd->Count() <= 0) return false;
    if (!isTip(nd)) return false;
    //if (nd->flags.GetBits(ARTERY) && nd->GetEdge(0)->r>max_radius_artery) return false;
    //if (nd->flags.GetBits(VEIN) && nd->GetEdge(0)->r>max_radius_vein) return false;
    return true;
  }

  void findTargets(const Int3 &pos, int neighborhood, DynArray<Int3> &parent_path)
  {
    if (neighborhood <= 0) return;
    for( int dir=0; dir<ld.NbCount(); ++dir )
    {
      const Int3 nbpos = ld.NbLattice(pos,dir);
      if( !ld.IsInsideLattice(nbpos) ) continue;
      int nbsite = ld.LatticeToSite(nbpos);
      if (visited.find(nbsite) != visited.end()) continue;
      visited.insert(nbsite);
      if( vl.FindVessel( pos,nbpos ) ) continue;
      //const VesselNode* nbnode =
      VesselNode* nbnode=NULL;
      Vessel* nbvessel=NULL;
      vl.FindAny(nbpos, nbvessel, nbnode);
      if (nbnode)
      {
        if (!canConnect(nbnode)) continue;
        if (nbnode->flags.GetBits(ARTERY | VEIN) != opposing_type) continue;
        paths.push_back(parent_path);
        paths.back().push_back(nbpos);
      }
      else if (nbvessel)
      {
        continue;
      }
      else if (neighborhood > 1)
      {
        DynArray<Int3> path = parent_path;
        path.push_back(nbpos);
        findTargets(nbpos, neighborhood-1, path);
      }
    }
  }

  void findTargets(const VesselNode* vc, int neighborhood)
  {
    if (!canConnect(vc)) return;
    source_type = vc->flags.GetBits(ARTERY | VEIN);
    myAssert(source_type == ARTERY || source_type == VEIN);
    opposing_type = (source_type==ARTERY) ? VEIN : ARTERY;
    DynArray<Int3> path;
    findTargets(vc->lpos,neighborhood, path);
  }
};
#endif


void Grower::AddCapillaries()
{
#if 0
  DynArray<const VesselNode*> nodelist(vl->GetNCount());
  for( int i=0; i<vl->GetNCount(); ++i )
    nodelist[i] = vl->GetNode(i);

  RandomPermute(rand, get_ptr(nodelist), nodelist.size());
  
  const LatticeData &ld = get_ld();

  for (int i=0; i<nodelist.size(); ++i)
  {
    const VesselNode* vc = nodelist[i];
    CapillaryGenerator capgen(*vl, max_sprout_radius_artery, max_sprout_radius_vein);
    capgen.findTargets(vc, 2);
    if (capgen.paths.size() <= 0) continue;
    int path_index = rand.Get(capgen.paths.size());
    const DynArray<Int3> &path = capgen.paths[path_index];
    Int3 pos = vc->lpos;
    for (int j=0; j<path.size(); ++j)
    {
      Vessel* v = vl->InsertVessel(pos, path[j]);
      pos = path[j];
      v->flags.AddBits( CAPILLARY );
      v->r = radius_capi;
    }
  }
#else
  DynArray<const VesselNode*> nodelist(vl->GetNCount());
  for( int i=0; i<vl->GetNCount(); ++i )
    nodelist[i] = vl->GetNode(i);
  while (nodelist.size() > 0)
  {
    const VesselNode* vc = GetSwapRemoved(nodelist, rand.Get(nodelist.size()));
    CapillaryGenerator capgen(*vl, max_sprout_radius_artery, max_sprout_radius_vein);
    capgen.findTargets(vc, generate_more_capillaries ? 6 : 1);
    if (capgen.paths.size() <= 0) continue;
    const DynArray<Int3> &path = capgen.paths[rand.Get(capgen.paths.size())];
    Int3 pos = vc->lpos;
    for (int j=0; j<path.size(); ++j)
    {
      Vessel* v = vl->InsertVessel(pos, path[j]);
      pos = path[j];
      v->flags.AddBits( CAPILLARY | vc->flags.GetBits(ARTERY | VEIN));
      v->r = radius_capi;
      if (j < path.size()-1 && generate_more_capillaries) 
      {
        // Add each new capillary nodes to the list of node from which to start capillaries.
        // This should generate a denser capillary plexus.
        VesselNode* vc2 = vl->FindNode(path[j]);
        vc2->flags.AddBits(CAPILLARY | vc->flags.GetBits(ARTERY | VEIN));
        nodelist.push_back(vc2);
      }
    }
  }
#endif
}

void Grower::FlushCapillaries()
{
  int i = 0;
  while( i < vl->GetECount() )
  {
    Vessel* v = vl->GetEdge(i);
    if( v->flags.GetBits( CAPILLARY ) ) vl->DeleteVessel(v);
    else ++i;
  }
}

/*---------------------------------------------------------
  --------------------------------------------------------*/

bool Grower::AddBeginning(OpenEnds &ends, Int3 pos, int dir, uchar flags, int len )
{
    if( !get_ld().IsInsideLattice(pos) ) {
      throw std::invalid_argument(str(format("Error: Point %s is not in lattice!") % pos));
    }
    if(len<=0) ends.add(pos,dir,flags);
    if (vl->FindNode(pos))
    {
      // silently reject
      return false;
    }
    VesselNode* initnode = vl->InsertNode(pos); 
    initnode->flags = flags|BOUNDARY;
    initnode->press = (flags & ARTERY) ? 1.0f : 0.0f;
    while(len-->0)
    {      
      Element elem;
      MakeElement(pos, dir, 0, ELEM_I, elem );
      if(!TestFitElement( elem, flags & (ARTERY | VEIN) ))
      {
        // silently reject
         //cout << (str(format("Warning: Collision in starting vessels at pos %s") % pos)) << endl;
         return len>=1;
      }
      AddElement( elem, flags, &ends );
      pos = get_ld().NbLattice(pos,dir);
    }
    return true;
}


#define THROW_LINE_ERROR throw std::runtime_error(str(format("bad root position %s in file") % pos))


void Grower::GenerateRootElements(const ptree& pt, OpenEnds &ends)
{
  const LatticeData &ld = get_ld();
  FloatBBox3 wbb = ld.GetWorldBox();
  Float3 wcenter = 0.5*(wbb.max+wbb.min);
  const BBox3 bb = ld.Box();
  
  BOOST_FOREACH(const ptree::value_type &vv, pt.get_child("roots"))
  {
    const ptree &v = vv.second;

    Int3 pos = v.get<Int3>("p");
    char type = v.get<char>("t");
    int len = v.get<int>("l", 0);
    int flags = type=='a' ? ARTERY : VEIN;
    // compute direction of extension: find direction whose angle matches best to vector from center to pos
    int dir = -1;
    Float3 wpos = (ld.LatticeToWorld(pos) - wcenter).normalized();

    if (len > 1 && (bb.min[0]==pos[0] || bb.max[0]==pos[0])) // we restrict the extensions to vessels at the left/right boundary
    {
      dir = (bb.min[0]==pos[0]) ?
        ld.GetAxisDirLen(pos, pos+Int3(1,0,0)).dir :
        ld.GetAxisDirLen(pos, pos+Int3(-1,0,0)).dir;
    }
    else
    { // pick a direction taht points to the center
      myAssert(len <= 1);
      float best_cosa = 0;
      len = 0;
      for (int j=0; j<ld.NbCount(); ++j)
      {
        if (dir_is_forbidden[j]) continue; // in 2d skip configurations which are out of the z=0 plane
        float cosa = wdiraxes[j].dot(wpos);
        if (cosa < best_cosa) {
          best_cosa = cosa;
          dir = j;
        }
      }
    }
    if(dir<0)
    {
      THROW_LINE_ERROR;
    }
    if(!ld.IsInsideLattice(pos) || !(type=='a' || type=='v'))
    {
      THROW_LINE_ERROR;
    }
    bool ok = AddBeginning(ends, pos, dir, flags,len);
    if (ok) {
      myAssert(vl->FindNode(pos));
      tree_roots.insert(std::make_pair(pos, TreeRoot(pos, dir, flags, len)));
    }
  }
}


void Grower::AddRandomElement( const OpenEnd &end, OpenEnds *ends, bool useGfGradient, int typesAllowedFlags)
{
    ElemType et[3] = { ELEM_I, ELEM_Y, ELEM_T };
    double prob_et[3] = { 1., 2.,2. };//probability to occure for the elements declared above.
    //T elements are not allowed for the fcc type
    if (lattice_type_id == LATTICE_TYPE_FCC)
      prob_et[2] = 0.;
    if ((typesAllowedFlags & ELEM_I) == 0)
      prob_et[0] = 0.;
    if ((typesAllowedFlags & ELEM_Y) == 0)
      prob_et[1] = 0.;
    if ((typesAllowedFlags & ELEM_I) == 0)
      prob_et[2] = 0.;
    assert(typesAllowedFlags > 0 && (prob_et[0]>0 || prob_et[1]>0 || prob_et[2]>0));
    NormalizeSumNorm( prob_et, 3 );
    /* NOTE: The old code from 2009 first tried to attach a Y and only if this failed,
     * a single segment was added!
     */

    int ielem = RandomPick( rand, 3, prob_et );//pick random element
    
    const VesselNode* attachnode = vl->FindNode( end.p );
    myAssert( attachnode );
    const uchar avflags = attachnode->flags&(ARTERY|VEIN);//the node knows wether its a Artery or a vein

    DynArray<std::pair<int,int> > config;
    DynArray<double> probs;
    
    for (int dir=0; dir<get_ld().NbCount(); ++dir)//loop over all neighbours in the lattice
    {
      if (dir_is_forbidden[dir]) continue; // in 2d skip configurations which are out of the z=0 plane

      float prob = std::numeric_limits<float>::quiet_NaN();
      if (end.dir >= 0)
      {
        float cosangle = wdiraxes[dir].dot(wdiraxes[end.dir]);
        if (cosangle < 1.e-3) continue; // must be forward dir

        prob = cosangle + 1.;
#if GFFIELD_ENABLE
        if (useGfGradient)
        {
          Float3 gfgrad = GetGfGradient(get_ld().LatticeToWorld(end.p));
          float pg = 1 - gfgrad.dot(wdiraxes[dir]); // 0 for up hill, 2 for downhill
          prob += pg;
        }
#endif
      }
      else prob = 1.;
      
      int num_test_elems = GetNumTrialElements(dir, et[ielem]);
      for (int index=0; index<num_test_elems; ++index)
      {
        Element elem;
        MakeElement( end.p, dir, index, et[ielem], elem );
        if( TestFitElement(elem, avflags) )
        {
          config.push_back(std::pair<int,int>(dir, index));
          probs.push_back(prob); // pick chance depends on angle to forward dir
        }
      }
    }
    if(!config.empty())
    {
      NormalizeSumNorm(get_ptr(probs), probs.size());
      int dir, index;
      boost::tie(dir, index) = config[RandomPick(rand,probs.size(),get_ptr(probs))];
      Element elem;
      MakeElement(end.p, dir, index, et[ielem], elem);
      myAssert(TestFitElement(elem, avflags));
      AddElement(elem, avflags, ends);
    }
}


enum Remodel {
  DIE, GROW, IDLE
};


struct VesselTipStuff
{
  float f, r;
  int cnt;
  bool is_circulated;
};

inline VesselTipStuff getStuff(const VesselNode* vc)
{
  const int cnt = vc->Count();
  VesselTipStuff x;
  x.cnt = cnt;
  if (cnt > 0)
  {
    const Vessel *w,*v = vc->GetEdge(0);
    float f, r;
    if (cnt == 1)
    {
      f=v->f; //vc->residual;
      r=v->r;
    x.is_circulated = v->IsCirculated();
    }
    else // just consider the second edge, because no change can happen to node with 3 adjacent vessels or more
    {
      w = vc->GetEdge(1);
      f = (v->f+w->f)*0.5f;
      r = (v->r+w->r)*0.5f;
      x.is_circulated = v->IsCirculated() | w->IsCirculated();
    }
    x.f = f;
    x.r = r;
  }
  else
  {
    x.f = x.r = 0;
    x.is_circulated = false;
  }
  return x;
}



void Grower::DetermineRemodelingAction(DynArray<int> &action)
{
  float death_prob = 1.;
  int ncnt = vl->GetNCount();
  
  float median_tip_f = 0;
  my::Averaged<float> force_all_bounds;
  my::Averaged<float> force_tip_bounds;

  {
    DynArray<float> ref_f_arr(1024, ConsTags::RESERVE);
    for (int i=0; i<vl->GetECount(); ++i)
    {
      const Vessel* v = vl->GetEdge(i);
      if (!v->IsCirculated()) continue;
      force_all_bounds.Add(v->f);
      if (v->NodeA()->Count() != 1 && v->NodeB()->Count() != 1) continue;
      force_tip_bounds.Add(v->f);
      ref_f_arr.push_back(v->f);
    }
    std::sort(ref_f_arr.begin(), ref_f_arr.end());
    if(ref_f_arr.size()>0)
      median_tip_f = ref_f_arr[(ref_f_arr.size()-1)/2];
    //cout << format("%f %f %f %i") % ref_f_arr.front() % ref_f % ref_f_arr.back() % ref_f_arr.size() << endl;
  }

  //float log_force_min = std::log(force_tip_bounds.Min());
  //float log_force_max = std::log(force_tip_bounds.Max());
  float log_force_min = std::log(force_all_bounds.Min());
  float log_force_max = std::log(force_all_bounds.Max());

#if GFFIELD_ENABLE
  my::Averaged<double> gf_mm = fieldgf.valueStatistics();
#endif
  
  for (int i=0; i<ncnt; ++i)
  {
    VesselNode* vc = vl->GetNode(i);
    VesselTipStuff x = getStuff(vc);
    uchar a = IDLE;

    bool canGrow = vc->Count() <= (vc->IsBoundary() ? 0 : 2);
    bool canDie  = !vc->IsBoundary() && vc->Count()<=1;
    canGrow &= !(x.r>max_sprout_radius_artery && vc->flags.GetBits(ARTERY));
    canGrow &= !(x.r>max_sprout_radius_vein && vc->flags.GetBits(VEIN));
    
#if GFFIELD_ENABLE
    {

      float gf = GetGf(get_ld().LatticeToWorld(vc->lpos));
      gf /= gf_mm.Max()+1.e-12;

      float feps2 = 1.e-2;
      float multiplier = vc->flags.GetBits(ARTERY) ? 1 : 1;
      float xi = x.f*multiplier/force_all_bounds.Max();
      float f_signal = x.is_circulated ? (xi / (xi + feps2)) : 0.;
      
      if( !std::isfinite(f_signal) ) continue;
      
      double gain   = 2.;
      double offset = 0.1;
      double pGrow = std::pow((1.0 - offset) * f_signal+offset * 0.5, gain);
      double pDie  = std::pow((1.0 - offset) * (1.0-f_signal) +offset * 0.5, gain);
      double pIdle = 1.0 - pDie - pGrow;
      
      if(gf < 0.01)
      {
        pGrow = 0.5;
        pDie  = 0.5;
        pIdle = 0.; 
      }
      
      double pSum  = pGrow+pDie+pIdle;
      uchar  apick[3] = { DIE, GROW, IDLE };
      double ppick[3] = { pDie/pSum, pGrow/pSum, pIdle/pSum };
      a = apick[RandomPick( rand, 3, ppick )];
      if (!canDie && a==DIE) a=IDLE;
      if (!canGrow && a==GROW) a=IDLE;
      action[i] = a;
    }
#else
    {
      float feps2 = 1.e-2;
      float multiplier = vc->flags.GetBits(ARTERY) ? 1 : 1;
      float xi = x.f*multiplier/force_all_bounds.Max();
      float f_signal = x.is_circulated ? (xi / (xi + feps2)) : 0.;

      if( !std::isfinite(f_signal) ) continue;

      double gain   = 2.;
      double offset = 0.01;
      double pGrow = std::pow((1.0 - offset) * f_signal+offset, gain);
      double pDie  = std::pow((1.0 - offset) * (1.0-f_signal) +offset, gain);
      double pIdle = 1.0 - pDie - pGrow;

      if (!x.is_circulated)
      {
        pGrow = 0.4;
        pDie  = 0.6;
        pIdle = 0.;
      }

      double pSum  = pGrow+pDie+pIdle;
      uchar  apick[3] = { DIE, GROW, IDLE };
      double ppick[3] = { pDie/pSum, pGrow/pSum, pIdle/pSum };
      a = apick[RandomPick( rand, 3, ppick )];
    if (!canDie && a==DIE) a=IDLE;
    if (!canGrow && a==GROW) a=IDLE;
    action[i] = a;
  }
#endif
  }
}


void Grower::RemodelTrees()
{
#if GFFIELD_ENABLE
  UpdateGf();
#endif
  
  DynArray<float> nodalShearStress(vl->GetNCount());
  DynArray<float> nodalAttachedCapis(vl->GetNCount());
  
  {
  int ecnt = vl->GetECount();
  for (int i=0; i<ecnt; ++i)
  {
      Vessel* v = vl->GetEdge(i);
      if (!v->IsCapillary()) continue;
      
      nodalShearStress[v->NodeA()->Index()] += v->f;
      nodalShearStress[v->NodeB()->Index()] += v->f;
      nodalAttachedCapis[v->NodeA()->Index()] += 1.;
      nodalAttachedCapis[v->NodeB()->Index()] += 1.;
  }
  for (int i=0; i<vl->GetNCount(); ++i)
  {
    VesselNode* nd = vl->GetNode(i);
    if (nodalAttachedCapis[i] > 0)
    {
      // we just use the residual variable here since we don't need it otherwise
      nd->residual = nodalShearStress[i] / nodalAttachedCapis[i];
    }
    else
      nd->residual = 0.;
  }
  }
 
  FlushCapillaries();

  VESSGEN_MAXDBG(vl->IntegrityCheck();)
  
  int ncnt = vl->GetNCount();
#ifdef WRITE_REMODELING_ACTIONS
  last_remodeling_actions.clear();
#endif
  DynArray<int> action(ncnt, IDLE);
  DetermineRemodelingAction(action);
  DynArray<VesselNode*> toKill;
  OpenEnds toGrow;
  for( int i=0; i<ncnt; ++i )
  {
    VesselNode* vc = vl->GetNode(i);
    int a = action[i];
    if( a==GROW )
    {
      int dir = GetDirection(vc);
      toGrow.add(vc->lpos, dir, vc->flags & (ARTERY|VEIN|CIRCULATED));
    }
    else if( a==DIE )
    {
      toKill.push_back(vc);
    }
#ifdef WRITE_REMODELING_ACTIONS
    last_remodeling_actions[vc] = a;
#endif
  }

  for( int i=0; i<toKill.size(); ++i )
  {
    myAssert(toKill[i]->Count() == 0 || toKill[i]->Count() == 1);
    while (toKill[i]->Count() > 0)
      vl->DeleteVessel(toKill[i]->GetEdge(0), false);
    vl->DeleteUnusedNode(toKill[i]);
  }
  
  //RandomGrowth(toGrow, false);

// obtain new terminal branches for random growth  
  while( toGrow.size()>0 )
  {
    OpenEnd end = PickEnd( toGrow );
    AddRandomElement( end, NULL, true, (end.flags&CIRCULATED) ? ELEM_ANY : (ELEM_Y|ELEM_T));
  }

  VESSGEN_MAXDBG(vl->IntegrityCheck();)

  AddCapillaries();

  VESSGEN_MAXDBG(vl->IntegrityCheck();)

  ComputeCirculatedComponents(vl.get());
  CalcRadi();
  CalcFlow();
#if GFFIELD_ENABLE
  UpdateGf();
#endif
  iteration_number++;
}



/*---------------------------------------------------------
  --------------------------------------------------------*/

#if GFFIELD_ENABLE
void Grower::GetGfAtNodes(DynArray<float> &data) const
{
  data.resize(vl->GetNCount());
  for (int i=0; i<vl->GetNCount(); ++i)
  {
    const Int3 lpos = vl->GetNode(i)->lpos;
    data[i] = GetGf(get_ld().LatticeToWorld(lpos));
  }
}

#endif


void Grower::SanityCheck()
{
  for (int i=0; i<vl->GetECount(); ++i)
  {
    const Vessel* v = vl->GetEdge(i);
    myAssert(!v->IsCapillary() || v->r == 2.);
  }
}



// void Grower::UpdateTreeRoots() //TODO: check if this is actually needed
// {
//   DynArray<TreeRoot> new_roots(1024, ConsTags::RESERVE);
//   
//   for (int j=0; j<tree_roots.size(); ++j)
//   {
//     bool ok = false;
//     for (int i=0; i<vl->GetNCount(); ++i)
//     {
//       const VesselNode* nd = vl->GetNode(i);
//       if (nd->lpos == tree_roots[j].p && nd->IsBoundary()) { ok = true; break; }
//     }
//     if (ok)
//       new_roots.push_back(tree_roots[j]);
//   }
//   tree_roots = new_roots;
// }


void Grower::UpscaleTrees()
{
  // first get vessels on subdivided lattice
  std::auto_ptr<LatticeData> oldld(vl->Ld().Clone());
  //hope this wont end too bad
  //std::auto_ptr<LatticeData> oldld(vl->Ld()->get());
  std::auto_ptr<VesselList3d> vl_new = GetSubdivided(vl, 2, oldld->Scale(), 0);
  vl = vl_new;
  const LatticeData &ld = vl->Ld();

  VESSGEN_MAXDBG(vl->IntegrityCheck();)

  // extend roots to the lattice boundary
  int reverse_dir[32];
  GetReverseDir(ld, reverse_dir);

  TreeRootList new_roots; // with updated lattice index keys
  BOOST_FOREACH(auto it, tree_roots)
  {
    TreeRoot root = it.second;
    int d = reverse_dir[root.dir];
    myAssert(!dir_is_forbidden[d]);
    
    root.p = oldld->GetLatticeIndexOnRefinedGrid(root.p, 1); // the root position on refined lattice. This is the starting point to pick a different position closer the boundary of the lattice
    
    VesselNode* previousRoot = vl->FindNode(root.p);
    myAssert(previousRoot && (previousRoot->flags & BOUNDARY)); // there should be a node, because it is the original root site (just copied to the refined lattice)
    VesselNode* vcnew = nullptr;
    {
      VesselNode* vc = previousRoot;
      while (true)
      {
        Int3 q = ld.NbLattice(root.p, d); // move one bound
        if (!ld.IsInsideLattice(q)) break;
        if (vl->FindNode(q)) break;  
        // add a vascular segment
        vcnew = vl->InsertNode(q);
        Vessel* v = vl->InsertVessel(vc, vcnew);
        v->flags = vc->flags & (~BOUNDARY);
        vcnew->flags = vc->flags & (~BOUNDARY);
        root.p = q;
        vc = vcnew;
      }
    }
    if (vcnew) // did this actually do anything
    { 
      // set boundary condition to new node
      auto it = vl->GetBCMap().find(previousRoot);
      myAssert(it != vl->GetBCMap().end());

      vcnew->flags.AddBits(BOUNDARY);
      if (it != vl->GetBCMap().end())
      {
        vl->SetBC(vcnew, it->second);
        vl->ClearBC(previousRoot);
      }
      previousRoot->flags.DelBits(BOUNDARY);
    }
    
    new_roots.insert(std::make_pair(root.p, root));
  }
  tree_roots = new_roots;

  VESSGEN_MAXDBG(vl->IntegrityCheck();)
}


void Grower::HierarchicalGrowth()
{
  if (debug_output_every_configuration)
    DebugOutVessels(*this, str(format("with_capillaries_hit_%i") % hierarchy_level));
  FlushCapillaries();
  if (debug_output_every_configuration)
    DebugOutVessels(*this, str(format("without_capillaries_hit_%i") % hierarchy_level));
  
  UpscaleTrees();
  SplitSegmentsToOneLatticeBond(*vl);

  // obtain new terminal branches for random growth
  { OpenEnds ends;
    for (int i=0; i<vl->GetNCount(); ++i)
    {
      const VesselNode* vc = vl->GetNode(i);
      if (hasTooMuchNeighborsForAttachment(vc)) continue;
      if (vc->flags.GetBits(ARTERY) && getR(vc) > max_sprout_radius_artery) continue;
      if (vc->flags.GetBits(VEIN) && getR(vc) > max_sprout_radius_vein) continue;
      int dir = GetDirection(vc);
      ends.add(vc->lpos, dir, vc->flags & (ARTERY|VEIN));
    }
  
    RandomGrowth(ends, false);
  }
  if (debug_output_every_configuration)
    DebugOutVessels(*this, str(format("after_growth_hit_%i") % hierarchy_level));
  
  AddCapillaries();
  
  if (debug_output_every_configuration)
    DebugOutVessels(*this, str(format("after_capillaries_hit_%i")% hierarchy_level));
  
  ComputeCirculatedComponents(vl.get());
  CalcRadi();
  CalcFlow();
  
  if (debug_output_every_configuration)
    DebugOutVessels(*this, str(format("after_calcflow_hit_%i")% hierarchy_level));

#if GFFIELD_ENABLE
  InitGfDistrib();
#endif
}

void Grower::RandomGrowth(OpenEnds &ends, bool isInitial)
{
  int iteration_number_on_level = 0;
    OpenEnds ends_other;
    while (ends.size()>0)
    {
      while (ends.size()>0)
      {
        OpenEnd end = PickEnd( ends );
      AddRandomElement(end, &ends_other, !isInitial, isInitial ? (ELEM_Y | ELEM_T) : ELEM_ANY);
      }
      ends_other.swap(ends);

    if (debug_output_every_configuration && isInitial)
      DebugOutVessels(*this, str(format("growth_initial_%02i")  %iteration_number_on_level));
    if (debug_output_every_configuration && not isInitial)
      DebugOutVessels(*this, str(format("growth_hit_%1i_%02i") % hierarchy_level %iteration_number_on_level));

    VESSGEN_MAXDBG(vl->IntegrityCheck();)
    ++iteration_number_on_level;
}

}

/*---------------------------------------------------------
  --------------------------------------------------------*/

void Grower::Run(const ptree &settings, boost::function1<bool, const Grower&> callback)
{
  Init(settings);

  { // placing root nodes
  OpenEnds ends;
  GenerateRootElements(settings, ends);
  VESSGEN_MAXDBG(
    vl->IntegrityCheck();
  )
	/**
	 * adds elements until bounding box is reached
	 * in case the simulation edge is reached the OpenEnd end won't transfere to ends_other 
	 * and is therefore out of the iteration, 
	 * and size of ends will decreased after the next swap
	 * clever
	 */
  RandomGrowth(ends, true);
    }

  if (debug_output_every_configuration)
    DebugOutVessels(*this, "after_initial_growth");

  AddCapillaries();

  if (debug_output_every_configuration)
    DebugOutVessels(*this, "after_initial_capillaries");

  ComputeCirculatedComponents(vl.get());
  CalcRadi();
  CalcFlow();

  if (debug_output_every_configuration)
    DebugOutVessels(*this, "after_initial_calcflow");

  if (my::checkAbort()) return;
  if (max_hierarchy_level > 0)
  {
    for (hierarchy_level=0; hierarchy_level<max_hierarchy_level; ++hierarchy_level)
    {
      for (iteration_number_on_level = 0; iteration_number_on_level<max_num_iter; ++iteration_number_on_level)
      {
        RemodelTrees();

        if (debug_output_every_configuration)
          DebugOutVessels(*this, str(format("after_remodel_%02i_%05i") % hierarchy_level % iteration_number_on_level));

        if (!callback(boost::cref(*this))) break;
        if (my::checkAbort()) return;
      }

      HierarchicalGrowth();
      if (my::checkAbort()) return;
      
      cout << "growth " << hierarchy_level << " with vessel " << get_ld() << endl;
#if GFFIELD_ENABLE
      cout << "field " << get_field_ld() << endl;
#endif
    }
    //max_num_iter += iteration_number;
  }

  cout << "iterating ... " << endl;

  for (iteration_number_on_level = 0; iteration_number_on_level<max_num_iter; ++iteration_number_on_level)
  {
    RemodelTrees();

    if (debug_output_every_configuration)
      DebugOutVessels(*this, str(format("after_remodel_%02i_%05i") % hierarchy_level % iteration_number_on_level));

    if (!callback(boost::cref(*this))) break;
    if (my::checkAbort()) return;
  }

#if 0
  if (max_hierarchy_level > 0)
  {
    FlushCapillaries();
    UpscaleTrees();
    AddCapillaries();
    ComputeCirculatedComponents(vl.get());
    FlushCapillaries();
    CalcRadi();
    AddCapillaries();
    ComputeCirculatedComponents(vl.get());
    CalcRadi();
    CalcFlow();
    
    cout << "growth " << hierarchy_level << " with vessel " << get_ld() << endl;

    hierarchy_level = max_hierarchy_level;
    iteration_number_on_level = max_num_iter;
    iteration_number++;
    callback(boost::cref(*this));
  }
#endif
}//end grower




