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
#ifndef _REMODELER_H
#define _REMODELER_H

#include "../common.h"
#include "mwlib/math_ext.h"
#include "mwlib/ptree_ext.h"
#include "mwlib/any_tree.h"

#include "mwlib/helpers.h"
#include "../shared-objects.h"
#include "../vessels3d.h"
#include "mwlib/dynamicarray.h"
#include "mwlib/histogram.h"
#include "../growthfactor_model.h"
#include "common/calcflow.h"


#define VESSGEN_MAXDBG(x) // x
#define GFFIELD_ENABLE 1
#define WRITE_REMODELING_ACTIONS

enum ElemFlags
{
  EL_START=1,
  EL_END=2
};

struct OpenEnd
{
  OpenEnd() {}
  OpenEnd( const Int3 &p, int dir, uint flags) : p(p),dir(dir),flags(flags) {}
  Int3 p;
  int dir;
  uint flags;
};


class OpenEnds
{
  DynArray<OpenEnd> ends;
public:
  OpenEnds() {}
  
  OpenEnd pop(int i)
  {
    OpenEnd end = ends[i];
    ends[i] = ends.back();
    ends.pop_back();
    return end;
  }

  void add(const Int3 &p, int dir, uint flags)
  {
    ends.push_back(OpenEnd(p,dir,flags));
  }

  size_t size() const { return ends.size(); }
  const OpenEnd& operator[](int i) const { return ends[i]; }

  void swap(OpenEnds &other)
  {
    ends.swap(other.ends);
  }
};


/**
 * @brief An Edge has a starting and an ending node, plus a direction.
 */
struct ElemEdge
{
  ElemEdge( const Int3 &a, const Int3 &b, int dir, uint flags ) : a(a),b(b),dir(dir),flags(flags) {}
  Int3 a;
  Int3 b;
  int dir;
  my::Bitfield<uint> flags;
};
enum ElemType
{
  ELEM_I = 1,
  ELEM_T = 2,
  ELEM_Y = 4,
  ELEM_ANY = 255 // all bits on
};

typedef DynArray<ElemEdge> Element;

struct PrecElem // precomputed
{
  ElemType t;
  int dirs[3];
  PrecElem() { dirs[0]=dirs[1]=dirs[2]=-1; }
};


struct TreeRoot : public OpenEnd
{
  int len;
  TreeRoot() {}
  TreeRoot( const Int3 &p, int dir, uchar flags, int len) : OpenEnd(p,dir,flags), len(len) {}
};
typedef boost::unordered_map<Int3, TreeRoot, std::hash<Int3> > TreeRootList;

enum LatticeTypeId {
  LATTICE_TYPE_FCC,
  LATTICE_TYPE_QUAD,
};



class Grower
{
  typedef VesselList3d::LatticeData LatticeData;
  std::auto_ptr<VesselList3d> vl;
  Random rand;
  LatticeTypeId lattice_type_id;
  Float3 wdiraxes[32];
  bool dir_is_forbidden[32];
  float gf_range;
  float max_sprout_radius_vein, max_sprout_radius_artery;
  float radius_vein, radius_capi, radius_artery;
  double murray_alpha_vein, murray_alpha_artery;
  bool debug_output_every_configuration;
  bool generate_more_capillaries;
  int capillariesUntilLevel;
  BloodFlowParameters bloodFlowParameters;

#if GFFIELD_ENABLE
  // WARNING: GF is really actually like oxygen, in that it surrounds vessels
  GfModel gfmodel;
  Array3df fieldgf;
  DomainDecomposition mtboxes;
  ContinuumGrid grid;
#endif
  boost::unordered_map<std::pair<int,int>, DynArray<PrecElem> > precomputed_elements;
  
#ifdef WRITE_REMODELING_ACTIONS
  boost::unordered_map<const VesselNode*, uchar> last_remodeling_actions;
#endif
  
  TreeRootList tree_roots; // for output, not actually used in the remodeller
  
public:
  void Run(const ptree &settings, boost::function1<bool, const Grower&> callback);

  int iteration_number, // these are read only; i'm just too lazy to write getters.
      hierarchy_level,
      iteration_number_on_level,
      max_num_iter,
      max_hierarchy_level;
  const VesselList3d& get_vl() const { return *vl; }
  const LatticeData& get_ld() const { return vl->Ld(); }
  const TreeRootList& get_tree_roots() const { return tree_roots; }
  HemodynamicBounds hemodynBounds;
  int dim() const { return (Size(get_ld().Box())[2] == 1) ? 2 : 3; }
  
#if GFFIELD_ENABLE
  const LatticeDataQuad3d& get_field_ld() const { return grid.ld; }
  Array3d<float> ComputeGfSources() const;
  ConstArray3d<float> GetGf() const { return fieldgf; }
  void GetGfAtNodes(DynArray<float> &data) const;
#endif
#ifdef WRITE_REMODELING_ACTIONS
  uchar last_remodeling_action(const VesselNode* nd) const
  {
    auto it = last_remodeling_actions.find(nd);
    if (it != last_remodeling_actions.end())
      return it->second;
    else
      return (uchar)-1;
  }
#endif
  
private:
  void Init(const boost::property_tree::ptree &settings_);
  // computes a lookup table which maps from orientation and type to a list of lattice sites
  void PrecomputeElements();
  // calc radii
  void CalcRadi();  
  // random growth till full
  bool InitialRandomGrowth(const ptree& pt);
  void GenerateRootElements(const ptree& pt, OpenEnds& ends);
  bool AddBeginning(OpenEnds& ends, Int3 pos, int dir, uchar flags, int len );
  void RandomGrowth(OpenEnds &ends, bool isInitial);
  // main remodel iteration
  void RemodelTrees();    
  // adds capillaries
  void AddCapillaries();
  // removes them
  void FlushCapillaries();
  // calc flow
  void CalcFlow();
  // 
  void UpdateGf();
  //
  void SanityCheck();
  //
  void HierarchicalGrowth();
  void UpscaleTrees();


private:
#if GFFIELD_ENABLE
  float GetGf( const Float3 &pos ) const;
  Float3 GetGfGradient( const Float3 &pos ) const;
  // fill gfsource array with distribution around a source
  void InitGfDistrib();
#endif
  // create element structure
  int GetNumTrialElements(int dir, ElemType type) const;
  void MakeElement( const Int3 &pos, int dir, int index, ElemType type, Element &edges ) const;
  // check if there is space
  bool TestFitElement (const Element &elem, uchar avflags);
  // add vessels to the vessel list, no check if fits
  void AddElement( const Element &elem, uchar flags, OpenEnds *ends );
  // main routine to add a new element
  void AddRandomElement( const OpenEnd &end, OpenEnds *ends, bool useGfGradient, int typesAllowedFlags);
  // pop random element from list
  const OpenEnd PickEnd( OpenEnds &ends );
  // calcradi subroutine
  void CalcRadi( Vessel* vroot, VesselNode* nroot );
  // determines what happens at the tip nodes
  void DetermineRemodelingAction(DynArray<int> &action);
  // updates tree_roots. at the moment removes entries where roots existed before but now no more.
  void UpdateTreeRoots();
  // find a direction index for the direction of further growth. This needs special treatment since
  // any number of vessels can be attached and root nodes have their direction stored in tree_roots.
  // Will return -1 if no preferred direction can be found
  int GetDirection(const VesselNode* vc);
};

h5cpp::Group DebugOutVessels(const Grower &grower, const string &name);
typedef BasicHistogram1D<float> Histo;
double MeasureBloodVolume( const VesselList3d& vl, float dest_lattice_scale, int &n_sites, int &n_sites_occupied );
const my::Averaged<float> MeasureBranchLengths( const VesselList3d& vl, Histo &distrib, Histo &byrad );
const my::Averaged<float> MeasureRadiDistribution( const VesselList3d& vl, Histo &hrad );
void MeasureBranchNumbers( const VesselList3d &vl,  Histo &hbranch );
void MeasureRoot( const VesselList3d &vl, double &arad, double &aflow, double &arootcnt, double &vrad, double &vflow, double &vrootcnt );
const my::Averaged<float> MeasureShearDistribution( const VesselList3d& vl, Histo &hshear, Histo &byrad );

void DoOutput(h5cpp::Group root,
              const VesselList3d &vl,
              const TreeRootList &tree_roots);

#endif
