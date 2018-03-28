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
#define DOTIMING 0

#include "vesselmodel1.h"
#include "calcflow.h"
#include "mwlib/hdf_wrapper_ptree.h"
#include "mwlib/ptree_ext.h"

#ifdef DEBUG
#include <signal.h>
#define VESSEL_THREAD_CHUNK_SIZE 16
#else
#define VESSEL_THREAD_CHUNK_SIZE 1024
#endif

using namespace my;

namespace VesselModel1
{

/** @brief for a single node this gives an average 
 * of maturation
 * for the vessels linked to that node
 */
inline float GetMaturation( const VesselNode *node )
{
    float vcMat = 0;
    //myAssert(node->Count()>0);
    //T.F. I think this was a bug, maybe even related to the bad 
    //multithread behaviour during the adaption stuff????
    if(node->Count()>0)
    {
      for( int i=0; i<node->Count(); ++i )
      {
        vcMat += node->GetEdge(i)->maturation;
      }
      vcMat /= node->Count();
    }
    return vcMat;
}

/* calc timeInTumor from vessels to node */
inline float GetTimeInTumor( const VesselNode *node)
{
    float val = 0;
    for( int i=0; i<node->Count(); ++i )
    {
      val += node->GetEdge(i)->timeInTumor;
    }
    val /= node->Count();
    return val;  
}

template< class F, class T >
T GetJunctionMax( const VesselNode* node )
{
  F getter;
  bool bInit = false;
  T m=T();//either 0 or 0.0
  for( uint i=0; i<node->Count(); ++i )
  {
    const Vessel* v = node->GetEdge(i);
    T x = getter(v);
    if( !bInit || x>m ) m=x;
    bInit = true;
  }
  return m;
}

static Bitfield<uchar> CombineFlags( const VesselNode* node )
{
  Bitfield<uchar> combi = node->flags;
  for( uint i=0; i<node->Count(); ++i )
  {
    const Vessel* v = node->GetEdge(i);
    combi.AddBits( v->flags );
  }
  return combi;
}

/*-------------------------------------------------------------
// parameters
-------------------------------------------------------------*/

#define PT_ASSIGN(name) boost::property_tree::get(name, #name, pt)

void Params::assign(const ptree &pt)
{
  PT_ASSIGN(seed);
  PT_ASSIGN(pressMax);
  PT_ASSIGN(pressMin);
  PT_ASSIGN(gfVessProl);
  PT_ASSIGN(timeProlEcSprout);
  PT_ASSIGN(timeEcRadiusInflate);
  PT_ASSIGN(timeEcRadiusDeflate);
  PT_ASSIGN(timeProlEcSproutLifetime);
  PT_ASSIGN(sproutDelay);
  PT_ASSIGN(probCollapse);
  PT_ASSIGN(badaption_on_off);
  PT_ASSIGN(bRelativeShearforceCollapse);
  PT_ASSIGN(forceCollapse);
  PT_ASSIGN(forceEnlarge);
  PT_ASSIGN(radMax);
  PT_ASSIGN(radInit);
  PT_ASSIGN(distSproutMin);
  PT_ASSIGN(dematuration_rate);
  PT_ASSIGN(maturation_crit);
  PT_ASSIGN(bSproutModelSimple);
  PT_ASSIGN(vesselCompressionFactor);
  PT_ASSIGN(sproutMaxVesselWallThickness);
  PT_ASSIGN(sproutMaxVesselWallThicknessArterial);
  PT_ASSIGN(onlyReduceMaturationIfUnderperfused);
  PT_ASSIGN(isShearForceEffectLinear);
  PT_ASSIGN(bShearStressControlledDilatationAndRadiusDependentCollapse);
  //PT_ASSIGN(radiusShrinkingTendency);
  PT_ASSIGN(radMin);
  PT_ASSIGN(bRadiusSmoothing);
  PT_ASSIGN(radiusSmoothingDiffusionCoeff);
}

ptree Params::pt_defaults()
{
  return Params().as_ptree();
}

#define AS_PTREE(name) pt.put(#name, name);

ptree Params::as_ptree() const
{
  ptree pt;
  AS_PTREE(seed);
  AS_PTREE(pressMax);
  AS_PTREE(pressMin);
  AS_PTREE(gfVessProl);
  AS_PTREE(timeProlEcSprout);
  AS_PTREE(timeEcRadiusInflate);
  AS_PTREE(timeEcRadiusDeflate);
  AS_PTREE(timeProlEcSproutLifetime);
  AS_PTREE(sproutDelay);
  AS_PTREE(probCollapse);
  AS_PTREE(bRelativeShearforceCollapse);
  AS_PTREE(forceCollapse);
  AS_PTREE(forceEnlarge);
  AS_PTREE(radMax);
  AS_PTREE(radInit);
  AS_PTREE(distSproutMin);
  AS_PTREE(dematuration_rate);
  AS_PTREE(maturation_crit);
  AS_PTREE(bSproutModelSimple);
  AS_PTREE(badaption_on_off);
  AS_PTREE(vesselCompressionFactor);
  AS_PTREE(sproutMaxVesselWallThickness);
  AS_PTREE(sproutMaxVesselWallThicknessArterial);
  AS_PTREE(onlyReduceMaturationIfUnderperfused);
  AS_PTREE(isShearForceEffectLinear);
  AS_PTREE(bShearStressControlledDilatationAndRadiusDependentCollapse);
  //AS_PTREE(radiusShrinkingTendency);
  AS_PTREE(radMin);
  AS_PTREE(bRadiusSmoothing);
  AS_PTREE(radiusSmoothingDiffusionCoeff);
  return pt;
}
/** 
 * \param maturation_crit	Vessels with maturation above this threshold 
 * 				cannot collapse anymore
 * 				initial wall thickness should be around 4 micron, so this is maturation_crit,
 *				so capillaries should immediately collapse
 * 				note that also vessels within normal tissue can collapse if their maturation value is < maturation_crit
 * \param dematuration_rate	images in (holash et al.) indicate approximately 0.04,
 * 				based on the time it takes to destroy a 25micron *radius* vessel
 * 				maturation of 25 micron vessel ca 16, 5 micron vessel ca 4.5
 * \param pressMax		kPa, pressure must be larger in order to get reasonable flow/shearforce
 * \param timeProlEcSprout	in hours per added segment
 */

Params::Params()
{
  seed = 12348;
  pressMax = 13.0;
  pressMin = 0.0;
  gfVessProl = 0.0005;
  timeProlEcSprout  = 2;  // 2h per added segment
  timeEcRadiusInflate = 60; // after 360 hours, vessel radius should be ca 27micron
  timeEcRadiusDeflate = 60.;
  timeProlEcSproutLifetime = 50;  // the time until sprouts regress
  sproutDelay = 1;
  probCollapse = 0.5; //0.05f; //0.05; //0.01;
  bRelativeShearforceCollapse = false;
  forceCollapse = 1.0/1000.0;  // kPa
  forceEnlarge  = 10./1000.0;
  radMax    = 25.0;
  radInit   = 4.0;
  distSproutMin =  3;
  dematuration_rate = 0.05;
  maturation_crit = 3.;
  bSproutModelSimple = false;
  badaption_on_off = false;
  sproutMaxVesselWallThicknessArterial = 50.;
  sproutMaxVesselWallThickness = 50.; // infinity
  // radius under max pressure is = reference_r * vesselCompressionFactor
  vesselCompressionFactor = 0.5f;
  onlyReduceMaturationIfUnderperfused = false;
  isShearForceEffectLinear = true;
  bShearStressControlledDilatationAndRadiusDependentCollapse = false;
  //radiusShrinkingTendency = 0.8;
  radMin = 2.0;
  bRadiusSmoothing = false;
  radiusSmoothingDiffusionCoeff = 100./timeEcRadiusInflate;
}



/*-------------------------------------------------------------
// vessels
-------------------------------------------------------------*/


uint Model::GetThreadRandomSeed()
{
  main_mutex.lock();
  uint s = main_rnd.Get();
  main_mutex.unlock();
  return s;
}


void Model::Init(VesselList3d *vl_, const Params &params_, Callbacks &callbacks_)
{
  callbacks = callbacks_;
  params = params_;
  vl = vl_;
  //cache vessel lattice data
  m_ld = &vl->Ld(); 
  //we fill lattice_axis_directions with the world directions
  //this is needed to calculate dot products
  GetWorldDirectionAxes<LatticeData>(*m_ld, lattice_axis_directions, GET_WORLD_DIRECTION_AXES_NORMALIZE);
  
  myAssert(callbacks.getGf &&
           callbacks.getGfGrad &&
           callbacks.getPress &&
           callbacks.getTumorDens);

  main_rnd.Init( params.seed>0 ? params.seed : 12345 );

  //set initial values for vessels
  int ecnt = vl->GetECount();
  for( int i=0; i<ecnt; ++i )
  {
    // WARNING only model specific data is set. everyhting else is set by the composit model
    Vessel* v = vl->GetEdge(i);
    v->timeInTumor = 0;
    v->timeSprout = -1;
    v->reference_r = v->r;
    v->initial_f = v->f;
  }
  num_iteration = 0;
  time = 0;
}

void myprint(boost::property_tree::ptree const& pt)
{
    using boost::property_tree::ptree;
    ptree::const_iterator end = pt.end();
    for (ptree::const_iterator it = pt.begin(); it != end; ++it) {
        std::cout << it->first << ": " << it->second.get_value<std::string>() << std::endl;
        myprint(it->second);
    }
}

void Model::DoStep(double dt, const BloodFlowParameters *bfparams)
{
  this->dt = dt;
  myAssert(dt == 1.);
  
#ifndef USE_ADAPTION
  {
    GenerateSprouts();
#ifndef NDEBUG
    std::cout << "collapse called" << std::endl; std::cout.flush();
#endif
    CollapseVessels();
#ifndef NDEBUG
    std::cout << "collapse ened" << std::endl; std::cout.flush();
#endif
    if (IS_DEBUG) vl->IntegrityCheck();
#ifndef NDEBUG
    std::cout << "integrity ended" << std::endl; std::cout.flush();
#endif
    EnlargeVessels();
#ifndef NDEBUG
    std::cout << "enlarge ended" << std::endl; std::cout.flush();
#endif
    SmoothVessels();
#ifndef NDEBUG
    std::cout << "smooth ended" << std::endl; std::cout.flush();
#endif
    MaturateVessel();
#ifndef NDEBUG
    std::cout << "maturate ended" << std::endl; std::cout.flush();
#endif
  }
#endif
  if(num_iteration%10 == 0)
  {
#ifndef NDEBUG
    cout << "optimize called" << endl;
#endif
    Optimize(vl);
#ifndef NDEBUG
    cout << "optimize finished" << endl;
#endif
  }
  if (IS_DEBUG) vl->IntegrityCheck();
  ClassifyTumorVessels();
  //DebugOutVessels(*vl, "firstiter");
  ComputeCirculatedComponents(vl);
  ++num_iteration;
  time += dt;
  this->dt = 0.;
}

/*
 * if various checks on the two node point succeded
 * we can insert a new sprout. New sprouts are labeled
 * as capillaries.
 */
Vessel* Model::AddSproutVessel(VesselNode* vc, VesselNode* dst_vc)
{
  const Vessel* v = vc->GetEdge(0);
  const bool bSproutMigration = vc->Count()==1;

  myAssert( dst_vc==NULL || vl->FindVessel(vc->lpos,dst_vc->lpos)==NULL );

  uint bridgeFlags = vc->flags;
  if( dst_vc ) bridgeFlags |= dst_vc->flags;

  int timeSprout  = GetJunctionMax<VData::PropTimeSprout,int>(vc);
  int timeInTumor = GetJunctionMax<VData::PropTimeInTumor,int>(vc);
  float initial_f = GetJunctionMax<VData::PropInitialF,float>(vc);
  float initial_h = GetJunctionMax<VData::PropHematocrit, float>(vc);
  Bitfield<uchar> flags = CombineFlags(vc);
  flags.DelBits(ARTERY|VEIN);
  flags.AddBits(CAPILLARY);

  Vessel* vBridge = vl->InsertVessel(vc, dst_vc);
  vBridge->r = params.radInit;
  vBridge->flags = flags;
  vBridge->maturation = GetInitialThickness( vBridge->r );
  vBridge->initial_f = initial_f;
  vBridge->hematocrit = initial_h;
  vBridge->reference_r = vBridge->r;

  if( bSproutMigration ) // sprout migration from an existing sprout
  {
    vBridge->timeSprout = v->timeSprout;
    vBridge->timeInTumor = v->timeInTumor;
  }
  else // creation of a new sprout, possibly from an existing sprout
  {
    vBridge->timeSprout = timeSprout;
    vBridge->timeInTumor = timeInTumor;
  }
  // the following conditions is true if the sprout starts from a non sprouting vessel
  if( vBridge->timeSprout<0 ) vBridge->timeSprout=params.timeProlEcSproutLifetime;
  return vBridge;
}


/* 
 *this generates a single sprout at position &pos in direction &forward_dir
 */
Vessel* Model::GenerateSprout(Random &rnd, const Int3 &pos, const Float3 &forward_dir)
{
  Float3 gfgrad = GetGfGrad(pos);
  
  bool has_fwd = squaredNorm(forward_dir)>1.e-13;

  double max_score = -std::numeric_limits<double>::max();
  int best_dir = -1;
  double scores[128];
  
  /* 
   * for the given &pos we look in all possible directions
   * and 
   * judge this direction by growth factor field
   */
  for (int dir=0; dir<Ld().NbCount(); ++dir)
  {// probe directions
    //initialize
    scores[dir] = 0;
    //get neighbor pos
    const Int3 nbpos = Ld().NbLattice(pos, dir);
    //ignore neigbors outside the lattice and neigbors already occupied by vessels
    if (!Ld().IsInsideLattice(nbpos) || vl->FindVessel(pos,nbpos)) continue;
    
    bool ok_strict_fwd = !has_fwd || lattice_axis_directions[dir].dot(forward_dir) > 0.9;
    bool ok_fwd = !has_fwd || lattice_axis_directions[dir].dot(forward_dir) > -0.1;
    Vessel* dst_v;
    VesselNode* dst_vc;
    //check if a vessel or a node is at the neigbor pos 
    vl->FindAny(nbpos, dst_v, dst_vc);
    //if there if something
    if (dst_v || dst_vc)
    {
      //ask for maturation
      double maturation = dst_v ? dst_v->maturation : GetMaturation(dst_vc);
      Bitfield<uchar> flags  = dst_v ? dst_v->flags : CombineFlags(dst_vc);

      bool isArtery = flags.GetBits(ARTERY);
      //if the artery or vein is already mature enough as defined
      //defined by the model parameters this direction is ignored
      if ( ((maturation > params.sproutMaxVesselWallThicknessArterial) && isArtery) ||
          ((maturation > params.sproutMaxVesselWallThickness) && !isArtery)) continue;
      
      //ignore nodes with more than 3 connected vessels
      if (dst_vc && dst_vc->Count()>=3) continue;
    }
    
    // reject if not forward facing and no adjacent node is there
    if (!dst_vc && !ok_strict_fwd) continue;
    if (dst_vc && !ok_fwd) continue;
    
    // get a score from the gradient and save it to the array
    // up the gradient
    double score = std::max<double>(0, lattice_axis_directions[dir].dot(safe_normalized(gfgrad))); 
    if (dst_vc)
      score += 100.;  // make it go to the adjacent node
    scores[dir] = score;
    //remember max_score and best_dir
    if (score > max_score)
    {
       max_score = score;
       best_dir = dir;
    }
  }// end probe directions
  //now we know the best direction to go!!!
  //failsafe
  if (best_dir < 0 || max_score<=0.) return NULL;

  NormalizeSumNorm(scores, Ld().NbCount());
  //now we draw a random direction according to the scores
  best_dir = RandomPick(rnd, Ld().NbCount(), scores);
  //get lattice point in best_dir
  Int3 nbpos = Ld().NbLattice(pos, best_dir);
  
  VesselNode* src_vc=NULL, *dst_vc=NULL;
  Vessel* src_v=NULL, *src_v2=NULL, *dst_v=NULL, *dst_v2=NULL;
  // find the source vessel or node if there are some in vicinity
  vl->FindAny(pos, src_v, src_vc); 
  myAssert(src_v || src_vc);
  //if we found no node but a vessel only
  //than we need to split the vessl into two
  if (!src_vc && src_v)
  {
    vl->SplitVessel(src_v, FindPositionOnVessel(Ld(), src_v, pos), src_v2, src_vc);
  }
  //find vessel or node in the vicinity of the lattice position
  vl->FindAny(nbpos, dst_v, dst_vc);
  //again, if we found a vessel but no node, we need to split it
  if (dst_v && !dst_vc) {
    vl->SplitVessel(dst_v, FindPositionOnVessel(Ld(), dst_v, nbpos), dst_v2, dst_vc);
  }
  //if the destination node dst_vc in not yet end point of a vessel
  //we can use it for the end point of the new sprout
  if (!dst_vc)
    dst_vc = vl->InsertNode(nbpos);

  Vessel* sprout_v = AddSproutVessel(src_vc, dst_vc); // boom??!!
  return sprout_v;
}


bool Model::CheckCanSprout(Random &rnd, const Int3 &lpos, const VesselNode* nd, const Vessel* v)
{
  double maturation = v ? v->maturation : GetMaturation(nd);
  Bitfield<uchar> flags  = v ? v->flags : CombineFlags(nd);

  bool isCapillary = flags.GetBits(CAPILLARY);
  bool isArtery = flags.GetBits(ARTERY) && !isCapillary; // just for safety we exclude explicitly labeled capillaries.
  bool isVein   = flags.GetBits(VEIN) && !isCapillary;
  if ( (maturation > params.sproutMaxVesselWallThicknessArterial && isArtery) ||
       (maturation > params.sproutMaxVesselWallThickness && !isVein)) return false;

  //if ( maturation > params.sproutMaxVesselWallThickness ) return false;

  double timeInTum  = v ? v->timeInTumor : GetTimeInTumor(nd);
  if ( timeInTum > params.sproutDelay) return false;

// #ifdef DEBUG
//   printf("gf is: %f\n", GetGfConc(lpos));
//   printf("threshold is: %f\n", params.gfVessProl);
// #endif
  if( GetGfConc(lpos)<params.gfVessProl)
  {
    return false;
  }

  return true;
}

/*
 * supervising function wich controlles all single sprouts
 */
void Model::GenerateSprouts()
{
  FUNC_TIMING_START
  typedef LatticeData::SiteType SiteType;
  tbb::spin_mutex mutex;
  DynArray<SiteType> sitesSprout;
  DynArray<VesselNode*> vcExtendSprout;
  
  #pragma omp parallel
  {
    Random rnd(GetThreadRandomSeed());
    DynArray<SiteType> th_sitesSprout(1024, ConsTags::RESERVE);
    DynArray<VesselNode*> th_vcExtendSprout(1024, ConsTags::RESERVE);

    //edge stuff
    #pragma omp for schedule(dynamic, VESSEL_THREAD_CHUNK_SIZE)
    for( int i=0; i<vl->GetECount(); ++i )
    {
      Vessel* v = vl->GetEdge(i);
      //for all which are sprouts
      if( v->timeSprout>=0 )
      {
	//if all this we decrease the max sprout time!
        if(!v->IsCirculated() && v->flags.GetBits(CONNECTED) && v->timeSprout>=0)
          --v->timeSprout; // NOTE: decrease by dt instead???!!!!!
        else
          v->timeSprout = -1; 
      }
    }
    #pragma omp for schedule(dynamic, VESSEL_THREAD_CHUNK_SIZE)
    for(int i=0; i<vl->GetNCount(); ++i)
    {//node stuff
      VesselNode* nd = vl->GetNode(i);
      //ignore unconnected stuff
      if (nd->Count() <= 0) continue;
      //get lattice position
      Int3 lpos = nd->lpos;
      //if we have excatly two node, so it is embedded in a vessel
      if( nd->Count()==2 )
      {
        if (CheckCanSprout(rnd, lpos, nd, NULL))
          th_sitesSprout.push_back(Ld().LatticeToSite(lpos)); 
      }
      //if it has no link
      else if( nd->Count()==1 )
      {
	//if vessel is no sprout, maybe at boundary or so
        if (nd->GetEdge(0)->timeSprout <= 0 ) continue;
	//throw the dice on wether this sprout will grow
        if (rnd.Get01()>1.0/params.timeProlEcSprout ) continue;
        th_vcExtendSprout.push_back(nd);
      }
    }//end node stuff

    #pragma omp for schedule(dynamic, VESSEL_THREAD_CHUNK_SIZE)
    for( int i=0; i<vl->GetECount(); ++i )
    {//edge loop
      Vessel* v = vl->GetEdge(i);
      //inital lpos
      Int3 lpos = v->NodeA()->lpos;
      //walk along the vessel
      for( int i=1; i<v->len-1; ++i )
      {
	//next lpos for walk
        lpos = Ld().NbLattice(lpos,v->dir);
        if (CheckCanSprout(rnd, lpos, NULL, v))
          th_sitesSprout.push_back( Ld().LatticeToSite(lpos) );
      }
    }//end edge loop

    {
      mutex.lock();
      sitesSprout.insert(sitesSprout.end(), th_sitesSprout.begin(), th_sitesSprout.end());
      vcExtendSprout.insert(vcExtendSprout.end(), th_vcExtendSprout.begin(), th_vcExtendSprout.end());
      mutex.unlock();
      th_sitesSprout.remove_all();
      th_vcExtendSprout.remove_all();
    }
  }// end #pragma omp parallel

  Random rnd(GetThreadRandomSeed());
  //create new sproutings
  //const int numIterations = sitesSprout.size()/params.timeProlEc;
  const int numIterations = (int)sitesSprout.size();
  for( int i=0; i<numIterations; ++i )
  {
    const SiteType siteSprout  = sitesSprout[i];
    const Int3 lposSprout = Ld().SiteToLattice( siteSprout );

    Vessel* v=NULL; VesselNode* vc=NULL;
    vl->FindAny( lposSprout, v, vc );
    myAssert( (v==NULL)^(vc==NULL) );  // either there is a vessel or a node

    if( rnd.Get01()>1.0/params.timeProlEcSprout ) continue;
    int dist = DistanceToJunction( lposSprout, params.distSproutMin );
    
    if( dist<params.distSproutMin ) continue;
    //single new sprout at lposSprout
    GenerateSprout(rnd, lposSprout, Float3(0));
  }
 
  // extend existings sprouts
  for( int iteration=0; iteration<vcExtendSprout.size(); ++iteration )
  {
    VesselNode* vc = vcExtendSprout[ iteration ];
    //a sprouting node is not allowed to have more than one neigbor, 
    //and it is not allowed be at the boundary
    if( vc->Count()!=1 || vc->IsBoundary() ) continue;

    const Int3 lposSprout = vc->lpos;
    Vessel* v = vc->GetEdge(0);

    //this gets the world vector reprensenting to the vessel
    Float3 forward_dir = lattice_axis_directions[v->dir];
    //direction always points from NodeA to NodeB
    //so if NodeA == vc this es exactly opposite
    if (v->NodeA() == vc)
      forward_dir = -forward_dir;

    Vessel* sprout_v = GenerateSprout(rnd, lposSprout, forward_dir);
    //if GenerateSprout fails we reset sprouting time of this vessel
    if (!sprout_v)
      v->timeSprout = -1;

  }
  FUNC_TIMING_END_OS(my::log())
}


inline double Model::StabilityAmountDueToShearStress(const Vessel* v)
{
  double factor_fc = 0.f;
  if (params.isShearForceEffectLinear)
  {
    if( params.bRelativeShearforceCollapse )
    {
      factor_fc = Cut<float>((v->f/v->initial_f/params.forceCollapse),0.f,1.f);
      //std::printf("factor_fc: %f, v->f: %f, v->initial_f: %f, params.forceCollapse: %f\n", factor_fc,v->f,v->initial_f, params.forceCollapse);
    }
    else
    {
      factor_fc = Cut<float>( v->f/params.forceCollapse, 0.f, 1.f );
    }
  }
  else
  {
    if (params.bRelativeShearforceCollapse)
    {
      factor_fc = (v->f > (v->initial_f*params.forceCollapse)) ? 1. : 0.;
    }
    else
    {
      factor_fc = (v->f > params.forceCollapse) ? 1. : 0.;
    }
  }
  return factor_fc;
}


void Model::CollapseVessels()
{
  FUNC_TIMING_START

  DynArray<Vessel*> toKill;
  tbb::spin_mutex mutex;

  #pragma omp parallel
  {
    Random rnd(GetThreadRandomSeed());
    DynArray<Vessel*> th_toKill(1024, ConsTags::RESERVE);

    #pragma omp for schedule(dynamic, VESSEL_THREAD_CHUNK_SIZE)
    for( int i=0; i<vl->GetECount(); ++i )
    {
      Vessel* v = vl->GetEdge( i );
      //what does the b stand for???? --> bool
      bool bKill = false;
#if 0
      if (params.badaption_on_off)
      {//do adaption stuff
	//do not consider sprouts here
	if( v->timeSprout>=0 ) continue;
	//if to small or not circulated --> kill
	if (v->reference_r<params.radMin || !v->IsCirculated())
	{
	  //walk along the vessel until its length is reached
	  for(int i=0; i<v->len-1 && !bKill; ++i)
	  {
	    // always kill
	    bKill = true;
	  }
	}
      }
#endif
      //else{// the stuff mw implemented 
	if ((bool)params.bShearStressControlledDilatationAndRadiusDependentCollapse)
	{
	  //do not consider sprouts here, v->timeSprout initialized by -1
	  if( v->timeSprout>=0 ) continue;
	  //if to small or not circulated --> kill
	  if (v->reference_r<params.radMin || !v->IsCirculated())
	  {
	    //walk along the vessel until its length is reached
	    for(int i=0; i<v->len-1 && !bKill; ++i)
	    {
	      //get random between 0 and 1
	      double r = rnd.Get01();
	      // sometimes we kill
	      bKill |= r<params.probCollapse;//bKill = bKill | r<params.probCollapse
	    }
	  }
	}
	else
	{
	  //do not kill sprouts
	  if( v->timeSprout>=0 ) continue;
	  //do not kill well established vessels
	  if( v->maturation>params.maturation_crit ) continue;
// 	  std::printf("v->timeSprout: %f \n",v->timeSprout);
// 	  std::printf("v->maturation: %f \n",v->maturation);
// 	  std::printf("params.maturation_crit: %f \n",params.maturation_crit);
	  // get pos of node A
	  Int3 pos = v->LPosA();
	  // while the vessel has a length; pos = pos of neigbor in direction
	  for( int i=0; i<v->len-1 && !bKill; ++i, pos=Ld().NbLattice(pos,v->dir) )
	  {
	    //role the dice
	    double r = rnd.Get01();
	    if(!v->IsCirculated())
	    {
	      //role the dice, sometimes this vessel will be killed here
	      bKill |= r<params.probCollapse;
	    }
	    else
	    {
	      //even if it is circulated it will be killed here
	      double factor_fc = StabilityAmountDueToShearStress(v);
	      double pcoll = Lerp<double>( factor_fc, params.probCollapse, 0 );
	      //std::printf("pcoll: %f\n",pcoll);
	      bKill |= r<pcoll;
	    }
	  }
	}
     // }
      if(bKill) th_toKill.push_back(v);
    }
    //parallel kill ;-)
    mutex.lock();
    toKill.insert(toKill.end(), th_toKill.begin(), th_toKill.end());
    mutex.unlock();
    th_toKill.remove_all();
  }
  #pragma omp barrier
 
  for( int i=0; i<toKill.size(); ++i )
  {
    vl->DeleteVessel( toKill[i] );
  }

  if(toKill.size()>0)
  {
    std::cout<<"KILL: #" << toKill.size() << " vessels" <<endl;
  }
  else
  {
    std::cout<<"Nothing to KILL!"<<endl;
  }

  FUNC_TIMING_END_OS(my::log())
}


// enlarge for a length segment equivalent to one lattice bond length
// distribute enlarged area to length of whole vessel, increasing its radius
// several steps needed due to checks for GF
inline void Model::EnlargeRadiusSubstep(Random &rnd, Vessel* v, float fr)
{
  if (v->reference_r > params.radMax) return;
  const float areaEC = 10.0*10.0;
  const float segment_length = v->WorldLength(Ld());
  float numECs = (my::mconst::fpi2()*v->reference_r*Ld().Scale()) / areaEC;
  float dn = 0.;
  while (numECs > 0.)
  {
    if (rnd.Get01f() < numECs && rnd.Get01()<1.0/params.timeEcRadiusInflate)
    {
      dn += 1.;
    }
    numECs -= 1.;
  }
  const float dr = (dn*areaEC)/(my::mconst::fpi2()*segment_length);
  v->reference_r += fr * dr;
}


float Model::GetAverageGrowthfactorOverSegment(const Vessel* v)
{
  float gf = 0.;
  Int3 pos = v->LPosA();
  Int3 pos2 = Ld().NbLattice(pos, v->dir);
  for (int l=0; l<v->len-1; ++l )
  {
    Float3 wpos = 0.5*(Ld().LatticeToWorld(pos)+Ld().LatticeToWorld(pos2));
    float localGf = GetGfConc(wpos);
    gf += localGf;
  }
  gf *= 1./v->WorldLength(Ld());
  return gf;
}


void Model::EnlargeVessels()
{
  FUNC_TIMING_START
  //why is ln(2) needed here an 1/parmas.time everywhere else???
  // if params.bShearStressControlledDilatationAndRadiusDependentCollapse == true then,
  // inflation and deflation rates are interpreted as halflife and doubling time of endothelial cells.
  // ln(2) is required for conversion to rate constants. Because df/dt = k f. And we want f(T) = 2 f(0).
  const double infl_rate = my::mconst::ln2()/params.timeEcRadiusInflate;
  const double defl_rate = my::mconst::ln2()/params.timeEcRadiusDeflate;
  #pragma omp parallel
  {
    Random rnd(GetThreadRandomSeed());

    #pragma omp for schedule(dynamic, VESSEL_THREAD_CHUNK_SIZE)
    for( int iteration=0; iteration<vl->GetECount(); ++iteration )
    {
      Vessel *v = vl->GetEdge( iteration );

      if( v->flags.GetBits(WITHIN_TUMOR) ) ++v->timeInTumor; //since time is measure in integers we can do this

      if (v->timeInTumor>params.sproutDelay && v->timeSprout<0 && v->IsCirculated())
      {
        if (params.bShearStressControlledDilatationAndRadiusDependentCollapse)
        {
            //double gf_conc = GetAverageGrowthfactorOverSegment(v);
            //double gf_growth_signal = gf_conc > params.gfVessProl ? 0.5 : 0.;
          
            double r = v->reference_r;
	    double growth_signal; //either 1. or 0.
            growth_signal = v->f > (v->initial_f*params.forceEnlarge) ? 1. : 0.;
	    double shrink_signal; //either 1. or 0.
            shrink_signal = v->f < (v->initial_f*params.forceCollapse) ? 1. : 0.;
            //double b = my::smooth_heaviside_sin(v->f/v->initial_f-params.forceGrow, 0.1)*my::mconst::ln2()/params.timeProlEcEnlarge;
            //double b = my::smooth_heaviside_sin(v->f/1.e-1, 0.1)*my::mconst::ln2()/params.timeProlEcEnlarge;
            //double b = 0.;
            //double dr = -a*r*r+b*r;
            //double dr = r*a*(1 - (1.-f)*r/20.);
            //double dr = r*a*(f - params.radiusShrinkingTendency);
            double dr = r*(infl_rate*(growth_signal) - defl_rate*shrink_signal);
            if (v->reference_r > params.radMax && dr>0.) dr=0.;
            v->reference_r += dr*this->dt;
        }
        else
        {
          Int3 pos = v->LPosA();
          Int3 pos2 = Ld().NbLattice(pos, v->dir);
          for (int l=0; l<v->len-1; ++l )
          {
            Float3 wpos = 0.5*(Ld().LatticeToWorld(pos)+Ld().LatticeToWorld(pos2));
            if( GetGfConc(wpos)>params.gfVessProl )
            {
              EnlargeRadiusSubstep(rnd, v, 1.0 );
            }
            pos=pos2;
            pos2=Ld().NbLattice(pos2, v->dir);
          }
        }
      } 
      // allow circumferential growth
      if(params.vesselCompressionFactor<1.)
      {
        myAssert(params.vesselCompressionFactor > 0.01);
        Int3 pos = v->LPosA();
        Int3 pos2 = Ld().NbLattice(pos, v->dir);
        float tiss_press = 0.0f;
        int l=0;
        while( l<v->len )
        {
          tiss_press += GetTissPress(pos);
          ++l;
          pos=pos2;
          pos2=Ld().NbLattice(pos2,v->dir);
        }
        tiss_press *= (1.0f/v->len);
        v->r = Lerp<float,float>(Cut<float>(tiss_press,0.,1.),v->reference_r,v->reference_r*params.vesselCompressionFactor);
      }
      else // no vessel compression
      {
        v->r = v->reference_r;
      }
    }//end #pragma omp for schedule(dynamic, VESSEL_THREAD_CHUNK_SIZE)
  }//end #pragma omp parallel
  FUNC_TIMING_END_OS(my::log())
}


void Model::SmoothVessels()
{
  FUNC_TIMING_START
  if (params.bRadiusSmoothing)
  {
    int ecnt = vl->GetECount();
    DynArray<double> dr(ecnt, 0.);

    #pragma omp parallel for schedule(dynamic, VESSEL_THREAD_CHUNK_SIZE)
    for (int i=0; i<ecnt; ++i)
    {
      const Vessel *v = vl->GetEdge(i);
      if (v->timeSprout>0 || v->timeInTumor<=0 || !v->IsCirculated()) continue;

      double l = v->WorldLength(vl->Ld());
      double dr_v = 0.;
      for (int i=0; i<2; ++i)
      {
        const VesselNode* n = v->GetNode(i);
        for (int j=0; j<n->Count(); ++j)
        {
          const Vessel* w = n->GetEdge(j);
          if (w == v || !w->IsCirculated()) continue;
          double lw = w->WorldLength(vl->Ld());
          dr_v += (w->reference_r-v->reference_r)/(lw+l);
        }
      }
      dr[v->Index()] = dr_v*this->dt*params.radiusSmoothingDiffusionCoeff*2.; // factor 2 because (lw+l) counts the way to the neighbors twice.
      myAssert((v->reference_r+dr[v->Index()])>0.);
    }

    for (int i=0; i<ecnt; ++i)
    {
      Vessel *v = vl->GetEdge(i);
      v->reference_r += dr[v->Index()];
    }
  }
  FUNC_TIMING_END_OS(my::log())
}


/* maturation is probably only relevant for vessel stability and 
 * determination of the wall permeability to interstitial fluid. */
void Model::MaturateVessel()
{
  FUNC_TIMING_START
  int ecnt = vl->GetECount();
  
  #pragma omp parallel for schedule(dynamic, VESSEL_THREAD_CHUNK_SIZE)
  for( int i=0; i<ecnt; ++i )
  {
    Vessel *v = vl->GetEdge( i );
    if(v->flags.GetBits(WITHIN_TUMOR) && v->maturation>0)
    {
      if (params.onlyReduceMaturationIfUnderperfused)
      {
        double f = StabilityAmountDueToShearStress(v);
        myAssert(f>=0. && f<=1.);
        double rate = (1.-f)*params.dematuration_rate;
        v->maturation -= rate;
      }
      else
      {
        v->maturation -= params.dematuration_rate;
      }
      if( v->maturation<0 ) v->maturation=0;
    }

  }
  FUNC_TIMING_END_OS(my::log())
}



void Model::ClassifyTumorVessels()
{
  FUNC_TIMING_START

  #pragma omp parallel for schedule(dynamic, VESSEL_THREAD_CHUNK_SIZE)
  for( int i=0; i<vl->GetECount(); ++i )
  {
    Vessel* v = vl->GetEdge(i);
    myAssert(v->r > 1.);
    // vessels are only ever set to be within the tumor, never unset
    if( v->flags.GetBits(WITHIN_TUMOR) ) continue;
    // check neighbor sites of vessel
    uint ntc=0,ns=0;
    Int3 lpos = v->LPosA();
    float avg = 0.;
    for( int l=0; l<v->len; ++l,lpos=Ld().NbLattice(lpos,v->dir))
    {
      avg += GetTumorDens(lpos);
    }
    /*v->len is in units of the base lattice length
     * so this is an average per unit volume!
     */
    avg /= v->len;
    if (avg > 0.5)
    {
      v->flags.AddBits(WITHIN_TUMOR);
    }
  }
  FUNC_TIMING_END_OS(my::log())
}



int Model::DistanceToJunction( const Int3 &pos, int maxdist )
{
  Vessel* vstart=NULL; VesselNode* vcstart=NULL;
  vl->FindAny(pos,vstart,vcstart);

  myAssert( (vstart==NULL)^(vcstart==NULL) );

  int posOnVess = -1;
  if(vstart) posOnVess=FindPositionOnVessel(Ld(), vstart, pos );

  return FindDistanceToJunction(vstart,posOnVess,vcstart,maxdist);
}


/*-------------------------------------------------------------
-------------------------------------------------------------*/

}
