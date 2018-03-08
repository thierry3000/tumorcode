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
  //#define DOTIMING 1
#include "shared-objects.h"
#include "vessels3d.h"
#include "continuum-flow.h"
#include "calcflow.h"
#include "levelset.h"
#include "continuum-grid.h"


#include <boost/tuple/tuple.hpp>
#include <boost/foreach.hpp>
#include <boost/property_tree/info_parser.hpp>
#include <boost/program_options.hpp>

using namespace my;

#define VESSEL_INTEGRITY_CHECK_SWITCH(x)

/*-----------------------------------------------------
  -----------------------------------------------------*/
/** 
 * @brief Decides which tumor type is present
 * myexception is thrown once no type is recognized
 */
class myexception: public std::exception
{
  virtual const char* what() const throw()
  {
    return "tumor type not implemented on cpp side";
  }
};

TumorTypes determineTumorType(boost::optional<H5::Group> tumorgroup)
{
// #ifdef DEBUG
//     cout<<format("tumorgroup.attrs().exists(\"TYPE\"): %s\n") % (bool)(tumorgroup.attrs().exists("TYPE"));
//     cout<<format("tumorgroup.attrs().get<string>(\"TYPE\"): %s\n") % tumorgroup.attrs().get<string>("TYPE");
// #endif
  
  TumorTypes tumortype;
  if(tumorgroup)//means it is initialized
  {
    //std::string detailedTumorDescription = tumorgroup.attrs().get<std::string>("TYPE");
    /*
    try {  // to determine if the dataset exists in the group
        dataset = new DataSet( group->openDataSet( "Compressed_Data" ));
    }
    catch( GroupIException not_found_error ) {
        cout << " Dataset is not found." << endl;
    }
    cout << "dataset \"/Data/Compressed_Data\" is open" << endl;
    */
    
    try
    {
      string detailedTumorDescription;
      //now the BUG is here!
      readAttrFromH5(*tumorgroup, string("TYPE"), detailedTumorDescription);
      if( detailedTumorDescription == "faketumor" )
      {
        tumortype = TumorTypes::FAKE;
      }
      else if ( detailedTumorDescription == "BulkTissueFormat1")
      {
        tumortype = TumorTypes::BULKTISSUE;
      }
    }
    catch( H5::AttributeIException not_found_error)
    {
      cout << " Attribute TYPE not found" << endl;
    }
  }
//     if( tumorgroup->attrs().exists("TYPE") )
//     {
//       std::string detailedTumorDescription = tumorgroup->attrs().get<std::string>("TYPE");
//       if( detailedTumorDescription == "faketumor" )
//       {
//         tumortype = TumorTypes::FAKE;
//       }
//       else if ( detailedTumorDescription == "BulkTissueFormat1")
//       {
//         tumortype = TumorTypes::BULKTISSUE;
//       }
//     }
//   }
//   else
//   {
//     tumortype = TumorTypes::NONE;
//   }
//   try
//     {
//       std::string detailedTumorDescription = tumorgroup.attrs().get<std::string>("TYPE");
//       if( detailedTumorDescription == "faketumor" )
//       {
//         tumortype = TumorTypes::FAKE;
//       }
//       else if ( detailedTumorDescription == "BulkTissueFormat1")
//       {
//         tumortype = TumorTypes::BULKTISSUE;
//       }
//       else
//       {
// //         throw myexception();
//         tumortype = TumorTypes::NONE;
//       }
//     }
//     catch(std::exception& e)
//     {
//       cout << "reading tumor type from hdf failed because of: " << e.what() << '\n';
//     }
  return tumortype;
//   tumortype = TumorTypes::FAKE;
//   return TumorTypes::FAKE;
}

void SetupTissuePhases(TissuePhases &phases, const ContinuumGrid &grid, DomainDecomposition &mtboxes, boost::optional<H5::Group> tumorgroup)
{
  TumorTypes tumortype = determineTumorType(tumorgroup);
  if( tumortype == TumorTypes::NONE )
  {
    // if there is no tumorgroup provided, we fill everything with normal tissue
    // 1 means singel tissue phase
    phases = TissuePhases(1, grid.Box());
    phases.phase_arrays[TISSUE].fill(1.);
  }
  else if( tumortype == TumorTypes::FAKE)
  {
    //there is a faketumor group provided
    //double tumor_radius = tumorgroup->attrs().get<double>("TUMOR_RADIUS");
    double tumor_radius;
    readAttrFromH5(*tumorgroup, "TUMOR_RADIUS",tumor_radius);
    phases = TissuePhases(3, grid.Box());//consistent labeling problems Tissue is 2 in enumeration
    #pragma omp parallel
    {
      BOOST_FOREACH(const BBox3 &bbox, mtboxes.getCurrentThreadRange())
      {
        FOR_BBOX3(p, bbox)
        {
          Float3 wp =  grid.ld.LatticeToWorld(p);
          double r = wp.norm();
          double f = my::smooth_heaviside_sin(-r+tumor_radius, grid.Spacing()); // is = 1 within the tumor
          phases.phase_arrays[TISSUE](p) = 1.f-f;
          phases.phase_arrays[TCS](p) = f;
        }
      }
    }
  }
  else if (tumortype == TumorTypes::BULKTISSUE)
  {
    //there is a bulktissue tumor group provided
    H5::DataSet cell_vol_fraction_ds = tumorgroup->openDataSet("conc");
    H5::DataSet tumor_fraction_ds = tumorgroup->openDataSet("ptc");
    H5::DataSet necro_fraction_ds = tumorgroup->openDataSet("necro");
    string latticePath;
    readAttrFromH5(*tumorgroup, "LATTICE_PATH", latticePath);
    
    H5::H5File file(tumor_fraction_ds.getFileName(), H5F_ACC_RDONLY);
    H5::Group ldgroup = file.openGroup(latticePath);
    
    //H5::Group   ldgroup        = tumor_fraction_ds.get_file().root().open_group(tumor_fraction_ds.attrs().get<string>("LATTICE_PATH"));

    Array3df cell_vol_fraction;
    Array3df tumor_fraction;
    Array3df necro_fraction;
    LatticeDataQuad3d fieldld;
    ReadHdfLd(ldgroup, fieldld);
    ReadArray3D(cell_vol_fraction_ds, cell_vol_fraction);
    ReadArray3D(tumor_fraction_ds, tumor_fraction);
    ReadArray3D(necro_fraction_ds, necro_fraction);

    phases = TissuePhases(3, grid.Box());

    #pragma omp parallel
    {
      BOOST_FOREACH(const BBox3 &bbox, mtboxes.getCurrentThreadRange())
      {
        FOR_BBOX3(p, bbox)
        {
          //interpolate fields to point p on the lattice
          float this_cell_vol_fraction = FieldInterpolate::ValueAveraged(cell_vol_fraction, fieldld, FieldInterpolate::Const(0.f), grid.ld.LatticeToWorld(p));
          float this_tumor_fraction = FieldInterpolate::ValueAveraged(tumor_fraction, fieldld, FieldInterpolate::Const(0.f), grid.ld.LatticeToWorld(p));
          float this_necro_fraction = FieldInterpolate::ValueAveraged(necro_fraction, fieldld, FieldInterpolate::Const(0.f), grid.ld.LatticeToWorld(p));
    
          phases.phase_arrays[TISSUE](p) = (1.f-this_tumor_fraction)*this_cell_vol_fraction;
          phases.phase_arrays[TCS](p) = this_tumor_fraction*this_cell_vol_fraction;
          phases.phase_arrays[DEAD](p) = this_necro_fraction*this_cell_vol_fraction;
        }
      }
    }
  }
}

double EstimateTumorRadius(const ContinuumGrid& grid, const DomainDecomposition& mtboxes, const Levelset& ls)
{
  tbb::spin_mutex mutex;
  float max_rad = 0.;
  #pragma omp parallel
  {
    Float3 center = 0.5 * (grid.ld.GetWorldBox().min + grid.ld.GetWorldBox().max);
    float thread_max_rad_sqr = 0.;
    BOOST_FOREACH(const BBox3 bbox, mtboxes.getCurrentThreadRange())
    {
      FOR_BBOX3(p, bbox)
      {
        if (ls.phi(p) > 0.)
        {
          float rr = (center - grid.ld.LatticeToWorld(p)).squaredNorm();
          if (rr > thread_max_rad_sqr)
            thread_max_rad_sqr = rr;
        }
      }
    }
    tbb::spin_mutex::scoped_lock lock(mutex); // sync changes made to tumor_wbox
    max_rad = std::max(max_rad, std::sqrt(thread_max_rad_sqr));
  }
  return max_rad;
}



template<class CheckInside>
void ClassifyVesselsByRegion( const VesselList3d* list, const CheckInside &idfield, DynArray<Bitfield<uchar> > &vessid )
{
  const LatticeDataQuad3d &latData = list->Ld();
  vessid.resize( list->GetECount() );
  vessid.fill(TUMBIT_REGION_OUTSIDE);
  for( uint i=0; i<list->GetECount(); ++i )
  {
    const Vessel* v = list->GetEdge(i);
    Int3 pos(v->LPosA());
    for( uint i=0; i<v->len; ++i, pos=latData.NbLattice(pos,v->dir) )
    {
      vessid[v->Index()].AddBits( idfield(latData.LatticeToWorld(pos)) );
    }
  }
}

#if 0
template
void ClassifyVesselsByRegion<CheckInside_ContTumor3d>( const VesselList3d* list,
                                                       const CheckInside_ContTumor3d &idfield,
                                                       DynArray<Bitfield<uchar> > &vessid );

static void BuildRegionFieldQuad3d( Array3d<uchar> f, const BBox3 &tumbbox, float regwidth, bool bPeripherySeperate )
{
  // WARNING: THIS SHIT HAS CHANGED AND IS UNTESTED!!!!
  const BBox3 bbox = tumbbox;
  const BBox3 exbbox = Extend(bbox,Int3(0,0,1));
  LatticeDataQuad3d ld; ld.Init(bbox,1.0f);

  DynArray<Int3> s0(1024,ConsTags::RESERVE),s1(1024,ConsTags::RESERVE),sb(1024,ConsTags::RESERVE);
  FOR_BBOX3_BORDER(Int3,p,exbbox)
  {
    if(!ld.IsInsideLattice(p)) continue;
    if(f(p)==TUMBIT_REGION_PERIPHERY)
    {
      f(p) = TUMBIT_REGION_OUTSIDE;
      s0.push_back(p);
    }
  }

  while(!s0.empty())
  {
    for(int64 i=0; i<s0.size(); ++i)
    {
      const Int3 p = s0[i];
      for( int d=0; d<LatticeDataQuad3d::DIR_CNT; ++d )
      {
        const Int3 pp = ld.NbLattice(p,d);
        const int site = f.offset(pp);
        if( !bbox.Overlaps(pp) || f(pp)==TUMBIT_REGION_OUTSIDE ) continue;
        uchar fval = f(pp);
        if(fval==TUMBIT_REGION_CENTER)
        {
          sb.push_back(pp);
        }
        else
        {
          myAssert(fval==TUMBIT_REGION_PERIPHERY);
          s1.push_back(pp);
        }
        f(pp) = TUMBIT_REGION_OUTSIDE;
      }
    }
    swap(s0,s1);
    s1.remove_all();
  }

  clear_and_free_memory(s0);
  clear_and_free_memory(s1);

  FOR_BBOX3(p,bbox)
  {
    if(f(p)==TUMBIT_REGION_PERIPHERY) f(p)=TUMBIT_REGION_CENTER;
  }

  if(regwidth<0.5f) return;

  // around the invasive edge, mark sites as TUMBIT_REGION_BORDER within radius 'regwidth'
  for(int64 i=0; i<sb.size(); ++i )
  {
    const Int3    p0 = sb[i];
    const Float3 p0w = ld.LatticeToWorld(p0);
    int ext = my::iceil(regwidth/ld.scale);
    const BBox3 b = BBox3().Add(p0).Extend(ext).Clip(ld.Box());
    FOR_BBOX3(p,b)
    {
      const Float3 pw = ld.LatticeToWorld(p);
      if(DiffLenSqr(p0w,pw)<regwidth*regwidth)
      {
        const int site  = f.offset(p);
          if(!bPeripherySeperate || f[site]==TUMBIT_REGION_CENTER)
            f[site]=TUMBIT_REGION_PERIPHERY;
          else if(f[site]==TUMBIT_REGION_OUTSIDE)
            f[site]=TUMBIT_REGION_BORDER;
      }
    }
  }
}


CheckInside_ContTumor3d::CheckInside_ContTumor3d()
{
}

void CheckInside_ContTumor3d::Init( const ContinuumTumor &tl, float regwidth, bool bPeripherySeperate )
{
  FUNC_TIMING_START
  this->ld = tl.GetLatticeData();
  const Array3d<float> &tcs = tl.GetState().ptc;
  f.initFromBox(ld.Box(), TUMBIT_REGION_OUTSIDE);
  bbox = BBox3();
  const float cnt_thres = 0.25f;
  FOR_BBOX3(p,ld.Box())
  {
    if(tcs(p)>cnt_thres) {
      f(p) = TUMBIT_REGION_CENTER;
      bbox.Add(p);
    }
    else {
      f(p) = TUMBIT_REGION_PERIPHERY;
    }
  }
  BuildRegionFieldQuad3d(f,bbox,regwidth/ld.scale*0.5f,bPeripherySeperate);
  bbox = Clip(Extend(bbox,iceil(regwidth/ld.scale)),ld.Box());
  //{
  //  // debug
  //  Array3d<uchar> slice = Array3d<uchar>::copy(f.slice(RgNd::fromRg()()(f.shape()[2]/2)));
  //  Image img;
  //  Array3d<uchar>::Shape s = slice.shape();
  //  img.Init(s[0],s[1]);
  //  FOR_REG2V(p,0,s[0]-1,0,s[1]-1)
  //  {
  //    uchar g = slice(p);
  //    if(g==TUMBIT_REGION_OUTSIDE)
  //      img.SetPixel(p.x,p.y,128,128,128);
  //    else if(g==TUMBIT_REGION_CENTER)
  //      img.SetPixel(p.x,p.y,255,255,255);
  //    else if(g==TUMBIT_REGION_PERIPHERY)
  //      img.SetPixel(p.x,p.y,255,0,0);
  //  }
  //  img.Write("tumorinterior.png");
  //}
  FUNC_TIMING_END_OS(my::log())
}


Bitfield<uchar> CheckInside_ContTumor3d::operator()(const Float3 &pf) const
{
  Int3 p = ld.WorldToLattice(pf);
  if(!bbox.Overlaps(p)) return TUMBIT_REGION_OUTSIDE;
  else return f(p);
}
#endif

// //-----------------------------------
// //-----------------------------------
#if 0
void InitCoarseLD( LatticeDataQuad3d &ld_coarse, // result size includes border
                   const LatticeDataQuad3d &ld_fine,
                   const int border_fine,
                   const int border_coarse,
                   const int subdiv )
{
  myAssert( border_fine>=0 && border_coarse>=0 );
  Int3 s = ld_fine.Size()-Vec<int,3>(2*border_fine);
  Int3 sc;
  sc.x() = (s.x()+subdiv-1)/subdiv;
  sc.y() = (s.y()+subdiv-1)/subdiv;
  sc.z() = (s.z()+subdiv-1)/subdiv;
  myAssert(sc.x()*subdiv>=s.x());
  myAssert(sc.y()*subdiv>=s.y());
  myAssert(sc.z()*subdiv>=s.z());
  sc += Vec<int,3>(2*border_coarse);
  ld_coarse.Init(sc.x(),sc.y(),sc.z(),ld_fine.scale*subdiv);
}
#endif

//-----------------------------------
//-----------------------------------
static inline float EffectiveRadius( float r1, float r2, float l1, float l2 )
{
  return r1*r2*std::pow((l1+l2)/(my::sqr(my::sqr(r1))*l2+my::sqr(my::sqr(r2))*l1),0.25f);
}

//-----------------------------------
//-----------------------------------


int FindPositionOnVessel(const VesselList3d::LatticeData &ld, const Vessel* v, const Int3 &p )
{
//   const Int3 p1 = v->LPosA();
//   const Int3 p2 = v->LPosB();
//  int axis = v->dir/2;
//  myAssert(axis>=0);
//  myAssert(p1[(axis+1)%3]==p2[(axis+1)%3]);
//  myAssert(p1[(axis+2)%3]==p2[(axis+2)%3]);
//  int ret = std::abs(p[axis]-p1[axis]);
//  myAssert( ret<=std::abs(p2[axis]-p1[axis]) );
//  return ret;
  AxisDirLen x = ld.GetAxisDirLen(v->LPosA(), p);
  assert(x.isValid());
  return x.len;
}

static int DistanceToJunctionDirected( const VesselNode* vc, const VesselNode* last, int maxdist )
{
  myAssert( vc && last );
  int dist = 0;
  while( dist<=maxdist )
  {
    // get down edge
    const Vessel* v = vc->GetEdge(0);
    if( v->NodeA()==last || v->NodeB()==last ) {
      v = vc->GetEdge(1);
    }
    // get down node
    vc = v->GetOther(vc);
    // add up dist
    dist += v->len-1;
    // done
    if( vc->Count()>2 ) return dist;
    if( vc->Count()==1 ) return 10000;
  }
  return 10000;
}


int FindDistanceToJunction( const Vessel* vstart, int posOnVess, const VesselNode* vcstart, int maxdist )
{
  myAssert( (vstart==NULL)^(vcstart==NULL) );
  if( vcstart )
  {
    if( vcstart->Count()==2 ) {
      int d1 = DistanceToJunctionDirected( vcstart, vcstart->GetNode(0).node, maxdist );
      int d2 = DistanceToJunctionDirected( vcstart, vcstart->GetNode(1).node, maxdist );
      return std::min( d1, d2 );
    }
    else if( vcstart->Count()==1 ) {
      return DistanceToJunctionDirected( vcstart, NULL, maxdist );
    }
    else return 0;
  }
  else if( vstart )
  {
    myAssert( posOnVess>0 && posOnVess<vstart->len-1 );
    // distances to next junction at least distance to end of the vessel segment,
    // or infinite if the sprout ends there
    int d1 = (vstart->NodeA()->Count()==1) ? std::numeric_limits<int>::max() : (posOnVess);
    int d2 = (vstart->NodeB()->Count()==1) ? std::numeric_limits<int>::max() : (vstart->len-1-posOnVess);
    // if this is a section of the vesselgraph is a simple polyline,
    // then move along the line in both directions until a junction is found where more than 2 vessels meet
    if( vstart->NodeA()->Count()==2 )
      d1 += DistanceToJunctionDirected( vstart->NodeA(), vstart->NodeB(), maxdist-d1 );
    if( vstart->NodeB()->Count()==2 )
      d2 += DistanceToJunctionDirected( vstart->NodeB(), vstart->NodeA(), maxdist-d2);
    return std::min(d1,d2);
  }
  return std::numeric_limits<int>::max();
}


std::unique_ptr<VesselList3d> ReadVesselList3d(H5::Group &vesselgroup, const ptree &params)
{
  float grid_scale = params.get<float>("scale subdivide", -1.);
  bool filter_uncirculated = params.get<bool>("filter", false);
  
  
  std::unique_ptr<VesselList3d> vl;
  typedef polymorphic_latticedata::LatticeData LatticeData;
  string type_of_vessel_network;
  readAttrFromH5(vesselgroup, string("CLASS"), type_of_vessel_network);
//   // Create new string datatype for attribute
//   H5::StrType strdatatype(H5::PredType::C_S1, 256); // of length 256 characters
//   H5std_string strreadbuff("");
//   H5::Attribute myatt_out = vesselgroup.openAttribute("CLASS");
//   myatt_out.read(strdatatype,strreadbuff);
  //std::string myString = vesselgroup.openAttribute("CLASS");
  //h5cpp::Attributes a = vesselgroup.attrs();
/*
  if(vesselgroup.attrs().get<std::string>("CLASS") == "GRAPH")
    */
  if(type_of_vessel_network == "GRAPH")
  {
#ifdef DEBUG
    cout << "read vl from: \n " << type_of_vessel_network << endl;
#endif
    //we have a lattice struture->get it, could also produce an error
        /*
     * Access "lattice" group .
     */
    H5::Group ldgroup;
    ldgroup = vesselgroup.openGroup("lattice");
//     try {  // to determine if the lattice group exists
//       ldgroup = vesselgroup.openGroup("lattice");
//     }
//     catch( H5::GroupIException not_found_error ) {
//         cout << " lattice group is not found." << endl;
//     }
    cout << "group lattice is open" << endl;
    
//     if(!vesselgroup.exists("lattice"))
//     {
//       string latticeIOError = str(format("No lattice group in hdf: %s\n") % vesselgroup.get_file_name());
//       throw std::runtime_error(latticeIOError);
//     }
//     h5cpp::Group ldgroup = vesselgroup.open_group("lattice");//may not be there

    
    std::unique_ptr<polymorphic_latticedata::LatticeData> ldp = polymorphic_latticedata::ReadHdf(ldgroup);
#ifdef DEBUG
    cout << "ReadVesselList3d read " << endl;
    ldp->print(cout);
#endif
    //std::unique_ptr<VesselList3d> vl_local;
    //vl_local.reset(new VesselList3d(ldp));
    //vl_local->Init(*ldp);
    vl = std::unique_ptr<VesselList3d> (new VesselList3d());
    vl->Init(ldp);
  
    
#ifdef DEBUG
    VESSEL_INTEGRITY_CHECK_SWITCH(vl->IntegrityCheck();)
#endif

    ReadHdfGraph(vesselgroup, vl.get());
    //this magic can only be done on a lattice
    float original_grid_scale_override = params.get<float>("scale override", -1.);
    if (original_grid_scale_override>0.)
    {
      throw std::runtime_error("vessel scale override not implemented");
      //ldp->Scale() = original_grid_scale_override;
    }
    // subdivide the grid to other lattice spacing
    if (grid_scale > 0)
    {
      // ldp is now bound to vl,  use vl->Ld().Scale() instead
      double scale_of_vessel_data = vl->Ld().Scale();
      // check whether subdivision makes sense
      myAssert(scale_of_vessel_data/grid_scale - int(scale_of_vessel_data/grid_scale) < 1.e-3  && scale_of_vessel_data/grid_scale > 1.);
      //vl = std::move(GetSubdivided( vl, grid_scale));
      GetSubdivided( vl, grid_scale);

    #ifdef DEBUG
      VESSEL_INTEGRITY_CHECK_SWITCH(vl->IntegrityCheck();)
    #endif
    }
  }
  else
  {
    //null!!!! 
    //world no lattice data needed
    std::unique_ptr<LatticeData> ldp;
    std::unique_ptr<VesselList3d> vl_local(new VesselList3d());
    vl_local->Init(ldp);
    vl=std::move(vl_local);
    ReadHdfGraph(vesselgroup, vl.get());
  }

#if 1
#ifdef DEBUG
#ifndef TOTAL_SILENCE
//typedef boost::unordered_map<int, FlowBC> BcsMap;
    for(auto bc: vl->GetBCMap())
    {
      printf("first: %i, second: %f\n", bc.first->Index(), bc.second.val);
    }
#endif
#endif
#endif
  
  
  {//hdf failscope
    int ecnt = vl->GetECount();
    int vcnt = vl->GetNCount();
    
    for( int i=0; i<ecnt; ++i )
    {
      Vessel* v = vl->GetEdge(i);
      v->maturation = GetInitialThickness( v->r );
      v->hematocrit = -1; // since i use hematocrit as a system parameter now i want things to fail hard if this is not properly set somewhere else
      v->f = v->q = 0.;
    }

#if 0
    DynArray<decltype(VNodeData::press)> pressure;
    if (vesselgroup.exists("nodes/pressure"))
    {
      h5::read_dataset(vesselgroup.open_dataset("nodes/pressure"),pressure);
    }

    if (vesselgroup.exists("nodes/roots_pressure"))
    {
      pressure.resize(vcnt);
      DynArray<decltype(VNodeData::press)> roots_pressure;
      h5::read_dataset(vesselgroup.open_dataset("nodes/roots_pressure"),roots_pressure);
      DynArray<int> idx;
      h5::read_dataset(vesselgroup.open_dataset("nodes/roots"),idx);
      for (int i=0; i<idx.size(); ++i)
      {
        pressure[idx[i]] = roots_pressure[i];
      }
    }

    DynArray<decltype(VData::q)> flow;
    if (vesselgroup.exists("edges/flow"))
    {
      h5::read_dataset(vesselgroup.open_dataset("edges/flow"),flow);
    }
    
    DynArray<decltype(VData::f)> shearforce;
    if (vesselgroup.exists("edges/shearforce"))
    {
      h5::read_dataset(vesselgroup.open_dataset("edges/shearforce"), shearforce);
    }
    
    DynArray<int> flags;
    if (vesselgroup.exists("edges/flags"))
    {
      h5::read_dataset(vesselgroup.open_dataset("edges/flags"), flags);
    }

    DynArray<decltype(VData::hematocrit)> hematocrit; // however, read hematocrit if available
    if (vesselgroup.exists("edges/hematocrit"))
    {
      h5::read_dataset(vesselgroup.open_dataset("edges/hematocrit"),hematocrit);
    }
    //set back the read values
    for(int i=0;i<vl->GetECount();++i)
    {
      Vessel* v=vl->GetEdge(i);
      //I don't know why this is there, and therefore I keep it
      v->NodeA()->flags.AddBits( v->flags.GetBits(ARTERY|VEIN) );
      if (!flow.empty())
        v->q = flow[i];
      if (!hematocrit.empty())
        v->hematocrit = hematocrit[i];
      if (!shearforce.empty())
	v->f = shearforce[i];
      if (!flags.empty())
	v->flags = flags[i];
    }
    //also set back the boundary conditions, if found
    if(!vl->GetBCMap().empty())
    {
      //if we have the BCMap, we try to use it
      //ReadHdfGraph executed before should provide that
      for( auto bc: vl->GetBCMap())
      {
	const VesselNode* nodeFromBCMap = bc.first;
	FlowBC theBC = bc.second;
	switch(theBC.typeOfInstance)
	{
	  case FlowBC::PIN:
	    vl->GetNode(nodeFromBCMap->Index())->press=theBC.val;
	    break;
	  case FlowBC::CURRENT:
	    const Vessel* theEdgeToTheNode=nodeFromBCMap->GetEdge(0);
	    vl->GetEdge(theEdgeToTheNode->Index())->q = theBC.val;
	    //cout<<"not implemented yet"<<endl;
	    break;
	}
      }
    }
    else
    {
      cout<<"Warning vl->GetBCMap().empty() "<<endl;
      //go through complete list
      for(int i=0; i<vcnt; ++i)
      {
	VesselNode* vc = vl->GetNode(i);
	if (vc->flags.GetBits(BOUNDARY) && pressure.empty())
	{
	//guess pressure, so there is an reasonable value to compute pressure
	  vc->press = PressureRadiusRelation(vc->GetEdge(0)->r,(vc->GetEdge(0)->IsArtery()));
	}
	else if (vc->flags.GetBits(BOUNDARY) && !pressure.empty())
	{
	  vc->press = pressure[vc->Index()];
	}
      }
    }
#endif
  }//end hdf failscope
    

  if (filter_uncirculated)
  {
    { DynArray<Vessel*> tokill(100, ConsTags::RESERVE);
    for (int i=0; i<vl->GetECount(); ++i)
    {
      Vessel* v = vl->GetEdge(i);
      if (!v->IsCirculated()) tokill.push_back(v);
    }
    for (int i=0; i<tokill.size(); ++i) vl->DeleteVessel(tokill[i]);
    }

    { DynArray<VesselNode*> tokill(100, ConsTags::RESERVE);
    for (int i=0; i<vl->GetNCount(); ++i)
    {
      VesselNode* nd = vl->GetNode(i);
      if (nd->Count() == 0) tokill.push_back(nd);
    }
    for (int i=0; i<tokill.size(); ++i) vl->DeleteUnusedNode(tokill[i]);
    }
  }

  VESSEL_INTEGRITY_CHECK_SWITCH(vl->IntegrityCheck();)

  const auto calcflow = params.get_child_optional("calcflow");
  if (calcflow)
  {
    BloodFlowParameters params;
    params.assign(*calcflow);
    ComputeCirculatedComponents(vl.get());
    CalcFlow(*vl, params);
  }
#if 1
#ifdef DEBUG
#ifndef TOTAL_SILENCE
    for(auto bc: vl->GetBCMap())
    {
      printf("first: %i, second: %f\n", bc.first->Index(), bc.second.val);
    }
#endif
#endif
#endif

  return vl;
}


void WriteVesselList3d(const VesselList3d &vl, H5::Group &vesselgroup, const ptree &params)
{
  WriteHdfGraph(vesselgroup,vl);

  
  int ecnt = vl.GetECount();
  int ncnt = vl.GetNCount();
  DynArray<double> arrd;
  DynArray<float> arrf;

  arrd.resize(ncnt);

  bool w_all      = params.get<bool>("w_all", true);
  bool w_pressure = w_all || params.get<bool>("w_pressure", true);
  bool w_flow     = w_all || params.get<bool>("w_flow", true);
  bool w_adaption = w_all || params.get<bool>("w_adaption", true);

  
  if (w_pressure)
  {
    for(int i=0; i<ncnt; ++i) { arrd[i] = vl.GetNode(i)->press; }
    //h5::Dataset ds = h5::create_dataset(vesselgroup,"nodes/pressure",arrd);
    //H5::DataSet ds = vesselgroup.createDataSet("nodes/pressure");
    writeDataSetToGroup(vesselgroup, string("nodes/pressure"), arrd);
    //H5::DataSet ds = h5::create_dataset(vesselgroup,"nodes/pressure",arrd);
    //ds.attrs().set("MODE","linear");
    H5::DataSet press;
    try
    {
      press = vesselgroup.openDataSet(string("nodes/pressure"));
    }
    catch(H5::Exception e)
    {
      e.printError();
    }
    writeAttrToH5(press,string("MODE"), string("linear"));

    //for(int i=0; i<ncnt; ++i) { arrd[i] = vl.GetNode(i)->residual; }
    //ds = h5::create_dataset_range(vesselgroup,"nodes/residual",arrd);
  }
// #define _WRITE_STATE_FILL_ARRAY_EDGE(arr,var,name)\
// {\
//   for(int i=0; i<ecnt; ++i) { arr[i] = vl.GetEdge(i)->var; }\
//   h5::Dataset ds = h5::create_dataset(vesselgroup,name,arr);\
//   ds.attrs().set("MODE","const");\
// }
#define _WRITE_STATE_FILL_ARRAY_EDGE(arr,var,name)\
{\
  for(int i=0; i<ecnt; ++i) { arr[i] = vl.GetEdge(i)->var; }\
  H5::DataSet ds = writeDataSetToGroup(vesselgroup, string(name), arr);\
  writeAttrToH5(ds, string("MODE"), string("const"));\
}
  arrf.resize(ecnt);
  arrd.resize(ecnt);

  
  if( w_adaption )
  {
    _WRITE_STATE_FILL_ARRAY_EDGE(arrd,conductivitySignal,"edges/conductivitySignal");
    _WRITE_STATE_FILL_ARRAY_EDGE(arrd,metabolicSignal,"edges/metabolicSignal");
    _WRITE_STATE_FILL_ARRAY_EDGE(arrd,S_total,"edges/S_tot");
  }
  if (w_flow)
  {
    _WRITE_STATE_FILL_ARRAY_EDGE(arrd,f,"edges/shearforce");
    _WRITE_STATE_FILL_ARRAY_EDGE(arrd,q,"edges/flow");
    _WRITE_STATE_FILL_ARRAY_EDGE(arrd,hematocrit,"edges/hematocrit");
  }
  if (w_all)
  {
    _WRITE_STATE_FILL_ARRAY_EDGE(arrd,initial_f,"edges/initial_shearforce");
    _WRITE_STATE_FILL_ARRAY_EDGE(arrd,reference_r,"edges/initial_radius");
    _WRITE_STATE_FILL_ARRAY_EDGE(arrd,maturation,"edges/maturation");
    _WRITE_STATE_FILL_ARRAY_EDGE(arrd,timeInTumor,"edges/timeInTumor");
    _WRITE_STATE_FILL_ARRAY_EDGE(arrd,timeSprout,"edges/timeSprout");
  }
  DynArray<uchar> tmp2;
  tmp2.resize(vl.GetNCount());
  //DynArray<uchar> tmp2(vl.GetNCount());
  for (int i=0; i<vl.GetNCount(); ++i)
  {
    tmp2[i] = vl.GetNode(i)->flags;
  }
  H5::DataSet nodeflags = writeDataSetToGroup(vesselgroup, string("nodes/nodeflags"), tmp2);
  //h5::create_dataset(vesselgroup, "nodes/nodeflags", tmp2);
}


//-----------------------------------
//-----------------------------------

VesselVolumeGenerator::VesselVolumeGenerator(const VesselList3d &vl_, const LatticeDataQuad3d &fieldld_, int dim_, const ptree &params)
  : vl(vl_), ld(vl_.Ld()), fieldld(fieldld_), dim(dim_)
{
  use_smooth_delta = params.get<bool>("smooth", true);
  upper_volume_bound = params.get<float>("cut_at", 1.f);
  sampler.Init(fieldld.Scale(), params);
}



void VesselVolumeGenerator::Fill(const BBox3& bbox, Array3df vessel_volume, const DynArray< const Vessel* > vessels, int& box_sample_cnt)
{
  for (int i=0; i<vessels.size(); ++i)
  {
    //get vessel
    const Vessel* v = vessels[i];
    //get real world coordinate of starting and ending point
    my::eqpair<Float3> edge(ld.LatticeToWorld(v->LPosA()),
                            ld.LatticeToWorld(v->LPosB()));
    //get segment radius
    float rad = v->r;

    sampler.Set(&edge, &rad, 1);
    int cnt = sampler.GenerateVolumeSamples();

    for (int j=0; j<cnt; ++j)
    {
      CylinderNetworkSampler::Sample smpl = sampler.GetSample(j);
      if (use_smooth_delta) {
        AddSmoothDelta(vessel_volume, bbox, fieldld, dim, smpl.wpos, sampler.weight_per_volume);
      }
      else
      {
        Int3 lpos = fieldld.WorldToLattice(smpl.wpos);
        if (bbox.Overlaps(lpos))
          vessel_volume(lpos) += sampler.weight_per_volume;
      }
    }
    box_sample_cnt += cnt;
  }

  if (upper_volume_bound > 0.)
  {
    // unfortunately large vessel can produce values larger than 1 because
    // several segments will overlap within a cell
    FOR_BBOX3(p, bbox)
    {
      float &v = vessel_volume(p);
      v = std::min(v, upper_volume_bound);
    }
  }
}


/*
 * parameters:
 *  smooth_ : bool (true) - use smothended delta function with 4 cell support width
 *  cut_at : float (1.) - limit resulting values to this upper bound. The limiting is applied to the values in the supplied vessel_volume array, use < 0 to disable
 *  samples_per_cell : int (2) - the number of samples (approximately) taken within one cell volume of the array lattice
 *  seed : uint (12345) - random seed
 */
ptree AddVesselVolumeFraction(const VesselList3d &vl, const LatticeDataQuad3d &field_ld, int dim, Array3d<float> vessel_volume, const ptree &params)
{
  my::Time t__;

  DynArray<BBox3> boxes = MakeMtBoxGrid(field_ld.Box(), Int3(1<<31, 32, 32));
  DynArray< DynArray<const Vessel*> > grid;
  SortVesselsIntoMtBoxGrid(field_ld, vl, 2, boxes, grid); // this is single threaded and O(n^2) in the number of boxes! -> must be replaced

  volatile int total_num_samples = 0;

  cout << format("sort time=%f") % (my::Time() - t__).to_ms() << endl;
  
  #pragma omp parallel
  {
    int thread_sample_cnt = 0;
    VesselVolumeGenerator volumegen(vl, field_ld, dim, params);

    #pragma omp for schedule(dynamic, 1)
    for (int i=0; i<grid.size(); ++i)
    {
      volumegen.Fill(boxes[i], vessel_volume, grid[i], thread_sample_cnt);
    }

    #pragma omp atomic
    total_num_samples += thread_sample_cnt;
  }
  cout << format("add vessel volume fraction: time=%f, np=%i, npboxes=%i") % (my::Time()-t__).to_ms() % omp_get_max_threads() % boxes.size() << endl;
  return make_ptree("num_samples", total_num_samples);
}

//-----------------------------------
//-----------------------------------



static void CalcGoodSpacing(Float3 size, float &spacing, Int3 &num_cells)
{
#if 0
  Int3 pow2;
  for (int i=0; i<3; ++i)
  {
    //ok ??  maybe?!
    pow2[i] = 2;  // for chombo, the field size must be divisible by 2
    while (true) {
      float q = size[i]/(pow2[i]*spacing);
      if (q < 20.) break;
      pow2[i] *= 2;
    }
    int n = std::max(1, int(size[i]/(pow2[i]*spacing)));
    num_cells[i] = n * pow2[i];
//     num_cells[i] = std::max(2, int(size[i]/spacing));
//     num_cells[i] &= ~1;
  }
  spacing = Max(size.cwise() / num_cells.cast<float>());
#endif
  Float3 h(spacing);
  for (int i=0; i<3; ++i)
  {
    num_cells[i] = std::max<int>(1, my::round(size[i]/h[i]));
  }
}


void CenterVesselListAndSetupFieldLattice(VesselList3d &vl, int dim, float spacing, LatticeDataQuad3d &field_ld)
{
  const FloatBBox3 wbb = vl.Ld().GetWorldBox();
  const Float3 origin = vl.Ld().LatticeToWorld(Int3(0,0,0));
  const Float3 size = ::Size(wbb);
  Float3 org(0.);
  for (int i=0; i<dim; ++i) org[i] = -0.5*size[i]+(origin[i]-wbb.min[i]);
  vl.SetDomainOrigin(org);

  const float safety = 0.1*vl.Ld().Scale();
  SetupFieldLattice(vl.Ld().GetWorldBox(), dim, spacing, safety, field_ld);
}


void SetupFieldLattice(const FloatBBox3 &wbbox, int dim, float spacing, float safety_spacing, LatticeDataQuad3d &ld)
{
  Float3 vessel_domain_size = wbbox.max - wbbox.min;
  Float3 domain_size(0.); Bool3 cell_centering(false);
  Int3 num_cells(1);
  Float3 world_offset(0.);
  for (int i=0; i<dim; ++i) {
    domain_size[i] = vessel_domain_size[i] + 2.*safety_spacing;
    cell_centering[i] = true;
    num_cells[i] = std::max<int>(1, my::iceil(domain_size[i]/spacing));
    world_offset[i] = -num_cells[i]*spacing*0.5 + vessel_domain_size[i]*0.5 + wbbox.min[i];
  }
  ld.Init(num_cells, spacing);
  ld.SetOriginPosition(world_offset);
  ld.SetCellCentering(cell_centering);
}


//-----------------------------------
//-----------------------------------
namespace FieldInterpolate
{

static Int3 projectToWithinLattice(const LatticeDataQuad3d &ld, const Int3 &p_)
{
  /*
   * if p is outside of the lattice, then project it to the boundary
   * else it is the identity
   */
  const BBox3 bb = ld.Box();
  Int3 p(p_);
  for (int ax=0; ax<3; ++ax)
  {
    if (p[ax]<bb.min[ax]) p[ax] = bb.min[ax];
    if (p[ax]>bb.max[ax]) p[ax] = bb.max[ax];
  }
  return p;
}

  
template<class T>
T inline get(const ConstArray3d<T>& ff, const LatticeDataQuad3d &ld, const OutOfDomainExtrapolateT &ood, const Int3& ip)
{
  if (!ld.IsInsideLattice(ip)) return ff(projectToWithinLattice(ld,ip));
  return ff(ip);
}

template<class T>
T inline get(const ConstArray3d<T>& ff, const LatticeDataQuad3d &ld, const OutOfDomainConstT<T> &ood, const Int3& ip)
{
    if (!ld.IsInsideLattice(ip)) return ood.m_value;
    return ff(ip);
}


template<class T, class OOD, int dim>
struct Avg
{
  static inline T Quad(const ConstArray3d<T> &ff, const LatticeDataQuad3d &ld, const OOD &ood, const Int3 &pos)
  {
    // average recursively along an axis, use zero-gradient boundary condition
    T values = 0.;
    values += Avg<T, OOD, dim-1>::Quad(ff, ld, ood, pos);
    Int3 pos_nb = pos;
    ++pos_nb[dim];
    values += Avg<T, OOD, dim-1>::Quad(ff, ld, ood, pos_nb);
    return 0.5 * values;
  }
  static inline T Grad(const ConstArray3d<T> &ff, const LatticeDataQuad3d &ld, const OOD &ood, const Int3 &pos, int grad_dim)
  {
    T values = 0;
    values = (dim == grad_dim ? -1. : 1.) * Avg<T, OOD, dim-1>::Grad(ff, ld, ood, pos, grad_dim);
    Int3 pos_nb = pos; ++pos_nb[dim];
    values += Avg<T, OOD, dim-1>::Grad(ff, ld, ood, pos_nb, grad_dim);
    if (dim == grad_dim)
      return values * ld.ScaleInv();
    else
      return 0.5f * values;
  }
};

template<class T, class OOD>
struct Avg<T, OOD, -1>
{
  static inline T Quad(const ConstArray3d<T> &ff, const LatticeDataQuad3d &ld, const OOD &ood, const Int3 &pos)
  {
    return get(ff, ld, ood, pos);
  }

  static inline T Grad(const ConstArray3d<T> &ff, const LatticeDataQuad3d &ld, const OOD &ood, const Int3 &pos, int grad_dim)
  {
    return get(ff, ld, ood, pos);
  }
};



template<class T, class OOD>
T ValueAveraged(const ConstArray3d<T>& ff, const LatticeDataQuad3d &ld, const OOD &ood, const Float3& pos)
{
  Int3 ip; Float3 q;
  boost::tie(ip, q) = ld.WorldToLatticeCell(pos);

  T f[2][2][2];
  for (int x=0; x<2; ++x)
    for (int y=0; y<2; ++y)
      for (int z=0; z<2; ++z)
      {
        f[x][y][z] = Avg<T, OOD, 2>::Quad(ff, ld, ood, ip+Int3(x-1,y-1,z-1));
      }
  return Lerp( q[2], Lerp2D( q[0],q[1], f[0][0][0], f[1][0][0], f[0][1][0], f[1][1][0] ),
                    Lerp2D( q[0],q[1], f[0][0][1], f[1][0][1], f[0][1][1], f[1][1][1] ));
}


template<class T, class OOD>
Vec<T,3> Gradient(const ConstArray3d<T>& ff, const LatticeDataQuad3d &ld, const OOD &ood, const Float3 &pos)
{
  Int3 ip; Float3 q;
  Vec<T,3> r;
  boost::tie(ip, q) = ld.WorldToLatticeCell(pos);

  for (int axis=0; axis<3; ++axis)
  {
    T f[2][2][2]; // has the gradient here
    for (int x=0; x<2; ++x)
      for (int y=0; y<2; ++y)
        for (int z=0; z<2; ++z)
        {
          f[x][y][z] = Avg<T, OOD, 2>::Grad(ff, ld, ood, ip+Int3(x-1,y-1,z-1), axis);
        }
    r[axis] = Lerp( q[2], Lerp2D( q[0],q[1], f[0][0][0], f[1][0][0], f[0][1][0], f[1][1][0] ),
                    Lerp2D( q[0],q[1], f[0][0][1], f[1][0][1], f[0][1][1], f[1][1][1] ));
  }
  return r;
}


template<class T, class OOD>
T Value(const ConstArray3d<T>& ff, const LatticeDataQuad3d &ld, const OOD &ood, const Int3& ip)
{
  return get(ff, ld, ood, ip);
}


template<class T, class OOD>
T Value(const ConstArray3d<T>& ff, const LatticeDataQuad3d &ld, const OOD &ood, const Float3& pos)
{
  Int3 ip; Float3 q;
  boost::tie(ip, q) = ld.WorldToFractionalCoordinate(pos);

  T f[2][2][2];
  for (int x=0; x<2; ++x)
    for (int y=0; y<2; ++y)
      for (int z=0; z<2; ++z)
      {
        f[x][y][z] = get(ff, ld, ood, ip+Int3(x,y,z));
      }
  return Lerp( q[2], Lerp2D( q[0],q[1], f[0][0][0], f[1][0][0], f[0][1][0], f[1][1][0] ),
                    Lerp2D( q[0],q[1], f[0][0][1], f[1][0][1], f[0][1][1], f[1][1][1] ));
}


#define INSTANTIATE_INTERPOLATE_FIELD(T, OOD)\
  template Vec<T,3> Gradient(const ConstArray3d<T>& ff, const LatticeDataQuad3d &ld, const OOD &ood, const Float3 &pos); \
  template T ValueAveraged(const ConstArray3d<T>& ff, const LatticeDataQuad3d &ld, const OOD &ood, const Float3& pos); \
  template T Value(const ConstArray3d<T>& ff, const LatticeDataQuad3d &ld, const OOD &ood, const Float3& pos); \
  template T Value(const ConstArray3d<T>& ff, const LatticeDataQuad3d &ld, const OOD &ood, const Int3& pos);

INSTANTIATE_INTERPOLATE_FIELD(float, OutOfDomainConstT<float>)
INSTANTIATE_INTERPOLATE_FIELD(float, OutOfDomainExtrapolateT)
INSTANTIATE_INTERPOLATE_FIELD(double, OutOfDomainConstT<double>)
INSTANTIATE_INTERPOLATE_FIELD(double, OutOfDomainExtrapolateT)

}

//-----------------------------------
//-----------------------------------

template<class Container>
static void SortVesselsIntoMtBoxGrid_(const LatticeDataQuad3d &ld, const VesselList3d &vl, int boxborder, const Container &boxes, DynArray< DynArray<const Vessel*> > &grid)
{
  my::Time t_;
  if(boxes.size()<=0)
    throw std::runtime_error("SortVesselsIntoMtBoxGrid got empty mtboxes list");

  DynArray<FloatBBox3> wboxes(boxes.size());
  for (int i=0; i<boxes.size(); ++i)
  {
    BBox3 bb = Extend(boxes[i], boxborder);
    wboxes[i] = FloatBBox3(ld.LatticeToWorld(bb.min) - Float3(0.5*ld.Scale()),
                           ld.LatticeToWorld(bb.max) + Float3(0.5*ld.Scale()));
  }

  grid.resize(boxes.size());
  for (int vi=0; vi<vl.GetECount(); ++vi)
  {
    const Vessel* v = vl.GetEdge(vi);
    FloatBBox3 vwbb = FloatBBox3().Add(vl.Ld().LatticeToWorld(v->LPosA())).Add(vl.Ld().LatticeToWorld(v->LPosB()));
    vwbb.Extend(v->r);
    for (int i=0; i<boxes.size(); ++i)
    {
      if (!(vwbb.Overlaps(wboxes[i]))) continue;
      grid[i].push_back(v);
    }
  }
  cout << format("sort vessel into mt boxes: %s") % (my::Time()-t_) << endl;

#if 0
  int total_in_grids = 0;
  for (int i=0; i<boxes.size(); ++i)
  {
    total_in_grids += grid[i].size();
    cout << format("%i has %s with %i vessels") % i % boxes[i] % grid[i].size() << endl;
  }
  cout << format("overhead: %f%% of %i vessels") % (float(total_in_grids)/vl.GetECount()) % vl.GetECount() << endl;
#endif
}

void SortVesselsIntoMtBoxGrid(const LatticeDataQuad3d &ld, const VesselList3d &vl, int boxborder, const DynArray<BBox3> &boxes, DynArray< DynArray<const Vessel*> > &grid)
{
  SortVesselsIntoMtBoxGrid_(ld, vl, boxborder, boxes, grid);
}

void SortVesselsIntoMtBoxGrid(const LatticeDataQuad3d &ld, const VesselList3d &vl, int boxborder, const DomainDecomposition &mtboxes, VesselsInBoxes &grid)
{
  SortVesselsIntoMtBoxGrid_(ld, vl, boxborder, mtboxes, grid);
}



//-----------------------------------
//-----------------------------------


void CylinderNetworkSampler::Init(double spacing_, const ptree &params)
{
  spacing = spacing_;
  seed = params.get<uint>("seed", 12345);
  samples_per_cell = params.get<int>("samples_per_cell", 2);
  rnd.Init(seed);
  buffer.reserve(1024*16);
  Set(NULL,NULL,0);
  Restart();
}


void CylinderNetworkSampler::Set(const WorldEdge *edges_, const float *radii_, int count_)
{
  edges = edges_;
  radii = radii_;
  edge_count = count_;
  index = 0;
}

void CylinderNetworkSampler::Set(const Float3 &posa, const Float3 &posb, float rad)
{
  this->single_edge.first = posa;
  this->single_edge.second = posb;
  this->single_rad = rad;
  Set(&single_edge, &single_rad, 1);
}

void CylinderNetworkSampler::Restart()
{
  index = 0;
}


int CylinderNetworkSampler::GenerateLineSamples()
{
  myAssert(spacing>0 && samples_per_cell>0);
  buffer.remove_all();
  if (index >= edge_count) return 0;

  float rad = radii[index];
  Float3 posa = edges[index].first,
         posb = edges[index].second;

  float len = norm(posa-posb);
  float q = len/spacing;
  int num_samples = int(q);
  q -= num_samples;
  if (q > rnd.Get01f())
    ++num_samples;

  if (num_samples > 0)
  {
    weight = len / num_samples;
    weight_per_volume = weight / my::cubed(spacing);

    HaltonSequence<1> qrnd;

    Sample s;
    for (int i=0; i<num_samples; ++i)
    {
      float r = GetJittered01f<1>(qrnd, rnd, num_samples)[0];
      float f = (i + 0.25 + 0.5*r)/num_samples;
      s.fpos[2] = f;
      s.wpos = (1-f)*posa + f*posb;
      buffer.push_back(s);
    }
  }
  ++index;
  return buffer.size();
}


int CylinderNetworkSampler::GenerateSurfaceSamples()
{
  myAssert(spacing>0 && samples_per_cell>0);
  buffer.remove_all();
  if (index >= edge_count) return 0;

  float rad = radii[index];
  Float3 posa = edges[index].first,
         posb = edges[index].second;

  float len2 = squaredNorm(posb-posa);
  float len = std::sqrt(len2);
  float area = len*rad*mconst::fpi2();
  float q = samples_per_cell*(area/sqr(spacing));
  int num_samples = int(q);
  q -= num_samples;
  if (q > rnd.Get01f())
    ++num_samples;
  
  if (num_samples > 0)
  {
    weight = area/num_samples;
    weight_per_volume = weight / my::cubed(spacing);

    HaltonSequence<2> qrnd;

    Float3x4 m = OrthogonalSystem(posb-posa);
    m.X *= rad;
    m.Y *= rad;
    m.Z *= len;
    m.T = posa;

    Sample s;

    for (int i=0; i<num_samples; ++i)
    {
    // volume sample of cylinder
      Float2 rv = GetJittered01f<2>(qrnd, rnd, num_samples);
      float u = mconst::fpi2()*rv[0];
      Float3 p(cosf(u),sinf(u),rv[1]);
      s.fpos = p;
      s.wpos = m*p;
      buffer.push_back(s);
    }
  }

  ++index;
  return buffer.size();
}



int CylinderNetworkSampler::GenerateVolumeSamples()
{
  myAssert(spacing>0 && samples_per_cell>0);
  buffer.remove_all();
  if (index >= edge_count) return 0;

  float rad = radii[index];
  Float3 posa = edges[index].first,
         posb = edges[index].second;

  float cellvol = my::cubed(spacing);
  float len = norm(posb-posa);
  float volume = my::sqr(rad)*mconst::fpi() * len; // added
  float q;
  if (rad < spacing) // 1d
     q = samples_per_cell*(len/spacing);
  else
     q = samples_per_cell*(volume/cellvol); // changed
  int num_samples = int(q);
  q -= num_samples;
  if (q > rnd.Get01f())
    ++num_samples;

  if (num_samples > 0)
  {
    weight = volume/num_samples;
    weight_per_volume = weight / my::cubed(spacing);

    HaltonSequence<3> qrnd;

    Float3x4 m = OrthogonalSystem(posb-posa);
    m.T = posa; // changed

    Sample s;
    for (int i=0; i<num_samples; ++i)
    {
      // volume sample of cylinder
      Float3 rv = GetJittered01f<3>(qrnd, rnd, num_samples);
      // transform to cylinder
      float u = mconst::fpi2()*rv[0];
      float r = std::sqrt(rv[1]);
      float t = rv[2];

      Float3 f(r*std::cos(u),r*std::sin(u),t);
      Float3 p(rad*f[0],rad*f[1],len*f[2]);

      s.fpos = f;
      s.wpos = m*p;

      buffer.push_back(s);
    }
  }

  ++index;
  return buffer.size();
}

void Set(CylinderNetworkSampler &sampler, const Vessel* v, const VesselList3d::LatticeData &ld)
{
  sampler.Set(ld.LatticeToWorld(v->LPosA()),
              ld.LatticeToWorld(v->LPosB()),
              v->r);
}





boost::optional<ptree> HandleSimulationProgramArguments(const ptree &default_params, int argc, char **argv, const ptree &vararg_params)
{
  namespace po = boost::program_options;
  po::options_description desc("Allowed Options");
  desc.add_options()
    ("help", "produce help message")
    ("parameter-file,p", po::value<string>(), "file name for a boost .info file which contains parameters (use 'show' to see the default parameters)")
    ("parameter-string,s", po::value<string>(), "parameters supplied by a string the content of which represents a boost .info file");
  po::variables_map vm;
  try
  {
    po::store(po::parse_command_line(argc, argv, desc), vm);
  }
  catch (const po::error &e)
  {
    cout << e.what() << endl;
    cout << desc << "\n";
    return boost::optional<ptree>();
  }
  
  po::notify(vm);

  if (vm.count("help")) {
      cout << desc << "\n";
      return boost::optional<ptree>();
  }

  ptree read_params;
  if (vm.count("parameter-file"))
  {
    string filename = vm["parameter-file"].as<string>();
    if (filename == "show")
    {
      cout << "parameters:" << endl;
      boost::property_tree::write_info(cout, default_params);
      return boost::optional<ptree>();
    }
    else if (filename == "stdin")
    {
      boost::property_tree::read_info(std::cin, read_params);
    }
    else
    {
      std::ifstream f(filename.c_str());
      boost::property_tree::read_info(f, read_params);
    }
  }

  if (vm.count("parameter-string"))
  {
    string parameter_file_string = vm["parameter-string"].as<string>();
    std::istringstream f(parameter_file_string);
    ptree p;
    boost::property_tree::read_info(f, p);
    boost::property_tree::update(read_params, p);
  }

  ptree unknown_params = boost::property_tree::subtract(read_params, default_params);
  unknown_params = boost::property_tree::remove(unknown_params, vararg_params);

  if (unknown_params.begin() != unknown_params.end())
  {
    std::cerr << "--- WARNING: Unknown input parameters detected ---" << endl;
    boost::property_tree::write_info(std::cerr, unknown_params);
    std::cerr << "-------------------------------------------------" << endl;
    return boost::optional<ptree>();
  }

  return read_params;
}
void isVesselListGood(VesselList3d &vl)
{
  for ( auto bc : vl.GetBCMap())
  {
    const VesselNode* vc = bc.first;
    if (vc->press<=0)
    {
#ifndef NDEBUG
      cout<<format("borders of vessellist not so nice!\n");
#endif
      exit;
    }
  }
}

void printPtree(ptree const &pt)
{
    using boost::property_tree::ptree;
    ptree::const_iterator end = pt.end();
    for (ptree::const_iterator it = pt.begin(); it != end; ++it) {
        std::cout << it->first << ": " << it->second.get_value<std::string>() << std::endl;
        printPtree(it->second);
    }
}

