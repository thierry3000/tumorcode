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

#include <fstream>
#include <stdexcept>

#include <boost/foreach.hpp>
#include <boost/property_tree/info_parser.hpp>
#include <boost/optional.hpp>

#include "../hdfio.h"
#include "remodeler.h"
#include "mwlib/helpers-vec.h"
#include "mwlib/any_tree.h"

using boost::str;
using boost::format;


void WriteHdfHistogram( H5::Group g, const string &id, const BasicHistogram1D<float> &h )
{
  Array3d<float> a;
  bool ue = h.UseError();
  a.init(Int3(ue ? 3 : 2, h.info.num_buckets, 1));
  for(int i=0; i<h.info.num_buckets; ++i)
  {
    a(0,i,0) = h.GetXC(i);
    float f = h.GetAvg(i);
    //if(!IsFinite(f)) f = -1.;
    a(1,i,0) = f;
    if(ue) {
      f = h.GetVar(i);
      //if(!IsFinite(f)) f = -1.;
      a(2,i,0) = sqrt( std::max<float>( 0, f ) );
    }
  }
  H5::DataSet ds = WriteArray3D(g,id,a);
  //h5::Dataset ds = WriteArray3D(f,id,a);
  writeAttrToH5(ds, string("TITLE"), h.name);
  //writeAttrToH5<decltype(h.info.size_bucket)>(ds, "BUCKED_SIZE", h.info.size_bucket );
  writeAttrToH5(ds, "BUCKED_SIZE", h.info.size_bucket );
  writeAttrToH5(ds, "FIELD_0_NAME", h.info.name);
  writeAttrToH5(ds,"FIELD_1_NAME", h.name );
  if( ue)
  {
    writeAttrToH5(ds,"FIELD_2_NAME",string("rmse") );
  }
  //h5::Attributes attrs = ds.attrs();
//   attrs.set("TITLE",h.name);
//   attrs.set("BUCKED_SIZE", h.info.size_bucket);
//   attrs.set("FIELD_0_NAME",h.info.name);
//   attrs.set("FIELD_1_NAME",h.name);
//   if(ue) {
//     attrs.set("FIELD_2_NAME","rmse");
//   }
}


H5::Group DebugOutVessels(const Grower &grower, const string &name)
{
  const VesselList3d &vl = grower.get_vl();
  static int number = 0;
  H5::H5File f;
  if ( number==0)
  {
    f = H5::H5File("dbgvessels.h5", H5F_ACC_TRUNC);
  }
  else
  {
    f = H5::H5File("dbgvessels.h5", H5F_ACC_RDWR);
  }
  //h5::File f("dbgvessels.h5", number==0 ? "w" : "a");
  H5::Group grp = f.createGroup(str(format("%s") % name));
  H5::Group vesselgrp = grp.createGroup("vessels");
  H5::Group h5_nodes = vesselgrp.openGroup("nodes");
  //since world coordinates this is also needed for proper hdf output
  //vesselgrp.attrs().set<std::string>("CLASS","GRAPH");
  writeAttrToH5(vesselgrp, string("CLASS"), string("GRAPH"));
  WriteVesselList3d(vl, vesselgrp, make_ptree("w_all",false)("w_pressure",true)("w_flow",true));
#if GFFIELD_ENABLE
  {
    DynArray<float> gf;
    grower.GetGfAtNodes(gf);
    //h5::create_dataset(vesselgrp.open_group("nodes"), "gf", gf);
    
    writeDataSetToGroup(h5_nodes, string("gf"), gf);
    //h5::create_dataset(vesselgrp.open_group("nodes"), "gf", gf);
  }
  H5::Group field_ld_grp = grp.createGroup("field_ld");
  //WriteHdfLd(field_ld_grp, grower.get_field_ld());
  //grower.get_field_ld().Lattice2Hdf(field_ld_grp);
  grower.get_field_ld().WriteHdfLd(field_ld_grp);
  WriteScalarField(grp, "gf", grower.GetGf(), grower.get_field_ld(), field_ld_grp);
#endif
  {
    DynArray<uchar> tmp2(vl.GetNCount());
    for (int i=0; i<vl.GetNCount(); ++i)
    {
      tmp2[i] = vl.GetNode(i)->flags;
    }
    //h5::create_dataset(vesselgrp.open_group("nodes"), "nodeflags", tmp2);
    writeDataSetToGroup(h5_nodes, string("nodeflags"), tmp2);
#ifdef WRITE_REMODELING_ACTIONS    
    for (int i=0; i<vl.GetNCount(); ++i)
    {
      tmp2[i] = grower.last_remodeling_action(vl.GetNode(i));
    }
    //h5::create_dataset(vesselgrp.open_group("nodes"), "action", tmp2);
    //writeDataSetToGroup<DynArray<uchar>>(vesselgrp.openGroup("nodes"), string("action"), tmp2);
    writeDataSetToGroup(h5_nodes, string("action"), tmp2);
#endif
  }
  ++number;
  return vesselgrp;
}


void DoOutput(H5::Group &root,
              const VesselList3d &vl,
              const TreeRootList &tree_roots)
{
#ifndef NDEBUG
  std::cout << "DoOutput for vl in vessgen_support" << std::endl;
#endif
  
  const VesselList3d::LatticeData &ld = vl.Ld();

  int num_sites,num_occ_sites;
  float bloodVolume = MeasureBloodVolume(vl, 10.0f, num_sites, num_occ_sites );
  BasicHistogram1D<float> plen,hrad,hbranch,hlenbyrad;
  my::Averaged<float> mlen = MeasureBranchLengths(vl, plen, hlenbyrad );
  my::Averaged<float> mrad = MeasureRadiDistribution(vl, hrad );
  MeasureBranchNumbers(vl, hbranch );
  double arad,aflow,arootcnt,vrad,vflow,vrootcnt;
  MeasureRoot(vl,arad,aflow,arootcnt,vrad,vflow,vrootcnt);

  cout << "ouput -> " << root.getFileName() << ":" << getH5GroupName(root) << endl;

  
  // vessels and stuff
   
  H5::Group vesselgrp;
  H5::Group h5_nodes;
  H5::Group h5_edges;
  H5::Group g;
  H5::Group gg;
  H5::DataSet h5_node_flags;
  
  try
  {
    //if this is created for the first time!
    H5::Exception::dontPrint();
    vesselgrp = root.openGroup("vessels");
    // since world coordinates are introduced this is needed for proper hdf output
    
    //writeAttrToH5(vesselgrp, string("CLASS"), string("GRAPH"));
  }
  catch(H5::Exception &e)
  {
    e.dontPrint();
    vesselgrp = root.createGroup("vessels");
    cout<<"catched vesselgroup" << endl;
    //if I put it like this, it happens only once.
    //writeAttrToH5(vesselgrp, string("CLASS"), string("GRAPH"));
  }
  cout<<"created vesselgroup" << endl;
  
  WriteVesselList3d(vl, vesselgrp, make_ptree("w_all",false)("w_pressure",true)("w_flow",true));
//   try
//   {
//     WriteVesselList3d(vl, vesselgrp, make_ptree("w_all",false)("w_pressure",true)("w_flow",true));
//   }
//   catch(H5::Exception &e)
//   {
//     e.printErrorStack();
//   }
  //07.19.2018 maybe stack allocation is not enought --> I switch to heap
  //07.22.2018 transfered issue to the writeDataSet routine
  DynArray<uchar> tmp2(vl.GetNCount());
  cout<<"allowcated tmp2" << endl;
  //fill tmp with flags
  for (int i=0; i<vl.GetNCount(); ++i)
  {
    tmp2[i] = vl.GetNode(i)->flags;
  }

  
  try
  {
    h5_nodes = vesselgrp.openGroup("nodes");
  }
  catch( H5::Exception &e)
  {
    //e.printErrorStack();
    h5_nodes = vesselgrp.createGroup("nodes");
    cout << "catched nodes" << endl;
    e.dontPrint();
  }
    
  try
  {
    h5_node_flags = h5_nodes.openDataSet("nodeflags");
  }
  catch( H5::Exception &e)
  {
    //only write, if it is not there!
    writeDataSetToGroup(h5_nodes, "nodeflags", tmp2);
    e.dontPrint();
    cout << "catched flags" << endl;
  }
  cout<<"deleted tmp2" << endl;
  
  DynArray<int> tmp3(vl.GetECount());
  
  cout<<"allowcated tmp3" << endl;
  
  for (int i=0; i<vl.GetECount(); ++i)
  {
    tmp3[i] = vl.GetEdge(i)->timeSprout;
  }
  
  
  try
  {
    h5_edges = root.openGroup("vessels/edges");
  }
  catch( H5::Exception &e)
  {
    h5_edges = root.createGroup("vessels/edges");
    e.dontPrint();
    cout << "catched edges" << endl;
  }
    
  writeDataSetToGroup(h5_edges, string("level"), tmp3);
  cout<<"deleted tmp3" << endl;
  
  MemUsage memusage = GetMemoryUsage();
  writeAttrToH5(root, string("mem_vsize"),(int)memusage.vmem_peak );
  writeAttrToH5(root, string("mem_rss"),(int)memusage.rss_peak );

    
  // measurement
  cout <<" data" << endl;
  try
  {
    g = root.openGroup("data");
  }
  catch( H5::Exception &e)
  {
    cout << "catch data" << endl;
    g = root.createGroup("data");
    e.printErrorStack();
  }
  
  cout <<" after data" << endl;
  WriteHdfHistogram(g,"lengths_prob",plen);
  WriteHdfHistogram(g,"lengths_by_rad",hlenbyrad);
  WriteHdfHistogram(g,"radii_prob",hrad);
  WriteHdfHistogram(g,"num_branches_by_rad",hbranch);
  writeAttrToH5(g, string("rBV"), bloodVolume);
  writeAttrToH5(g, string("SITES_OCCUPIED"), num_occ_sites);
  writeAttrToH5(g, string("SITES_TOTAL"), num_sites);
  writeAttrToH5(g, string("MEAN_BRANCH_LENGTH"), mlen.Avg());
  writeAttrToH5(g, string("MEAN_BRANCH_RADIUS"), mrad.Avg());
  writeAttrToH5(g, string("ROOT_A_RADIUS"), arad);
  writeAttrToH5(g, string("ROOT_A_FLOW"), aflow);
  writeAttrToH5(g, string("ROOT_A_COUNT"), arootcnt);
  writeAttrToH5(g, string("ROOT_V_RADIUS"), vrad);
  writeAttrToH5(g, string("ROOT_V_FLOW"), vflow);
  writeAttrToH5(g, string("ROOT_V_COUNT"), vrootcnt);
  
  // roots
  cout << " roots " << endl;
  int N = tree_roots.size();
  DynArray<int64> pos(N);
  DynArray<int> len(N);
  DynArray<char> dir(N);
  DynArray<uchar> flags(N);
  
  auto it = tree_roots.begin();
  for (int i=0; i<N; ++i, ++it)
  {
    const TreeRoot &e = it->second;
    pos[i] = ld.LatticeToSite(e.p);
    len[i] = e.len;
    dir[i] = e.dir;
    flags[i] = e.flags;
  }
  
  try
  {
    gg = g.openGroup("roots");
  }
  catch(H5::Exception &e)
  {
    gg = g.createGroup("roots");
    e.printErrorStack();
    cout << "catch roots" << endl;
  }
  writeDataSetToGroup(gg, string("lattice_pos"), pos);
  writeDataSetToGroup(gg, string("flags"), flags);
  writeDataSetToGroup(gg, string("len"), len);
  writeDataSetToGroup(gg, string("dir"), dir);
  
  vesselgrp.close();
  h5_nodes.close();
  h5_edges.close();
  g.close();
  gg.close();
  h5_node_flags.close();
  cout<< "Error 5 DoOutput" << std::endl;
}


void DoOutput(H5::H5File &file,
              const VesselList3d &vl,
              const Grower &grower,
              const boost::property_tree::atree &additional_data,
              const ptree &input_pt
             )
{
#ifndef NDEBUG
  std::cout << "DoOutput for Grower in vessgen_support" << std::endl;
#endif
  const VesselList3d::LatticeData &ld = vl.Ld();
  H5::Group root = file.openGroup("/");
  
  DoOutput(root, vl, grower.get_tree_roots());

      
        
#if GFFIELD_ENABLE
  {
    //DynArray<float> tmp;
    DynArray<float> *tmp = new DynArray<float>(vl.GetNCount());
    grower.GetGfAtNodes(*tmp);
    //h5::create_dataset(root.open_group("vessels/nodes"), "gf", tmp);
    //h5::create_dataset(root.open_group("vessels/nodes"), "gf", tmp);
    H5::Group h5_nodes= root.openGroup("vessels/nodes");
    writeDataSetToGroup(h5_nodes, string("gf"), *tmp);
    //tmp.clear();
    delete tmp;
  }
#endif
  
  writeAttrToH5(root, string("real_time"), (my::Time() - additional_data.get<my::Time>("real_start_time")).to_s() );
  
  // measurement
  H5::Group h5_data;
  try
  {
    h5_data = root.openGroup("data");
  }
  catch(H5::Exception &e)
  {
    e.dontPrint();
    h5_data = root.createGroup("data");
  }
  {
    const int n = additional_data.get<int>("num_iter");
    Array3d<float> t(Int3(3,n,1));
    int i = 0;
    BOOST_FOREACH(const boost::property_tree::atree::value_type &vv, additional_data.get_child("iters"))
    {
      t(0,i,0) = vv.second.get<int>("iter");
      t(1,i,0) = vv.second.get<float>("rBV");
      t(2,i,0) = vv.second.get<int>("sites");
      ++i;
    }
    H5::DataSet ds = WriteArray3D(h5_data,string("rBV_by_iter"),t);
    writeAttrToH5(ds, string("FIELD_0_NAME"), string("iter"));
    writeAttrToH5(ds, string("FIELD_1_NAME"), string("rBV"));
  }
  
  H5::Group h5_param;
  try
  {
    h5_param = root.openGroup("parameters");
  }
  catch(H5::Exception &e)
  {
    e.dontPrint();
    h5_param = root.createGroup("parameters");
  }
        
  ptree parameter_pt = input_pt;// removes definitions since we already got those
        
        
        
        //parameter_pt.erase("roots");
        //h5::Dataset::create_scalar(g, "SEED", input_pt.get<uint>("seed"));
        //h5::Dataset::create_scalar(g, "IN_FILENAME", input_pt.get<string>("input_fn", ""));
        //h5::Dataset::create_scalar(g, "MESSAGE", input_pt.get<string>("message",""));
        //h5::Dataset::create_scalar(g, "ENSEMBLE_INDEX", input_pt.get<int>("ensemble_index", 0));
  boost::optional<int> seed = input_pt.get_optional<int>("seed");
  if(seed)
  {
    writeAttrToH5(h5_param, string("SEED"), seed.get());
    parameter_pt.erase("seed"); 
  }
  
  boost::optional<string> in_filename = input_pt.get_optional<string>("input_fn");
  if( in_filename)
  {
    writeAttrToH5(h5_param, string("IN_FILENAME"), in_filename.get());
    parameter_pt.erase("input_fn");
  }
  else
  {
    writeAttrToH5(h5_param, string("IN_FILENAME"), string(""));
  }
  
  boost::optional<string> message = input_pt.get_optional<string>("message");
  if( message)
  {
    writeAttrToH5(h5_param, string("MESSAGE"), message.get());
    parameter_pt.erase("message");
  }
  else
  {
    writeAttrToH5(h5_param, string("MESSAGE"), string(""));
  }
  boost::optional<int> ensemble_index = input_pt.get_optional<int>("ensemble_index");
  if( ensemble_index)
  {
    writeAttrToH5(h5_param, string( "ENSEMBLE_INDEX"), ensemble_index.get());
    parameter_pt.erase("ensemble_index");
  }
  else
  {
    writeAttrToH5(h5_param, string( "ENSEMBLE_INDEX"), 0);
  }
	
//         g.attrs().set("SEED", input_pt.get<uint>("seed"));
//         g.attrs().set("IN_FILENAME", input_pt.get<string>("input_fn", ""));
//         g.attrs().set("MESSAGE", input_pt.get<string>("message",""));
//         g.attrs().set( "ENSEMBLE_INDEX", input_pt.get<int>("ensemble_index", 0));
        
  WriteHdfPtree(h5_param, parameter_pt, HDF_WRITE_PTREE_AS_ATTRIBUTE);
  std::cout << "after write Ptree" << std::endl;std::cout.flush();
#if GFFIELD_ENABLE
        H5::Group field_ld_grp = root.createGroup("field_ld");
        //WriteHdfLd(field_ld_grp, grower.get_field_ld());
	grower.get_field_ld().WriteHdfLd(field_ld_grp);
        WriteScalarField(root, "gf", grower.GetGf(), grower.get_field_ld(), field_ld_grp);
        WriteScalarField(root, "gfsources", grower.ComputeGfSources(), grower.get_field_ld(), field_ld_grp);
#endif
  
#ifndef NDEBUG
  std::cout << "exit DoOutput" << std::endl;
#endif
}


/*---------------------------------------------------------
  --------------------------------------------------------*/

/* Apparently this measures BLOOD VOLUME NOT MVD??!!!
 */
double MeasureBloodVolume( const VesselList3d& vl, float dest_lattice_scale, int &n_sites, int &n_sites_occupied )
{
  n_sites = Volume(vl.Ld().Box());
  n_sites_occupied = 0;
  for (int i=0; i<vl.GetNCount(); ++i)
  {
    const VesselNode* vc = vl.GetNode(i);
    bool is_circ = false;
    for (int j=0; j<vc->Count(); ++j)
      if (vc->GetEdge(j)->IsCapillary() && vc->GetEdge(j)->IsCirculated()) is_circ |= true;
    if (is_circ)
      n_sites_occupied++;
  }

#if 0
  float multi = vl.Ld().Scale()/dest_lattice_scale;
  double cnt_occupied = vl.GetNCount();
  for( int i=0; i<vl.GetECount(); ++i )
  {
    const Vessel* v = vl.GetEdge(i);
    if( !v->IsCirculated() ) continue;
    cnt_occupied += (v->len-1)*(multi-1);
    n_sites_occupied += (v->len-2);
  }

  double mvd = cnt_occupied/product_of_components((Size(vl.Ld().Box())-Int3(1))*multi+Int3(1));
#endif

  FloatBBox3 bb = vl.Ld().GetWorldBox();
  double dom_vol = 1.;
  for (int i=0; i<3; ++i) dom_vol *= (bb.max[i]>bb.min[i]) ? (bb.max[i]-bb.min[i]) : 1.;
  //double len = 0.;
  double vol = 0.;
  for (int i=0; i<vl.GetECount(); ++i)
  {
    const Vessel* v = vl.GetEdge(i);
    if (!v->IsCirculated()) continue;
    double l = v->WorldLength(vl.Ld());
    vol += l * v->r*v->r * my::mconst::pi();
  }
  double bloodVolume = vol / dom_vol;
  /** might also break for uncirculated network */
  //myAssert(bloodVolume>0.);
  return bloodVolume;
}


struct BranchDat
{
  float r; // sum of the radii encountered. Will be used to compute the average radius of this branch.
  float l; // world length, in microns
  int   n; // number of encountered segments
};

/* Recursion procedure for MeasureBranchLengths.
 * Will visit all adjacent vessels that are attached in a chain, 
 * and will stop at junctions with more than two attached vessels. */
BranchDat BranchLengthRecurse( const Vessel *v, const VesselList3d::LatticeData &ld, uchar* marker )
{
  BranchDat b;
  b.l = v->WorldLength(ld);
  b.n = 1;
  b.r = v->r;
  marker[v->Index()] = 1;  // mark as visited
  for( int i=0; i<2; ++i )  // look at both end points
  {
    if( v->GetNode(i)->Count()==2 )  // is there one and only one other vessels?
    {
      const Vessel* w = v->GetNode(i)->GetEdge(0);
      if( w==v ) w = v->GetNode(i)->GetEdge(1);  // not going to go where we came from.
      if( !marker[w->Index()] )   // protection against closed vascular loops. If this was not there, the algorithm would go in circles around the loop forever.
      {
        BranchDat bb = BranchLengthRecurse(w,ld,marker);
        b.l += bb.l;
        b.r += bb.r;
        b.n += bb.n;
      }
    }
  }
  return b;
}

std::vector<BranchDat> MeasureBranchLengths( const VesselList3d& vl, const std::vector<bool> wanted)
{
  my::Averaged<float> av; 
  std::vector<uchar> marker( vl.GetECount() );
  std::vector<BranchDat> lengths; 
  //lengths.reserve(1024);
  for( int i=0; i<vl.GetECount(); ++i )
  {
    const Vessel* v = vl.GetEdge(i);
    if( marker[v->Index()] or !wanted[i]) continue;
    BranchDat b = BranchLengthRecurse( v, vl.Ld(), &marker[0] );
    lengths.push_back( b );
    av.Add( b.l );
  }

  return lengths;
}

const my::Averaged<float> MeasureBranchLengths( const VesselList3d& vl, Histo &distrib, Histo &byrad )
{
  typedef Histo HT;
  my::Averaged<float> av; // average length (?)
  distrib.Clear();
  std::vector<uchar> marker( vl.GetECount() );  // edges that we already visited
  std::vector<BranchDat> lengths; lengths.reserve(1024);  // data for each branch, i.e. a section composed of several vessel segments in a chain
  for( int i=0; i<vl.GetECount(); ++i )
  {
    const Vessel* v = vl.GetEdge(i);
    if( marker[v->Index()] ) continue;  // visited within  BranchLengthRecurse, so just continue
    BranchDat b = BranchLengthRecurse( v, vl.Ld(), &marker[0] );  // look at another branch
    lengths.push_back( b );  // store result
    av.Add( b.l );
  }

  // compute histograms
  // @thierry: might as well do this in python. Much easier!
  distrib.Init( HT::HInfo(av.Min(),av.Max(),20u,false, "len" ), HistoBase::Normalized, "prob" );
  for( int i=0; i<lengths.size(); ++i )
  {
    distrib.Fill( lengths[i].l, 1.0 );
  }

  byrad.Init( HT::HInfo(0,200,40u,false,"rad"), HistoBase::SquareSum, "branchlen" );
  byrad.ResetCountMode( true );
  for( int i=0; i<lengths.size(); ++i )
  {
    byrad.Fill( lengths[i].r/lengths[i].n, lengths[i].l );
  }
  return av;
}


const my::Averaged<float> MeasureRadiDistribution( const VesselList3d& vl, Histo &hrad)
{
  typedef Histo HT;
  my::Averaged<float> av;
  hrad.Init( HT::HInfo(0,200,40u,false,"rad"), HistoBase::Normalized, "prob" );
  for( int i=0; i<vl.GetECount(); ++i )
  {
    const Vessel* v = vl.GetEdge(i);
    hrad.Fill( v->r, 1.0f );
    av.Add( v->r );
  }
  return av;
}

void MeasureBranchNumbers( const VesselList3d &vl,  Histo &hbranch )
{
  typedef Histo HT;
  hbranch.Init( HT::HInfo(0,200,40u,false,"rad"), HistoBase::SquareSum, "count" );
  hbranch.ResetCountMode( true );
  for( int i=0; i<vl.GetNCount(); ++i )
  {
    const VesselNode* vc = vl.GetNode(i);
    my::MinMax<float> r;
    for( int j=0; j<vc->Count(); ++j ) r.add(vc->GetEdge(j)->r);
    hbranch.Fill( r.max, vc->Count() );
  }
}

void MeasureRoot( const VesselList3d& vl, double& arad, double& aflow, double& arootcnt, double& vrad, double& vflow, double& vrootcnt )
{
  arad = aflow = arootcnt = vrad = vflow = vrootcnt = 0.;
  for( int i=0; i<vl.GetNCount(); ++i )
  {
    const VesselNode* root = vl.GetNode(i);
    if( !root->flags.GetBits( BOUNDARY ) ) continue;
    for( int k=0; k<root->Count(); ++k )
    {
      const Vessel* v = root->GetEdge(k);
      if(v->flags.GetBits(ARTERY))
      {
        arad += v->r;
        aflow += v->q;
        arootcnt += 1.;
      }
      else if(v->flags.GetBits(VEIN))
      {
        vrad += v->r;
        vflow += v->q;
        vrootcnt += 1.;
      }
    }
  }
  arad /= arootcnt;
  aflow /= arootcnt;
  vrad /= vrootcnt;
  vflow /= vrootcnt;
}


