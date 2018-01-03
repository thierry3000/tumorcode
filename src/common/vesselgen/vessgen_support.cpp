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
    writeAttrToH5(ds,"FIELD_2_NAME","rmse" );
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
  //since world coordinates this is also needed for proper hdf output
  //vesselgrp.attrs().set<std::string>("CLASS","GRAPH");
  writeAttrToH5(vesselgrp, string("CLASS"), string("GRAPH"));
  WriteVesselList3d(vl, vesselgrp, make_ptree("w_all",false)("w_pressure",true)("w_flow",true));
#if GFFIELD_ENABLE
  {
    DynArray<float> gf;
    grower.GetGfAtNodes(gf);
    //h5::create_dataset(vesselgrp.open_group("nodes"), "gf", gf);
    writeDataSetToGroup(vesselgrp.openGroup("nodes"), string("gf"), gf);
    //h5::create_dataset(vesselgrp.open_group("nodes"), "gf", gf);
  }
  H5::Group field_ld_grp = grp.createGroup("field_ld");
  WriteHdfLd(field_ld_grp, grower.get_field_ld());
  WriteScalarField(grp, "gf", grower.GetGf(), grower.get_field_ld(), field_ld_grp);
#endif
  {
    DynArray<uchar> tmp2(vl.GetNCount());
    for (int i=0; i<vl.GetNCount(); ++i)
    {
      tmp2[i] = vl.GetNode(i)->flags;
    }
    //h5::create_dataset(vesselgrp.open_group("nodes"), "nodeflags", tmp2);
#ifdef WRITE_REMODELING_ACTIONS    
    for (int i=0; i<vl.GetNCount(); ++i)
    {
      tmp2[i] = grower.last_remodeling_action(vl.GetNode(i));
    }
    //h5::create_dataset(vesselgrp.open_group("nodes"), "action", tmp2);
    writeDataSetToGroup<DynArray<uchar>>(vesselgrp.openGroup("nodes"), string("action"), tmp2);
#endif
  }
  ++number;
  return vesselgrp;
}


void DoOutput(H5::Group root,
              const VesselList3d &vl,
              const TreeRootList &tree_roots)
{
  const VesselList3d::LatticeData &ld = vl.Ld();

  int num_sites,num_occ_sites;
  float bloodVolume = MeasureBloodVolume(vl, 10.0f, num_sites, num_occ_sites );
  BasicHistogram1D<float> plen,hrad,hbranch,hlenbyrad;
  my::Averaged<float> mlen = MeasureBranchLengths(vl, plen, hlenbyrad );
  my::Averaged<float> mrad = MeasureRadiDistribution(vl, hrad );
  MeasureBranchNumbers(vl, hbranch );
  double arad,aflow,arootcnt,vrad,vflow,vrootcnt;
  MeasureRoot(vl,arad,aflow,arootcnt,vrad,vflow,vrootcnt);

  //cout << "ouput -> " << root.get_file().get_file_name() << ":" << root.get_name() << endl;
  cout << "ouput -> " << root.getFileName() << ":" << getH5Name(root) << endl;

  {
    //h5::Attributes a;
    // vessels and stuff
    H5::Group vesselgrp = root.createGroup("vessels");
    // since world coordinates are introduced this is needed for proper hdf output
    writeAttrToH5(vesselgrp, string("CLASS"), string("GRAPH"));
    //vesselgrp.attrs().set<std::string>("CLASS","GRAPH");
    WriteVesselList3d(vl, vesselgrp, make_ptree("w_all",false)("w_pressure",true)("w_flow",true));
    {
//       DynArray<uchar> tmp2(vl.GetNCount());
//       for (int i=0; i<vl.GetNCount(); ++i)
//       {
//         tmp2[i] = vl.GetNode(i)->flags;
//       }
//       h5::create_dataset(root.open_group("vessels/nodes"), "nodeflags", tmp2);
      DynArray<int> tmp3(vl.GetECount());
      for (int i=0; i<vl.GetECount(); ++i)
      {
        tmp3[i] = vl.GetEdge(i)->timeSprout;
      }
      //h5::create_dataset(root.open_group("vessels/edges"), "level", tmp3);
      writeDataSetToGroup(root.openGroup("vessels/edges"), string("level"), tmp3);
    }
    {
      MemUsage memusage = GetMemoryUsage();
      writeAttrToH5(root, string("mem_vsize"),(int)memusage.vmem_peak );
      writeAttrToH5(root, string("mem_rss"),(int)memusage.rss_peak );
      //a = root.attrs();
//       a.set<uint64>("mem_vsize", memusage.vmem_peak);
//       a.set<uint64>("mem_rss", memusage.rss_peak);
    }
    // measurement
    H5::Group g = root.createGroup("data");
    WriteHdfHistogram(g,"lengths_prob",plen);
    WriteHdfHistogram(g,"lengths_by_rad",hlenbyrad);
    WriteHdfHistogram(g,"radii_prob",hrad);
    WriteHdfHistogram(g,"num_branches_by_rad",hbranch);
    writeAttrToH5(g, string("rBV"), bloodVolume);
//     a = g.attrs();
//     a.set("rBV",bloodVolume);
//     a.set("SITES_OCCUPIED",num_occ_sites);
//     a.set("SITES_TOTAL",num_sites);
//     a.set("MEAN_BRANCH_LENGTH",mlen.Avg());
//     a.set("MEAN_BRANCH_RADIUS",mrad.Avg());
//     a.set("ROOT_A_RADIUS",arad);
//     a.set("ROOT_A_FLOW",aflow);
//     a.set("ROOT_A_COUNT",arootcnt);
//     a.set("ROOT_V_RADIUS",vrad);
//     a.set("ROOT_V_FLOW",vflow);
//     a.set("ROOT_V_COUNT",vrootcnt);
    // roots
    {
      int N = tree_roots.size();
      DynArray<int> pos(N), len(N);
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
      H5::Group gg = g.createGroup("roots");
      writeDataSetToGroup<DynArray<int>>(gg, string("lattice_pos"), pos);
      writeDataSetToGroup<DynArray<uchar>>(gg, string("flags"), flags);
      writeDataSetToGroup<DynArray<int>>(gg, string("len"), len);
      writeDataSetToGroup<DynArray<char>>(gg, string("dir"), dir);
//       h5::create_dataset(gg, "lattice_pos", pos); int
//       h5::create_dataset(gg, "flags", flags); uchar
//       h5::create_dataset(gg, "len", len); int
//       h5::create_dataset(gg, "dir", dir); char
    }
  }
}


void DoOutput(H5::H5File &file,
              const VesselList3d &vl,
              const Grower &grower,
              const boost::property_tree::atree &additional_data,
              const ptree &input_pt
             )
{
      const VesselList3d::LatticeData &ld = vl.Ld();

      DoOutput(file.openGroup("/"), vl, grower.get_tree_roots());

      {
        //h5::Attributes a;
        H5::Group root = file.openGroup("/");
#if GFFIELD_ENABLE
        {
          DynArray<float> tmp;
          grower.GetGfAtNodes(tmp);
          //h5::create_dataset(root.open_group("vessels/nodes"), "gf", tmp);
	  //h5::create_dataset(root.open_group("vessels/nodes"), "gf", tmp);
	  writeDataSetToGroup<DynArray<float>>(root.openGroup("vessels/nodes"), string("gf"), tmp);
          tmp.clear();
        }
#endif
        {
          //a = root.attrs();
	  writeAttrToH5(root, string("real_time"), (my::Time() - additional_data.get<my::Time>("real_start_time")).to_s() );
          //a.set<double>("real_time", (my::Time() - additional_data.get<my::Time>("real_start_time")).to_s());
        }
        // measurement
        H5::Group g = root.createGroup("data");
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
          H5::DataSet ds = WriteArray3D(g,"rBV_by_iter",t);
          //a = ds.attrs();
	  writeAttrToH5(ds, string("FIELD_0_NAME"), string("iter"));
	  writeAttrToH5(ds, string("FIELD_1_NAME"), string("rBV"));
//           a.set<string>("FIELD_0_NAME","iter");
//           a.set<string>("FIELD_1_NAME","rBV");
        }

        // parameters
        //g = root.require_group("parameters");
	g = root.createGroup("parameters");
        //h5::Dataset::create_scalar(g, "SEED", input_pt.get<uint>("seed"));
        //h5::Dataset::create_scalar(g, "IN_FILENAME", input_pt.get<string>("input_fn", ""));
        //h5::Dataset::create_scalar(g, "MESSAGE", input_pt.get<string>("message",""));
        //h5::Dataset::create_scalar(g, "ENSEMBLE_INDEX", input_pt.get<int>("ensemble_index", 0));
	writeAttrToH5(g, string("SEED"), input_pt.get<int>("seed"));
	writeAttrToH5(g, string("IN_FILENAME"), input_pt.get<string>("input_fn", ""));
	writeAttrToH5(g, string("MESSAGE"), input_pt.get<string>("message",""));
	writeAttrToH5(g, string( "ENSEMBLE_INDEX"), input_pt.get<int>("ensemble_index", 0));
	
//         g.attrs().set("SEED", input_pt.get<uint>("seed"));
//         g.attrs().set("IN_FILENAME", input_pt.get<string>("input_fn", ""));
//         g.attrs().set("MESSAGE", input_pt.get<string>("message",""));
//         g.attrs().set( "ENSEMBLE_INDEX", input_pt.get<int>("ensemble_index", 0));
        ptree parameter_pt = input_pt;
        parameter_pt.erase("seed"); // removes definitions since we already got those
        parameter_pt.erase("input_fn");
        parameter_pt.erase("message");
        parameter_pt.erase("ensemble_index");
        parameter_pt.erase("roots"); 
        WriteHdfPtree(g, parameter_pt, HDF_WRITE_PTREE_AS_ATTRIBUTE);
#if GFFIELD_ENABLE
        H5::Group field_ld_grp = root.createGroup("field_ld");
        WriteHdfLd(field_ld_grp, grower.get_field_ld());
        WriteScalarField(root, "gf", grower.GetGf(), grower.get_field_ld(), field_ld_grp);
        WriteScalarField(root, "gfsources", grower.ComputeGfSources(), grower.get_field_ld(), field_ld_grp);
#endif
      }
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


