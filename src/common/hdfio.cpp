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
#include "hdfio.h"

#include "mwlib/lattice-data.h"
#include "vessels3d.h"
#include "continuum-grid.h"

#include <boost/foreach.hpp>
#include <boost/multi_array.hpp>
#include <stdexcept>

template<class T>
H5::DataType getH5TypeFromCpp()
{
  H5::DataType thisWritingType;
  if(typeid(T) == typeid(float))
  {
    thisWritingType = H5::PredType::NATIVE_FLOAT;
  }
  else if(typeid(T) == typeid(Float3))
  {
    thisWritingType = H5::PredType::NATIVE_FLOAT;
  }
  else if(typeid(T) == typeid(Float2))
  {
    thisWritingType = H5::PredType::NATIVE_FLOAT;
  }
  else if(typeid(T) == typeid(double))
  {
    thisWritingType = H5::PredType::NATIVE_DOUBLE;
  }
  else if(typeid(T) == typeid(Double3))
  {
    thisWritingType = H5::PredType::NATIVE_DOUBLE;
  }
  else if(typeid(T) == typeid(Bool3))
  {
    thisWritingType = H5::PredType::NATIVE_CHAR;
  }
  else if(typeid(T) == typeid(int))
  {
    thisWritingType = H5::PredType::NATIVE_INT;
  }
  else if(typeid(T) == typeid(Int3))
  {
    thisWritingType = H5::PredType::NATIVE_INT;
  }
  else if(typeid(T) == typeid(Int6))
  {
    thisWritingType = H5::PredType::NATIVE_INT;
  }
  else if(typeid(T) == typeid(bool))
  {
    thisWritingType = H5::PredType::NATIVE_CHAR;
  }
  else if(typeid(T) == typeid(char))
  {
    thisWritingType = H5::PredType::NATIVE_CHAR;
  }
  else if(typeid(T) == typeid(uchar))
  {
    thisWritingType = H5::PredType::NATIVE_UCHAR;
  }
  else if(typeid(T) == typeid(long))
  {
    thisWritingType = H5::PredType::NATIVE_LONG;
  }
  else if(typeid(T) == typeid(int64))
  {
    thisWritingType = H5::PredType::NATIVE_INT64;
  }
  else if(typeid(T) == typeid(unsigned long))
  {
    thisWritingType = H5::PredType::NATIVE_UINT64;
  }
  else
  {
    cout << "unsupported Template type in writeDataSetToGroup!" << endl;
    exit(1);
  }
  return thisWritingType;
}




// template<class U,class T>
// void readAttrFromH5(U &g, const string &attr_name, T &output_buffer)
// {
//   H5::Attribute att_to_read;
//   H5::DataType type; 
//   try
//   {
//     att_to_read = g.openAttribute(attr_name);
//     type = att_to_read.getDataType();
//     att_to_read.read(type, &output_buffer);
//   }
//   catch(H5::Exception error)
//   {
//     error.printErrorStack();
//   }
// }
template<class T>
void readAttrFromH5(H5::H5Location &g, const string &attr_name, T &output_buffer)
{
  H5::Attribute att_to_read;
  H5::DataType type; 
  try
  {
    att_to_read = g.openAttribute(attr_name);
    type = att_to_read.getDataType();
    att_to_read.read(type, &output_buffer);
  }
  catch(H5::Exception error)
  {
    error.printErrorStack();
  }
}


// template<>
// void readAttrFromH5<H5::DataSet, string>(H5::DataSet &g, const string &attr_name, string &output_buffer)
// { 
//   H5::Attribute att_to_read = g.openAttribute(attr_name);
//   //H5::DataType type = att_to_read.getDataType();
//   
//   // Create new dataspace for attribute
//   H5::DataSpace attr_dataspace = H5::DataSpace(H5S_SCALAR);
// 
//   // Create new string datatype for attribute
//   H5::StrType strdatatype(H5::PredType::C_S1, H5T_VARIABLE); // of length 256 characters
// 
//   // Set up read buffer for attribute
//   H5std_string strreadbuf ("");
// 
//   // Create attribute and write to it
//   try
//   {
//     att_to_read.read(strdatatype, strreadbuf);
//   }
//   catch(H5::Exception error)
//   {
//     error.printErrorStack();
//   }
//   output_buffer = strreadbuf;
// }
// template<>
// void readAttrFromH5<H5::Group, string>(H5::Group &g, const string &attr_name, string &output_buffer)
// { 
//   H5::Attribute att_to_read = g.openAttribute(attr_name);
//   //H5::DataType type = att_to_read.getDataType();
//   
//   // Create new dataspace for attribute
//   H5::DataSpace attr_dataspace = H5::DataSpace(H5S_SCALAR);
// 
//   // Create new string datatype for attribute
//   H5::StrType strdatatype(H5::PredType::C_S1, H5T_VARIABLE); // of length 256 characters
// 
//   // Set up read buffer for attribute
//   H5std_string strreadbuf ("");
// 
//   // Create attribute and write to it
//   try
//   {
//     att_to_read.read(strdatatype, strreadbuf);
//   }
//   catch(H5::Exception error)
//   {
//     error.printErrorStack();
//   }
//   output_buffer = strreadbuf;
// }

template<>
void readAttrFromH5<string>(H5::H5Location &g, const string &attr_name, string &output_buffer)
{ 
  H5::Attribute att_to_read = g.openAttribute(attr_name);
  //H5::DataType type = att_to_read.getDataType();
  
  // Create new dataspace for attribute
  H5::DataSpace attr_dataspace = H5::DataSpace(H5S_SCALAR);

  // Create new string datatype for attribute
  H5::StrType strdatatype(H5::PredType::C_S1, H5T_VARIABLE); // of length 256 characters

  // Set up read buffer for attribute
  H5std_string strreadbuf ("");

  // Create attribute and write to it
  try
  {
    att_to_read.read(strdatatype, strreadbuf);
  }
  catch(H5::Exception error)
  {
    error.printErrorStack();
  }
  output_buffer = strreadbuf;
}

template <class T>
void writeAttrToH5(H5::H5Location &h, const string &attr_name,  const T &value)
{ 
  H5::DataType thisType = getH5TypeFromCpp<T>();
  const int rank = 2;
  hsize_t dims[rank];
  dims[0] = 1;
  if(typeid(T) == typeid(Float3) or typeid(T) == typeid(Int3) or typeid(T) == typeid(Bool3))
  {
    dims[1] = 3;
  }
  else
  {
    dims[1] = 1;
  }
  if(typeid(T) == typeid(Int6))
  {
    dims[1] = 6;
  }
  H5::DataSpace mspace = H5::DataSpace( rank, dims);
  H5::Attribute attr_out;
  try{
    attr_out = h.createAttribute(attr_name, thisType, mspace);
  }
  catch(H5::Exception e)
  {
    e.printErrorStack();
  }
  attr_out.write(thisType, &value);
};
template<>
void writeAttrToH5<string>(H5::H5Location &h, const string &attr_name, const string &value)
{ 
  // Create new dataspace for attribute
  H5::DataSpace attr_dataspace = H5::DataSpace(H5S_SCALAR);
  // Create new string datatype for attribute
  H5::StrType strdatatype(H5::PredType::C_S1, H5T_VARIABLE); // of length 256 characters
  // Set up write buffer for attribute
  const H5std_string strwritebuf (value);
  try
  {
    H5::Attribute myatt_in = h.createAttribute(attr_name, strdatatype, attr_dataspace);
    myatt_in.write(strdatatype, strwritebuf);
  }
  catch(H5::Exception error)
  {
    error.printErrorStack();
  }
};
// template<>
// void writeAttrToH5<H5::DataSet,string>(H5::DataSet &h, const string &attr_name, const string &value)
// { 
//   {
//     // Create new dataspace for attribute
//     H5::DataSpace attr_dataspace = H5::DataSpace(H5S_SCALAR);
//     // Create new string datatype for attribute
//     H5::StrType strdatatype(H5::PredType::C_S1, H5T_VARIABLE); // of length 256 characters
//     // Set up write buffer for attribute
//     const H5std_string strwritebuf (value);
//     try
//     {
//       H5::Attribute myatt_in = h.createAttribute(attr_name, strdatatype, attr_dataspace);
//       myatt_in.write(strdatatype, strwritebuf);
//     }
//     catch(H5::Exception error)
//     {
//       error.printErrorStack();
//     }
//   }
// };

H5::Group RequireLatticeDataGroup(H5::H5File &f, const string &name, const LatticeDataQuad3d &ld)
{
  bool had_ld_group = false;
  H5::Group ld_group;
  try 
  {  // to determine if the dataset exists in the file
    //dataset = new DataSet(file->openDataSet( "/Data_new/Compressed_Data" ));
    //group = new Group(file->openGroup("Data"));
    ld_group = H5::Group(f.openGroup("name"));
    //ld_group =H5Gopen2(g.getLocId(), name, H5P_DEFAULT);
  }
  catch( H5::GroupIException not_found_error )
  {
      cout << " Lattice group not found in hdf, we create it." << endl;
      //WriteHdfLd(ld_group, ld);
      ld.WriteHdfLd(ld_group);
  }
//   H5::Group ld_group = g.createGroup()
//   h5cpp::Group ld_group = g.require_group(name, &had_ld_group);
//   if (!had_ld_group)
//     WriteHdfLd(ld_group, ld);
  return ld_group;
}

H5::Group RequireLatticeDataGroup(H5::Group &g, const LatticeDataQuad3d &ld)
{
  //"field_ld"
  bool had_ld_group = false;
  H5::Group ld_group;
  try 
  {  // to determine if the dataset exists in the file
    //dataset = new DataSet(file->openDataSet( "/Data_new/Compressed_Data" ));
    //group = new Group(file->openGroup("Data"));
    ld_group = g.openGroup("field_ld");
    //ld_group =H5Gopen2(g.getLocId(), name, H5P_DEFAULT);
  }
  catch( H5::GroupIException not_found_error )
  {
      cout << " Lattice group not found in hdf, we create it." << endl;
      //WriteHdfLd(ld_group, ld);
      ld.WriteHdfLd(ld_group);
  }
//   H5::Group ld_group = g.createGroup()
//   h5cpp::Group ld_group = g.require_group(name, &had_ld_group);
//   if (!had_ld_group)
//     WriteHdfLd(ld_group, ld);
  return ld_group;
}



void WriteHdfGraph( H5::Group &g, const VesselList3d &vl )
{
#ifdef DEBUG
  printf("Starting to write hdf\n");
#endif
  int ncnt = vl.GetNCount();
  int ecnt = vl.GetECount();
  myAssert(ncnt>0 && ecnt>0);
  
  H5::Group h5_nodes;
  try
  {
    h5_nodes = g.createGroup("nodes");
  }
  catch(H5::Exception e)
  {
    //e.printErrorStack();
    h5_nodes = g.openGroup("nodes");
  }
  writeAttrToH5(h5_nodes, string("COUNT"), ncnt);
  //gg.attrs().set("COUNT",ncnt);
  //h5cpp::Attributes attrs = g.attrs();
//   H5::Group gg = g.createGroup("nodes");
//   gg.createAttribute("COUNT");
//   gg.attrs().set("COUNT",ncnt);
//   h5cpp::Attributes attrs = g.attrs();
  if(vl.HasLattice())//LATTICE IS THERE
  {
#ifndef NDEBUG
    std::cout<<"lattice found" << std::endl;
#endif
    //reading data in a write function does not make sense.
    writeAttrToH5(g,string("CLASS"), string("GRAPH"));
    
    // CHECK THIS!!!
//     try
//     {
//       string graph_type;
//       readAttrFromH5(g,string("CLASS"),graph_type);
//     }
//     catch(H5::Exception e)
//     {
//       //e.dontPrint();
//       writeAttrToH5(g,string("CLASS"), string("GRAPH"));
//     }
    
    //lattice stuff is writen to hdf
    H5::Group lattice_group = g.createGroup("lattice");
    //vl.Ld().WriteHdf(g.create_group("lattice"));
    //vl.Ld().WriteHdf(lattice_group);
    const VesselList3d::LatticeData &ld = vl.Ld();
    ld.Lattice2Hdf(lattice_group);
    //ld.writeToHdf(lattice_group);
    //vl.Ld().WriteHdf(lattice_group);
    DynArray<int> a(ncnt);
    for(int i=0; i<ncnt; ++i) 
    {
      a[i] = ld.LatticeToSite(vl.GetNode(i)->lpos);
    }
    H5::DataSet dataset = writeDataSetToGroup(h5_nodes, string("lattice_pos"), a);
    //h5cpp::Dataset ds = h5cpp::create_dataset<int>(gg,  "lattice_pos", a);
//       hsize_t dims[2] = {1,ncnt};
//       int RANK = 2;
//       H5::DataSpace dataspace = H5::DataSpace(RANK, dims); // create new dspace
//       H5::DSetCreatPropList ds_creatplist;
//       H5::DataSet dataset = H5::DataSet(gg.createDataSet("lattice_pos", H5::PredType::NATIVE_INT,dataspace, ds_creatplist ));
//       dataset.write(&a, H5::PredType::NATIVE_INT, dataspace);
    writeAttrToH5(dataset, "MODE", string("linear"));
    //ds.attrs().set("MODE","linear");
  }
  else//no lattice pressent
  {
    try
    {
      string theType;
      readAttrFromH5(g, string("CLASS"), theType);
      //myAssert(theType == "REALWORLD");
    }
    catch( H5::Exception not_found_error )
    {
      not_found_error.dontPrint();
      cout << " CLASS type not found" << endl;
      writeAttrToH5(g,string("CLASS"), string("REALWORLD"));
    }  
    // write the special world stuff
    //needs old school linearizing
// 	      DynArray<float> a(3*ncnt);
// 	      for(int i=0; i<ncnt; ++i) 
// 	      {
// 		a[3*i+0]= vl.GetNode(i)->worldpos[0];
// 		a[3*i+1] = vl.GetNode(i)->worldpos[1];
// 		a[3*i+2] = vl.GetNode(i)->worldpos[2];
// 	      }
    DynArray<Float3> a;
    a.resize(ncnt);
    for(int i=0; i<ncnt; ++i) 
    {
      Float3 abuffer=Float3(vl.GetNode(i)->worldpos[0],
			    vl.GetNode(i)->worldpos[1],
			    vl.GetNode(i)->worldpos[2]);
      a[i] = abuffer;
    }
    // due to creation of h5 files with python this could already be pressent
    
    try
    {
      //DataSet dataset = file.openDataSet(DATASET_NAME);
      H5::DataSet ds = h5_nodes.openDataSet("world_pos");
      //grp=h5file.openGroup("A");
    } 
    catch(H5::Exception& e)
    {
      /* group does not exists, create it */
      //grp=h5file.createGroup("A");
      const int	 RANK = 2;
      hsize_t dims[2] = {3, ncnt};
      H5::DataSpace dataspace(RANK, dims);
      H5::DataSet ds = h5_nodes.createDataSet("world_pos", H5::PredType::NATIVE_FLOAT, dataspace);
      ds.write(&a, H5::PredType::NATIVE_FLOAT);
    }
//     if(!gg.exists("world_pos"))
//     {
//       h5cpp::Dataset ds = h5cpp::create_dataset<float>(gg, "world_pos", h5cpp::Dataspace::simple_dims(a.size()/3,3), &a[0]);
//     }
  }//end no lattice pressent
  
  {//interrupt begin
  //write roots
  DynArray<int> roots(16,ConsTags::RESERVE);
  for(int i=0; i<ncnt; ++i) 
  {
    if(vl.GetNode(i)->IsBoundary())
      roots.push_back(i);
  }
  H5::DataSet h5_dataset_roots;
  try
  {
    //h5_dataset_roots = h5_nodes.openDataSet("roots");
    h5_dataset_roots = writeDataSetToGroup(h5_nodes, string("roots"), roots);
  } 
  catch(H5::Exception& e)
  {
    /* group does not exists, create it */
    e.printErrorStack();
  }

  {
    DynArray <int> bc_node_index;
    DynArray <int> bctyp_index;
    DynArray <float> values_of_bcs;
    DynArray <float> bc_conductivity_value;
    if(!vl.GetBCMap().empty())
    {
      for (auto it = vl.GetBCMap().begin(); it != vl.GetBCMap().end(); ++it)
      {
        bc_node_index.push_back(it->first->Index());
        bctyp_index.push_back(it->second.typeOfInstance);
        values_of_bcs.push_back(it->second.val);
        bc_conductivity_value.push_back(it->second.w);
      }
    }
    else
    {
      //the map has not been used up to now,
      //for the adaption it is good to have any way, so we write it.
      //we assume all pressure boundary conditions
      for( auto const &value: roots)
      {
	const VesselNode *nd = vl.GetNode(value);
	bc_node_index.push_back(value);
	bctyp_index.push_back(FlowBC::PIN);
	values_of_bcs.push_back(nd->press);
	bc_conductivity_value.push_back(0);
      }
    }
    writeDataSetToGroup(h5_nodes, string("bc_node_index"), bc_node_index);
    writeDataSetToGroup(h5_nodes, string("bc_type"), bctyp_index);
    writeDataSetToGroup(h5_nodes, string("bc_value"), values_of_bcs);
    writeDataSetToGroup(h5_nodes, string("bc_conductivity_value"), bc_conductivity_value);

  }
  }//for interrupt 

  H5::Group h5_edges = g.createGroup("edges");
  writeAttrToH5(h5_edges,string("COUNT"), ecnt);
  //gg.attrs().set("COUNT",ecnt);
  
  {//Write edge stuff
    DynArray<int> va(ecnt),vb(ecnt),flags(ecnt);
    DynArray<float> float_radius(ecnt);
    for(int i=0; i<ecnt; ++i) 
    { 
      const Vessel* v = vl.GetEdge(i);
      va[i] = v->NodeA()->Index();
      vb[i] = v->NodeB()->Index();
      flags[i] = v->flags; 
      float_radius[i] = v->r; 
    }
    writeDataSetToGroup(h5_edges, string("node_a_index"), va);
    writeDataSetToGroup(h5_edges, string("node_b_index"), vb);
    writeDataSetToGroup(h5_edges, string("flags"), flags);
    writeDataSetToGroup(h5_edges, string("radius"), float_radius);
    writeAttrToH5(h5_edges, string("MODE"), string("const"));
    
   
//     h5cpp::create_dataset<int>( gg, "node_a_index", va);
//     h5cpp::create_dataset<int>( gg, "node_b_index", vb );
//     h5cpp::create_dataset<int>( gg, "flags", flags );
//     h5cpp::Dataset ds = h5cpp::create_dataset<float>( gg, "radius", float_radius );
//     ds.attrs().set("MODE","const");
  }
}

void ReadHdfGraph( H5::Group &g, VesselList3d *vl )
{
  H5::Group gnodes;
  H5::Group gedges;
  try
  {
    gnodes = g.openGroup("nodes");
    gedges = g.openGroup("edges");
  }
  catch(H5::Exception &e)
  {
    e.printError();
  }
  
  int ecnt=0,ncnt=0;
  readAttrFromH5(gnodes, string("COUNT"), ncnt);
  readAttrFromH5(gedges, string("COUNT"), ecnt);

  if(vl->HasLattice())//LATTICE IS THERE
  {//read lattice stuffe
    const VesselList3d::LatticeData &ld = vl->Ld();
    {//begin interrupt
      
      //node stuff
//       std::vector<VesselList3d::SiteType> a;
//       h5cpp::read_dataset<VesselList3d::SiteType>(gnodes.open_dataset("lattice_pos"),a);
      DynArray<VesselList3d::SiteType> a;
      readDataSetFromGroup(gnodes,string("lattice_pos"),a);
      for(int i=0; i<ncnt; ++i)
      {
        vl->InsertNode(ld.SiteToLattice(a[i]));
      }
    }//end interrupt
  }
  else//no lattice pressent
  {
    {//begin interupt
      // WARNING: UNTESTET
//       DynArray<float> a;
//       h5cpp::Dataset ds = gnodes.open_dataset("world_pos");
//       h5cpp::Dataspace sp = ds.get_dataspace();
//       hsize_t dims[H5S_MAX_RANK];
//       int rank = sp.get_dims(dims);
//       h5cpp::read_dataset<float>(ds,a);
//       for(int i=0; i<ncnt; ++i)
//       {
// 	vl->InsertNode(
// 	  Float3(
// 	    a[3*i+0],a[3*i+1],a[3*i+2]
// 	  )
// 	);
//       }
      DynArray<Float3> a;
      readDataSetFromGroup(gnodes,string("world_pos"), a);
      for(int i=0; i<ncnt; ++i)
      {
	vl->InsertNode(a[i]); //now a[i] is a Float3
      }
    }//end interupt
  }
  /* Read all other stuff which is the same*/
  {
    //edge stuff
    DynArray<int> va,vb;
    va.resize(ncnt);
    vb.resize(ncnt);
//     h5cpp::read_dataset<int>(gedges.open_dataset("node_a_index"),va);
//     h5cpp::read_dataset<int>(gedges.open_dataset("node_b_index"),vb);
    readDataSetFromGroup(gedges, string("node_a_index"), va);
    readDataSetFromGroup(gedges, string("node_b_index"), vb);
    for(int i=0; i<ecnt; ++i)
    {
#ifndef NDEBUG
      std::cout << "va[" << i << "]: " << va[i] << std::endl;
      std::cout << "vb[" << i << "]: " << vb[i] << std::endl;
#endif
      Vessel* v = vl->InsertVessel(vl->GetNode(va[i]),vl->GetNode(vb[i]));
    }
  }
  {
    DynArray<int> root_indices;
    //h5cpp::read_dataset<int>(gnodes.open_dataset("roots"),root_indices);
    readDataSetFromGroup(gnodes, string("roots"), root_indices);
    for(int i=0; i<root_indices.size(); ++i)
      vl->GetNode(root_indices[i])->flags.AddBits(BOUNDARY);
  }
    
//   if(gnodes.exists("bc_node_index") && 
//     gnodes.exists("bc_type") && 
//     gnodes.exists("bc_value") && 
//     gnodes.exists("bc_conductivity_value"))
//   {
//     DynArray <int> bc_node_index;
//     DynArray <int> bctyp_index;
//     DynArray <float> values_of_bcs;
//     DynArray <float> bc_conductivity_value;
//     h5cpp::read_dataset<int>(gnodes.open_dataset("bc_node_index")     ,bc_node_index);
//     h5cpp::read_dataset<int>(gnodes.open_dataset("bc_type")           ,bctyp_index);
//     h5cpp::read_dataset<float>(gnodes.open_dataset("bc_value")          ,values_of_bcs);
//     h5cpp::read_dataset<float>(gnodes.open_dataset("bc_conductivity_value") ,bc_conductivity_value);
//     for(int i=0; i<bc_node_index.size(); ++i)
//     {
//       #ifdef DEBUG
//       cerr<<format("root index %i of #%i, bctype %i, value: %f, conductivity: %f\n") % bc_node_index[i] % bc_node_index.size() % bctyp_index[i] % values_of_bcs[i] % bc_conductivity_value[i];
//       #endif
//       VesselNode* nd = vl->GetNode(bc_node_index[i]);
//       nd->flags.AddBits(BOUNDARY);
//       FlowBC bc; 
//       switch (bctyp_index[i])
//       {
//         case FlowBC::PIN:
//           bc = FlowBC(FlowBC::PIN, values_of_bcs[i]);
//           break;
//         case FlowBC::CURRENT:
//           bc = FlowBC(FlowBC::CURRENT, values_of_bcs[i]);
//           break;
//         case FlowBC::RESIST:
//           bc = FlowBC(FlowBC::RESIST, bc_conductivity_value[i], values_of_bcs[i]);
//           break;
//       }
//       vl->SetBC(nd, bc);
//     }
//   }

//   if(gnodes.exists("bc_node_index") && 
//     gnodes.exists("bc_type") && 
//     gnodes.exists("bc_value") && 
//     gnodes.exists("bc_conductivity_value"))
//   {
//     DynArray <int> bc_node_index;
//     DynArray <int> bctyp_index;
//     DynArray <float> values_of_bcs;
//     DynArray <float> bc_conductivity_value;
//     h5cpp::read_dataset<int>(gnodes.open_dataset("bc_node_index")     ,bc_node_index);
//     h5cpp::read_dataset<int>(gnodes.open_dataset("bc_type")           ,bctyp_index);
//     h5cpp::read_dataset<float>(gnodes.open_dataset("bc_value")          ,values_of_bcs);
//     h5cpp::read_dataset<float>(gnodes.open_dataset("bc_conductivity_value") ,bc_conductivity_value);
//     for(int i=0; i<bc_node_index.size(); ++i)
//     {
//       #ifdef DEBUG
//       cerr<<format("root index %i of #%i, bctype %i, value: %f, conductivity: %f\n") % bc_node_index[i] % bc_node_index.size() % bctyp_index[i] % values_of_bcs[i] % bc_conductivity_value[i];
//       #endif
//       VesselNode* nd = vl->GetNode(bc_node_index[i]);
//       nd->flags.AddBits(BOUNDARY);
//       FlowBC bc; 
//       switch (bctyp_index[i])
//       {
//         case FlowBC::PIN:
//           bc = FlowBC(FlowBC::PIN, values_of_bcs[i]);
//           break;
//         case FlowBC::CURRENT:
//           bc = FlowBC(FlowBC::CURRENT, values_of_bcs[i]);
//           break;
//         case FlowBC::RESIST:
//           bc = FlowBC(FlowBC::RESIST, bc_conductivity_value[i], values_of_bcs[i]);
//           break;
//       }
//       vl->SetBC(nd, bc);
//     }
//   }

/**
 * Error handling is now done inside the hdf5 interface
 */
  DynArray <int> bc_node_index;
  DynArray <int> bctyp_index;
  DynArray <float> values_of_bcs;
  DynArray <float> bc_conductivity_value;
  /* for old files the boundary arrays are not present
    */
  readDataSetFromGroup(gnodes, string("bc_node_index"),bc_node_index );
  readDataSetFromGroup(gnodes, string("bc_type"),bctyp_index );
  readDataSetFromGroup(gnodes, string("bc_value"),values_of_bcs );
  readDataSetFromGroup(gnodes, string("bc_conductivity_value"),bc_conductivity_value );
  for(int i=0; i<bc_node_index.size(); ++i)
  {
    #ifdef DEBUG
    cerr<<format("root index %i of #%i, bctype %i, value: %f, conductivity: %f\n") % bc_node_index[i] % bc_node_index.size() % bctyp_index[i] % values_of_bcs[i] % bc_conductivity_value[i];
    #endif
    VesselNode* nd = vl->GetNode(bc_node_index[i]);
    nd->flags.AddBits(BOUNDARY);
    FlowBC bc; 
    switch (bctyp_index[i])
    {
      case FlowBC::PIN:
        bc = FlowBC(FlowBC::PIN, values_of_bcs[i]);
        break;
      case FlowBC::CURRENT:
        bc = FlowBC(FlowBC::CURRENT, values_of_bcs[i]);
        break;
      case FlowBC::RESIST:
        bc = FlowBC(FlowBC::RESIST, bc_conductivity_value[i], values_of_bcs[i]);
        break;
    }
    vl->SetBC(nd, bc);
  }
    

  DynArray<int> flags;
  DynArray<float> aflt;
//   h5cpp::read_dataset<int>(gedges.open_dataset("flags"),flags);
//   h5cpp::read_dataset<float>(gedges.open_dataset("radius"),aflt);
  readDataSetFromGroup(gedges, string("flags"), flags);
  readDataSetFromGroup(gedges, string("radius"), aflt);
  for(int i=0; i<ecnt; ++i)
  {
    Vessel* v = vl->GetEdge(i);
    v->flags = flags[i];
    v->r = aflt[i];
  }
}

template<class T>
H5::DataSet WriteScalarField(H5::Group &g, const string &name, ConstArray3d<T> arr, const LatticeDataQuad3d &ld, const H5::Group &ldgroup)
{
  H5::DataSet ds = WriteArray3D<T>(g, name, arr);
  writeAttrToH5(ds, string("TYPE"), string("FIELD_QUAD3D"));
  writeAttrToH5(ds, string("LATTICE_PATH"), ldgroup.getObjName());
  return ds;
}

template<class Vector>
H5::DataSet WriteVectorField(H5::Group &g, const string &name, ConstArray3d<Vector> arr, const LatticeDataQuad3d &ld, const H5::Group &ldgroup)
{
  arr = arr[ld.Box()];
  H5::DataSet ds = WriteVectorArray3D<Vector>(g, name, arr);
  //h5cpp::Attributes a = ds.attrs();
  writeAttrToH5(ds, string("TYPE"), string("FIELD_QUAD3D"));
  writeAttrToH5(ds, string("LATTICE_PATH"), ldgroup.getObjName());
//   a.set("TYPE", "FIELD_QUAD3D");
//   a.set("LATTICE_PATH", ldgroup.getObjName());
  return ds;
}


template<class T, int mydim>
H5::DataSet WriteAveragedFaceVariables(H5::Group &file, const string &id, const ConstArray3d<T> *face_fields)
{
  BBox3 bb = face_fields[0].getBox(); --bb.max[0];
  Array3d<Vec<T, mydim> > tmp(bb);
  FOR_BBOX3(p, bb)
  {
    Vec<T, mydim> u;
    for (int axis=0; axis<mydim; ++axis)
    {
      Int3 pp = p; ++pp[axis];
      u[axis] = 0.5 * (face_fields[axis](p) + face_fields[axis](pp));
    }
    tmp(p) = u;
  }
  return WriteVectorArray3D(file, id, tmp[bb]);
}


template<class T>
H5::DataSet WriteAveragedFaceVariables(H5::Group &file, const string &id, int dim, const ConstArray3d<T> *face_fields)
{
  switch(dim)
  {
    case 1: return WriteAveragedFaceVariables<T,1>(file, id, face_fields);
    case 2: return WriteAveragedFaceVariables<T,2>(file, id, face_fields);
    case 3: return WriteAveragedFaceVariables<T,3>(file, id, face_fields);
    default:
      throw std::runtime_error("WriteAveragedFaceVariables: dim must be in [1,3]");
  }
}

template<class T>
H5::DataSet WriteAveragedFaceVariableField(H5::Group &file, const string &id, int dim, const ConstArray3d<T> *face_fields, const LatticeDataQuad3d &ld, const H5::Group &ldgroup)
{
  LatticeIndexRanges ir = LatticeIndexRanges::FromCellRange(ld.Box(), dim);
  ConstArray3d<T> my_fields[3];
  for (int j=0; j<dim; ++j)
    my_fields[j] = face_fields[j][ir.faces[j]];
  H5::DataSet ds = WriteAveragedFaceVariables(file, id, dim, my_fields);
  writeAttrToH5(ds,string("TYPE"), string("FIELD_QUAD3D"));
  writeAttrToH5(ds,string("LATTICE_PATH"), ldgroup.getObjName());
//   H5::Attributes a = ds.attrs();
//   a.set("TYPE", "FIELD_QUAD3D");
//   a.set("LATTICE_PATH", ldgroup.get_name());
  return ds;
}

template<class T>
string getH5GroupName(T g)
{
  size_t len = H5Iget_name(g.getId(),NULL,0);
  char buffer[len];
  H5Iget_name(g.getId(),buffer,len+1);
  std::string n = buffer;
  return n;
}
/* direct instantiation, since we are in a library setting */
template string getH5GroupName<H5::Group>(H5::Group g);
template string getH5GroupName<H5::DataSet>(H5::DataSet g);

// template<typename T>
// T readAttrFromGroup(H5::Group g, string attr_name)
// {
//   if (typeid(T) == typeid(string()))
//   {
//     H5::Attribute myatt_out = g.openAttribute(attr_name);
//     H5::StrType strdatatype(H5::PredType::C_S1, 256); // of length 256 characters
//     H5std_string strreadbuff("");
//     myatt_out.read(strdatatype,strreadbuff);
//     return string(strreadbuff);
//   }
//   if(typeid(T) == typeid(float()))
//   {
//     H5::Attribute myatt_out = g.openAttribute(attr_name);
//     float test = 0.0;
//     H5::DataType type = myatt_out.getDataType();
//     myatt_out.read(type, &test);
//     return test;
//   }
//   if(typeid(T) == typeid(double()))
//   {
//     H5::Attribute myatt_out = g.openAttribute(attr_name);
//     double test = 0.0;
//     H5::DataType type = myatt_out.getDataType();
//     myatt_out.read(type, &test);
//     return test;
//   }
//   
//   if(typeid(T) == typeid(Float3))
//   {
//     try
//     {
//        /*
//         * Turn off the auto-printing when failure occurs so that we can
//         * handle the errors appropriately
//         */
//       H5::Exception::dontPrint();
//       /*
//        * Open attribute space
//        */
//       H5::Attribute myatt_out = g.openAttribute(attr_name);
//       H5::DataSpace attributeFileSpace = myatt_out.getSpace();
//       H5::DataType attributeFileType = myatt_out.getDataType();
//       //int rank = attributeSpace.getSimpleExtentNdims();
//       /*
//        * Get the dimensions sizes of the file dataspace
//        */
//       hsize_t dims[2];
//       int rank = attributeFileSpace.getSimpleExtentDims(dims);
//       myAssert(dims[1] == 3); // need to be a Float3
//       
//       /* 
//        * allocate memory space to read dataset
//        */
//       H5::DataSpace mspace( rank, dims);
// 
//       Float3 avec;
//       //myatt_out.read(&avec, H5::PredType::NATIVE_FLOAT, mspace, attributeFileSpace);
//       myatt_out.read(attributeFileType,&avec);
//     return avec;
//       
//     }
//     // catch failure caused by the H5File operations
//     catch( H5::FileIException error )
//     {
//       error.printErrorStack();
//     }
//     // catch failure caused by the DataSet operations
//     catch( H5::DataSetIException error )
//     {
//       error.printErrorStack();
//     }
//     // catch failure caused by the DataSpace operations
//     catch( H5::DataSpaceIException error )
//     {
//       error.printErrorStack();
//     }
//     
//   }//end if
// }


// float readAttrFromGroup(H5::Group g, string attr_name)
// {
//   {
//     H5::Attribute myatt_out = g.openAttribute(attr_name);
//     float test = 0.0;
//     H5::DataType type = myatt_out.getDataType();
//     myatt_out.read(type, &test);
//     return test;
//   } 
// }

template <class T>
void readDataSetFromGroup(H5::Group &g, string ds_name, T *placeToStore)
{
  H5::DataSet dset = g.openDataSet(ds_name);
  H5::DataSpace dspace = dset.getSpace();
  /*
  * Get the dimensions sizes of the file dataspace
  */
  auto sizearray = placeToStore.size();
  int rank_of_cpp = sizearray.size();
  hsize_t dims[rank_of_cpp];
  
  int rank = dspace.getSimpleExtentDims(dims);
  myAssert( rank_of_cpp == rank );
  //void read( void* buf, const DataType& mem_type, const DataSpace& mem_space = DataSpace::ALL, const DataSpace& file_space = DataSpace::ALL, const DSetMemXferPropList& xfer_plist = DSetMemXferPropList::DEFAULT ) const;
  dset.read(placeToStore, H5::PredType::NATIVE_FLOAT);
}




template<class T>
H5::DataSet writeDataSetToGroup(H5::Group &g, const string &dataset_name, DynArray<T> &value)
{
  int sizeOfDynArray = value.size();
  T continousMemoryArrary[sizeOfDynArray];
  for( int i=0; i<sizeOfDynArray;++i)
  {
    continousMemoryArrary[i] = T(value[i]);
  }
  int rank = 2;
  hsize_t dims[rank];
  dims[0]=sizeOfDynArray;
  H5::DataType thisWritingType = getH5TypeFromCpp<T>();
  if(typeid(T) == typeid(Float3) or typeid(T) == typeid(Int3) or typeid(T) == typeid(Bool3))
  {
    dims[1] = 3;
  }
  else if( typeid(T) == typeid(Float2))
  {
    dims[1] = 2;
  }
  else
  {
    dims[1] = 1;
  }
#ifndef NDEBUG
  cout << "we are writting data of size (" << dims[0] << ", " << dims[1] << "!" <<endl;
#endif
  H5::DataSet ds;
  try 
  {
    H5::DataSpace dataspace( rank, dims);
    ds = g.createDataSet(dataset_name, thisWritingType, dataspace);
    ds.write(&continousMemoryArrary, thisWritingType, dataspace);
  }
  catch( H5::Exception e)
  {
    e.printErrorStack();
  }
  return ds;
}

template <class T>
void readDataSetFromGroup(H5::Group &g, const string &dataset_name, DynArray<T> &readIn)
{
  H5::DataSet dset;
  H5::DataSpace dataspace;
  
  try{
    dset = g.openDataSet(dataset_name);
    dataspace = dset.getSpace();
    hsize_t rank = dataspace.getSimpleExtentNdims();
    hsize_t dims[rank];
    hsize_t max_dims[rank];
    hsize_t again = dataspace.getSimpleExtentDims(dims,max_dims);
#ifndef NDEBUG
    cout << "rank: " << rank << endl;
    for(int i=0; i<rank; i++)
    {
      cout << "dims[" << i << "]" << " = " << dims[i] << endl;
    }
#endif
  
    boost::multi_array<T,1> arr_data(boost::extents[dims[0]]);
    H5::DataType thisWritingType = getH5TypeFromCpp<T>();

    dset.read(arr_data.data(), thisWritingType);
    readIn.resize(dims[0]);
    
    for( int i=0;i<dims[0];++i)
    {
      readIn[i] = arr_data[i];
    }
  }
  catch(H5::Exception e)
  {
    e.printErrorStack();
  }
}

////// from hdf_wrapper 
// template<class LD>
// void WriteHdfLdGenericPart_(H5::Group g, const LD &ld)
// {
//   writeAttrToH5(g,"SIZEX", ld.Size()[0]);
//   writeAttrToH5(g,"SIZEY", ld.Size()[1]);
//   writeAttrToH5(g,"SIZEZ", ld.Size()[2]);
//   writeAttrToH5(g, "SIZE", ld.Size());
//   //H5::Attribute attrs = g.attrs();
// //   attrs.set("SIZEX",ld.Size()[0]);
// //   attrs.set("SIZEY",ld.Size()[1]);
// //   attrs.set("SIZEZ",ld.Size()[2]);
// //   set_array(attrs, "SIZE", ld.Size());
//   BBox3 bb = ld.Box();
//   // box is stored in a 6 component vector, xxyyzz, must match python code
//   Int6 bv;
//   //Int6 bv;
//   for (int i=0; i<3; ++i)
//   {
//     bv[i*2  ] = bb.min[i];
//     bv[i*2+1] = bb.max[i];
//   }
//   writeAttrToH5(g, string("BOX"), bv);
//   writeAttrToH5(g, string("WORLD_OFFSET"),ld.GetOriginPosition() );
//   //set_array(attrs, "BOX", bv);
//   //attrs.set("SCALE",ld.Scale());
//   //set_array(attrs, "WORLD_OFFSET", ld.GetOriginPosition());
// }

// void WriteHdfLd( H5::Group &g, const LatticeDataQuad3d &ld )
// {
//   //writeAttrToH5<Bool3>(g, "CENTERING", ld.GetCellCentering().cast<int>());
//   Bool3 centering = ld.GetCellCentering();
//   writeAttrToH5(g, "CENTERING", centering);
//   writeAttrToH5(g, "TYPE", string("QUAD3D"));
// //   h5cpp::Attributes attrs = g.attrs();
// //   attrs.set("TYPE","QUAD3D");
//   //auto buffer = ld.GetCellCentering().cast<int>();
//   //set_array<int,3>(attrs, "CENTERING", ld.GetCellCentering().cast<int>());
//   WriteHdfLdGenericPart_(g, ld);
// }

// void WriteHdfLd( H5::Group &g, const LatticeDataFCC &ld )
// {
// //   h5cpp::Attributes attrs = g.attrs();
// //   attrs.set("TYPE","FCC");
//   writeAttrToH5(g, "TYPE", string("FCC"));
//   WriteHdfLdGenericPart_(g, ld);
// }




template<class LD>
void ReadHdfLdGenericPart_(H5::Group &g, LD &ld)
{
  //h5cpp::Attributes attrs = g.attrs();
  BBox3 bb; float scale;
  try
  {
    // box is stored in a 6 component vector, xxyyzz
    Int6 bv;
    readAttrFromH5(g, string("BOX"), bv);
    //Vec<int, 6> bv = get_array<int, 6>(attrs, "BOX");
    for (int i=0; i<3; ++i)
    {
      // this must match the python code!!
      bb.min[i] = bv[i*2  ];
      bb.max[i] = bv[i*2+1];
    }
  }
  catch (const H5::AttributeIException &e)
  {
    // legacy code :[
    readAttrFromH5(g,string("SIZEX"),bb.max[0]);
    readAttrFromH5(g,string("SIZEY"),bb.max[1]);
    readAttrFromH5(g,string("SIZEZ"),bb.max[2]);
//     attrs.get("SIZEX", bb.max[0]);
//     attrs.get("SIZEY", bb.max[1]);
//     attrs.get("SIZEZ", bb.max[2]);
    for (int i=0; i<3; ++i)
    {
      bb.min[i] = 0;
      bb.max[i] -= 1;
    }
  }
  readAttrFromH5(g, string("SCALE"), scale);
  //attrs.get("SCALE", scale);
  ld.Init(bb, (double)scale);
  try 
  {
    Float3 world_offset;
    readAttrFromH5(g, "WORLD_OFFSET", world_offset);
    //ld.SetOriginPosition(get_array<float,3>(attrs, "WORLD_OFFSET"));
    ld.SetOriginPosition(world_offset);
  } 
  catch (const H5::AttributeIException &e) {}
}


// void ReadHdfLd( h5cpp::Group g, LatticeDataQuad3d &ld )
// {
//   h5cpp::Attributes attrs = g.attrs();
//   if(attrs.get<string>("TYPE")!="QUAD3D") throw std::runtime_error("LatticeDataQuad3d from hdf5 mismatch");
//   ReadHdfLdGenericPart_(g, ld);
//   try {
//     ld.SetCellCentering(get_array<int,3>(attrs, "CENTERING").cast<bool>());
//   } catch (const h5cpp::NameLookupError &e) {}
// }
// 
// 
// void ReadHdfLd( h5cpp::Group g, LatticeDataFCC &ld )
// {
//   h5cpp::Attributes attrs = g.attrs();
//   if(attrs.get<string>("TYPE")!="FCC") throw std::runtime_error("LatticeDataFCC from hdf5 mismatch");
//   ReadHdfLdGenericPart_(g, ld);
// }

void ReadHdfLd( H5::Group &g, LatticeDataQuad3d &ld )
{
  //h5cpp::Attributes attrs = g.attrs();
  
  string type;
  readAttrFromH5(g, string("TYPE"), type);
  if(type!=string("QUAD3D")) throw std::runtime_error("LatticeDataQuad3d from hdf5 mismatch");
  ReadHdfLdGenericPart_(g, ld);
  //do this properly!!!
//   try {
//     ld.SetCellCentering(get_array<int,3>(attrs, "CENTERING").cast<bool>());
//   } catch (const h5cpp::NameLookupError &e) {}
}


void ReadHdfLd( H5::Group &g, LatticeDataFCC &ld )
{
  //h5cpp::Attributes attrs = g.attrs();
  // __cxx1112basic_stringIcSt11char_traitsIcESaIcEEEEvT_RKS7_RT0_
  // __cxx1112basic_stringIcSt11char_traitsIcESaIcEEERT0_
  // IN2H55GroupENSt7
  // IN2H55GroupEiEvT_RKNSt7
  // WARNING HERE IS A BUG! --> partly solved
  // maybe from inherence of differnt string types e.g. H5string std::string etc...
  string type;
  readAttrFromH5(g, string("TYPE"), type);
  if(type!="FCC") throw std::runtime_error("LatticeDataFCC from hdf5 mismatch");
  ReadHdfLdGenericPart_(g, ld);
}
////// end hdf_wrapper 
#define INSTANTIATE(T)\
  template H5::DataSet WriteScalarField<T>(H5::Group &g, const string &name, ConstArray3d<T> arr, const LatticeDataQuad3d &ld, const H5::Group &ldgroup);
INSTANTIATE(float)
INSTANTIATE(double)
INSTANTIATE(int)
INSTANTIATE(char)
#undef INSTANTIATE

#define INSTANTIATE2(T)\
  template H5::DataSet writeDataSetToGroup<T>(H5::Group &g, const string &name, DynArray<T> &value);
INSTANTIATE2(float)
INSTANTIATE2(double)
INSTANTIATE2(int)
INSTANTIATE2(bool)
INSTANTIATE2(Float3)
INSTANTIATE2(Float2)
INSTANTIATE2(Int3)
INSTANTIATE2(Bool3)
INSTANTIATE2(char)
INSTANTIATE2(uchar)
INSTANTIATE2(long)
INSTANTIATE2(string)

#undef INSTANTIATE2

// #define INSTANTIATE_H5Cpp_read(U,T)\
//   template void readAttrFromH5<U, T>(U &g, const string &name, T &output_buffer);
// INSTANTIATE_H5Cpp_read(H5::Group, float)
// INSTANTIATE_H5Cpp_read(H5::Group, Float3)
// INSTANTIATE_H5Cpp_read(H5::Group, double)
// INSTANTIATE_H5Cpp_read(H5::Group, Double3)
// INSTANTIATE_H5Cpp_read(H5::Group, int)
// INSTANTIATE_H5Cpp_read(H5::Group, Int3)
// INSTANTIATE_H5Cpp_read(H5::Group, Int6)
// INSTANTIATE_H5Cpp_read(H5::Group, bool)
// INSTANTIATE_H5Cpp_read(H5::Group, uchar)
// INSTANTIATE_H5Cpp_read(H5::Group, Bool3)
// 
// 
// INSTANTIATE_H5Cpp_read(H5::DataSet, float)
// INSTANTIATE_H5Cpp_read(H5::DataSet, Float3)
// INSTANTIATE_H5Cpp_read(H5::DataSet, double)
// INSTANTIATE_H5Cpp_read(H5::DataSet, Double3)
// INSTANTIATE_H5Cpp_read(H5::DataSet, int)
// INSTANTIATE_H5Cpp_read(H5::DataSet, Int3)
// INSTANTIATE_H5Cpp_read(H5::DataSet, Int6)
// INSTANTIATE_H5Cpp_read(H5::DataSet, bool)
// INSTANTIATE_H5Cpp_read(H5::DataSet, uchar)
// INSTANTIATE_H5Cpp_read(H5::DataSet, Bool3)
// #undef INSTANTIATE_H5Cpp_read

#define INSTANTIATE_H5Cpp_read(T)\
  template void readAttrFromH5<T>(H5::H5Location &g, const string &name, T &output_buffer);
INSTANTIATE_H5Cpp_read(float)
INSTANTIATE_H5Cpp_read(Float3)
INSTANTIATE_H5Cpp_read(double)
INSTANTIATE_H5Cpp_read(Double3)
INSTANTIATE_H5Cpp_read(int)
INSTANTIATE_H5Cpp_read(Int3)
INSTANTIATE_H5Cpp_read(Int6)
INSTANTIATE_H5Cpp_read(bool)
INSTANTIATE_H5Cpp_read(uchar)
INSTANTIATE_H5Cpp_read(Bool3)
#undef INSTANTIATE_H5Cpp_read

// #define INSTANTIATE_H5Cpp1_write(U,T)\
//   template void writeAttrToH5<U,T>(U &h, const string &name, const T &output_buffer);
// INSTANTIATE_H5Cpp1_write(H5::Group, float)
// INSTANTIATE_H5Cpp1_write(H5::Group, Float3)
// INSTANTIATE_H5Cpp1_write(H5::Group, double)
// INSTANTIATE_H5Cpp1_write(H5::Group, Double3)
// INSTANTIATE_H5Cpp1_write(H5::Group, int)
// INSTANTIATE_H5Cpp1_write(H5::Group, Int3)
// INSTANTIATE_H5Cpp1_write(H5::Group, Int6)
// INSTANTIATE_H5Cpp1_write(H5::Group, bool)
// INSTANTIATE_H5Cpp1_write(H5::Group, uchar)
// INSTANTIATE_H5Cpp1_write(H5::Group, Bool3)
// INSTANTIATE_H5Cpp1_write(H5::Group, unsigned long)
// //INSTANTIATE_H5Cpp1_write(H5::Group, string)
// 
// 
// INSTANTIATE_H5Cpp1_write(H5::DataSet, float)
// INSTANTIATE_H5Cpp1_write(H5::DataSet, Float3)
// INSTANTIATE_H5Cpp1_write(H5::DataSet, double)
// INSTANTIATE_H5Cpp1_write(H5::DataSet, Double3)
// INSTANTIATE_H5Cpp1_write(H5::DataSet, int)
// INSTANTIATE_H5Cpp1_write(H5::DataSet, Int3)
// INSTANTIATE_H5Cpp1_write(H5::DataSet, Int6)
// INSTANTIATE_H5Cpp1_write(H5::DataSet, bool)
// INSTANTIATE_H5Cpp1_write(H5::DataSet, uchar)
// INSTANTIATE_H5Cpp1_write(H5::DataSet, Bool3)
// INSTANTIATE_H5Cpp1_write(H5::DataSet, unsigned long)
// //INSTANTIATE_H5Cpp1_write(H5::DataSet, string)
// #undef INSTANTIATE_H5Cpp1_write


#define INSTANTIATE_H5Cpp1_write(T)\
  template void writeAttrToH5<T>(H5::H5Location &h, const string &name, const T &output_buffer);
INSTANTIATE_H5Cpp1_write(float)
INSTANTIATE_H5Cpp1_write(Float3)
INSTANTIATE_H5Cpp1_write(double)
INSTANTIATE_H5Cpp1_write(Double3)
INSTANTIATE_H5Cpp1_write(int)
INSTANTIATE_H5Cpp1_write(Int3)
INSTANTIATE_H5Cpp1_write(Int6)
INSTANTIATE_H5Cpp1_write(bool)
INSTANTIATE_H5Cpp1_write(uchar)
INSTANTIATE_H5Cpp1_write(Bool3)
INSTANTIATE_H5Cpp1_write(unsigned long)
//INSTANTIATE_H5Cpp1_write(H5::Group, string)
#undef INSTANTIATE_H5Cpp1_write

#define INSTANTIATE_VEC(T)\
  template H5::DataSet WriteAveragedFaceVariableField<T>(H5::Group &file, const string &id, int dim, const ConstArray3d<T> *face_fields, const LatticeDataQuad3d &ld, const H5::Group &ldgroup);\
  template H5::DataSet WriteVectorField<Vec<T,3> >(H5::Group &g, const string &name, ConstArray3d<Vec<T,3> > arr, const LatticeDataQuad3d &ld, const H5::Group &ldgroup);\
  template H5::DataSet WriteVectorField<Vec<T,2> >(H5::Group &g, const string &name, ConstArray3d<Vec<T,2> > arr, const LatticeDataQuad3d &ld, const H5::Group &ldgroup);\
  template H5::DataSet WriteVectorField<Vec<T,1> >(H5::Group &g, const string &name, ConstArray3d<Vec<T,1> > arr, const LatticeDataQuad3d &ld, const H5::Group &ldgroup);\
  template H5::DataSet WriteAveragedFaceVariables<T>(H5::Group &file, const string &id, int dim, const ConstArray3d<T> *face_fields);

INSTANTIATE_VEC(float)
INSTANTIATE_VEC(double)
#undef INSTANTIATE_VEC

/** from hdf_wrapper_array3d.cpp
 */


//arr3d.data = malloc(arr3d.a * arr3d.b * arr3d.c * sizeof *arr3d.data);

//arr3d[r][c][d]
// becomes:
//arr3d.data[r * (arr3d.b * arr3d.c) + c * arr3d.c + d];

template<class T>
H5::DataSet WriteArray3D(H5::Group &file, const std::string &DATASET_NAME, const ConstArray3d<T> &a)
{
#ifdef NDEBUG
  //std::cout << " in write Array 3D" << std::endl;std::cout.flush();
#endif
  const Int3 s = a.size();
  const int rank = 3;
  hsize_t     dims[rank];       // dataset dimensions
  dims[0] = s[0];
  dims[1] = s[1];
  dims[2] = s[2];
  // allocate array using the "extents" helper. 
  // This makes it easier to see how big the array is
  // this works fine with hdf5
  boost::multi_array<T, 3>  arr_3d_data(boost::extents[s[0]][s[1]][s[2]]);
  
  H5::DataSpace dspace = H5::DataSpace( rank, dims);
  H5::DataSet dataset = file.createDataSet(DATASET_NAME, H5::PredType::NATIVE_FLOAT,dspace);
 
  /** 
   * need to allocate array on heap, stack might be too small
   * I messed up the heap allocation, by hand. 
   * Now I delegate the messy work to boost.
   */
//   T ***ptr3D = NULL;
//   
//   ptr3D = new T**[s[0]];
//   for(int i = 0; i< s[0];i++)
//   {
//     ptr3D[i] = new T*[s[1]];
//     for(int ii=0; ii< s[1]; ii++)
//     {
//       ptr3D[i][ii]= new T[s[2]];
//       
//       for(int iii=0; iii<s[2]; iii++)
//       {
//         ptr3D[i][ii][iii] = a(i,ii,iii);
//       }
//     }
//   }

  for(int i = 0; i< s[0];i++)
    for(int ii=0; ii< s[1]; ii++)
      for(int iii=0; iii<s[2]; iii++)
        arr_3d_data[i][ii][iii] = a(i,ii,iii);
    
  // Write the data to the dataset using default memory space, file
	// space, and transfer properties.
  try{
    dataset.write(arr_3d_data.data(), H5::PredType::NATIVE_FLOAT);
  }
  catch(H5::Exception e)
  {
    e.printErrorStack();
  }
  return dataset;
}


template<class Vector>
H5::DataSet WriteVectorArray3D(H5::Group  &file,const std::string &id, const ConstArray3d<Vector> &a)
{
  const Int3 s = a.size();
  const int rank = 4;
  hsize_t     dims[rank];       // dataset dimensions
  dims[0] = s[0];
  dims[1] = s[1];
  dims[2] = s[2];
  dims[3] = 3;
  std::cout << " Warning: untested!" << std::endl;
  boost::multi_array<typename Vector::Scalar, 4>  arr_3d_data(boost::extents[s[0]][s[1]][s[2]][3]);
  H5::DataSpace dspace = H5::DataSpace( rank, dims);
  H5::DataSet dataset = file.createDataSet(id, H5::PredType::NATIVE_FLOAT,dspace);
  
  for(int i = 0; i< s[0];i++)
    for(int ii=0; ii< s[1]; ii++)
      for(int iii=0; iii<s[2]; iii++)
        for(int iiii=0; iiii<3; iiii++)
          arr_3d_data[i][ii][iii][iiii] = a(i,ii,iii)[iiii];
	

  try{
    //dataset.write(&tmp, H5::PredType::NATIVE_FLOAT);
    dataset.write(arr_3d_data.data(), H5::PredType::NATIVE_FLOAT);
  }
  catch(H5::Exception e)
  {
    e.printErrorStack();
  }
  
  return dataset;
}



template<class T>
void ReadArray3D(H5::DataSet &ds, Array3d<T> &a)
{
//   h5cpp::Dataspace sp = ds.get_dataspace();
//   hsize_t dims[H5S_MAX_RANK];
//   int r = sp.get_dims(dims);
//   if (r != 3) throw h5cpp::Exception("ReadArray3d expected Dataset of Rank 3!");
//   Array3d<T> tmp(Int3(dims[2], dims[1], dims[0]));
//   tmp.swapAxes(0,2);
//   ds.read<>(tmp.getPtr());
//   a = Array3d<T>(tmp.size());
//   a.fill(tmp);
}


template<class Vector>
void ReadVectorArray3D(H5::DataSet &ds, Array3d<Vector> &a)
{
//   h5cpp::Dataspace sp = ds.get_dataspace();
//   hsize_t dims[H5S_MAX_RANK];
//   int r = sp.get_dims(dims);
//   if (r != 4) throw h5cpp::Exception("ReadVectorArray3d expected Dataset of Rank 4.");
//   if (dims[3] != Vector::SizeAtCompileTime) throw h5cpp::Exception("ReadVectorArray3d expected Dataset dimension 3 to be the same size as the Vector type.");
//   Array3d<Vector> tmp(Int3(dims[2], dims[1], dims[0]));
//   tmp.swapAxes(0,2);
//   ds.read<>((typename Vector::value_type*)tmp.getPtr());
//   a = Array3d<Vector>(tmp.size());
//   a.fill(tmp);
}

/*
#define INSTANTIATE(T)\
  template h5cpp::Dataset WriteArray3D<T>(h5cpp::Group file, const std::string &id, const ConstArray3d<T> &a, const h5cpp::Datatype &disktype);\
  template void ReadArray3D<T>(h5cpp::Dataset ds, Array3d<T> &a);

#define INSTANTIATE_VEC(T)\
  template h5cpp::Dataset WriteVectorArray3D<Vec<T,3> >(h5cpp::Group  file,const std::string &id, const ConstArray3d<Vec<T,3> > &a);\
  template h5cpp::Dataset WriteVectorArray3D<Vec<T,2> >(h5cpp::Group  file,const std::string &id, const ConstArray3d<Vec<T,2> > &a);\
  template h5cpp::Dataset WriteVectorArray3D<Vec<T,1> >(h5cpp::Group  file,const std::string &id, const ConstArray3d<Vec<T,1> > &a);\
  template void ReadVectorArray3D<Vec<T,3> >(h5cpp::Dataset ds, Array3d<Vec<T,3> > &a);

INSTANTIATE(float)
INSTANTIATE(double)
INSTANTIATE(int)
INSTANTIATE(char)
INSTANTIATE_VEC(float)
INSTANTIATE_VEC(double)
*/

#define INSTANTIATE(T)\
  template H5::DataSet WriteArray3D<T>(H5::Group &file, const std::string &id, const ConstArray3d<T> &a);\
  template void ReadArray3D<T>(H5::DataSet &ds, Array3d<T> &a);

#define INSTANTIATE_VEC(T)\
  template H5::DataSet WriteVectorArray3D<Vec<T,3> >(H5::Group  &file,const std::string &id, const ConstArray3d<Vec<T,3> > &a);\
  template H5::DataSet WriteVectorArray3D<Vec<T,2> >(H5::Group  &file,const std::string &id, const ConstArray3d<Vec<T,2> > &a);\
  template H5::DataSet WriteVectorArray3D<Vec<T,1> >(H5::Group  &file,const std::string &id, const ConstArray3d<Vec<T,1> > &a);\
  template void ReadVectorArray3D<Vec<T,3> >(H5::DataSet &ds, Array3d<Vec<T,3> > &a);

INSTANTIATE(float)
INSTANTIATE(double)
INSTANTIATE(int)
INSTANTIATE(char)
INSTANTIATE_VEC(float)
INSTANTIATE_VEC(double)
