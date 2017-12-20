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

#include <stdexcept>

H5::Group RequireLatticeDataGroup(H5::H5File f, const string name, const LatticeDataQuad3d &ld)
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
      WriteHdfLd(ld_group, ld);
  }
//   H5::Group ld_group = g.createGroup()
//   h5cpp::Group ld_group = g.require_group(name, &had_ld_group);
//   if (!had_ld_group)
//     WriteHdfLd(ld_group, ld);
  return ld_group;
}

H5::Group RequireLatticeDataGroup(H5::Group g, const LatticeDataQuad3d &ld)
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
      WriteHdfLd(ld_group, ld);
  }
//   H5::Group ld_group = g.createGroup()
//   h5cpp::Group ld_group = g.require_group(name, &had_ld_group);
//   if (!had_ld_group)
//     WriteHdfLd(ld_group, ld);
  return ld_group;
}



void WriteHdfGraph( H5::Group g, const VesselList3d &vl )
{
#ifdef DEBUG
  printf("Starting to write hdf\n");
#endif
  const int ncnt = vl.GetNCount();
  const int ecnt = vl.GetECount();
  myAssert(ncnt>0 && ecnt>0);
  
  H5::Group gg = g.createGroup("nodes");
  writeAttrToGroup<int>(gg, "COUNT", ncnt);
  //gg.attrs().set("COUNT",ncnt);
  //h5cpp::Attributes attrs = g.attrs();
//   H5::Group gg = g.createGroup("nodes");
//   gg.createAttribute("COUNT");
//   gg.attrs().set("COUNT",ncnt);
//   h5cpp::Attributes attrs = g.attrs();
  if(vl.HasLattice())//LATTICE IS THERE
  {
    try
    {
      string theType = readAttrFromGroup<string>(g, string("CLASS"));
      myAssert(theType == "GRAPH");
    }
    catch( H5::AttributeIException not_found_error )
    {
        cout << " Graph type not found" << endl;
	writeAttrToGroup<string>(g,string("CLASS"), string("GRAPH"));
    }
//     if(attrs.exists("CLASS"))//correct hdf attribute is checke
//     {
//       myAssert(attrs.get<std::string>("CLASS") == "GRAPH");
//     }
//     else
//     {
//       attrs.create<std::string>("CLASS", "GRAPH");
//     }
    //lattice stuff is writen to hdf
    H5::Group lattice_group = g.createGroup("lattice");
    //vl.Ld().WriteHdf(g.create_group("lattice"));
    vl.Ld().WriteHdf(lattice_group);
    const VesselList3d::LatticeData &ld = vl.Ld();
    DynArray<int> a(ncnt);
    for(int i=0; i<ncnt; ++i) 
    {
      a[i] = ld.LatticeToSite(vl.GetNode(i)->lpos);
    }
    //h5cpp::Dataset ds = h5cpp::create_dataset<int>(gg,  "lattice_pos", a);
    hsize_t dims[2] = {1,ncnt};
    int RANK = 2;
    H5::DataSpace dataspace = H5::DataSpace(RANK, dims); // create new dspace
    H5::DSetCreatPropList ds_creatplist;
    H5::DataSet dataset = H5::DataSet(gg.createDataSet("lattice_pos", H5::PredType::NATIVE_INT,dataspace, ds_creatplist ));
    dataset.write(&a, H5::PredType::NATIVE_INT, dataspace);
    writeAttrToDataset<string>(dataset, "MODE", "linear");
    //ds.attrs().set("MODE","linear");
  }
  else//no lattice pressent
  {
    try
    {
      string theType = readAttrFromGroup<string>(g, "CLASS");
      myAssert(theType == "REALWORLD");
    }
    catch( H5::AttributeIException not_found_error )
    {
        cout << " CLASS type not found" << endl;
	writeAttrToGroup<string>(g,"CLASS", string("REALWORLD"));
    }  
//     if(attrs.exists("CLASS"))
//     {
//       myAssert(attrs.get<std::string>("CLASS") == "REALWORLD");
//     }
//     else
//     {
//       attrs.create<std::string>("CLASS", "REALWORLD");
//     }
    // write the special world stuff
    //needs old school linearizing
    DynArray<float> a(3*ncnt);
    for(int i=0; i<ncnt; ++i) 
    {
      a[3*i+0]= vl.GetNode(i)->worldpos[0];
      a[3*i+1] = vl.GetNode(i)->worldpos[1];
      a[3*i+2] = vl.GetNode(i)->worldpos[2];
    }
    // due to creation of h5 files with python this could already be pressent
    
    try{
      //DataSet dataset = file.openDataSet(DATASET_NAME);
      H5::DataSet ds = gg.openDataSet("world_pos");
      //grp=h5file.openGroup("A");
    } catch(H5::Exception& e){
      /* group does not exists, create it */
      //grp=h5file.createGroup("A");
      const int	 RANK = 2;
      hsize_t dims[2] = {3, ncnt};
      H5::DataSpace dataspace(RANK, dims);
      H5::DataSet ds = gg.createDataSet("world_pos", H5::PredType::NATIVE_FLOAT, dataspace);
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
  try{
      //DataSet dataset = file.openDataSet(DATASET_NAME);
      H5::DataSet ds = gg.openDataSet("roots");
      //grp=h5file.openGroup("A");
    } catch(H5::Exception& e){
      /* group does not exists, create it */
      //grp=h5file.createGroup("A");
      const int	 RANK = 2;
      hsize_t dims[2] = {1, 16};
      H5::DataSpace dataspace(RANK, dims);
      H5::DataSet ds = gg.createDataSet("roots", H5::PredType::NATIVE_INT, dataspace);
      ds.write(&roots, H5::PredType::NATIVE_INT);
      writeDataSetToGroup(gg,"roots", roots);
    }
//   if(!gg.exists("roots"))
//   {
//     h5cpp::create_dataset<int>(gg,  "roots", roots);
//   }
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
    writeDataSetToGroup<DynArray<int>>(gg, string("bc_node_index"), bc_node_index);
    writeDataSetToGroup<DynArray<int>>(gg, string("bc_type"), bctyp_index);
    writeDataSetToGroup<DynArray<float>>(gg, string("bc_value"), values_of_bcs);
    writeDataSetToGroup<DynArray<float>>(gg, string("bc_conductivity_value"), bc_conductivity_value);
//     h5cpp::create_dataset<int>(gg,  "bc_node_index", bc_node_index);
//     h5cpp::create_dataset<int>(gg,  "bc_type", bctyp_index);
//     h5cpp::create_dataset<float>(gg,  "bc_value", values_of_bcs);
//     h5cpp::create_dataset<float>(gg,  "bc_conductivity_value", bc_conductivity_value);
  }
  }//for interrupt 

  gg = g.createGroup("edges");
  writeAttrToGroup<int>(gg,string("COUNT"), ecnt);
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
    writeDataSetToGroup<DynArray<int>>(gg, string("node_a_index"), va);
    writeDataSetToGroup<DynArray<int>>(gg, string("node_b_index"), vb);
    writeDataSetToGroup<DynArray<int>>(gg, string("flags"), flags);
    writeDataSetToGroup<DynArray<float>>(gg, string("radius"), float_radius);
    writeAttrToGroup<string>(gg, string("MODE"), string("const"));
    
   
//     h5cpp::create_dataset<int>( gg, "node_a_index", va);
//     h5cpp::create_dataset<int>( gg, "node_b_index", vb );
//     h5cpp::create_dataset<int>( gg, "flags", flags );
//     h5cpp::Dataset ds = h5cpp::create_dataset<float>( gg, "radius", float_radius );
//     ds.attrs().set("MODE","const");
  }
}

void ReadHdfGraph( H5::Group g, VesselList3d *vl )
{
  H5::Group gnodes = g.openGroup("nodes");
  H5::Group gedges = g.openGroup("edges");
  
  int ecnt=0,ncnt=0;
  //H5::StrType strdatatype(H5::PredType::C_S1, 256); // of length 256 characters
  //H5std_string strreadbuff("");
  H5::Attribute att_gnodes = gnodes.openAttribute("COUNT");
  H5::IntType intdatatype(H5::PredType::NATIVE_INT); // of length 256 characters
  att_gnodes.read(intdatatype,&ncnt);
  //gedges.openAttribute("COUNT").read(H5::IntType, &ecnt);
  //gnodes.attrs().get("COUNT",ncnt);
  //gedges.attrs().get("COUNT",ecnt);

#if 0
  if(vl->HasLattice())//LATTICE IS THERE
  {//read lattice stuffe
    const VesselList3d::LatticeData &ld = vl->Ld();
    {//begin interrupt
      
      //node stuff
      std::vector<VesselList3d::SiteType> a;
      h5cpp::read_dataset<VesselList3d::SiteType>(gnodes.open_dataset("lattice_pos"),a);
      for(int i=0; i<ncnt; ++i)
      {
	vl->InsertNode(ld.SiteToLattice(a[i]));
      }
    }//end interrupt
  }
  else//no lattice pressent
  {
    {//begin interupt
      DynArray<float> a;
      h5cpp::Dataset ds = gnodes.open_dataset("world_pos");
      h5cpp::Dataspace sp = ds.get_dataspace();
      hsize_t dims[H5S_MAX_RANK];
      int rank = sp.get_dims(dims);
      h5cpp::read_dataset<float>(ds,a);
      for(int i=0; i<ncnt; ++i)
      {
	vl->InsertNode(
	  Float3(
	    a[3*i+0],a[3*i+1],a[3*i+2]
	  )
	);
      }
    }//end interupt
  }
  /* Read all other stuff which is the same*/
  {
    //edge stuff
    DynArray<int> va,vb;
    h5cpp::read_dataset<int>(gedges.open_dataset("node_a_index"),va);
    h5cpp::read_dataset<int>(gedges.open_dataset("node_b_index"),vb);
    for(int i=0; i<ecnt; ++i)
    {
      Vessel* v = vl->InsertVessel(vl->GetNode(va[i]),vl->GetNode(vb[i]));
    }
  }
  {
    DynArray<int> root_indices;
    h5cpp::read_dataset<int>(gnodes.open_dataset("roots"),root_indices);
    for(int i=0; i<root_indices.size(); ++i)
      vl->GetNode(root_indices[i])->flags.AddBits(BOUNDARY);
  }
    
  if(gnodes.exists("bc_node_index") && 
    gnodes.exists("bc_type") && 
    gnodes.exists("bc_value") && 
    gnodes.exists("bc_conductivity_value"))
  {
    DynArray <int> bc_node_index;
    DynArray <int> bctyp_index;
    DynArray <float> values_of_bcs;
    DynArray <float> bc_conductivity_value;
    h5cpp::read_dataset<int>(gnodes.open_dataset("bc_node_index")     ,bc_node_index);
    h5cpp::read_dataset<int>(gnodes.open_dataset("bc_type")           ,bctyp_index);
    h5cpp::read_dataset<float>(gnodes.open_dataset("bc_value")          ,values_of_bcs);
    h5cpp::read_dataset<float>(gnodes.open_dataset("bc_conductivity_value") ,bc_conductivity_value);
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
      //vl.GetBCMap().emplace(nd,bc);
    }
  }
  
  DynArray<int> flags;
  DynArray<float> aflt;
  h5cpp::read_dataset<int>(gedges.open_dataset("flags"),flags);
  h5cpp::read_dataset<float>(gedges.open_dataset("radius"),aflt);
  for(int i=0; i<ecnt; ++i)
  {
    Vessel* v = vl->GetEdge(i);
    v->flags = flags[i];
    v->r = aflt[i];
  }
#endif
}

template<class T>
H5::DataSet WriteScalarField(H5::Group g, const string &name, ConstArray3d<T> arr, const LatticeDataQuad3d &ld, const H5::Group ldgroup)
{
  arr = arr[ld.Box()];
  //h5cpp::Dataset ds = WriteArray3D<T>(g, name, arr, disktype);
  //H5::DataSet ds
  H5::DataSet ds = WriteArray3D<T>(g, name, arr);
//   h5cpp::Attributes a = ds.attrs();
  writeAttrToDataset<string>(ds, string("TYPE"), string("FIELD_QUAD3D"));
  writeAttrToDataset<string>(ds, string("LATTICE_PATH"), ldgroup.getObjName());
//   a.set("TYPE", "FIELD_QUAD3D");
//   a.set("LATTICE_PATH", ldgroup.get_name());
  return ds;
}

template<class Vector>
H5::DataSet WriteVectorField(H5::Group g, const string &name, ConstArray3d<Vector> arr, const LatticeDataQuad3d &ld, const H5::Group ldgroup)
{
  arr = arr[ld.Box()];
  H5::DataSet ds = WriteVectorArray3D<Vector>(g, name, arr);
  //h5cpp::Attributes a = ds.attrs();
  writeAttrToDataset<string>(ds, string("TYPE"), string("FIELD_QUAD3D"));
  writeAttrToDataset<string>(ds, string("LATTICE_PATH"), ldgroup.getObjName());
//   a.set("TYPE", "FIELD_QUAD3D");
//   a.set("LATTICE_PATH", ldgroup.getObjName());
  return ds;
}


template<class T, int mydim>
H5::DataSet WriteAveragedFaceVariables(H5::Group file, const string &id, const ConstArray3d<T> *face_fields)
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
H5::DataSet WriteAveragedFaceVariables(H5::Group file, const string &id, int dim, const ConstArray3d<T> *face_fields)
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
H5::DataSet WriteAveragedFaceVariableField(H5::Group file, const string &id, int dim, const ConstArray3d<T> *face_fields, const LatticeDataQuad3d &ld, const H5::Group ldgroup)
{
  LatticeIndexRanges ir = LatticeIndexRanges::FromCellRange(ld.Box(), dim);
  ConstArray3d<T> my_fields[3];
  for (int j=0; j<dim; ++j)
    my_fields[j] = face_fields[j][ir.faces[j]];
  H5::DataSet ds = WriteAveragedFaceVariables(file, id, dim, my_fields);
  writeAttrToDataset<string>(ds,string("TYPE"), string("FIELD_QUAD3D"));
  writeAttrToDataset<string>(ds,string("LATTICE_PATH"), ldgroup.getObjName());
//   H5::Attributes a = ds.attrs();
//   a.set("TYPE", "FIELD_QUAD3D");
//   a.set("LATTICE_PATH", ldgroup.get_name());
  return ds;
}

template<class T>
string getH5Name(T g)
{
  size_t len = H5Iget_name(g.getId(),NULL,0);
  char buffer[len];
  H5Iget_name(g.getId(),buffer,len+1);
  std::string n = buffer;
  return n;
}

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
//       error.printError();
//     }
//     // catch failure caused by the DataSet operations
//     catch( H5::DataSetIException error )
//     {
//       error.printError();
//     }
//     // catch failure caused by the DataSpace operations
//     catch( H5::DataSpaceIException error )
//     {
//       error.printError();
//     }
//     
//   }//end if
// }


float readAttrFromGroup(H5::Group g, string attr_name)
{
  {
    H5::Attribute myatt_out = g.openAttribute(attr_name);
    float test = 0.0;
    H5::DataType type = myatt_out.getDataType();
    myatt_out.read(type, &test);
    return test;
  } 
}

template <class T>
void readDataSetFromGroup(H5::Group g, string ds_name, T *placeToStore)
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
// 
// string readAttrFromGroup(H5::Group g, string attr_name)
// {
//   //if (typeid(T) == typeid(string()))
//   {
//     H5::Attribute myatt_out = g.openAttribute(attr_name);
//     H5::StrType strdatatype(H5::PredType::C_S1, 256); // of length 256 characters
//     H5std_string strreadbuff("");
//     myatt_out.read(strdatatype,strreadbuff);
//     return string(strreadbuff);
//   } 
// }

template <class T>
void writeAttrToGroup(H5::Group g, string attr_name, T value);

void writeAttrToGroup(H5::Group g, string attr_name, string value)
{ 
  {
    // Create new dataspace for attribute
    H5::DataSpace attr_dataspace = H5::DataSpace(H5S_SCALAR);
    // Create new string datatype for attribute
    H5::StrType strdatatype(H5::PredType::C_S1, 256); // of length 256 characters
    // Set up write buffer for attribute
    //const H5::H5std_string strwritebuf (value);
    //const string strwritebuf (value);
    H5::Attribute myatt_in = g.createAttribute(attr_name, strdatatype, attr_dataspace);
    myatt_in.write(strdatatype, &value);
  }
};

template<class T>
void writeAttrToDataset(H5::DataSet g, string attr_name, T value)
{
  if (typeid(value) == typeid(string()))
  {
    // Create new dataspace for attribute
    H5::DataSpace attr_dataspace = H5::DataSpace(H5S_SCALAR);
    // Create new string datatype for attribute
    H5::StrType strdatatype(H5::PredType::C_S1, 256); // of length 256 characters
    // Set up write buffer for attribute
    //const H5::H5std_string strwritebuf (value);
    const std::string strwritebuf (value);
    H5::Attribute myatt_in = g.createAttribute(attr_name, strdatatype, attr_dataspace);
    myatt_in.write(strdatatype, strwritebuf);
  }
  if (typeid(value) == typeid(int()))
  {
    // Create new dataspace for attribute
    H5::DataSpace attr_dataspace = H5::DataSpace(H5S_SCALAR);
    // Create new string datatype for attribute
  }
}

template<class T>
H5::DataSet writeDataSetToGroup(H5::Group g, string attr_name, T value)
{
  if( typeid(value) == typeid(Float3()))
  {
    hsize_t dims[2] = {1,3};
    int rank = 2;
    myAssert(dims[1] == 3); // need to be a Float3
      
    /* 
      * allocate memory space to read dataset
      */
    H5::DataSpace mspace( rank, dims);
  }
  H5::DataSet ds = H5::DataSet();
  return ds;
}

// #define INSTANTIATE2(T)\
//   template H5::DataSet WriteDataSetToGroup<T>(H5::Group g, string attr_name, T value);
// 
// INSTANTIATE2(DynArray<double>)
  
#define INSTANTIATE(T)\
  template H5::DataSet WriteScalarField<T>(H5::Group g, const string &name, ConstArray3d<T> arr, const LatticeDataQuad3d &ld, const H5::Group ldgroup);

INSTANTIATE(float)
INSTANTIATE(double)
INSTANTIATE(int)
INSTANTIATE(char)

// #define INSTANTIATE_H5Cpp(T)\
//   template<T> readAttrFromGroup<T>(H5::Group g, string &name);
//template float readAttrFromGroup<float>(H5::Group g, string &name);
//template string readAttrFromGroup<string>(H5::Group g, string &name);

// INSTANTIATE_H5Cpp(float)
// INSTANTIATE_H5Cpp(string)

#define INSTANTIATE_VEC(T)\
  template H5::DataSet WriteAveragedFaceVariableField<T>(H5::Group file, const string &id, int dim, const ConstArray3d<T> *face_fields, const LatticeDataQuad3d &ld, const H5::Group ldgroup);\
  template H5::DataSet WriteVectorField<Vec<T,3> >(H5::Group g, const string &name, ConstArray3d<Vec<T,3> > arr, const LatticeDataQuad3d &ld, const H5::Group ldgroup);\
  template H5::DataSet WriteVectorField<Vec<T,2> >(H5::Group g, const string &name, ConstArray3d<Vec<T,2> > arr, const LatticeDataQuad3d &ld, const H5::Group ldgroup);\
  template H5::DataSet WriteVectorField<Vec<T,1> >(H5::Group g, const string &name, ConstArray3d<Vec<T,1> > arr, const LatticeDataQuad3d &ld, const H5::Group ldgroup);\
  template H5::DataSet WriteAveragedFaceVariables<T>(H5::Group file, const string &id, int dim, const ConstArray3d<T> *face_fields);

INSTANTIATE_VEC(float)
INSTANTIATE_VEC(double)
