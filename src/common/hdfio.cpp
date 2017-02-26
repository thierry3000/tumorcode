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


h5cpp::Group RequireLatticeDataGroup(h5cpp::Group g, const string &name, const LatticeDataQuad3d &ld)
{
  bool had_ld_group = false;
  h5cpp::Group ld_group = g.require_group(name, &had_ld_group);
  if (!had_ld_group)
    WriteHdfLd(ld_group, ld);
  return ld_group;
}




void WriteHdfGraph( h5cpp::Group g, const VesselList3d &vl )
{
#ifdef DEBUG
  printf("Starting to write hdf\n");
#endif
  const int ncnt = vl.GetNCount();
  const int ecnt = vl.GetECount();
  myAssert(ncnt>0 && ecnt>0);
  
  h5cpp::Group gg = g.create_group("nodes");
  gg.attrs().set("COUNT",ncnt);
  h5cpp::Attributes attrs = g.attrs();
  if(vl.HasLattice())//LATTICE IS THERE
  {
    if(attrs.exists("CLASS"))//correct hdf attribute is checke
    {
      myAssert(attrs.get<std::string>("CLASS") == "GRAPH");
    }
    else
    {
      attrs.create<std::string>("CLASS", "GRAPH");
    }
    //lattice stuff is writen to hdf
    vl.Ld().WriteHdf(g.create_group("lattice"));
    const VesselList3d::LatticeData &ld = vl.Ld();
    DynArray<int> a(ncnt);
    for(int i=0; i<ncnt; ++i) 
    {
      a[i] = ld.LatticeToSite(vl.GetNode(i)->lpos);
    }
    h5cpp::Dataset ds = h5cpp::create_dataset<int>(gg,  "lattice_pos", a);
    ds.attrs().set("MODE","linear");
  }
  else//no lattice pressent
  {
    if(attrs.exists("CLASS"))
    {
      myAssert(attrs.get<std::string>("CLASS") == "REALWORLD");
    }
    else
    {
      attrs.create<std::string>("CLASS", "REALWORLD");
    }
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
    if(!gg.exists("world_pos"))
    {
      h5cpp::Dataset ds = h5cpp::create_dataset<float>(gg, "world_pos", h5cpp::Dataspace::simple_dims(a.size()/3,3), &a[0]);
    }
  }//end no lattice pressent
  
  {//interrupt begin
  //write roots
  DynArray<int> roots(16,ConsTags::RESERVE);
  for(int i=0; i<ncnt; ++i) 
  {
    if(vl.GetNode(i)->IsBoundary())
      roots.push_back(i);
  }
  if(!gg.exists("roots"))
  {
    h5cpp::create_dataset<int>(gg,  "roots", roots);
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
	bctyp_index.push_back((int)it->second.type);
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
    h5cpp::create_dataset<int>(gg,  "bc_node_index", bc_node_index);
    h5cpp::create_dataset<int>(gg,  "bc_type", bctyp_index);
    h5cpp::create_dataset<float>(gg,  "bc_value", values_of_bcs);
    h5cpp::create_dataset<float>(gg,  "bc_conductivity_value", bc_conductivity_value);
  }
  }//for interrupt 

  gg = g.create_group("edges");
  gg.attrs().set("COUNT",ecnt);
  
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
    h5cpp::create_dataset<int>( gg, "node_a_index", va);
    h5cpp::create_dataset<int>( gg, "node_b_index", vb );
    h5cpp::create_dataset<int>( gg, "flags", flags );
    h5cpp::Dataset ds = h5cpp::create_dataset<float>( gg, "radius", float_radius );
    ds.attrs().set("MODE","const");
  }
}


void ReadHdfGraph( h5cpp::Group g, VesselList3d *vl )
{
  h5cpp::Group gnodes = g.open_group("nodes");
  h5cpp::Group gedges = g.open_group("edges");
  
  int ecnt=0,ncnt=0;
  gnodes.attrs().get("COUNT",ncnt);
  gedges.attrs().get("COUNT",ecnt);
  
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
}

template<class T>
h5cpp::Dataset WriteScalarField(h5cpp::Group g, const string &name, ConstArray3d<T> arr, const LatticeDataQuad3d &ld, const h5cpp::Group ldgroup, const h5cpp::Datatype &disktype)
{
  arr = arr[ld.Box()];
  h5cpp::Dataset ds = WriteArray3D<T>(g, name, arr, disktype);
  h5cpp::Attributes a = ds.attrs();
  a.set("TYPE", "FIELD_QUAD3D");
  a.set("LATTICE_PATH", ldgroup.get_name());
  return ds;
}

template<class Vector>
h5cpp::Dataset WriteVectorField(h5cpp::Group g, const string &name, ConstArray3d<Vector> arr, const LatticeDataQuad3d &ld, const h5cpp::Group ldgroup)
{
  arr = arr[ld.Box()];
  h5cpp::Dataset ds = WriteVectorArray3D<Vector>(g, name, arr);
  h5cpp::Attributes a = ds.attrs();
  a.set("TYPE", "FIELD_QUAD3D");
  a.set("LATTICE_PATH", ldgroup.get_name());
  return ds;
}


template<class T, int mydim>
h5cpp::Dataset WriteAveragedFaceVariables(h5cpp::Group file, const std::string &id, const ConstArray3d<T> *face_fields)
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
h5cpp::Dataset WriteAveragedFaceVariables(h5cpp::Group file, const std::string &id, int dim, const ConstArray3d<T> *face_fields)
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
h5cpp::Dataset WriteAveragedFaceVariableField(h5cpp::Group file, const std::string &id, int dim, const ConstArray3d<T> *face_fields, const LatticeDataQuad3d &ld, const h5cpp::Group ldgroup)
{
  LatticeIndexRanges ir = LatticeIndexRanges::FromCellRange(ld.Box(), dim);
  ConstArray3d<T> my_fields[3];
  for (int j=0; j<dim; ++j)
    my_fields[j] = face_fields[j][ir.faces[j]];
  h5cpp::Dataset ds = WriteAveragedFaceVariables(file, id, dim, my_fields);
  h5cpp::Attributes a = ds.attrs();
  a.set("TYPE", "FIELD_QUAD3D");
  a.set("LATTICE_PATH", ldgroup.get_name());
  return ds;
}


#define INSTANTIATE(T)\
  template h5cpp::Dataset WriteScalarField<T>(h5cpp::Group g, const string &name, ConstArray3d<T> arr, const LatticeDataQuad3d &ld, const h5cpp::Group ldgroup, const h5cpp::Datatype &disktype);

INSTANTIATE(float)
INSTANTIATE(double)
INSTANTIATE(int)
INSTANTIATE(char)

#define INSTANTIATE_VEC(T)\
  template h5cpp::Dataset WriteAveragedFaceVariableField<T>(h5cpp::Group file, const std::string &id, int dim, const ConstArray3d<T> *face_fields, const LatticeDataQuad3d &ld, const h5cpp::Group ldgroup);\
  template h5cpp::Dataset WriteVectorField<Vec<T,3> >(h5cpp::Group g, const string &name, ConstArray3d<Vec<T,3> > arr, const LatticeDataQuad3d &ld, const h5cpp::Group ldgroup);\
  template h5cpp::Dataset WriteVectorField<Vec<T,2> >(h5cpp::Group g, const string &name, ConstArray3d<Vec<T,2> > arr, const LatticeDataQuad3d &ld, const h5cpp::Group ldgroup);\
  template h5cpp::Dataset WriteVectorField<Vec<T,1> >(h5cpp::Group g, const string &name, ConstArray3d<Vec<T,1> > arr, const LatticeDataQuad3d &ld, const h5cpp::Group ldgroup);\
  template h5cpp::Dataset WriteAveragedFaceVariables<T>(h5cpp::Group file, const std::string &id, int dim, const ConstArray3d<T> *face_fields);

INSTANTIATE_VEC(float)
INSTANTIATE_VEC(double)
