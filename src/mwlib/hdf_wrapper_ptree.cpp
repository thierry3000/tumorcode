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

#include <boost/optional.hpp>

#include "hdf_wrapper_ptree.h"
#include "hdf_wrapper_vec.h"
#include "../common/hdfio.h"
#include "dynamicarray.h"
#include <iostream>

using boost::property_tree::ptree;

void WriteHdfPtree(H5::Group &f, const ptree &pt, HdfWritePtreeAs storage_mode)
{
  typedef ptree::const_iterator I;
  for (I it = pt.begin(); it != pt.end(); ++it)
  {
    // get all values as string. There is no other way, because there is no type information stored in ptree
    const string k = it->first;
    const ptree &v = it->second;
    if (!v.data().empty())
    {
      if (storage_mode == HDF_WRITE_PTREE_AS_ATTRIBUTE)
      {
// 	std::cout << f.getFileName() << std::endl;
// 	std::cout << "print k: "<< k << "v: " << v.data() << std::endl; 
	writeAttrToH5(f, k, v.data());
      }
      else
      {
	DynArray<string> this_data;
	for( auto entry: v.data())
	{
	  this_data.push_back(&entry);
	}
	writeDataSetToGroup(f,k, this_data);
      }
//         f.attrs().set(k, v.data());
/*	
      else
        h5cpp::create_dataset_scalar<string>(f, k, v.data());*/
    }
    // recurse
    if (v.begin() != v.end())
    {
      H5::Group h5_group = f.createGroup(k);
      WriteHdfPtree(h5_group, v, storage_mode);
    }
  }
}
void printObjectName(const H5::Group &g)
{
  try
  {
    std::cout << g.getObjName() << std::endl;
  }
  catch(H5::Exception e)
  {
    e.printError();
  }
}

void printObjectName( H5::H5Location &g, const H5std_string attr_name, void* operator_data)
{
  try
  {
    std::cout << attr_name << std::endl;
  }
  catch(H5::Exception e)
  {
    e.printError();
  }
}

void add_all_attributes_to_ptree( H5::H5Location &g, const H5std_string attr_name, void* p_to_user_data)
{
  boost::property_tree::ptree *pt = static_cast<boost::property_tree::ptree *>(p_to_user_data);
  try
  {
    std::string output_buffer; // I hope all ptree values are also strings
    readAttrFromH5(g, attr_name, output_buffer);
#ifndef NDEBUG
    std::cout << "adding " << attr_name << " : " << output_buffer << " to ptree." << std::endl;
#endif
    pt->put(attr_name, output_buffer);
  }
  catch(H5::Exception e)
  {
    e.printError();
  }
}

void ReadHdfPtree(boost::property_tree::ptree &pt, H5::Group &f, HdfWritePtreeAs storage_mode)
{
#ifndef NDEBUG
  try
  {
    std::cout << "read: " << f.getFileName() << std::endl;
    std::cout << "at:   " << f.getObjName() << std::endl;
  }
  catch(H5::Exception e)
  {
    e.printError();
  }
#endif

  //f.iterateAttrs(&printObjectName, NULL, NULL);
  f.iterateAttrs(&add_all_attributes_to_ptree, NULL, &pt);
  //h5_vessel_parameters = h5_parameters.createGroup("vessels");
  //f.itera
#ifndef NDEBUG
  std::cout << "content of ptree: " << std::endl;
  for( ptree::const_iterator it = pt.begin(); it !=pt.end(); ++it)
  {
    // get all values as string. There is no other way, because there is no type information stored in ptree
    const string k = it->first;
    const ptree &v = it->second;
    std::cout << "key: content " << k << " : " << v.data() << std::endl; 
  }
#endif
  /*
  typedef ptree::const_iterator I;
  for (I it = pt.begin(); it != pt.end(); ++it)
  {
    // get all values as string. There is no other way, because there is no type information stored in ptree
    const string k = it->first;
    const ptree &v = it->second;
    if (!v.data().empty())
    {
      if (storage_mode == HDF_WRITE_PTREE_AS_ATTRIBUTE)
      {
// 	std::cout << f.getFileName() << std::endl;
// 	std::cout << "print k: "<< k << "v: " << v.data() << std::endl; 
	writeAttrToH5(f, k, v.data());
      }
      else
      {
	DynArray<string> this_data;
	for( auto entry: v.data())
	{
	  this_data.push_back(&entry);
	}
	writeDataSetToGroup(f,k, this_data);
      }
//         f.attrs().set(k, v.data());

//      else
//        h5cpp::create_dataset_scalar<string>(f, k, v.data());
    }
    // recurse
    if (v.begin() != v.end())
    {
      H5::Group h5_group = f.createGroup(k);
      WriteHdfPtree(h5_group, v, storage_mode);
    }
  }
  */
}
