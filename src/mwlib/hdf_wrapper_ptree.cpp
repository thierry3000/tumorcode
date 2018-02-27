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
	std::cout << "print k: "<< k << "v: " << v.data() << std::endl; 
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
