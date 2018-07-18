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
#ifndef HDF_WRAPPER_PTREE_H
#define HDF_WRAPPER_PTREE_H

#include <boost/property_tree/ptree.hpp>

#include <H5Cpp.h>

enum HdfWritePtreeAs {
  HDF_WRITE_PTREE_AS_ATTRIBUTE,
  HDF_WRITE_PTREE_AS_DATASETS
};

void WriteHdfPtree(H5::Group &f, const boost::property_tree::ptree &pt, HdfWritePtreeAs storage_mode = HDF_WRITE_PTREE_AS_ATTRIBUTE);
void ReadHdfPtree(boost::property_tree::ptree &pt, H5::Group &f, HdfWritePtreeAs storage_mode = HDF_WRITE_PTREE_AS_ATTRIBUTE);


#endif
