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
#include "any_tree.h"
#include "ptree_ext.h"
#include <boost/format.hpp>

namespace boost { namespace property_tree {


static void asptree_(ptree &dst, const atree &src)
{
  boost::any q = src.data();
  if (q.empty())
  { /* do nothing */ }
  else if (q.type() == typeid(int))
    dst.put_value(boost::any_cast<int>(q));
  else if (q.type() == typeid(float))
    dst.put_value(boost::any_cast<float>(q));
  else if (q.type() == typeid(double))
    dst.put_value(boost::any_cast<double>(q));
  else if (q.type() == typeid(std::string))
    dst.put_value(boost::any_cast<std::string>(q));
  else
    throw std::invalid_argument(boost::str(boost::format("asptree: extraction of type %s not implemented") % q.type().name()));

  typedef atree::const_iterator I;
  for (I it = src.begin(); it != src.end(); ++it)
  {
    const std::string name = it->first;
    // swapped string() woth ptree()
    ptree &dstchild = dst.push_back(std::make_pair(name, ptree()))->second;
    asptree_(dstchild, it->second);
  }
}

ptree asptree(const atree &src)
{
  ptree pt;
  asptree_(pt, src);
  return pt;
}

}}