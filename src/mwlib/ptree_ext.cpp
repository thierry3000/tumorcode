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
#include <string>
#include <boost/property_tree/ptree.hpp>

using std::string;

namespace boost { namespace property_tree {

void update(ptree &dst, const ptree &src)
{
  typedef ptree::const_iterator I;
  for (I it = src.begin(); it != src.end(); ++it)
  {
    const string name = it->first;
    if (name.empty() || dst.find(name) == dst.not_found())
    {
      // if the source branch name is empty, or the destination tree does not have this branch,
      // then the whole branch is copied over to the destination
      dst.push_back(*it);
    }
    else
    {
      // in this case the destination tree has a matching branch
      // thus the data is copied and then this function is called recursively for the branch
      ptree &c = dst.to_iterator(dst.find(name))->second;
      c.data() = it->second.data();
      update(c, it->second);  
    }  
  }
}


static void subtract_(ptree &pta, const ptree &ptb)
{
  typedef ptree::const_iterator I;
  for (I it = ptb.begin(); it != ptb.end(); ++it) // go through b and look for matching nodes in a
  {
    const string name = it->first;
    if (name.empty() || pta.find(name) == pta.not_found()) continue; // pta does not have it, or the node name in ptb is empty
    ptree &c = pta.to_iterator(pta.find(name))->second; // find the node in pta
    subtract_(c, it->second); // recurse
    if (c.begin() == c.end()) // if the recursion left this node with no childs, then it is okay to delete it as well
      pta.erase(name);
  }
}

ptree subtract(const ptree &pta, const ptree &ptb)
{
  ptree res(pta);
  subtract_(res, ptb);
  return res;
}



static void remove_(ptree &pta, const ptree &ptb)
{
  typedef ptree::const_iterator I;
  for (I it = ptb.begin(); it != ptb.end(); ++it) // go through b and look for matching nodes in a
  {
    const string name = it->first;
    if (name.empty() || pta.find(name) == pta.not_found()) continue; // pta does not have it, or the node name in ptb is empty
    ptree &c = pta.to_iterator(pta.find(name))->second; // find the node in pta
    if (it->second.begin() == it->second.end()) // remove node from a if b has no children here
    {
      pta.erase(name);
    }
    else
    {
      remove_(c, it->second); // recurse
    }
  }
}

ptree remove(const ptree &pta, const ptree &ptb)
{
  ptree res(pta);
  remove_(res, ptb);
  return res;
}

}}