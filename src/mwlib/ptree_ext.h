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
#ifndef PTREE_EXT_H
#define PTREE_EXT_H

#include <boost/property_tree/ptree.hpp>
#include <boost/any.hpp>
#include <boost/optional.hpp>

class T;
namespace boost { namespace property_tree {

ptree parse_program_options(int argc, char* argv[], ptree implicit_opts);


// Merge src into dst, overwriting existing items
// note:
// -- branches which have no name are copied (even if the content of the branch already exists in an unnamed branch of the destination. This function cannot now that!)
// -- branches which do not exsist are created
// -- if the branch names are not unique, one of the non-unique branches is picked to be updated. Which one is not specified. The others are unchanged.
void update(boost::property_tree::ptree &dst, const boost::property_tree::ptree &src);

ptree subtract(const ptree &pta, const ptree &ptb);
ptree remove(const ptree &pta, const ptree &ptb);

// Get a value from a property tree. This doesn't require the template argument.
template<class T>
inline void get(T &val, const std::string &name, const boost::property_tree::ptree &pt)
{
  val = pt.get<T>(name);
}

template<class T>
inline T get(const ptree &pt, const std::string &name, const T &default_)
{
  boost::optional<T> r = pt.get_optional<T>(name);
  if (r) return *r;
  else return default_;
}


class make_ptree_
{
  ptree pt;
public:
  make_ptree_() {}
  make_ptree_(const ptree &other) : pt(other) {}
  
  template<class T>
  make_ptree_& operator()(const std::string &key, const T &value)
  {
    pt.put(key, value);
    return *this;
  }

  template<class T>
  make_ptree_& operator()(const ptree &other)
  {
    update(pt, other);
    return *this;
  }
  
  operator ptree() const { return pt; }
};

template<class T>
inline make_ptree_ make_ptree(const std::string &key, const T &value)
{
  return make_ptree_()(key, value);
}

inline make_ptree_ make_ptree(const ptree &other)
{
  return make_ptree_(other);
}


template<class Datatype>
inline basic_ptree<std::string, Datatype>& require_child(basic_ptree<std::string, Datatype> &pt, const std::string &name)
{
  typedef basic_ptree<std::string, Datatype> TreeType;
  boost::optional<TreeType &> c = pt.get_child_optional(name);
  if (c)
    return *c;
  else
    return pt.add_child(name, TreeType());
}



#if 0
#include <boost/foreach.hpp>
#include <vector>
// this works!
//template<class Key, class Value>
class basic_ptree_walker
{
  //typedef basic_ptree<Key, Value> treetype;
  typedef ptree treetype;
  typedef treetype::iterator iterator;
  std::vector<std::pair<iterator, iterator> > stack;
  
public:
  basic_ptree_walker(ptree &pt)
  {
    stack.push_back(std::make_pair(pt.begin(), pt.end()));
  }
  void increment()
  {
    if (stack.empty()) return; // done
    while (true)
    {
      iterator it  =  stack.back().first;
      if (it == stack.back().second)
      {
        while (it == stack.back().second)
        {
          stack.pop_back();
          if (stack.empty()) return;
          ++stack.back().first;
          it = stack.back().first;
        }
        return;
      }
      stack.push_back(std::make_pair(it->second.begin(), it->second.end()));
      if (!it->second.empty()) return;
    }
  }
  treetype::value_type& dereference() const { return *stack.back().first; }

  bool is_valid() const { return !stack.empty(); }
};
#endif

}}


#endif
