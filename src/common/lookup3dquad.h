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
#ifndef _LOOKUP3D_QUAD_H
#define _LOOKUP3D_QUAD_H

#include "mwlib/lattice-data.h"
#include "lattice-data-polymorphic.h"

#include<boost/unordered_map.hpp>

#define DEBUG_BondLookup(x)

//template<class LD, class V_>
class BondLookup
{
  typedef polymorphic_latticedata::LatticeData LD;
  typedef LD::SiteType S;
  typedef LD::LatticeIndexType I;
  typedef Vessel* V;
//   typedef typename LD::SiteType S;  // site type
//   typedef typename LD::LatticeIndexType I;
  //typedef V_ V; // value type
  typedef std::pair<S,S> K;
  typedef boost::unordered_map<K, V> Map;
  const LD *m_ld;
  Map map;
  V default_value;

  inline K make_key(S a, S b) const
  {
    if (a < b) std::swap(a,b);
    return K(a,b);
  }
  
public:
  BondLookup() {}

  void Init(const LD& _ld, V default_value_ = V())
  {
    map.clear();
    m_ld = &_ld;
    default_value = default_value_;
  }

  const LD& ld() const { return *m_ld; }

  void Clear() 
  {
    m_ld = NULL;
    map.clear();
  }

  void Insert(S site, int dir, int len, V x);

  V Remove(S site, int dir, int len)
  {
    V ret = default_value;
    S nbsite;
DEBUG_BondLookup(
    bool found = false;
)
    for( int j=0; j<len; ++j)
    {
      nbsite = ld().NbSite(site, dir);
      Map::iterator it = map.find(make_key(site, nbsite));
      if (it == map.end()) continue;
      V val = it->second;
      map.quick_erase(it);
DEBUG_BondLookup(
      myAssert(found == false || val == ret);
      found = true;
)     
      ret = val;
      site = nbsite;
    }
    return ret;
  }

  V Find(S site, int dir, int len)
  {
    V ret = default_value;
    S nbsite;
DEBUG_BondLookup(
    bool found = false;
)
    for( int j=0; j<len; ++j)
    {
      nbsite = ld().NbSite(site, dir);
      Map::iterator it = map.find(make_key(site, nbsite));
      if (it == map.end()) continue;
      V val = it->second;
DEBUG_BondLookup(
      myAssert(found == false || val == ret);
      found = true;
)
      ret = val;
      site = nbsite;
    }
    return ret;
  }

  void PrintContents(std::ostream &os) const;

//   std::size_t estimateMemoryUsage() const
//   {
//     return sizeof(*this)-sizeof(Map)+estimateMemoryUsage(map);
//   }
};


//template<class LD, class V>
class SiteLookup
{
  typedef polymorphic_latticedata::LatticeData LD;
  typedef LD::SiteType S;
  typedef VesselNode* V;
  
  typedef S K; // key type
  typedef boost::unordered_map<K, V> Map;
  Map map;
  const LD* m_ld;
  //boost::shared_ptr<LD> m_ld;
  V default_value;

public:
  SiteLookup() : default_value(),m_ld(NULL) {}
  const LD& ld() const { return *m_ld; }
  //const LD& ld() const { return *m_ld.get(); }
  S Count() const { return map.size(); }
  void Init( const LD &_ld, V default_value_ = V())
  {
    map.clear();
    m_ld = &_ld;
    //boost::shared_ptr<LD> p(new (_ld.Clone()));
    //m_ld = _ld.Clone();
    default_value = default_value_;
  }

  void Clear()
  {
    map.clear();
    m_ld = NULL;
  }

  void InsertAtSite(S s, V x )
  {
    map[s] = x;
  }

  V RemoveAtSite(S s )
  {
    Map::iterator i = map.find(s);
    if(i == map.end()) return default_value;
    V q = i->second;
    map.quick_erase(i);
    return q;
  }
  V FindAtSite(S s )
  {
    Map::iterator i = map.find(s);
    if (i == map.end()) return default_value;
    else return i->second;
  }

  void PrintContents(std::ostream &os) const;

//   std::size_t estimateMemoryUsage() const
//   {
//     return sizeof(*this)-sizeof(Map)+estimateMemoryUsage(map);
//   }
};


#endif
