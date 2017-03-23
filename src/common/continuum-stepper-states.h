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
#ifndef _CONTINUUM_STEPPER_STATES_H_
#define _CONTINUUM_STEPPER_STATES_H_

#include "mwlib/myAssert.h"

#include <boost/any.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/unordered_map.hpp>

#include <boost/function.hpp>
#include <boost/foreach.hpp>

/*
 * The idea is to increment a counter every time a state is changed via
 * the time stepper interface. The counter value could be used in a lookup table
 * for additional data.
 */
class IncrementingIntOps
{
public:
  typedef int state_type;
  IncrementingIntOps() {}

  void addScaled(double fa, state_type &u, double fb, const state_type &v)
  {
    if (fa == 1. && fb == 0.) return;
    ++u;
  }

  void initFrom(state_type &u, const state_type &other, ConsMode mode)
  {
    if (mode == ConsMode::AS_COPY)
      u = other;
    else
      u = 0; // clean
  }
};


/*------------------------------------------------------
------------------------------------------------------*/

#if 0
class Container
{
  typedef boost::unordered_map<std::string, boost::any> unordered_any_map;
  unordered_any_map map;
public:

  template<class T>
  T& make_insert(const std::string &key)
  {
    boost::shared_ptr<T> p(boost::make_shared<T>());
    map[key] = boost::any(p);
    return *p;
  }

  template<class T, class A1>
  T& make_insert(const std::string &key, A1 const & a1)
  {
    boost::shared_ptr<T> p(boost::make_shared<T>(a1));
    map[key] = boost::any(p);
    return *p;
  }

  template<class T, class A1, class A2>
  T& make_insert(const std::string &key, A1 const & a1, A2 const &a2)
  {
    boost::shared_ptr<T> p(boost::make_shared<T>(a1, a2));
    map[key] = boost::any(p);
    return *p;
  }

  template<class T>
  T& get(const std::string &key)
  {
    myAssert(map.find(key) != map.end());
    typedef boost::shared_ptr<T> ptr_t;
    ptr_t p = boost::any_cast<ptr_t>(map[key]);
    return *p;
  }

  template<class T>
  const T& get(const std::string &key) const
  {
    return const_cast<Container*>(this)->get<T>(key);
  }

  bool has(const std::string &key) const
  {
    return map.find(key) != map.end();
  }
  
  void insert_shared_into(Container &other, const std::string &key) const
  {
    // the const cast seems ok. The shared element can be manipulated through
    // 'other' but this function does not change 'this'
    other.map[key] = const_cast<Container*>(this)->map[key];
  }

  void clear()
  {
    map.clear();
  }
};







namespace ContainerOpsInternal
{

template<class Ops>
static void addScaled(const std::string &key, const Ops &ops, double fu, Container &u, double fv, const Container &v)
{
  typedef typename Ops::state_type X;
  ops.addScaled(fu, u.get<X>(key), fv, v.get<X>(key));
}

template<class Ops>
static void initFrom(const std::string &key, const Ops &ops, Container &u, const Container &v, ConsMode mode)
{
  typedef typename Ops::state_type X;
  ops.initFrom(u.make_insert<X>(key), v.get<X>(key), mode);
}


static void addScaled_nop(double fu, Container &u, double fv, const Container &v)
{
}

static void initFrom_shared(const std::string &key, Container &u, const Container &v, ConsMode mode)
{
  v.insert_shared_into(u, key);
}


}


template<class Ops>
static OpsFuncs<Container> make_opsfuncs(const std::string &key, const Ops &ops)
{
  OpsFuncs<Container> fcops;
  fcops.addScaled = boost::bind(ContainerOpsInternal::addScaled<Ops>, key, boost::cref(ops), _1, _2, _3, _4);
  fcops.initFrom  = boost::bind(ContainerOpsInternal::initFrom<Ops>, key, boost::cref(ops), _1, _2, _3);
  return fcops;
}

static OpsFuncs<Container> make_opsfuncs_so_initFrom_shares(const std::string &key)
{
  OpsFuncs<Container> fcops;
  fcops.addScaled = ContainerOpsInternal::addScaled_nop;
  fcops.initFrom  = boost::bind(ContainerOpsInternal::initFrom_shared, key, _1, _2, _3);
  return fcops;
}


class ContainerOps
{
  std::vector<OpsFuncs<Container> > data;
  
public:
  typedef Container state_type;
  
  void push_back(const OpsFuncs<Container> &ops)
  {
    data.push_back(ops);
  }
  
  void addScaled(double fu, state_type &u, double fv, const state_type &v)
  {
    BOOST_FOREACH(const OpsFuncs<Container> &o, data)
    {
      o.addScaled(fu, u, fv, v);
    }
  }

  void initFrom(state_type &u, const state_type &v,  ConsMode mode)
  {
    BOOST_FOREACH(const OpsFuncs<Container> &o, data)
    {
      o.initFrom(u, v, mode);
    }
  }
};
#endif

#endif