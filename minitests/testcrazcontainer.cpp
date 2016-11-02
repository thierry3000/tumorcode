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
#include <boost/any.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/unordered_map.hpp>
#include <boost/property_tree/detail/ptree_utils.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/fusion/algorithm/iteration/for_each.hpp>
#include <boost/fusion/tuple.hpp>
//#include <boost/tuple/tuple.hpp>
#include <iostream>
 

class Container
{
  typedef boost::unordered_map<std::string, boost::any> unordered_any_map;
  unordered_any_map map;
public:

  template<class T>
  void insert_make(const std::string &key)
  {
    map[key] = boost::any(boost::make_shared<T>());
  }
  
  template<class T, class A1>
  void insert_make(const std::string &key, A1 const & a1)
  {
    map[key] = boost::any(boost::make_shared<T, A1>(a1));
  }
  
  template<class T, class A1, class A2>
  void insert_make(const std::string &key, A1 const & a1, A2 const &a2)
  {
    map[key] = boost::any(boost::make_shared<T, A1, A2>(a1, a2));
  }
  
  template<class T>
  T& get(const std::string &key) 
  {
    typedef boost::shared_ptr<T> ptr_t;
    ptr_t p = boost::any_cast<ptr_t>(map[key]);
    return *p;
  }
  
  template<class T>
  const T& get(const std::string &key) const
  {
    return const_cast<Container*>(this)->get<T>(key);
  }
  
  void add_shared_to(Container &other, const std::string &key)
  {
    other.map[key] = map[key];
  }
  
  void clear()
  {
    map.clear();
  }
};

#define DEFINE_FUNCTOR1_START(name, t1, a1) \
  struct name { \
    t1 a1; \
    name(t1 a1) : a1(a1) {}

#define DEFINE_FUNCTOR1_END };


#define DEFINE_FUNCTOR4_START(name, t1, a1, t2, a2, t3, a3, t4, a4) \
  struct name { \
    t1 a1; t2 a2; t3 a3; t4 a4; \
    name(t1 a1, t2 a2, t3 a3, t4 a4) : a1(a1), a2(a2), a3(a3), a4(a4) {}

#define DEFINE_FUNCTOR4_END };

DEFINE_FUNCTOR1_START(Fubar, double, fu)
DEFINE_FUNCTOR1_END

DEFINE_FUNCTOR4_START(AddScaled, double, fu, double, fv, Container &, u, const Container &, v)
  template<typename T>
  void operator()(T &key) const
  {
    // t is a pair as above
    u.get<double>(key) *= fu;
    u.get<double>(key) += fv * v.get<double>(key);
  }
DEFINE_FUNCTOR4_END



void test()
{
  boost::fusion::tuple<const char*, const char*> keys("k1","k2");
  Container c;
  c.insert_make<double>("k1",5);
  c.insert_make<double>("k2",5.);
  
  Container c2;
  c2.insert_make<double>("k1",6);
  c2.insert_make<double>("k2",6.);
  
  boost::fusion::for_each(keys, AddScaled(100., 1., c, c2));
  std::cout << c.get<double>("k1") << std::endl;
  std::cout << c.get<double>("k2") << std::endl;
}


 
int main(int argc, char **argv)
{
  test();
//   Container c;
//   c.insert_make<int>("k1",5);
//   c.insert_make<std::pair<double,int> >("k2",5., 1.);
//   std::cout << c.get<int>("k1") << std::endl;
//   std::pair<double,int> x1 = c.get<std::pair<double, int> >("k2");
//   std::cout << x1.first << " " << x1.second << std::endl;
//   //c.get<double>("k1"); // fails
//   
//   
//   
//   Container c2;
//   c.add_shared_to(c2,"k1");
//   c.get<int>("k1") = 9001;
//   c.clear();
//   std::cout << c2.get<int>("k1") << std::endl;
}