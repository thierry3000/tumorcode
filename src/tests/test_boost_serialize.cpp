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

/* Our class managing the network is inheret from 
 * boost::noncopyable.
 * The pagmo extension relys heavily on
 * boost::serialization. So lets see if we can serialize
 * a noncopyable object???
 * 
 * example from https://theboostcpplibraries.com/boost.serialization-pointers-and-references
 * https://theboostcpplibraries.com/boost.serialization-pointers-and-references
 * 
 * we see that we can serialize a boost shared_ptr and classes inheret from boost::noncopyable
 */
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <iostream>
#include <sstream>

using namespace boost::archive;
std::stringstream ss;

class animal : boost::noncopyable
{
public:
  animal() = default;
  animal(int legs) : legs_{legs} {}
  int legs() const { return legs_; }

private:
  friend class boost::serialization::access;
  template <typename Archive>
  void serialize(Archive &ar, const unsigned int version) { ar & legs_; }
  int legs_;
};
class bird : public animal
{
public:
  bird() = default;
  bird(int legs, bool can_fly) :
    animal{legs}, can_fly_{can_fly} {}
  bool can_fly() const { return can_fly_; }

private:
  friend class boost::serialization::access;
  template <typename Archive>
  void serialize(Archive &ar, const unsigned int version)
  {
    ar & boost::serialization::base_object<animal>(*this);
    ar & can_fly_;
  }
  bool can_fly_;
};

void save()
{
  text_oarchive oa{ss};
  boost::shared_ptr<animal> a{new animal{4}};
  oa << a;
  boost::shared_ptr<bird> b{new bird{2, false}};
  oa << b;
}

void load()
{
  text_iarchive ia{ss};
  boost::shared_ptr<animal> a;
  ia >>a;
  std::cout << a->legs() << '\n';
  boost::shared_ptr<bird> penguin;
  ia >>penguin;
  std::cout << penguin->legs() << " and " << std::boolalpha << penguin->can_fly() << '\n';
}

int main()
{
  save();
  load();
}