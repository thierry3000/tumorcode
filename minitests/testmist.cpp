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

#include <iostream>
#include <vector>
#include <auto_ptr.h>
#include <boost/scoped_array.hpp>

using std::cout;
using std::endl;

namespace Aspace
{

struct A
{
  int x;
};

void f(A &a) { cout << a.x << endl; }

}


namespace Bspace
{

struct B : public Aspace::A
{
  int y;
  void foo() { f((Aspace::A&)*this); }
};

void f(B &b) { cout << b.y << endl; }

}

template <class Y>
struct Args {
  typedef Y value_type;
};

template <class A>
struct Base : public A
{
};

template <class A>
struct Derived : public Base<A>
{
  void fun(typename A::value_type &x)
  {
    x = 0;
  }
};



#if 0
void test2()
{
  Smallest<double> a(2.);
  //a = 5;
  //a = MinimumDouble(10.);
  a |= Smallest<double>(0.1);
  a |= 3.;
  double q = a;
  std::cout << q << std::endl;
}
#endif

void test()
{
  Aspace::A a;
  Bspace::B b;
  f(a);
  f(b);
  b.foo();
  Derived<Args<double> > d;
  std::auto_ptr<int> p(new int(3));
  boost::scoped_array<std::auto_ptr<int> > v(new std::auto_ptr<int>[2]);
  v[0] = p;
  cout << *v[0];
};

int main(int argc, char **argv)
{
  test();
  return 0;
}


