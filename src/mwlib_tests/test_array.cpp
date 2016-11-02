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
#include <boost/tuple/tuple.hpp>
#include <boost/format.hpp>

#include "mwlib/field.h"

using std::cout;
using std::endl;
using boost::format;

// void bla()
// {
//   std::vector<BBox3> par_boxes = MpBoxes(bb);
//   #pragma omp parallel for
//   for (int i=0; i<par_boxes.size(); ++i) FOR_BBOX3(p, par_boxes[i])
//   {
//   }
// }

#define FOR_BBOX3(p,bb)\
  for(std::pair<const BBox3, bool> __forbb(bb,true); __forbb.second; __forbb.second=!__forbb.second)\
    for(Int3 p=__forbb.first.min; p[2]<=__forbb.first.max[2]; ++p[2]) for(p[1]=__forbb.first.min[1]; p[1]<=__forbb.first.max[1]; ++p[1]) for(p[0]=__forbb.first.min[0]; p[0]<=__forbb.first.max[0]; ++p[0])


      
void TestMoveFuncCall(ConstArray3d<int> arr)
{
  cout << "const arr has ref count " << arr.refCount() << endl;
}

void TestMoveFuncCallNonConst(Array3d<int> arr)
{
  cout << "arr has ref count " << arr.refCount() << endl;
}


void TestArray()
{
  const BBox3 bb(0,0,0, 10, 10, 0);
  const BBox3 bbext = Extend(bb, Int3(2,2,0));
  const BBox3 b1(-1,-1, 0, 2, 2, 0);
  
  Array3d<int> a(bbext);
  a.fill(-1);
  int k=0;
  FOR_BBOX3(p, bb)
  {
    a(p) = k++;
  }
  cout << "initial array:" << endl;
  PrintArray(a, cout);
  a.setBox(bb);
  cout << "center array:" << endl;
  PrintArray(a, cout);

  Array3d<int> a_slice = a[b1];
  cout << format("a slice by %s:") % b1 << endl;
  PrintArray(a_slice, cout);

  a_slice = a[bbext];
  cout << format("a slice by %s:") % bbext << endl;
  PrintArray(a_slice, cout);
  
  Array3d<int> a_moved = a;
  Int3 m(10, 10, 0);
  a_moved.move(m);
  cout << format("a moved by %s:") % m << endl;
  PrintArray(a_moved, cout);


  Array3d<int> a_cloned; a_cloned.initDeepCopy(a_moved);
  cout << "a deep copied" << endl;
  PrintArray(a_cloned, cout);

//   CopyBorder(a, 2, 2);
//   cout << "a with copied border" << endl;
//   PrintArray(a[bbext]);

  Array3d<int> keeper = a;
  cout << "---- move cons ----" << endl;
  cout << "before a refcount is " << a.refCount() << endl;
  Array3d<int> moved(a.release());
  cout << "after a refcount is" << a.refCount() <<  endl;
  cout << "move target refCount is" << moved.refCount() << endl;

  cout << "---- move assign ----" << endl;
  a = keeper;
  cout << "before a refcount is " << a.refCount() << endl;
  moved = a.release();
  cout << "after a refcount is" << a.refCount() <<  endl;
  cout << "move target refCount is" << moved.refCount() << endl;

  cout << "---- non const function call ---" << endl;
  cout << "refcount outside of the function " << moved.refCount() << " and is empty " << moved.empty() << endl;
  TestMoveFuncCallNonConst(moved.release());
  //TestMoveFuncCallNonConst(ConstArray3d<int>(moved).release()); // compile error
  cout << "empty outside of the function " << moved.empty() << endl;

  cout << "---- const function call ---" << endl;
  ConstArray3d<int> const_moved(a_cloned);
  a_cloned.clear();
  cout << "refcount outside of the function " << const_moved.refCount() << " and is empty " << const_moved.empty() << endl;
  TestMoveFuncCall(const_moved.release());
  cout << "empty outside of the function " << const_moved.empty() << endl;
}

#if 0
// doesn't compile!
#define CRAZY_LOOP(x, arr, n)\
  #pragma omp parallel for\
  for (int i=0, b=1; i<n; ++i, b=1)\
    for(x = arr[i]; b; b=0)
  int at[5] = { 1, 2, 3 ,4, 5 };
  CRAZY_LOOP(int x, at, 5)
  {
    printf("parallel %i\n", x);
  }
#endif 

int main(int argc, char **argv)
{ 
  TestArray();
}
