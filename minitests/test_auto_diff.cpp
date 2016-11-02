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
 
//#include <Sacado_Fad_SFad.hpp>
#include "Sacado.hpp"
#include <iostream>


bool g_sw = true;

template<class T>
T func(const T &a, const T &b, const T &c) {
  T q = a*a;
  T r = b*b;
  T s = c*c;
  T ret = 10. * q + r + s;
  if (g_sw)
    ret += a*b*c;
  return ret;
  //return a*a + b*b + a*b*c + c*c;
}


int main(int argc, char **argv)
{
  g_sw = argc>1 && argv[1][0]=='1';
  double a = 1, b = 2, c = 3;
  
  double h = 1.e-6;
  double r = func(a, b, c);
  double drda = (func(a+h, b, c)-func(a-h, b, c))*0.5/h;
  double drdb = (func(a, b+h, c)-func(a, b-h, c))*0.5/h;
  std::cout << r << " " << drda << " " << drdb << std::endl;  

#if 0
  Sacado::Fad::DFad<double> afad(2, 0, a);
  Sacado::Fad::DFad<double> bfad(2, 1, b);
  Sacado::Fad::DFad<double> cfad(c);
  Sacado::Fad::DFad<double> rfad;
  rfad = func(afad, bfad, cfad);
  
  double r_ad = rfad.val();
  double drda_ad = rfad.dx(0);
  double drdb_ad = rfad.dx(1);
#else
  Sacado::Fad::SFad<double,2> afad(2, 0, a);
  Sacado::Fad::SFad<double,2> bfad(2, 1, b);
  Sacado::Fad::SFad<double,2> cfad(c);
  Sacado::Fad::SFad<double,2> rfad;
  
  rfad = func(afad, bfad, cfad);
  double r_ad = rfad.val();
  double drda_ad = rfad.dx(0);
  double drdb_ad = rfad.dx(1);
#endif
  std::cout << r_ad << " " << drda_ad << " " << drdb_ad << std::endl;
  return 0;
}