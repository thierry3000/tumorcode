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
#include <iostream>
 
 
class PrintCons
{
public:
  PrintCons()
  {
    std::cout << "cons" << std::endl;
  }
  PrintCons(const PrintCons &)
  {
    std::cout << "copy cons" << std::endl;
  }
  
  ~PrintCons()
  {
    std::cout << "dstor" << std::endl;
  }
};
 
 
int main(int argc, char **argv)
{
  PrintCons x;
  std::cout << "-----------" << std::endl;
  boost::any a(x);
  std::cout << "----- addr ------" << std::endl;
  boost::any_cast<PrintCons>(&a);
  std::cout << "----- val ------" << std::endl;
  boost::any_cast<PrintCons>(a);
}