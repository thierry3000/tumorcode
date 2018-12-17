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

#include "lattice-data-polymorphic.h"
#include <string>
#include <iostream>

namespace polymorphic_latticedata
{

std::unique_ptr< LatticeData > Make_ld(const char* ldtype, const BBox3& bb, float scale)
{
  //if (strcmp(ldtype, "quad")==0)
  if (strcmp(ldtype, "QUAD3D")==0)
    return std::unique_ptr<LatticeData>(new Derived<LatticeDataQuad3d>(bb, scale));
  //if (strcmp(ldtype, "fcc")==0)
  if (strcmp(ldtype, "FCC")==0)
    return std::unique_ptr<LatticeData>(new Derived<LatticeDataFCC>(bb, scale));
  throw std::invalid_argument(boost::str(boost::format("LatticeData::Make got ldtype %s") % ldtype));
}


template<class LD>
Int3 WorldToLatticeWrapper(const LD &ld, const Float3 &p)
{
  return ld.WorldToLattice(p);
}

template<>
Int3 WorldToLatticeWrapper<LatticeDataFCC>(const LatticeDataFCC &ld, const Float3 &p)
{
  std::runtime_error("WorldToLattice not implemented for FCC lattice");
  return Int3();
}

template Int3 WorldToLatticeWrapper(const LatticeDataQuad3d &ld, const Float3 &p);
template Int3 WorldToLatticeWrapper(const LatticeDataFCC &ld, const Float3 &p); // needs explicit instantiation. It doesn't work automatically for some reason.



template<class LD>
std::unique_ptr<LatticeData> ReadHdfLdGeneric(H5::Group &g)
{
  LD ld;
  ReadHdfLd(g, ld);
  //auto returner = std::unique_ptr<Derived<LD>>(new Derived<LD>(ld));// works on g++ > 5
  std::unique_ptr<LatticeData> returner = std::unique_ptr<LatticeData>(new Derived<LD>(ld));
  returner->print(std::cout);
  return returner;
}


std::unique_ptr<LatticeData> ReadHdf(H5::Group &g)
{
  string type;
  readAttrFromH5(g, string("TYPE"), type);
  if (type == "QUAD3D")
    return ReadHdfLdGeneric<LatticeDataQuad3d>(g);
  else if (type == "FCC")
    return ReadHdfLdGeneric<LatticeDataFCC>(g);
  else
    throw std::runtime_error(boost::str(boost::format("unknown lattice data type %s in hdf file") % type));
}

template <class Ld>
void Derived<Ld>::Lattice2Hdf(H5::Group& g) const
{
  ld.WriteHdfLd(g);
}

template <class Ld>
Derived<Ld>::Derived(const Ld& _ld): ld(_ld)
{
  ld.print(std::cout);
}

template <class Ld>
string Derived<Ld>::GetType() const
{
  return string("FCC");
}


}//polymorphic_latticedata
