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

//#include <mwlib/histogram.h>
#include <H5Cpp.h>

struct BranchDat
{
  float r;
  float l;
  int   n;
};

std::vector<BranchDat> MeasureBranchLengths( const VesselList3d& vl, const std::vector<bool> wanted);
void WriteHdfHistogram( H5::Group f, const string &id, const BasicHistogram1D<float> &h );
const my::Averaged<float> MeasureBranchLengths( const VesselList3d& vl, Histo &distrib, Histo &byrad );