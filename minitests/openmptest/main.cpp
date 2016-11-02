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
#include <omp.h>
#include <stdio.h>
#include "timer.h"
#include <string>
#include <vector>
#include <math.h>
#include <fstream>

using std::string;
using std::cout;
using std::endl;

#include <boost/random/linear_congruential.hpp>
//#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/format.hpp>

#include "Epetra_SerialComm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Version.h"
#include "Epetra_CrsMatrix.h"

#include "timer.h"

typedef boost::rand48 RNG;
typedef boost::variate_generator<RNG&, boost::uniform_01<double> > Uniform01Generator;

enum TimingOp {
  Construction,
  Random,
  Update_xb,
  Update_xbya,
  Norm2,
  Dot,
  Abs,
  Reciprocal,
  CrsMatrixBuilding,
  CrsMatrixMult,
  PutScalar,
  NUM_OP
};
#define DECL_OPSTR(s) #s

static const char* opstr[] =
{
  DECL_OPSTR(Construction),
  DECL_OPSTR(Random),
  DECL_OPSTR(Update_xb),
  DECL_OPSTR(Update_xbya),
  DECL_OPSTR(Norm2),
  DECL_OPSTR(Dot),
  DECL_OPSTR(Abs),
  DECL_OPSTR(Reciprocal),
  DECL_OPSTR(CrsMatrixBuilding),
  DECL_OPSTR(CrsMatrixMult),
  DECL_OPSTR(PutScalar)
};

#define TSTART(id) t_ = my::Time();
#define TSTOP(id) timings[id] += (my::Time()-t_).to_ms();

std::vector<double> do_it(const string &prefix, bool manual, const int N)
{
  double res = 0.;
  //const int N = 100000;
  const int M = 10;

  std::vector<double> timings(NUM_OP, 0.);
  my::Time t_;

  double time_per_iter = 0;
  for (int j = 0; j<M; ++j) 
  {
    Epetra_SerialComm Comm;
    Epetra_Map Map(N, 0, Comm);

    TSTART(Construction)
    Epetra_Vector values(Map, false);
    Epetra_Vector copy(Map, false);
    TSTOP(Construction)

    TSTART(Random)
    values.Random();
    TSTOP(Random)

    TSTART(CrsMatrixBuilding)
    Epetra_CrsMatrix A(Copy, Map, 3);
    {
    // Add  rows one-at-a-time
    double negOne = -1.0;
    double posTwo = 2.0;
    for (int i=0; i<N; i++) {
      int GlobalRow = A.GRID(i); int RowLess1 = GlobalRow - 1; int RowPlus1 = GlobalRow + 1;

      if (RowLess1!=-1) A.InsertGlobalValues(GlobalRow, 1, &negOne, &RowLess1);
      if (RowPlus1!=N) A.InsertGlobalValues(GlobalRow, 1, &negOne, &RowPlus1);
      A.InsertGlobalValues(GlobalRow, 1, &posTwo, &GlobalRow);
    }
    // Finish up
    A.FillComplete();
    }
    TSTOP(CrsMatrixBuilding)

    TSTART(PutScalar)
    values.PutScalar(1.);
    TSTOP(PutScalar)

    TSTART(Update_xb)
    copy.Update(2., values, 0.0);
    TSTOP(Update_xb)

    TSTART(Update_xbya)
    values.Update(3., copy, 5.);
    TSTOP(Update_xbya)

    double vnorm = 0., cnorm = 0.;
    TSTART(Norm2)
    if (manual)
    {
      #pragma omp parallel for reduction(+:cnorm) shared(copy)
      for (int i=0; i<N; ++i)
      {
        cnorm += copy[i] * copy[i];
      }
      cnorm = std::sqrt(cnorm);
      #pragma omp parallel for reduction(+:vnorm) shared(values)
      for (int i=0; i<N; ++i)
      {
        vnorm += values[i] * values[i];
      }
      vnorm = std::sqrt(vnorm);
    }
    else
    {
      copy.Norm2(&cnorm);
      values.Norm2(&vnorm);
    }
    TSTOP(Norm2)
 
//     TSTART(Abs)
//     values.Abs(copy);
//     TSTOP(Abs)

    TSTART(Reciprocal)
    copy.Reciprocal(values);
    TSTOP(Reciprocal)

    TSTART(CrsMatrixMult) 
    A.Multiply(false, copy, values);
    TSTOP(CrsMatrixMult)

      double dot = 0; 
    TSTART(Dot)
    values.Dot(copy, &dot);
    TSTOP(Dot)

    double avg = 0.5 * ( vnorm + cnorm + dot);
    if (j > M - 2)
      printf("some result %s: %f\n",prefix.c_str(), avg);
  }
  for (int o=0; o<timings.size(); ++o) {
    if (timings[o]>0)
      timings[o] /= M;
    else
      timings[o] = 1.;
  }
  return timings;
}


int main(int arc, char **argv)
{
  int sizes[] = {
    1000000,
     800000,
     500000,
     200000,
     100000,
      80000,
      50000,
      20000,
      10000,
       8000,
       5000,
       2000,
       1000,
        800,
        500,
        200,
        100
  };
  std::ofstream of("perdata2.dat");
  for (int j=0; j<sizeof(sizes)/sizeof(int); ++j)
  {
    omp_set_num_threads(1); 
    std::vector<double> t1 = do_it("", true, sizes[j]);
    omp_set_num_threads(3); 
    std::vector<double> t2 = do_it("", true, sizes[j]);
    for (int i=0; i<NUM_OP; ++i)
      printf("%4i: np1: %.2f, np3: %.2f, (%.2fx)\n", sizes[j], t1[i], t2[i], t1[i]/t2[i]);
    printf("-------");
    of << sizes[j] << " " << (t1[CrsMatrixMult]/t2[CrsMatrixMult]) << " " << t2[CrsMatrixMult] << " " << t1[CrsMatrixMult] << endl;
  }
  return 1;
}