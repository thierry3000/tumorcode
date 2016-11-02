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


typedef boost::rand48 RNG;
typedef boost::variate_generator<RNG&, boost::uniform_01<double> > Uniform01Generator;

double do_it(const string &prefix)
{
  double res = 0.;
  const int N = 100000;
  const int M = 100;

  double time_per_iter = 0;
  for (int j = 0; j<M; ++j)
  {
    Epetra_SerialComm Comm;
    Epetra_Map Map(N, 0, Comm);

    Epetra_Vector values(Map);
    Epetra_Vector copy(Map);
    
    my::Time _t;
//     #pragma omp parallel shared(values)
//     {
//       boost::uniform_01<double> dist;
//       RNG gen(23466 * omp_get_thread_num());
//       Uniform01Generator gen01(gen, dist);
// 
//       #pragma omp for
//       for (int i=0; i<N; ++i)
//       {
//         double localval = gen01();
//         values[i] = localval;
//       }
//     }
    values.Random();

//     #pragma omp parallel for shared(values, copy)
//     for (int i=0; i<N; ++i)
//     {
//       copy[i] = values[i];
//     }
    copy.Update(2., values, 0.0);

    //double avg = 0.;
//     #pragma omp parallel for reduction(+:avg) shared(copy)
//     for (int i=0; i<N; ++i)
//     {
//       avg += copy[i];
//     }
    double vnorm, cnorm;
    copy.Norm2(&cnorm);
    values.Norm2(&vnorm);
    double avg = 0.5 * ( vnorm + cnorm);
    //avg /= N;
    time_per_iter += (my::Time() - _t).to_ms()*1000. / N;
    if (j > M-5) {
      res = time_per_iter/(j+1);
      printf("%s - %i: %f, (%f us)\n",prefix.c_str(), j, avg, res);
    }
  }
  return res;
}


int main(int arc, char **argv)
{
  //do_it(boost::str(boost::format("np=%i") % omp_get_max_threads()));
  omp_set_num_threads(1);
  double a = do_it("1 thread ");
  omp_set_num_threads(4);
  double b = do_it("4 threads");
  printf("speedup: %f\n", a/b);
  return 1;
}