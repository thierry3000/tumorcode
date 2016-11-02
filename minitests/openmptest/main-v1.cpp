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

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/uniform_int.hpp>


typedef boost::mt19937 RNG;
typedef boost::variate_generator<RNG&, boost::uniform_01<double> > Uniform01Generator;

int main(int arc, char **argv)
{
  const int N = 1000000;
  const int M = 10000;
  double time_per_iter = 0;
  for (int j = 0; j<M; ++j)
  {
    double avg = 0.;
    my::Time _t;  
    #pragma omp parallel  shared(avg)
    {
      boost::uniform_01<double> dist;
      RNG gen(23466 * omp_get_thread_num());
      Uniform01Generator gen01(gen, dist);
      #pragma omp for reduction(+:avg) // reduction scales good
      for (int i=0; i<N; ++i)
      {
        double localval = gen01();
        //#pragma omp atomic // slows stuff down
        avg += localval;
      }
    }
    avg /= N;
    time_per_iter += (my::Time() - _t).to_ms()*1000. / N;
    printf("%i: %f, (%f us)\n",j, avg, time_per_iter/(j+1));
  }
}