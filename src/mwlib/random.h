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
#ifndef RANDOM_H
#define RANDOM_H


#include "helpers-defs.h"
#include "math_ext.h"

#if 0
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/uniform_int.hpp>

class Random
{
public:
  //typedef boost::mt19937 RNG;
  typedef boost::mt11213b RNG;
  typedef boost::variate_generator<RNG&, boost::uniform_01<double> > Uniform01Generator;
  RNG rng;
private:
  boost::uniform_01<double> dist_uniform01;
  Uniform01Generator gen01;
public:
  Random(uint seed = 123456) : rng(seed), gen01(rng, dist_uniform01) {}
  void Init(uint seed) { rng = RNG(seed); }

  double Get01() { return gen01(); }
  float Get01f() { return float(Get01()); }
  double Get11() { return Get01() * 2. - 1.; }
  float Get11f() { return float(Get11()); }
  int Get(int max)
  {
    boost::uniform_int<> dist(0, max-1);
    return boost::variate_generator<RNG&, boost::uniform_int<> >(rng, dist)();
  }
  uint Get() { return rng(); }

  template <class Map>
  boost::variate_generator<RNG&, Map> var_gen(const Map &map)
  {
    return boost::variate_generator<RNG&, Map>(rng);
  }
};
#endif

// Mersenne Twister
class Random
{
public:
  Random( uint seed = 0 ) { Init(seed); }
  void Init( uint seed );
  inline float Get01f()
  {
    uint y;
    if (--left == 0) NextState();
    //y = *next++;
    y = state[inext++];
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);
    return (float)y * (1.0f/4294967295.0f);
  }
  inline double Get01()
  {
    uint y;
    if (--left == 0) NextState();
    //y = *next++;
    y = state[inext++];
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);
    return (double)y * (1.0/4294967295.0);
  }

  inline uint Get()
  {
    uint y;
    if (--left == 0) NextState();
    //y = *next++;
    y = state[inext++];
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);
    return y;
  }
  
  inline double Get11()
  {
    return Get01()*2.0-1.0;
  }
  inline float Get11f()
  {
    return Get01f()*2.0f-1.0f;
  }
  inline int Get( int max ) // apparently not including "max"
  {
    uint y;
    if (--left == 0) NextState();
    y = state[inext++];
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);
    return int(((double)y * (1.0/4294967296.0))*max);
  }

private:
  void NextState();
  enum {
    N = 624
  };
  int left;
  int inext;
  uint state[N];
};


class Halton
{
  int prime;
  float last;
  float b;
  int i;
public:
  Halton() {}

  Halton(int p, int i = 0)
  {
    prime = p;
    b = 1.0f/prime;
    //i = _i;
    last = Instance(i,p);
    this->i = 0;
  }
  float operator()()
  {
//     float h,hh,r;
//     r = 1.0f-1.0e-10f-last;
//     if( b<r )
//       last += b;
//     else
//     {
//       h = b;
//       do {
//         hh = h;
//         h *= b;
//       }
//       while( h>=r );
//       last += hh + h -1.0f;
//     }
//     return last;
    return Instance(i++, prime);
  }
  static float Instance(int i, int prime)
  {
    i++;
    int j = i;
    int d = 1;
    int n = 0;
    while (j) {
      n = n * prime + j % prime;
      d = d * prime;
      j = j / prime;
    }
    return (n/float(d));
  }
};


template<int dim>
class HaltonSequence
{
    Halton seq[dim];
  public:
    HaltonSequence()
    {
      for (int i=0; i<dim; ++i)
        seq[i] = Halton(my::mconst::primes[i], 0);
    }

    Vec<float,dim> Get01f()
    {
      Vec<float,dim> r;
      for (int i=0; i<dim; ++i)
        r[i] = seq[i]();
      return r;
    }

    Vec<float,dim> Get11f()
    {
      Vec<float,dim> r;
      for (int i=0; i<dim; ++i)
        r[i] = 2.f*(seq[i]()-0.5f);
      return r;
    }
};


template<int dim>
inline Vec<float, dim> GetJittered01f(HaltonSequence<dim> &qrnd, Random &rnd, int num_samples)
{
  float ampl = 0.1 * pow(float(num_samples), float(-1./dim));
  Vec<float, dim> q(qrnd.Get01f());
  for (int i=0; i<dim; ++i) {
    q[i] = my::cut<float>(q[i] + ampl*rnd.Get11f(), 0., 1.);
  }
  return q;
}




#if 0
struct TMSNet
{
  typedef unsigned int uint;
  static const float DENOM;
  static inline float RI_vdC(uint bits, uint r = 0)
  {
    bits = ( bits               << 16) | ( bits               >> 16);
    bits = ((bits & 0x00ff00ff) <<  8) | ((bits & 0xff00ff00) >>  8);
    bits = ((bits & 0x0f0f0f0f) <<  4) | ((bits & 0xf0f0f0f0) >>  4);
    bits = ((bits & 0x33333333) <<  2) | ((bits & 0xcccccccc) >>  2);
    bits = ((bits & 0x55555555) <<  1) | ((bits & 0xaaaaaaaa) >>  1);
    bits ^= r;
    return (float)bits * DENOM;
  }

  static inline float RI_S(uint i, uint r = 0)
  {
    for (uint v = (unsigned int)1<<31; i; i >>= 1, v ^= v>>1)
      if (i & 1)
  r ^= v;
    return (float)r * DENOM;
  }

  static inline float RI_LP(uint i, uint r = 0)
  {
    for (uint v = (unsigned int)1<<31; i; i >>= 1, v |= v>>1)
      if (i & 1)
  r ^= v;
    return (float)r * DENOM;
  }
};
#endif


#endif
