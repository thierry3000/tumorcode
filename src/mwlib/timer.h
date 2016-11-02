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
#ifndef TIMER_H
#define TIMER_H

#include <sys/time.h>
#include <iostream>

namespace my
{

class Time
{
    timeval t;
  public:
    Time() {
      gettimeofday(&t, NULL);
    }
    Time(const Time &t) : t(t.t) {
    }
    double to_ms() const { return (t.tv_sec * 1000.0) + (t.tv_usec * 0.001); }
    double to_s() const {  return (t.tv_sec * 1.) + (t.tv_usec * 0.001 * 0.001); }
    friend Time operator-(const Time &a, const Time &b) {
      Time r(a);
      r.t.tv_sec -= b.t.tv_sec;
      r.t.tv_usec -= b.t.tv_usec;
      if (r.t.tv_usec < 0) { 
        r.t.tv_usec += 1000000;
        r.t.tv_sec -= 1;
      }
      return r;
    }
    //operator double () const { return to_ms(); }
    void print(std::ostream &os) const { os << to_ms() << "ms"; }
    static Time zero() {
      Time t;
      t.t.tv_sec = 0;
      t.t.tv_usec = 0;
      return t;
    }
};

class Timer
{
    Time tstart;
  public:
    void restart() { tstart = Time(); }
    Time elapsed() { return Time()-tstart; }
};


inline std::ostream& operator<<(std::ostream &os, const Time &t) { t.print(os); return os; }


}

#endif
