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
#ifndef _OPERATORS_H
#define _OPERATORS_H

#include<functional>
#include "math_ext.h"
#include "helpers-mem.h"

namespace Func
{

template<class T, class S>
struct mult_add : std::binary_function<T,T,void> {
  const S s;
  mult_add(const S &s) : s(s) {}
  void operator()(T& x, const T& y) const {
    x += s * y;
  }
};

template<class T>
struct abs : std::unary_function<T,T> {
  //typedef T result_type;
  T operator()(const T &x) const { return x>=T(0) ? x : -x; }
};

template<class T>
struct cut : std::unary_function<T,T> {
  const T a, b;
  cut(const T &a, const T &b) : a(a), b(b), std::unary_function<T,T>() {}
  T operator()(const T &x) const { return my::clamp<T>(x, a, b); }
};


}

namespace Ops
{
  template<class T, class F>
  struct _Bind2nd {
    _Bind2nd(F f, const T c) : c(c),f(f) {}
    void inplace(T &a) const { f.inplace(a,c); }
    const T operator()(const T &a) const { return f(a,c); }
  private:
    F f;
    const T c;
  };

  template<class T, class F>
  const _Bind2nd<T,F> Bind2nd(F f, const T c) { return _Bind2nd<T,F>(f,c); }

  template<class T, class F, class G>
  struct _Chain {
    _Chain(F f, G g) : f(f),g(g) {}
    void inplace(T &a, const T &b) const { f.inplace(a, g(b)); }
    const T operator()(const T &a, const T &b) const { return f(a, g(b)); }
  private:
    F f; G g;
  };

  template<class T, class F, class G>
  const _Chain<T,F,G> Chain(F f, G g) { return _Chain<T,F,G>(f,g); }

  template<class T>
  struct Construct {
    void inplace(T &a, const T&b) const { ::Construct<T>(&a,b); }
  };

  template<class T>
  struct DefaultConstruct {
    void inplace(T &a) const { ::Construct<T>(&a); }
  };

  template<class T>
  struct Assign {
    void inplace(T &a, const T&b) const { a = b; }    
  };

  template<class T>
  struct Add {
     void inplace(T &a, const T&b) const { a += b; }
     const T operator()(const T &a, const T &b) const { return a + b; }
  };

  template<class T>
  struct Subtract {
    void inplace(T &a, const T&b) const { a -= b; }
    const T operator()(const T &a, const T &b) const { return a - b; }
  };

  template<class T>
  struct Negate {
    void inplace(T &a) const { a = -a; }
    const T operator()(const T &a) const { return -a; }
  };

  template<class T>
  struct Mult {
    void inplace(T &a, const T &b) const { a *= b; }
    const T operator()(const T &a, const T &b) const { return a * b; }
  };
  
  template<class T>
  struct Less {
    bool operator()(const T &a, const T &b) const { return a<b; }    
  };
  
  template<class T>
  struct Greater {
    bool operator()(const T &a, const T &b) const { return a>b; }
  };

  template<class T>
  struct MultAdd {
   MultAdd(const T c) : c(c) {}
   void inplace(T &a, const T&b) const { a += c*b; }
   const T operator()(const T &a, const T &b) const { return a+c*b; }
  private:
   const T c;
  };
  
  template<class T>
  struct MultAddScaleSelf {
   MultAddScaleSelf(const T c, const T c_self) : c(c), c_self(c_self) {}
   void inplace(T &a, const T&b) const { a = a*c_self + c*b; }
   const T operator()(const T &a, const T &b) const { return a*c_self+c*b; }
  private:
   const T c, c_self;
  };
}


#endif
