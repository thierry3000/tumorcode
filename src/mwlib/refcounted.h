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
#ifndef REFCOUNTED_H
#define REFCOUNTED_H

#include <tbb/spin_mutex.h>
#include <boost/intrusive_ptr.hpp>

#include "myAssert.h"



template<class T>
class StorageRef
{
  T* o;
public:
  StorageRef() : o(NULL) {}
  explicit StorageRef(T &o_) : o(&o_) {}
  operator T& () { return *o; }
};


template<class T>
class StorageCopy
{
  T o;
public:
  StorageCopy() {}
  explicit StorageCopy(T &o_) : o(o) {}
  operator T& () { return o; }
};



class RefCounted
{
  mutable int ref_cnt;
  mutable tbb::spin_mutex mutex;
public:
  RefCounted() : ref_cnt(0) {}
  virtual ~RefCounted() { /*myAssert(cnt==0);*/ }
  const RefCounted& operator=( const RefCounted& x ) { return *this; } // don't overwrite count
  RefCounted( const RefCounted& x ) : ref_cnt(0) {}

  static void incRef(const RefCounted* h)
  {
    if (!h) return;
    h->mutex.lock();
    h->ref_cnt += 1;
    h->mutex.unlock();
  }

  static void decRef(const RefCounted* h)
  {
    if (!h) return;
    h->mutex.lock();
    h->ref_cnt -= 1;
    if(h->ref_cnt<=0)
    {
      h->mutex.unlock();
      delete h;
    }
    else
      h->mutex.unlock();
  }

  static int refCount(const RefCounted* h) { return h->ref_cnt; }
};

namespace boost
{

inline void intrusive_ptr_add_ref(RefCounted *p)
{
  RefCounted::incRef(p);
}

inline void intrusive_ptr_release(RefCounted *p)
{
  RefCounted::decRef(p);
}

}


template<class T>
struct StorageSmartRef
{
  boost::intrusive_ptr<T> o;
public:
  StorageSmartRef() {}
  StorageSmartRef(T &o_, bool manage_ref) : o(&o_)
  {
    if(!manage_ref) {
      RefCounted::incRef(o.get());
      myAssert(RefCounted::refCount(o.get())==2); // the pointer is never freed
    }
  }
  StorageSmartRef(T *o_, bool manage_ref) : o(o_)
  {
    if(!manage_ref) {
      RefCounted::incRef(o.get());
      myAssert(RefCounted::refCount(o.get())==2); // the pointer is never freed
    }
  }
  operator T& () { return *o; }
};



#endif