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
#ifndef ARRAY_REFCOUNTED_H
#define ARRAY_REFCOUNTED_H

#include "myAssert.h"
#include "helpers-defs.h"
#include "helpers-mem.h"

#include <tbb/spin_mutex.h>
#include <tbb/scalable_allocator.h>
#include <tbb/cache_aligned_allocator.h>

#ifdef DEBUG
  #define ARRAYND_DEBUG
#endif

#ifdef ARRAYND_DEBUG
#  define ARRAYND_ASSERT(exp) myAssert(exp)
#else
#  define ARRAYND_ASSERT(exp) ((void)0)
#endif


struct ArrayNdRef
{
  typedef void (*DestructElements)(void *storage, size_t cnt);

  struct Header {
    size_t cnt, itemsize, cnt_constructed;
    DestructElements destructfunc;
    void *storage;
   
    mutable size_t ref_cnt;
    mutable tbb::spin_mutex mutex;
    tbb::cache_aligned_allocator<char> allocator;

    Header(size_t _cnt, size_t _itemsize, DestructElements _destructFunc, void *extstorage) : cnt(_cnt),itemsize(_itemsize),cnt_constructed(_cnt),destructfunc(_destructFunc),storage(0),ref_cnt(0)
    {
      /*
        if you dont supply the destruct func then the memory will not be freed. The destruct func should call destructors of the classes in the array.
      */
      if (extstorage)
        storage = extstorage;
      else {
        assert(destructfunc != NULL);
        storage = allocator.allocate(cnt*itemsize);
      }
      //cout << "header cons @ " << this << endl;
    }

    ~Header()
    {
      //cout << "header dest @ " << this << endl;
      if(destructfunc)
      {
        destructfunc(storage, cnt_constructed);
        allocator.deallocate((char*)storage, cnt*itemsize);
      }
    }

    static void incRef(const Header* h)
    {
      if (!h) return;
      h->mutex.lock();
      h->ref_cnt += 1;
      h->mutex.unlock();
    }

    static void decRef(const Header* h)
    {
      if (!h) return;
      myAssert(h->ref_cnt>0);
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
  };

  ArrayNdRef() : _header(NULL) {}

  ArrayNdRef(const ArrayNdRef &a) : _header(a._header)
  {
    Header::incRef(_header);
  }

  ArrayNdRef& operator=(const ArrayNdRef &a)
  {
    take(const_cast<ArrayNdRef&>(a).header());
    return *this;
  }

  ~ArrayNdRef()
  {
    Header::decRef(_header);
    _header = NULL;
  }

  void clear()
  {
    Header::decRef(_header);
    _header = NULL;
  }

  void swap(ArrayNdRef &a)
  {
    ::std::swap(a._header,_header);
  }

  void alloc(size_t cnt, size_t itemsize, DestructElements destructFunc)
  {
    Header *h = new Header(cnt, itemsize, destructFunc, NULL);
    take(h);
  }

  void take(Header *h)
  {
    Header::decRef(_header);
    Header::incRef(h);
    _header = h;
  }

  void take(void* mem, size_t cnt, size_t itemsize, DestructElements destructFunc)
  {
    Header *h = new Header(cnt, itemsize, destructFunc, mem);
    take(h);
  }

  void* storage() { return _header ? _header->storage : NULL; }
  const void* storage() const { return _header ? _header->storage : NULL; }

  const Header* header() const { return _header; }
  Header* header() { return _header; }

  // the number of elements the memory allocation can hold
  size_t allocation_size() const { return header() ? header()->cnt : 0; }

  // see if the two reference the same memory
  friend bool operator==(const ArrayNdRef &a, const ArrayNdRef &b)
  {
    bool ok = a.header() == b.header();
    myAssert(!ok || a.storage() == b.storage());
    return ok;
  }

private:
  Header* _header;
};


#endif