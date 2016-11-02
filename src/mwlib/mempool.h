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
#ifndef _MEMPOOL_H
#define _MEMPOOL_H


#include "helpers-mem.h"
//#include "helpers-defs.h"
#include "helpers-containers.h"

class Image;

/*---------------------------------------*/
// simple segregated storage
/*---------------------------------------*/

class SimpleSegregatedStorage
{
  struct EmptyHeader
  {
    EmptyHeader() : next(NULL) {}
    EmptyHeader* next;
  };
  struct BlockHeader : public DblLLNode<BlockHeader>
  {
    BlockHeader( size_t s, size_t num, size_t binsize ) : 
      DblLLNode<BlockHeader>(), size(s),num(num),nalloc(0),empty_begin(NULL),empty_end(NULL) 
    {
      EmptyHeader *h = (EmptyHeader*)this->Get();
      empty_begin = h;
      for( int i=0; i<size-1; ++i )
      {
        h->next = (EmptyHeader*)(((uchar*)h)+binsize);
        h = h->next;
      }
      h->next = NULL;
      empty_end = h;
    }

    size_t size;
    size_t num;
    size_t nalloc;
    EmptyHeader* empty_begin;
    EmptyHeader* empty_end;
    uchar* Get() { return (uchar*)(this+1); }
  };

  enum {
    INIT_BLOCKSIZE = 8
  };

  BlockHeader* blocks_begin;
  BlockHeader* blocks_end;  
  size_t nblocks;
  size_t binsize;

  BlockHeader* FindBlock( void* p ) const
  {
    BlockHeader *bl = blocks_begin;
    while( bl )
    {
      if( p>=bl->Get() && p<bl->Get()+(bl->size*binsize) ) 
      {
        return bl;
      }
      bl = bl->next;
    }
    return NULL;
  }

  void* GetFromBlock( BlockHeader* bl )
  {
    myAssert( bl->nalloc<bl->size );
    void* ret = bl->empty_begin;
    bl->empty_begin = bl->empty_begin->next;
    if( !bl->empty_begin ) bl->empty_end = NULL;
    bl->nalloc++;
    return ret;
  }

  void ReturnToBlock( BlockHeader* bl, void* _p )
  {
      EmptyHeader* p = (EmptyHeader*)_p;
      if( bl->empty_end ) bl->empty_end->next = p;
      bl->empty_end = p;
      p->next = NULL;
      if( !bl->empty_begin ) bl->empty_begin = p;
      --bl->nalloc;
  }

  size_t SizeFunc( size_t n, size_t s )
  {
    s *= binsize;
    return size_t(s<(1<<26) ? s*2 : s*1.5)/binsize;
  }

  BlockHeader* MakeNewBlock()
  {
    BlockHeader* bl = blocks_end;
    const size_t s = bl ? SizeFunc(bl->num,bl->size) : std::size_t(INIT_BLOCKSIZE);
    bl = new(malloc(s*binsize+sizeof(BlockHeader))) BlockHeader(s,bl ? bl->num+1 : 0, binsize);
    if( !bl ) { myAssert(false); throw std::bad_alloc(); }
    if( !bl ) return NULL;
    bl->InsertLast(blocks_begin,blocks_end);
    ++nblocks;
    return bl;
  }

  void FreeBlock( BlockHeader *bl )
  {
    myAssert(bl->nalloc==0);
    bl->Remove(blocks_begin,blocks_end);
    free( bl );
    --nblocks;
  }


public:

  SimpleSegregatedStorage( size_t allocsize ) : blocks_begin(NULL),blocks_end(NULL),nblocks(0),binsize(allocsize)
  {
    if( binsize < sizeof(EmptyHeader) ) binsize = sizeof(EmptyHeader);
  }

  ~SimpleSegregatedStorage()
  {
    Flush();
  }

  void* Alloc()
  {
    BlockHeader* bl = blocks_begin;
    while( bl )
    {
      if( bl->empty_begin ) break;
      bl = bl->next;
    }  
    if( !bl || !bl->empty_begin )
    {
      bl = MakeNewBlock();      
    }
    return GetFromBlock(bl);
  }

  void Free( void* p )
  {    
    BlockHeader *bl = FindBlock( p );
    myAssert( bl );
    if( bl )
    {
      ReturnToBlock( bl, p );
      if( bl->nalloc<=0 )
      {
        FreeBlock(bl);
      }
    }
  }

  void Flush( int allocsize = -1 );

  size_t allocationSize() const;
  void GetStats( size_t &bins_total, size_t &bins_alloc, size_t &num_blocks );
  void DebugDraw( Image &img, const vector<void*> &allocated );
  void DebugDrawSmall( Image &img );
};




/*----------------------------------------------------------
  Heaplayers Fixed Size Small Object Allocator
  ---------------------------------------------------------*/

#define USE_FIXED_SIZE_ALLOCATOR 0

#if USE_FIXED_SIZE_ALLOCATOR
template<class T>
class FixedSizeAllocator : public SimpleSegregatedStorage
{
  FixedSizeAllocator(const FixedSizeAllocator &);
  FixedSizeAllocator& operator=(const FixedSizeAllocator&);
public:
  enum { SIZE = sizeof(T) };
  typedef SimpleSegregatedStorage SUPER;
public:
  FixedSizeAllocator() : SimpleSegregatedStorage(SIZE) {}
  
  inline T* Alloc() { void*p=SUPER::Alloc(); return new(p) T; }
  inline void Free( T* p ) { if(p) { p->~T(); SUPER::Free(p); } }

  inline void* malloc(size_t s) { myAssert(s==sizeof(T)); return SUPER::Alloc(); }
  inline void free(void* p) { SUPER::Free(p); }
};
#else

template<class T>
class FixedSizeAllocator
{
  FixedSizeAllocator(const FixedSizeAllocator &);
  FixedSizeAllocator& operator=(const FixedSizeAllocator&);
public:
  enum { SIZE = sizeof(T) };
  typedef SimpleSegregatedStorage SUPER;
public:
  FixedSizeAllocator() {}
  
  inline T* Alloc() { return new T; }
  inline void Free( T* p ) { delete(p); }

  inline void* malloc(size_t s) { return ::malloc(s); }
  inline void free(void* p) { return ::free(p); }

  friend void Swap(FixedSizeAllocator<T> &a, FixedSizeAllocator<T> &b) {}
};

#endif

template<class T>
class FixedSizeAllocated 
{
public:
#if USE_FIXED_SIZE_ALLOCATOR
  typedef FixedSizeAllocator<T> Allocator;
  static Allocator* GetAllocator() 
  {    
    static char buf[sizeof(Allocator)];
    static Allocator* ptr = new(buf) Allocator;    
    return ptr;  
  }
  void* operator new( size_t s ) { myAssert(Allocator::SIZE==s);  return GetAllocator()->malloc(s); }
  void operator delete( void* p ) { GetAllocator()->free(p); }
#else
  void* operator new( size_t s ) { return malloc(s); }
  void operator delete( void* p ) { free(p); }
#endif
};

#endif
