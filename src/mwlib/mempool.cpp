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

#include <iostream>
#include <boost/format.hpp>

#include "mempool.h"
#include "drawing.h"
#include "helpers-vec.h"
#include "random.h"

typedef SimpleSegregatedStorage SSG;


void SimpleSegregatedStorage::Flush( int allocsize )
{
  while( blocks_begin )
  {
    BlockHeader* next = blocks_begin->Next();
    free( (void*)blocks_begin );
    blocks_begin = next;
  }
  blocks_begin = blocks_end = NULL;
  nblocks = 0;
  if( allocsize>0 ) {
    binsize = std::max<size_t>(sizeof(EmptyHeader),allocsize);
  }
}


void SimpleSegregatedStorage::GetStats( size_t &bins_total, size_t &bins_alloc, size_t &_num_blocks )
{
  bins_total = bins_alloc = 0;
  BlockHeader* bl = blocks_begin;
  for(;bl;bl=bl->Next())
  {
    bins_total += bl->size;
    bins_alloc += bl->nalloc;
  }
  _num_blocks = nblocks;
}

size_t SimpleSegregatedStorage::allocationSize() const
{
  size_t ret = sizeof(SimpleSegregatedStorage);
  BlockHeader* bl = blocks_begin;
  for(;bl;bl=bl->Next())
  {
    ret += bl->size*binsize+sizeof(BlockHeader);
  }
  return ret;
}

void SimpleSegregatedStorage::DebugDraw( Image &img, const vector<void*> &allocated )
  {
    size_t s = 0;
    BlockHeader* bl = blocks_begin;
    while( bl )
    {
      s += bl->size;
      bl = bl->next;
    }
    const int cellw = 16;
    size_t w = cellw*s+nblocks*4;
    if(w>MAX_IMG_SIZE) return;
    img.Init( int(w), cellw );
    bl = blocks_begin;
    int xoff = 2;
    while( bl )
    {
      img.SetColor( 200,100,100 );
      img.DrawRect( xoff, 2, int(xoff+bl->size*cellw), cellw-2 );            
      EmptyHeader* h = bl->empty_begin;
      int i=0;
      while( h )
      {
        int off = int((((uchar*)h)-bl->Get())/binsize)*cellw;
        img.SetColor( 0,0,0 );
        img.DrawRect( xoff+off+1, 3, xoff+off+cellw-1, cellw-3 );
        img.SetColor( 100,100,100 );
        img.DrawText( xoff+off+2, 4, boost::str(boost::format("%i") % i).c_str() );
        h = h->next;
        ++i;
      }

      for( int i=0; i<allocated.size(); ++i )
      {
        if( allocated[i]==NULL ) continue;
        BlockHeader* bbl = FindBlock(allocated[i]);
        if( bbl != bl ) continue;
        int off = int((((uchar*)allocated[i])-bl->Get())/binsize)*cellw;
        img.SetColor( 255,0,0 );
        img.DrawRect( xoff+off+2, 4, xoff+off+6, 8 );
      }

      xoff += int(bl->size*cellw)+4;
      bl = bl->next;
    }    
  }

void SimpleSegregatedStorage::DebugDrawSmall( Image &img )
{
  uint bcnt = 0;
  size_t cnt = 0;
  BlockHeader* bl = blocks_begin;
  for(;bl;bl=bl->Next())
  {
    ++bcnt;
    cnt += bl->size;
  }
  size_t ss = 2;
  size_t tot_cnt = cnt+4*bcnt;
  while( (ss)*(ss) <= tot_cnt ) ++ss;
  if(ss>MAX_IMG_SIZE) return;
  int s = int(ss);
  Random rnd;
  img.Init(s,s);
  uint i=0;
  bl = blocks_begin;
  img.Fill(128,128,128);
  for(;bl;bl=bl->Next())
  {
    img.SetColor(0,255,0);    
    img.DrawPixel(i%s,i/s);
    ++i;
    const Float3 col = Vec<float,3>(0.5f)+0.5f*Vec<float,3>(rnd.Get01(),rnd.Get01(),rnd.Get01());
    img.SetColor(col);
    for( uint x=i; x<i+bl->size; ++x )
    {
      myAssert(x<tot_cnt);
      img.DrawPixel(x%s,x/s);
    }
    EmptyHeader *h = bl->empty_begin;
    img.SetColor(col-Vec<float,3>(0.25));
    for(;h;h=h->next)
    {
      myAssert((uchar*)h>(uchar*)bl && (uchar*)h<((uchar*)bl)+bl->size*binsize);
      uint x = i+uint(((uchar*)h-bl->Get())/binsize);
      myAssert(x<tot_cnt);
      img.DrawPixel(x%s,x/s);
    }
    i += uint(bl->size);
    img.SetColor(0,255,0);
    img.DrawPixel(i%s,i/s);
    ++i;
  }
}
