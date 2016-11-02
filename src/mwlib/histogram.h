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
#ifndef HISTOGRAMS_H
#define HISTOGRAMS_H

#include "helpers.h"
#include "math_ext.h"

//------------------------------------------------------------------------------


template< class T >
inline T LogTrafo( T x, T eps ) { return log(x+eps); }
template< class T >
inline T ExpTrafo( T x, T eps ) { return exp(x)-eps; }


template< class T >
struct HistogramInfo
{
  typedef T ValType;
  my::MinMax<T> limits;
  uint num_buckets;
  T size_bucket;
  T size_bucket_inv;
  T eps;
  bool bLogScale;
  string name;

  HistogramInfo();  
  HistogramInfo( const my::MinMax<T> &_limits, uint _num_buckets, bool _bLogscale, 
                 const string &_name = string(), T _eps = std::numeric_limits<T>::epsilon() ) { Init(_limits,_num_buckets,_bLogscale,_name,_eps); }
  HistogramInfo( T min, T max, uint _num_buckets, bool _bLogscale, 
                 const string &_name = string(), T _eps = std::numeric_limits<T>::epsilon() ) { Init(my::MinMax<T>(min,max),_num_buckets,_bLogscale,_name,_eps); }
  HistogramInfo( T min, T max, T _size_bucket, bool _bLogscale, 
                 const string &_name = string(), T _eps = std::numeric_limits<T>::epsilon() ) { Init(my::MinMax<T>(min,max),my::iceil((max-min)/_size_bucket),_bLogscale,_name,_eps); }
  void Init( const my::MinMax<T> &limits, uint num_buckets, bool bLogscale, const string &name = string() , T eps = std::numeric_limits<T>::epsilon() );
  inline T GetBucketCenter( uint i ) const;
  inline T GetBucketBegin( uint i ) const;
  inline std::pair<int,bool> GetBucket( T x ) const;
  inline bool GetBucketFract( T x, int &ix, float &fx ) const;
};

template<class T>
bool EqualAxis( const HistogramInfo<T> &a, const HistogramInfo<T> &b );


//------------------------------------------------------------------------------

class HistoBase
{
public:
  enum ErrorMode {
    NoError = 0,
    Normalized = 1,
    VarianceSum = 2,
    SquareSum = 4
  };
  HistoBase() : id(0),emode(NoError) {}
  inline bool UseError() const { return emode&(VarianceSum|SquareSum); }
  inline ErrorMode GetErrorMode() const { return emode; }
  inline bool IsNormalized() const { return emode==Normalized; }

  uint            id;
  string     name;
  ErrorMode       emode;
};

template< class T >
class BasicHistogram : public HistoBase
{
public:
  typedef T ValType;
  typedef uint CntType;

  BasicHistogram();
  void Clear();
  void Init( uint cnt, ErrorMode emode );
  void ResetCountMode( bool bCount, const BasicHistogram<T>* cntref = NULL );

  void Add( const BasicHistogram<T> &src );
  void FillWithAverage( const BasicHistogram<T> &src );  
  
  int FillBucket( uint bucket, ValType y, ValType yv=0.0 );  
  
  inline ValType GetAvg( uint i ) const;
  inline ValType GetVar( uint i ) const;
  inline CntType GetCnt( uint i ) const;  
   
  vector<ValType> w;
  vector<ValType> w2;
  ValType         sum_w;
  ValType         sum_w2;

  const BasicHistogram<T>*   cntref;
  vector<CntType> cnt;
};


template< class T >
class BasicHistogram1D : public BasicHistogram<T>
{ 
public:
  typedef double AxisType;  
  typedef HistogramInfo<AxisType> HInfo;
	typedef T ValType;
  BasicHistogram1D();
  void Init( const HistogramInfo<AxisType> &src, HistoBase::ErrorMode emode, const string &name = string() );
  void Init( uint num_buckets, HistoBase::ErrorMode emode );  
  int Fill( AxisType x, ValType y, ValType yv=0.0 );
  inline AxisType GetXC( uint i ) const { return info.GetBucketCenter(i); }
  inline AxisType GetXB( uint i ) const { return info.GetBucketBegin(i); }
  inline uint Size() const { return info.num_buckets; } 
  
  HistogramInfo<AxisType> info;
};



template< class T >
class BasicHistogram2D : public BasicHistogram<T>
{
public:
  typedef double AxisType;
	typedef T ValType;
  typedef HistogramInfo<AxisType> HInfo;
	typedef typename BasicHistogram<T>::CntType CntType;

  BasicHistogram2D();
  void Init( const HistogramInfo<AxisType> &srcx, const HistogramInfo<AxisType> &srcy, HistoBase::ErrorMode emode, const string &name = string() );

  int Fill( AxisType x, AxisType y, ValType z, ValType zv=ValType(0) );
  inline AxisType GetXC( uint i ) const { return infox.GetBucketCenter(i); }
  inline AxisType GetXB( uint i ) const { return infox.GetBucketBegin(i);  }
  inline AxisType GetYC( uint i ) const { return infoy.GetBucketCenter(i); }
  inline AxisType GetYB( uint i ) const { return infoy.GetBucketBegin(i);  }
  inline uint SizeX() const { return infox.num_buckets; } 
  inline uint SizeY() const { return infoy.num_buckets; } 
  inline uint ToSite( uint x, uint y ) const { return y*infox.num_buckets+x; }

  inline ValType GetAvg( uint x, uint y ) const { return BasicHistogram<T>::GetAvg(ToSite(x,y)); }
  inline ValType GetVar( uint x, uint y ) const { return BasicHistogram<T>::GetVar(ToSite(x,y)); }
  inline CntType GetCnt( uint x, uint y ) const { return BasicHistogram<T>::GetCnt(ToSite(x,y)); }
  inline ValType GetAvg( uint i ) const { return BasicHistogram<T>::GetAvg(i); }
  inline ValType GetVar( uint i ) const { return BasicHistogram<T>::GetVar(i); }
  inline CntType GetCnt( uint i ) const { return BasicHistogram<T>::GetCnt(i); }  


  HistogramInfo<AxisType> infox;
  HistogramInfo<AxisType> infoy;
};

template< class T >
void PrintHistograms1D( std::ostream &os, 
                        const BasicHistogram1D<T> *h,
                        const uint cnt,
                        const string &prefix = string (),
                        const string &prefixlabel = string() );

template< class T >
void PrintHistograms2D( std::ostream &os, 
                        const BasicHistogram2D<T> *h,
                        const uint cnt );

//------------------------------------------------------------------------------

template< class T >
class BasicHistogram3D : public BasicHistogram<T>
{
public:
  typedef double AxisType;
	typedef T ValType;
  typedef HistogramInfo<AxisType> HInfo;
	typedef typename BasicHistogram<T>::CntType CntType;

  BasicHistogram3D();
  void Init( const HistogramInfo<AxisType> &srcx, const HistogramInfo<AxisType> &srcy, const HistogramInfo<AxisType> &srcz, HistoBase::ErrorMode emode, const string &name = string() );

  int Fill( AxisType x, AxisType y, AxisType z, ValType vv, ValType zv=0.0 );
  inline AxisType GetXC( uint i ) const { return infox.GetBucketCenter(i); }
  inline AxisType GetXB( uint i ) const { return infox.GetBucketBegin(i);  }
  inline AxisType GetYC( uint i ) const { return infoy.GetBucketCenter(i); }
  inline AxisType GetYB( uint i ) const { return infoy.GetBucketBegin(i);  }
  inline AxisType GetZC( uint i ) const { return infoz.GetBucketCenter(i); }
  inline AxisType GetZB( uint i ) const { return infoz.GetBucketBegin(i);  }
  inline uint SizeX() const { return infox.num_buckets; } 
  inline uint SizeY() const { return infoy.num_buckets; } 
  inline uint SizeZ() const { return infoz.num_buckets; }
  inline uint ToSite( uint x, uint y, uint z ) const { return z*infox.num_buckets*infoy.num_buckets+y*infox.num_buckets+x; }

  inline ValType GetAvg( uint x, uint y, uint z ) const { return BasicHistogram<T>::GetAvg(ToSite(x,y,z)); }
  inline ValType GetVar( uint x, uint y, uint z ) const { return BasicHistogram<T>::GetVar(ToSite(x,y,z)); }
  inline CntType GetCnt( uint x, uint y, uint z ) const { return BasicHistogram<T>::GetCnt(ToSite(x,y,z)); }
  inline ValType GetAvg( uint i ) const { return BasicHistogram<T>::GetAvg(i); }
  inline ValType GetVar( uint i ) const { return BasicHistogram<T>::GetVar(i); }
  inline CntType GetCnt( uint i ) const { return BasicHistogram<T>::GetCnt(i); }  


  HistogramInfo<AxisType> infox;
  HistogramInfo<AxisType> infoy;
  HistogramInfo<AxisType> infoz;
};


//------------------------------------------------------------------------------


template< class T >
T HistogramInfo<T>::GetBucketCenter( uint i ) const
{ 
  T x = limits.min + size_bucket*(i+0.5f); 
  if( bLogScale ) x = ExpTrafo<T>( x, eps );
  return x;
}

template< class T >
T HistogramInfo<T>::GetBucketBegin( uint i ) const
{
  T x = limits.min + size_bucket*(i); 
  if( bLogScale ) x = ExpTrafo<T>( x, eps );
  return x;
}

template< class T >
std::pair<int,bool> HistogramInfo<T>::GetBucket( T x ) const
{  
  if( bLogScale ) x = LogTrafo<T>( x, eps );
  const T f = my::normalize<T>(x, limits.min, limits.max );
  const int indx = (int)(f*num_buckets);
  return std::make_pair(indx,((uint)indx)<num_buckets);
}

template< class T >
bool HistogramInfo<T>::GetBucketFract( T x, int &ix, float &fx  ) const
{  
  if( bLogScale ) x = LogTrafo<T>( x, eps );
  const T f = my::normalize<T>(x, limits.min, limits.max )*num_buckets;
  ix = std::min<int>((int)(f),num_buckets-1);
  fx = float(f-ix);
  return f>=T(0) && f<=T(num_buckets);
}

//------------------------------------------------------------------------------

template< class T >
typename BasicHistogram<T>::CntType BasicHistogram<T>::GetCnt( uint i ) const
{
  if( cntref || !cnt.empty() )
  {
    myAssert( !IsNormalized() );
    return cntref ? cntref->cnt[i] : cnt[i];
  }
  if( IsNormalized() )
  {
    // what the fuck is this?????
    //myAssert(false);
    //return 1;
    return sum_w;
  }
  return 1;
}

template< class T >
typename BasicHistogram<T>::ValType BasicHistogram<T>::GetVar( uint i ) const
{  
  if( UseError() )
  {
    const CntType cnt = GetCnt(i);
    if( emode==SquareSum ) return (w2[i] - w[i]*w[i]/cnt)/cnt;
    else return (w2[i]/(cnt*cnt));
  }
  else return 0;
}

template< class T >
typename BasicHistogram<T>::ValType BasicHistogram<T>::GetAvg( uint i ) const 
{ 
  //if( emode )
  //{
    const CntType cnt = GetCnt(i);
    return (1.0f/cnt)*w[i]; 
  //}
  //else return w[i];
}


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

template< class T, class R >
struct LookupTable2D
{
  vector<T> data;
  HistogramInfo<R> xAxis;
  HistogramInfo<R> yAxis;
  uint nx,ny;
  void Init( R xmin, R xmax, uint nx, bool bLogScaleX,
             R ymin, R ymax, uint ny, bool bLogScaleY );
  T DoLookup( R x, R y ) const;
  T &operator()( uint i, uint j ) { return data[j*nx+i]; }
  T operator()( uint i, uint j ) const { return data[j*nx+i]; }
};

template< class T, class R >
void LookupTable2D<T,R>::Init( R xmin, R xmax, uint nx, bool bLogScaleX,
                                 R ymin, R ymax, uint ny, bool bLogScaleY )
{
  data.clear();
  this->nx = nx;
  this->ny = ny;
  xAxis.Init( my::MinMax<R>(xmin,xmax), nx-1, bLogScaleX );
  yAxis.Init( my::MinMax<R>(ymin,ymax), ny-1, bLogScaleY );
  data.resize( nx*ny );
}

template< class T, class R >
T LookupTable2D<T,R>::DoLookup( R x, R y ) const
{
  T vx0,vx1;
  int iy,ix;
  R fx,fy;

  if( xAxis.bLogScale ) x = LogTrafo<T>( x, xAxis.eps );
  if( yAxis.bLogScale ) y = LogTrafo<T>( y, yAxis.eps );
  x = (x-xAxis.limits.min)*xAxis.size_bucket_inv;
  y = (y-yAxis.limits.min)*yAxis.size_bucket_inv;

  if( x<xAxis.num_buckets && x>=0.0 )
  {
    ix = (int)x;
    fx = x-ix;
  }
  else
  {
    if( x < 0.0 ) {
      ix = 0;
      fx = 0;
    }
    else {
      ix = xAxis.num_buckets-1;
      fx = 1;
    }
  }

  if( y<yAxis.num_buckets && y>=0.0 )
  {
    iy = (int)y;
    fy = y-iy;
  }
  else
  {
    if( y<0.0 ) {
      iy = 0;
      fy = 0;
    }
    else {
      iy = yAxis.num_buckets-1;
      fy = 1;
    }
  }

  vx0 = data[ (iy+0)*nx + ix ]*(1.0-fx) + (fx)*data[ (iy+0)*nx + ix+1 ];
  vx1 = data[ (iy+1)*nx + ix ]*(1.0-fx) + (fx)*data[ (iy+1)*nx + ix+1 ];
  return vx0*(1.0-fy)+fy*vx1;
}



#endif
