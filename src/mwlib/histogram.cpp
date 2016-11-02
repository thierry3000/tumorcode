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
#include "histogram.h"
//#include "more_helpers.h"
#include <algorithm>

template< class T >
BasicHistogram<T>::BasicHistogram() : HistoBase(),sum_w(0),sum_w2(0),cntref(NULL)
{
}

template< class T >
void BasicHistogram<T>::Init( uint _cnt, ErrorMode _emode )
{
  Clear();
  emode  = _emode;  
  cntref = NULL;
  w.resize( _cnt );
  if( UseError() ) 
  {
    w2.resize( _cnt );
    ResetCountMode( true, NULL );
  }  
  else 
  {
    ResetCountMode( false, NULL );
  }
  sum_w = 0;
  sum_w2 = 0;
}

template< class T >
void BasicHistogram<T>::ResetCountMode( bool _bCount, const BasicHistogram<T>* _cntref )
{
  if( _bCount )
  {
    cntref = _cntref;
    cnt.clear();
    if( !cntref )  cnt.resize( w.size() );
  }
  else
  {
    //myAssert( emode==NoError );
    cnt.clear();
    cntref = NULL;
  }
}


template< class T >
void BasicHistogram<T>::Clear()
{
  w.clear();
  w2.clear();
  cnt.clear();
  sum_w = 0;
  sum_w2 = 0;  
}

template< class T >
void BasicHistogram<T>::Add( const BasicHistogram<T> &src ) // just add up the values
{
  myAssert( src.w.size() == w.size() );
  myAssert( src.emode == emode );
  for( uint i=0; i<w.size(); ++i ) 
  {
    w[i] += src.w[i];
    sum_w += src.w[i];
    if( emode==SquareSum )  // add summed sample values and their squares
    {      
      //myAssert( src.emode==NoError || src.emode==SquareSum );
      //const ValType c = (src.emode==SquareSum) ? src.w2[i] : src.w[i]*src.w[i];
      const ValType c = src.w2[i];
      w2[i] += c;
      sum_w2 += c;
    }
    else if( emode==VarianceSum ) // add summed sample values from src together with their variances
    {
      //myAssert( src.emode==VarianceSum );
      const ValType c = src.w2[i];
      w2[i] += c;
      sum_w2 += c;
    }
  }
  if( !cnt.empty() ) for( uint i=0; i<cnt.size(); ++i ) {
    cnt[i] += src.GetCnt(i);
  }
}

template< class T >
void BasicHistogram<T>::FillWithAverage( const BasicHistogram<T> &src )
{
  myAssert( src.w.size() == w.size() );
  const bool bUseError = UseError();
  const bool bSrcError = src.UseError();

  ValType sum_var_add = 0;
  ValType sum_w_add = 0;
  for( uint i=0; i<w.size(); ++i ) 
  {
    if( src.GetCnt(i)<=0 ) continue;
    
    if( !cnt.empty() ) cnt[i] += 1;
    
    const ValType x = src.GetAvg( i );
    sum_w_add += x;
    w[i] += x;
    //sum_w += x;
    if( emode==SquareSum ) // add a sample wich is the average of samples in src, but don't add the variance from src
    {
      w2[i] += x*x;
      //sum_w2 += x*x;
    }
    else if( emode==VarianceSum ) // add a sample wich is the average of samples in src including the sample variance
    {
      myAssert( src.emode!=NoError );
      const ValType c = src.GetVar( i );
      w2[i] += c;
      //sum_w2 += c*c;
    }
  }

  sum_w += sum_w_add;
  if( emode==SquareSum ) {    
    sum_w2 += sum_w_add*sum_w_add;
  }
  else if( emode==VarianceSum ) {
    sum_w2 += sum_var_add;
  }
}


template< class T >
int BasicHistogram<T>::FillBucket( uint bucket, ValType y, ValType yv )
{
  w[bucket] += y;
  sum_w += y;
  if( emode==SquareSum ) {
    myAssert( yv>=0.0 );
    w2[bucket] += y*y;
    sum_w2 += y*y;
  }
  else if( emode==VarianceSum ) {
    myAssert( yv>0.0 );
    w2[bucket] += yv;
    sum_w += yv;
  }
  if( !cnt.empty() ) cnt[bucket]+=1;
  return bucket;
}


//------------------------------------------------------------------------------

template< class T >
BasicHistogram1D<T>::BasicHistogram1D() : BasicHistogram<T>(),info()
{
}

template< class T >
void BasicHistogram1D<T>::Init( const HistogramInfo<AxisType> &src, HistoBase::ErrorMode _emode, const string &_name )
{
  if( (&src)!=(&info) ) info = src;
	this->name = _name;
  BasicHistogram<ValType>::Init( info.num_buckets, _emode );
}

template< class T >
void BasicHistogram1D<T>::Init( uint _num_buckets, HistoBase::ErrorMode _emode )
{
  Init( HistogramInfo<AxisType>(my::MinMax<AxisType>( 0, _num_buckets ),_num_buckets, false ), _emode );
}

template< class T >
int BasicHistogram1D<T>::Fill( AxisType x, ValType y, ValType yv )
{
  const std::pair<int,bool> b = info.GetBucket( x );
  if( !b.second ) return -1;
  return BasicHistogram<ValType>::FillBucket( b.first, y, yv );
}


//------------------------------------------------------------------------------


template< class T >
BasicHistogram2D<T>::BasicHistogram2D() : infox(),infoy(),BasicHistogram<T>() 
{
}

template< class T >
void BasicHistogram2D<T>::Init( const HistogramInfo<AxisType> &srcx, const HistogramInfo<AxisType> &srcy, HistoBase::ErrorMode _emode, const string &_name  )
{
  if( &srcx!=&infox ) infox = srcx;
  if( &srcy!=&infoy ) infoy = srcy;
  this->name = _name;
  BasicHistogram<ValType>::Init( infox.num_buckets*infoy.num_buckets, _emode );  
}

template< class T >
int BasicHistogram2D<T>::Fill( AxisType x, AxisType y, ValType z, ValType zv )
{
  const std::pair<int,bool> bx = infox.GetBucket( x );
  const std::pair<int,bool> by = infoy.GetBucket( y );
  if( !bx.second || !by.second ) return -1;
  return BasicHistogram<ValType>::FillBucket( bx.first+by.first*infox.num_buckets, z, zv );
}

//------------------------------------------------------------------------------


template< class T >
BasicHistogram3D<T>::BasicHistogram3D() : infox(),infoy(),infoz(),BasicHistogram<T>() 
{
}

template< class T >
void BasicHistogram3D<T>::Init( const HistogramInfo<AxisType> &srcx, 
                                const HistogramInfo<AxisType> &srcy, 
                                const HistogramInfo<AxisType> &srcz, 
                                HistoBase::ErrorMode _emode, const string &_name  )
{
  if( &srcx!=&infox ) infox = srcx;
  if( &srcy!=&infoy ) infoy = srcy;
  if( &srcz!=&infoz ) infoz = srcz;
  this->name = _name;
  BasicHistogram<ValType>::Init( infox.num_buckets*infoy.num_buckets*infoz.num_buckets, _emode );  
}

template< class T >
int BasicHistogram3D<T>::Fill( AxisType x, AxisType y, AxisType z, ValType vv, ValType zv )
{
  const std::pair<int,bool> bx = infox.GetBucket( x );
  const std::pair<int,bool> by = infoy.GetBucket( y );
  const std::pair<int,bool> bz = infoz.GetBucket( z );
  if( !bx.second || !by.second || !bz.second ) return -1;
  return BasicHistogram<ValType>::FillBucket( ToSite(bx.first,by.first,bz.first), vv, zv );
}


//------------------------------------------------------------------------------

template< class T >
HistogramInfo<T>::HistogramInfo() : num_buckets(0),size_bucket(0),eps( std::numeric_limits<T>::epsilon() ), limits()
{
}


template< class T >
void HistogramInfo<T>::Init( const my::MinMax<T> &_limits, uint _num_buckets, bool _bLogScale, const string &_name, T _eps )
{
  bLogScale = _bLogScale;
  limits = _limits;
  eps = _eps;
  num_buckets = _num_buckets;
  name = _name;

  if( bLogScale ) {
    limits.min = LogTrafo( limits.min, eps );
    limits.max = LogTrafo( limits.max, eps );
  }  
  size_bucket = (limits.max-limits.min)/num_buckets;  
  size_bucket_inv = T(1.0)/size_bucket;
}



template<class T>
bool EqualAxis( const HistogramInfo<T> &a, const HistogramInfo<T> &b )
{
  return a.bLogScale==b.bLogScale &&
          a.limits   ==b.limits &&
          a.num_buckets == b.num_buckets &&
          a.eps      ==b.eps;
}

//------------------------------------------------------------------------------

#define STR_QUOTE( s ) "\"" s "\""

template< class HistogramType >
void PrintHeaderMulti( std::ostream &os, 
                        const HistogramType *h,                           
                        const string *labels,
                        const uint cnt )
{
  using namespace std;
  ios_width w(15);

  os << "#" << w;
  if( !labels[0].empty() )  {
    os << quoted(labels[0]) << " ";
  }
  else os << STR_QUOTE("x") << " ";

  for( uint i=1; i<=cnt; ++i ) 
  {
    const HistogramType* hh = h+i-1;
    os << w;
    if( labels && !labels[i].empty() )
      os << quoted(labels[i]) << " ";
    else if( !hh->name.empty() )
      os << quoted(hh->name) << " ";
    else os << boost::format("\"y%i\"") % i << " ";
    if( hh->UseError() ) 
      os << w << STR_QUOTE("RmsDev") << " ";
  }
  os << endl;
}


template< class NamedType > 
struct PrintName
{
  std::ostream &os;
  PrintName( std::ostream &os ) : os(os) {}
  std::ostream& operator()( const NamedType &h ) { 
    os << " " << h.name; 
    if( h.UseError() ) {
      os << " " << h.name << "-rmsd";
    }
    return os; 
  }
};


template< class T >
void PrintHistograms1D( std::ostream &os, 
                        const BasicHistogram1D<T> *h,
                        const uint cnt,
                        const string &prefix,
                        const string &prefixlabel
                         )
{
  using namespace std;
  if( cnt==0 ) return;
  const bool bUseError = h[0].UseError();
  ios_width w(15);
  // basic sanity check
#ifdef DEBUG
  for( uint i=1; i<cnt; ++i ) {
    myAssert( EqualAxis(h[i].info,h[0].info) );
  }
#endif

  os << "#" << prefixlabel;
  os << " " << h[0].info.name;
  std::for_each( h, h+cnt, PrintName<BasicHistogram1D<T> >(os) );
  os << endl;

  // data
  for( uint x=0; x<h[0].info.num_buckets; ++x )
  {
    os << prefix << " " << w << h[0].GetXC( x ) << " ";
    for( uint fi=0; fi<cnt; ++fi )
    {            
      if( h[fi].GetCnt(x)>0 )
      {
        double avg = h[fi].GetAvg( x );
        os << w << avg << " ";
        if( h[fi].UseError() )  {
          double dev = sqrt( std::max<T>( 0.0, h[fi].GetVar( x ) ) );
          os << w << dev  << " ";
        }
      }
      else 
      {
        os << w << "-" << " ";
        if( h[fi].UseError() ) 
          os << w << "-"  << " ";
      }
    }
    os << endl;
  }
}


template< class T >
void PrintHistograms2D( std::ostream &os, 
                        const BasicHistogram2D<T> *h,
                        const uint cnt
                         )
{
  using namespace std;
  if( cnt==0 ) return;
  const bool bUseError = h[0].UseError();
  ios_width w(15);
  // basic sanity check
#ifdef DEBUG
  for( uint i=1; i<cnt; ++i ) {
    myAssert( EqualAxis(h[i].infox,h[0].infox) );
    myAssert( EqualAxis(h[i].infoy,h[0].infoy) );
  }
#endif

  os << "#";
  os << " " << h[0].infox.name;
  os << " " << h[0].infoy.name;
  std::for_each( h, h+cnt, PrintName<BasicHistogram2D<T> >(os) );
  os << endl;

  // data with one extra row/col for correct gp map plots
  for( uint _y=0; _y<h[0].infoy.num_buckets+1; ++_y )
  {
    uint y = std::min(_y,h[0].infoy.num_buckets-1);
    for( uint _x=0; _x<h[0].infox.num_buckets+1; ++_x )
    {
      uint x = std::min(_x,h[0].infox.num_buckets-1);
      os << " " << w << h[0].GetYB( _y ) << " " << w << h[0].GetXB( _x ) << " ";
      for( uint fi=0; fi<cnt; ++fi )
      {            
        if( h[fi].GetCnt(x,y)>0 )
        {
          double avg = h[fi].GetAvg( x, y );
          os << w << avg << " ";
          if( h[fi].UseError() )  {
            double dev = sqrt( std::max<T>( 0.0, h[fi].GetVar( x, y ) ) );
            os << w << dev  << " ";
          }
        }
        else 
        {
          os << w << "0" << " ";
          if( h[fi].UseError() ) 
            os << w << "0"  << " ";
        }
      }
      os << endl;
    }
    os << endl;
  }
}


//------------------------------------------------------------------------------

template struct HistogramInfo<float>;
template class BasicHistogram<float>;
template class BasicHistogram1D<float>;
template class BasicHistogram2D<float>;
template class BasicHistogram3D<float>;

template struct HistogramInfo<double>;
template class BasicHistogram<double>;
template class BasicHistogram1D<double>;
template class BasicHistogram2D<double>;
template class BasicHistogram3D<double>;

template
void PrintHistograms1D<double>( std::ostream &os, 
                        const BasicHistogram1D<double> *h,
                        const uint cnt,
                        const string &prefix,
                        const string &prefixlabel);
template
void PrintHistograms1D<float>( std::ostream &os, 
                        const BasicHistogram1D<float> *h,
                        const uint cnt,
                        const string &prefix,
                        const string &prefixlabel);

template
void PrintHistograms2D<float>( std::ostream &os, 
                        const BasicHistogram2D<float> *h,
                        const uint cnt
                         );
template
void PrintHistograms2D<double>( std::ostream &os, 
                        const BasicHistogram2D<double> *h,
                        const uint cnt
                         );
