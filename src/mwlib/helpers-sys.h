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
#ifndef HELPERS_SYS_H
#define HELPERS_SYS_H

#include <iostream>
#include <string>

#if 1
#include "helpers-defs.h"
#include <stdlib.h>
#include <string.h>
#include <string>
#include <sstream>
#include <typeinfo>
#include <exception>
#include <stdexcept>

using std::cerr;
using std::cout;
using std::endl;
using std::ostream;
using std::istream;
using std::string;


uint _GetTimer();
bool _OpenFileChooser( string &sfile );
void _System(const char *command);
void _Sleep(const uint milliseconds);
bool _HasDebugger();
// enum MemUsageMode {
//   MEMUSAGE_CURRENT,
//   MEMUSAGE_PEAK,
//   MEMUSAGE_VMPEAK,
//   MEMUSAGE_VM,
// };

struct MemUsage
{
  uint64 vmem, vmem_peak, rss, rss_peak; // in bytes
};

MemUsage GetMemoryUsage_();

// size_t _GetMemoryUsage(MemUsageMode);
//bool _CopyFile( const string &src, const string &dst );
//bool _MoveFile( const std::string &src, const std::string &dst );
//bool _RemoveFile( const string &fn );

#ifndef INTERNAL

  inline uint GetTimer() { return _GetTimer(); }
  inline bool OpenFileChooser( string &sfile ) { return _OpenFileChooser( sfile ); }
  inline void System(const char *command) { _System(command); }
  inline void Sleep(const uint milliseconds) { _Sleep( milliseconds ); }
  //inline size_t GetMemoryUsage(MemUsageMode m) { return _GetMemoryUsage(m); }
  inline MemUsage GetMemoryUsage() { return GetMemoryUsage_(); }
  inline bool HasDebugger() { return _HasDebugger(); }
//  inline bool CopyFile( const string &src, const string &dst ) { return _CopyFile(src,dst); }
//  inline bool MoveFile( const std::string &src, const std::string &dst ) { return _MoveFile( src,dst); }
//  inline bool RemoveFile( const string &fn ) { return _RemoveFile(fn); }

#endif

bool GetLineFile(void *adr, int &line, char* &file);
bool GetCallstack( void** stack, size_t n );
void* MemDebugAlloc( size_t n );
void MemDebugFree( void *p );

const string FormatMemSize(size_t size);
const string FormatTimeSec( double s );
const string FormatTimeMs( uint ms );
const string strprintf( const char* format, ... );
void DebugOutString( const char* format, ... );
const string RemoveExtension( const string &s );
const string RemoveAllExtensions( const string &s );
const string GetExtension( const string &s );
const string GetPath( const string &s );
const string RemovePath( const string &s );
const string JoinPaths( const string &a, const string &b );
bool SplitPathBack( const string &s, string &a, string &b );

bool   FileExists( const string &s );
bool   PathExists( const string &s, bool bCheckForDir=false );
bool   IsFileGZip( const string &s );
bool   CreatePath( const string &s );


template<class T> inline const string toStr(const T &x)
{
  std::ostringstream os; // much more reasonable than the old implementation 
  os << x;
  return os.str();
}
template<> inline const string toStr<char>( const char &x ) { char b[2] = {x,0}; return string(b); }
template<> inline const string toStr<uchar>( const uchar &x ) { return strprintf("%i",int(x)); }
template<> inline const string toStr<int>( const int &x ) { return strprintf("%i",x); }
template<> inline const string toStr<uint>( const uint &x ) { return strprintf("%u",x); }
template<> inline const string toStr<float>( const float &x ) { return strprintf("%f",x); }
template<> inline const string toStr<double>( const double &x ) { return strprintf("%lf",x); }
template<> inline const string toStr<string>( const string &x ) { return x; }
template<> inline const string toStr<bool>( const bool &x ) { return x ? "true" : "false"; }

template<class T> inline T strTo(const string &s)
{
  std::istringstream is(s);
  T dst; is >> dst;
  if (is.fail()) throw std::invalid_argument(s);
  return dst;
}

//bool strToInt( const char* s, int &x );
// template<class T>
// inline bool FromString( const string& arg, T &dst )
// {
//   std::istringstream is(arg);
//   is >> dst;
//   return !is.fail();
// }
// 
// 
// template<>
// inline bool FromString<bool>( const string& arg, bool &dst )
// {
//   if( 0==arg.compare("true") || 0==arg.compare("1") ) {
//     dst = true;
//     return true;
//   }
//   else if( 0==arg.compare("false") || 0==arg.compare("0") ) {
//     dst = false;
//     return true;
//   }
//   else
//     return false;
// }
// 
// template<>
// inline bool FromString<string>( const string &arg, string &dst )
// {
//   dst.assign(arg);
//   return true;
// }

#endif

struct ios_width {
  unsigned int w;
  ios_width( unsigned int w ) : w(w) {};
};
inline std::ostream& operator<<( std::ostream& os, const ios_width &width ) { os.width(width.w); return os; }

inline const std::string quoted(const std::string &s,const char* left = "\"", const char* right = "\"" ) {
  return std::string(left)+s+std::string(right);
}

#endif
