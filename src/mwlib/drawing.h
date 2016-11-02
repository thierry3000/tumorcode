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
#ifndef DRAWING_H_
#define DRAWING_H_

//#include "helpers.h"
#include <string>
#include <vector>


#if (((defined __GNUC__) && (defined __LP64__)) || (defined _WIN64))
  #define IMG_PRIVATE_DATA_SIZE 32
#else
  #define IMG_PRIVATE_DATA_SIZE 24
#endif

#define MAX_IMG_SIZE 16384 //2^14

class Image
{
public:
  typedef unsigned char uchar;
  
  Image();
  Image( const Image &img );
  Image( int w, int h );

  void Init( int w, int h );
  
  ~Image();
  void Clear() { Init(0,0); }

  int Width() const;
  int Height() const;
  bool IsEmpty() const { return Width()<=0 || Height()<=0; }

  void Fill( uchar r, uchar g, uchar b );
  void Rescale( int w, int h, int interpolation=1 ) { MyRescale(w,h,interpolation); }

  void DrawLine( int x1, int y1, int x2, int y2 );
  void DrawLine( int x1, int y1, int x2, int y2, int width );
  void DrawRect( int x1, int y1, int x2, int y2 );
  void DrawRectOutline( int x1, int y1, int x2, int y2 );
  void DrawPixel( int x, int y );
  void DrawPoint( int x, int y, int width );
  void DrawImage( const Image &img, int x, int y );
  void DrawText( int x, int y, const char* text, const uchar* colbg = NULL )         { MyDrawText(x,y,text,colbg); }
  void DrawText( int x, int y, const std::string &text, const uchar* colbg = NULL )  { MyDrawText(x,y,text.c_str(),colbg); }

  void SetColor( uchar r, uchar g, uchar b );
  void SetColorF( double r, double g, double b )                       { SetColor((uchar)(255.99*r), (uchar)(255.99*g), (uchar)(255.99*b)); }
  void SetOpacity( float a );

  void GetPixel( int x, int y, uchar &r, uchar &g, uchar &b) const;
  void SetPixel( int x, int y, uchar r, uchar g, uchar b )            { SetColor(r, g, b); DrawPixel(x, y); }

  bool WritePNM( const std::string &filename ) const;
  bool WritePng( const std::string &filename ) const;
  bool Write( const std::string &filename ) const;

  bool Read( const std::string &filename );

  Image& operator*=( float x );
  Image& operator+=( Image &img );
  Image& operator=( const Image &img );

  uchar* GetDataAddress( int x, int y, int c );

  template<class Array>  void DrawPixel( const Array &p )                    { DrawPixel(p[0],p[1]); }
  template<class Array>  void SetColor( const Array &v )                     { SetColor(v[0], v[1], v[2]); }
  template<class Array>  void SetColorF( const Array &v )                     { SetColorF(v[0], v[1], v[2]); }
  template<class Array, class Array2>  void SetPixel(const Array &p, const Array2 &v) { SetColor(v); DrawPixel(p); }
  template<class Array, class Array2>  void SetPixelF(const Array &p, const Array2 &v) { SetColorF(v); DrawPixel(p); }
  template<class Array, class Array2>  void GetPixel(const Array &p, Array2 &v) { GetPixel(p[0], p[1], v[0], v[1], v[2]); }
  
private:
  // hack around stupid macros in windows headers that redefine DrawText and Rescale
  void MyDrawText( int x, int y, const char* text, const uchar* colbg );
  void MyRescale( int w, int h, int interpolation );

  struct CImgInstanceHolder
  {
    uchar buffer[IMG_PRIVATE_DATA_SIZE];
  } __attribute__((__may_alias__)) priv;
  
  uchar col[3];
  float opacity;
};


void DrawImageGrid(Image& dst, const std::vector< Image >& images, int dir = 0);
void DrawImageGrid(Image &dst, const std::vector<std::vector<Image> > &images);



#endif
