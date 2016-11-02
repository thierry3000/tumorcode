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

#include "field.h"
#define cimg_use_png 1
#define cimg_debug 0
#define cimg_OS 1  // linux
#define cimg_display 0

#include "CImg.h"
#include "drawing.h"
#include "field_impl.h"


template<class T>
Array3d<T> resizedArray3d(ConstArray3d<T> arr, const Int3 &nz, float filterwidth)
{
  if(!arr.isContiguous())
  {
    arr = Array3d<T>(arr,Cons::COPY);
  }
  // a shared image taking the address of the array data
  cimg_library::CImg<T> t(const_cast<T*>(arr.getPtr()),arr.size()[0],arr.size()[1],arr.size()[2],1,true);
  // blur
  if(filterwidth>0.)
  {
    Float3 s(0.5 * filterwidth);
    for(int i=0; i<3; ++i) 
    {
      if(nz[i]<0) 
      {
        s[i] *= -100./nz[0];
      }
      else 
      {
        s[i] *= float(arr.size()[i])/nz[i];
      }
    }
    t.assign(t.get_blur(s[0],s[1],s[2],true));
  }
  // get resized image
  cimg_library::CImg<T> a = t.get_resize(nz[0],nz[1],nz[2],-100,5,0);
  // prevent freeing of memory
  a._is_shared = true;
  Int3 newl(a.width(),a.height(),a.depth());
  return Array3d<T>(newl, Array3d<T>::Lattice::CalcStrides(newl), a.data(), newl.prod(), true);
}


template
Array3d<float> resizedArray3d<float>(ConstArray3d<float> arr, const Int3 &nz, float filterwidth);

template
Array3d<double> resizedArray3d<double>(ConstArray3d<double> arr, const Int3 &nz, float filterwidth);

template
Array3d<uchar> resizedArray3d<uchar>(ConstArray3d<uchar> arr, const Int3 &nz, float filterwidth);



//template<class T>
//inline T transferPixel(uchar val) { return val*(1./255.); }
//
//template<>
//inline uchar transferPixel<uchar>(uchar val) { return val; }


template<class T>
inline double pixelTransferFactor() { return 255.; }

template<>
inline double pixelTransferFactor<uchar>() { return 1.; }


template<class T>
Array3d<T> ReadArrayFromImage(const string &fn, int channel, int zsize)
{
  Image img;
  img.Read(fn.c_str());
  Array3d<T> res(Int3(img.Width(),img.Height(),zsize));
  FOR_REG2V2(p,res.size2())
  {
    uchar cc[3];
    img.GetPixel(p, cc);
    for(int z=0; z<zsize; ++z)
    {
      res(p.x(),p.y(),z) = cc[channel]/pixelTransferFactor<T>();
    }
  }
  return res;
}


template Array3d<float> ReadArrayFromImage<float>(const string &fn, int channel, int zsize);
template Array3d<double> ReadArrayFromImage<double>(const string &fn, int channel, int zsize);
template Array3d<uchar> ReadArrayFromImage<uchar>(const string &fn, int channel, int zsize);


template<class T>
void DrawArrayInt(Image &img, const ConstArray3d<T> &arr, bool normalized, double val_scale, double val_offset, bool bText, const string &title)
{
  if (arr.getBox().min != Int3(0))
    throw std::runtime_error("Array Box must be at (0,0,0) for drawing");
  Int3 s(arr.size());
  img.Init(s.x(),s.y());
  img.Fill(128,128,128);
  const int z = s.z()/2;

  my::MinMax<double> mm;
  if(normalized || bText)
  {
    FOR_BBOX3(p, arr.getBox())
    {
      mm.add(arr(p));
    }
  }
  
  const BBox3 bb = arr.getBox().Set(2,z);  
  FOR_BBOX3(p,bb)
  {
    double v = arr(p);
    if(normalized) v = my::normalize<double>(v,mm.min,mm.max);
    else v = val_offset + val_scale * v;
    uchar g = pixelTransferFactor<T>() * v;
    img.SetPixel(p.x(), p.y(), g, g, g);
  }
  
  if(bText || !title.empty())
  {
    int lineheight = 12, charwidth = 6;
    int lines = (bText ? 2 : 0) + (title.empty() ? 0 : 1);
    string s1, s2;
    if (bText)
    {
      s1 = boost::str(boost::format("min: %f") %  mm.min);
      s2 = boost::str(boost::format("max: %f") %  mm.max);
    }
    int chars = std::max(std::max(s1.length(), s2.length()), title.length());
    Image bigimg(std::max(img.Width(), 2+chars*charwidth), img.Height()+lines*lineheight+2);
    bigimg.Fill(128,50,50);
    bigimg.DrawImage(img, 0, lineheight*lines + 2);
    bigimg.SetColor( 255, 255, 255 );
    uchar c[3] = { 0, 0, 0 };
    int y = 2;
    if (!title.empty())
    {
      bigimg.DrawText( 2, y, title, c);
      y += lineheight;
    }
    if (bText)
    {
      bigimg.DrawText( 2, y, s1 + "\n" + s2, c);
    }
    img = bigimg;
  }
}


template void DrawArrayInt<uchar>(Image &img, const ConstArray3d<uchar> &arr, bool normalized, double val_scale, double val_offset, bool bText, const string &title);
template void DrawArrayInt<float>(Image &img, const ConstArray3d<float> &arr, bool normalized, double val_scale, double val_offset, bool bText, const string &title);
template void DrawArrayInt<double>(Image &img, const ConstArray3d<double> &arr, bool normalized, double val_scale, double val_offset, bool bText, const string &title);
template void DrawArrayInt<char>(Image &img, const ConstArray3d<char> &arr, bool normalized, double val_scale, double val_offset, bool bText, const string &title);
template void DrawArrayInt<int>(Image &img, const ConstArray3d<int> &arr, bool normalized, double val_scale, double val_offset, bool bText, const string &title);

#define INSTANTIATE(T) \
  template void PrintArray<T>(const ConstArray3d<T> a, std::ostream &os); \
  template DynArray<T> ComputeConvolutionValues(ConstArray3d<T> &stencil, ConstArray3d<bool> mask);

INSTANTIATE(float)
INSTANTIATE(double)
INSTANTIATE(int)
//INSTANTIATE(char)
INSTANTIATE(Float3)

DynArray<int> ComputeConvolutionOffsets(const BBox3 &stencil_box, const Int3 &strides, ConstArray3d<bool> mask)
{
  int n = Volume(stencil_box);
  DynArray<int> offsets(n);
  int k = 0;
  FOR_BBOX3(p, stencil_box)
  {
    if (!mask.empty() && !mask(p)) continue;
    offsets[k] = strides.dot(p);
    ++k;
  }
  offsets.resize(k);
  return offsets;
}



