#!/usr/bin/env python2
# -*- coding: utf-8 -*-
'''
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
'''

"""
  Here i collect some utilities. Inclusion of this module actually modifies the content of other
  modules such as numpy. For example it replaces numpys asarray, ma.array, asanyarray by
  function which can read data directly from hdf5 datasets. This is its main use. Other
  stuff is rarely used if at all.
"""


def check_and_set_augmentation_need(module_name):
  try:
    module = __import__(module_name)
  except ImportError:
    return None, False
  try:
    getattr(module, '__mw_extension_hook_installed__')
  except:
    setattr(module, '__mw_extension_hook_installed__', True)
    return module, True
  else:
    return module, False

math, _augment_math = check_and_set_augmentation_need('math')
np  , _augment_np   = check_and_set_augmentation_need('np')
gsl , _augment_gsl  = check_and_set_augmentation_need('pygsl')
pil , _augment_pil  = check_and_set_augmentation_need('PIL')
h5py, _augment_h5py = check_and_set_augmentation_need('h5py')

if _augment_math:
  def sgn(x):
      return 1. if x>=0 else -1.

  def cut(x,a,b):
      return a if x<a else (b if x>b else x)

  def nextPowerOf2(x):
      if x==0:
          r = 0.
      else:
          e = math.log(x,2)
          e = math.ceil(e)
          r = math.pow(2.,e)
      return int(r) if isinstance(x,int) else r

  math.sgn = sgn
  math.abs = abs
  math.cut = cut
  math.nextPowerOf2 = nextPowerOf2

if _augment_np:
    def resize2(a,t):
        s = a.shape
        if s==t: return a
        b = np.zeros(t,dtype=a.dtype)
        q = tuple( slice(min(ls,lt)) for ls,lt in zip(s,t) )
        b[q] = a
        return b
    np.resize2 = resize2

    def derivative(x,y):
        n, = y.shape
        assert (n,)==x.shape, "x and y must have equal shape and be one dimensional"
        def bounds1(v):
            vm2 = v[n-2:n-1]
            vm1 = v[n-1:n]
            vp1 = v[1:2]
            vp0 = v[0:1]
            v = np.append(v,2.*vm1-vm2,axis=0)
            v = np.append(2.*vp0-vp1,v,axis=0)
            return v
        y = bounds1(y)
        x = bounds1(x)
        r = (y[2:2+n]-y[0:n])/(x[2:2+n]-x[0:n])
        return r
    np.derivative = derivative

    def vstackDot(a,b):
        return np.sum(a*b,axis=1)
    np.vstackDot = vstackDot

    def vstackNorm(a):
        return np.sqrt(vstackDot(a,a))
    np.vstackNorm = vstackNorm

    def vstackScalarMult(s,a):
        assert len(a.shape)==2
        assert len(s.shape)==1
        assert s.shape[0] == a.shape[0]
        return np.column_stack([s]*a.shape[1])*a
    np.vstackScalarMult = vstackScalarMult

    def vstackComponentMult(v,a):
        v = np.asarray(v)
        assert len(v.shape)==1
        assert len(a.shape)==2
        assert a.shape[1] == v.shape[0]
        return np.tile(v,(a.shape[0],1))*a
    np.vstackComponentMult = vstackComponentMult

if _augment_np and _augment_gsl:
    def linearSimple(x,y,w=None):
        if len(x)<1:
            raise ArgumentError('x-array must have at least two entries')
        if len(x)<>len(y)<>len(w):
            raise ArgumentError('input arrays must have equal size')
        if w is None:
            _t = pygsl.fit.linear(x,y)
        else:
            _t = pygsl.fit.wlinear(x,w,y)
        c, m, c00, c01, c11, chisqr = _t
        try:
            c11 = math.sqrt(c11)
        except:
            c11 = 0.
        try:
            c00 = math.sqrt(c00)
        except:
            c00 = 0.
        return m,c,c11,c00

    import pygsl.fit
    pygsl.fit.linearSimple = linearSimple


if _augment_np:
  def linearFit(x,y):
    import numpy as np
    beta, alpha = np.polyfit(x,y,1,None,False)
    n = len(x)
    eps = y - alpha - beta * x
    a = np.sum(np.power(eps,2.))/(n-2)
    b = np.sum(np.power(x - np.average(x),2))
    sbeta = np.sqrt(a/b)
    return beta, alpha, sbeta, 0.
  np.linearFit = linearFit



if _augment_np and _augment_h5py:
    np.h5py = h5py
    np._orig_as_array = np.asarray
    np._orig_ma_array = np.ma.array
    np._orig_asanyarray = np.asanyarray

    def asarray_(ds, dtype=None, order=None):
      import h5py
      import numpy as np
      if isinstance(ds, h5py.Dataset):
  #      a = np.empty(ds.shape, ds.dtype, 'C')
  #      ds.read_direct(a)
  #      ds = a
        ds = np.array(ds)
      return np._orig_as_array(ds, dtype, order)

    def ma_array_(ds, **kwargs):
      import h5py
      import numpy as np
      if isinstance(ds, h5py.Dataset) and 'CLASS' in ds.attrs and ds.attrs['CLASS']=='MASKED_ARRAY':
        ds = np._orig_ma_array(asarray_(ds['_data']), mask=asarray_(ds['_mask']))
      return np._orig_ma_array(ds, **kwargs)

    def asanyarray_(ds, dtype=None, order=None):
      import h5py
      import numpy as np
      if isinstance(ds, h5py.Dataset):
        if ds.attrs.get('CLASS','')=='MASKED_ARRAY':
          ds = np._orig_ma_array(asarray_(ds['_data']), mask=asarray_(ds['_mask']))
        else:
          ds = asarray_(ds)
      return np._orig_asanyarray(ds, dtype, order)

    np.asarray = asarray_
    np.ma.array = ma_array_
    np.asanyarray = asanyarray_





if _augment_h5py:
#  def recreate_group_(f, name):
#    try:
#      del f[name]
#    except KeyError:
#      pass
#    return f.create_group(name)

#  h5py.highlevel.Group.recreate_group = recreate_group_
#  h5py.highlevel.File.recreate_group = recreate_group_

  #h5py.highlevel.File.__enter__ = lambda self: self
  #h5py.highlevel.File.__exit__()

  def recreate_(func):
    def f(self, *args, **kwargs):
      name = args[0]
      if name in self:
        del self[name]
      return func(self, *args, **kwargs)
    return f

  h5py.highlevel.Group.recreate_group = recreate_(h5py.highlevel.Group.create_group)
  h5py.highlevel.File.recreate_group = recreate_(h5py.highlevel.File.create_group)
  h5py.highlevel.Group.recreate_dataset = recreate_(h5py.highlevel.Group.create_dataset)
  h5py.highlevel.File.recreate_dataset = recreate_(h5py.highlevel.File.create_dataset)

  def try_del(self, *names):
    for name in names:
      try:
        del self[name]
      except KeyError:
        pass
  h5py.highlevel.Group.try_del = try_del
  h5py.highlevel.File.try_del = try_del
  #h5py.highlevel.Attributes.try_del = try_del

  def new_get_item(func):
    def f(self, name):
      try:
        return func(self, name)
      except KeyError, e:
        fn = self.file.filename
        path = self.name
        raise KeyError(r'"%s:%s/%s": %s' % (fn, path, name, str(e)))
    return f

  h5py.highlevel.Group.__getitem__ = new_get_item(h5py.highlevel.Group.__getitem__)
  h5py.highlevel.File.__getitem__ = new_get_item(h5py.highlevel.File.__getitem__)


  def new_exception_handling(func):
    def f(self, filename, *args, **kwargs):
      try:
        return func(self, filename, *args, **kwargs)
      except Exception, e:
        raise type(e)(r'%s in %s' % (str(e),filename))
    return f

  h5py.highlevel.File.__init__ = new_exception_handling(h5py.highlevel.File.__init__)


  def test_h5py():
    with h5py.File('memfile2', 'w', 'core', backing_store = False) as f:
      f.create_dataset('fubar', data = [1,2,3])
      f.recreate_dataset('fubar', data = [4, 5, 6])
      print np.asarray(f['fubar'])
      f.try_del('fubar')
      f.try_del('bar')
      print f.keys()



if _augment_pil and _augment_np:
    import PIL.Image as Image

    def image2array(im):
        import numpy as np
        if im.mode not in ("L", "F","RGB"):
            raise ValueError, "can only convert L,F and RGB images"
        if im.mode == "RGB":
            r,g,b = tuple(np.fromstring(q.tostring(), np.ubyte) for q in im.split())
            a = np.concatenate((r,g,b),axis=2)
            a.shape = im.size[1], im.size[0], 3
        elif im.mode == "L":
            a = np.fromstring(im.tostring(), np.ubyte)
            a.shape = im.size[1], im.size[0]
        else:
            a = np.fromstring(im.tostring(), np.float32)
            a.shape = im.size[1], im.size[0]
        return a
    Image.toArray = image2array

    def array2image(a):
        import numpy as np
        if len(a.shape)==2:
            h,w = a.shape
            if a.dtype == np.ubyte:
                mode = "L"
            elif a.dtype == np.float32:
                mode = "F"
            else:
                raise ValueError, "unsupported image mode"
        elif len(a.shape)==3 and a.shape[2]==3:
            h,w,c = a.shape
            if a.dtype == np.ubyte:
                mode = "RGB"
            else:
                raise ValueError, "unsupported image mode"
        else:
            raise ValueError,"unsupported shape"
        return Image.fromstring(mode, (w,h), a.tostring())
    Image.fromArray = array2image
