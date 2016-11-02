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
import sys, os
from math import *

hasNumpy = True
try:
    import numpy as np
except ImportError:
    hasNumpy = False

class File:
  def __init__(self,f,*items):
    self.file = f
    self.__indent = 0
    self.write(*items)
  def include(self,name):
    self.writeln( '#include "%s"'%name )
    self.writeln()
  def indent(self):
    self.__indent += 1
  def dedent(self):
    self.__indent -= 1
    assert self.__indent >= 0
  def block_begin(self):
    self.writeln( "{" )
    self.indent()
  def block_end(self):
    self.dedent()
    self.writeln( "}" )
    if self.__indent == 0:
      # blank line if this is a top level end
      self.writeln( )
  def write(self,*items):
    for item in items:
      if isinstance(item,str):
        lines = item.splitlines()
        for l in lines:
            self.writeln(l)
      else:
        item.write(self)
  def writeln(self,s=""):
    #print "  "*self.__indent+s
    self.file.write("  "*self.__indent+s+'\n')

class Vector:
  def __init__(self,*args):
    if len(args) == 1:
      self.v = args[0]
    else:
      self.v = args
  def __str__(self):
    return "<%s>"%(", ".join([str(x)for x in self.v]))
  def __repr__(self):
    return "Vector(%s)"%self.v
  def __mul__(self,other):
    return Vector( [r*other for r in self.v] )
  def __rmul__(self,other):
    return Vector( [r*other for r in self.v] )

def _isFloat(x):
    try:
        float(x)
        return True
    except:
        return False

def _isVector(q):
    return hasattr(q,'__iter__') and len(q) in (3,4) and all(_isFloat(c) for c in q)

def _toStr(v):
    return str(Vector(*v)) if _isVector(v) else str(v)

class ItemBase:
    def write(self, file):
        raise RuntimeError("write not implemented")

class Item(ItemBase):
  def __init__(self,name,args=[],opts=[],**kwargs):
    self.name = name
    args=list(args)
    for i in range(len(args)):
      if type(args[i]) == tuple or type(args[i]) == list:
        args[i] = Vector(args[i])
    self.args = args
    self.opts = opts
    self.kwargs=kwargs
  def append(self, item):
    self.opts.append( item )
  def write(self, file):
    file.writeln( self.name )
    file.block_begin()
    if self.args:
      file.writeln( ", ".join([_toStr(arg) for arg in self.args]) )
    for opt in self.opts:
      if isinstance(opt,ItemBase): #hasattr(opt,"write"):
        opt.write(file)
      elif hasattr(opt,'__iter__'):
        file.writeln( "%s %s"%(opt[0],', '.join(_toStr(v) for v in opt[1:])) )
      elif opt is not None:
        file.writeln( str(opt) )
    for key,val in self.kwargs.items():
      if _isVector(val):
        val = Vector(*val)
        file.writeln( "%s %s"%(key,val) )
      elif hasattr(val,'__iter__'):
        val = ', '.join(_toStr(v) for v in val)
        file.writeln( "%s %s"%(key,val) )
      else:
        file.writeln( "%s %s"%(key,val) )
    file.block_end()
  def __setattr__(self,name,val):
    self.__dict__[name]=val
    if name not in ["kwargs","args","opts","name"]:
      self.__dict__["kwargs"][name]=val
  def __setitem__(self,i,val):
    if i < len(self.args):
      self.args[i] = val
    else:
      i += len(args)
      if i < len(self.opts):
        self.opts[i] = val
  def __getitem__(self,i,val):
    if i < len(self.args):
      return self.args[i]
    else:
      i += len(args)
      if i < len(self.opts):
        return self.opts[i]

class Interior(Item):
  def __init__(self,*opts,**kwargs):
    Item.__init__(self,"interior",(),opts,**kwargs)

class Material(Item):
  def __init__(self,*opts,**kwargs):
    Item.__init__(self,"material",(),opts,**kwargs)

class Texture(Item):
  def __init__(self,*opts,**kwargs):
    Item.__init__(self,"texture",(),opts,**kwargs)

class Pigment(Item):
  def __init__(self,*opts,**kwargs):
    Item.__init__(self,"pigment",(),opts,**kwargs)

class Finish(Item):
  def __init__(self,*opts,**kwargs):
    Item.__init__(self,"finish",(),opts,**kwargs)

class Normal(Item):
  def __init__(self,*opts,**kwargs):
    Item.__init__(self,"normal",(),opts,**kwargs)

class Camera(Item):
  def __init__(self,*opts,**kwargs):
    Item.__init__(self,"camera",(),opts,**kwargs)

class LightSource(Item):
  def __init__(self,v,*opts,**kwargs):
    Item.__init__(self,"light_source",(Vector(v),),opts,**kwargs)

class Background(Item):
  def __init__(self,*opts,**kwargs):
    Item.__init__(self,"background",(),opts,**kwargs)

class Box(Item):
  def __init__(self,v1,v2,*opts,**kwargs):
    #self.v1 = Vector(v1)
    #self.v2 = Vector(v2)
    Item.__init__(self,"box",(v1,v2),opts,**kwargs)

class Cylinder(Item):
  def __init__(self,v1,v2,r,*opts,**kwargs):
    " opts: open "
    Item.__init__(self,"cylinder",(v1,v2,r),opts,**kwargs)

class Plane(Item):
  def __init__(self,v,r,*opts,**kwargs):
    Item.__init__(self,"plane",(v,r),opts,**kwargs)

class Torus(Item):
  def __init__(self,r1,r2,*opts,**kwargs):
    Item.__init__(self,"torus",(r1,r2),opts,**kwargs)

class Cone(Item):
  def __init__(self,v1,r1,v2,r2,*opts,**kwargs):
    " opts: open "
    Item.__init__(self,"cone", (v1,r1,v2,r2),opts,**kwargs)

class Polygon(Item):
    def __init__(self,pointlist,*opts,**kwargs):
        q = []
        q.append(str(len(pointlist)))
        for p in pointlist: q.append(str(Vector(p)))
        Item.__init__(self,"polygon",(','.join(q),),opts,**kwargs)

class Sphere(Item):
  def __init__(self,v,r,*opts,**kwargs):
    Item.__init__(self,"sphere",(v,r),opts,**kwargs)

class Union(Item):
  def __init__(self,*opts,**kwargs):
    Item.__init__(self,"union",(),opts,**kwargs)

class Intersection(Item):
  def __init__(self,*opts,**kwargs):
    Item.__init__(self,"intersection",(),opts,**kwargs)

class Difference(Item):
  def __init__(self,*opts,**kwargs):
    Item.__init__(self,"difference",(),opts,**kwargs)

class Merge(Item):
  def __init__(self,*opts,**kwargs):
    Item.__init__(self,"merge",(),opts,**kwargs)

class ClippedBy(Item):
    def __init__(self,*opts,**kwargs):
        Item.__init__(self,"clipped_by",(),opts,**kwargs)

class ContainedBy(Item):
    def __init__(self,*opts,**kwargs):
        Item.__init__(self,"contained_by",(),opts,**kwargs)

class Function(Item):
    def __init__(self,*opts,**kwargs):
        Item.__init__(self,"function",(),opts,**kwargs)

class Isosurface(Item):
    def __init__(self,*opts,**kwargs):
        Item.__init__(self,"isosurface",(),opts,**kwargs)

class Pattern(Item):
    def __init__(self,*opts,**kwargs):
        Item.__init__(self,"pattern",(),opts,**kwargs)

class Object(Item):
    def __init__(self,*opts,**kwargs):
        Item.__init__(self,"object",(),opts,**kwargs)

class Fog(Item):
    def __init__(self,*opts,**kwargs):
        Item.__init__(self,"fog",(),opts,**kwargs)

class Media(Item):
    def __init__(self,*opts,**kwargs):
        Item.__init__(self,"media",(),opts,**kwargs)

class Include(ItemBase):
    def __init__(self,filename):
        self.filename = filename
    def write(self, file):
        file.writeln( '#include "%s"' % (self.filename) )

class Colormap(ItemBase):
    def __init__(self, *opts):
        self.map = []
        #err = ArgumentError("colormap arguments must be of the form (value1, color1, value2, color2, ...")
        if not len(opts)%2==0:
            raise ValueError("colormap arguments must be of the form (value1, color1, value2, color2, ...")
        for i in xrange(0,len(opts),2):
            val = opts[i]
            col = opts[i+1]
            if not (_isFloat(val) and _isVector(col)):
                raise ValueError("colormap arguments must be of the form (value1, color1, value2, color2, ...")
            self.map.append((val, Vector(*col)))
    def write(self, file):
        file.writeln("color_map")
        file.block_begin()
        for val, col in self.map:
            file.writeln("[%s color %s]" % (str(val),str(col)))
        file.block_end()

class Declare(ItemBase):
    def __init__(self,name,item):
        self.item = item
        self.name = name
    def write(self, file):
        file.writeln( '#declare %s =' % self.name )
        self.item.write(file)

class String(ItemBase):
    def __init__(self,thestring):
        self.thestring = thestring
    def write(self, file):
        file.writeln(self.thestring)


if hasNumpy:
    max_uint32 = np.uint32(1)<<32
    max_uint32_halve = np.uint32(1)<<31

    def f_(dtype):
      x = np.array(1., dtype = dtype)
      while True:
        if np.uint32(x * max_uint32) > max_uint32_halve: return x
        else: x *= 0.9999999
    max_float_which_converts_to_less_or_equal_than_max_uint32 = f_


    def writeArrayAsDensityFile(filename, array, value_bounds):
        """ write an array as povray density file, always normalized to the range min..max """
        import struct
        array = np.atleast_3d(array)
        x,y,z = array.shape
        a, b = value_bounds

        array = np.clip(((array-a)/(b-a))*max_uint32, 0, max_float_which_converts_to_less_or_equal_than_max_uint32)
        array = np.asarray(array, dtype=np.uint32)

#        import matplotlib.pyplot as pyplot
#        pyplot.imshow(array[:,:,0])
#        pyplot.show()

        array = array.swapaxes(0,2)
        array.byteswap(True)
        file = open(filename,'wb')
        file.write(struct.pack('>hhh',x,y,z))
        array.tofile(file)
        file.close()


