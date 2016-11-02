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
#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
  Tools to convert nested dicts to .info files from boost::property_tree.
"""

import cStringIO
import numpy

__all__ = ['Vec', 'dicttoinfo']

class Vec(numpy.ndarray):
  """
    This wrapper around numpy arrays has the only purpose of making it
    easy to put a representation of my c++ VecNd class in .info files.
    
    Therefore the str conversion is overloaded so it matches the c++ side.
  """
  def __new__(cls, input_array):
    obj = numpy.asarray(input_array).view(cls)
    return obj
  
  def __str__(self):
    return "<%s>" % (','.join(str(q) for q in self))


def print_dicttoinfo(write, name, obj, indent):
  """
    writes a .info file from boost::property_tree based on nested python dicts
    writer = a callable, callable with a string argument
    arg = a dictionary, can have nested dictornaries
  """
  indentstr = " "*indent*2
  not_toplevel = indent>=0
  if isinstance(obj, dict):
    if not_toplevel: write("%s%s {\n" % (indentstr, name,))
    for k, v in obj.iteritems():
      print_dicttoinfo(write, k, v, indent+1)
    if not_toplevel: write("%s}\n" % (indentstr))
  elif type(obj) in (tuple, list, numpy.ndarray):
    if not_toplevel: write("%s%s {\n" % (indentstr, name,))
    for q in obj:
      print_dicttoinfo(write, '""', q, indent+1)
    if not_toplevel: write("%s}\n" % indentstr)  
  else:
    if isinstance(obj, (str,unicode)):
      obj = '"%s"' % obj
    elif isinstance(obj, bool):
      obj = "true" if obj else "false"
    write("%s%s %s\n" % (indentstr, name,obj)) 


def dicttoinfo(args):
  """
    convert nested dicts to a string which contains a .info file from boost::property_tree
    args = a dictionary, can have nested dictornaries
  """ 
  f = cStringIO.StringIO()
  print_dicttoinfo(f.write, "", args, -1)
  return f.getvalue()
  


if __name__ == '__main__':
  import numpy as np
  d = dict(
    a = 5,
    b = 'str',
    c = np.concatenate((np.arange(0., 1., 0.25), np.arange(1., 3., 1.), np.arange(3., 97., 3.))),
    d = Vec((1,2,3,4))
  )
  print dicttoinfo(d)