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
#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import h5py as h5
import os,sys
import getopt
import h5py.h5a
import h5py.h5l
import numpy as np
import posixpath
from mystruct import Struct


def exitmsg(msg):
    print msg
    sys.exit(0)

def H5PrintOpt():
  return Struct(
    recursive = False,
    print_full_dataset = 0,
    print_attributes = False,
    print_eval = None
  )


class Printer:
    def __init__(self, h5opt = H5PrintOpt()):
        self.h5opt = h5opt
        self.indent = 0
        self.maxlevel = 1000000 if h5opt.recursive else 1
    def writeln(self, s):
        print ('  '*self.indent)+s
    def begin(self):
        self.indent += 1
    def end(self):
        self.indent -= 1
    def writeattrs(self, o, level):
      #for k in o.attrs.iterkeys():
        #a = h5py.h5a.open(o.attrs.get_id(), k)
        #print k, a.shape, a.get_type(),
        #if isinstance(a.get_type(), h5py.h5t.TypeArrayID):
          #print a.get_type().get_array_dims()
        #else:
          #print ''
      for k, v in o.attrs.iteritems():
          self.writeln('| %s = %s' % (k, str(v)))
    def writeds(self, g, level):
      s = '%s (%s) (%i A) %s %s' % (
        posixpath.basename(g.name),
        g.__class__.__name__,
        len(g.attrs),
        g.shape,
        g.dtype)
      self.writeln(s)
      if self.h5opt.print_full_dataset == 1:
        print g[...]
      elif self.h5opt.print_full_dataset == 2:
        # read it
        arr = np.zeros(g.shape, dtype = g.dtype)
        g.read_direct(arr)
        for index, x in np.ndenumerate(arr):
          print '%s = %s' % (index, x)
      elif self.h5opt.print_eval:
        data = np.zeros(g.shape, dtype = g.dtype)
        g.read_direct(data)
        globals_ = dict(np.__dict__)
        globals_.update(data = data)
        print eval(self.h5opt.print_eval, globals_)
    
    def writeobject(self,g,level=0):
      if isinstance(g,h5.Dataset):
        self.writeds(g, level)
      else:
        self.writeln('%s (%s) (%i A)' % (posixpath.basename(g.name), g.__class__.__name__, len(g.attrs)))
      self.begin()
      if self.h5opt.print_attributes:
        self.writeattrs(g, level)
      if isinstance(g,h5.Group) and level<self.maxlevel:
          for oname in g.keys():
              if g.id.links.get_info(oname).type in (h5py.h5l.TYPE_SOFT, h5py.h5l.TYPE_EXTERNAL):
                self.writeln('%s (Link) %s' % (oname, str(g.id.links.get_val(oname))))
              else:
                self.writeobject(g[oname],level+1)
      self.end()

if __name__ == '__main__':
  usage = """
  %s [-r | -f | -a | --fullds | -e expr] hdf5file [path-in-hdf5-file]
  """
  try:
    opt, args = getopt.getopt(sys.argv[1:],'rfae:',['help','fullds'])
  except getopt.GetoptError, e:
    exitmsg(usage % sys.argv[0])
  help = False
  h5opt = H5PrintOpt()
  for o, a in opt:
    if o == '-r':
      h5opt.recursive = True
    elif o == '-f':
      h5opt.print_full_dataset = max(1, h5opt.print_full_dataset)
    elif o == '--fullds':
      h5opt.print_full_dataset = 2
    elif o == '-a':
      h5opt.print_attributes = True
    elif o == '-e':
      h5opt.print_eval = a
    elif o == '--help':
      help = True
  if help:
    exitmsg(usage % sys.argv[0])
  fn = args[0]
  if not os.path.isfile(fn):
      exitmsg('%s must exist' % fn)
  if len(args)>1:
      path = args[1]
  else:
      path = '/'
  try:
      f = h5.File(fn,'r')
  except:
      exitmsg('%s must be a valid hdf5 file' % fn)
  try:
      g = f[path]
  except:
      exitmsg('%s must be a valid path in the hdf file' % path)
  p  = Printer(h5opt)
  p.writeobject(g)
