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
if __name__ == '__main__':
  import os.path, sys
  sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'..'))

import h5py as h5
import pprint
import vtk

import myutils
import extractVtkFields
import hdfvessels2vtk


def writeFields_(f, search_groups, outfn):
  e = extractVtkFields.Extractor(f, search_groups, recursive = True) 
  print 'found field datasets:'
  pprint.pprint(e.getDatasetPaths())
  e.write(outfn)
  del e

def writeVessels_(f, search_groups, outfn):
  vessels = hdfvessels2vtk.ConvertMyHdfVesselsToVTKPolydata(f)
  writer = vtk.vtkPolyDataWriter()
  writer.SetInput(vessels)
  writer.SetFileName(outfn)
  writer.Write()
  del vessels


def hdftumor2vtk(fn, search_groups=['.'], outfn=None, dstdir='', writeVessels=True, writeFields=True):
  f = h5.File(fn, 'r')  
  a = f['/'].attrs
  if outfn:
    outfn = "%s-%%s.vtk" % (outfn)
  elif ('OUTPUT_NUM' in a and 'OUTPUT_NAME' in a):
    outfn = "%s-%%s-%04i.vtk" % (outfn, a['OUTPUT_NUM'])
  else:
    outfn = "%s-%%s.vtk" % (os.path.splitext(os.path.basename(fn))[0])
  outfn = os.path.join(dstdir, outfn)
  if writeFields:
    writeFields_(f, search_groups, outfn % 'fields')
  if writeVessels:
    writeVessels_(f, search_groups, outfn % 'vessels')


if __name__ == '__main__':
  fn, search_groups = myutils.splitH5PathsFromFilename(sys.argv[1])

    

