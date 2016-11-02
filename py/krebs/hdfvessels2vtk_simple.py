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
#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os, sys
from os.path import join, basename, dirname, splitext, normpath, abspath
sys.path.append(join(dirname(__file__),'..'))
import krebsutils
import numpy as np
import h5py
from vtk import *
from vtkcommon import *
import extensions
import itertools
import myutils
import fnmatch


def ConvertMyHdfVesselsToVTKPolydata(h5file, grpname):
  import posixpath
  #edges, pos, rad, flags = krebsutils.read_vessels_from_hdf(h5file[grpname]['vessels'], ('position', 'radius','flags'))
  edges, pos, rad, flags, hema, flow= krebsutils.read_vessels_from_hdf(h5file, ('position',
                                                                     'radius',
                                                                     'flags',
                                                                     'hematocrit',
                                                                     'flow'))
  
  num_verts = len(pos)
  num_edges = len(edges)
  print "v %i e %i" % (num_verts, num_edges)

  pts = vtkPoints()
  pts.SetNumberOfPoints(len(pos))
  for i,(x,y,z) in enumerate(pos):
    pts.SetPoint(i, (x,y,z))

  cells = asVtkCellArray(edges)

  polydata = vtkPolyData()
  polydata.SetPoints(pts)
  polydata.SetLines(cells)

  if 1:
    #radius
    polydata.GetCellData().AddArray(asVtkArray(rad, "radius", vtkFloatArray))
    #nodes
    node_rad = np.zeros(num_verts)
    node_n   = np.zeros(num_verts)
    for r,(a,b) in itertools.izip(rad,edges): # itertools, or it will blow up in yr face cuz ultra slow python lists and shit will be created
      node_rad[a] += r
      node_rad[b] += r
      node_n[a] += 1.
      node_n[b] += 1.
    node_rad /= np.maximum(1, node_n)
    polydata.GetPointData().AddArray(asVtkArray(node_rad, "point_radius", vtkFloatArray))

  if 1:
    #flags
    iscirc = np.asarray(np.bitwise_and(np.uint64(flags),np.uint64(krebsutils.CIRCULATED)) > 0, dtype=np.int32)
    polydata.GetCellData().AddArray(asVtkArray(iscirc, "isCirculated", vtkIntArray))
    
  if 1:
    #flags
    isartery = np.asarray(np.bitwise_and(np.uint64(flags),np.uint64(krebsutils.ARTERY)) > 0, dtype=np.int32)
    polydata.GetCellData().AddArray(asVtkArray(isartery, "isArtery", vtkIntArray))
    
  if 1:
    #flags
    isvein = np.asarray(np.bitwise_and(np.uint64(flags),np.uint64(krebsutils.VEIN)) > 0, dtype=np.int32)
    polydata.GetCellData().AddArray(asVtkArray(isvein, "isVein", vtkIntArray))

  if 1:
    #flags
    iscap = np.asarray(np.bitwise_and(np.uint64(flags),np.uint64(krebsutils.CAPILLARY)) > 0, dtype=np.int32)
    polydata.GetCellData().AddArray(asVtkArray(iscap, "isCapillary", vtkIntArray))
    
  if 1:
    #flags
    isboundary = np.asarray(np.bitwise_and(np.uint64(flags),np.uint64(krebsutils.BOUNDARY)) > 0, dtype=np.int32)
    polydata.GetCellData().AddArray(asVtkArray(isboundary, "isBoundary", vtkIntArray))
    
  if 1:
    #edge data
    #indices of edges
    print len(edges)
    polydata.GetCellData().AddArray(asVtkArray(np.arange(0,len(edges)), "IndexOfEdges2", vtkIntArray))
    polydata.GetCellData().AddArray(asVtkArray(edges, "IndexOfEdges", vtkIntArray))
    
    polydata.GetCellData().AddArray(asVtkArray(hema, "hematocrit", vtkFloatArray))
    polydata.GetCellData().AddArray(asVtkArray(flow, "flow", vtkFloatArray))

  if 0:
    #vesselgroup = h5file[posixpath.join(grpname, 'vessels')]
    vesselgroup = h5file['/']

    #everything else
    ignore = ['lattice_pos','roots','node_a_index','node_b_index','radius']
    g = vesselgroup['edges'].values() + vesselgroup['nodes'].values()
    for item in g:
      short_name = item.name.split('/')[-1]
      if short_name in ignore: continue
      n = item.shape
      if len(n) <> 1:
        print 'ignoring item %s because if unsupported shape %s' % (item.name, n)
        continue
      n = n[0]
      if n == num_verts:
        print 'vert data %s [%s,%s]' % (item.name, item.dtype, str(item.shape))
        polydata.GetPointData().AddArray(asVtkArray(np.asarray(item), short_name))
      elif n == num_edges:
        print 'edge data %s [%s,%s]' % (item.name, item.dtype, str(item.shape))
        polydata.GetCellData().AddArray(asVtkArray(np.asarray(item), short_name))
      else:
        print 'ignoring item %s because array size does not match' % (item.name)
  return polydata




def removeUncirculatedVessels(dataset):
  flagarr = fromVtkArray(dataset.GetCellData().GetArray("flags"))
  ids = [ i for i, f in enumerate(flagarr) if (f & krebs.CIRCULATED) ]
  return copyPolyDataCells(dataset, ids)



#def walkh5(g, pattern):
#  res = []
#  if not fnmatch.fnmatch(g.name, pattern):
#    if isinstance(g, h5py.Group):
#      for gchild in g.itervalues():
#        res += walkh5(gchild, pattern)
#  else:
#    res.append(g.name)
#  return res


if __name__ == '__main__':
  from os.path import basename, dirname, join, splitext

  fn = sys.argv[1]
  fn, _ = myutils.splitH5PathsFromFilename(fn)
  grpnames = [ sys.argv[2] ]
  f = h5py.File(fn, 'r')
  if grpnames[0]:
    dirs = myutils.walkh5(f['/'], grpnames[0])
  else:
    dirs = ['']
  print 'found items %s to search pattern %s' % (str(dirs), grpnames[0])
  for d in dirs:
    polydata = ConvertMyHdfVesselsToVTKPolydata(f[join('/',d)], '.');
    writer = vtkPolyDataWriter()
    writer.SetInputData(polydata)
    writer.SetFileName("%s-vessels-%s.vtk" % (splitext(basename(fn))[0], d.replace('/','-')))
    writer.Write()

  #args = sys.argv[1:]
  #fn, (grpname, ) = myutils.splitH5PathsFromFilename(args[0])
  #f = h5.File(fn,'r')
  #a = f.attrs

  #out = basename(fn);
  #if out.endswith('.h5'):
    #out = splitext(out)[0]
  #if grpname and grpname <> '.':
    #tmp = grpname.replace('/','-')
    #if not tmp.startswith('-'):
      #tmp = '-'+tmp
    #out += tmp
  #out += '-vessels.vtk'

  #print grpname

  #polydata = ConvertMyHdfVesselsToVTKPolydata(f, grpname)
  #writer = vtkPolyDataWriter()
  #writer.SetInput(polydata)
  #writer.SetFileName(out)
  #writer.Write()
