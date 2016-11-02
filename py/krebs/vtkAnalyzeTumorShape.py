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
# -*- coding: utf-8 -*-

import os, sys
import numpy as np
import h5py
import math
import pprint
import vtk
import vtkcommon
import posixpath
from mystruct import Struct


def analyzeTumorShape(dataset,
                      do_curvature = True,
                      do_integrals = True,
                      do_linecomps = False,
                      contourspec = ('S', 0.)):
  out = Struct()
  #generate tumorsurface
  contour_name, contour_value = contourspec
  assert contour_value in (0., 0.5)
  dataset.GetPointData().SetActiveScalars(contour_name)
  iso = vtkcommon.vtkContour(dataset, contour_value)
  # measurement
  if do_linecomps:
    lc, pts = vtkcommon.vtkGetLineComponents(iso)
    out['lines'] = lc
    out['points'] = pts
  if do_curvature:
    # curvature
    tmp = vtkcommon.vtkCurvature(iso, 9)
    out['curvature'] = np.asarray(vtkcommon.fromVtkArray(tmp.GetPointData().GetArray("curvature")), dtype=np.float32)
    del tmp
  if do_integrals:
    # integrate volume (area, length)
    _, _, out['iso_area'] = vtkcommon.vtkIntegrateData(iso)
    # integrated data
    cd, pd, volume = vtkcommon.vtkIntegrateData(dataset)
    cd = dict((dataset.GetCellData().GetArray(i).GetName(), cd[i]) for i in xrange(dataset.GetCellData().GetNumberOfArrays()))
    pd = dict((dataset.GetPointData().GetArray(i).GetName(), pd[i]) for i in xrange(dataset.GetPointData().GetNumberOfArrays()))
    data = cd.copy()
    data.update(pd)
    out['tumor_volume'] = (data[contour_name] + volume)*0.5 if contour_value==0 else data[contour_name]
    out['radius_estimate'] = math.sqrt(out.tumor_volume/math.pi)
    out['area_estimate'] = 2.*math.pi*out.radius_estimate
  return out


def writeShapeDataToH5(m, dataset, data):
    curv = data.pop('curvature', None)
    m.attrs['version'] = 1
    for k, v in data.iteritems():
      m.attrs[k] = v
    if curv is not None:
      #m.create_dataset('curvatures', data = curv, compression = 'gzip', compression_opts = 9)
      if 'curvatures' in m:
        del m['curvatures']
      ds = m.require_dataset('curvatures', curv.shape, curv.dtype, compression = 'gzip', compression_opts = 9)
      ds.write_direct(curv)

  
def analyzeDealiiTumorShapeAndSaveToH5(filenames, dt, dstpath):
  #re dt: dt = 2 for older dealii simulations, dt = 100 for the ones where the tumor can die
  filenames.sort()
  for i, filename in enumerate(filenames):
    print "processing %i: %s" % (i, filename)
    dataset, outfilename = vtkcommon.ZippedRead(vtk.vtkXMLUnstructuredGridReader, filename)
    # create file and store data    
    i = int(outfilename[-4:])
    f = h5py.File(os.path.join(dstpath, outfilename+".h5"), 'w')
    g = f['/']
    g.attrs['time'] = i * float(dt)
    g.attrs['out_num'] = i
    m = g.create_group("tumor")
    ds = m.create_dataset("vtkfile", data = np.fromstring(vtkcommon.ZippedOpen(filename).read(), dtype=np.uint8), compression = 'gzip', compression_opts = 9)
    ds.attrs["TYPE"] = "VTU_FILE"
    # do the measurement and save
    data = analyzeTumorShape(dataset)
    writeShapeDataToH5(g.require_group("measure_tumglob"), dataset, data)
    print "-> %s" % (f.filename)
    f.close()


def analyzeBulkTissueTumorShapeAndSaveToH5(filename, dstfilename):
  from hdffields2vtk import ConvertMyHdfFieldsToVtkImage
  f = h5py.File(filename, "r")
  fdst = h5py.File(dstfilename, "a")
  # read write mode, get the tumor data, analyze and save
  groups = [ g for k, g in f['/'].iteritems() if k.startswith('out') ]
  groups.sort(key = lambda g: g.attrs['time'])

  for g in groups:
      print 't = %f' % g.attrs['time']
      gdst = fdst.require_group(posixpath.join(
                                #'measurements', 
                                posixpath.basename(g.name), 
                                'measure_tumglob'))
      dataset = ConvertMyHdfFieldsToVtkImage(g, ['ptc'], f['field_ld'])
      dataset = vtkcommon.vtkCellDataToPointData(dataset)
      data = analyzeTumorShape(dataset,
                               contourspec = ('ptc', 0.5),
                               #do_curvature = (g is groups[-1])
                               )
      pprint.pprint(data)
      writeShapeDataToH5(gdst, dataset, data)


if __name__ == '__main__':
  analyzeBulkTissueTumorShapeAndSaveToH5(sys.argv[1], sys.argv[2])