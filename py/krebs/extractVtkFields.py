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
if __name__ == '__main__':
  import os.path, sys
  sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'..'))

import sys, os
import posixpath
import numpy as np
import h5py
import collections
import extensions # for asarray with h5py support
import identifycluster
if identifycluster.getname()=='snowden':
  import pyvtk as vtk
else:
  import vtk
import vtkcommon
import myutils
import krebsutils as ku


class VtkFiles(object):
  pass


class Extractor(object):
  """
    Extract a vtk dataset from a h5file
    
    Arguments:
      obj -- hdf5 file or group
      paths -- list of strings representing the paths to the datasets
      lattice_paths -- path to a group which describes the lattice on which the data lives (optional)
    
    paths must be either
      - one item, pointing to a char dataset, which holds the content of a vtk file
      - several paths, pointing to nd-float datasets, which all have a 'LATTICE_PATH' attribute
        if lattice_path is not given. Currently my LatticeDataQuad3d type is the only supported
        lattice.
  """

#  gridtype2vtureader = {
#    'vtkImageData' : vtk.vtkXMLImageDataReader,
#    'vtkPolyData' : vtk.vtkXMLPolyDataReader,
#    'vtkUnstructuredGrid' : vtk.vtkXMLUnstructuredGridReader,
#  }
  
  def __init__(self, obj, paths = None, lattice_path = None, recursive = False):
    assert isinstance(obj, (h5py.File, h5py.Group, h5py.Dataset))
    assert isinstance(obj, h5py.Dataset) or isinstance(paths, (list,tuple))
    self.file = f = obj if isinstance(obj, h5py.File) else obj.file
    # data is a mapping from lattices to lists of datasets
    self.data = collections.defaultdict(list)
    self.unidentified_data = []
    self.lattice_path = lattice_path
    self.recursive = recursive
    if paths is not None:
      for p in paths:
        self.identifyObject(obj[p], permissive = recursive)
    else:
      self.identifyObject(obj, permissive = recursive)
    # check stuff
    if self.unidentified_data and not recursive:
      raise RuntimeError('got unidentified data: %s' % (self.unidentified_data))

  def identifyObject(self, g, permissive):
    if isinstance(g, h5py.Dataset):
      # if it is a vtk file stored in a string
      if 'TYPE' in g.attrs and g.attrs['TYPE'] in ('VTU_FILE', 'VTK_FILE'):
        assert g.dtype == np.uint8
        self.data[VtkFiles].append(g)
        return
      else:
        ldpath = g.attrs.get('LATTICE_PATH', self.lattice_path)
        if not ldpath:
          if permissive: return
          raise RuntimeError('Error: %s needs a LATTICE_PATH  attribute' % (g.name))
        if not ldpath.startswith('/'):
          ldpath = g[ldpath].name # absolute path
        if not ldpath in self.file:
          if permissive: return
          raise RuntimeError('Error: lattice path %s of %s does not exist in %s' % (ldpath, g.name, f.filename))
        self.data[ldpath].append(g)
        return
    elif isinstance(g, h5py.Group) and permissive:
      for subg in g.itervalues():
        self.identifyObject(subg, permissive)
      return
    # if we get here we don't know what to do with it
    self.unidentified_data.append(g.name)

  def write(self, filename):
    if len(self.data) > 1 or (VtkFiles in self.data and len(self.data[VtkFiles])>1):
      raise RuntimeError('cannot save data on multiple grids')
    k, g = self.data.popitem()
    if k == VtkFiles:
      g = g.pop()
      ext = '.vtk' if g.attrs['TYPE'] == 'VTK_FILE' else '.vtu'
      f = open(os.path.splitext(filename)[0]+ext, 'wb')
      f.write(np.asarray(g).tostring())
      f.close()
    else:
      #ds = vtkcommon.vtkImageDataFromLd(self.file[k].attrs)
      #ld = ku.read_lattice_data_from_hdf(self.file[k])
      fn=str(self.file[k].file.filename)
      path=str(self.file[k].name)
      Pyld = ku.read_lattice_data_from_hdf_by_filename(fn, path)
      
      ds = vtkcommon.vtkImageDataFromLd(Pyld)
      for q in g:
        # iterate over hdf datasets and add them to the image data
        try:
          vtkcommon.vtkImageDataAddData(ds, q, 'CellData', posixpath.basename(q.name))
        except RuntimeError, e:
          print 'Warning: cannot add data %s' % q.name
          print '  Exception reads "%s"' % str(e)
          pass
      writer = vtk.vtkDataSetWriter()
      if int(vtk.vtkVersion.GetVTKSourceVersion()[12])>5:
        writer.SetInputData(ds)
      else:
        writer.SetInput(ds)
      writer.SetFileName(os.path.splitext(filename)[0]+'.vtk')
      writer.Write()
      del ds

  def asVtkDataSets(self, return_ld_path = False):
    import krebsutils
    result = {}
    for k, gg in self.data.iteritems():
      if k == VtkFiles:
        for g in gg:
#          ds = vtkcommon.vtkStringToDataset(
#            np.asarray(g).tostring(),
#            vtk.vtkDataSetReader if g.attrs['TYPE'] == 'VTK_FILE' else self.gridtype2vtureader[g.attrs.get('VTK_DATASET_TYPE', 'vtkUnstructuredGrid')]
#          )
          ds = vtkcommon.vtkDatasetFromHdf5(g)
          result[g.name] = ds
      else:
        #ds = vtkcommon.vtkImageDataFromLd(self.file[k].attrs)
        ld = krebsutils.read_lattice_data_from_hdf(self.file[k])
        ds= vtkcommon.vtkImageDataFromLd(ld)
        name = posixpath.commonprefix([q.name for q in gg])
        for q in gg:
          # iterate over hdf datasets and add them to the image data
          vtkcommon.vtkImageDataAddData(ds, q, 'CellData', posixpath.basename(q.name))
        if not return_ld_path:
          result[name] = ds
        else:
          result[name] = (ds, k)
    return result.values()
    
  def getDatasetPaths(self):
    def it(x0):
      for x1 in x0:
        for x2 in x1:
          yield x2.name
    return list(it(self.data.itervalues()))


if __name__=="__main__":
  from os.path import basename, dirname, join, splitext
  import pprint
  args = sys.argv[1:]
  
  # process
  fn, (grpname, ) = myutils.splitH5PathsFromFilename(args[0])
  f = h5py.File(fn,'r')
  
  e = Extractor(f, [grpname], recursive = True)

  print 'found datasets:'
  pprint.pprint(e.getDatasetPaths())
  
  out = basename(fn);
  if out.endswith('.h5'):
    out = splitext(out)[0]
  if grpname and grpname <> '.':
    tmp = grpname.replace('/','-')
    if not tmp.startswith('-'):
      tmp = '-'+tmp
    out += tmp
  out += '-fields.vtk'  
  
  e.write(out)
