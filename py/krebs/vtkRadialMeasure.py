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
#!/usr/bin/env python

"""
This does radial averaging of .vtu and .vtk files.
Meant for stuff i did with the deal.ii lib!
"""


import sys, os
import numpy as np
from vtk import *
from vtkcommon import *
from hdfvessels2vtk import removeUncirculatedVessels
import math



def computeRadius(center, bounds):
  bounds = bounds.copy()
  for q in bounds:
    q -= center
  bounds = np.max(np.abs(bounds), axis=0)
  return np.max(bounds) #np.linalg.norm(bounds)


def computeCellBinAndWeight(c, center, bin_size):
  bounds = np.asarray(c.GetBounds()).reshape((3,2)).transpose()
  size = bounds[1]-bounds[0]
  mp = np.average(bounds, axis=0)
  dim = c.GetCellDimension()
  if dim == 0:
    w = 1.
  elif dim == 1:
    w = np.linalg.norm(size)
  elif dim == 2:
    w = 0.33*(size[0]*size[1]+size[1]*size[2]+size[0]*size[2])
  elif dim == 3:
    w = size[0]*size[1]*size[2]
  mp -= center
  r = int(np.linalg.norm(mp)/bin_size)
  #print mp, size, dim, w, r
  return r, w


def computeBinVolumes(center, bounds, shell_width, num_bins):
  bounds = bounds.copy()
  for q in bounds:
    q -= center
  cellsize_base = shell_width*0.25
  num_cells = np.zeros(center.shape, dtype=np.int)
  cellsize = np.zeros(center.shape, dtype=np.float)
  for i, (a, b) in enumerate(bounds.transpose()):
    num_cells[i] = (b-a)/cellsize_base
    cellsize[i] = (b-a)/num_cells[i]
  arr = np.zeros(num_bins, dtype=np.double)
  for idx in np.ndindex(*num_cells):
    wpos = cellsize * (np.asarray(idx, dtype=np.float)+0.5) + bounds[0]
    r = int(np.linalg.norm(wpos)/shell_width)
    if r < num_bins:
      arr[r] += 1.
  return np.prod(cellsize)*arr


def radialMeasureCellData(dataset, center, bounds, bin_size):
  """compute average data in concentric shells"""
  num_cells = dataset.GetNumberOfCells()
  num_points = dataset.GetNumberOfPoints()

  #print dataset

  max_radius = computeRadius(center, bounds)
  num_bins = max(1, int(max_radius / bin_size)+1)

  celldata = dataset.GetCellData()
  num_data = celldata.GetNumberOfArrays()

  print 'Bins, center=%s, bins=%i, radius=%f, num cells=%i, bin size=%f' % (center, num_bins, max_radius, num_cells, bin_size)
  # precompute weights and bins for each cell
  rad_weights = np.zeros(num_bins) # normalization constants, one for each bin
  cell_weights = np.zeros(num_cells) # integration weight attached to each cell
  cell_bin     = -np.ones(num_cells, dtype=np.int) # the bin each cell falls into
  for i in xrange(num_cells):
    c = dataset.GetCell(i)
    r, w = computeCellBinAndWeight(c, center, bin_size)
    if r >= num_bins:
      continue
    cell_weights[i] = w
    cell_bin[i] = r
    rad_weights[r] += w
  max_w = np.max(cell_weights)
  cell_weights *= (1./max_w)
  rad_weights *= (1./max_w)

  res_data = []
  res_names = []
  for a in vtkIterArrays(celldata):
    # go through all cell data arrays
    # check if it makes sense to do the measurement
    nc = a.GetNumberOfComponents()
    if nc > 1:
      print "skipping data '%s' because number of components %i > 1" % (a.GetName(), nc)
      continue
    if a.GetDataTypeAsString() not in ('float', 'double'):
      print "skipping data '%s' because datatype is '%s'" % (a.GetName(), a.GetDataTypeAsString())
      continue
    res_names.append(a.GetName())
    print "measuring", a.GetName()
    # allocate result array
    loc_res = np.zeros((2, num_bins))
    getfunc = getattr(a, "GetValue")
    # iterate over all cells, get the data and insert it into the bins
    for cellnum in xrange(num_cells):
      bin = cell_bin[cellnum]
      if bin < 0:
        continue
      w = cell_weights[cellnum]
      val = getfunc(cellnum)
      loc_res[0, bin] += w * val
      loc_res[1, bin] += w * val * val  # <val^2> for std deviation
    # divide by normalization constants
    for q in loc_res:
      np.divide(q, rad_weights, q)
    # compute std dev
    loc_res = np.nan_to_num(loc_res)
    loc_res[1] = np.sqrt(loc_res[1] - np.power(loc_res[0], 2.))
    res_data.append(loc_res)

  # x coordinates of the bins
  bins = (np.arange(num_bins)+0.5) * bin_size

  # unnormalize -> lengh, area, or volume (depending on the cell dimension) per bin
  rad_weights *= max_w
  
  return bins, rad_weights, res_names, res_data



def vtkFileToRadialHdf(in_filename, h5_group, shell_width = 100.):
  fn = in_filename
  
  reader = vtkXMLUnstructuredGridReader() if fn.endswith(".vtu") else vtkDataSetReader()
  reader.SetFileName(fn)
  reader.Update();
  dataset = reader.GetOutput()
  del reader
  dataset.ComputeBounds()

  ptc_algo = vtkPointDataToCellData()
  ptc_algo.PassPointDataOn()
  ptc_algo.SetInput(dataset)
  ptc_algo.Update()
  dataset = ptc_algo.GetOutput()

  #print dataset

  center = np.asarray((0.,0.,0.))

  bounds = vtkGetDataSetBounds(dataset)

  if 'vessel' in os.path.basename(fn):
    print 'vessel data detected -> filtering uncirculated'
    dataset = removeUncirculatedVessels(dataset)

  xarr, dens, names, data = radialMeasureCellData(dataset, center, bounds, shell_width)

  print "computing bin volumes"
  binvol = computeBinVolumes(center, bounds, shell_width, len(xarr))

  for i in xrange(len(dens)):
    vol = binvol[i]
    #print vol, dens[i], dens[i]/vol
    dens[i] /= vol
  data.append((dens, np.zeros_like(dens)))
  names.append('grid_density')

  #f = open(outfn+".dat",'w')
  #for i, (name, (yarr, yyarr)) in enumerate(zip(names, data)):
    #print "data @index %i = %s" % (i, name)
    #print >>f, '# %s' % name
    #for q in zip(xarr, yarr, yyarr):
      #print >>f, '%f %f %f' % q
    #print >>f, '\n'
  #f.close()

  #f = h5py.File(outfn+'.h5', outmode)
  #f.attrs['SOURCE'] = fn
  #g = f['/']
  #if outgroup:
    #g = g.require_group(outgroup)
  g = h5_group
  for i, (name, (yarr, yyarr)) in enumerate(zip(names, data)):
    ds = g.create_dataset(name, data = np.asarray((xarr, yarr, yyarr)))
    ds.attrs['TYPE'] = 'PLOT_XYE'
    ds.attrs['XLABEL'] = 'r'



if __name__ == '__main__':
  import h5py
  import getopt
  
  def exitmsg(msg):
      print msg
      sys.exit(0)
  usage = """
  %s [-o | -a filename] [-g groupname] input
  """ % sys.argv[0]
  
  try:
    opt, args = getopt.getopt(sys.argv[1:],'o:a:g:',['help'])
  except getopt.GetoptError, e:
    exitmsg(usage)

  help = False
  outmode = 'w'  # append or write
  outfn = None
  outgroup = ''
  for o, a in opt:
    if o == '-o':
      outmode = 'w'
      outfn = a
    elif o == '-a':
      outmode = 'a'
      outfn = a
    elif o == '-g':
      outgroup = a
    elif o == '--help':
      help = True
  if help:
    exitmsg(usage)

  if len(args) > 1:
    exitmsg(usage)

  fn = args[0]

  if outfn is None:
    if not (fn.endswith('.vtk') or fn.endswith('.vtu')):
      exitmsg('input filename must end with .vtk or .vtu')
    outfn = os.path.basename(fn)+"-rad"
  else:
    if not any((outfn.endswith('.dat'), outfn.endswith('.h5'))):
      exitmsg('out filename must end with .dat or .h5') 
  outfn = os.path.splitext(outfn)[0]

  f = h5py.File(outfn+'.h5', outmode)
  f.attrs['SOURCE'] = fn
  g = f['/']
  if outgroup:
    g = g.require_group(outgroup)

  vtkFileToRadialHdf(fn, g)
