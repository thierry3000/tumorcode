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
'''
very nice code which calculates a local mvd field
'''

import os, sys
import numpy as np
import h5py
import extensions

from matplotlib import pyplot

import krebsutils
import myutils
import mpl_utils

if __name__ == '__main__':
  filename = sys.argv[1]
  grouppath = sys.argv[2]
  
  with h5py.File(filename, 'r') as file:
    vesselgroup = file[grouppath]
    ldvessels = krebsutils.read_lattice_data_from_hdf(vesselgroup['lattice'])
    wbbox     = ldvessels.worldBox
    
    graph = krebsutils.read_vesselgraph(vesselgroup, ['position', 'flags'])
    graph = graph.get_filtered(myutils.bbitwise_and(graph['flags'], krebsutils.CIRCULATED))

  print 'vessel ld:'
  print ldvessels
  ''' this splits splits space in lattices of 300.
      the second 300 adds a 'safety layer of 100. mum
      so we do not consider the outermost data for calculating the
      actual mvd
  '''
  fieldld = krebsutils.SetupFieldLattice(wbbox, 3, 20., 100.)
  wbbox     = fieldld.worldBox
  print 'field ld:'
  print fieldld
  z = fieldld.shape[2]/2

  weights = krebsutils.sample_edges_weights(graph.nodes['position'], graph.edgelist, 30.)
  positions = krebsutils.sample_edges(graph.nodes['position'], graph.edgelist, graph.nodes['position'], 30., krebsutils.VesselSamplingFlags.DATA_PER_NODE | krebsutils.VesselSamplingFlags.DATA_LINEAR)

  eps = 1.0-1.e-15
  x0, x1, y0, y1, z0, z1 = wbbox
  ranges = [np.arange(x0, x1, fieldld.scale*eps),
            np.arange(y0, y1, fieldld.scale*eps),
            np.arange(z0, z1, fieldld.scale*eps),
            ]
  print 'histogram bin ends:',map(lambda r: (r.min(), r.max()), ranges)
  mvd, _ = np.histogramdd(positions, bins = ranges, weights = weights)
  mvd *= 1.e6/(fieldld.scale**3)
  print 'result shape:',mvd.shape
  print('average mvd')
  print(np.mean(mvd[1:-1,1:-1,1:-1]))
  
  fig, ax = pyplot.subplots(1)
  plt = ax.imshow(mvd[:,:,z], interpolation = 'none')
  ax.set(title = 'MVD')
  divider = mpl_utils.make_axes_locatable(ax)
  cax = divider.append_axes("right", size = "5%", pad = 0.05)
  fig.colorbar(plt, cax = cax)    
  
  pyplot.show()