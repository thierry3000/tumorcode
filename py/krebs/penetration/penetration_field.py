#!/usr/bin/env python2
# -*- coding: utf-8 -*-
'''
This file is part of tumorcode project.
(http://www.uni-saarland.de/fak7/rieger/homepage/research/tumor/tumor.html)

Copyright (C) 2017  Michael Welter and Thierry Fredrich

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
from os.path import join, basename, dirname, splitext
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'..'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../..'))

import os, sys
import numpy as np
import h5py
#import extensions

from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot

import krebsutils
import myutils
import mpl_utils

def sample_vessel_system(goodArguments):
  filename = goodArguments.vesselFileNames;
  grouppath = goodArguments.grp_pattern
      
  
  with h5py.File(filename, 'r') as file:
    if( 'vessels' in grouppath):
      print('found vessels!')
      vesselgroup = file[grouppath]
    else:
      if( 'out' in grouppath):
        outgroup = file[grouppath]
        print('found tumor of type: %s' % str(outgroup['tumor'].attrs.get('TYPE')))
        vesselgroup = file[grouppath+'/vessels']
      else:
        print("unknown data structure!")
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
  sampling_lattice_spacing = goodArguments.sampling_lattice_spacing
  fieldld = krebsutils.SetupFieldLattice(wbbox, 3, sampling_lattice_spacing, 100.)
  wbbox     = fieldld.worldBox
  print 'field ld:'
  print fieldld
  z = fieldld.shape[2]/2

  longitudinal_sampling_distance = goodArguments.longitudinal;
  weights = krebsutils.sample_edges_weights(graph.nodes['position'], graph.edgelist, longitudinal_sampling_distance)
  positions = krebsutils.sample_edges(graph.nodes['position'], graph.edgelist, graph.nodes['position'], longitudinal_sampling_distance, krebsutils.VesselSamplingFlags.DATA_PER_NODE | krebsutils.VesselSamplingFlags.DATA_LINEAR)


  
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
  
  ''' new stuff '''
  from scipy import ndimage
   
  bool_field = mvd>0;
  bool_field = np.logical_not(bool_field)
  distance_map = ndimage.morphology.distance_transform_edt(bool_field)
  distance_map = distance_map*sampling_lattice_spacing
  
#  fig, ax = pyplot.subplots(1)
#  plt = ax.imshow(mvd[:,:,z], interpolation = 'none')
#  ax.set(title = 'MVD')
#  divider = mpl_utils.make_axes_locatable(ax)
#  cax = divider.append_axes("right", size = "5%", pad = 0.05)
#  fig.colorbar(plt, cax = cax)
  
  fig, ax = pyplot.subplots(1)
  plt = ax.imshow(distance_map[:,:,z], interpolation = 'none')
  #ax.set(title = 'Distance Map \n group: %s, file: %s' %(grouppath, filename))
  ax.set(title = 'Distance Map \n smp_logitudinal: %s, lattice_const: %s' %(longitudinal_sampling_distance,sampling_lattice_spacing))
  divider = mpl_utils.make_axes_locatable(ax)
  cax = divider.append_axes("right", size = "5%", pad = 0.05)
  fig.colorbar(plt, cax = cax)
  basename(filename)
  with PdfPages('distmap_' + basename(filename)+'_' + grouppath + '.pdf') as pdf:
    pdf.savefig(fig)
  #pyplot.show()
  
if __name__ == "__main__":
  import argparse
  parser = argparse.ArgumentParser(description='Plot/ Analyze infos penetration distances.',
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('vesselFileNames',type=str, help='Vessel file to calculate')  
  parser.add_argument('grp_pattern',type=str, help='Where to find the vessel group in the file')    
  parser.add_argument("-l","--longitudinal", type=float, help="logitudinal sampling distance", default=5.)
  parser.add_argument("-s","--sampling_lattice_spacing", type=float, help="sampling distance of lattice", default=20.)  
  goodArguments, otherArguments = parser.parse_known_args()
  
  dirs = set()
  #for fn in filenames: # only one file, so far
  with h5py.File(goodArguments.vesselFileNames, 'r') as f:
    d = myutils.walkh5(f, goodArguments.grp_pattern)
    assert len(d), 'your pattern "%s" is not found in "%s"!' % ( goodArguments.grp_pattern, goodArguments.vesselFileNames)
    dirs =set.union(dirs, d)
  print 'and resolved groups therein: %s' % ','.join(dirs)
  
  for aDir in dirs:
    goodArguments.grp_pattern = str(aDir)
    sample_vessel_system(goodArguments)