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
if __name__ == '__main__':
  import os.path, sys
  sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'..'))

import os,sys
from os.path import join, basename, dirname, splitext
import time
import krebsutils
import h5py
import posixpath
import numpy as np
import matplotlib.cm
import myutils
from copy import deepcopy
import math

#import krebs.povrayRenderVessels
from krebs.povrayRenderVessels import  addVesselTree, cm_redblue
from krebs.povrayEasy import *


colors = [(0., 0.2, 0., 0.), (1., 0., 1.0, 0.) ]
cm_gf = matplotlib.colors.LinearSegmentedColormap('', {
  'red'   : [(x, r, r) for (x,r,g,b) in colors],
  'green' : [(x, g, g) for (x,r,g,b) in colors],
  'blue'  : [(x, b, b) for (x,r,g,b) in colors],
}, N = 256, gamma = 1.)


def colorfactory(vesselgraph):
  is_set = lambda flags_,flag: np.asarray(np.bitwise_and(flags_, flag), np.bool)
  
  edges = vesselgraph.edgelist
  num_nodes = len(vesselgraph.nodes['position'])
  flags = vesselgraph.edges['flags']
  data = vesselgraph.nodes['pressure']
  nflags = vesselgraph.nodes['nodeflags']
  #action = vesselgraph.nodes['action']
  #data = vesselgraph.edges['shearforce']
  #data = krebsutils.edge_to_node_property(num_nodes, edges, data, 'avg')  
  #nflags = krebsutils.edge_to_node_property(num_nodes, edges, flags, 'and')
  #ncirc  = krebsutils.edge_to_node_property(num_nodes, edges, np.bitwise_and(flags, krebsutils.CIRCULATED), 'or')
  #nflags = np.bitwise_or(nflags, ncirc)

  gray = np.asarray((0.1,0.1,0.1))
  ltgray = np.asarray((0.0, 0.8, 0.0))
  boundarycol = np.asarray((1., 1., 0))

  circulated = is_set(flags,krebsutils.CIRCULATED)
  has_circ = circulated.max()==True
  if 0:
    DIE = 0
    GROW = 1
    IDLE = 2
    edgecolors = np.zeros((len(edges),3), dtype = np.float)
    nodecolors = np.zeros((num_nodes,3), dtype = np.float)
    edgecolors[:,:] = gray
    nodecolors[:,:] = gray
    nodecolors[action == DIE] = np.asarray((1.,0.,0.))
    nodecolors[action == GROW] = np.asarray((0.,1.,0.))
    nodecolors[action == IDLE] = np.asarray((0.,0.,1.))
    
  elif has_circ:
    ncirculated = is_set(nflags,krebsutils.CIRCULATED)
    p0 = np.amin(data[ncirculated])
    p1 = np.amax(data[ncirculated])
    cm = matplotlib.cm.ScalarMappable(cmap=cm_redblue)
    cm.set_clim(p0, p1)

    edgedata = np.average((data[edges[:,0]], data[edges[:,1]]), axis=0)
    edgecolors = cm.to_rgba(edgedata)[:,:3]
    edgecolors[np.where(circulated==0)] = gray

    nodecolors = cm.to_rgba(data)[:,:3]
    nodecolors[np.where(ncirculated==0)] = gray
  else:
    edgecolors = np.zeros((len(edges),3), dtype = np.float)
    nodecolors = np.zeros((num_nodes,3), dtype = np.float)
    edgecolors[is_set(flags,krebsutils.ARTERY),0] = 1
    nodecolors[is_set(nflags,krebsutils.ARTERY),0] = 1
    edgecolors[is_set(flags,krebsutils.VEIN),2] = 1
    nodecolors[is_set(nflags,krebsutils.VEIN),2] = 1
  edgecolors[is_set(flags,krebsutils.CAPILLARY)] = ltgray
  nodecolors[is_set(nflags,krebsutils.CAPILLARY)] = ltgray
  nodecolors[is_set(nflags,krebsutils.BOUNDARY)] = boundarycol
  edgecolors[is_set(flags,krebsutils.BOUNDARY)] = boundarycol

  vesselgraph.edges['colors'] = edgecolors
  vesselgraph.nodes['colors'] = nodecolors


def renderSliceWithDistribution(vesselgroup, imagefn, **kwargs):
  vessel_ld = krebsutils.read_lattice_data_from_hdf(vesselgroup['vessels/lattice'])
  vessel_graph = krebsutils.read_vessels_from_hdf(vesselgroup['vessels'], ['position', 'flags', 'radius', 'pressure', 'shearforce', 'nodeflags'], return_graph=True)
  vessel_graph.edges['radius'] *= 4.

  kwargs = deepcopy(kwargs)
  wbbox = vessel_ld.worldBox
  trafo = calc_centering_normalization_trafo(wbbox)
  height = (wbbox[5]-wbbox[4])*trafo.w

  print 'Vessel BBox:', vessel_ld.worldBox
  print vessel_ld
  print 'Post Trafo Ld BBox:', transform_ld(trafo, vessel_ld).worldBox

  hasGfField = 'field_ld' in vesselgroup and vessel_ld.shape[2] == 1
  #hasGfField = False
  if hasGfField:
    volume_ld = krebsutils.read_lattice_data_from_hdf(vesselgroup['field_ld'])
    print 'Volume BBox:', volume_ld.worldBox
    print volume_ld
    volumedata = np.asarray(vesselgroup['gf'])
    #print volumedata.shape, volumedata.min(), volumedata.max()
    volume_ld = transform_ld(trafo, volume_ld)
    print 'Post Trafo Volume BBox:', volume_ld.worldBox
    print 'volume data bounds:', volumedata.min(), volumedata.max()

  colorfactory(vessel_graph)

  with EasyPovRayRender(**kwargs) as epv:
    epv.setBackground(kwargs.pop('background',0.0))

    cam_fov = 60.
    cam_distance_factor = ComputeCameraDistanceFactor(cam_fov, kwargs['res'], wbbox)

    epv.setCamera((0,0,cam_distance_factor*1.05), lookat = (0,0,0), fov = cam_fov, up = 'y')

    epv.addLight(10.*Vec3(1,0.5,2), 1.2)

    cm = matplotlib.cm.ScalarMappable(cmap = cm_gf)
    cm.set_clim(-0.01,1.01)
    pvcm = matplotlibColormapToPovray('DATACOLORMAP', cm)
    epv.declareColorMap(pvcm)

    if hasGfField:
      volumedata = epv.declareVolumeData(volumedata, volume_ld.worldBox, volume_ld)
      epv.addVolumeDataSlice(volumedata, (0,0,0), (0, 0, 1.), pvcm)
    addVesselTree(epv, vessel_graph, trafo = trafo, **kwargs)

    imagefn = epv.render(imagefn)




if __name__ == '__main__':
  sys.path.append(join(dirname(__file__),'..'))

if (__name__ == '__main__'):
  fn = sys.argv[1]
  pattern = sys.argv[2]

  f = h5py.File(fn,'r')
  dirs = myutils.walkh5(f['/'], pattern) #'*vessels')

  #for group in f.itervalues():
  for grpname in dirs:
    vesselgrp = f[grpname]
    renderSliceWithDistribution(vesselgrp, splitext(basename(fn))[0]+posixpath.basename(vesselgrp.name),
      res=(1024,1024), aa=3,
      #clip_pie=0,
      #clip_zslice=(0.4, 0.6),
      #cam='corner',
      cam='topdown',
      colored_slice=False,
      out_alpha=False,
      num_threads=4,
      debug=True,
      colorfactory = colorfactory,
      #keep_files = True,
      #temp_file_dir = '.',
      )
