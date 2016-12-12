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
"""
Created on Thu Aug 11 16:09:05 2016

@author: mwelter
"""
'''
NOTE:
since synchronizing with the repository
is to tedious when developing, we decided
to spare rendering settings in the
repository.

Handle:
copy povrayRenderSettingsTemplate.py
to
povrayRenderSettings.py

and make changes as you wish.
'''
import numpy as np
if __name__ == '__main__':
  import os.path, sys
  sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'..'))
import krebsutils
import matplotlib
from krebs.povrayRenderVessels import  cm_redblue

def colorfactory_vessel_debug(vesselgraph):
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

image = dict(
      #res=(3600,3600), # brauche ca. 3600 fuer 150mm bei 600 dpi
      #res = (1800, 1800), # 75 mm 600 dpi
      #aa=4,
      res=(1024, 1024), aa=1,
      background = 1.,
      out_alpha=False,
      num_threads=3,
)

tumor = dict(
      #cam='topdown_slice',
      #cam='corner',
      #clip_zslice=(0.4, 0.6),
      cam_distance_multiplier = 1.0,
      cam='topdown_slice',
      debug=False,
      colored_slice = True,
)

vessels = dict(
      #clip_pie=0,
      #clip_zslice=(0.4, 0.6),
      #cam='corner',
      cam='topdown_slice',
      #cam = 'topdown',
      colored_slice = True,
      background = 1.,
      #ambient_color = (0.5, 0.5, 0.5),
      ambient_color = (0.3, 0.3, 0.3),
)
dbg_vessels = dict(
      res=(1024,1024),
      aa=3,
      #clip_pie=0,
      #clip_zslice=(0.4, 0.6),
      #cam='corner',
      cam='topdown',
      colored_slice=False,
      out_alpha=False,
      num_threads=7,
      debug=True,
      background = 0.0,
      colorfactory = colorfactory_vessel_debug,
      overlay = False,
      #colorfactory = colorfactory,
      #keep_files = True,
      #temp_file_dir = '.',
)