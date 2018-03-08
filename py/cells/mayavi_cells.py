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

import numpy as np
import mayavi.mlab as mlab
mlab.options.backend = 'envisage'
import h5py
import krebsutils as ku
import itertools
'''
possible color maps:

accent       flag          hot      pubu     set2
autumn       gist_earth    hsv      pubugn   set3
black-white  gist_gray     jet      puor     spectral
blue-red     gist_heat     oranges  purd     spring
blues        gist_ncar     orrd     purples  summer
bone         gist_rainbow  paired   rdbu     winter
brbg         gist_stern    pastel1  rdgy     ylgnbu
bugn         gist_yarg     pastel2  rdpu     ylgn
bupu         gnbu          pink     rdylbu   ylorbr
cool         gray          piyg     rdylgn   ylorrd
copper       greens        prgn     reds
dark2        greys         prism    set1
'''
def test_points3d():
    t = numpy.linspace(0, 4 * numpy.pi, 20)
    cos = numpy.cos
    sin = numpy.sin

    x = sin(2 * t)
    y = cos(t)
    z = cos(2 * t)
    s = 2 + sin(t)

    return mlab.points3d(x, y, z, s, colormap="copper", scale_factor=.25)
    
@mlab.show
def plot_cells(goodArguments):
  with h5py.File(goodArguments.vbl_simulation_output_filename, 'r') as f:
    h5_cells_grp = f[goodArguments.output_grp_name + "/cells"]
    pos = h5_cells_grp['cell_center_pos']
    pos = np.asarray(pos)
    rad = h5_cells_grp['cell_radii']
    rad = np.asarray(rad)
    o2 = h5_cells_grp['o2']
    o2 = np.asarray(o2)
    x = pos[:,0]
    y = pos[:,1]
    z = pos[:,2]
    s = rad[:,0]
    pts = mlab.quiver3d(x,y,z, s,s,s, scalars=o2, colormap="blue-red", scale_factor=2, mode='sphere')
    pts.glyph.color_mode = 'color_by_scalar'
    pts.glyph.glyph_source.glyph_source.center = [0,0,0]
    
    
  ''' CELLS READY '''
  #this ready the graph info from hdf
  with h5py.File(goodArguments.vbl_simulation_output_filename, 'r') as f:
    h5_vessel_grp = f[goodArguments.output_grp_name + "/vessels"]
    graph = ku.read_vessels_from_hdf(h5_vessel_grp, ['position', 'flags', 'radius', 'nodeflags'], return_graph=True)
    
  
  node_pos=np.asarray(graph['position'])
  myEdges = graph.edgelist
  ''' create radii at the vertices '''
  num_verts = len(node_pos)
  node_rad = np.zeros(num_verts)
  node_n   = np.zeros(num_verts)
  for r,(a,b) in itertools.izip(np.asarray(graph['radius']),np.asarray(myEdges)): # itertools, or it will blow up in yr face cuz ultra slow python lists and shit will be created
    node_rad[a] += r
    node_rad[b] += r
    node_n[a] += 1.
    node_n[b] += 1.
  node_rad /= np.maximum(1, node_n)
  
  ''' create points data set, similar to vtk points data set
      we take the radius as associated scalar
  '''
  pts = mlab.points3d(np.asarray(node_pos[:,0]), np.asarray(node_pos[:,1]), np.asarray(node_pos[:,2]), node_rad, scale_factor=1.0)
  #create lines
  pts.mlab_source.dataset.lines = np.array(graph.edgelist)
  #same as in the constructor! could be changed here
  pts.mlab_source.scalars = node_rad
  #set up a tube pipeline
  tube = mlab.pipeline.tube(pts)
  #I use this construction to the the possible options
  if 0:
    alloptions = dir(tube.filter)
    for option in alloptions:
      print(option)
  #tube.filter.set_input_data=node_rad
  #tube.filter.radius= 42
  #tube.filter.radius_factor = 1.
      
  tube.filter.vary_radius = 'vary_radius_by_absolute_scalar'
  
  ''' colors '''
  #mlab.pipeline.surface(tube, color=(0.8, 0.8, 0))
  mlab.pipeline.surface(tube, colormap="blue-red")
  
  mlab.colorbar()
  #mlab.savefig("test.eps")
  
  

if __name__ == "__main__":
  import argparse
  parser = argparse.ArgumentParser(description='Plot MTS with mayavi')  
  parser.add_argument('--s',dest='vbl_simulation_output_filename', type=str, default='check_types.h5', help='output file name in hdf5 format')
  
  interactive = True;
  goodArguments, otherArguments = parser.parse_known_args()
  
  goodArguments.output_grp_name = 'out0100'
  #goodArguments.output_grp_name = 'out0458'
  
  plot_cells(goodArguments)
  
  #plot_vessels(goodArguments)