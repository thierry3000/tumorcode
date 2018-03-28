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
    t = np.linspace(0, 4 * np.pi, 20)
    cos = np.cos
    sin = np.sin

    x = sin(2 * t)
    y = cos(t)
    z = cos(2 * t)
    s = 2 + sin(t)

    return mlab.points3d(x, y, z, s, colormap="copper", scale_factor=.25)
    
def test_points3d_gauss():
    n=20
    pts = np.random.randn(n,3)
    s = np.sqrt(np.sum(pts**2,1))
    
    no_grid_points = 10
    x_min = np.min(pts[:,0])
    x_max = np.max(pts[:,0])
    x_dist=(x_max-x_min)/no_grid_points
    
    y_min = np.min(pts[:,1])
    y_max = np.max(pts[:,1])
    y_dist=(y_max-y_min)/no_grid_points
    
    z_min = np.min(pts[:,2])
    z_max = np.max(pts[:,2])
    z_dist=(z_max-z_min)/no_grid_points
    
    
    X, Y, Z = np.meshgrid(np.arange(x_min,x_max,x_dist), np.arange(y_min,y_max,y_dist), np.arange(z_min,z_max,z_dist))
    
    from scipy.interpolate import griddata
    
    grid1 = griddata(pts,s,(X,Y,Z), method='nearest')
    obj = mlab.contour3d(grid1, contours=1, transparent=True)
    return obj


def test_contour3d():
    n=20
    pts = np.random.randn(3,20)
    s = np.sum(pts**2,0)
    

    obj = mlab.contour3d(pts[0,:],pts[1,:],pts[2,:], s, contours=0.5, transparent=True)
    return obj
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
  
@mlab.show
def plot_contour_cells(goodArguments):
  with h5py.File(goodArguments.vbl_simulation_output_filename, 'r') as f:
    h5_cells_grp = f[goodArguments.output_grp_name + "/cells"]
    pos = h5_cells_grp['cell_center_pos']
    pos = np.asarray(pos)
    rad = h5_cells_grp['cell_radii']
    rad = np.asarray(rad)
    o2 = h5_cells_grp['o2']
    o2 = np.asarray(o2)
    o2 = o2[:,0]
    x = pos[:,0]
    y = pos[:,1]
    z = pos[:,2]
    s = rad[:,0]
    pts = mlab.quiver3d(x,y,z, s,s,s, scalars=o2, colormap="blue-red", scale_factor=2, mode='sphere')
    pts.glyph.color_mode = 'color_by_scalar'
    #pts.glyph.glyph_source.glyph_source.center = [0,0,0]
    
    
    ''' interpolate to grid '''
    
    no_grid_points = 100
    x_min = np.min(pos[:,0])
    x_max = np.max(pos[:,0])
    x_dist=(x_max-x_min)/no_grid_points
    
    y_min = np.min(pos[:,1])
    y_max = np.max(pos[:,1])
    y_dist=(y_max-y_min)/no_grid_points
    
    z_min = np.min(pos[:,2])
    z_max = np.max(pos[:,2])
    z_dist=(z_max-z_min)/no_grid_points
    
    
    X, Y, Z = np.meshgrid(np.arange(x_min,x_max,x_dist), np.arange(y_min,y_max,y_dist), np.arange(z_min,z_max,z_dist))
    
    from scipy.interpolate import griddata
    
    grid1 = griddata(pos,o2,(X,Y,Z), method='linear')
    isosurface = (np.max(o2)-np.min(o2))/2
    print("isosurface at: %f" % isosurface)
    obj = mlab.contour3d(grid1, contours=[isosurface], transparent=True)
    

if __name__ == "__main__":
  import argparse
  parser = argparse.ArgumentParser(description='Plot MTS with mayavi')  
  parser.add_argument('--s',dest='vbl_simulation_output_filename', type=str, default='check_types.h5', help='output file name in hdf5 format')
  
  interactive = True;
  goodArguments, otherArguments = parser.parse_known_args()
  
  goodArguments.output_grp_name = 'out0100'
  #goodArguments.output_grp_name = 'out0458'
  
  #test_points3d_gauss()
  #test_contour3d()
  
  #plot_cells(goodArguments)
  plot_contour_cells(goodArguments)
  
  #plot_vessels(goodArguments)
  
  mlab.show()