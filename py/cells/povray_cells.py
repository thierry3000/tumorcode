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
import matplotlib
from matplotlib import cm
import h5py
import krebsutils as ku
import itertools
import vapory
def povray_test():
  """ Just a purple sphere """
  scene = vapory.Scene(  vapory.Camera('location',  [0.0, 0.5, -4.0],
                         'direction', [0,0,1.5],
                         'look_at',  [0, 0, 0]),
  
                  objects = [
  
                      vapory.Background("color", [0.85, 0.75, 0.75]),
  
                      vapory.LightSource([0, 0, 0],
                                    'color',[1, 1, 1],
                                    'translate', [-5, 5, -5]),
  
                      vapory.LightSource ([0, 0, 0],
                                      'color', [0.25, 0.25, 0.25],
                                      'translate', [600, -600, -600]),
  
  
                      vapory.Box([-0.5, -0.5, -0.5], [0.5, 0.5, 0.5],
                           vapory.Texture( vapory.Pigment( 'color', [1,0,0]),
                                    vapory.Finish('specular', 0.6),
                                    vapory.Normal('agate', 0.25, 'scale', 0.5)),
                          'rotate', [45,46,47])
                 ]
  )
  # We use antialiasing. Remove this option for faster rendering.
  scene.render("cube.png", width=300, height=300, antialiasing=0.001)
def povray_cells(goodArguments):
  with h5py.File(goodArguments.vbl_simulation_output_filename, 'r') as f:
    h5_cells_grp = f[goodArguments.grp_pattern + "/cells"]
    pos = h5_cells_grp['cell_center_pos']
    pos = np.asarray(pos)
    rad = h5_cells_grp['cell_radii']
    rad = np.asarray(rad)
    o2 = h5_cells_grp['o2']
    o2 = np.asarray(o2)
    x_min= np.min(pos[:,0])
    x_max = np.max(pos[:,0])
    center_x = x_min+0.5*(x_max-x_min)
    y_min= np.min(pos[:,1])
    y_max = np.max(pos[:,1])
    center_y = x_min+0.5*(y_max-y_min)
    z_min= np.min(pos[:,2])
    z_max = np.max(pos[:,2])
    center_z = z_min+0.5*(z_max-z_min)
    print('x: [%f,%f]' % (x_min, x_max))
    print('y: [%f,%f]' % (y_min, y_max))
    print('z: [%f,%f]' % (z_min, z_max))
    print('%f, %f, %f' %(center_x,center_y,center_z))
    #o2 = o2/np.max(o2)
#    x = pos[:,0]
#    y = pos[:,1]
#    z = pos[:,2]
#    s = rad[:,0]
  camera = vapory.Camera('location', [700,700,-700], 'look_at', [0,0,0])
  light = vapory.LightSource([1000,-1000,-1000], 'color', [1, 1, 1])
  light2 = vapory.LightSource([0,0,0], 'color',[1, 1, 1], 'translate', [1000,-1000,-1000] )
  light3 = vapory.LightSource([500,-1000,500], 'color', [1, 1, 1] )
  myObjectList = []
  myObjectList.append(light)
  myObjectList.append(light2)
  myObjectList.append(light3)
  
  cuttingY = vapory.Plane([0,1,0], 0,)
  cuttingX = vapory.Plane([1,0,0], -1,)
  max_rad = np.max(rad)
  max_o2 = np.max(o2)
  n= 10000
  for (aPosition, aRadius, aO2Value) in zip(pos[0:n], rad[0:n], o2[0:n]):
    thisSphere = vapory.Sphere( aPosition, aRadius[0])
    color = matplotlib.cm.hsv(aO2Value[0]/max_o2)
    #print(color[0:3])
    #cuttedSphere = vapory.Intersection(thisSphere, cuttingY, vapory.Texture( vapory.Pigment( 'color', color[0:3]  )))    
    #cuttedSphere = vapory.Intersection(thisSphere, cuttingY, cuttingX)    
    #cuttedSphere = thisSphere  
    #myObjectList.append(cuttedSphere)
    #myObjectList.append(thisSphere)
   # myObjectList.append(vapory.Sphere( aPosition, aRadius[0], vapory.Texture( vapory.Pigment( 'color', matplotlib.cm.Blues(aO2Value[0]/max_o2) ))))
    myObjectList.append(vapory.Sphere( aPosition, aRadius[0], vapory.Texture( vapory.Pigment( 'color', [1,0,0] ))))
    
  scene = vapory.Scene( camera, objects= myObjectList,  defaults = [vapory.Finish( 'ambient', 1.5)],)
  scene.render("purple_sphere.png", width=400, height=300,  antialiasing=0.01, remove_temp=True)
    
  
  

if __name__ == "__main__":
  import argparse
  parser = argparse.ArgumentParser(description='Plot MTS with mayavi')  
  parser.add_argument('vbl_simulation_output_filename', type=str, default='check_types.h5', help='output file name in hdf5 format')
  parser.add_argument('grp_pattern')
  
  interactive = True;
  goodArguments, otherArguments = parser.parse_known_args()
  
  goodArguments.grp_pattern = 'out0490'
  #goodArguments.output_grp_name = 'out0458'
  
  povray_cells(goodArguments)
  povray_test()
  #plot_vessels(goodArguments)