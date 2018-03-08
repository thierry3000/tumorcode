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

def povray_cells(goodArguments):
  with h5py.File(goodArguments.vbl_simulation_output_filename, 'r') as f:
    h5_cells_grp = f[goodArguments.output_grp_name + "/cells"]
    pos = h5_cells_grp['cell_center_pos']
    pos = np.asarray(pos)
    rad = h5_cells_grp['cell_radii']
    rad = np.asarray(rad)
    o2 = h5_cells_grp['o2']
    o2 = np.asarray(o2)
    #o2 = o2/np.max(o2)
    x = pos[:,0]
    y = pos[:,1]
    z = pos[:,2]
    s = rad[:,0]
  camera = vapory.Camera('location', [40,40,-80], 'look_at', [0,0,0])
  light = vapory.LightSource( [20,20,-50], 'color', 1 )
  myObjectList = []
  myObjectList.append(light)
  
  cuttingY = vapory.Plane([0,1,0], 0,)
  cuttingX = vapory.Plane([1,0,0], -1,)
  max_rad = np.max(rad)
  max_o2 = np.max(o2)
  for (aPosition, aRadius, aO2Value) in zip(pos, rad, o2):
    thisSphere = vapory.Sphere( aPosition, aRadius[0])
    color = matplotlib.cm.hsv(aO2Value[0]/max_o2)
    print(color[0:3])
    cuttedSphere = vapory.Intersection(thisSphere, cuttingY, vapory.Texture( vapory.Pigment( 'color', color[0:3]  )))    
    #cuttedSphere = vapory.Intersection(thisSphere, cuttingY, cuttingX)    
    #cuttedSphere = thisSphere  
    myObjectList.append(cuttedSphere)
    #myObjectList.append(vapory.Sphere( aPosition, aRadius[0], vapory.Texture( vapory.Pigment( 'color', matplotlib.cm.Blues(aO2Value[0]) ))))

  scene = vapory.Scene( camera, objects= myObjectList)
  scene.render("purple_sphere.png", width=400, height=300, antialiasing=0.001, remove_temp=False)
    
  
  

if __name__ == "__main__":
  import argparse
  parser = argparse.ArgumentParser(description='Plot MTS with mayavi')  
  parser.add_argument('--s',dest='vbl_simulation_output_filename', type=str, default='check_types.h5', help='output file name in hdf5 format')
  
  interactive = True;
  goodArguments, otherArguments = parser.parse_known_args()
  
  goodArguments.output_grp_name = 'out0300'
  #goodArguments.output_grp_name = 'out0458'
  
  povray_cells(goodArguments)
  
  #plot_vessels(goodArguments)