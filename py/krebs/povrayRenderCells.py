#!/usr/bin/env python2
# -*- coding: utf-8 -*-
'''
This file is part of tumorcode project.
(http://www.uni-saarland.de/fak7/rieger/homepage/research/tumor/tumor.html)

Copyright (C) 2018  Thierry Fredrich

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
a nice rendering command
submitPovrayRender ../../fakeTumMTS-default-typeI-sample00-vbl_3486978.h5 out0* --cam=pie_only_cells --cells --noLabel --background=black --fontcolor=white -f --numTreads=4 --cellsProperty=o2 --cellsColorLimits 6 24
pH_ex 5.8 7.5

@artery
python2 /localdisk/thierry/tc_install/py/krebsjobs/submitPovrayRender.py fakeTumMTS-default-typeI-sample00-vbl_3484331.h5 out0540 --clip_box 240 400 -310 800 800 800 --cam_distance_multiplier=0.8 -f --cells --tumor_clip 240 400 -310
pwd: /localdisk/thierry/output/choose_position_2/artery/finished
python2 /localdisk/thierry/tc_install/py/krebsjobs/submitPovrayRender.py fakeTumMTS-default-typeI-sample00-vbl_3484331.h5 out0340 --clip_ball 240 400 -310 300 --cam_distance_multiplier=1.0 -f --cells --tumor_clip 240 400 -310
NOTE: something with the cam_distance_multiplier is still going wrong.
ffmpeg -r 20 -f image2 -s 1920x1080 -i "fakeTum-gero_3d_8mm-typeI-sample00-gero_3month_to_5mmd-out0%03d.png_slice.png" -vcodec libx264 -crf 25 -pix_fmt yuv420p growth.mp4
'''

if __name__ == '__main__':
  import os.path, sys
  sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'..'))


import os,sys
from os.path import basename, splitext

import krebsutils
import numpy as np
import povrayRenderVessels
#from povrayRenderVessels import *
import povrayEasy
import myutils
import math
import copy
import matplotlib

def paraview_color_space(scheme='jet'):
  colors = eval(scheme);
  N = len(colors);
  outfile_name = scheme + '.xml'
  with open(outfile_name, 'w') as f:
    f.write('<ColorMap name="%s" space="HSV">\n'%scheme)
    for i in range(N):
      #x = [(i-1)/(N-1); colors(i,:)'];
      f.write('  <Point x="%f" o="1" r="%f" g="%f" b="%f"/>\n' % x)
    f.write('</ColorMap>');

def convert_to_mmHg(o2_in_pg, cell_radii):
  cell_o2_concentration = o2_in_pg/ (4/float(3)* np.pi*np.power(cell_radii,3))
  #cell_o2_concentration = cell_o2_mass
  volume_o2_ml = cell_o2_concentration/(1.429*1e9)
  #volume_o2_ml = cell_o2_mass/1.429
  ''' o2 density 1.429 g/L --> 1.429*10^9 pg/ml
      1cm^3 = 10^12 (mum^3)
  '''
  solubility = 3.1e-3 #ml O2/cm^3 mmHg
  #solubility = 2.8e-3 #ml O2/cm^3 mmHg
  #solubility = 1.1e-4 #ml O2/cm^3 mmHg
  solubility = solubility*1e-12 #ml O2/mum^3 mmHg
  #volume_density = 1.429e9 #pg/ml
  #x = cell_o2_concentration/volume_density # ml / mum^3
  
  #cell_po2 = x/solubility
  
  return volume_o2_ml/solubility

#def createColormapForCells(data_color):
#  ''' convert data to color here'''
#  data_color = data_color
#  print(data_color[0:100])
#  cNorm = matplotlib.colors.Normalize(vmin=np.min(data_color), vmax=np.max(data_color))
#  scalar_map = matplotlib.cm.ScalarMappable(norm=cNorm, cmap='terrain')
#  return scalar_map
#class nlcmap(matplotlib.colors.LinearSegmentedColormap):
#    """A nonlinear colormap"""
#
#    name = 'nlcmap'
#
#    def __init__(self, cmap, levels):
#        self.cmap = cmap
#        self.monochrome = self.cmap.monochrome
#        self.levels = np.asarray(levels, dtype='float64')
#        self._x = self.levels/ self.levels.max()
#        self.levmax = self.levels.max()
#        self.levmin = self.levels.min()
#        self._y = np.linspace(self.levmin, self.levmax, len(self.levels))
#
#    def __call__(self, xi, alpha=1.0, **kw):
#        yi = np.interp(xi, self._x, self._y)
#        return self.cmap(yi/self.levmax, alpha)
      
def addVBLCells(epv, quantity , cell_hdf_group, options):
    print('path to vbl group is: %s' % cell_hdf_group.name)
    print('adding quantity: %s' % quantity)
    
    position_of_cells= np.asarray(cell_hdf_group['cell_center_pos'])
    radii_of_cells= np.asarray(cell_hdf_group['cell_radii'])
    
    ''' convert data to color here'''
    data_of_cells = np.asarray(cell_hdf_group[quantity])
    if quantity == 'o2':
      data_of_cells= convert_to_mmHg(data_of_cells, radii_of_cells )

    cNorm = matplotlib.colors.Normalize(vmin=np.min(data_of_cells[:,0]), vmax=np.max(data_of_cells[:,0]))
    if quantity == 'o2':
      cells_cm = matplotlib.cm.ScalarMappable(norm=cNorm, cmap='terrain')
    elif quantity == 'cell_phase':
      colors = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]  # R -> G -> B
      cmap_name = 'cell_phase_colors'
      cells_cm = matplotlib.colors.LinearSegmentedColormap.from_list(cmap_name, colors, N=len(colors))
      levels = [1.0, 4.0, 6.0]
      cells_cm = matplotlib.cm.ScalarMappable(norm=cNorm, cmap=cells_cm)
    else:
      cells_cm = matplotlib.cm.ScalarMappable(norm=cNorm, cmap='jet')
    #cells_cm = createColormapForCells(data_of_cells[:,0])
   
    
    
    #data_of_cells= np.asarray(cell_hdf_group[quantity])
    data_color = data_of_cells[:,0]
    
#    print(data_color[0:100])
#    cNorm = matplotlib.colors.Normalize(vmin=np.min(data_color), vmax=np.max(data_color))
#    scalar_map = matplotlib.cm.ScalarMappable(norm=cNorm, cmap='terrain')
#    scalar_map.get_cmap()
    if quantity == 'cell_phase':
      the_colors = cells_cm.to_rgba(data_color)
      rgb_colors_of_cells = copy.deepcopy(the_colors)
      for (i, data) in enumerate(data_color):
        if data == 1 or data == 2:
          #G1m, G1p
          #green
          rgb_colors_of_cells[i] = (0,1,0,0)
        if data == 3:
          #S
          #yellow
          rgb_colors_of_cells[i] = (1,1,0,0)
        if data == 4 or data ==5:
          #G2, M
          #yellow
          rgb_colors_of_cells[i] = (1,0,0,0)
        if data == 6:
          #dead
          if options.background == 1.0:
            rgb_colors_of_cells[i] = (0.,0.,0.,0.)
          elif options.background == 0.0:
            rgb_colors_of_cells[i] = (0.75,0.75,0.75,0)
          else:
            inv_val = 1.0 - options.background
            rgb_colors_of_cells[i] = (inv_val,inv_val,inv_val,0)
            
      rgb_colors_of_cells = rgb_colors_of_cells[:,:3]
      wbbox = options.wbbox
      trafo = povrayEasy.calc_centering_normalization_trafo(wbbox)
    
      epv.addVBLCells(trafo, position_of_cells, radii_of_cells, rgb_colors_of_cells, options)
    
      return None
      
    else:
      print('cells_cm.get_clim()')
      print(cells_cm.get_clim())
      print('data_color.shape')
      print(data_color.shape)
      data_color = cells_cm.to_rgba(data_color)
      print(data_color.shape)
      data_color=data_color[:,:3]
      print(data_color.shape)
      rgb_colors_of_cells=data_color
    
      wbbox = options.wbbox
      trafo = povrayEasy.calc_centering_normalization_trafo(wbbox)
      
      epv.addVBLCells(trafo, position_of_cells, radii_of_cells, rgb_colors_of_cells, options)
      
      return cells_cm
    
#    ld = krebsutils.read_lattice_data_from_hdf_by_filename(str(tumorgroup.file.filename), str(tumorgroup['conc'].attrs['LATTICE_PATH']))
#    ld = transform_ld(trafo, ld)
#
#    ds_necro    = np.asarray(tumorgroup['necro'])
#    data        = np.clip(np.asarray(tumorgroup['conc']), 0., 1.) # - np.asanyarray(tumorgroup['necro']), 0., 1.)
#
#    ds_levelset = -np.minimum(np.asarray(tumorgroup['ls']), 0.4 - ds_necro)
#    ds_levelset = krebsutils.distancemap(ds_levelset) * ld.scale
#
##    import matplotlib
##    matplotlib.use('Qt4Agg')
##    import matplotlib.pyplot as pyplot
##    pyplot.imshow(ds_levelset[:,:,8])
##    pyplot.contour(ds_levelset[:,:,8],[0.])
##    pyplot.show()
#    if 'tumor_clip' in options:
#      clip = clipFactory(options.tumor_clip)
#    else:
#      clip = clipFactory('None')
#
#    voldata_ls    = epv.declareVolumeData(ds_levelset, ld.GetWorldBox())
#    voldata_cells = epv.declareVolumeData(data, ld.GetWorldBox())
#
#    value_bounds = voldata_cells.value_bounds
#    style = """
#      texture {
#        pigment {
#          function { %f + %f*%s(x,y,z) }
#          color_map {
#            [0.0  color <0.3,0,0>]
#            [0.5  color <1,0.8, 0.3>]
#            [0.8  color <1,1,0.1>]
#          }
#        }
#        finish { 
#          specular 0.3
#        }
#      }""" % (value_bounds[0], (value_bounds[1]-value_bounds[0]), voldata_cells.name)
#    #style = " texture { pigment { color rgb<1,0.8,0.3> }  finish { specular 0.3 }}"
#    epv.addIsosurface(voldata_ls, 0., lambda : style, clip, style)

if __name__ == '__main__':
  paraview_color_space()