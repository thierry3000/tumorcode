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
#!/usr/bin/env python2
# -*- coding: utf-8 -*-
from os.path import basename, dirname, join, splitext
import os,sys
#if __name__ == '__main__':
#  import os.path, sys
#  sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','..'))
''' import of this must come in front of import of detailedo2 libs because some required initialization stuff on the c++ side (see mainboost.cpp)
    e.g. mulithreading initializer,...
    '''
import qsub
import h5files
import numpy as np
from copy import deepcopy
from dicttoinfo import Vec

import krebsutils 
if sys.flags.debug:
  tumor_cpp = __import__('libtumors_d', globals(), locals())
else:
  tumor_cpp = __import__('libtumors_', globals(), locals())
  
def vess_size_to_tum_size(fn, tum_lattice_const):
  vesselgroup = h5files.open(fn, 'r')['/vessels']
  if __debug__:
    print("vesselgroup.id.id: %i" % vesselgroup.id.id)
  ld = krebsutils.read_lattice_data_from_hdf(krebsutils.find_lattice_group_(vesselgroup))
  ld_vessels = krebsutils.read_lattice_data_from_hdf(vesselgroup['lattice'])
  bbox = ld_vessels.worldBox
  boxLengths = ((bbox[1]-bbox[0]),(bbox[3]-bbox[2]),(bbox[5]-bbox[4]))
  number_of_numerical_grid_points = np.ceil((np.asarray(boxLengths)+.5)/tum_lattice_const)
  for (i,x) in enumerate(number_of_numerical_grid_points):
    if x==1:
      number_of_numerical_grid_points[i] = 2
  
  del vesselgroup
  l = ld.shape
  s = int(ld.scale)
  """
  mapping from: latticesize initial vesselsnetwork 
  to numical lattice of tissue lattice
  first number: size x
  second number: size z
  third number: lattice constant
  """
  vs = {
    (9,1,130) : (10,10), #div runtime test
    (100,6,80) : (228, 10), # 100x5.80
    (100,62,80) : (196, 128), #100x50.80
    (200,6,80) : (384, 10), #200x5.80
    (60,75,80) : (128, 128), #60x60.80
    # now hierarchical builds
    (167,17,80) : (300, 34),
    (103,121,80) : (150, 150),
    # with large lattice spacing
    (59, 69, 150) : (160, 160),
    (63, 65, 150) : (160, 160),
    #(103,121,80) : (200, 200), # a special extra large configuration
    # 7x7 * 2^3 hierarchical biulds
    (127, 129, 80) : (150, 150),  # actual size could be 270**3
    # twodimensional system
    (160, 1, 80) : (350, 1),
    # quasi twodimensional
    (163, 9, 80) : (350, 24),
    (59, 9, 150) : (250, 32),
    # small system
    (31, 5, 150) : (80, 16),
    # 16 cm flat
    (111, 9, 150) : (267, 32),
    # 8mm P10 3d 
    #(55,61, 130) : (80,90),
    #(55,61, 130) : (155,155),
    (55,61, 130) : (100,100),
    # 8mm P10 q2d
    (59,5, 130) : (100, 15),
    (47,9, 130) : (50, 20),
    # for testing 3d_mini_mini
    (14, 18, 130) : (40, 40),
    # trastuzumab calibration growth
    (53, 65, 80) : (100, 100),
    # 5mm P11 q2d, 5mm x 5mm x 0.3mm
    (59,5, 90) : (int(4000./tum_lattice_const), int(300./tum_lattice_const+0.5)),
    # swine h1
    (17,19, 75) : (int(1200./tum_lattice_const), int(1200./tum_lattice_const)),
    # swine h2
    (35,37, 75) : (int(1400./tum_lattice_const), int(1400./tum_lattice_const)),
    # swine h2 -big (67, 77, 75
    (67,77, 75) : (int(5000./tum_lattice_const), int(5000./tum_lattice_const)),
    # 4mm x 0.3 P? ??
    (51,57, 88) : (int(4000./tum_lattice_const), int(4000./tum_lattice_const)),
    #the vessBigTras
    (119, 141, 88): (int(10000./tum_lattice_const), int(10000./tum_lattice_const)),
  }#[(l[0],l[2],s)]
  theKey = (l[0],l[2],s)
  if theKey in vs:
    return (vs[theKey][0],vs[theKey][0],vs[theKey][1])
  else:
    print('Warning: Guessed tum grid size')
    return number_of_numerical_grid_points.astype(int)
  #return (vs[0], vs[0], vs[1])
def set_lattice_size(c, vesselfilename):
  # set the lattice size based on translation table for vessel data sets
  c = deepcopy(c)
  sx, sy, sz = vess_size_to_tum_size(vesselfilename, c['lattice_scale'])
  #print('sx: %i, sy: %i, sz: %i')% (sx,sy,sz)
  if sz * 3 < sx:
    c['tumor']['init_shape'] = 'cylinder'
  else:
    c['tumor']['init_shape'] = 'sphere'
  c['lattice_size'] = Vec((sx,sy,sz))
  return c

def run_faketum_mts(configstr):
  if qsub.is_client:
    qsub.printClientInfo()
  return tumor_cpp.run_faketum_mts_(configstr)

def run_faketum(configstr):
  if qsub.is_client:
    qsub.printClientInfo()
  return tumor_cpp.run_faketum_(configstr)

def run_bulktissue_no_vessels(configstr):
  if qsub.is_client:
    qsub.printClientInfo()
  return tumor_cpp.run_bulktissue_no_vessels_(configstr)  

def run_bulktissue_w_vessels(configstr):
  if qsub.is_client:
    qsub.printClientInfo()
  return tumor_cpp.run_bulktissue_with_vessels_(configstr)
