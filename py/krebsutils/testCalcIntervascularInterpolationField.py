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
import os, sys
import numpy as np
import h5py
import extensions

from matplotlib import pyplot

import krebsutils
import myutils
import mpl_utils

# blurring is used to actually compute the local average with the help of a box filter.
SHOW_BLURRED_IMAGES = False

if SHOW_BLURRED_IMAGES:
  import scipy.signal as sig
import scipy.ndimage as ndimage

if __name__ == '__main__':
  filename = sys.argv[1]
  grouppath = sys.argv[2]
  
  with h5py.File(filename, 'r') as file:
    vesselgroup = file[grouppath]
    ldvessels = krebsutils.read_lattice_data_from_hdf(vesselgroup['lattice'])
    wbbox     = ldvessels.worldBox
    
    graph = krebsutils.read_vesselgraph(vesselgroup, ['position', 'radius', 'flags', 'pressure'])
    graph = graph.get_filtered(myutils.bbitwise_and(graph['flags'], krebsutils.CIRCULATED))
    
  fieldld = krebsutils.SetupFieldLattice(wbbox, 3, 50., 0.)
  print fieldld
  
  edgevalues = graph['pressure']
  edgevalues = edgevalues[graph.edgelist]
  
  thefield = krebsutils.CalcIntervascularInterpolationField(graph.edgelist, graph['radius'], graph['position'], edgevalues, fieldld, 1.)    
  volfraction = krebsutils.make_vessel_volume_fraction_field(graph['position'], graph.edgelist, graph['radius'], fieldld, 5)
  thegrad = ndimage.gaussian_filter(thefield, 1.0, 0, mode='nearest')
  gradfield_ = krebsutils.field_gradient(thegrad)
  thegrad = np.sum([np.square(g) for g in gradfield_], axis = 0)
  thegrad = np.sqrt(thegrad)
  del gradfield_

  if SHOW_BLURRED_IMAGES:
    kernel = (1.0/64.0)*np.ones((4,4,4), dtype = np.float32)
    thefield = sig.fftconvolve(thefield,kernel)
    volfraction = sig.fftconvolve(volfraction,kernel)
    thegrad =  sig.fftconvolve(thegrad,kernel)

  z = fieldld.shape[2]/2

  fig, axes = pyplot.subplots(2,2)
  ax = axes[0,0]
  plt = ax.imshow(thefield[:,:,z])
  ax.set(title = 'pressure')
  divider = mpl_utils.make_axes_locatable(ax)
  cax = divider.append_axes("right", size = "5%", pad = 0.05)
  fig.colorbar(plt, cax = cax)
  
  ax = axes[0, 1]
  plt = ax.imshow(volfraction[:,:,z])
  ax.set(title = 'volume fraction')
  divider = mpl_utils.make_axes_locatable(ax)
  cax = divider.append_axes("right", size = "5%", pad = 0.05)
  fig.colorbar(plt, cax = cax)  
  
  ax = axes[1,0]
  plt = ax.imshow(thegrad[:,:,z])
  ax.set(title = 'gradient magnitude')
  divider = mpl_utils.make_axes_locatable(ax)
  cax = divider.append_axes("right", size = "5%", pad = 0.05)
  fig.colorbar(plt, cax = cax)    

  del edgevalues, thefield, volfraction, thegrad

  pyplot.show()