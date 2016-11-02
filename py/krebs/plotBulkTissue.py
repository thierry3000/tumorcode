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
if __name__ == '__main__':
  import os.path, sys
  sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'..'))

import numpy as np
import h5py
import extensions
from mystruct import Struct
import myutils
import os
import mpl_utils
import collections

#import krebsutils

import matplotlib

##############################################################################
### mostly helper functions for plotting
##############################################################################

def commonOutputName(filenames):
  dst = os.getcwd()
  fnmeasure = os.path.commonprefix(filenames)
  fnmeasure = os.path.basename(fnmeasure)
  if fnmeasure.endswith('.h5'): fnmeasure = fnmeasure[:-len('.h5')]
  fnmeasure = os.path.join(dst, fnmeasure)
  return fnmeasure


def imshow(axes, a, ld, **kwargs):
  '''
    uses lattice data;
    add new keywords:
      crange = one of ((min,max), 'zero-centered', or None==auto-computed)
      vscale = value scaling. crange refers to unscaled data values.
  '''
  crange = kwargs.pop('crange', None)
  if crange is None:
    vmin = kwargs.pop('vmin', None)
    vmax = kwargs.pop('vmax', None)
    if vmin is None: vmin = a.min()
    if vmax is None: vmax = a.max()
    crange = (vmin, vmax)
  elif crange == 'zero-centered':
    q = np.abs(a).max()
    crange = (-q,q)
  vscale = kwargs.pop('vscale', None)
  if vscale is not None:
    a = a*vscale
  else:
    vscale = 1.
  kwargs['vmin'] = vscale*crange[0]
  kwargs['vmax'] = vscale*crange[1]
  return axes.imshow(
    a,
    extent = ld.worldBox[:4] * kwargs.pop('worldscaling',1.),
    origin = 'lower',
    interpolation = kwargs.pop('interpolation','nearest'),
    **kwargs)



def imslice(a):
  if len(a.shape)==3:
    a = a[:,:,a.shape[2]/2]
  return np.asarray(a).transpose()


def contour(axes, a, ld, **kwargs):
  return axes.contour(
    a,
    extent = ld.worldBox[:4] * kwargs.pop('worldscaling',1.),
    linewidths = kwargs.pop('linewidth', 0.25),
    colors = kwargs.pop('colors', 'w'),
    **kwargs)


def colorbar(fig, axes, theplot, height = "50%", **kwargs):
#  axins = mpl_utils.inset_axes(
#    axes,
#    width = "3%",
#    height = height,
#    loc=2
#  )
  divider = mpl_utils.make_axes_locatable(axes)
  cax = divider.append_axes("right", size = "5%", pad = 0.05)
  return fig.colorbar(theplot, cax = cax, **kwargs), cax



def with_cb(func):
  '''Calls the wrapped function and generates a colorbar for the
     plot (axes Object, which has to be the first argument of the function)'''
  def wrapper(ax, *args, **kwargs):
    ret = func(ax, *args, **kwargs)
    colorbar(ax.get_figure(), ax, ret)
    return ret
  return wrapper


def with_contour_factory(levelset, ld, **contour_kwargs):
  def decorator(func):
    def wrapper(ax, *args, **kwargs):
      ret = func(ax, *args, **kwargs)
      contour(ax, levelset, ld, levels = [0.], **contour_kwargs)
      return ret
    return wrapper
  return decorator



def logmap(arr, Nmin = -2, Nmax = 2):
  bb = np.power(10., np.arange(Nmin, Nmax))
  N = len(bb)
  def f(x):
    if x < bb[0]:
      return x/bb[0]/N
    for i in xrange(N-1):
      b1 = bb[i+1]
      if x > b1: continue
      b0 = bb[i]
      return ((x - b0)/(b1-b0) + i+1)/N
    return 1.

  uf = np.frompyfunc(f, 1, 1)
  arr = np.asarray(uf(arr), dtype=np.float)
  return arr
logticks = [ '0', '$10^{-2}$', '$10^{-1}$', '$10^{0}$', '$10^1$']

if 0: # my mapping vs standard logarithmic mapping
  x = np.arange(1.e-4, 10., 1.e-5)
  y = logmap(x)
  pyplot.plot(x,y)
  pyplot.plot(x, (np.log10(x)-1.)/(2.-1)*0.25+1.)
  pyplot.setp(pyplot.gca(), xscale = 'log')
  pyplot.show()

def getSortedTimeGroupNames(f):
  gg = myutils.getTimeSortedGroups(f['.'], "out")
  return [ g.name for g in gg ]

from quantities import Prettyfier
