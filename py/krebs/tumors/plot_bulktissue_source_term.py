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

"""
      Algorithmic design
      This plot is not for evaluating data, but for 
      visualize some algorithmic aspects.
"""
if __name__ == '__main__':
  import os.path, sys
  sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'..'))

import os,sys
import time
import h5py
import numpy as np
import extensions # for asarray with h5py support
import krebsutils
import math
from myutils import f2s

from dicttoinfo import dicttoinfo, Vec

import matplotlib
#matplotlib.use('Qt4Agg')
import matplotlib.pyplot as pyplot
import matplotlib.patches
from mpl_toolkits.mplot3d import Axes3D
import mpl_utils

rc = matplotlib.rc
rc('figure', figsize=(5, 5))
rc('font', size = 8.)
rc('axes', titlesize = 10., labelsize = 8.)
rc('pdf', compression = 6, fonttype = 42)


def generate_data(params, default_args, **kwargs):
  """
    tensor product style sampling of source term
  """
  arg_dict = default_args.copy()
  arg_dict.update(kwargs)
  param_string = dicttoinfo(params)
  args = tuple(np.atleast_1d(np.asarray(arg_dict[q])) for q in 'theta ncells nother nnecro o2'.split())
  shape = tuple(len(a) for a in args)
  arglist = []
  for ii in np.ndindex(shape):
    arglist.append(tuple(args[j][i] for j,i in enumerate(ii)))
  q = krebsutils.calcBulkTissueSourceTerm(arglist, param_string)
  for i, ii in enumerate(np.ndindex(shape)):
    q[i] /= args[1][ii[1]]
  return q


if __name__ == '__main__':
  
  #print args
  #d = generate_data(*args, params=dict())
  #print d
  
  with mpl_utils.PdfWriter('bulktissuesources_new.pdf') as pdfpages:

    params = dict(source_model_version=3,
                  o2_necro_threshold=0.1,
                  o2_prol_threshold=0.4, 
                  time_necrosis=48,
                  time_death_tumor=240000.,
                  use_necrotic_regions=False)  
    if 1:
      o2values = [ 0., 0.11, 0.041, 1. ]
      
      fig = pyplot.figure()
      
      plot = fig.add_subplot(111, xlabel='$\phi$',
                             title = r'proliferation rate coefficient $\Gamma_\phi(\phi,O_2) / \phi$')
      plot.grid(linestyle=':', linewidth=0.5, color=(0.7,0.7,0.7))
      plot.plot([0., 1.], [0., 0.], 'k-')
      
      # full o2, no tumor
      args = dict(ncells = np.linspace(0.01, 1., 200), nother=0.2, o2=1., nnecro=0., theta=0.)
      
      def mydata(**kwargs):
        return generate_data(params, args, **kwargs)

      plot.plot(args['ncells'], mydata(), 'g', label = r'$\theta = 0$ (normal)', linewidth=2.)
      # varying o2
      for o2 in o2values:
        q = mydata(o2=o2)
        plot.plot(args['ncells'], q, 'g--')
        
      # intermediate 1/2 tumor
      q = mydata(theta=0.5)
      plot.plot(args['ncells'], q, 'b', label = r'$\theta = 1/2$', linewidth=2.)
      
      # full o2, tumor
      q = mydata(theta=1.)
      plot.plot(args['ncells'], q, 'r', label = r'$\theta = 1$ (tumor)', linewidth=2.)
      # varying o2
      for o2 in o2values:
        q = mydata(theta=1., o2=o2)
        plot.plot(args['ncells'], q, 'r--')
        plot.text(0.1, mydata(o2=o2, ncells=0.1)[0]+0.001, r'$O_2 = %s$' % f2s(o2))
     
      # finished
      plot.legend()  

      levels = [ -1./48, -1./238, 0., 1./238, 1./48, 1./24 ]
      levellabels = [ '-1/48', '-1/240', '0', '1/240', '1/48', '1/24']

      plot.set(yticks = levels,
               yticklabels = levellabels)
      
      pdfpages.savefig(fig)
      
      # plot vs o2
      fig = pyplot.figure()
      plot = fig.add_subplot(111, xlabel='$O_2$',
                                  title = r'proliferation rate coefficient $\Gamma_\phi(\phi,O_2) / \phi$')
      plot.grid(linestyle=':', linewidth=0.5, color=(0.7,0.7,0.7))
      plot.plot([0., 1.], [0., 0.], 'k-')
      
      args = dict(ncells = 0., nother=0.2, o2=np.linspace(0., 1., 500), nnecro=0., theta=0.)
      
      def mydata(**kwargs):
        return generate_data(params, args, **kwargs)
      
      for ncells in [ 0.2, 0.8 ]:
        style_ls, style_lw = ('--',2.) if ncells==0.2 else (':',1.)
        
        q = mydata(ncells = ncells, theta=0.)
        plot.plot(args['o2'], q, style_ls+'g', label = r'$\theta = 0$ (normal)', linewidth=style_lw)  
        q = mydata(ncells = ncells, theta=0.5)
        plot.plot(args['o2'], q, style_ls+'b', label = r'$\theta = 1/2$', linewidth=style_lw)  
        q = mydata(ncells = ncells, theta=1.)
        plot.plot(args['o2'], q, style_ls+'r', label = r'$\theta = 1$ (tumor)', linewidth=style_lw)
      
      plot.set(yticks = levels,
               yticklabels = levellabels)
      pdfpages.savefig(fig)
    
    if 1:
      ncells_values = np.linspace(0.01, 1., 100)
      o2_values = np.linspace(0., 1., 100)
      data = generate_data(params, dict(o2=o2_values, ncells=ncells_values, nother=0., theta=0., nnecro=0.))
      X, Y = np.meshgrid(ncells_values, o2_values)
      Z = data.reshape(X.shape)
      fig = pyplot.figure()
      plot = fig.add_subplot(111, xlabel=r'$O_2$', ylabel=r'$\phi$', title=r'$\Gamma_\phi(\phi,O_2) / \phi$')
      p = plot.imshow(Z, cmap=matplotlib.cm.binary_r, origin='lower', extent=(ncells_values[0], ncells_values[-1], o2_values[0], o2_values[-1]))
      #fig.colorbar(p, shrink=0.5)
      p = plot.contour(X, Y, Z, 
                   levels = levels,
                   colors = [ ('r' if l<0 else ('g' if l > 0. else 'w')) for l in levels ]
                   )
      #fig.colorbar(p, shrink=0.5)
      pdfpages.savefig(fig)


    if 1:
      def make3d(nc, no2, **kwargs):
        #grid setup
        ncells_values = np.linspace(0.01, 1., nc)
        o2_values = np.linspace(0., 1., no2)
        data = generate_data(params, dict(o2=o2_values, ncells=ncells_values, nother=0., theta=0., nnecro=0.), **kwargs)
        X, Y = np.meshgrid(o2_values, ncells_values)
        Z = data.reshape(X.shape)
        return X, Y, Z

      fig = pyplot.figure()
      plot = fig.add_subplot(111, 
                             projection='3d', 
                             xlabel=r'$O_2$', ylabel=r'$\phi$', 
                             title=r'$\Gamma_\phi(\phi,O_2) / \phi$')

      X, Y, Z = make3d(50, 100)      

      plot.plot_wireframe(X, Y, Z, cstride = 20, rstride = 20, color='g')
#      plot.plot_surface(X, Y, Z,
#                        rstride=1, cstride=1, 
#                        cmap=matplotlib.cm.binary_r,
#                        linewidth = 0.,
#                        antialiased = False,
#                        )
#      levels = np.asarray([ -1./24, -1./48, -1./238, 0., 1./238, 1./48, 1./27 ])

      p = plot.contour(X, Y, Z, 
                       levels = levels,
                       colors = [ ('g' if l<-0.0001 else ('g' if l > 0.0001 else 'k')) for l in levels ]
                       )
      plot.set(zticks = levels,
            zticklabels = levellabels)
      plot.view_init(30, 180-45)
#      p = plot.contour(X, Y, Z, 
#                   levels = levels,
#                   colors = [ ('r' if l<0 else ('g' if l > 0. else 'w')) for l in levels ]
#                   )
#      fig.colorbar(p, shrink=0.5)
      pdfpages.savefig(fig)
