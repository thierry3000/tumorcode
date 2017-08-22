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
'''this file is still buggy

However it is intended to be a modern interface
to the older IFF data published by
1.Welter, M. & Rieger, H. Interstitial Fluid Flow and Drug Delivery in Vascularized Tumors: A Computational Model. PLoS ONE 8, e70395 (2013).
'''
import os, sys
from os.path import join, basename, dirname, splitext
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'..'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../..'))
import myutils
import h5py
import numpy as np
import sys
import h5files
import matplotlib.pyplot as pyplot

from os.path import basename

import mpl_utils

from plotIff import DataTissue, DataGlobalIff, DataRadialIff, DataTumorBloodFlow, averaged

if __name__ == '__main__':
  
    
  if 0:#this is read in  
    filenames = sys.argv[1:]
    files = [ h5py.File(fn, 'r') for fn in filenames ]
    fmeasure = h5files.open('analyzeIFF.h5', 'a')
    dataman = myutils.DataManager(100, [ DataTissue(), DataGlobalIff(), DataRadialIff(), DataTumorBloodFlow()])
    #data = dataman('iff_radial', files[0], 'vs_dr', 'iff_pressure') 
    xlim = -1.5, 1.0

    array_of_plateau = []
    array_of_v_max = []
    
    for f in files:
      print('Running file: %s' % f.filename)
      iff_pressure = dataman('iff_radial', f, 'vs_dr', 'iff_pressure')
      iff_velocity_mag = dataman('iff_radial', f, 'vs_dr', 'iff_velocity_mag')
      bins =(1./1000.) * dataman('iff_radial', f, 'vs_dr', 'bins')
      mask = np.logical_and(bins < xlim[1], bins >= xlim[0])
      nice_iff_pressure = iff_pressure.avg[mask]
      nice_iff_velocity_mag = iff_velocity_mag.avg[mask]
      
      plateau = 1/4.*(np.sum(nice_iff_pressure[0:4]))
      max_velocity = np.max(nice_iff_velocity_mag)
      print('max IFP: %f, max IFV: %f' % (plateau, max_velocity))
      array_of_plateau.append(plateau)
      array_of_v_max.append(max_velocity)
      
      common_base=basename(f.filename)
      fmeasure.create_group(common_base)
      fmeasure[common_base].attrs.create('plateau',plateau)
      fmeasure[common_base].attrs.create('max_velocity', max_velocity)
      f.close()
    
    fmeasure.create_dataset('all_plateaus', data=array_of_plateau)
    fmeasure.create_dataset('all_v_max', data=array_of_v_max)

  if 1:
    fmeasure = h5files.open('/localdisk/thierry/review_data/analyzeIFF.h5', 'a')
    
    a = np.asarray(fmeasure['/all_plateaus'])
    b = np.asarray(fmeasure['/all_v_max'])
    fig1 = pyplot.figure()
    ax = pyplot.gca()
    ax.scatter(a,b)
    pyplot.xlabel('IFP/$kPa$')
    ax.set_yscale('log')
    pyplot.ylabel(r'max IFV/ $\frac{\mu m}{s}$')
    pyplot.grid()
    pyplot.savefig('IFP_vs_max_IFV.png')    
    pyplot.show()
    
    
    if 0:
      radialgetter = lambda f,name: dataman('iff_radial', f, 'vs_dr', name)
      radial = averaged(radialgetter,files,
                          ['iff_pressure', 'ivp', 'iff_velocity_mag'])
  
      bins = (1./1000.) * dataman('iff_radial', files[0], 'vs_dr', 'bins')
      xlim = -1.5, 1.0
      mask = np.logical_and(bins < xlim[1], bins >= xlim[0])
      
      def plot(name, **kwargs):
        f = kwargs.pop('value_prefactor', 1.)
        avg, std = radial[name]
        avg, std = avg[mask], std[mask]
        kwargs.update(every = 5, marker = kwargs.pop('marker', 's'), color = kwargs.pop('color', 'k'))
        ret = mpl_utils.errorbar(ax, bins[mask], f*avg, yerr = f*std, **kwargs)
        return ret
      name='iff_pressure'
      avg, std = radial[name]
      avg, std = avg[mask], std[mask]

      
      fig, axes = pyplot.subplots(3, 1)
      axes = axes.ravel()
      mpl_utils.subplots_adjust_abs(fig, left = mpl_utils.mm_to_inch*20,
                                    right = -mpl_utils.mm_to_inch*10,
                                    top  = -mpl_utils.mm_to_inch*5,
                                    bottom = mpl_utils.mm_to_inch*10,
                                    hspace = mpl_utils.mm_to_inch*30,)
    
      ax = axes[0]
      ax.set(xticklabels=[])
      ax.set(ylabel = r'[kPa]')
      plot('iff_pressure', label = r'$p_i$', marker='o', color = 'r')
      #plot('ivp', label = '$p_v$', marker = 's', color = 'k')
      #plot('ivp_minus_ifp', label = '$p_v - p_i$', marker='>', color = 'b')
      #end pressure
      ax = axes[1]
      ax.set(ylabel = ur'[\u03BCm/s]') #, xlabel = r'$\theta$ [mm]', title = 'velocity')
      ax.set(xticklabels=[])
  #    plot('iff_velocity_out', label = r'$v_{||}$', marker = 'o', color = 'r')
  #    plot('iff_velocity_mag', label = r'$|v|$', marker = 's', color = 'b')
    
      ax.legend()    
      gridcolor = (0.7,0.7,0.7)
      for ax in axes:
        ax.set_xlim(*xlim)
        ax.grid(linestyle=':', linewidth=0.5, color=gridcolor)
        mpl_utils.add_crosshair(ax, (0,0), color = gridcolor)
  
      
      pyplot.show()