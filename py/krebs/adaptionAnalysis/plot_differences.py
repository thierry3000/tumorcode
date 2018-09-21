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
Created on Wed Aug  5 12:31:50 2015

@author: thierry

"""
import h5py
import matplotlib.pyplot as plt
import itertools
import numpy as np
import os
import sys
from matplotlib.backends.backend_pdf import PdfPages
if __name__ == '__main__':
  sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../..'))
import myutils

def plot_differences(grp1,grp2, quantity,pp,interactive_output=True):
  if quantity == 'pressure':
    quantity_1 = grp1['nodes/pressure']
    quantity_2 = grp2['nodes/pressure']
    fig2 = plt.figure()
    plt.title('Differences of quantity: '+quantity)
    plt.plot(quantity_1,'*r')
    plt.plot(quantity_2,'*b')
    plt.legend(['before','after'])
    plt.grid()
    
  else:
    quantity_1 = grp1['edges/'+quantity]
    quantity_2 = grp2['edges/'+quantity]
    fig2 = plt.figure()
    plt.title('Differences of quantity: '+quantity)
    plt.semilogy(quantity_1,'*r')
    plt.semilogy(quantity_2,'*b')
    plt.legend(['before','after'])
    plt.grid()

#    fig3 = plt.subplot(2,1,2)
#    foobar = []
#    for aValue, bValue in itertools.izip(quantity_1,quantity_2):
#        foobar.append((float(aValue)-float(bValue))*60/1000000)
#    plt.plot(np.abs(foobar),'*')
  pp.savefig()
  if interactive_output:
    plt.show()
  plt.close()
  
def plot_hists(grp1,grp2, quantity,pp,interactive_output=True):
  no_bins=30
  fig2 = plt.figure()
  fig2.suptitle(quantity)
  ax1 = fig2.add_subplot(211)
  ax2 = fig2.add_subplot(212)
  if quantity == 'pressure':
    quantity_1 = grp1['nodes/pressure']
    quantity_2 = grp2['nodes/pressure']
    ax1.hist(quantity_1,no_bins)
    ax1.set_title('before')
    ax2.hist(quantity_2,no_bins)
    ax2.set_title('after')
    
    
  else:
    quantity_1 = grp1['edges/'+quantity]
    quantity_2 = grp2['edges/'+quantity]
    ax1.hist(quantity_1,no_bins)
    ax1.set_title('before')
    ax2.hist(quantity_2,no_bins)
    ax2.set_title('after')

#    fig3 = plt.subplot(2,1,2)
#    foobar = []
#    for aValue, bValue in itertools.izip(quantity_1,quantity_2):
#        foobar.append((float(aValue)-float(bValue))*60/1000000)
#    plt.plot(np.abs(foobar),'*')
  pp.savefig()
  if interactive_output:
    plt.show()
  plt.close()
    
if __name__ == '__main__':
  filenames = sys.argv[1:]
  #filename= 'mesentry_secomb_adption_p_mesentry_subset.h5'
  for fn in filenames:
    f = h5py.File(fn)
    vesselgrp_before = f['/adaption/recomputed']
    vesselgrp_after = f['/adaption/vessels_after_adaption']
    common_filename = os.path.splitext(os.path.basename(fn))[0]
    pp = PdfPages('diff_summary' + '_' + common_filename + '.pdf')
    for quan in 'radius flow pressure'.split():    
      
      #plot_differences(vesselgrp_before,vesselgrp_after,quan,pp, False)
      plot_hists(vesselgrp_before,vesselgrp_after,quan,pp, False)
      
    f.close
    pp.close()