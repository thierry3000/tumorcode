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
from os.path import join, basename, dirname, splitext
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'..'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../..'))

import h5files
import numpy as np
import krebsutils as ku
import os
import matplotlib.pyplot as plt
#import adaption
import mpl_utils

def DoGetMurrayForSingleFileGENERAL(fn,pattern, pdfpages):
  afile = h5files.open(fn, 'r')
  vesselgrp = afile[pattern]
  print(vesselgrp)
  result_venous, result_a = ku.get_Murray(vesselgrp)
  #result = ku.get_Murray(vesselgrp, alpha)
  result_venous = removeZeros(result_venous)
  result_a = removeZeros(result_a)
 
  #x=np.arange(0,max_,0.1)
  #y=x
  #fig1, (ax1, ax2) = plt.subplots(2,1,sharex=True, sharey=True)
  fig1, (ax1, ax2) = plt.subplots(2,1)
  ### plot murrays law with exponent
  alphas = np.arange(2,6,0.5)
  max_ = np.max(result_venous[2,:])
  min_ = np.min(result_venous[0,:])
  x = np.arange(min_,max_, 0.1)
  font = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 4,
        'rotation': 30,
        'verticalalignment': 'bottom',
        'horizontalalignment': 'left',
        }
  for alpha in alphas:
    ax1.text(np.power(2 * x[-1]**alpha, 1/float(alpha)), x[-1], r'$\alpha = %0.1f $' % alpha, fontdict=font)
    ax1.plot(np.power(2 * x**alpha, 1/float(alpha)) , x, c='k', linewidth=0.1)
  
  ### plot data
  ax1.set_title('Murray for file:\n%s\n%i venous points, %i arterial points' % (os.path.basename(afile.filename),result_venous.shape[1], result_a.shape[1]), size=6)
  pointsize = 0.1
  venous = ax1.scatter(np.power(result_venous[0,:]**3+result_venous[1,:]**3,1/float(3)),result_venous[2,:],c='b', s=pointsize, label='bldk')
  arterial = ax1.scatter(np.power(result_a[0,:]**3+result_a[1,:]**3,1/float(3)),result_a[2,:],c='r',s=pointsize)
  ax1.set_xlabel(r'$\sqrt[3]{r_{daughter_1}^3 + r_{daughter_2}^3}$')
  ax1.set_ylabel(r'$r_{mother}$')
  ax1.legend((venous,arterial),('Venous', 'Arterial'),loc='lower right')
  ax1.grid()
  
  
  ### histogramm deviation of daughters
  hist, bin_edges = np.histogram(result_a[1,:]-result_a[0,:], bins=50) #note. they are sorted in on c++ side
  width = 0.7*(bin_edges[1]-bin_edges[0])
  centers = (bin_edges[:-1]+bin_edges[1:])/2
  ax2.bar(centers,hist/float(len(result_a[0,:])), align='center', width=width, color='red')
  ax2.set_xlabel(r'$\|r_{daughter 1}- r_{daughter 2}\|$')
  ax2.set_ylabel(r'$p$')
  fig1.tight_layout()
  pdfpages.savefig(fig1) 

def DoGetMurrayForSingleFile(fn,pattern):
    vesselgrp_before = fn['adaption/recomputed']
    result = ku.get_Murray(ku.find_lattice_group_(vesselgrp_before),vesselgrp_before)
    vesselgrp_after = fn['adaption/vessels_after_adaption']
    result_after = ku.get_Murray(ku.find_lattice_group_(vesselgrp_after),vesselgrp_after)
    #print result
    max_before = np.max(result[0,:])
    max_after = np.max(result_after[0,:])
    xmax = np.ceil(np.fmax(max_before,max_after))
    x=np.arange(0,xmax,0.1)
    y=x
    fig1, (ax1, ax2) = plt.subplots(2,1,sharex=True, sharey=True)

    #fig1.suptitle('this is the figure title', fontsize=12)
    ax1.set_title('Murray for \n%s\n%i points considered' % (os.path.basename(afile.filename),result.shape[1]))
    #plt.subplot(2,1,1)
    ax1.scatter(result[0,:],np.sqrt(result[1,:]**2+result[2,:]**2))
    ax1.set_xlabel('radius mother/ daughter vessel')
    ax1.set_ylabel(r"$\sqrt{r_a^2+r_b^2}$")
    #plt.title('Murray for \n%s\n%i points considered' % (os.path.basename(afile.filename),result.shape[1]))
    ax1.plot(x,y,'k')
    ax1.grid()
    ax1.legend(['murray','before_data'], loc=4)
#      plt.subplot(2,1,2)
    ax2.scatter(result_after[0,:],np.sqrt(result_after[1,:]**2+result_after[2,:]**2))
    #ax2.set_xlabel('radius mother/ daughter vessel')
    ax2.set_ylabel(r"$\sqrt{r_a^2+r_b^2}$")
    #plt.title('Murray for \n%s\n%i points considered' % (os.path.basename(afile.filename),result.shape[1]))
    ax2.plot(x,y,'k')
    ax2.legend(['murray','after_data'], loc=4)
    ax2.grid()
    plt.savefig('murray_%s.png' % os.path.basename(afile.filename))

def removeZeros(atwoarray):
    first_line = atwoarray[0,:]
    good_indeces_first_line = first_line>0
    second_line = atwoarray[1,:]
    good_indeces_second_line = second_line>0
    all_good_indeces = np.bitwise_and(good_indeces_first_line,good_indeces_second_line)
    return atwoarray[:,all_good_indeces]

def printMurray(dataman, f_measure, filenames, options, pdfpages):
    alpha = 3
    filenames = adaption.get_files_with_successful_adaption(filenames)
    files = [h5files.open(fn, 'r') for fn in filenames]

def DoGetMurray(filenames, pdfpages):
    alpha = 3
    filenames = adaption.get_files_with_successful_adaption(filenames)
    files = [h5files.open(fn, 'r') for fn in filenames]
    for afile in files:
      print("consider file: %s " % afile.filename)
      if('adaption/recomputed' in afile ):
        vesselgrp_before = afile['adaption/recomputed']
      if('vessels/recomputed_flow' in afile ):
        vesselgrp_before = afile['vessels/recomputed_flow']
      if(not 'vesselgrp_before' in locals()):
        print("bad")
      if ('results' in locals()):
        results = np.hstack((results, ku.get_Murray(vesselgrp_before, alpha)))
      else:
        results = ku.get_Murray(vesselgrp_before, alpha)
      vesselgrp_after = afile['adaption/vessels_after_adaption']
      if ('results_after' in locals()):
        results_after = np.hstack((results_after, ku.get_Murray(vesselgrp_after, alpha)))
      else:
        results_after = ku.get_Murray(vesselgrp_after, alpha)
    
    results = removeZeros(results)
    results_after = removeZeros(results_after)
      #print result
    max_before = np.amax(results)
    max_after = np.amax(results_after)
    xmax = np.ceil(np.fmax(max_before,max_after))
    x=np.arange(0,xmax,0.1)
    y=x
    fig1, (ax1, ax2) = plt.subplots(2,1,sharex=True, sharey=True)
    fig2, (ax3, ax4) = plt.subplots(2,1,sharex=True, sharey=True)
    #fig2 = plt.hist(results[0,:]-results[1,:],20)
    #fig1.suptitle('this is the figure title', fontsize=12)
    ax1.set_title('Murray for \n%s\n%i points considered' % (os.path.basename(afile.filename),results.shape[1]))
    #plt.subplot(2,1,1)
    ax1.scatter(results[0,:],results[1,:])
    #ax1.set_xlabel('daughter vessel')
    ax1.set_ylabel('mother')
    #plt.title('Murray for \n%s\n%i points considered' % (os.path.basename(afile.filename),result.shape[1]))
    ax1.plot(x,y,'k')
    ax1.grid()
    ax1.legend(['murray','before_data'], loc=4)
#      plt.subplot(2,1,2)
    ax2.scatter(results_after[0,:],results_after[1,:])
    ax2.set_xlabel(r"daughter $\sqrt{\sum_ir_i^2}$")
    #ax2.set_ylabel(r"$\sqrt{r_a^2+r_b^2}$")
    #plt.title('Murray for \n%s\n%i points considered' % (os.path.basename(afile.filename),result.shape[1]))
    ax2.plot(x,y,'k')
    ax2.legend(['murray','after_data'], loc=4)
    ax2.grid()
    pdfpages.savefig(fig1)
    #fig1.savefig('murray_%s.png' % os.path.basename(afile.filename))
    
    hist, bin_edges = np.histogram(results[0,:]-results[1,:], bins=50)
    width = 0.7*(bin_edges[1]-bin_edges[0])
    centers = (bin_edges[:-1]+bin_edges[1:])/2
    ax3.bar(centers,hist, align='center', width=width)
    
    hist_after, bin_edges = np.histogram(results_after[0,:]-results_after[1,:],bin_edges)
    ax4.bar(centers,hist_after, align='center', width=width)
    #fig2.savefig('dens_murray_%s.png' % os.path.basename(afile.filename))
    pdfpages.savefig(fig2)  
    #n,bins,patches = plt.hist(results[0,:]-results[1,:],20, normed=1)
    
    #ax_inset = fig1.add_axes([0.6,0.7,0.3,0.3])
    #ax3 = plt.axes([])
    #ax4    
    #n,bins,patches = plt.hist(results[0,:]-results[1,:],20, normed=1)
    
    #ax_inset2 = fig1.add_axes([0.0,0.3,0.3,0.3])
    #n2,bins2,patches2 = plt.hist(results_after[0,:]-results_after[1,:],20, normed=1)
    #ax3 = plt.hist(results[0,:]-results[1,:],20)
    #plt.grid()
    #plt.show()

if __name__ == "__main__":
  import optparse
  parser = optparse.OptionParser()
  parser.add_option("-O","--with-o2", dest="with_o2", help="look at detailed o2 data", default=False, action="store_true")
  options, args = parser.parse_args()
  filenames = args[:]
  
#either specify group, like in tumor simuations or take all from the ensemble
  if 0:    
    for afilename in filenames:
        print("got file: %s" % afilename)
    if not filenames:
        print("python2 got no filenames")    
    with mpl_utils.PdfWriter('murray' + '.pdf') as pdfpages:
      DoGetMurray(filenames,pdfpages)
  else:
    pattern =  args[-1]
    print(pattern)
    print(filenames[0])
    with mpl_utils.PdfWriter('murray' + '.pdf') as pdfpages:
      DoGetMurrayForSingleFileGENERAL(filenames[0],pattern, pdfpages)
    