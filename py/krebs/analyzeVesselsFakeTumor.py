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
import os,sys
from os.path import join, basename, dirname, splitext, normpath, abspath
sys.path.append(join(dirname(__file__),'..'))

from os.path import basename, dirname, join, splitext
import h5py
import h5files
import numpy as np
import extensions # for hdf5 support in np.asarray
import krebsutils as krebs
import myutils
import posixpath
import math
import collections
import matplotlib.pyplot as plt

from collections import defaultdict
from pprint import pprint

from analyzeGeneral import DataVesselLength, DataVesselTumor

import mpl_utils


import myutils

def create_mean_plot_length():
    my_res_file = h5files.open('analyzeVesselTumor.h5','r')
    count_hists = 0
    the_keys= my_res_file['/global'].keys()
    print(the_keys[0])
    the_density_in = np.zeros_like(my_res_file['/global/%s/length_dist/inside_density' %the_keys[0]])
    the_density_out = np.zeros_like(my_res_file['/global/%s/length_dist/outside_density' %the_keys[0]])
    a_bin_in = my_res_file['/global/%s/length_dist/inside_bins' %the_keys[0]]    
    a_bin_out = my_res_file['/global/%s/length_dist/outside_bins' %the_keys[0]]    
    for key in my_res_file['/global'].keys():
        print('working with %s' % key)
        count_hists=count_hists+1
        the_density_in = the_density_in+my_res_file['/global/%s/length_dist/inside_density' %key]
        the_density_out = the_density_out+my_res_file['/global/%s/length_dist/outside_density' %key]
    print('we have %i histograms' % count_hists)
    the_density_in = np.divide(the_density_in,count_hists)
    the_density_out = np.divide(the_density_out,count_hists)
    mu1 = np.sum(the_density_in*np.diff(a_bin_in)*(0.5*(a_bin_in[1:]+a_bin_in[:-1])))
    quad1 = np.sum(the_density_in*np.diff(a_bin_in)*np.square((0.5*(a_bin_in[1:]+a_bin_in[:-1]))))
    std1 = np.sqrt(quad1-mu1**2)    
    mu2 = np.sum(the_density_out*np.diff(a_bin_out)*(0.5*(a_bin_out[1:]+a_bin_out[:-1])))
    quad2 = np.sum(the_density_out*np.diff(a_bin_out)*np.square((0.5*(a_bin_out[1:]+a_bin_out[:-1]))))
    std2 = np.sqrt(quad2-mu2**2)    
    plt.title('Length distribution')    
    ax = plt.subplot(211)
    ax.bar(a_bin_out[:-1],the_density_out,np.diff(a_bin_out))
    ax.text(300,0.035,'$< l_{out} > = %f \pm %f$' % (mu2, std2))
    plt.legend(['Outside'])
    plt.ylabel('Propability')    
    ax = plt.subplot(212)    
    plt.bar(a_bin_in[:-1],the_density_in,np.diff(a_bin_in))
    ax.text(300,0.008,'$< l_{in} > = %f \pm %f$' % (mu1, std1))
    plt.xlabel('length of vessel segments/ $\mu m$')
    plt.ylabel('Propability')
    
    plt.legend(['Inside'])
    plt.savefig('Length_distribution.png')
    
    #plt.show()

def create_mean_plot_radii():
    my_res_file = h5files.open('analyzeVesselTumor.h5','r')
    count_hists = 0
    the_keys= my_res_file['/global'].keys()
    print(the_keys[0])
    no_bins=55
    the_density_in = np.zeros(no_bins)
    the_density_out = np.zeros(no_bins)
    for key in my_res_file['/global'].keys():
        print('working with %s' % key)
        count_hists=count_hists+1
        count_in = my_res_file['/global/%s/radii_within_tumor/radii_inside' %key]
        count_out = my_res_file['/global/%s/radii_within_tumor/radii_outside' %key]
        density_in, bin_in_array = np.histogram(count_in,no_bins, density=True)
        density_out, bin_out_array = np.histogram(count_out,no_bins, density=True)
        the_density_in = the_density_in+density_in
        the_density_out = the_density_out+density_out
    print('we have %i histograms' % count_hists)
    the_density_in = np.divide(the_density_in,count_hists)
    the_density_out = np.divide(the_density_out,count_hists)
    mu1 = np.sum(the_density_in*np.diff(bin_in_array)*(0.5*(bin_in_array[1:]+bin_in_array[:-1])))
    quad1 = np.sum(the_density_in*np.diff(bin_in_array)*np.square((0.5*(bin_in_array[1:]+bin_in_array[:-1]))))
    std1 = np.sqrt(quad1-mu1**2)    
    mu2 = np.sum(the_density_out*np.diff(bin_out_array)*(0.5*(bin_out_array[1:]+bin_out_array[:-1])))
    quad2 = np.sum(the_density_out*np.diff(bin_out_array)*np.square((0.5*(bin_out_array[1:]+bin_out_array[:-1]))))
    std2 = np.sqrt(quad2-mu2**2)    
    plt.title('Radii distribution')    
    ax = plt.subplot(211)
    ax.bar(bin_out_array[:-1],the_density_out,np.diff(bin_out_array))
    ax.text(30,0.2,'$< r_{out} > = %f \pm %f$' % (mu2, std2))
    plt.legend(['Outside'])
    plt.ylabel('Propability')    
    ax = plt.subplot(212)    
    plt.bar(bin_in_array[:-1],the_density_in,np.diff(bin_in_array))
    ax.text(20,0.3,'$< r_{in} > = %f \pm %f$' % (mu1, std1))
    plt.xlabel('radius of vessel segments/ $\mu m$')
    plt.ylabel('Propability')
    
    plt.legend(['Inside'])
    plt.savefig('Radii_distribution.png')
    
    plt.show()
    
def create_mean_plot_flows():
    my_res_file = h5files.open('analyzeVesselTumor.h5','r')
    count_hists = 0
    the_keys= my_res_file['/global'].keys()
    print(the_keys[0])
    no_bins=55
    the_density_in = np.zeros(no_bins-1)
    the_density_out = np.zeros(no_bins-1)
    for key in my_res_file['/global'].keys():
        print('working with %s' % key)
        count_hists=count_hists+1
        count_in = my_res_file['/global/%s/flows_within_tumor/flows_inside' %key]
        count_out = my_res_file['/global/%s/flows_within_tumor/flows_outside' %key]
        count_in = np.asarray(count_in)
        count_out= np.asarray(count_out)        
        density_in, bin_in_array = np.histogram(count_in,np.logspace(1.,6,no_bins),density=True)
        density_out, bin_out_array = np.histogram(count_out,np.logspace(1.,6,no_bins),density=True)
        the_density_in = the_density_in+density_in
        the_density_out = the_density_out+density_out
    print('we have %i histograms' % count_hists)
    the_density_in = np.divide(the_density_in,count_hists)
    the_density_out = np.divide(the_density_out,count_hists)
    mu1 = np.sum(the_density_in*np.diff(bin_in_array)*(0.5*(bin_in_array[1:]+bin_in_array[:-1])))
    quad1 = np.sum(the_density_in*np.diff(bin_in_array)*np.square((0.5*(bin_in_array[1:]+bin_in_array[:-1]))))
    std1 = np.sqrt(quad1-mu1**2)    
    mu2 = np.sum(the_density_out*np.diff(bin_out_array)*(0.5*(bin_out_array[1:]+bin_out_array[:-1])))
    quad2 = np.sum(the_density_out*np.diff(bin_out_array)*np.square((0.5*(bin_out_array[1:]+bin_out_array[:-1]))))
    std2 = np.sqrt(quad2-mu2**2)    
    plt.title('Flow distribution')    
    ax = plt.subplot(211)
    ax.bar(bin_out_array[:-1],the_density_out,np.diff(bin_out_array))
    ax.set_xscale('log')
    ax.text(10**4,0.0001,'$< q_{out} > = %f \pm %f$' % (mu2, std2))
    plt.legend(['Outside'])
    plt.ylabel('Propability')    
    ax = plt.subplot(212)    
    ax.bar(bin_in_array[:-1],the_density_in,np.diff(bin_in_array))
    ax.set_xscale('log')
    ax.text(10**4,0.0004,'$< q_{in} > = %f \pm %f$' % (mu1, std1))
    plt.xlabel('Volumeflow through vessel segments/ $\mu m$')
    plt.ylabel('Propability')
    
    plt.legend(['Inside'])
    plt.savefig('Flow_distribution.png')
    
    plt.show()

def create_mean_plot_velocities():
    my_res_file = h5files.open('analyzeVesselTumor.h5','r')
    count_hists = 0
    the_keys= my_res_file['/global'].keys()
    print(the_keys[0])
    no_bins=55
    the_density_in = np.zeros(no_bins-1)
    the_density_out = np.zeros(no_bins-1)
    for key in my_res_file['/global'].keys():
        print('working with %s' % key)
        count_hists=count_hists+1
        count_in = my_res_file['/global/%s/flows_within_tumor/flows_inside' %key]
        count_out = my_res_file['/global/%s/flows_within_tumor/flows_outside' %key]
        radii_in = np.asarray(my_res_file['/global/%s/radii_within_tumor/radii_inside' %key])
        radii_out = np.asarray(my_res_file['/global/%s/radii_within_tumor/radii_outside' %key])        
        count_in = np.asarray(count_in)/(np.pi*np.square(radii_in))
        count_out= np.asarray(count_out)/(np.pi*np.square(radii_out))
        density_in, bin_in_array = np.histogram(count_in,np.logspace(1,4,no_bins),density=True)
        density_out, bin_out_array = np.histogram(count_out,np.logspace(1,4,no_bins),density=True)
        the_density_in = the_density_in+density_in
        the_density_out = the_density_out+density_out
    print('we have %i histograms' % count_hists)
    the_density_in = np.divide(the_density_in,count_hists)
    the_density_out = np.divide(the_density_out,count_hists)
    mu1 = np.sum(the_density_in*np.diff(bin_in_array)*(0.5*(bin_in_array[1:]+bin_in_array[:-1])))
    quad1 = np.sum(the_density_in*np.diff(bin_in_array)*np.square((0.5*(bin_in_array[1:]+bin_in_array[:-1]))))
    std1 = np.sqrt(quad1-mu1**2)    
    mu2 = np.sum(the_density_out*np.diff(bin_out_array)*(0.5*(bin_out_array[1:]+bin_out_array[:-1])))
    quad2 = np.sum(the_density_out*np.diff(bin_out_array)*np.square((0.5*(bin_out_array[1:]+bin_out_array[:-1]))))
    std2 = np.sqrt(quad2-mu2**2)    
    plt.title('Flow distribution')    
    ax = plt.subplot(211)
    ax.bar(bin_out_array[:-1],the_density_out,np.diff(bin_out_array))
    ax.set_xscale('log')
    ax.text(10**3,0.002,'$< v_{out} > = %f \pm %f$' % (mu2, std2))
    plt.legend(['Outside'])
    plt.ylabel('Propability')    
    ax = plt.subplot(212)    
    ax.bar(bin_in_array[:-1],the_density_in,np.diff(bin_in_array))
    ax.set_xscale('log')
    ax.text(10**2,0.006,'$< v_{in} > = %f \pm %f$' % (mu1, std1))
    plt.xlabel('Volumeflow through vessel segments/ $\mu m$')
    plt.ylabel('Propability')
    
    plt.legend(['Inside'])
    plt.savefig('Flow_distribution.png')
    
    plt.show()

def create_mean_plot_pressures():
    my_res_file = h5files.open('analyzeVesselTumor.h5','r')
    count_hists = 0
    the_keys= my_res_file['/global'].keys()
    print(the_keys[0])
    no_bins=55
    the_density_in = np.zeros(no_bins)
    the_density_out = np.zeros(no_bins)
    for key in my_res_file['/global'].keys():
        print('working with %s' % key)
        count_hists=count_hists+1
        count_in = my_res_file['/global/%s/pressures_within_tumor/pressures_inside' %key]
        count_out = my_res_file['/global/%s/pressures_within_tumor/pressures_outside' %key]
        density_in, bin_in_array = np.histogram(count_in,no_bins, density=True)
        density_out, bin_out_array = np.histogram(count_out,no_bins, density=True)
        the_density_in = the_density_in+density_in
        the_density_out = the_density_out+density_out
    print('we have %i histograms' % count_hists)
    the_density_in = np.divide(the_density_in,count_hists)
    the_density_out = np.divide(the_density_out,count_hists)
    mu1 = np.sum(the_density_in*np.diff(bin_in_array)*(0.5*(bin_in_array[1:]+bin_in_array[:-1])))
    quad1 = np.sum(the_density_in*np.diff(bin_in_array)*np.square((0.5*(bin_in_array[1:]+bin_in_array[:-1]))))
    std1 = np.sqrt(quad1-mu1**2)    
    mu2 = np.sum(the_density_out*np.diff(bin_out_array)*(0.5*(bin_out_array[1:]+bin_out_array[:-1])))
    quad2 = np.sum(the_density_out*np.diff(bin_out_array)*np.square((0.5*(bin_out_array[1:]+bin_out_array[:-1]))))
    std2 = np.sqrt(quad2-mu2**2)    
    plt.title('Pressure distribution')    
    ax = plt.subplot(211)
    ax.bar(bin_out_array[:-1],the_density_out,np.diff(bin_out_array))
    ax.text(6,0.2,'$< p_{out} > = %f \pm %f$' % (mu2, std2))
    plt.legend(['Outside'])
    plt.ylabel('Propability')    
    ax = plt.subplot(212)    
    plt.bar(bin_in_array[:-1],the_density_in,np.diff(bin_in_array))
    ax.text(5,0.3,'$< p_{in} > = %f \pm %f$' % (mu1, std1))
    plt.xlabel('pressure within vessel segments/ $kPa$')
    plt.ylabel('Propability')
    
    plt.legend(['Inside'])
    plt.savefig('Pressure_distribution.png')
    
    plt.show()

def create_mean_plot_pressures_av():
    my_res_file = h5files.open('analyzeVesselTumor.h5','r')
    count_hists = 0
    the_keys= my_res_file['/global'].keys()
    print(the_keys[0])
    no_bins=55
    the_density_in_a = np.zeros(no_bins)
    the_density_in_v = np.zeros(no_bins)
    the_density_out_a = np.zeros(no_bins)
    the_density_out_v = np.zeros(no_bins)
    for key in my_res_file['/global'].keys():
        print('working with %s' % key)
        count_hists=count_hists+1
        count_in_a = my_res_file['/global/%s/pressures_within_tumor/pressures_inside_a' %key]
        count_in_v = my_res_file['/global/%s/pressures_within_tumor/pressures_inside_v' %key]        
        count_out_a = my_res_file['/global/%s/pressures_within_tumor/pressures_outside_a' %key]
        count_out_v = my_res_file['/global/%s/pressures_within_tumor/pressures_outside_v' %key]
        density_in_a, bin_in_a_array = np.histogram(count_in_a,no_bins, density=True)
        density_in_v, bin_in_v_array = np.histogram(count_in_v,no_bins, density=True)
        density_out_a, bin_out_a_array = np.histogram(count_out_a,no_bins, density=True)
        density_out_v, bin_out_v_array = np.histogram(count_out_v,no_bins, density=True)
        the_density_in_a = the_density_in_a+density_in_a
        the_density_in_v = the_density_in_v+density_in_v
        the_density_out_a = the_density_out_a+density_out_a
        the_density_out_v = the_density_out_v+density_out_v
    print('we have %i histograms' % count_hists)
    the_density_in_a = np.divide(the_density_in_a,count_hists)
    the_density_in_v = np.divide(the_density_in_v,count_hists)
    the_density_out_a = np.divide(the_density_out_a,count_hists)
    the_density_out_v = np.divide(the_density_out_v,count_hists)
    mu_in_a = np.sum(the_density_in_a*np.diff(bin_in_a_array)*(0.5*(bin_in_a_array[1:]+bin_in_a_array[:-1])))
    quad_in_a = np.sum(the_density_in_a*np.diff(bin_in_a_array)*np.square((0.5*(bin_in_a_array[1:]+bin_in_a_array[:-1]))))
    std_in_a = np.sqrt(quad_in_a-mu_in_a**2)
    mu_in_v = np.sum(the_density_in_v*np.diff(bin_in_v_array)*(0.5*(bin_in_v_array[1:]+bin_in_v_array[:-1])))
    quad_in_v = np.sum(the_density_in_v*np.diff(bin_in_v_array)*np.square((0.5*(bin_in_v_array[1:]+bin_in_v_array[:-1]))))
    std_in_v = np.sqrt(quad_in_v-mu_in_v**2)
    
    mu_out_a = np.sum(the_density_out_a*np.diff(bin_out_a_array)*(0.5*(bin_out_a_array[1:]+bin_out_a_array[:-1])))
    quad_out_a = np.sum(the_density_out_a*np.diff(bin_out_a_array)*np.square((0.5*(bin_out_a_array[1:]+bin_out_a_array[:-1]))))
    std_out_a = np.sqrt(quad_out_a-mu_out_a**2)
    mu_out_v = np.sum(the_density_out_v*np.diff(bin_out_v_array)*(0.5*(bin_out_v_array[1:]+bin_out_v_array[:-1])))
    quad_out_v = np.sum(the_density_out_v*np.diff(bin_out_v_array)*np.square((0.5*(bin_out_v_array[1:]+bin_out_v_array[:-1]))))
    std_out_v = np.sqrt(quad_out_v-mu_out_v**2)    
    
    fig1 = plt.figure()    
    ax = plt.subplot(211)
    ax.xaxis.grid(True)
    ax.bar(bin_out_a_array[:-1],the_density_out_a,np.diff(bin_out_a_array),alpha=0.3,color='r')
    ax.text(5.2,0.3,'$< p_a{out} > = %f \pm %f$' % (mu_out_a, std_out_a))
    ax.bar(bin_out_v_array[:-1],the_density_out_v,np.diff(bin_out_v_array),alpha=0.3,color='b')
    ax.text(5.2,0.2,'$< p_v{out} > = %f \pm %f$' % (mu_out_v, std_out_v))
    ax.legend(['Outside tumor artery','Outside tumor vein'])
    ax.set_ylabel('Propability')
    
    plt.gca().xaxis.set_major_locator(plt.NullLocator())
    #ax.setp(ax.get_xticklabels(),visible=False)  
    ax = plt.subplot(212)
    ax.xaxis.grid(True)    
    ax.bar(bin_in_a_array[:-1],the_density_in_a,np.diff(bin_in_a_array),alpha=0.3,color='r')
    ax.bar(bin_in_v_array[:-1],the_density_in_v,np.diff(bin_in_v_array),alpha=0.3,color='b')
    ax.text(5,0.3,'$< p_a{in} > = %f \pm %f$' % (mu_in_a, std_in_a))
    ax.text(5,0.2,'$< p_v{in} > = %f \pm %f$' % (mu_in_v, std_in_v))
    plt.xlabel('pressure within vessel segments/ $kPa$')
    plt.ylabel('Propability')
    fig1.suptitle('Pressure distribution av')
    plt.legend(['Inside tumor artery', 'Inside tumor vein'])
    plt.savefig('Pressure_distribution_av.png')
    
#    ax = plt.subplot(211)
#    plt.title('Pressure distribution veins')    
#    ax.bar(bin_out_v_array[:-1],the_density_out_v,np.diff(bin_out_v_array),alpha=0.3,color='b')
#    ax.text(6,0.2,'$< p_{out} > = %f \pm %f$' % (mu_out_v, std_out_v))
#    plt.legend(['Outside'])
#    plt.ylabel('Propability')    
#    ax = plt.subplot(212)    
#    plt.bar(bin_in_v_array[:-1],the_density_in_v,np.diff(bin_in_v_array))
#    ax.text(5,0.3,'$< p_{in} > = %f \pm %f$' % (mu_in_v, std_in_v))
#    plt.xlabel('pressure within vessel segments/ $kPa$')
#    plt.ylabel('Propability')
#    
#    plt.legend(['Inside'])
#    plt.savefig('Pressure_distribution_vein.png')
    plt.show()

def create_nice_file_containing_all_the_data():
  filenames = sys.argv[1:]
  files = [ h5files.open(fn, 'r') for fn in filenames ]
  fmeasure = h5files.open('analyzeVesselTumor.h5', 'a')
  dataman = myutils.DataManager(10, [DataVesselLength(),DataVesselTumor()])

  allgroups = defaultdict(list)
  for f in files:
    keys = filter(lambda k: k.startswith('out'), f.keys())
    for k in keys:
      allgroups[k].append(f[k])
  allgroups = [ (k, allgroups[k]) for k in sorted(allgroups.keys()) ]
  groupslist = [allgroups[-1]]
  outfn = 'Length-%s.pdf' % splitext(basename(filenames[0]))[0]
  pprint(groupslist)
  print '-> %s' % (outfn)
  
  for key,groups in groupslist:
    for group in groups:
      cachelocation = (fmeasure, '%s_file_%s_FILEID%s' % (group.name, group.file.filename, myutils.checksum(group.file.filename)))
      data = dataman.obtain_data('length_dist','length_dist', group['vessels'],group['tumor'], cachelocation)
      data = dataman.obtain_data('within_fake_tumor','radii_within_tumor', group['vessels'],group['tumor'], cachelocation)
      data = dataman.obtain_data('within_fake_tumor','pressures_within_tumor', group['vessels'],group['tumor'], cachelocation)
      data = dataman.obtain_data('within_fake_tumor','flows_within_tumor', group['vessels'],group['tumor'], cachelocation)
  print('you just created a big nice file')
  
if __name__== '__main__':
#first: we build a nice cache file where all
#the relevant data is stored
#    create_nice_file_containing_all_the_data()

#in case this file has been created
#we can now use it to plot the intended data
#    create_mean_plot_length()
#    create_mean_plot_radii()
#    create_mean_plot_pressures()
    create_mean_plot_pressures_av()
#    create_mean_plot_flows()
#    create_mean_plot_velocities()
    