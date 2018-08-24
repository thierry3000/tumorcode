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

import h5py
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.mplot3d import Axes3D

no_bins = 30

def distances_to_vessels(out_grp_name,pp):
  with h5py.File(goodArguments.vbl_simulation_output_filename, 'r') as h5_f:
    h5_out_grp = h5_f[goodArguments.output_grp_name]
    distances_to_nearest_vessel = np.asarray(h5_out_grp['cells/distance_to_nearest_vessel'])
    counts, bin_edges = np.histogram(distances_to_nearest_vessel, density=True)
    width = bin_edges[1]-bin_edges[0]
    centers = (bin_edges[:-1] +bin_edges[1:])/2
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax1.plot(centers, counts*width)
    ax1.set_xlabel(r'distance from nearest vessel')
    ax1.set_ylabel(r'probability')
    ax1.set_title('file: %s \n at %s' % (goodArguments.vbl_simulation_output_filename, goodArguments.output_grp_name))
    ax1.grid()
    if interactive:
      plt.show()
    else:
      pp.savefig()
      
def scatter_cell_endity_vs_distances_to_next_vessel(endity, out_grp_name,pp):
  with h5py.File(goodArguments.vbl_simulation_output_filename, 'r') as h5_f:
    h5_out_grp = h5_f[goodArguments.output_grp_name]
    distances_to_nearest_vessel = np.asarray(h5_out_grp['cells/distance_to_nearest_vessel'])
    endity_value_of_cells = np.asarray(h5_out_grp['cells/' + endity])
    
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax1.scatter(distances_to_nearest_vessel, endity_value_of_cells)
    plt.xlabel(r'distance from nearest vessel')
    plt.ylabel(r' %s of cell' % endity)
    plt.title('file: %s \n at %s' % (goodArguments.vbl_simulation_output_filename, goodArguments.output_grp_name))
    plt.grid()
    if interactive:
      plt.show()
    else:
      pp.savefig()

def hist_cell_endity_vs_distances_to_next_vessel(endity, out_grp_name,pp):
  with h5py.File(goodArguments.vbl_simulation_output_filename, 'r') as h5_f:
    h5_out_grp = h5_f[goodArguments.output_grp_name]
    distances_to_nearest_vessel = np.asarray(h5_out_grp['cells/distance_to_nearest_vessel'])
    distances_to_nearest_vessel = distances_to_nearest_vessel[:,0]
    endity_value_of_cells = np.asarray(h5_out_grp['cells/' + endity])
    endity_value_of_cells = endity_value_of_cells[:,0]
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    plt.hist2d(distances_to_nearest_vessel, endity_value_of_cells, bins = no_bins,norm=matplotlib.colors.LogNorm())
    
    plt.xlabel(r'distance from nearest vessel')
    plt.ylabel(r' %s of cell' % endity)
    plt.title('file: %s \n at %s' % (goodArguments.vbl_simulation_output_filename, goodArguments.output_grp_name))
    plt.grid()
    plt.colorbar()
    if interactive:
      plt.show()
    else:
      pp.savefig()

def plot_cell_endity_vs_distances_to_next_vessel(endity, out_grp_name,pp):
  with h5py.File(goodArguments.vbl_simulation_output_filename, 'r') as h5_f:
    h5_out_grp = h5_f[goodArguments.output_grp_name]
    distances_to_nearest_vessel = np.asarray(h5_out_grp['cells/distance_to_nearest_vessel'])
    distances_to_nearest_vessel = distances_to_nearest_vessel[:,0]
    min_distance_to_nearest_vessel = np.min(distances_to_nearest_vessel)
    max_distance_to_nearest_vessel = np.max(distances_to_nearest_vessel)
    print('min: %f, max: %f' %(min_distance_to_nearest_vessel,max_distance_to_nearest_vessel))
    
#    spread = np.random.rand(50) * 100
#    center = np.ones(25) * 50
#    flier_high = np.random.rand(10) * 100 + 100
#    flier_low = np.random.rand(10) * -100
#    data = np.concatenate((spread, center, flier_high, flier_low), 0)
#    
#    # basic plot
#    plt.boxplot(data)    
    
    endity_value_of_cells = np.asarray(h5_out_grp['cells/' + endity])
    endity_value_of_cells = endity_value_of_cells[:,0]
    volume_o2_ml = endity_value_of_cells*1e-9/1.429
    ''' o2 density 1.429 g/L --> 1.429*10^9 pg/ml
    '''
    solubility = 2.8e-4 #ml O2/cm^3 mmHg
    solubility = solubility*1e-12 #ml O2/mum^3 mmHg
    endity_value_of_cells = volume_o2_ml/solubility
#    max_endity_value_of_cells = np.max(endity_value_of_cells)
#    min_endity_value_of_cells = np.min(endity_value_of_cells)
    
    max_endity_value_of_cells = 200
    min_endity_value_of_cells = 10
    ''' this is not showing something -> I try it the other way round'''
    
    bins =[]    
    no_of_bins = 10
    goWithDistance = False
    if goWithDistance:
      diff = (max_distance_to_nearest_vessel - min_distance_to_nearest_vessel)/no_of_bins
      print('diff: %f' % diff)
      for i in range(no_of_bins):
        bins.append(min_distance_to_nearest_vessel+i*diff)
      print('bins:')
      print(bins)
      big_data = list()
      for (i, a_lower_bound) in enumerate(bins):
        upper_bound = a_lower_bound+diff
        print('lower: %f, upper: %f' %(a_lower_bound, upper_bound))
        good_indexes = np.where(np.logical_and(distances_to_nearest_vessel<upper_bound, distances_to_nearest_vessel > a_lower_bound))
  #    good_indexes = np.where(endity_value_of_cells>=min_distance_to_nearest_vessel)    
        print('found %i indeces for %f' % (len(good_indexes[0]), a_lower_bound))
        data_on_this = endity_value_of_cells[good_indexes]
        print('min: %f, max: %f' % (np.min(data_on_this),np.max(data_on_this)))
        big_data.append(endity_value_of_cells[good_indexes])
    
    
    else:
      diff = (max_endity_value_of_cells - min_endity_value_of_cells)/no_of_bins
      print('diff: %f' % diff)
      for i in range(no_of_bins):
        bins.append(min_endity_value_of_cells+i*diff)
      print('bins:')
      print(bins)
      big_data = list()
      for (i, a_lower_bound) in enumerate(bins):
        upper_bound = a_lower_bound+diff
        print('lower: %f, upper: %f' %(a_lower_bound, upper_bound))
        good_indexes = np.where(np.logical_and(endity_value_of_cells<upper_bound, endity_value_of_cells > a_lower_bound))
  #    good_indexes = np.where(endity_value_of_cells>=min_distance_to_nearest_vessel)    
        print('found %i indeces for %f' % (len(good_indexes[0]), a_lower_bound))
        data_on_this = distances_to_nearest_vessel[good_indexes]
        print('min: %f, max: %f' % (np.min(data_on_this),np.max(data_on_this)))
        big_data.append(endity_value_of_cells[good_indexes])
      #plt.boxplot(big_data)
      #plt.show()
      
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    #ax1.scatter(distances_to_nearest_vessel[0:1000], endity_value_of_cells[0:1000])
    ax1.boxplot(big_data)
    if goWithDistance:
      plt.ylabel(r'pO2 / mmHg')
      plt.xlabel(r' distance to nearest vessle/ $\mu m$')
    else:
      plt.xlabel(r'pO2 / mmHg')
      plt.ylabel(r' distance to nearest vessle/ $\mu m$')
    plt.title('file: %s \n at %s' % (goodArguments.vbl_simulation_output_filename, goodArguments.output_grp_name))
    plt.grid()
    if interactive:
      plt.show()
    else:
      pp.savefig()
      
    '''WHAT ABOUT THE CELLS ON THE SURFACE AND NOT SURFACE?'''
    
    
if __name__ == '__main__':
  import argparse
  parser = argparse.ArgumentParser(description='Analyze MTS distances')  
  parser.add_argument('--s',dest='vbl_simulation_output_filename', type=str, default='safe.h5', help='output file name in hdf5 format')
  parser.add_argument('--g',dest='output_grp_name', type=str, default='out0001', help='output group withing hdf5 file')
  interactive = False;
  goodArguments, otherArguments = parser.parse_known_args()
  
  with PdfPages('analysisMTS_%s_%s.pdf' % (goodArguments.vbl_simulation_output_filename, goodArguments.output_grp_name)) as pp:
    pp.attach_note(r"$\beta$ ", positionRect=[-100,-100,0,0])
    pp.attach_note(r'$\beta$ ')
    pp.attach_note("klsdfjal")
    
    distances_to_vessels(goodArguments, pp);
#    cell_endities = ['o2', 'pH_ex', 'glucose_ex','cell_age','cell_no_neigh', 'cell_o2_consumption_rate', 'cell_phase', 'cell_phase_age','cell_radii']
    cell_endities = ['o2']    
    for cell_endity in cell_endities:
      #scatter_cell_endity_vs_distances_to_next_vessel(cell_endity, goodArguments, pp)
      #hist_cell_endity_vs_distances_to_next_vessel(cell_endity, goodArguments, pp)
      plot_cell_endity_vs_distances_to_next_vessel(cell_endity, goodArguments, pp)