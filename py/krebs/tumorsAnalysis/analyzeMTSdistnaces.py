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

if __name__ == '__main__':
  import os.path, sys
  sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'..'))
  sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../..'))

import povrayRenderCells
import myutils

from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

no_bins = 50

def example_plot():
    x_edges = np.arange(6)
    y_edges = np.arange(6)
    data = np.random.rand(340,2)*5
    
    ### using numpy.histogram2d
    bin_values,_,__ = np.histogram2d(data[:,0],data[:,1],bins=(x_edges, y_edges) )
    X, Y = np.meshgrid(x_edges,y_edges)
    
    fig, (ax,ax2) = plt.subplots(ncols=2)
    ax.set_title("numpy.histogram2d \n + plt.pcolormesh")
    ax.pcolormesh(X, Y, bin_values.T)
    
    ### using plt.hist2d
    ax2.set_title("plt.hist2d")
    ax2.hist2d(data[:,0],data[:,1],bins=(x_edges, y_edges))
    
    
    plt.show()
    
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
    ax1.set_xlabel(r'distance from nearest vessel / $\mu m')
    ax1.set_ylabel(r'probability')
    if infos:
      ax1.set_title('file: %s \n at %s' % (goodArguments.vbl_simulation_output_filename, goodArguments.output_grp_name))
    ax1.grid()
    if interactive:
      plt.show()
    else:
      pp.savefig()
      
#def scatter_cell_endity_vs_distances_to_next_vessel(endity, out_grp_name,pp):
#  with h5py.File(goodArguments.vbl_simulation_output_filename, 'r') as h5_f:
#    h5_out_grp = h5_f[goodArguments.output_grp_name]
#    distances_to_nearest_vessel = np.asarray(h5_out_grp['cells/distance_to_nearest_vessel'])
#    endity_value_of_cells = np.asarray(h5_out_grp['cells/' + endity])
#    
#    fig1 = plt.figure()
#    ax1 = fig1.add_subplot(111)
#    ax1.scatter(distances_to_nearest_vessel, endity_value_of_cells)
#    plt.xlabel(r'distance from nearest vessel')
#    plt.ylabel(r' %s of cell' % endity)
#    plt.title('file: %s \n at %s' % (goodArguments.vbl_simulation_output_filename, goodArguments.output_grp_name))
#    plt.grid()
#    if interactive:
#      plt.show()
#    else:
#      pp.savefig()


class Hist2d_data(object):
  keywords = ['hist_2d_data']
  
  def obtain_data(self, dataman, dataname, *args):
    if dataname == 'hist_2d_data':
      print(args)
      endity = args[0]
      this_out_grp_name =args[1]
      no_of_bins = args[2]
      def read(hdf_cache_grp, data_name):
        print('read data at: %s' %hdf_cache_grp.name)
        print('data_name: %s ' % data_name)
        h = np.asarray(hdf_cache_grp[data_name + '/' + 'h'])
        xedges = np.asarray(hdf_cache_grp[data_name + '/' + 'xedges'])
        yedges = np.asarray(hdf_cache_grp[data_name + '/' + 'yedges'])
        return (h, xedges, yedges)
      def write(hdf_cache_grp, data_name):
        this_out_grp = hdf_cache_grp.create_group(data_name)
        (h, xedges, yedges) = create_2d_histo(endity)
        #(average_value, errors, distances) = sample_line_general('o2', this_out_grp_name, vein_parallel_p1, vein_parallel_p2)
        #group_of_single_timepoint = hdf_cache_grp.create_group(this_out_grp_name)
        this_out_grp.create_dataset('h', data=h)
        this_out_grp.create_dataset('xedges', data=xedges)
        this_out_grp.create_dataset('yedges', data=yedges)
        print('created data at: %s' % this_out_grp.name)
#        print(hdf_cache_grp)
#        print(data_name)
#        print('before create')
#        group_of_single_timepoint.create_dataset('average_po2', data=average_value)
#        group_of_single_timepoint.create_dataset('average_po2_error', data=errors)
#        group_of_single_timepoint.create_dataset('distances', data=distances)
      possible_hdf_group_name = 'hist_2d_data_bins_%s/' % no_of_bins
      possible_hdf_group_name = possible_hdf_group_name+'/' + endity
      if not possible_hdf_group_name in f_cache:
        f_cache.create_group(possible_hdf_group_name)
    return myutils.hdf_data_caching(read, write, f_cache, possible_hdf_group_name)
def create_2d_histo(endity):
  with h5py.File(goodArguments.vbl_simulation_output_filename, 'r') as h5_f:
    h5_out_grp = h5_f[goodArguments.output_grp_name]
    distances_to_nearest_vessel = np.asarray(h5_out_grp['cells/distance_to_nearest_vessel'])
    distances_to_nearest_vessel = distances_to_nearest_vessel[:,0]
    endity_value_of_cells = np.asarray(h5_out_grp['cells/' + endity])
    endity_value_of_cells = endity_value_of_cells[:,0]
    cell_radii = np.asarray(h5_out_grp['cells/cell_radii'])
    cell_radii=cell_radii[:,0]
    
    if(endity == 'o2'):
      endity_value_of_cells = povrayRenderCells.convert_to_mmHg(endity_value_of_cells, cell_radii)
    #h, xedges, yedges, image) = plt.hist2d(distances_to_nearest_vessel, endity_value_of_cells, bins = no_bins,norm=matplotlib.colors.LogNorm())
    hist_return = np.histogram2d(distances_to_nearest_vessel, endity_value_of_cells, bins = no_bins,normed=True)
    h = hist_return[0]
    xedges = hist_return[1]
    yedges = hist_return[2]
    return (h, xedges, yedges)

def hist_cell_endity_vs_distances_to_next_vessel(endity, out_grp_name, no_of_bins,pp):
  h, xedges, yedges = dataman.obtain_data('hist_2d_data', endity, out_grp_name, no_of_bins)
    
  fig1 = plt.figure()
  ax1 = fig1.add_subplot(111)
  
  X,Y = np.meshgrid(xedges,yedges)
  im = ax1.pcolormesh(X,Y, h.T,norm=matplotlib.colors.LogNorm())
  fig1.colorbar(im, ax=ax1)
  ax1.set_xlabel(r'distance from nearest vessel')
  if endity == 'o2':
      ax1.set_ylabel(r'$pO_2 / mmHg$')
  else:
    ax1.set_ylabel(r' %s of cell' % endity)
  if infos:
    fig1.set(title = 'file: %s \n at %s' % (goodArguments.vbl_simulation_output_filename, goodArguments.output_grp_name))
  
  if interactive:
    plt.show()
  else:
    pp.savefig()

class Boxplot_data(object):
  keywords = ['boxplot_data']
  
  def obtain_data(self, dataman, dataname, *args):
    if dataname == 'boxplot_data':
      print(args)
      endity = args[0]
      this_out_grp_name =args[1]
      no_of_bins = args[2]
      def read(hdf_cache_grp, data_name):
        print('read data at: %s' %hdf_cache_grp.name)
        print('data_name: %s ' % data_name)
        big_data=list()
        x_labels = list()
        for the_key in hdf_cache_grp[data_name].keys():
          big_data.append(np.asarray(hdf_cache_grp[data_name][the_key]))
          x_labels.append(str(the_key))
        return (big_data,x_labels)
      def write(hdf_cache_grp, data_name):
        this_out_grp = hdf_cache_grp.create_group(data_name)
        (big_data, my_x_labels) = create_box_plot(this_out_grp_name, endity, no_of_bins)
        #(average_value, errors, distances) = sample_line_general('o2', this_out_grp_name, vein_parallel_p1, vein_parallel_p2)
        #group_of_single_timepoint = hdf_cache_grp.create_group(this_out_grp_name)
        for (data_around, center_pos) in zip(big_data, my_x_labels):
          proper_string = '%03i' % (int)(center_pos)
          this_out_grp.create_dataset(proper_string, data=data_around)
        print('created data at: %s' % this_out_grp.name)

      possible_hdf_group_name = 'box_plot_data_bins_%s/' % no_of_bins
      possible_hdf_group_name = possible_hdf_group_name+'/' + endity
      if not possible_hdf_group_name in f_cache:
        f_cache.create_group(possible_hdf_group_name)
          
      #return myutils.hdf_data_caching(read, write, f_cache[possible_hdf_group_name], this_out_grp_name)
      return myutils.hdf_data_caching(read, write, f_cache, possible_hdf_group_name)
    
def create_box_plot(this_out_grp_name, endity, no_of_bins):
  with h5py.File(goodArguments.vbl_simulation_output_filename, 'r') as h5_f:
    h5_out_grp = h5_f[this_out_grp_name]
    distances_to_nearest_vessel = np.asarray(h5_out_grp['cells/distance_to_nearest_vessel'])
    distances_to_nearest_vessel = distances_to_nearest_vessel[:,0]
    min_distance_to_nearest_vessel = np.min(distances_to_nearest_vessel)
    max_distance_to_nearest_vessel = np.max(distances_to_nearest_vessel)
    print('min: %f, max: %f' %(min_distance_to_nearest_vessel,max_distance_to_nearest_vessel))
    
    endity_value_of_cells = np.asarray(h5_out_grp['cells/' + endity])
    endity_value_of_cells = endity_value_of_cells[:,0]
    cell_radii = np.asarray(h5_out_grp['cells/cell_radii'])
    cell_radii=cell_radii[:,0]
    
  if(endity == 'o2'):
    endity_value_of_cells = povrayRenderCells.convert_to_mmHg(endity_value_of_cells, cell_radii)
  
  max_endity_value_of_cells = np.max(endity_value_of_cells)
  min_endity_value_of_cells = np.min(endity_value_of_cells)
  
#    max_endity_value_of_cells = 200
#    min_endity_value_of_cells = 20
  ''' this is not showing something -> I try it the other way round'''
  bins =[]    
  
  big_data = list()
  my_x_labels = list()
  
  ''' this will plot with distance '''
  
  diff = (max_distance_to_nearest_vessel - min_distance_to_nearest_vessel)/no_of_bins
  print('diff: %f' % diff)
  for i in range(no_of_bins):
    bins.append(min_distance_to_nearest_vessel+i*diff)
  print('bins:')
  print(bins)
  
  for (i, a_lower_bound) in enumerate(bins):
    upper_bound = a_lower_bound+diff
    print('lower: %f, upper: %f' %(a_lower_bound, upper_bound))
    good_indexes = np.where(np.logical_and(distances_to_nearest_vessel<upper_bound, distances_to_nearest_vessel > a_lower_bound))
#    good_indexes = np.where(endity_value_of_cells>=min_distance_to_nearest_vessel)    
    print('found %i indeces for %f' % (len(good_indexes[0]), a_lower_bound))
    data_on_this = endity_value_of_cells[good_indexes]
    print('min: %f, max: %f' % (np.min(data_on_this),np.max(data_on_this)))
    big_data.append(endity_value_of_cells[good_indexes])
    my_x_labels.append('%1.0f' % float(a_lower_bound+0.5*diff))
  return (big_data, my_x_labels)


def plot_cell_endity_vs_distances_to_next_vessel(out_grp_name,endity, no_of_bins,pp):
  big_data, x_ticks_labels = dataman.obtain_data('boxplot_data', out_grp_name, endity, no_of_bins)
  fig1 = plt.figure()
  ax1 = fig1.add_subplot(111)
  #ax1.scatter(distances_to_nearest_vessel[0:1000], endity_value_of_cells[0:1000])
  ax1.boxplot(big_data)
  ax1.set_xticklabels(x_ticks_labels,rotation=75)
  ax1.set_xticks(np.arange(len(x_ticks_labels))+1)
  if endity == 'o2':
    ax1.set_ylabel(r'pO2 / mmHg')
    ax1.set_xlabel(r' distance to nearest vessle/ $\mu m$')
    
  if endity == 'pH_ex':
    ax1.set_ylabel(r'pH')
    ax1.set_xlabel(r' distance to nearest vessle/ $\mu m$')
  ax1.set(title='file: %s \n at %s' % (goodArguments.vbl_simulation_output_filename, goodArguments.output_grp_name))
  ax1.grid(color='k', linestyle=':', linewidth=0.5)
  pp.savefig()    
  
  
  ''' this will plot the data on the absissis '''
#  bins =[]
#  big_data = list()
#  my_x_labels = list()
#  diff = (max_endity_value_of_cells - min_endity_value_of_cells)/no_of_bins
#  print('diff: %f' % diff)
#  for i in range(no_of_bins):
#    bins.append(min_endity_value_of_cells+i*diff)
#  print('bins:')
#  print(bins)
#  for (i, a_lower_bound) in enumerate(bins):
#    upper_bound = a_lower_bound+diff
#    print('lower: %f, upper: %f' %(a_lower_bound, upper_bound))
#    good_indexes = np.where(np.logical_and(endity_value_of_cells<upper_bound, endity_value_of_cells > a_lower_bound))
##    good_indexes = np.where(endity_value_of_cells>=min_distance_to_nearest_vessel)    
#    print('found %i indeces for %f' % (len(good_indexes[0]), a_lower_bound))
#    data_on_this = distances_to_nearest_vessel[good_indexes]
#    print('min: %f, max: %f' % (np.min(data_on_this),np.max(data_on_this)))
#    big_data.append(data_on_this)
#    #endity_value_of_cells[good_indexes
#    if(endity == 'pH_ex'):
#      my_x_labels.append('%.2f' % float(a_lower_bound+0.5*diff))
#    else:
#      my_x_labels.append('%1.0f' % float(a_lower_bound+0.5*diff))
#  fig2 = plt.figure()
#  ax2 = fig2.add_subplot(111)
#  #ax1.scatter(distances_to_nearest_vessel[0:1000], endity_value_of_cells[0:1000])
#  ax2.boxplot(big_data)
#  ax2.set_xticklabels(my_x_labels,rotation=75)
#  ax2.set_xticks(np.arange(len(my_x_labels))+1)
#  if endity == 'o2':
#    ax2.set_xlabel(r'pO2 / mmHg')
#    ax2.set_ylabel(r' distance to nearest vessle/ $\mu m$')
#    
#  if endity == 'pH_ex':
#    ax2.set_xlabel(r'pH')
#    ax2.set_ylabel(r' distance to nearest vessle/ $\mu m$')
#  ax2.set(title='file: %s \n at %s' % (goodArguments.vbl_simulation_output_filename, goodArguments.output_grp_name))
#  ax2.grid(color='k', linestyle=':', linewidth=0.5)
#  pp.savefig()  
      
  '''WHAT ABOUT THE CELLS ON THE SURFACE AND NOT SURFACE?'''
    
if __name__ == '__main__':
  import argparse
  parser = argparse.ArgumentParser(description='Analyze MTS distances')  
  parser.add_argument('--s',dest='vbl_simulation_output_filename', type=str, default='safe.h5', help='output file name in hdf5 format')
  parser.add_argument('--g',dest='output_grp_name', type=str, default='out0001', help='output group withing hdf5 file')
  interactive = False;
  infos = False;
  goodArguments, otherArguments = parser.parse_known_args()
  
  
  ''' begin of code '''
  '''register a clases at data manager'''
  with h5py.File('cache_'+ goodArguments.vbl_simulation_output_filename, 'a') as f_cache:
    dataman = myutils.DataManager(20, [Hist2d_data(), Boxplot_data()])
    
    
  #  pp.attach_note(r"$\beta$ ", positionRect=[-100,-100,0,0])
  #  pp.attach_note(r'$\beta$ ')
  #  pp.attach_note("klsdfjal")
    
    with PdfPages('analysisMTS_prop_%s_%s.pdf' % (goodArguments.vbl_simulation_output_filename, goodArguments.output_grp_name)) as pp:
      distances_to_vessels(goodArguments, pp);
    cell_endities = ['o2','pH_ex','cell_radii']
    #cell_endities = ['o2'] 
    for cell_endity in cell_endities:
      with PdfPages('analysisMTS_hist_%s_%s_%s.pdf' % (goodArguments.vbl_simulation_output_filename[:-3], cell_endity, goodArguments.output_grp_name)) as pp:
        hist_cell_endity_vs_distances_to_next_vessel(cell_endity, goodArguments.output_grp_name, no_bins,pp)
        plot_cell_endity_vs_distances_to_next_vessel(cell_endity, goodArguments.output_grp_name, no_bins,pp)