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

no_bins = 30
max_distance_to_vessel = 120
intervall_for_plot_numbers=20

my_markers = ['o', 'v', 's','h', 'x', 'd']
my_colors   = ['b', 'g', 'c', 'm', 'k', 'k']
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
    counts, bin_edges = np.histogram(distances_to_nearest_vessel, bins=no_bins,density=True)
    width = bin_edges[1]-bin_edges[0]
    centers = (bin_edges[:-1] +bin_edges[1:])/2
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax1.plot(centers, counts*width, '*')
    ax1.set_xlabel(r'distance from nearest vessel / $\mu m')
    ax1.set_ylabel(r'probability')
    if infos:
      ax1.set_title('file: %s \n at %s' % (goodArguments.vbl_simulation_output_filename, goodArguments.output_grp_name))
    ax1.grid()
    if interactive:
      plt.show()
    else:
      pp.savefig()

def distances_to_vessels_for_multiple_group(out_grp_name,pp):
  with h5py.File(goodArguments.vbl_simulation_output_filename, 'r') as h5_f:
    
    centers_list =[]
    for aOutGroup in out_groups_to_consider:
      h5_out_grp = h5_f[ aOutGroup]
      distances_to_nearest_vessel = np.asarray(h5_out_grp['cells/distance_to_nearest_vessel'])
      counts, bin_edges = np.histogram(distances_to_nearest_vessel, bins=no_bins,density=True)
      width = bin_edges[1]-bin_edges[0]
      centers = (bin_edges[:-1] +bin_edges[1:])/2
      centers_list.append(centers)
      
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    abc=np.asarray(centers_list)
    ax1.plot(abc.transpose(), counts*width, '*')
    
    ax1.legend(out_groups_to_hours(out_groups_to_consider))
    if infos:
      ax1.set_xlabel(r'distance from nearest vessel / $\mu m$')
      ax1.set_ylabel(r'probability')
      #ax1.set_title('file: %s \n at %s' % (goodArguments.vbl_simulation_output_filename, goodArguments.output_grp_name))
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
      possible_hdf_group_name = '%s/hist_2d_data_bins_%s/' % (this_out_grp_name, no_of_bins)
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
    if 'max_distance_to_vessel' in globals():
      ymin=np.min(endity_value_of_cells)
      ymax=np.max(endity_value_of_cells)
      hist_return = np.histogram2d(distances_to_nearest_vessel, endity_value_of_cells, bins = no_bins,normed=True, range=[[0, max_distance_to_vessel],[ymin,ymax]])
    else:
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
  keywords = ['boxplot_data', 'boxplot_data_multiple']
  
  def obtain_data(self, dataman, dataname, *args):
    if dataname == 'boxplot_data':
      print(args)
      this_out_grp_name =args[0]
      endity = args[1]
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
          print(proper_string)
          this_out_grp.create_dataset(proper_string, data=data_around)
        print('created data at: %s' % this_out_grp.name)

      possible_hdf_group_name = '%s/box_plot_data_bins_%s/' % (this_out_grp_name,no_of_bins)
      possible_hdf_group_name = possible_hdf_group_name+'/' + endity
      if not possible_hdf_group_name in f_cache:
        f_cache.create_group(possible_hdf_group_name)
          
      #return myutils.hdf_data_caching(read, write, f_cache[possible_hdf_group_name], this_out_grp_name)
      return myutils.hdf_data_caching(read, write, f_cache, possible_hdf_group_name)
    if dataname == 'boxplot_data_multiple':
      print(args)
      multiple_out_grp_names =args[0]
      endity = args[1]
      no_of_bins = args[2]
      def read(hdf_cache_grp, data_name):
        print('read data at: %s' %hdf_cache_grp.name)
        print('data_name: %s ' % data_name)
        big_big_data=dict()
        
        for aOutGroup in multiple_out_grp_names:
          big_big_data[aOutGroup] = list()
          x_labels = list()
          for the_key in hdf_cache_grp[data_name+'/'+aOutGroup].keys():
            big_big_data[aOutGroup].append(np.asarray(hdf_cache_grp[data_name+'/'+aOutGroup][the_key]))
            x_labels.append(str(the_key))
        return (big_big_data,x_labels)
      def write(hdf_cache_grp, data_name):
        for this_out_grp_name in multiple_out_grp_names:
          this_out_grp = hdf_cache_grp.create_group(data_name+'/'+this_out_grp_name)
          (big_data, my_x_labels) = create_box_plot(this_out_grp_name, endity, no_of_bins)
        #(average_value, errors, distances) = sample_line_general('o2', this_out_grp_name, vein_parallel_p1, vein_parallel_p2)
        #group_of_single_timepoint = hdf_cache_grp.create_group(this_out_grp_name)
          for (data_around, center_pos) in zip(big_data, my_x_labels):
            proper_string = '%03i' % (int)(center_pos)
            #print(proper_string)
            this_out_grp.create_dataset(proper_string, data=data_around)
        print('created data at: %s' % this_out_grp.name)

      possible_hdf_group_name = 'box_plot_data_bins_%s_multiple/' % (no_of_bins)
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
    
    if 'max_distance_to_vessel' in globals():
      min_distance_to_nearest_vessel= 0.
      max_distance_to_nearest_vessel= max_distance_to_vessel
      
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
    if(len(data_on_this)>0):
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
  if infos:
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
def plot_cell_endity_vs_distances_to_next_vessel_multiple_times(multiple_out_grp_names, endity, no_of_bins,pp):
  big_data, x_ticks_labels = dataman.obtain_data('boxplot_data_multiple', multiple_out_grp_names, endity, no_of_bins)
  fig1 = plt.figure()
  ax1 = fig1.add_subplot(111)
  #ax1.scatter(distances_to_nearest_vessel[0:1000], endity_value_of_cells[0:1000])
  ax1.boxplot(big_data[multiple_out_grp_names[-1]])
  ax1.set_color_cycle(my_colors)
  #ax1.set_marker_cycle(my_markers)
  for (i,outGrpName) in enumerate(multiple_out_grp_names[0:-1]):
    list_of_avg=list()
    patchList = list()
    for aDataEntry in big_data[outGrpName]:
      this_avg= np.average(aDataEntry)
      print(this_avg)
      list_of_avg.append(this_avg)
    ax1.plot(list_of_avg, linestyle='None', marker=my_markers[i], label='%s h' % outGrpName[4:])
  ax1.set_xticklabels(x_ticks_labels,rotation=75)
  ax1.set_xticks(np.arange(len(x_ticks_labels))+1)
  #legendlabel = ['%s hours' % entry[4:] for entry in multiple_out_grp_names]
  ax1.legend()
  if infos:
    ax1.set_xlabel(r' distance to nearest vessle/ $\mu m$')
    if endity == 'o2':
      ax1.set_ylabel(r'pO2 / $mmHg$')
      
    if endity == 'pH_ex':
      ax1.set_ylabel(r'pH')
    
    if endity == 'cell_radii':
      ax1.set_ylabel(r'radius of cells/ $\mu m$')
      
    
    #ax1.set(title='file: %s \n at %s' % (goodArguments.vbl_simulation_output_filename, goodArguments.output_grp_name))
  ax1.grid(color='k', linestyle=':', linewidth=0.5)
  pp.savefig()  


class Dev_from_sphere_data(object):
  keywords = ['dev_from_sphere']
  
  def obtain_data(self, dataman, dataname, *args):
    if dataname == 'dev_from_sphere':
      print(args)
      this_out_grp_name =args[0]
      no_of_bins = args[1]
      rangeMin = args[2]
      rangeMax = args[3]
      def read(hdf_cache_grp, data_name):
        print('read data at: %s' %hdf_cache_grp.name)
        print('data_name: %s ' % data_name)
        hist_data = np.asarray(hdf_cache_grp[data_name + '/hist_data'])
        hist_edges= np.asarray(hdf_cache_grp[data_name + '/hist_edges'])
        return (hist_data,hist_edges)
      def write(hdf_cache_grp, data_name):
        with h5py.File(goodArguments.vbl_simulation_output_filename, 'r') as h5_f:
          isonAS = np.asarray(h5_f[this_out_grp_name + '/vbl/isonAS'], dtype=bool)
          isonAS = isonAS[:,0]
          pos = np.asarray(h5_f[this_out_grp_name +'/cells/cell_center_pos'])
          dist_to_center = np.sqrt(np.sum(np.power(pos,2),1))
          hist, bin_edges = np.histogram(dist_to_center[isonAS], bins=no_of_bins, range=[rangeMin, rangeMax])
          #print(pos)

        this_out_grp = hdf_cache_grp.create_group(data_name)
        this_out_grp.create_dataset('hist_data', data=hist)
        this_out_grp.create_dataset('hist_edges', data=bin_edges)
        print('created data at: %s' % this_out_grp.name)

      possible_hdf_group_name = '%s/dev_from_sphere_%s/' % (this_out_grp_name,no_of_bins)
      #possible_hdf_group_name = possible_hdf_group_name+'/' + endity
      if not possible_hdf_group_name in f_cache:
        f_cache.create_group(possible_hdf_group_name)
          
      #return myutils.hdf_data_caching(read, write, f_cache[possible_hdf_group_name], this_out_grp_name)
      return myutils.hdf_data_caching(read, write, f_cache, possible_hdf_group_name)
    

      possible_hdf_group_name = 'box_plot_data_bins_%s_multiple/' % (no_of_bins)
      possible_hdf_group_name = possible_hdf_group_name+'/' + endity
      if not possible_hdf_group_name in f_cache:
        f_cache.create_group(possible_hdf_group_name)
          
      #return myutils.hdf_data_caching(read, write, f_cache[possible_hdf_group_name], this_out_grp_name)
      return myutils.hdf_data_caching(read, write, f_cache, possible_hdf_group_name)
def plot_dev_from_sphere_multiple_times(multiple_out_grp_names, no_of_bins,pp):
  rangeMin = 0.0
  rangeMax = 400.5
  fig1 = plt.figure()
  #multiple_out_grp_names = [multiple_out_grp_names[-3]]
  for aOutGroupName in multiple_out_grp_names:
  #aOutGroupName = multiple_out_grp_names[0] 
    hist_data, bin_edges = dataman.obtain_data('dev_from_sphere', aOutGroupName, 1000, rangeMin, rangeMax)
    bin_centers = bin_edges[0:-1]+0.5*(bin_edges[1:]-bin_edges[:-1])
  
    ax1 = fig1.add_subplot(111)
  #ax1.scatter(distances_to_nearest_vessel[0:1000], endity_value_of_cells[0:1000])
  #ax1.boxplot(big_data[multiple_out_grp_names[-1]])
  #ax1.set_color_cycle(my_colors)
  #ax1.set_marker_cycle(my_markers)
  
    ax1.plot(bin_centers,hist_data, label='%s h' % aOutGroupName[4:])
#  for (i,outGrpName) in enumerate(multiple_out_grp_names[0:-1]):
#    list_of_avg=list()
#    patchList = list()
#    for aDataEntry in big_data[outGrpName]:
#      this_avg= np.average(aDataEntry)
#      print(this_avg)
#      list_of_avg.append(this_avg)
#    ax1.plot(list_of_avg, linestyle='None', marker=my_markers[i], label='%s h' % outGrpName[4:])
    
  #ax1.set_xticklabels(x_ticks_labels,rotation=75)
  #ax1.set_xticks(np.arange(len(x_ticks_labels))+1)
  #legendlabel = ['%s h' % entry[4:] for entry in multiple_out_grp_names]
  ax1.legend()
  if infos:
    ax1.set_xlabel(r' distance from center/ $\mu m$')
    ax1.set_ylabel(r' number of found cells')
      
    
    #ax1.set(title='file: %s \n at %s' % (goodArguments.vbl_simulation_output_filename, goodArguments.output_grp_name))
  ax1.grid(color='k', linestyle=':', linewidth=0.5)
  pp.savefig()
def create_file_name_multiple_out(out_groups_to_consider):
  string_to_append= ''
  for aGroup in out_groups_to_consider:
    string_to_append=string_to_append+'_%s' % aGroup[-4:]
  return string_to_append[1:]
def out_groups_to_hours(outGroups):
  return ['%s h' % aGroup[-3:] for aGroup in outGroups]
  
def get_dict_for_intervalls_by_nearest_vessel(filename, intervalls):
  with h5py.File(filename, 'r') as h5_f:
    h5_out_grp = h5_f[goodArguments.output_grp_name]
    distances_to_nearest_vessel = np.asarray(h5_out_grp['cells/distance_to_nearest_vessel'])
    distances_to_nearest_vessel = distances_to_nearest_vessel[:,0] 
    cell_phases= np.asarray(h5_out_grp['cells/cell_phase'])
    cell_phases = cell_phases[:,0]
    counts_per_intervall = list()
    for (i,aIntervall) in enumerate(intervalls):
      if i>0:
        good_indexes_in_intervall = np.where(np.logical_and(distances_to_nearest_vessel<aIntervall, distances_to_nearest_vessel > intervalls[i-1]))
      
        good_indexes_in_intervall = good_indexes_in_intervall[0]
        #cell_phases_in_this_intervall = cell_phases[good_indexes_in_intervall]
        array_of_1 = np.where(cell_phases[good_indexes_in_intervall] == 1)
        array_of_2 = np.where(cell_phases[good_indexes_in_intervall] == 2)
        array_of_3 = np.where(cell_phases[good_indexes_in_intervall] == 3)
        array_of_4 = np.where(cell_phases[good_indexes_in_intervall] == 4)
        array_of_5= np.where(cell_phases[good_indexes_in_intervall] == 5)
        array_of_6= np.where(cell_phases[good_indexes_in_intervall] == 6)
        dict_of_this_intervall = dict()
        dict_of_this_intervall['S'] = len(array_of_1[0]) + len(array_of_2[0])
        dict_of_this_intervall['G'] = len(array_of_3[0]) + len(array_of_4[0]) + len(array_of_5[0])
        dict_of_this_intervall['dead'] = len(array_of_6[0]) 
        counts_per_intervall.append(dict_of_this_intervall)
  return counts_per_intervall

def get_dict_for_intervalls_by_center(filename, intervalls):
  with h5py.File(filename, 'r') as h5_f:
    h5_out_grp = h5_f[goodArguments.output_grp_name]
    cell_centers = np.asarray(h5_out_grp['cells/cell_center_pos'])
    distances_to_center = np.sqrt(np.sum(np.power(cell_centers,2),1))
    #distances_to_nearest_vessel = np.asarray(h5_out_grp['cells/distance_to_nearest_vessel'])
    #distances_to_nearest_vessel = distances_to_nearest_vessel[:,0] 
    cell_phases= np.asarray(h5_out_grp['cells/cell_phase'])
    cell_phases = cell_phases[:,0]
    counts_per_intervall = list()
    for (i,aIntervall) in enumerate(intervalls):
      if i>0:
        good_indexes_in_intervall = np.where(np.logical_and(distances_to_center<aIntervall, distances_to_center > intervalls[i-1]))
      
        good_indexes_in_intervall = good_indexes_in_intervall[0]
        #cell_phases_in_this_intervall = cell_phases[good_indexes_in_intervall]
        array_of_1 = np.where(cell_phases[good_indexes_in_intervall] == 1)
        array_of_2 = np.where(cell_phases[good_indexes_in_intervall] == 2)
        array_of_3 = np.where(cell_phases[good_indexes_in_intervall] == 3)
        array_of_4 = np.where(cell_phases[good_indexes_in_intervall] == 4)
        array_of_5= np.where(cell_phases[good_indexes_in_intervall] == 5)
        array_of_6= np.where(cell_phases[good_indexes_in_intervall] == 6)
        dict_of_this_intervall = dict()
        dict_of_this_intervall['S'] = len(array_of_1[0]) + len(array_of_2[0])
        dict_of_this_intervall['G'] = len(array_of_3[0]) + len(array_of_4[0]) + len(array_of_5[0])
        dict_of_this_intervall['dead'] = len(array_of_6[0]) 
        counts_per_intervall.append(dict_of_this_intervall)
  return counts_per_intervall
  
def plot_cell_phases(s_by_intervall,g_by_intervall,dead_by_intervall,intervalls,pp):
  fig1 = plt.figure()
  ax1 = fig1.add_subplot(111)
  #ax1.scatter(distances_to_nearest_vessel[0:1000], endity_value_of_cells[0:1000])
  def draw_plot(data, edge_color, fill_color):
    bp = ax1.boxplot(data, patch_artist=True)

    for element in ['boxes', 'whiskers', 'fliers', 'means', 'medians', 'caps']:
        plt.setp(bp[element], color=edge_color)

    for patch in bp['boxes']:
        patch.set(facecolor=fill_color)
#  ax1.boxplot(s_by_intervall.transpose())
#  ax1.boxplot(g_by_intervall.transpose())
#  ax1.boxplot(dead_by_intervall.transpose())
  s_phase_color = '#FFD700'
  g_phase_color = 'red'
  dead_phase_color = 'green'
  draw_plot(s_by_intervall.transpose(), 'k' , s_phase_color)
  draw_plot(g_by_intervall.transpose(), g_phase_color, g_phase_color)
  draw_plot(dead_by_intervall.transpose(), dead_phase_color, dead_phase_color)
  labels = []
  for (i,inter) in enumerate(intervalls):
    if i < len(intervalls[:-1]):
      labels.append('%i - %i' %(inter, intervalls[i+1]))
  ax1.set_xticklabels(labels,rotation=75)
  print(len(s_by_intervall.transpose()))
  s_patch = matplotlib.patches.Patch(color=s_phase_color,       label=r's- phase # $%i \pm %i$' % (np.mean(np.sum(s_by_intervall,0)),np.std(np.sum(s_by_intervall,0))))
  g_patch = matplotlib.patches.Patch(color=g_phase_color,       label='g- phase # $%i \pm %i$' % (np.mean(np.sum(g_by_intervall,0)),np.std(np.sum(g_by_intervall,0))))
  dead_patch = matplotlib.patches.Patch(color=dead_phase_color, label='dead       # $%i \pm %i$' % (np.mean(np.sum(dead_by_intervall,0)),np.std(np.sum(dead_by_intervall,0))))
  ax1.legend(handles=[s_patch, g_patch, dead_patch])
  pp.savefig(fig1)
  
def create_matrix_nearest_vessels(list_of_filenames, max_distance_to_vessel=max_distance_to_vessel):
  list_by_file_name=dict()  
  '''find max distance 
    only if max_distance_to_vessel == 0.0
    otherwise we use the give value
  '''
  if max_distance_to_vessel == 0.0:
    for aFilename in list_of_filenames:
      with h5py.File(goodArguments.vbl_simulation_output_filename, 'r') as h5_f:
        h5_out_grp = h5_f[goodArguments.output_grp_name]
        distances_to_nearest_vessel = np.asarray(h5_out_grp['cells/distance_to_nearest_vessel'])
        distances_to_nearest_vessel = distances_to_nearest_vessel[:,0]
        
        if(np.max(distances_to_nearest_vessel)>max_distance_to_vessel):
          max_distance_to_vessel = np.max(distances_to_nearest_vessel)
          
  intervalls = np.arange(0,max_distance_to_vessel,intervall_for_plot_numbers )


  for aFilename in list_of_filenames:
    list_for_that_file = get_dict_for_intervalls_by_nearest_vessel(aFilename, intervalls)
    
    list_by_file_name[aFilename] = list_for_that_file
  
  s_by_intervall = np.zeros([len(intervalls), len(list_of_filenames)])
  g_by_intervall = np.zeros([len(intervalls), len(list_of_filenames)])
  dead_by_intervall = np.zeros([len(intervalls), len(list_of_filenames)])
  
  k=0
  for aFilename in list_of_filenames:
    i=0
    for aIntervall in intervalls[0:-1]:
      s_by_intervall[i,k] = list_by_file_name[aFilename][i]['S']
      g_by_intervall[i,k] = list_by_file_name[aFilename][i]['G']
      dead_by_intervall[i,k] = list_by_file_name[aFilename][i]['dead']
      i=i+1
    k=k+1
  return s_by_intervall, g_by_intervall, dead_by_intervall, intervalls

def create_matrix_center(list_of_filenames, max_distance_to_center=None):
  list_by_file_name=dict()  
  '''find max distance 
    only if max_distance_to_vessel == 0.0
    otherwise we use the give value
  '''
  if not max_distance_to_center:
    max_distance_to_center = 0.0
  for aFilename in list_of_filenames:
    with h5py.File(goodArguments.vbl_simulation_output_filename, 'r') as h5_f:
      h5_out_grp = h5_f[goodArguments.output_grp_name]
      cell_centers = np.asarray(h5_out_grp['cells/cell_center_pos'])
      distances_to_center = np.sqrt(np.sum(np.power(cell_centers,2),1))
      
      if(np.max(distances_to_center)>max_distance_to_center):
        max_distance_to_center = np.max(distances_to_center)
          
  intervalls = np.arange(0,max_distance_to_center,intervall_for_plot_numbers )
  print('max_distance_to_center')
  print(max_distance_to_center)

  for aFilename in list_of_filenames:
    list_for_that_file = get_dict_for_intervalls_by_center(aFilename, intervalls)
    
    list_by_file_name[aFilename] = list_for_that_file
  
  s_by_intervall = np.zeros([len(intervalls), len(list_of_filenames)])
  g_by_intervall = np.zeros([len(intervalls), len(list_of_filenames)])
  dead_by_intervall = np.zeros([len(intervalls), len(list_of_filenames)])
  
  k=0
  for aFilename in list_of_filenames:
    i=0
    for aIntervall in intervalls[0:-1]:
      s_by_intervall[i,k] = list_by_file_name[aFilename][i]['S']
      g_by_intervall[i,k] = list_by_file_name[aFilename][i]['G']
      dead_by_intervall[i,k] = list_by_file_name[aFilename][i]['dead']
      i=i+1
    k=k+1
  return s_by_intervall, g_by_intervall, dead_by_intervall, intervalls

if __name__ == '__main__':
  import argparse
  parser = argparse.ArgumentParser(description='Analyze MTS distances')  
  parser.add_argument('--s',dest='vbl_simulation_output_filename', type=str, default='None', help='output file name in hdf5 format')
  parser.add_argument('--g',dest='output_grp_name', type=str, default='out0001', help='output group withing hdf5 file')
  interactive = False;
  infos = True;
  goodArguments, otherArguments = parser.parse_known_args()
  
  if goodArguments.vbl_simulation_output_filename == 'None':
    print('no --s option found')
    list_of_filenames=['/localdisk/output/center_with_AS_latest/3489294/fakeTumMTS-default-typeI-sample00-vbl_safe_1.h5','/localdisk/output/center_with_AS_latest/3489295/fakeTumMTS-default-typeI-sample00-vbl_safe_1.h5']
  else:
    list_of_filenames=[goodArguments.vbl_simulation_output_filename]
  
  
  
    
  #print('doing for upper bount %s' % aIntervall)
  ''' begin of code '''
  '''register a clases at data manager'''
  print(os.path.basename(goodArguments.vbl_simulation_output_filename))
  
    
  with PdfPages('analysisCell_phase_dist_from_nearestvessel_%s.pdf' % goodArguments.output_grp_name ) as pp:
    s_by_intervall, g_by_intervall, dead_by_intervall, intervalls = create_matrix_nearest_vessels(list_of_filenames)
    plot_cell_phases(s_by_intervall,g_by_intervall,dead_by_intervall,intervalls,pp)
  with PdfPages('analysisCell_phase_dist_from_center_%s.pdf' % goodArguments.output_grp_name ) as pp:
    s_by_intervall, g_by_intervall, dead_by_intervall, intervalls = create_matrix_center(list_of_filenames, max_distance_to_center=422.0)
    plot_cell_phases(s_by_intervall,g_by_intervall,dead_by_intervall,intervalls,pp)  
 