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
from __future__ import print_function  
import os, sys
from os.path import join, basename, dirname, splitext
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../..'))

import krebsutils
#import extractVtkFields
#import pprint
import numpy as np
import h5py
#import identifycluster
import myutils
import multiprocessing
import scoop
import matplotlib
matplotlib.use('agg')
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

eps_tube = 5


''' ARTERY
initial seed at: 240,400,-310
vessels at state 520
             299
           *    *
edge:492  *      * edge: 1215
         *        *
        *          *
        340        1357  
'''
artery_p1 = np.asarray([216,459,-324])
artery_p2 = np.asarray([223,383,-370])
artery_p1 = np.asarray([169,477,-306])
artery_p2 = np.asarray([189,352,-372])

''' VEIN
initial seed at: 360,180,-310
vessels at state 605
             364
           *    *
edge:553  *      * edge: 552
         *        *
        *          *
        3026      22  
'''

vein_parallel_p1 = np.asarray([350,150,-370])
vein_parallel_p2 = np.asarray([344,191,-262])
vein_parallel_p1 = np.asarray([349,157,-350])
vein_parallel_p2 = np.asarray([344,185,-281])

vein_ortho_p1 = np.asarray([350,248,-327])
vein_ortho_p2 = np.asarray([327,78,-283])

metadata_dict = dict()
metadata_dict['norm'] = 'maximum norm'
metadata_dict['eps_tube'] = eps_tube
metadata_dict['number_of_sampling_points'] = 40
'''
      this creates sample points along the line
      from starting_pos to end_pos
      returns the sampling points and their distance along the line
  '''
def get_sample_points_along_line(pos_a, pos_b, number_of_samples=2):
  metadata_dict['starting_pos'] = pos_a
  metadata_dict['end_pos'] = pos_b
  positions=[]
  factors = []
  distances = []
#  print(pos_a)
#  print(pos_b)
#  print(pos_b-pos_a)
#  print(number_of_samples)
  for i in range(number_of_samples):
#    print(i)
    factor= i/(float(number_of_samples)-1)
#    print("factor: %f" %factor)
    factors.append(factor)
    new_sample_pos = pos_a+ factor*(pos_b-pos_a)
    positions.append(new_sample_pos);
    distances.append(np.linalg.norm(np.asarray(new_sample_pos)- pos_a))
#  print(positions)
  return (np.asarray(positions), np.asarray(distances))
'''
  checks which indices are within eps_tube around the line
  
  returns the indeces an their distances
'''
def is_index_good(sample_pos):
  array_of_cell_centers=scoop.shared.getConst('cell_center_pos_')
  box_length_of_max_norm=scoop.shared.getConst('box_length_of_max_norm')
  list_of_good_indeces=[]
  #print(array_of_cell_centers.shape)
  ''' euclidean norm'''
#  for (i,aPos) in enumerate(array_of_cell_centers[:,0:3]):
#    #print('aPos')
#    #print(aPos)
#    distance_from_sample_pos_to_a_cell_center_pos_i = np.linalg.norm(np.asarray(aPos) -sample_pos, axis=0)
#    if(distance_from_sample_pos_to_a_cell_center_pos_i<eps_tube):
#      list_of_good_indeces.append(i)
#  return list_of_good_indeces
  ''' max norm'''
  for (i,aPos) in enumerate(array_of_cell_centers[:,0:3]):
    #print('aPos')
    #print(aPos)
    #distance_from_sample_pos_to_a_cell_center_pos_i = np.linalg.norm(np.asarray(aPos) -sample_pos, axis=0)
    distance_from_sample_pos_to_a_cell_center_pos_i = np.linalg.norm(np.asarray(aPos) -sample_pos, np.inf, axis=0)
    if(distance_from_sample_pos_to_a_cell_center_pos_i<eps_tube):
      list_of_good_indeces.append(i)
  return list_of_good_indeces
      
''' this is my development function
'''
def sample_line_artery(**keywords):
  starting_pos=artery_p1
  end_pos=artery_p2
  
  print(goodArguments.vbl_simulation_file_name)
  print(goodArguments.grp_pattern)
  
  print("sample along line from:")
  print("%s to %s" % (starting_pos, end_pos))
  #f = h5py.File(goodArguments.vbl_simulation_file_name, 'r')
  #vesselgroup = f[os.path.join(goodArguments.grp_pattern, 'vessels')]
  #graph = krebsutils.read_vessels_from_hdf(vesselgroup, ['position', 'radius', 'hematocrit', 'pressure', 'flow', 'flags','shearforce'] + datalist, return_graph=True)
  ''' these are the interesting vessels for the artery'''
  interesting_vessels = [492, 1215]
  edges = graph.edgelist
  pos = graph['position']
  for vessel in interesting_vessels:
    print("edge: %i nodes: %s" %(vessel,edges[vessel]))
    this_edges = edges[vessel]
    print("node: %i, pos: %s" %(this_edges[0],pos[this_edges[0]]))
    print("node: %i, pos: %s" %(this_edges[1],pos[this_edges[1]]))
  print("pos 252: %s" % pos[252])
  print("pos 293: %s" % pos[293])
  print("pos 1259: %s" % pos[1259])
  
  sample_pos, distances = get_sample_points_along_line(starting_pos,end_pos,metadata_dict['number_of_sampling_points'])
  box_length_of_max_norm = distances[1]*0.5
  print(distances)
  scoop.shared.setConst(box_length_of_max_norm=box_length_of_max_norm)
  print('box_length_of_max_norm %f' % box_length_of_max_norm )
  #cell_center_pos = np.asarray(f[os.path.join(goodArguments.grp_pattern, 'cells/cell_center_pos')])
  #print('cell_center_pos')
  #scoop.shared.setConst(cell_center_pos_=cell_center_pos)
  lists = list(scoop.futures.map(is_index_good, sample_pos))
  for (aList, aSamplePos) in zip(lists,sample_pos):
    print('for pos: %s found %i points in eps range' % (aSamplePos,len(aList)))
    
  ''' average '''
  average_value = []
  errors = []
  
  for (aList, aSamplePos) in zip(lists,sample_pos):
    this_avg = np.average(quantity_to_average[aList])
    average_value.append(this_avg)
    #errors.append(np.sqrt(1/float(len(aList)))*this_avg)
    errors.append(np.std(quantity_to_average[aList]))
      
  print('finished sample line artery')
  return (np.asarray(average_value), np.asarray(errors), np.asarray(distances))

def sample_line_vein_parallel_bifurcation(**keywords):
  starting_pos=vein_parallel_p1
  end_pos=vein_parallel_p2
  
  edges = graph.edgelist
  sample_pos, factors = get_sample_points_along_line(starting_pos,end_pos,metadata_dict['number_of_sampling_points'])
  print('found sample positions')
  lists = list(scoop.futures.map(is_index_good, sample_pos))
  for (aList, aSamplePos) in zip(lists,sample_pos):
    print('for pos: %s found %i points in eps range' % (aSamplePos,len(aList)))
    
  ''' average '''
  average_value = []
  errors = []
  for (aList, aSamplePos) in zip(lists,sample_pos):
    this_avg = np.average(quantity_to_average[aList])
    average_value.append(this_avg)
    errors.append(np.sqrt(1/float(len(aList)))*this_avg)
      
  print('finished sample line artery')
  return (np.asarray(average_value), np.asarray(errors), np.asarray(factors))

def sample_line_vein_orthogonal_bifurcation(**keywords):
  starting_pos=vein_ortho_p1
  end_pos=vein_ortho_p2
  
  sample_pos, factors = get_sample_points_along_line(starting_pos,end_pos,metadata_dict['number_of_sampling_points'])
  print('found sample positions')

  lists = list(scoop.futures.map(is_index_good, sample_pos))
  for (aList, aSamplePos) in zip(lists,sample_pos):
    print('for pos: %s found %i points in eps range' % (aSamplePos,len(aList)))
    
  ''' average '''
  average_value = []
  errors = []
  for (aList, aSamplePos) in zip(lists,sample_pos):
    this_avg = np.average(quantity_to_average[aList])
    average_value.append(this_avg)
    errors.append(np.sqrt(1/float(len(aList)))*this_avg)
      
  print('finished sample line artery')
  return (np.asarray(average_value), np.asarray(errors), np.asarray(factors))

if __name__ == '__main__':
  import argparse
  parser = argparse.ArgumentParser(description='plot cell o2 along a line')  
  parser.add_argument('vbl_simulation_file_name')
  parser.add_argument('grp_pattern')
  parser.add_argument('type', help="seed of initial spheroid: v: vein or a: artery")
  parser.add_argument("-d","--data", dest="datalist", help="which data (pressure, flow, shearforce, hematocrit) as comma separated list", default='pressure', action="store")
  parser.add_argument("-f","--filter-uncirculated", dest="filteruncirculated", help="filter uncirculated vessels", default=False, action="store_true")
  parser.add_argument("--filter-radius-high-pass", dest="filterradiushighpass", type=float, default=-1.0)
  parser.add_argument("--no-overlay", dest="overlay", default = True, action="store_false")
  parser.add_argument("--dpi", dest="dpi", default=None, action="store")
  parser.add_argument("--format", dest="format", default=None, action="store")
  parser.add_argument("-n", dest="out_nl", default = False, action="store") 
  parser.add_argument("--outFilename", dest="outfn", default= None, type=str)
  '''this option come from the tumor side'''
  parser.add_argument("--writeVessels", help="when doing the tumor, export vesesls", default=True, action="store_true")  
  parser.add_argument("--writeFields", help="when doing the tumor, export Fields", default=None)  
  goodArguments, otherArguments = parser.parse_known_args()
  
  '''check input'''
  try:
    dirs = set()
    if not os.path.isfile(goodArguments.vbl_simulation_file_name):
      raise AssertionError('The file %s is not present!'%goodArguments.vbl_simulation_file_name)
    with h5py.File(goodArguments.vbl_simulation_file_name, 'r') as f:
      d = myutils.walkh5(f, goodArguments.grp_pattern)
      if not len(d)>0:
        raise AssertionError('pattern "%s" not found in "%s"!' % (goodArguments.grp_pattern, goodArguments.vbl_simulation_file_name))
      else:
        dirs = set.union(dirs,d)
  except Exception, e:
    print(e.message)
    sys.exit(-1)
    
    
  ''' begin of code '''
  f = h5py.File(goodArguments.vbl_simulation_file_name, 'r')
  datalist = map(lambda s: s, map(str.strip, goodArguments.datalist.split(',')))
  vesselgroup = f[os.path.join(goodArguments.grp_pattern, 'vessels')]
  graph = krebsutils.read_vessels_from_hdf(vesselgroup, ['position', 'radius', 'hematocrit', 'pressure', 'flow', 'flags','shearforce'] + datalist, return_graph=True)
  cell_center_pos = np.asarray(f[os.path.join(goodArguments.grp_pattern, 'cells/cell_center_pos')])
  cell_o2_mass = np.asarray(f[os.path.join(goodArguments.grp_pattern, 'cells/o2')])
  cell_radii = np.asarray(f[os.path.join(goodArguments.grp_pattern, 'cells/cell_radii')])
  # pg/ mum^3
  cell_o2_concentration = cell_o2_mass/ (4/float(3)* np.pi*np.power(cell_radii,3))
  #cell_o2_concentration = cell_o2_mass/ np.power(eps_tube,3)
  volume_o2_ml = cell_o2_mass*1e-9/1.429
  ''' o2 density 1.429 g/L --> 1.429*10^9 pg/ml
  '''
  solubility = 2.8e-4 #ml O2/cm^3 mmHg
  solubility = solubility*1e-12 #ml O2/mum^3 mmHg
  volume_density = 1.429e9 #pg/ml
  x = cell_o2_concentration/volume_density # ml / mum^3
  
  cell_po2 = x/solubility
  
  #quantity_to_average = cell_o2_concentration
  quantity_to_average = volume_o2_ml/solubility
  #quantity_to_average = cell_o2_mass
  print('cell_center_pos')
  scoop.shared.setConst(cell_center_pos_=cell_center_pos)
  ''' pdf output '''
  meta_data_fig = plt.figure(figsize=(11.69,8.27))
  meta_data_fig.clf()
  result_string = ''

  infos = False  
  
  if goodArguments.type == 'a':
    avg_values, errors_avg, distances = sample_line_artery()
    fig1, ax1 = plt.subplots(1)
    ax1.errorbar(distances,avg_values, yerr=errors_avg)
    if infos:
      ax1.set(title = 'Cell based oxygen along parallel line \n at arterial bifurcation')
    ax1.set_xlabel(r'distance along line/ $\mu m$')
    ax1.set_ylabel('pO2/mmHg')
    ax1.grid(color='k', linestyle=':', linewidth=0.5)
    ax1.set_xlim([-10, 150])    
    with PdfPages('arterial.pdf') as pp:
      for entry in metadata_dict:
        #print(entry)
        result_string+= '%s:\t%s\n' % (entry, metadata_dict[entry])
      meta_data_fig.text(0.1,0.1,result_string, transform=meta_data_fig.transFigure, size=14, ha="left")
      pp.savefig(fig1)
      if infos:
        pp.savefig(meta_data_fig)
    
  if goodArguments.type == 'v':
    with PdfPages('venous_parallel.pdf') as pp:
      avg_values, errors_avg, distances = sample_line_vein_parallel_bifurcation()
      fig1, ax1 = plt.subplots(1)
      ax1.errorbar(distances,avg_values, yerr=errors_avg)
      if infos:
        ax1.set(title = 'Cell based oxygen along parallel line \n at venous bifurcation')
      ax1.set_xlabel(r'distance along line/ $\mu m$')
      ax1.set_ylabel('pO2/mmHg')
      ax1.grid(color='k', linestyle=':', linewidth=0.5)
      
      for entry in metadata_dict:
        #print(entry)
        result_string+= '%s:\t%s\n' % (entry, metadata_dict[entry])
      meta_data_fig.text(0.1,0.1,result_string, transform=meta_data_fig.transFigure, size=14, ha="left")
      pp.savefig(fig1)
      if infos:
        pp.savefig(meta_data_fig)
      meta_data_fig.clf()
    
    
    with PdfPages('venous_ortho.pdf') as pp:
      avg_values, errors_avg, distances = sample_line_vein_orthogonal_bifurcation()
      fig2, ax2 = plt.subplots(1)
      ax2.errorbar(distances,avg_values, yerr=errors_avg)
      if infos:
        ax2.set(title = 'Cell based oxygen along orthogonal line \n at venous bifurcation')
      ax2.set_xlabel(r'distance along line/ $\mu m$')
      ax2.set_ylabel('pO2/mmHg')
      ax2.grid(color='k', linestyle=':', linewidth=0.5)
      for entry in metadata_dict:
        #print(entry)
        result_string+= '%s:\t%s\n' % (entry, metadata_dict[entry])
      meta_data_fig.text(0.1,0.1,result_string, transform=meta_data_fig.transFigure, size=14, ha="left")
      pp.savefig(fig2)
      if infos:
        pp.savefig(meta_data_fig)
    
  #print(errors_avg)
  #plt.scatter(sample_factor, avg_values)
  
  #plt.title("Cell based oxygen along line at bifurcation")
  #plt.xlabel(r"distance along line/ $\mu m$")
  #plt.ylabel("")
  #plt.show()