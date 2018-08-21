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

import matplotlib.pyplot as plt

eps_tube = 20

def get_sample_points_along_line(pos_a, pos_b, number_of_samples=2):
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
def is_index_good(sample_pos):
  array_of_cell_centers=scoop.shared.getConst('cell_center_pos_')
  list_of_good_indeces=[]
  print(array_of_cell_centers.shape)
  for (i,aPos) in enumerate(array_of_cell_centers[:,0:3]):
    #print('aPos')
    #print(aPos)
    distance_from_sample_pos_to_a_cell_center_pos_i = np.linalg.norm(np.asarray(aPos) -sample_pos, axis=0)
    if(distance_from_sample_pos_to_a_cell_center_pos_i<eps_tube):
      list_of_good_indeces.append(i)
  return list_of_good_indeces
      
def sample_line_artery(**keywords):
  artery_p1 = np.asarray([216,459,-324])
  artery_p2 = np.asarray([223,383,-370])
  
  starting_pos=artery_p1
  end_pos=artery_p2
  
  print(goodArguments.vbl_simulation_file_name)
  print(goodArguments.grp_pattern)
  
  print("sample along line from:")
  print("%s to %s" % (starting_pos, end_pos))
  f = h5py.File(goodArguments.vbl_simulation_file_name, 'r')
  vesselgroup = f[os.path.join(goodArguments.grp_pattern, 'vessels')]
  graph = krebsutils.read_vessels_from_hdf(vesselgroup, ['position', 'radius', 'hematocrit', 'pressure', 'flow', 'flags','shearforce'] + datalist, return_graph=True)
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
  sample_pos, factors = get_sample_points_along_line(artery_p1,artery_p2,number_of_samples=30)
  print('found sample positions')
  cell_center_pos = np.asarray(f[os.path.join(goodArguments.grp_pattern, 'cells/cell_center_pos')])
  print('cell_center_pos')
  scoop.shared.setConst(cell_center_pos_=cell_center_pos)
  lists = list(scoop.futures.map(is_index_good, sample_pos))
  for (aList, aSamplePos) in zip(lists,sample_pos):
    print('for pos: %s found %i points in eps range' % (aSamplePos,len(aList)))
    
  ''' average '''
  average_value = []
  errors = []
  for (aList, aSamplePos) in zip(lists,sample_pos):
    this_avg = np.average(aList)
    average_value.append(this_avg)
    errors.append(np.sqrt(1/float(len(aList)))*this_avg)
    
#  print('average values:')
#  print(average_value)
  
  
  print('finished sample line artery')
  return (np.asarray(average_value), np.asarray(errors), np.asarray(factors))

if __name__ == '__main__':
  import argparse
  parser = argparse.ArgumentParser(description='plot cell o2 along a line')  
  parser.add_argument('vbl_simulation_file_name')
  parser.add_argument('grp_pattern')
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
  datalist = map(lambda s: s, map(str.strip, goodArguments.datalist.split(',')))
  
  avg_values, errors_avg, distances = sample_line_artery()
  print(errors_avg)
  #plt.scatter(sample_factor, avg_values)
  plt.errorbar(distances,avg_values, yerr=errors_avg)
  plt.show()