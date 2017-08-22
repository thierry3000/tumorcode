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
import os, sys
from os.path import join, basename, dirname, splitext
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'..'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../..'))
import h5py
import matplotlib
import mpl_utils
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy as np

import os, sys
from os.path import join, basename, dirname, splitext
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../..'))

from krebsjobs import parameters
from parameters import parameterSetsAdaption
from krebs import adaption as _adap

#f = h5py.File('deap_results.h5','r')

def convergenz_plot(f,pdf):
  group_labels = []
  group_bests = []
  group_fitnesses = []
  param_id_label = []
  param_id = []
  doLogPlot = True
  for group in f.keys():
    print(group)
    group_labels.append(group)
    print( f[group + '/best'] )
    group_bests.append(f[group + '/best'])
    print( f[group + '/fitness values'] )
    fitness_value = np.asarray(f[group + '/fitness values'])[0]
    if fitness_value<0:
      group_fitnesses.append(-1*fitness_value)
      doLogPlot = False
    else:
      group_fitnesses.append(fitness_value)
    print(f[group].attrs.get('paramListIndex'))
    param_id.append(int(f[group].attrs.get('paramListIndex')))
    param_id_label.append('%i %s' % (f[group].attrs.get('paramListIndex') , group[0:5]) )
  group_fitnesses_mapped = np.zeros(len(group_fitnesses))
  for i, _ in enumerate(group_fitnesses):
    group_fitnesses_mapped[i] = group_fitnesses[i]
  print(group_fitnesses_mapped)
  print(param_id)
  plt.scatter(param_id, group_fitnesses_mapped)
  print("finshed that")
  plt.xticks(param_id, param_id_label, rotation='vertical')
  plt.tight_layout()
  ax = plt.gca()
  if doLogPlot:
    ax.set_yscale('log')
  pdf.savefig()
  print("finished convergenz plot")
  #plt.show()
  
def find_convergent_groups(f):
  convergent_groups=[]
  variance_stuff = True
  for group in f.keys():
    if variance_stuff:
      convergent_groups.append(group)
    else:
      if float(np.asarray(f[group + '/fitness values'])) < 1e15:
        convergent_groups.append(group)
  return convergent_groups

def redo_adaption_for_convergent_sets(f):
  print("starting redo...")
  convergent_groups = find_convergent_groups(f)
  for group in convergent_groups:
    print("Found group: %s" % group)
    inputFileName = f[group].attrs.get('vfile')
    vessel_grp = f[group].attrs.get('vessel_grp')
    setName = f[group].attrs.get('params')
    paramSetIndex = f[group].attrs.get('paramListIndex')
    print("redo:")
    print("file: %s,\t group: %s,\t set: %s,\t index: %i\n"% (inputFileName,vessel_grp,setName,paramSetIndex))
    factory = getattr(parameterSetsAdaption, setName)
    factory = factory[paramSetIndex]
    bests = np.asarray(f[group + '/best'])
    print("bests found: %s " % bests)
    factory['adaption'].update(
        k_c = bests[0],
        k_m = bests[1],
        k_s = bests[2],
        outputFileName = 'redoAdaption_set_%s_%i.h5' % (setName,paramSetIndex) ,
        )
    _adap.do_simple_adaption(inputFileName,vessel_grp, factory)
if __name__ == '__main__':
  import argparse
  parser = argparse.ArgumentParser(description='redo adaption from deap optimization', formatter_class=argparse.ArgumentDefaultsHelpFormatter)  
  #parser.add_argument('AdaptionParamSet')
  #parser.add_argument('--listindex', type=int, help="index of value list" )
  #parser.add_argument('--outputFileFolder', type=str, help="where to store the output, default is working directory")
  parser.add_argument('--fileName', type=str,help="deap output")
  parser.add_argument('--redo', default=False,help="create the optimized network",action="store_true")
  goodArguments, otherArguments = parser.parse_known_args()
  #print("running with %s" % goodArguments.AdaptionParamSet)
  
  
  f = h5py.File(goodArguments.fileName,'r')
  basenameOfFile = basename(goodArguments.fileName)
  with PdfPages('out_deap_results_%s.pdf' % basenameOfFile) as pdf:
    rc = matplotlib.rc
    rc('font', size = 8.)
    rc('axes', titlesize = 10., labelsize = 8.)
    if 1:
      convergenz_plot(f,pdf)
  if goodArguments.redo:
    redo_adaption_for_convergent_sets(f)