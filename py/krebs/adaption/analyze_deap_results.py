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

import h5py
import matplotlib.pyplot as plt
import numpy as np


from krebsjobs.parameters import parameterSetsAdaption
from krebs import adaption as _adap

f = h5py.File('deap_results.h5','r')

def convergenz_plot(f):
  group_labels = []
  group_bests = []
  group_fitnesses = []
  param_id_label = []
  param_id = []
  for group in f.keys():
    print(group)
    group_labels.append(group)
    print( f[group + '/best'] )
    group_bests.append(f[group + '/best'])
    print( f[group + '/fitness values'] )
    group_fitnesses.append(np.asarray(f[group + '/fitness values'])[0])
    print(f[group].attrs.get('paramListIndex'))
    param_id.append(int(f[group].attrs.get('paramListIndex')))
    param_id_label.append('%i %s' % (f[group].attrs.get('paramListIndex') , group[0:5]) )
  
  group_fitnesses_mapped = np.zeros(len(group_fitnesses))
  for i, _ in enumerate(group_fitnesses):
    group_fitnesses_mapped[i] = group_fitnesses[i]
  plt.scatter(param_id, group_fitnesses_mapped)
  plt.xticks(param_id, param_id_label, rotation='vertical')
  plt.tight_layout()
  ax = plt.gca()
  ax.set_yscale('log')
  plt.show()
  
def find_convergent_groups(f):
  convergent_groups=[]
  for group in f.keys():
    if float(np.asarray(f[group + '/fitness values'])) < 1000:
      convergent_groups.append(group)
  return convergent_groups

def redo_adaption_for_convergent_sets(f):
  print("starting redo...")
  convergent_groups = find_convergent_groups(f)
  for group in convergent_groups:
    print("Found group: %s" % group)
    inputFileName = f[group].attrs.get('vfile')
    vessel_grp = f[group].attrs.get('vessel_grp')
    setName = "value_list"
    paramSetIndex = f[group].attrs.get('paramListIndex')
    factory = getattr(parameterSetsAdaption, setName)
    factory = factory[paramSetIndex]
    bests = np.asarray(f[group + '/best'])
    print("bests found: %s " % bests)
    factory['adaption'].update(
        k_c = bests[0],
        k_m = bests[1],
        k_s = bests[2],
        outputFileName = 'redoAdaption_value_list_%i.h5' % paramSetIndex ,
        )
    _adap.do_simple_adaption(inputFileName,vessel_grp, factory)
if __name__ == '__main__':
  f = h5py.File('deap_results.h5','r')
  if 1:
    convergenz_plot(f)
  if 0:
    redo_adaption_for_convergent_sets(f)