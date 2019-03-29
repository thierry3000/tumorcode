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

def whiskers_plot_parameters(f,pdf):
  print('do whiskers_plot')
  all_bests = []
  convergent_groups = find_convergent_groups(f)
  for group in convergent_groups:
    inputFileName = f[group].attrs.get('vfile')
    vessel_grp = f[group].attrs.get('vessel_grp')
    setName = f[group].attrs.get('params')
    paramSetIndex = f[group].attrs.get('paramListIndex')
    #print("redo:")
    #print("file: %s,\t group: %s,\t set: %s,\t index: %i\n"% (inputFileName,vessel_grp,setName,paramSetIndex))
    factory = getattr(parameterSetsAdaption, setName)
    factory = factory[paramSetIndex]
    all_bests.append(np.asarray(f[group + '/best']))
  all_bests= np.asarray(all_bests)
  average_k_c = np.median(all_bests[:,0])
  average_k_m = np.median(all_bests[:,1])
  average_k_s = np.median(all_bests[:,2])
  print("averge kc: %f, average k_m: %f, average k_s %f" % (average_k_c, average_k_m,average_k_s))
  #f['/'].attrs['average_kc'] = average_k_c
  #f['/'].attrs['average_km'] = average_k_m
  #f['/'].attrs['average_ks'] = average_k_s
  ''' whiskers plot '''
  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.boxplot(all_bests)
  plt.tick_params(axis='both', which='major', labelsize=tickfontsize)
  #plt.tick_params(axis='both', which='minor', labelsize=8)
  plt.setp(ax, xticklabels=[r'$k_c$', r'$k_m$', r'$k_s$'])
  ax.text(1.,average_k_c, "%0.2f" % average_k_c,horizontalalignment='center',verticalalignment='bottom',fontsize=tickfontsize)
  ax.text(2.,average_k_m, "%0.2f" % average_k_m,horizontalalignment='center',verticalalignment='bottom',fontsize=tickfontsize)
  ax.text(3.,average_k_s, "%0.2f" % average_k_s,horizontalalignment='center',verticalalignment='bottom',fontsize=tickfontsize)
  total_groups = len(f.keys())
  ax.text(1.3,3.8,'%i convergent boundary conditions of\n%i tested' % (len(all_bests), total_groups),fontsize=12)
  #ax.xaxis.set_ticks_([r'$k_c$', r'$k_m$', r'$S_5$'])
  #box1 = plt.boxplot(all_bests)
  #plt.xticks([1,2,3],[r'$k_c$', r'$k_m$', r'$S_5$'])
  pdf.savefig()
  
def convergenz_phase_diagram(f,pdf):
  group_labels = []
  group_bests = []
  group_fitnesses = []
  param_id_label = []
  param_id = []
  nameOfParamSet = str(f['0'].attrs.get('params'))
  parameterDict = getattr(parameterSetsAdaption, nameOfParamSet)
  allpressures = []
  allflows = [] 
  for aParamDict in parameterDict:
    allpressures.append(aParamDict['adaption']['a_pressure'])
    allflows.append(aParamDict['adaption']['a_flow'])
  flows = np.sort(list(set(allflows)))
  flows = -1 * flows[::-1]
  pressures = np.sort(list(set(allpressures)))
  FLOW,PRESSURE = np.meshgrid(flows,pressures)
  #create map
  flow_to_id = dict()
  for (i,flow) in enumerate(flows):
    flow_to_id[flow] = i
  pressure_to_id = dict()
  for (i,pressure) in enumerate(pressures):
    pressure_to_id[pressure] = i
  no_of_flow = len(flows)
  no_of_pressures = len(pressures)
  Z = np.zeros([no_of_flow,no_of_pressures])
  for group in f.keys():
    print(group)
    group_labels.append(group)
    print( f[group + '/best'] )
    group_bests.append(np.asarray(f[group + '/best']))
    print( f[group + '/fitness values'] )
    fitness_value = np.asarray(f[group + '/fitness values'])[0]
    if fitness_value<0:
      group_fitnesses.append(-1*fitness_value)
      doLogPlot = False
    else:
      group_fitnesses.append(fitness_value)
    paramsetListIndex = int(f[group].attrs.get('paramListIndex'))
    current_flow=parameterDict[paramsetListIndex]['adaption']['a_flow']
    current_flow= -1*current_flow
    current_pressure= parameterDict[paramsetListIndex]['adaption']['a_pressure']
    Z[flow_to_id[current_flow],pressure_to_id[current_pressure]] = fitness_value
    print(paramsetListIndex)
    param_id.append(paramsetListIndex)
    param_id_label.append('%i %s' % (paramsetListIndex , group[0:5]) )
  group_fitnesses_mapped = np.zeros(len(group_fitnesses))
  for i, _ in enumerate(group_fitnesses):
    group_fitnesses_mapped[i] = group_fitnesses[i]
  #print(group_fitnesses_mapped)
  #print(param_id)
  '''contour plot'''
  plt.figure()
  max_Z=np.max(Z)
  bad_indeces = Z==max_Z
  sorted_fitnesses = np.sort(Z.flatten())
  sorted_fitnesses_smaller_max = sorted_fitnesses[sorted_fitnesses<max_Z]
  largest_value_except_false = sorted_fitnesses_smaller_max[-1]
  Z[bad_indeces] = largest_value_except_false
  CS = plt.contour(FLOW,PRESSURE,np.log10(np.transpose(Z)))
  CS.ax.tick_params(labelsize=tickfontsize)
  #CS = plt.contour(np.log10(FLOW),PRESSURE,np.log10(np.transpose(Z)))
  plt.clabel(CS, inline=1, fontsize=10)
  #plt.title('Convergent boundary conditions for adaption')
  #plt.xlabel(r'$\log_{10}$ (Flow through arterial inlet) / $\mu m^3/s$')
  plt.xlabel(r'Flow through arterial inlet / $\mu$m$^3$/s',fontsize=tickfontsize)
  plt.ylabel(r'Pressure at venous outlet/ kPa',fontsize=tickfontsize)
  labels = [r'log$_{10}$ of fitness']
  for i in range(len(labels)):
    CS.collections[i].set_label(labels[i])
  plt.legend(loc='lower right',fontsize=20)
  pdf.savefig()
  ''' smothed image'''
  if 0:
    plt.figure()
    im = plt.imshow(np.log10(np.transpose(Z)),interpolation='bilinear')
    plt.colorbar()
    #print("finshed that")
    #plt.xticks(param_id, param_id_label, rotation='vertical')
    #plt.tight_layout()
    #ax = plt.gca()
    #if doLogPlot:
    #  ax.set_yscale('log')
    pdf.savefig()
    print("finished phase diagram")
  
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
  variance_stuff = False
  for group in f.keys():
    if variance_stuff:
      convergent_groups.append(group)
    else:
      #if float(np.asarray(f[group + '/fitness values'])) < 1e15:
      if float(np.asarray(f[group + '/fitness values'])) < np.power(10,4.8):
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
    returnState, mean, varOfMean, total_surface = _adap.do_simple_adaption(inputFileName,vessel_grp, factory)
    #f[group].attrs['returnState'] = returnState
    #f[group].attrs['mean'] = mean
    #f[group].attrs['varOfMean'] = varOfMean
    #f[group].attrs['total_surface'] = total_surface
    #note: this is not flushing the file, so wait until program finished!
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
  
  basenameOfFile = basename(goodArguments.fileName)
  tickfontsize=14
  rc = matplotlib.rc
  rc('font', size = 8.)
  rc('axes', titlesize = 10., labelsize = 8.)
  
  with h5py.File(goodArguments.fileName,'r') as f:
    if 0:
      with PdfPages('out_deap_results_convergenz_plot_%s.pdf' % basenameOfFile) as pdf:
        convergenz_plot(f,pdf)
    if 1:
      with PdfPages('out_deap_results_phase_dia_%s.pdf' % basenameOfFile) as pdf:
        convergenz_phase_diagram(f,pdf)
    if 1:
      with PdfPages('out_deap_results_whiskers_params_%s.pdf' % basenameOfFile) as pdf:
        whiskers_plot_parameters(f,pdf)
      
    if goodArguments.redo:
      redo_adaption_for_convergent_sets(f)