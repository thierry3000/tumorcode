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

"""
        visualize runtime behaviour of MTS simulation
"""

import h5py
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def plot_vbl_runtime(goodArguments,pp):
  f = open(goodArguments.vbl_timing_filename,'r')
  lines = f.readlines()
  diff = lines[0::3]
  cellEvents = lines[1::3]
  geometry = lines[2::3]
  diff = [aDiff[aDiff.find('for ')+3:aDiff.find('ms')] for aDiff in diff]
  cellEvents = [aCellEvent[aCellEvent.find('for ')+3:aCellEvent.find('ms')] for aCellEvent in cellEvents]
  geometry = [aGeomentryEvent[aGeomentryEvent.find('for ')+3:aGeomentryEvent.find('ms')] for aGeomentryEvent in geometry]
  
  diff = np.asarray(diff,dtype=float)
  cellEvents = np.asarray(cellEvents,dtype=float)
  geometry = np.asarray(geometry,dtype=float)
  
  diff = [adiff/1000000.0 for adiff in diff]
  cellEvents = [acellEvents/1000000.0 for acellEvents in cellEvents]
  geometry = [ageomentry/1000000.0 for ageomentry in geometry]
  
  fig, ax = plt.subplots()
  plt.plot(diff, label='diff call')
  #plt.plot(cellEvents, label='cell Events')
  plt.plot(geometry, label='geometry call')
  plt.ylabel('runtime/ s')
  plt.xlabel('calls at day 3')
  #ax.set_yscale("log", nonposy='clip')
  
  # Now add the legend with some customizations.
  legend = ax.legend(loc='upper center', shadow=True)
  
  if interactive:
    plt.show()
  else:
    pp.savefig()

def read_data(filename):
  timePoints = []
  noOfCells=[]
  with h5py.File(filename, 'r') as f:
    initialTime = f['out0000'].attrs.get('secondsSinceEpoch')
    print("Found inital time: %i\n" % initialTime)
    for key in f.keys():
      if 'out' in key:
        ex = key+'/cells/cell_radii' in f
        timePoints.append(f[key].attrs.get('secondsSinceEpoch'))
        if ex:
          noOfCells.append(len(f[key+'/cells/cell_radii']))
        else:
          noOfCells.append(0)
  timePoints = np.asarray(timePoints)
  timePoints = timePoints - initialTime
  timePoints = timePoints/3600
  noOfCells = np.asarray(noOfCells)
  return timePoints,noOfCells
def plot_no_cells_over_time_since_epoch(goodArguments,pp):
  timePointsInHours_28, noOfCells_28 = read_data('safe_28.h5')
  timePointsInHours_14, noOfCells_14 = read_data('safe_14.h5')
  timePointsInHours_7, noOfCells_7 = read_data('safe_7.h5')
  
  fig, ax = plt.subplots()
  ax.plot(timePointsInHours_28, noOfCells_28,label='28 threads')
  ax.plot(timePointsInHours_14, noOfCells_14,label='14 threads')
  ax.plot(timePointsInHours_7, noOfCells_7,label='7 threads')
  plt.xlabel('time/h')
  plt.ylabel('#cells')
  legend = ax.legend(loc='lower right', shadow=True)
  
  if interactive:
    plt.show()
  else:
    pp.savefig()
def plot_runtime_from_h5(goodArguments, pp):
  myH5Keys = [
              #'run_ann',
              'run_doMilottiStep',
              #'run_doStep',
              'run_o2',
              #'run_vbl_cellEvents',
              'run_vbl_diff',
              'run_vbl_diff_loop_1',
              'run_vbl_diff_loop_2',
              'run_vbl_diff_loop_3',
              'run_vbl_dynamics',
              'run_vbl_geometry',
              #'run_vbl_bico_call',
              #'run_vbl_writeToFile',
              'total_time', # NOTE: this is the integrated time, all others are runtimes!
              'number_of_cells',
              ]
  myDictOfData = dict()
  for key in myH5Keys:
    myDictOfData[key] = []
    
  no_of_threads = 42
  print("filename: %s" % str(goodArguments.vbl_simulation_output_filename))
  with h5py.File(str(goodArguments.vbl_simulation_output_filename), 'r') as f:
    initialTime = f['out0000/timing'].attrs.get('secondsSinceEpoch')
    no_of_threads=f.attrs.get('detectedNumberOfThreads')
    lastKey = 'out0000'
#    def isKeyGood(key):
#      if 'out' in key:
#        return str(key)
    goodKeys = [str(x) for x in f.keys() if 'out' in x]
    for mykey in myH5Keys:
      for (i,key) in enumerate(goodKeys):
        if mykey == 'total_time':
          if i ==0:
            totalTime = 0
          else:
            totalTime = int(f[key+'/timing'].attrs.get('secondsSinceEpoch')) - int(f[goodKeys[i-1]+'/timing'].attrs.get('secondsSinceEpoch'))
          myDictOfData[mykey].append(totalTime)
        elif mykey == 'number_of_cells':
          if 'out' in key:
            ex = key+'/cells/cell_radii' in f
            if ex:
              arrayOfRadii = np.asarray(f[key+'/cells/cell_radii'])
              myDictOfData[mykey].append(len(np.concatenate(arrayOfRadii, axis=0)))
            else:
              myDictOfData[mykey].append(0)
        else:
          #myDictOfData[mykey].append(f[key+'/timing'].attrs.get(mykey)/1e6)
          timeForKey = f[key+'/timing'].attrs.get(mykey)
          myDictOfData[mykey].append(timeForKey[0][0])
    b_sytem_parameters_exist = '/parameters/system' in f
    if b_sytem_parameters_exist:
      no_of_threads = f['/parameters/system'].attrs.get('num_threads')
  
  plotVsCells = True;
  fig, ax = plt.subplots()
  for mykey in myH5Keys:
    print("mykey: %s" % mykey)
    if not mykey == 'number_of_cells':
      if plotVsCells:
          ax.plot(myDictOfData['number_of_cells'],myDictOfData[mykey][0:],label=mykey)
      else:
          ax.plot(myDictOfData[mykey][0:],label=mykey)
          
  
  if plotVsCells:
    ax.set_xlabel('#cells')
    ax.set_ylabel('runtime/ s')
  else:
    ax.set_xlabel('output group/ one unit is 1h of simulated time')
    ax.set_ylabel('runtime/ s')
  
  
  ax.set_title('cluster run of faketumor with vbl \n run on snowden with # %i threads' % int(no_of_threads))
  
  legend = ax.legend(loc='upper left', shadow=True)
  if interactive:
    plt.show()
  else:
    pp.savefig()
    
def plot_memory_from_h5(goodArguments, pp):
  myH5Keys = [
              #'run_ann',
              #'run_doMilottiStep',
              #'run_doStep',
              #'run_o2',
              #'run_vbl_cellEvents',
              #'run_vbl_diff',
              #'run_vbl_diff_loop_1',
              #'run_vbl_diff_loop_2',
              #'run_vbl_diff_loop_3',
              #'run_vbl_dynamics',
              #'run_vbl_geometry',
              #'run_vbl_bico_call',
              #'run_vbl_writeToFile',
              #'total_time', # NOTE: this is the integrated time, all others are runtimes!
              'number_of_cells',
              ]
  myDictOfData = dict()
  for key in myH5Keys:
    myDictOfData[key] = []
  
  h5_memory_keys = [
        'rss',
        'rss_peak',
        'vmem',
        'vmem_peak'
        ]
  myDictOfDataForMemory = dict()
  for key in h5_memory_keys:
    myDictOfDataForMemory[key] =[]
    
  no_of_threads=0
  print("filename: %s" % str(goodArguments.vbl_simulation_output_filename))
  with h5py.File(str(goodArguments.vbl_simulation_output_filename), 'r') as f:
    initialTime = f['out0000/timing'].attrs.get('secondsSinceEpoch')
    no_of_threads=f.attrs.get('detectedNumberOfThreads')
    lastKey = 'out0000'
#    def isKeyGood(key):
#      if 'out' in key:
#        return str(key)
    goodKeys = [str(x) for x in f.keys() if 'out' in x]
    for mykey in myH5Keys:
      for (i,key) in enumerate(goodKeys):
        if mykey == 'total_time':
          if i ==0:
            totalTime = 0
          else:
            totalTime = int(f[key+'/timing'].attrs.get('secondsSinceEpoch')) - int(f[goodKeys[i-1]+'/timing'].attrs.get('secondsSinceEpoch'))
          myDictOfData[mykey].append(totalTime)
        elif mykey == 'number_of_cells':
          if 'out' in key:
            ex = key+'/cells/cell_radii' in f
            if ex:
              arrayOfRadii = np.asarray(f[key+'/cells/cell_radii'])
              myDictOfData[mykey].append(len(np.concatenate(arrayOfRadii, axis=0)))
            else:
              myDictOfData[mykey].append(0)
        else:
          #myDictOfData[mykey].append(f[key+'/timing'].attrs.get(mykey)/1e6)
          timeForKey = f[key+'/timing'].attrs.get(mykey)
          myDictOfData[mykey].append(timeForKey[0][0])
    ''' memory '''
    for mykey in h5_memory_keys:
      for (i,key) in enumerate(goodKeys):
        ex = key +'/memory' in f;
        if ex:
          value = float(f[key+'/memory'].attrs.get(mykey))
          ''' value is in bytes'''
          value = value/1000 #KB
          value = value/1000 #MB
          value = value/1000 #GB
          myDictOfDataForMemory[mykey].append(value)
    exist_system_parameters_in_hdf = '/parameters/system' in f
    if exist_system_parameters_in_hdf:
      no_of_threads = f['/parameters/system'].attrs.get('num_threads')
    else:
      no_of_threads = 42
  
  plotVsCells = False;
  fig, ax = plt.subplots()
  for mykey in myDictOfDataForMemory:
    print("mykey: %s" % mykey)
    if not mykey == 'number_of_cells':
      if plotVsCells:
          ax.plot(myDictOfData['number_of_cells'],myDictOfDataForMemory[mykey][0:],label=mykey)
      else:
          ax.plot(myDictOfDataForMemory[mykey][0:],label=mykey)
        
  
  if plotVsCells:
    ax.set_xlabel('#cells')
    ax.set_ylabel('use Memory/ GB')
  else:
    ax.set_xlabel('output group/ one unit is 1h of simulated time')
    ax.set_ylabel('use Memory/ GB')
    
  ax.set_title('cluster run of faketumor with vbl \n run on snowden with # %i threads' % int(no_of_threads))
  
  legend = ax.legend(loc='center', shadow=True)
  if interactive:
    plt.show()
  else:
    pp.savefig()   
    
if __name__ == '__main__':
  import argparse
  parser = argparse.ArgumentParser(description='Analyze MTS runtime')  
  #parser.add_argument('tumParamSet', help='Valid configuration are found in /py/krebsjobs/parameters/fparameterSetsFakeTumor.py')
  parser.add_argument('--v',dest='vbl_timing_filename', type=str, default='timing_of_runMainLoop_in_run.txt', help='timing file to evalutate')
  parser.add_argument('--s',dest='vbl_simulation_output_filename', type=str, default='safe.h5', help='output file name in hdf5 format')
  
  interactive = False;
  goodArguments, otherArguments = parser.parse_known_args()
  #qsub.parse_args(otherArguments)
  if 0:
    with PdfPages('runtime_analysis.pdf') as pp:
      plot_vbl_runtime(goodArguments,pp)
      plot_no_cells_over_time_since_epoch(goodArguments,pp)
  if 1:
    with PdfPages('runtime_analysis_file_%s.pdf' % str(goodArguments.vbl_simulation_output_filename)) as pp:
      plot_runtime_from_h5(goodArguments,pp)
      plot_memory_from_h5(goodArguments,pp)
