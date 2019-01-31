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

if __name__ == '__main__':
  import os.path, sys
  sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'..'))
  
import os, sys
from os.path import join, basename, dirname, splitext, abspath
sys.path.append(abspath(join(dirname(__file__),'..')))
import qsub
import glob
import h5py
import numpy as np
import krebsutils
import krebs
import krebs.tumors
from copy import deepcopy
from dicttoinfo import dicttoinfo, Vec
from myutils import f2s

#exe = krebsutils.get_full_tumor_executable_path('tum-only-vessels')
dstdir = os.getcwd()

#note: according to doeme et al. the tumor radius growth should be around 0.7 mu m / h

import krebsjobs
import krebsjobs.parameters.parameterSetsMTS as parameterSets
from krebsjobs.parameters import parameterSetsO2

#from submitAdaption import create_auto_dicts

from krebsutils import typelist


def getDomainSizeFromVesselFile(fn):
  with h5py.File(fn, 'r') as f:
    ld = krebsutils.read_lattice_data_from_hdf(krebsutils.find_lattice_group_(f['vessels']))
    size = np.amax(ld.GetWorldSize())
  # longest axis times the lattice spacing
  return size


def MakeVesselFilenamePart(fn):
  with h5py.File(fn, mode='a') as f:
    if 'parameters' in f:
      if 'MESSAGE' in f['parameters'].attrs:
        msg = f['parameters'].attrs['MESSAGE']
        ensemble_index = f['parameters'].attrs['ENSEMBLE_INDEX']
        if msg.startswith('vessels-'): msg=msg[len('vessels-'):]
    if 'msg' not in locals():
      msg = "hm"
      ensemble_index = 1
      f['parameters'].attrs['MESSAGE'] = msg 
      f['parameters'].attrs['ENSEMBLE_INDEX'] = ensemble_index
    
  name = '%s-sample%02i' % (msg, ensemble_index)
  return name

''' returns estimated runtime in hours'''
def estimateMTS_runtime(runTimeIndex):
#    with h5py.File(fn_of_previous_run, 'r') as previousFile:
#      last_out_num = int(previousFile['last_state'].attrs['OUTPUT_NUM'][0])
#      print(last_out_num)
#    if last_out_num < 100:
#      return 1
#    if last_out_num >=100 and last_out_num < 200:
#      return 2
  if runTimeIndex < 100:
    return 2./60.
  if runTimeIndex < 200:
    return 5./60. # 5min
  if runTimeIndex < 300:
    return 0.5
  if runTimeIndex < 400:
    return 1
  if runTimeIndex < 500:
    return 6
  if runTimeIndex < 600:
    return 12
  print('estimate hours for %i: %d' % (runTimeIndex, hours))
  #return hours
''' returns estimated memory in MB'''
def estimateMTS_memory(runTimeIndex):

  if runTimeIndex < 100:
    return '1GB'
  if runTimeIndex < 200:
    return '2GB'
  if runTimeIndex < 300:
    return '4GB'
  if runTimeIndex < 400:
    return '60GB'
  if runTimeIndex < 500:
    return '90GB'
  if runTimeIndex < 600:
    return '110GB'

def estimateMTS_nProc(runTimeIndex):
  if runTimeIndex < 2:
    return 16
  if runTimeIndex < 100:
    return 16
  if runTimeIndex < 200:
    return 16
  if runTimeIndex < 300:
    return 16
  if runTimeIndex < 400:
    return 20
  if runTimeIndex < 600:
    return 28
def checkIfFileExists(fn):
  print('checking if file exists...')
  if (os.path.isfile(fn)):
    print('your are about to destroy the file: %s' % fn)
    sys.exit()
    
def run_pipeline(vessel_fn, name, paramSet, mem, days, pipelineLength):
  ''' initial job'''
  name, paramSet = krebsjobs.PrepareConfigurationForSubmission(vessel_fn, name, 'fakeTumMTS', paramSet)
  print("name: %s" %name)
  configstr = dicttoinfo(paramSet)
  config_file_name = '%s.info' % paramSet['fn_out']
  
  checkIfFileExists(paramSet['fn_out'])
  
  with open(config_file_name, 'w') as f:
    f.write(configstr)
  jobID_to_continue=qsub.submit(qsub.func(krebs.tumors.run_faketum_mts, config_file_name),
                            name = 'job_'+name,
                            mem = mem,
                            hours = 1,
                            change_cwd = True)
  '''now we store in the same folder --> better for cluster'''
#  if qsub.determine_submission_program_() == 'run_locally':
#    fn_of_previous_run = name+'.h5'
#  else:
#    fn_of_previous_run = '%i/%s.h5' % (jobID_to_continue, name)
  
  for i in range(pipelineLength):
    jobID_to_continue = qsub.submit(qsub.func(krebs.tumors.rerun_faketum_mts, name+'.h5'),
                            name = 'job_'+ name,
                            num_cpus = estimateMTS_nProc(i),
                            mem = estimateMTS_memory(i),
                            hours = estimateMTS_runtime(i),
                            change_cwd = True,
                            dependOn = jobID_to_continue)
  
  
def run(vessel_fn, name, paramSet, mem, days):
  
  name, paramSet = krebsjobs.PrepareConfigurationForSubmission(vessel_fn, name, 'fakeTumMTS', paramSet)
  configstr = dicttoinfo(paramSet)
  config_file_name = '%s.info' % paramSet['fn_out']

  checkIfFileExists(paramSet['fn_out'])  
  
  with open(config_file_name, 'w') as f:
    f.write(configstr)
    
  #o2params = getattr(parameterSetsO2, "breastv3")
  qsub.submit(qsub.func(krebs.tumors.run_faketum_mts, config_file_name),
                            name = 'job_'+name,
                            mem = mem,
                            days = days,
                            #num_cpus = paramSet['num_threads'],
			    #set by slurm
                            change_cwd = True)

def rerun(fn_of_previous_run, job_name, mem, days):
  #name, c = krebsjobs.PrepareConfigurationForSubmission(vessel_fn, name, 'fakeTum', paramSet)
  #configstr = dicttoinfo(c)
  #config_file_name = '%s.info' % c['fn_out']
  #with open(config_file_name, 'w') as f:
  #  f.write(configstr)
  
  estimateMTS_runtime(fn_of_previous_run)
  qsub.submit(qsub.func(krebs.tumors.rerun_faketum_mts, fn_of_previous_run),
                            name = 'job_'+ job_name,
                            mem = mem,
                            hours = estimateMTS_runtime(fn_of_previous_run),
			    #days = days,
                            #num_cpus = 2,
			    #set by slurm
                            change_cwd = True)

if __name__ == '__main__':
  import argparse
  parser = argparse.ArgumentParser(description='Compute Fake tumor. Vessels are needed for that')  
  parser.add_argument('tumParamSet', help='Valid configuration are found in /py/krebsjobs/parameters/fparameterSetsFakeTumor.py')
  parser.add_argument('vesselFileNames', nargs='*', type=argparse.FileType('r'), default=sys.stdin, help='Vessel file to calculate')
  parser.add_argument('--rerun', help='allows to rerun a simulation', default=False, action='store_true')
  parser.add_argument('--pipelineLength', help= 'number of successsive pipeline calls', type=int ,default = 0)
  
  goodArgumentsMTS, otherArgumentsMTS = parser.parse_known_args()
  qsub.parse_args(otherArgumentsMTS)
  
  if(not goodArgumentsMTS.rerun):
    tumorParameterName = goodArgumentsMTS.tumParamSet
    #create filename due to former standards
    filenames=[]
    for fn in goodArgumentsMTS.vesselFileNames:
      filenames.append(fn.name)
      
    try:
      if not (tumorParameterName in dir(parameterSets)) and (not 'auto' in tumorParameterName):
          raise AssertionError('Unknown parameter set %s!' % tumorParameterName)
      for fn in filenames:
          if not os.path.isfile(fn):
              raise AssertionError('The file %s is not present!'%fn)
    except Exception, e:
      print e.message
      sys.exit(-1)
  
    
    factory = getattr(parameterSets, tumorParameterName)
    for fn in filenames:
      #run(fn, tumorParameterName, factory, '4GB', 2.)
      #run(fn, tumorParameterName, factory, '2GB', 5.)
      run_pipeline(fn, tumorParameterName, factory, '2GB', 5., goodArgumentsMTS.pipelineLength)
  else:
    print('starting rerun with file: %s' % goodArgumentsMTS.vesselFileNames[0].name)
    if not os.path.isfile(goodArgumentsMTS.vesselFileNames[0].name):
      raise AssertionError('The file %s is not present!'%goodArgumentsMTS.vesselFileNames[0].name)
    string_to_provide = str(goodArgumentsMTS.vesselFileNames[0].name)
    goodArgumentsMTS.vesselFileNames[0].close()
    rerun(string_to_provide, 'rerun_of_', '2GB', 5.)
  
        
