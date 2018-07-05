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

import krebsjobs

#exe = krebsutils.get_full_tumor_executable_path('tum-only-vessels')
dstdir = os.getcwd()

#note: according to doeme et al. the tumor radius growth should be around 0.7 mu m / h

import krebsjobs.parameters.parameterSetsFakeTumor as parameterSets
#from submitAdaption import create_auto_dicts

from krebsutils import typelist


def getDomainSizeFromVesselFile(fn):
  with h5files.open(fn, 'r') as f:
    ld = krebsutils.read_lattice_data_from_hdf(krebsutils.find_lattice_group_(f['vessels']))
    size = np.amax(ld.GetWorldSize())
  # longest axis times the lattice spacing
  return size


#def MakeVesselFilenamePart(fn):
#  with h5files.open(fn, mode='a') as f:
#    if 'parameters' in f:
#      if 'MESSAGE' in f['parameters'].attrs:
#        msg = f['parameters'].attrs['MESSAGE']
#        ensemble_index = f['parameters'].attrs['ENSEMBLE_INDEX']
#        if msg.startswith('vessels-'): msg=msg[len('vessels-'):]
#    if 'msg' not in locals():
#      msg = "hm"
#      ensemble_index = 1
#      f['parameters'].attrs['MESSAGE'] = msg 
#      f['parameters'].attrs['ENSEMBLE_INDEX'] = ensemble_index
#    
#  name = '%s-sample%02i' % (msg, ensemble_index)
#  return name


#def PrepareConfigurationForSubmission(vessel_fn, name, prepend, config_):
#  if vessel_fn is not None: 
#    dstdir = os.getcwd()  
#    c = deepcopy(config_)
#    vessel_fn_part = MakeVesselFilenamePart(vessel_fn)
#    out_fn = join(dstdir, '%s-%s-%s.h5' % (prepend,vessel_fn_part, name))
#    print 'generating tumor run with'
#    print '  vessel:', vessel_fn
#    print '  output:', out_fn
#    print ' paramset:', name
#    c['fn_out'] = out_fn
#    c['fn_vessel'] = vessel_fn
#    c['paramset_name'] = name
#    name = splitext(basename(out_fn))[0]
#  else:
#    dstdir = os.getcwd()  
#    c = deepcopy(config_)
#    #vessel_fn_part = MakeVesselFilenamePart(vessel_fn)
#    out_fn = join(dstdir, '%s-%s.h5' % (prepend,name))
#    print 'generating tumor run no vessels'
#    print '  output:', out_fn
#    print ' paramset:', name
#    c['fn_out'] = out_fn
#    c['paramset_name'] = name
#    name = splitext(basename(out_fn))[0]
#  return name, c
  
#def worker_on_client(fn, grp_pattern, tumParams, num_threads):
#  print('Fake tum on %s / %s / param: %s' % (fn, grp_pattern, tumParams['name']))
#  h5files.search_paths = [dirname(fn)] # so the plotting and measurement scripts can find the original tumor files using the stored basename alone
#  krebsutils.set_num_threads(num_threads)
#  
#  
#  fake_tum_refs = krebs.tumors.run_faketum(fn, grp_pattern, tumParams)
#  
#  h5files.closeall() # just to be sure

def run(vessel_fn, name, paramSet, mem, days):
  name, c = krebsjobs.PrepareConfigurationForSubmission(vessel_fn, name, 'fakeTum', paramSet)
  configstr = dicttoinfo(c)
  config_file_name = '%s.info' % c['fn_out']
  with open(config_file_name, 'w') as f:
    f.write(configstr)
  qsub.submit(qsub.func(krebs.tumors.run_faketum, config_file_name),
                            name = 'job_'+name,
                            mem = mem,
                            days = days,
                            num_cpus = c['num_threads'],
                            change_cwd = True)

def rerun(fn_of_previous_run, job_name, mem, days):
  #name, c = krebsjobs.PrepareConfigurationForSubmission(vessel_fn, name, 'fakeTum', paramSet)
  #configstr = dicttoinfo(c)
  #config_file_name = '%s.info' % c['fn_out']
  #with open(config_file_name, 'w') as f:
  #  f.write(configstr)
  qsub.submit(qsub.func(krebs.tumors.rerun_faketum, fn_of_previous_run),
                            name = 'job_'+ job_name,
                            mem = mem,
                            days = days,
                            #num_cpus = c['num_threads'],
                            change_cwd = True)

if __name__ == '__main__':
  import argparse
  parser = argparse.ArgumentParser(description='Compute Fake tumor. Vessels are needed for that')  
  parser.add_argument('tumParamSet', help='Valid configuration are found in /py/krebsjobs/parameters/fparameterSetsFakeTumor.py')
  parser.add_argument('vesselFileNames', nargs='*', type=argparse.FileType('r'), default=sys.stdin, help='Vessel file to calculate')
  parser.add_argument("--rerun", help=" ", default=False, action="store_true")
  
  goodArguments, otherArguments = parser.parse_known_args()
  qsub.parse_args(otherArguments)
  
  if not goodArguments.rerun:
      tumorParameterName = goodArguments.tumParamSet
      #create filename due to former standards
      filenames=[]
      for fn in goodArguments.vesselFileNames:
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
    
    #  if not 'auto' in tumorParameterName:
      factory = getattr(parameterSets, tumorParameterName)
      if type(factory).__name__ == 'function':
        configs = factory(len(filenames))
        for fn, cfg in zip(filenames, configs):
          run(fn, factory.name, cfg, '4GB', 2.)
      else:
        for fn in filenames:
          #run(fn, tumorParameterName, factory, '4GB', 2.)
          run(fn, tumorParameterName, factory, '2GB', 5.)
  else:
    print('starting rerun with file: %s' % goodArguments.vesselFileNames[0].name)
    if not os.path.isfile(goodArguments.vesselFileNames[0].name):
      raise AssertionError('The file %s is not present!'%goodArguments.vesselFileNames[0].name)
    string_to_provide = str(goodArguments.vesselFileNames[0].name)
    goodArguments.vesselFileNames[0].close()
    rerun(string_to_provide, 'rerun_of_', '2GB', 5.)
    