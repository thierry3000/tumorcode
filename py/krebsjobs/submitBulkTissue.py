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
  

from os.path import join, dirname
import qsub
from dicttoinfo import dicttoinfo, Vec
from copy import deepcopy
import numpy as np
import myutils
#import krebsutils #---> creates multiprocessor environment,also in __init__ tumor
import krebs
import krebs.tumors

import krebsjobs.parameters.tumorParameters as parameterSets
from krebsjobs.submitFakeTum import PrepareConfigurationForSubmission


def f2s(x):
  return myutils.f2s(x, exponential=False, latex=False)

def mkdir(outfn):
  p = dirname(outfn)
  if p and not os.path.isdir(p):
    print 'creating', p
    os.makedirs(p)



def run_with_vessels(vessel_fn, name, config_, mem, days):
  name, config_ = PrepareConfigurationForSubmission(vessel_fn, name, 'tumBulk', config_)
  config_ = krebs.tumors.set_lattice_size(config_, vessel_fn)
  
  sx, sy, sz = config_['lattice_size']
  print('cont size: %i, %i, %i' % (sx, sy, sz))
  #c = PrepConfig_new_python_dict(c)  
  configstr = dicttoinfo(config_)
  qsub.submit(qsub.func(krebs.tumors.run_bulktissue_w_vessels,configstr),
                        name = 'job_'+name,
                        mem = mem,
                        days = days,
                        num_cpus = config_['num_threads'],
                        change_cwd = True)
def run_no_vessels(name, config_, mem, days):
  name, config_ = PrepareConfigurationForSubmission(None, name, 'tumBulk', config_)
  #config_ = set_lattice_size(config_, None)
  if 1: #HACK
    print('Warning: lattice_size not set properly')
    config_['lattice_size'] = Vec((20,20,20))
  try:
    if not 'lattice_size' in config_.keys():
      raise AssertionError('No lattice_size found in configuration %s' % getattr(config_))
    if type(config_['lattice_size'])==str:
      if config_['lattice_size']=='set me to match the vessel domain':
        raise AssertionError('Find better lattice size' )
  except Exception, e:
    print e.message
    sys.exit(-1) 
  sx, sy, sz = config_['lattice_size']
  print('cont size: %i, %i, %i' % (sx, sy, sz))
  #c = PrepConfig_new_python_dict(c)  
  configstr = dicttoinfo(config_)
  qsub.submit(qsub.func(krebs.tumors.run_bulktissue_no_vessels,configstr),
                        name = 'job_'+name,
                        mem = mem,
                        days = days,
                        num_cpus = config_['num_threads'],
                        change_cwd = True)

if not qsub.is_client and __name__ == '__main__':

  import argparse
  parser = argparse.ArgumentParser(description='Compute BulkTissue tumor. Either with or without vessels')  
  parser.add_argument('tumParamSet', help='Valid configuration are found in /py/krebsjobs/parameters/tumorParameters.py')
  #this is not needed in the case without vessels
  parser.add_argument('vesselFileNames', nargs='*', type=argparse.FileType('r'), default=sys.stdin, help='Vessel file to calculate')
  parser.add_argument('--no_vessel', help = 'compute the continuum model of tumor cells, no vessels needed for that', default=False, action='store_true')

  goodArguments, otherArguments = parser.parse_known_args()
  qsub.parse_args(otherArguments)
  
  try:
    if not goodArguments.tumParamSet in dir(parameterSets):
      raise AssertionError('Unknown parameter set %s!' % goodArguments.tumParamSet)
  except Exception, e:
    print e.message
    sys.exit(-1)

  factory = getattr(parameterSets, goodArguments.tumParamSet)
  #create filename due to former standards
  filenames=[]
  for fn in goodArguments.vesselFileNames:
    filenames.append(fn.name)
  
  if goodArguments.no_vessel:
    run_no_vessels(goodArguments.tumParamSet, factory, '1GB', 2.)
  else:
    for fn in filenames:
      run_with_vessels(fn, goodArguments.tumParamSet, factory, '5GB', 28.)
