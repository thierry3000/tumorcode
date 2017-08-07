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
  

from os.path import join, dirname
import qsub
from dicttoinfo import dicttoinfo, Vec
from copy import deepcopy
import numpy as np
import myutils
import krebsutils #---> creates multiprocessor environment,also in __init__ tumor
import krebs
import krebs.tumors

import krebsjobs.parameters.parameterSetsBulkTissueTumor as parameterSets
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
#parser = argparse.ArgumentParser(parents=[general_group])
#subparsers=parser.add_subparsers(dest='action')
#subparsers.add_parser('Restart',parents=[general_group, second_group])
#subparsers.add_parser('Start', parents=[general_group])
  import argparse
  parent_parser = argparse.ArgumentParser(add_help=False,description='Compute BulkTissue tumor. Either with or without vessels')
  parent_parser.add_argument('--no_vessel', help = 'compute the continuum model of tumor cells, no vessels needed for that', default=False, action='store_true')
  
  known, other = parent_parser.parse_known_args()

  
  parser = argparse.ArgumentParser(add_help=False) 
  subparsers = parser.add_subparsers()  
  
  if known.no_vessel:
    #s = p2.add_subparsers()
    parser_no_vessels = subparsers.add_parser('no', parents=[parent_parser])
    parser_no_vessels.add_argument('tumParamSet', help='Valid configuration are found in /py/krebsjobs/parameters/parameterSetsBulkTissueTumor.py')
    goodArguments, otherArguments = parser_no_vessels.parse_known_args()
    qsub.parse_args(otherArguments)
    try:
      if not goodArguments.tumParamSet in dir(parameterSets):
        raise AssertionError('Unknown parameter set %s!' % goodArguments.tumParamSet)
    except Exception, e:
      print e.message
      sys.exit(-1)
    factory = getattr(parameterSets, goodArguments.tumParamSet)
    run_no_vessels(goodArguments.tumParamSet, factory, '1GB', 2.)
  #p1 = argparse.ArgumentParser( parents = [ p2 ] )
  #s = p1.add_subparsers()
  #p = s.add_parser( 'group' )
  #p.set_defaults( group=True )  
  if not known.no_vessel:
  #p1 = argparse.ArgumentParser(description='Compute BulkTissue tumor. Either with or without vessels')  
    #s2 = p2.add_subparsers()
    parser_vessels = subparsers.add_parser('with', parents=[parent_parser])    
    parser_vessels.add_argument('tumParamSet', help='Valid configuration are found in /py/krebsjobs/parameters/tumorParameters.py')
  #this is not needed in the case without vessels
    parser_vessels.add_argument('vesselFileNames', nargs='*', type=argparse.FileType('r'), default=sys.stdin, help='Vessel file to calculate')
  #parser.add_argument('--no_vessel', help = 'compute the continuum model of tumor cells, no vessels needed for that', default=False, action='store_true')
    #parser_vessels.add_argument('no_vessel', help = 'compute the continuum model of tumor cells, no vessels needed for that', default=False, action='store_true')
  
    goodArguments, otherArguments = parser_vessels.parse_known_args()
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
    

    for fn in filenames:
      run_with_vessels(fn, goodArguments.tumParamSet, factory, '5GB', 28.)
