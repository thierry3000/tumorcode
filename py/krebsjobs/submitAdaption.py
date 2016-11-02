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
import identifycluster
if identifycluster.name == 'snowden':
  import matplotlib
  matplotlib.use('Agg') 
import qsub
import dicttoinfo
import krebsutils
import myutils
import h5py
import h5files
import itertools
import time
import string

import krebs.adaption

from krebsjobs.parameters import parameterSetsAdaption

#globals
typelist = 'typeA typeB typeC typeD typeE typeF typeG typeH typeI'.split()

def create_auto_dicts(param_group):
  cadidatesList = []
  type_to_highest_number = dict()
  for t in typelist:
    type_to_highest_number[t] = 0
  for adict in dir(parameterSets_tum):
    if param_group in adict:
     for t in typelist:
       if t in adict:
         latest_vary_of_type =0
         index_start_type = string.find(adict,t+'_vary')
         #print adict[index_start_type+10:]
         zahl = int(adict[index_start_type+10:])
         if 0:#this gives no value for non present configs
           if( zahl > type_to_highest_number[t]):
             type_to_highest_number[t] = zahl
         if 1:
           type_to_highest_number[t] = zahl
  type_to_parameterset = dict()
  for t in typelist:
    astring = param_group+t+'_vary'+str(type_to_highest_number[t])
    #check if exists
    if astring in dir(parameterSets_tum):
      type_to_parameterset[t] = astring
    else:
      print('Warning no paramset found for: %s' % astring)
  return type_to_parameterset    
         
         
       

def worker_on_client(fn, grp_pattern, adaptionParams, num_threads):
  print('Adaption on %s / %s / param: %s' % (fn, grp_pattern, adaptionParams['name']))
  h5files.search_paths = [dirname(fn)] # so the plotting and measurement scripts can find the original tumor files using the stored basename alone
  krebsutils.set_num_threads(num_threads)
  
  #params['name'] = parameter_set_name
  adaption_refs = krebs.adaption.doit(fn, grp_pattern, adaptionParams)
  
#  for ref in o2_refs:
#    po2group = h5files.open(ref.fn)[ref.path]
#    detailedo2Analysis.WriteSamplesToDisk(po2group)
#  povrayRenderOxygenDetailed.doit(o2_refs[-1].fn, o2_refs[-1].path)
  h5files.closeall() # just to be sure

''' I tried to assing parameters due to the RC type here '''
def run(parameter_set_name, filenames, grp_pattern):
  print 'submitting ...', parameter_set_name
  usetumorparams = True
  if (not 'auto' in parameter_set_name):
    if usetumorparams:
      print dicttoinfo.dicttoinfo(getattr(parameterSets_tum, parameter_set_name))
      print 'for files', filenames
    else:
      print dicttoinfo.dicttoinfo(getattr(parameterSetsAdaption, parameter_set_name))
      print 'for files', filenames
 
  

  dirs = set()
  for fn in filenames:
    with h5py.File(fn, 'r') as f:
      d = myutils.walkh5(f, grp_pattern)
      assert len(d), 'you fucked up, pattern "%s" not found in "%s"!' % (grp_pattern, fn)
      dirs =set.union(dirs, d)
  print 'and resolved groups therein: %s' % ','.join(dirs)

  num_threads = 1
  if (not 'auto' in parameter_set_name):
    if usetumorparams:  
      adaptionParams = getattr(parameterSets_tum, parameter_set_name)
      adaptionParams['name'] = parameter_set_name
      if callable(adaptionParams):
        adaptionParamsList = adaptionParams['adaption'](len(filenames))
      else:
        if 'num_threads' in adaptionParams:
          if adaptionParams['num_threads']>num_threads:
            num_threads = adaptionParams['num_threads']
        adaptionParamsList = itertools.repeat(adaptionParams)
    else:
      adaptionParams = getattr(parameterSetsAdaption, parameter_set_name)
      adaptionParams['name'] = parameter_set_name
      if callable(adaptionParams):
        adaptionParamsList = adaptionParams(len(filenames))
      else:
        adaptionParamsList = itertools.repeat(adaptionParams)
  else:
    adaptionParamsList = []
    ### change here for different types
    parameter_set_name_from_pipe = parameter_set_name #begins with auto_
    parameter_set_name =  parameter_set_name_from_pipe[5:]
    print('Found param identifiyer: %s' % parameter_set_name)
    type_to_paramset = create_auto_dicts(parameter_set_name+'_')
    for fn in filenames:
      for t in typelist:
        if t in fn:
          adaptionParams = getattr(parameterSets_tum, type_to_paramset[t]) 
          adaptionParams['name'] = type_to_paramset[t]
          if 'num_threads' in adaptionParams:
            if adaptionParams['num_threads']>num_threads:
              num_threads = adaptionParams['num_threads']
          adaptionParamsList.append(adaptionParams)
#print(type_to_paramset['typeB'])
#        
#        if t in fn:
#          adaptionParams = getattr(parameterSets_tum, 'p2d_11layer_'+t)
#          adaptionParams['name'] = 'p2d_11layer_'+t
#          if 'num_threads' in adaptionParams:
#            if adaptionParams['num_threads']>num_threads:
#              num_threads = adaptionParams['num_threads']
#          adaptionParamsList.append(adaptionParams)

    
  for (adaptionParams, fn) in zip(adaptionParamsList, filenames):
    qsub.submit(qsub.func(worker_on_client, fn, grp_pattern, adaptionParams, num_threads),
                  name = 'job_adaption_'+parameter_set_name+'_'+basename(fn),
                  num_cpus = num_threads,
                  days = 4.,
                  mem = '3500MB',
                  change_cwd = True)
    
def run2(parameter_set, filenames, grp_pattern):
  print 'submitting ...', parameter_set['name']
 

  num_threads = 1
  if 'num_threads' in parameter_set:
    num_threads = parameter_set['num_threads']
    
  for fn in filenames:
    qsub.submit(qsub.func(worker_on_client, fn, grp_pattern, parameter_set, num_threads),
                  name = 'job_adaption_'+parameter_set['name']+'_'+basename(fn),
                  num_cpus = num_threads,
                  days = 4.,
                  mem = '3500MB',
                  change_cwd = True)


if not qsub.is_client and __name__=='__main__':
  import argparse
  parser = argparse.ArgumentParser(description='Compute adaption see Secomb model')  
  parser.add_argument('AdaptionParamSet')
  #this is not needed in the case without vessels
  parser.add_argument('vesselFileNames', nargs='*', type=argparse.FileType('r'), default=sys.stdin)
  parser.add_argument('grp_pattern')
  parser.add_argument('-t','--tumorParams', help='by explicitly enable this you can use tumor parameters for the adaption as well', action='store_true')

  goodArguments, otherArguments = parser.parse_known_args()
  qsub.parse_args(otherArguments)
  
  try:
    if not goodArguments.AdaptionParamSet in dir(parameterSetsAdaption):
      raise AssertionError('Unknown parameter set %s!' % goodArguments.AdaptionParamSet)
    dirs = set()
    for fn in goodArguments.vesselFileNames:
      if not os.path.isfile(fn.name):
        raise AssertionError('The file %s is not present!'%fn)
      with h5py.File(fn.name, 'r') as f:
        d = myutils.walkh5(f, goodArguments.grp_pattern)
        if not len(d)>0:
          raise AssertionError('pattern "%s" not found in "%s"!' % (grp_pattern, fn))
        else:
          dirs = set.union(dirs,d)
  except Exception, e:
    print e.message
    sys.exit(-1)
  
  print('Resolved groups: %s' % ','.join(dirs))
  
  #create filename due to former standards
  filenames=[]
  for fn in goodArguments.vesselFileNames:
    filenames.append(fn.name)
  
  factory = getattr(parameterSetsAdaption, goodArguments.AdaptionParamSet)
  factory['name'] = goodArguments.AdaptionParamSet
  run2(factory, filenames, goodArguments.grp_pattern)
  
      