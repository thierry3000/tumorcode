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
import cProfile

import krebs
import krebs.adaption

import krebsjobs.parameters
from krebsjobs.parameters import parameterSetsAdaption

from krebsutils import typelist 

def worker_on_client(fn, grp_pattern, adaptionParams, num_threads):
  print('Adaption on %s / %s / param: %s' % (fn, grp_pattern, adaptionParams['name']))
  qsub.printClientInfo()  
  h5files.search_paths = [dirname(fn)] # so the plotting and measurement scripts can find the original tumor files using the stored basename alone
  krebsutils.set_num_threads(num_threads)
  
  #params['name'] = parameter_set_name
  adaptionParams['adaption'].update(
      vesselFileName = fn,
      vesselGroupName = grp_pattern,
      )

  krebs.adaption.doit( adaptionParams)
  
  #h5files.closeall() # just to be sure

''' not we cannot use the whole goodArguments namespace here,
    because these object need to be pickeled on the 
    cluster!
'''
def worker_on_client_optimize(fn, grp_pattern, adaptionParams, num_threads, timing):
  print('Adaption on %s / %s / param: %s' % (fn, grp_pattern, adaptionParams['name']))
  h5files.search_paths = [dirname(fn)] # so the plotting and measurement scripts can find the original tumor files using the stored basename alone
  krebsutils.set_num_threads(num_threads)
  
  if timing:
    prof = cProfile.Profile()
    
    adaption_refs_optimize = prof.runcall(krebs.adaption.doit_optimize, fn, adaptionParams['adaption'],adaptionParams['calcflow'])
    #prof.dump_stats(file)
    prof.print_stats()

  else:
    adaption_refs_optimize = krebs.adaption.doit_optimize(fn, adaptionParams['adaption'],adaptionParams['calcflow'])
  print("got back: ")
  print(adaption_refs_optimize)
  h5files.closeall() # just to be sure

''' I tried to assing parameters due to the RC type here 
    which was kind of stupid
'''

    
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

def run_optimize(parameter_set, filenames, grp_pattern, timing):
  print 'submitting ...', parameter_set['name']
 

  num_threads = 1
  if 'num_threads' in parameter_set:
    num_threads = parameter_set['num_threads']
    
  for fn in filenames:
    qsub.submit(qsub.func(worker_on_client_optimize, fn, grp_pattern, parameter_set, num_threads, timing),
                  name = 'job_adaption_'+parameter_set['name']+'_'+basename(fn),
                  num_cpus = num_threads,
                  days = 4.,
                  mem = '3500MB',
                  change_cwd = True)

def do_optimize_adaption_deap():
  print("need pso here")


if not qsub.is_client and __name__=='__main__':
  import argparse
  parser = argparse.ArgumentParser(description='Compute adaption see Secomb model', formatter_class=argparse.ArgumentDefaultsHelpFormatter)  
  parser.add_argument('AdaptionParamSet', default="deap_test")
  #this is not needed in the case without vessels
  parser.add_argument('vesselFileNames', nargs='*', type=argparse.FileType('r'), default=sys.stdin)
  parser.add_argument('grp_pattern', default="adaption/vessels_after_adaption")
  parser.add_argument('-t','--tumorParams', help='by explicitly enable this you can use tumor parameters for the adaption as well', action='store_true')
  parser.add_argument('--time', default=False, help='Show time profile', action='store_true')
  parser.add_argument('--ks', help='ks')
  parser.add_argument('--kc', help='kc')
  parser.add_argument('--km', help='km')
  
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
          raise AssertionError('pattern "%s" not found in "%s"!' % (goodArguments.grp_pattern, fn))
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
  if factory.__class__ == list:
    factory=factory[7]
    print("warning: you are using several parameter sets")
  #single parameter set chosen  
  if factory.__class__ == dict:
    factory['name'] = goodArguments.AdaptionParamSet
    #run_optimize(factory, filenames, goodArguments.grp_pattern, goodArguments.time)
    if goodArguments.kc is not None:
      print("setting kc = %f" % float(goodArguments.kc))
      factory['adaption']['k_c'] = float(goodArguments.kc)
    if goodArguments.ks is not None:
      print("setting ks = %f" % float(goodArguments.ks))
      factory['adaption']['k_s'] = float(goodArguments.ks)
    if goodArguments.km is not None:
      print("setting km = %f" % float(goodArguments.km))
      factory['adaption']['k_m'] = float(goodArguments.km)
  else:
    print("WARNING: factory is not a dict")

  returnState, mean, varOfMean, total_surface = krebs.adaption.do_simple_adaption(filenames[0],goodArguments.grp_pattern, factory)
