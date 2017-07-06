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
import qsub
import dicttoinfo
import krebsutils
import myutils
import h5py
import h5files
import itertools
import copy

from krebs import detailedo2
from krebsjobs.parameters import parameterSetsO2
#from krebs import povrayRenderOxygenDetailed
from krebs import detailedo2Analysis


def worker_on_client(fn, pattern, o2params):
  print 'detailedo2 on %s / %s' % (fn, pattern)
  h5files.search_paths = [dirname(fn)] # so the plotting and measurement scripts can find the original tumor files using the stored basename alone
  num_threads = o2params.pop('num_threads')
  krebsutils.set_num_threads(num_threads)
  o2_refs = detailedo2.doit(fn, pattern, (o2params, o2params['name']))
  if 0: #this is for data analysis on the clusters
    for ref in o2_refs:
      po2group = h5files.open(ref.fn)[ref.path]
      detailedo2Analysis.WriteSamplesToDisk(po2group)
  h5files.closeall() # just to be sure


def worker_plots_for_paper(filenames, pattern):
  from krebs.detailedo2Analysis import plotsForPaper
  import h5files
  h5files.search_paths.append('..')
  plotsForPaper.doit(filenames, pattern)

def worker_render_images(fn, pattern, num_threads):
  from krebs import povrayRenderOxygenDetailed
  import h5files
  h5files.search_paths.append('..')
  rendersettings = dict(
    format = 'png',
    res=(1024, 1024),
    aa=2,
    num_threads=num_threads,
  )
  povrayRenderOxygenDetailed.doit(fn, pattern, rendersettings)


def prepareParametersWithNumThreads(o2params, numThreadsOverride):
  ''' insert num_threads in o2params dict if not already present or
      override is used. Returns tuple of provided o2params and num_threads used.'''
  if numThreadsOverride:
    o2params['num_threads'] = num_threads = numThreadsOverride
  else:
    try:
      num_threads = o2params['num_threads']
    except KeyError:
      o2params['num_threads'] = num_threads = 6
  return o2params, num_threads


def run(parameter_set_name, filenames, grp_pattern, systemsize):
  print 'submitting ...', parameter_set_name
  print dicttoinfo.dicttoinfo(getattr(parameterSetsO2, parameter_set_name))
  print 'for files', filenames

  dirs = set()
  for fn in filenames:
    with h5py.File(fn, 'r') as f:
      d = myutils.walkh5(f, grp_pattern)
      assert len(d), 'you fucked up, pattern "%s" not found in "%s"!' % (grp_pattern, fn)
      dirs =set.union(dirs, d)
  print 'and resolved groups therein: %s' % ','.join(dirs)

  o2params = getattr(parameterSetsO2, parameter_set_name)
  o2params['name'] = parameter_set_name
  if callable(o2params):
    o2paramsList = o2params(len(filenames))
  else:
    o2paramsList = itertools.repeat(o2params)

  for (o2params, fn) in zip(o2paramsList, filenames):
    o2params, num_threads = prepareParametersWithNumThreads(copy.deepcopy(o2params), systemsize)
    qsub.submit(qsub.func(worker_on_client, fn, grp_pattern, o2params),
                  name = 'job_o2_'+parameter_set_name+'_'+basename(fn),
                  num_cpus = num_threads,
                  days = 5,  # about one day per thread, the idea being that number of threads is proportional to systems size and runtime is 
                  mem = '%iMB' % (2000*num_threads),
                  change_cwd = True)


if not qsub.is_client and __name__=='__main__':
  import argparse
  parser = argparse.ArgumentParser(description='Compute po2 distributions and analyze them. In normal mode it takes arguments: parameter set name, filenames, h5 group pattern. In analysis mode -a, it takes arguments: filenames, h5 group pattern.')  
  parser.add_argument('o2params', help = 'choose the parameter for the simulation. found at /py/krebsjobs/parameters/parameterSetsO2.py')  
  parser.add_argument('vesselFileNames', nargs='*', type=argparse.FileType('r'), default=sys.stdin, help='Vessel file to calculate')   
  parser.add_argument('grp_pattern',help='Where to find the vessel group in the file')  
  parser.add_argument('-a', '--analyze', help = 'loop through all files analyze data and make plot', default=False, action='store_true')
  parser.add_argument('-r', '--render',  help = 'launch rendering of overview images', default = False, action = 'store_true')
  parser.add_argument('--systemsize'   , help = 'num threads = 1 * systemsize, memory = 2 Gb * systemsize (for o2 computation; analysis uses different scaling (see source)', default = None, action = 'store')
  goodArguments, otherArguments = parser.parse_known_args()
  #print('before: %i' % qsub.goodArgumentsQueue.days)
  qsub.parse_args(otherArguments)
  #print('after: %i' % qsub.defaultDays)
  if qsub.goodArgumentsQueue.days:
    print('qsub, days used: set to %i' % qsub.goodArgumentsQueue.days)

  #create filename due to former standards
  filenames=[]
  for fn in goodArguments.vesselFileNames:
    filenames.append(fn.name)
  
  if goodArguments.analyze:
    filenames   = goodArguments.vesselFileNames
    grp_pattern = goodArguments.grp_pattern
    try:
      systemsize = int(goodArguments.systemsize)
    except:
      print 'no valid --systemsize given, using 2 as default'
      systemsize = 2
    qsub.submit(qsub.func(worker_plots_for_paper, filenames, grp_pattern),
                  name = 'job_o2_analysis',
                  num_cpus = 1, # single threaded only, but scale memory usage with num_threads
                  days = 2,
                  mem = '%iMB' % (1000*systemsize),
                  change_cwd = True)
  if goodArguments.render:
    #filenames   = parseResult.args[:-1]
    grp_pattern = goodArguments.grp_pattern    
    try:
      systemsize = int(goodArguments.systemsize)
    except:
      print 'no valid --systemsize given, using 2 as default'
      systemsize = 2
    for fn in filenames:
      qsub.submit(qsub.func(worker_render_images, fn, grp_pattern, systemsize),
                  name = 'job_o2_render',
                  num_cpus = systemsize,
                  days = 1,
                  mem = '%iMB' % (1000*systemsize),
                  change_cwd = True)
  if not goodArguments.analyze and not goodArguments.render:
    try:
      systemsize = goodArguments.systemsize
      if not goodArguments.o2params in dir(parameterSetsO2):
        raise AssertionError('Unknown parameter set %s!' % goodArguments.o2params)
      for fn in filenames:
          if not os.path.isfile(fn):
              raise AssertionError('The file %s is not present!'%fn)
    except Exception, e:
      print e.message
      sys.exit(-1)
    run(goodArguments.o2params, filenames, goodArguments.grp_pattern, systemsize)
