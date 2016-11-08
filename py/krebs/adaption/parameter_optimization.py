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
#! /usr/bin/env python2
# -*- coding: utf-8 -*-
if __name__ == '__main__':
  import os.path, sys
  sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../..'))
import os, sys
from os.path import  basename, dirname
  
import identifycluster
if identifycluster.name == 'snowden':
  import matplotlib
  matplotlib.use('Agg') 
import qsub
import dicttoinfo
import myutils
import h5py
import h5files
import itertools
import time
import string
import numpy as np

#sys.path.append(join(dirname(__file__),'/localdisk/thierry/tc_install/py/krebs/adaptionAnalysis'))
from krebs.adaption import getVesselTypes
from pso import *
from os.path import basename


import krebs.adaption

from krebsjobs.parameters import parameterSetsAdaption
#import krebsjobs.parameters.tumorParameters as parameterSets_tum

import krebsutils
if sys.flags.debug: # library search path made known by initialization of krebsutils
  adaption_cpp = __import__('libadaption_d', globals(), locals())
else:
  adaption_cpp = __import__('libadaption_', globals(), locals())

#globals
from krebsutils import typelist
calculated_mean_cap_flow = 0.0
use_initial_guess = False

def con(x, *args):
  #global calculated_mean_cap_flow
  return np.square(calculated_mean_cap_flow-10e4)
def con2(x, *args):
  #global calculated_mean_cap_flow
  return np.square(calculated_mean_cap_flow-10e5)
def con3(x, *args):
  #global calculated_mean_cap_flow
  return np.square(calculated_mean_cap_flow-10e6)

def do_pSO_for_file(x,dst_file, vesselgroup, parameters):
  parameters['counter'] = parameters['counter'] + 1
  alabel = 'data_%i' % parameters['counter']
  
  #hack
  parameters['adaption']['avgRootNodeConductivity'] = 0.0
  
  #pso stuff
  #write current parametr variation to parametres object
  #so this goes to h5 file and goes to the adaption
  parameters['adaption']['k_m'] = x[0]
  parameters['adaption']['k_c'] = x[1]
  parameters['adaption']['k_s'] = x[2]
  parameters['adaption']['cond_length'] = x[3]
  #parameters['adaption']['Q_refdot'] = x[4]
  # running adaption
  print('creating dataset: %s' % alabel)
  gdst = dst_file.create_group(alabel)
  
  if parameters['counter']%50>0:
    parameters['adaption']['write2File'] = False
  else:
    parameters['adaption']['write2File'] = True
    myutils.hdf_write_dict_hierarchy(gdst, 'parameters', parameters)
  #print('We run: vesselgroup: %s, parameters[\'adaption\']: %s, parameters[\'calcflow\'] %s, gdst %s\n' % (vesselgroup, parameters['adaption'], parameters['calcflow'], gdst))
  
  return_state, mean_flow_cap,std_flow_cap = adaption_cpp.computeAdaption(vesselgroup, parameters['adaption'], parameters['calcflow'], gdst)
  
  # creating opt criterion
  
  if return_state == 0:
    #it was convergent
#    veins, arteries, capillaries = getVesselTypes(gdst['vessels_after_adaption'])
#    edges, flows = krebsutils.read_vessels_from_hdf(gdst['vessels_after_adaption'],['flow'])
#    this_std = np.std(flows[capillaries])
#    #gdst.file.flush()
#    del gdst['recomputed']
#    del gdst['vessels_after_adaption']
    print('Adaption convergent!')
    print('mean(flow): %f, std(flow): %f ' % (mean_flow_cap,std_flow_cap))
    global calculated_mean_cap_flow
    calculated_mean_cap_flow = mean_flow_cap
    return std_flow_cap
  else:
    print('Adaption crashed!') #we dont like that, big penalty
    return sys.float_info.max

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
  
  vesselgroup = h5files.open(fn, 'r', search = False)[grp_pattern]
  f_opt_data = h5files.open('PSO_data_%s.h5' % basename(fn), 'a', search = False)
  
  adaptionParams['counter'] = 0
  #for default_pso use only 20, hope this speed up
  #for all other simulations 50 was used here
  adaptionParams['adaption']['max_nun_iterations'] = 100
  adaptionParams['adaption']['Q_refdot'] = 40
  adaptionParams['adaption']['S_0'] = 20
  # initialize parameters
  x_0= [1.5, 1.5, 1, 2500.]
  lb= [0.5, 0.5, 0.5, 500]
  ub= [4, 4, 4, 3500]
#  x_0= [1.5, 1.5, 1, 2500.,40]
#  lb= [0.5, 0.5, 0.5, 500,10]
#  ub= [4, 4, 4, 3500,80]
  x_0[0] =adaptionParams['adaption']['k_m']
  x_0[1] =adaptionParams['adaption']['k_c']
  x_0[2] =adaptionParams['adaption']['k_s']
  x_0[3] =adaptionParams['adaption']['cond_length']
  #x_0[4] =adaptionParams['adaption']['Q_refdot']
  
  if 0:#hope_final
    swarmsize = 50
    maxiter = 10
    minstep = 1e-8  #default 1e-8
    minfunc = 0.001  #default 1e-8
    omega = 0.5     #default 0.5
  if 0:#hope_final_phase
    adaptionParams['calcflow']['includePhaseSeparationEffect'] = True
    swarmsize = 50
    maxiter = 10
    minstep = 1e-8  #default 1e-8
    minfunc = 0.001  #default 1e-8
    omega = 0.5     #default 0.5
  if 0:#big_swarm
    swarmsize = 100
    maxiter = 20
    minstep = 1e-8  #default 1e-8
    minfunc = 0.001  #default 1e-8
    omega = 0.5     #default 0.5
    adaptionParams['starting_radii'] = 10.,
  if 0:#big_swarm2
    swarmsize = 500
    maxiter = 20
    minstep = 1e-8  #default 1e-8
    minfunc = 0.001  #default 1e-8
    omega = 0.5     #default 0.5
  '''
  omega : scalar
      Particle velocity scaling factor (Default: 0.5)
  phip : scalar
      Scaling factor to search away from the particle's best known position
      (Default: 0.5)
  phig : scalar
      Scaling factor to search away from the swarm's best known position
      (Default: 0.5)
  '''
  if 1:#big_swarm3
    swarmsize = 1000
    maxiter = 20
    minstep = 1e-8  #default 1e-8
    minfunc = 0.001  #default 1e-8
    omega = 0.5     #default 0.5
  if 0:#big_swarm3_phase
    swarmsize = 1000
    maxiter = 20
    minstep = 1e-8  #default 1e-8
    minfunc = 0.001  #default 1e-8
    omega = 0.5     #default 0.5
    adaptionParams['calcflow']['includePhaseSeparationEffect'] = True
  if 0:#big_swarm_phase
    adaptionParams['calcflow']['includePhaseSeparationEffect'] = True
    swarmsize = 100
    maxiter = 20
    minstep = 1e-8  #default 1e-8
    minfunc = 0.001  #default 1e-8
    omega = 0.5     #default 0.5
    adaptionParams['starting_radii'] = 10.,
  if 0:#bigbig
    adaptionParams['calcflow']['includePhaseSeparationEffect'] = True
    swarmsize = 250
    maxiter = 30
    minstep = 1e-8  #default 1e-8
    minfunc = 0.001  #default 1e-8
    omega = 0.5     #default 0.5
  
#  swarmParams1 = dict()
#  swarmParams1['name'] = 'compare_phase_big'
#  swarmParams1['swarmsize'] = 100
#  swarmParams1['maxiter'] = 100
#  swarmParams1['minstep'] = 1e-8  #default 1e-8
#  swarmParams1['minfunc'] = 0.001 #default 1e-8
#  swarmParams1['omega'] = 0.5     #default 0.
#  maxiter = 20
  processes = 1 #somehow only 1 supported, pickels problem
  adaptionParams['num_threads'] = 8
  
  args=(f_opt_data, vesselgroup, adaptionParams)
  xopt, fopt = pso(do_pSO_for_file,
                   x_0,
                   lb,
                   ub,
                   ieqcons=[con3],
                   debug=True,
                   swarmsize= swarmsize,
                   maxiter= maxiter,
                   minstep=minstep,
                   minfunc=minfunc,
                   omega=omega,
                   processes=processes,
                   use_initial_guess=use_initial_guess,
                   phig=0.25,
                   args=args)
  
  f_opt_data.attrs.create('xopt', data = xopt)
  f_opt_data.attrs.create('fopt', data = fopt)
  

  h5files.closeall() # just to be sure
  
def run(parameter_set_name, filenames, grp_pattern):
  print 'submitting ...', parameter_set_name
  usetumorparams = False
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
  parser = argparse.ArgumentParser(description='runs particle swarm optimization for all files in file in filenames')  
  parser.add_argument('AdaptionParamSet')  
  parser.add_argument('vesselFileNames', nargs='*', type=argparse.FileType('r'), default=sys.stdin, help='Vessel file to calculate')  
  parser.add_argument('grp_pattern',help='Where to find the vessel group in the file')  
  parser.add_argument('-a', '--analyze', help = 'loop through all files analyze data and make plot', default=False, action='store_true')
  
  goodArguments, otherArguments = parser.parse_known_args()
  qsub.parse_args(otherArguments)
  
  #create filename due to former standards
  filenames=[]
  for fn in goodArguments.vesselFileNames:
    filenames.append(fn.name)
    
  try:
    if not goodArguments.AdaptionParamSet in dir(parameterSetsAdaption):
      raise AssertionError('Unknown parameter set %s!' % goodArguments.AdaptionParamSet)
  except Exception, e:
      print e.message
      sys.exit(-1)
      
  factory = getattr(parameterSetsAdaption, goodArguments.AdaptionParamSet)
  factory['name'] = goodArguments.AdaptionParamSet
  run2(factory, filenames, goodArguments.grp_pattern)
      