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
from os.path import basename
import qsub
import krebsutils
from krebs.vesselgenerator import *
import krebsjobs.parameters.vesselGenParams as vParams

from krebsutils import typelist
import identifycluster

def run_vesselgen_client(configstring_file, workdir, vesselfilename):
  qsub.printClientInfo()
  os.chdir(workdir)
  with open(configstring_file, 'r') as f:
    configstring = f.read()
  krebsutils.run_vesselgen(configstring)

## if this is imported on a cluster node we are done!
#if not qsub.is_client:
import copy

#qsub.parse_args(sys.argv)
 
computer_relative_speed = os.environ.get('RELATIVE_COMPUTATIONAL_SPEED_TO_HOME', 1)
runtime_per_iterdof = 5.e-3 / computer_relative_speed

if identifycluster.name == 'sleipnir':
  runtime_per_iterdof = 50.e-3
elif identifycluster.name == 'durga':
  runtime_per_iterdof = 5.e-3
  num_threads = 6
elif identifycluster.name == 'snowden':
  runtime_per_iterdof = 2.5e-3
  num_threads = 8
else:
  runtime_per_iterdof = 8.e-6
  num_threads = 4


def runtime_sec(vd):
  """
    compute runtime as machine dependent factor x lattice sites x number of iterations
    returns seconds
  """
  #ndofs = vd.shape[0]*vd.shape[1]*vd.shape[2]
  #return float(runtime_per_iterdof) * ndofs * vd.num_iter * (1+vd.num_hierarchical_iterations)
  res = 0.
  for i in range(vd.num_hierarchical_iterations):
    volume = (vd.shape[0]*vd.shape[1]*vd.shape[2])*((2**3)**i)
    niter = max(q for q  in vd.shape)*100
    res += niter * volume
  return float(runtime_per_iterdof) * res


def run_config_samples(configfactory, num_samples, client_worker_function):
  vdcopy = configfactory(0)
  vdcopy.roots = []
  print vdcopy.generate_info_string()
  for sample_num in num_samples:
    vd = configfactory(sample_num)
    vd.outfilename += '-sample%02i' % sample_num
    vd.num_threads = num_threads
    print 'submitting %s, estimated runtime %f h, %i iters' % (vd.outfilename, runtime_sec(vd)/60/60, vd.num_iter)
    fn = vd.outfilename + '.info'
    s  = vd.generate_info_string()
    with open(fn, 'w') as f:
      f.write(s)
    qsub.submit(qsub.func(client_worker_function, fn, os.getcwd(), vd.outfilename+'.h5'),
                name = 'job_'+basename(vd.outfilename),
                days = (runtime_sec(vd)/60/60)/24.,
                mem="2000MB", num_cpus=num_threads)


def write_config(configfactory):
  vd = configfactory(0)
  s = vd.generate_info_string()
  fn = vd.outfilename + '.info'
  f = open(fn, 'w')
  f.write(s)
  f.close()
  return fn


configuration_type_factories = dict([
  ('typeA', lambda vd, i: generate_onesided_alternating_config(vd)),
  ('typeB', lambda vd, i: generate_alternating_config(vd, 1)),
  ('typeC', lambda vd, i: generate_alternating_config(vd, 2)),
  ('typeD', lambda vd, i: generate_handmade1_config(vd, configs1[0])),
  ('typeE', lambda vd, i: generate_handmade1_config(vd, configs1[1])),
  ('typeF', lambda vd, i: generate_handmade1_config(vd, configs1[6])),
  ('typeG', lambda vd, i: generate_handmade1_config_w_stems(vd, configs1dict['baumlecfg11'], vd.shape[0]*0.5, vd.shape[0]*0.8)),
  ('typeH', lambda vd, i: generate_handmade1_config_w_stems(vd, configs1dict['baumlecfg12'], vd.shape[0]*0.5, vd.shape[0]*0.8)),
  ('typeI', lambda vd, i: generate_alternating_diluted_config(vd))
])


  ##==============================================================================
  ##                       main function
  ##==============================================================================


if (not qsub.is_client) and __name__=="__main__":
  import argparse
  parser = argparse.ArgumentParser(prog='submitVesselgeneration',description='Compute an artificial blood vessel network.')
  #default is False
  parser.add_argument('-t','--type', nargs='+', help='Type of root node configurations, range from 0 to 8', default=[8], type=int)
  parser.add_argument('-p','--VesselParamSet', help='specify Parameterset for creation, possible configs are found in /krebsjobs/parameters/vesselGenParams.py')
  parser.add_argument('-w','--width_cube', help='Size of vesselcube you want to create', default=1000, type=float)
  parser.add_argument('--twoD', help='Creates only 2D network, default is 3D', default=False, action='store_true')  
  parser.add_argument('-i','--iter_h', help='Number of hieracial iterations', default=1, type=int)  
  parser.add_argument('-e','--ensemble_size', help='Number of realizations', default=1, type=int)  
  
  
  goodArguments, otherArguments = parser.parse_known_args()
  qsub.parse_args(otherArguments)

  # first define a function which does something with a configuration
  # either run a job or write a config file (for debug)
  def doit(name_, type, shape, hiter, options_):
   def real_factory(i):
     options = copy.deepcopy(options_)
     name = '%s-%s' % (name_, type)
     vd = VD(shape=shape, scale=options.pop('scale'), latticetype='fcc', num_hierarchical_iterations=hiter, name=name, ensemble_index=i, **options)
     vd = fix_worldsize(vd)
     vd.outfilename = 'vessels-'+make_fn(vd)
     vd.num_threads = num_threads
     vd = configuration_type_factories[type](vd, i)
     return vd
   run_config_samples(real_factory, index_range, run_vesselgen_client)

  
  
  index_range = range(0, goodArguments.ensemble_size)
  if goodArguments.VesselParamSet:
    #check if specified paramSet exists
    try:
      if not goodArguments.VesselParamSet in dir(vParams):
        raise AssertionError('Unknown parameter set %s!' % goodArguments.VesselParamSet)
      else:
        factory = getattr(vParams, goodArguments.VesselParamSet)
        nums_points = int(goodArguments.width_cube / (2**goodArguments.iter_h * factory['scale']) + 1)
      if nums_points<5:
        print('nums_points: %i'% nums_points)
        raise AssertionError('This width might be a little small!')
    except Exception, e:
      print e.message
      sys.exit(-1)
    
    #for t in 'typeA typeB typeC typeD typeE typeF typeG typeH typeI'.split():
    for type_index in goodArguments.type:
      #default nums_points for 3D defined here
      if goodArguments.twoD: # 2D flag enabled
        print('2D forced!')
        doit(goodArguments.VesselParamSet, typelist[type_index], (nums_points,nums_points,1), goodArguments.iter_h, factory)
      else:
        doit(goodArguments.VesselParamSet, typelist[type_index], (nums_points,nums_points,nums_points), goodArguments.iter_h, factory)
  else: #this was used befor the release and is here for compatibilty reasons
    print('Falling back to default configs from .py file!')
    print('No VesselParamSet provided')
    if 0:
      doit('4mm-P3-q2d', 'typeA', (7,7,2), 2, vParams.paramset3)
  
    if 0:
      doit('q2d-30mm-P2', 'typeA', (13, 13, 2), 4, vParams.paramset2)
      doit('q2d-30mm-P2', 'typeB', (13, 13, 2), 4, vParams.paramset2)
      doit('q2d-30mm-P2', 'typeC', (13, 13, 2), 4, vParams.paramset2)
      doit('q2d-30mm-P2', 'typeD', (13, 13, 2), 4, vParams.paramset2)
      doit('q2d-30mm-P2', 'typeE', (13, 13, 2), 4, vParams.paramset2)
      doit('q2d-30mm-P2', 'typeF', (13, 13, 2), 4, vParams.paramset2)
  
    if 0:
      for t in 'typeA typeB typeC typeD typeE typeF typeG typeH'.split():
        doit('3d-8mm-P2', t, (14, 14, 14), 2, vParams.paramset2)
  
    if 0:
      for t in 'typeA typeB typeC typeD typeE typeF typeG typeH'.split():
        doit('q2d-8mm-P3', t, (14, 14, 2), 2, vParams.paramset3)
        doit('3d-8mm-P3', t, (14, 14, 14), 2, vParams.paramset3)
  
    if 0:
      for t in 'typeA typeB typeC typeD typeE typeF typeG typeH'.split():
        doit('q2d-8mm-P4', t, (14, 14, 2), 2, vParams.paramset4)
        doit('3d-8mm-P4', t, (14, 14, 14), 2, vParams.paramset4)
  
    if 0:
      for t in 'typeA typeB typeC typeD typeE typeF typeG typeH'.split():
        doit('q2d-15mm-P4', t, (8, 8, 1), 3, vParams.paramset11)
  
    if 0:
      for t in 'typeA typeB typeC typeD typeE typeF typeG typeH'.split():
        #doit('q2d-8mm-P6', t, (13, 13, 2), 2, paramset6)
        doit('3d-8mm-P7', t, (13, 13, 13), 2, vParams.paramset7)
  
    if 0:
      #for t in 'typeD'.split():
      for t in 'typeA typeB typeC typeD typeE typeF typeG typeH typeI'.split():
        doit('r2d-P12-5mm', t, (21, 21, 1), 2, vParams.paramset12)
        #doit('3d-8mm-P10', t, (13, 13, 13), 2, paramset10)
    if 0:
      for t in 'typeI'.split():
        doit('r2d-8mm-P10', t, (13, 13, 1), 1, vParams.paramset12)
        
    if 0:
      #for t in 'typeA typeB typeC typeD typeE typeF typeG typeH typeI'.split():
      for t in 'typeA'.split():
        #doit('q2d-mini-2layer-H2-P10', t, (14, 14, 2), 2, paramset10)
        doit('3d-bigger-H2-P10', t, (17, 17, 17), 2, vParams.paramset10)
        
    if 0:
      for t in 'typeA typeB typeC typeD typeE typeF typeG typeH typeI'.split():
        doit('3d-mini-P10', t, (14, 14, 14), 1, vParams.paramset10)
    if 0:
      for t in 'typeA typeB typeC typeD typeE typeF typeG typeH typeI'.split():
        doit('3d-mini-mini-P10', t, (7, 7, 7), 1, vParams.paramset10)
    if 0:
      for t in 'typeA typeB typeC typeD typeE typeF typeG typeH typeI'.split():
      #for t in 'typeI'.split():
        #doit('3d-swine-H1-P21', t, (8, 8, 8), 1, paramset21)
        #doit('3d-swine-H2-P21', t, (8, 8, 8), 2, paramset21)
        doit('3d-big-swine-H2-P24', t, (16, 16, 16), 2, vParams.paramset24)
    if 0:# adaption configuration
      for t in 'typeA typeB typeC typeD typeE typeF typeG typeH typeI'.split():
        doit('3d-P11', t, (9, 9, 9), 2, vParams.paramset11)
    if 0:# adaption configuration
      for t in 'typeA typeB typeC typeD typeE typeF typeG typeH typeI'.split():
        doit('q2d-P11', t, (9, 9, 2), 2, vParams.paramset11)
    if 0:# trastuzumab configurations --> MVD \approx 100
      #from submitVesseltreeCalibration.py
      hiter = 2
      cube_width  = 4000.
      scale = 88 #see run2 results from submitVesseltreeCalibration.py
      paramset13['scale']=scale
      nums_points = int(cube_width / (2**hiter * scale) + 1)
      for t in 'typeA typeB typeC typeD typeE typeF typeG typeH typeI'.split():
        doit('vessTras-P13', t, (nums_points,nums_points,nums_points), hiter, vParams.paramset13)
    if 0:# trastuzumab configurations --> MVD \approx 100
      #from submitVesseltreeCalibration.py
      hiter = 2
      cube_width  = 10000.
      scale = 88 #see run2 results from submitVesseltreeCalibration.py
      paramset13['scale']=scale
      nums_points = int(cube_width / (2**hiter * scale) + 1)
      for t in 'typeA typeB typeC typeD typeE typeF typeG typeH typeI'.split():
        doit('vessBigTras-P13', t, (nums_points,nums_points,nums_points), hiter, vParams.paramset13)
    
    if 0:# mini_test
      #for t in 'typeA typeB typeC typeD typeE typeF typeG typeH typeI'.split():
      for t in 'typeI'.split():
        doit('q2d_mini_test_Pdefault', t, (7, 7, 3), 1, vParams.default)
#  else:
#  #this is the tutorial testing case quasi2D
#    index_range = range(0, 1)
#    for t in 'typeI'.split():
#      doit('tutorialVessels', t, (5, 5, 5), 1, vParams.default)