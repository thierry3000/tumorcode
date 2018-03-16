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
from __future__ import print_function

""" DESCRIPTION 
This tool is designed to debug the vessel generator.
You have to switch the flag 

full_debug_output = True,

in the parameterSetsVesselGen.py set to True in order to produce
a file called dbgvessels.h5
there all intermediate steps from the vessel creation process 
are stored.
"""

if __name__ == '__main__':
  import os.path, sys
  sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'..'))

import os,sys
from os.path import join, basename, dirname, splitext
import time
import krebsutils
import h5py
import posixpath
import numpy as np
import matplotlib.cm
import myutils
from copy import deepcopy
import math
import argparse
import qsub
import third_party
import third_party.ffmpy
import collections
import shutil

#import krebs.povrayRenderVessels
from krebs.povrayRenderVessels import  addVesselTree, cm_redblue
from krebs.povrayEasy import *
import krebs.povrayRenderSettings
from krebsjobs.submitPovrayRender import RenderJob  
from krebsjobs.submitPovrayRender import clientfunc

import datetime
import identifycluster
if(identifycluster.getname()=='snowden'):
  FFMPEG_BIN = "ffmpeg_latest"
else:
  FFMPEG_BIN = "ffmpeg" # on Linux ans Mac OS


colors = [(0., 0.2, 0., 0.), (1., 0., 1.0, 0.) ]
cm_gf = matplotlib.colors.LinearSegmentedColormap('', {
  'red'   : [(x, r, r) for (x,r,g,b) in colors],
  'green' : [(x, g, g) for (x,r,g,b) in colors],
  'blue'  : [(x, b, b) for (x,r,g,b) in colors],
}, N = 256, gamma = 1.)





def renderSliceWithDistribution(vesselgroup, imagefn, options):
  vessel_ld = krebsutils.read_lattice_data_from_hdf(vesselgroup['vessels/lattice'])
  vessel_graph = krebsutils.read_vessels_from_hdf(vesselgroup['vessels'], ['position', 'flags', 'radius', 'pressure', 'shearforce', 'nodeflags'], return_graph=True)
  vessel_graph.edges['radius'] *= 4.

  #kwargs = deepcopy(kwargs)
  wbbox = vessel_ld.worldBox
  trafo = calc_centering_normalization_trafo(wbbox)
  height = (wbbox[5]-wbbox[4])*trafo.w

  print('Vessel BBox:' + str(vessel_ld.worldBox ))
  print(vessel_ld)
  print('Post Trafo Ld BBox:' + str(transform_ld(trafo, vessel_ld).worldBox))

  hasGfField = 'field_ld' in vesselgroup and vessel_ld.shape[2] == 1
  #hasGfField = False
  if hasGfField:
    volume_ld = krebsutils.read_lattice_data_from_hdf(vesselgroup['field_ld'])
    print ('Volume BBox:' + str(volume_ld.worldBox))
    print (volume_ld)
    volumedata = np.asarray(vesselgroup['gf'])
    #print volumedata.shape, volumedata.min(), volumedata.max()
    volume_ld = transform_ld(trafo, volume_ld)
    print ('Post Trafo Volume BBox:'+ str(volume_ld.worldBox))
    print ('volume data bounds:'+ str(volumedata.min()), str(volumedata.max()))

  colorfactory(vessel_graph)

  with EasyPovRayRender(**options) as epv:
    epv.setBackground(options.pop('background',0.0))

    cam_fov = 60.
    cam_distance_factor = ComputeCameraDistanceFactor(cam_fov, options['res'], wbbox)

    epv.setCamera((0,0,cam_distance_factor*1.05), lookat = (0,0,0), fov = cam_fov, up = 'y')

    epv.addLight(10.*Vec3(1,0.5,2), 1.2)

    cm = matplotlib.cm.ScalarMappable(cmap = cm_gf)
    cm.set_clim(-0.01,1.01)
    pvcm = matplotlibColormapToPovray('DATACOLORMAP', cm)
    epv.declareColorMap(pvcm)

    if hasGfField:
      volumedata = epv.declareVolumeData(volumedata, volume_ld.worldBox, volume_ld)
      epv.addVolumeDataSlice(volumedata, (0,0,0), (0, 0, 1.), pvcm)
    addVesselTree(epv, vessel_graph, trafo = trafo, **options)

    imagefn = epv.render(imagefn)

def create_movie(options):
  vessel_file = h5py.File(options.filename.name,'r')
  extension = '-vessels_pressure'
  num_hier = int(vessel_file['parameters'].attrs['num_hierarchical_iterations'])
  if not os.path.isfile('dbgvessels_after_initial_capillaries_0%s.png' % extension):
    if os.path.isfile('dbgvessels_after_initial_capillaries%s.png'% extension):
      shutil.copyfile('dbgvessels_after_initial_capillaries%s.png'% extension,'dbgvessels_after_initial_capillaries_0%s.png'% extension)
  if not os.path.isfile('dbgvessels_after_initial_capillaries_1%s.png'% extension):
    if os.path.isfile('dbgvessels_after_initial_growth%s.png'% extension):
      shutil.copyfile('dbgvessels_after_initial_growth%s.png'% extension,'dbgvessels_after_initial_capillaries_1%s.png'% extension)
  print('num_hier read: %i' % num_hier)  
  stamp = datetime.datetime.now().time()
  stamp = 'a_stamp'
  fn_1 = 'out_initial_%s.webm' % stamp
  fn_3 = 'out_after_initial_capillaries_%s.webm' % stamp
  fn_4 = 'out_after_initial_calcflow_%s.webm' % stamp
  fn_5 = 'out_after_remodel_hit_%s_%s.webm'
  fn_6 = 'out_capillaries_hit_%s_%s.webm'
  fn_7 = 'out_growth_hit_%s_%s.webm'
  fn_8 = 'out_capillaries_after_growth_hit_%s_%s.webm'
  myframe_rate = str('35')
  mydebug = True
  
  #using ffmpy wrapper
  list_of_filenames=[]
  optionstring_input = '-framerate %s -an' % myframe_rate
  #optionstring_output = '-c:v libx264 -pix_fmt yuv420p'
  #optionstring_output = '-c copy -f mpegts'
  optionstring_output = None
  ff_1 = third_party.ffmpy.FFmpeg(inputs={'dbgvessels_growth_initial_%%02d%s.png'% extension:'-framerate 1 -an'},
                    outputs={fn_1:'-r 1 -y'})
  ff_1.run()
  list_of_filenames.append(fn_1)

  ff_3 = third_party.ffmpy.FFmpeg(inputs={'dbgvessels_after_initial_capillaries_%%1d%s.png'% extension:'-an -t 3 -loop 1'},
                    outputs={fn_3:'-r 5 -y'})
  ff_3.run()
  list_of_filenames.append(fn_3)
  
  ff_4 = third_party.ffmpy.FFmpeg(inputs={'dbgvessels_after_initial_calcflow%s.png'% extension:'-loop 1 -an -t 4'},
                    outputs={fn_4:'-an -r 10 -y'})
  ff_4.run()
  list_of_filenames.append(fn_4)
  isInitialStep = True
  for i in range(num_hier+1):
    print('running hit %i'%i)
    
    ff_5 = third_party.ffmpy.FFmpeg(inputs={'dbgvessels_after_remodel_0%i_%%05d%s.png'%(i,extension):'-framerate %s'%myframe_rate},
                    outputs={fn_5%(stamp,i):'-y'})
    if mydebug:
      ff_5.run()
    del ff_5
    list_of_filenames.append(fn_5%(stamp,i))
    
    if i < num_hier: # if there is a hierachical step and this continues
      #this is inital growth not hierachical
      if os.path.isfile('dbgvessels_with_capillaries_hit_%i%s.png' % (i,extension)):
        #os.rename('dbgvesselswith_capillaries_hit_%i.png' % i, 'capillaries_0_hit_%i.png' % i)
        shutil.copyfile('dbgvessels_with_capillaries_hit_%i%s.png' % (i,extension), 'capillaries_0_hit_%i%s.png' % (i,extension))
      if os.path.isfile('dbgvessels_without_capillaries_hit_%i%s.png' % (i,extension)):
        #os.rename('dbgvesselswithout_capillaries_hit_%i.png' % i, 'capillaries_1_hit_%i.png' % i)
        shutil.copyfile('dbgvessels_without_capillaries_hit_%i%s.png' % (i,extension), 'capillaries_1_hit_%i%s.png' % (i,extension))
      ff_6 = third_party.ffmpy.FFmpeg(inputs={'capillaries_%%1d_hit_%i%s.png' % (i,extension) : '-an -t 3 -loop 1'},
                      outputs={fn_6%(stamp,i):'-r 5 -y'})
      if mydebug:
        ff_6.run()
      list_of_filenames.append(fn_6%(stamp,i))
      del ff_6
      
      if os.path.isfile('dbgvessels_growth_hit_%i_00%s.png' % (i,extension)):
        ff_7 = third_party.ffmpy.FFmpeg(inputs={'dbgvessels_growth_hit_%i_%%2d%s.png' % (i,extension) : '-framerate 1 -an'},
                      outputs={fn_7%(stamp,i):'-r 5 -y'})
        ff_7.run()
        del ff_7
      list_of_filenames.append(fn_7%(stamp,i))
      
      if os.path.isfile('dbgvessels_after_growth_hit_%i%s.png' % (i,extension)):
        shutil.copyfile('dbgvessels_after_growth_hit_%i%s.png' % (i,extension), 'capillaries_after_growth_0_hit_%i%s.png' % (i,extension))
      if os.path.isfile('dbgvessels_after_capillaries_hit_%i%s.png' % (i,extension)):
        shutil.copyfile('dbgvessels_after_capillaries_hit_%i%s.png' % (i,extension), 'capillaries_after_growth_1_hit_%i%s.png' % (i,extension))
      if os.path.isfile('dbgvessels_growth_hit_%i_00%s.png' % (i,extension)):
        ff_8 = third_party.ffmpy.FFmpeg(inputs={'capillaries_after_growth_%%1d_hit_%i%s.png' % (i,extension) : '-an -t 3 -loop 1'},
                      outputs={fn_8%(stamp,i):'-r 5 -y'})
        ff_8.run()
        del ff_8
      list_of_filenames.append(fn_8%(stamp,i))
    
    
  with open('file.txt','w') as f:
    for fn in list_of_filenames:
      print('file %s' % fn, file=f)
    
  #merge
  ff = third_party.ffmpy.FFmpeg(inputs={'file.txt':'-f concat'},outputs={'construct_vessel_network.webm':'-codec copy -y'})
  ff.run()

if (__name__ == '__main__'):
  """ note: when using the submitPovray.py the
            -f filter options is not applicable for
            the initial stage since the pressure is not yet calculated
  """
  parser = argparse.ArgumentParser(description='Make snapshot of blood vessel network creation')  
  parser.add_argument('filename', type=argparse.FileType('r'), default=sys.stdin)
  parser.add_argument('debug_filename', type=argparse.FileType('r'), default=sys.stdin) 
  parser.add_argument('-g', '--grp_pattern', default='*/vessels')
  parser.add_argument('--noPov', default=False, action='store_true')
  parser.add_argument('-v', '--movie', default=False, action='store_true')

  goodArguments, otherArguments = parser.parse_known_args()
  qsub.parse_args(otherArguments)
  dbg_fn = goodArguments.debug_filename.name
  options = getattr(krebs.povrayRenderSettings,'dbg_vessels')
  f = h5py.File(dbg_fn,'r')
  if not goodArguments.noPov:
    jobs = []
    a = jobs.append
  
    #for fn in filenames:
    with h5files.open(dbg_fn,'r') as f:
      paths = myutils.walkh5(f, '*/vessels')
      for path in paths:
        j = RenderJob(f, path, '', options)
        a(j)
  
    for job in jobs:
      t, m = job.runtime_and_mem
      print('submit %s, %i mb, %f h' % (job.imageFilename, m, t))
      qsub.submit(
        qsub.func(clientfunc, job, os.getcwd()),
        name='job_render_'+basename(job.imageFilename),
        num_cpus=job.params['num_threads'],
        mem=('%iMB' % m),
        days=t/24.)
  
    
  if goodArguments.movie:
    time.sleep(30) #wait 30 seconds for the queing system to finish
    create_movie(goodArguments)
    
  
