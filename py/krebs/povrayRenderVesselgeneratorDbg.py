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

#import krebs.povrayRenderVessels
from krebs.povrayRenderVessels import  addVesselTree, cm_redblue
from krebs.povrayEasy import *
import krebs.povrayRenderSettings

import subprocess as sp
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


def colorfactory(vesselgraph):
  is_set = lambda flags_,flag: np.asarray(np.bitwise_and(flags_, flag), np.bool)
  
  edges = vesselgraph.edgelist
  num_nodes = len(vesselgraph.nodes['position'])
  flags = vesselgraph.edges['flags']
  data = vesselgraph.nodes['pressure']
  nflags = vesselgraph.nodes['nodeflags']
  #action = vesselgraph.nodes['action']
  #data = vesselgraph.edges['shearforce']
  #data = krebsutils.edge_to_node_property(num_nodes, edges, data, 'avg')  
  #nflags = krebsutils.edge_to_node_property(num_nodes, edges, flags, 'and')
  #ncirc  = krebsutils.edge_to_node_property(num_nodes, edges, np.bitwise_and(flags, krebsutils.CIRCULATED), 'or')
  #nflags = np.bitwise_or(nflags, ncirc)

  gray = np.asarray((0.1,0.1,0.1))
  ltgray = np.asarray((0.0, 0.8, 0.0))
  boundarycol = np.asarray((1., 1., 0))

  circulated = is_set(flags,krebsutils.CIRCULATED)
  has_circ = circulated.max()==True
  if 0:
    DIE = 0
    GROW = 1
    IDLE = 2
    edgecolors = np.zeros((len(edges),3), dtype = np.float)
    nodecolors = np.zeros((num_nodes,3), dtype = np.float)
    edgecolors[:,:] = gray
    nodecolors[:,:] = gray
    nodecolors[action == DIE] = np.asarray((1.,0.,0.))
    nodecolors[action == GROW] = np.asarray((0.,1.,0.))
    nodecolors[action == IDLE] = np.asarray((0.,0.,1.))
    
  elif has_circ:
    ncirculated = is_set(nflags,krebsutils.CIRCULATED)
    p0 = np.amin(data[ncirculated])
    p1 = np.amax(data[ncirculated])
    cm = matplotlib.cm.ScalarMappable(cmap=cm_redblue)
    cm.set_clim(p0, p1)

    edgedata = np.average((data[edges[:,0]], data[edges[:,1]]), axis=0)
    edgecolors = cm.to_rgba(edgedata)[:,:3]
    edgecolors[np.where(circulated==0)] = gray

    nodecolors = cm.to_rgba(data)[:,:3]
    nodecolors[np.where(ncirculated==0)] = gray
  else:
    edgecolors = np.zeros((len(edges),3), dtype = np.float)
    nodecolors = np.zeros((num_nodes,3), dtype = np.float)
    edgecolors[is_set(flags,krebsutils.ARTERY),0] = 1
    nodecolors[is_set(nflags,krebsutils.ARTERY),0] = 1
    edgecolors[is_set(flags,krebsutils.VEIN),2] = 1
    nodecolors[is_set(nflags,krebsutils.VEIN),2] = 1
  edgecolors[is_set(flags,krebsutils.CAPILLARY)] = ltgray
  nodecolors[is_set(nflags,krebsutils.CAPILLARY)] = ltgray
  nodecolors[is_set(nflags,krebsutils.BOUNDARY)] = boundarycol
  edgecolors[is_set(flags,krebsutils.BOUNDARY)] = boundarycol

  vesselgraph.edges['colors'] = edgecolors
  vesselgraph.nodes['colors'] = nodecolors


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
  num_hier = int(vessel_file['parameters'].attrs['num_hierarchical_iterations'])
  if os.path.isfile('dbgvesselsafter_initial_capillaries.png'):
    os.rename('dbgvesselsafter_initial_capillaries.png','dbgvesselsafter_initial_capillaries_0.png')
  if os.path.isfile('dbgvesselsafter_initial_growth.png'):  
    os.rename('dbgvesselsafter_initial_growth.png','dbgvesselsafter_initial_capillaries_1.png')
  
  print('num_hier read: %i' % num_hier)  
  stamp = datetime.datetime.now().time()
  stamp = 'a_stamp'
  fn_1 = 'out_initial_%s.webm' % stamp
  #frame  after_initial_growth
  #frame  after_initial_capillaries
  #frame  after_initial_calcflow
  #frame  after_remodel_ii_IIIII
  #fn_2 = 'out_after_initial_growth_%s.webm' % stamp
  fn_3 = 'out_after_initial_capillaries_%s.webm' % stamp
  fn_4 = 'out_after_initial_calcflow_%s.webm' % stamp
  fn_5 = 'out_after_remodel_hit_%s_%s.webm'
  fn_6 = 'out_capillaries_hit_%s_%s.webm'
  myframe_rate = str('2')
  if 0: #direct command line approach
    print("FFMPEG_BIN: %s" % FFMPEG_BIN)
    command = [ FFMPEG_BIN,
          '-y', # (optional) overwrite output file if it exists
          '-framerate', myframe_rate,
          '-i', 'dbgvesselsinitial_growth_%02d.png',
          '-an', # Tells FFMPEG not to expect any audio
          #'-c:v', 'libx264',
          #'-s', 'WxH',
          #'-r', '30',
          #'-pix_fmt', 'yuv420p'
          ]
    command.append(fn_initial_movie)
    sp.call(command)
  else: #using ffmpy wrapper
    if 1:
      list_of_filenames=[]
      optionstring_input = '-framerate %s -an' % myframe_rate
      #optionstring_output = '-c:v libx264 -pix_fmt yuv420p'
      #optionstring_output = '-c copy -f mpegts'
      optionstring_output = None
      ff_1 = third_party.ffmpy.FFmpeg(inputs={'dbgvesselsgrowth_initial_%02d.png':'-framerate 1 -an'},
                        outputs={fn_1:'-r 1 -y'})
      ff_1.run()
      list_of_filenames.append(fn_1)

      ff_3 = third_party.ffmpy.FFmpeg(inputs={'dbgvesselsafter_initial_capillaries_%1d.png':'-an -t 3 -loop 1'},
                        outputs={fn_3:'-r 5 -y'})
      ff_3.run()
      list_of_filenames.append(fn_3)
      
      ff_4 = third_party.ffmpy.FFmpeg(inputs={'dbgvesselsafter_initial_calcflow.png':'-loop 1 -an -t 4'},
                        outputs={fn_4:'-an -r 10 -y'})
      ff_4.run()
      list_of_filenames.append(fn_4)
      for i in range(num_hier+1):
        print('running hit %i'%i)
        
        ff_5 = third_party.ffmpy.FFmpeg(inputs={'dbgvesselsafter_remodel_0%i_%%05d.png'%i:'-framerate %s'%myframe_rate},
                        outputs={fn_5%(stamp,i):None})
        ff_5.run()
        del ff_5
        list_of_filenames.append(fn_5%(stamp,i))
        
        if i < num_hier:
          if os.path.isfile('dbgvesselswith_capillaries_hit_%i.png' % i):
            os.rename('dbgvesselswith_capillaries_hit_%i.png' % i, 'capillaries_0_hit_%i.png' % i)
          if os.path.isfile('dbgvesselswithout_capillaries_hit_%i.png' % i):
            os.rename('dbgvesselswithout_capillaries_hit_%i.png' % i, 'capillaries_1_hit_%i.png' % i)
          ff_6 = third_party.ffmpy.FFmpeg(inputs={'capillaries_%%1d_hit_%i.png' % i : '-an -t 3 -loop 1'},
                          outputs={fn_6%(stamp,i):'-r 5 -y'})
          ff_6.run()
          del ff_6
          list_of_filenames.append(fn_6%(stamp,i))
        
        
      with open('file.txt','w') as f:
        for fn in list_of_filenames:
          print('file %s' % fn, file=f)
      
    #merge
    if 1:
      ff = third_party.ffmpy.FFmpeg(inputs={'file.txt':'-f concat'},outputs={'construct_vessel_network.webm':'-codec copy -y'})
      ff.run()
    else:
      command = [ FFMPEG_BIN,
          '-y', # (optional) overwrite output file if it exists
          '-i', 'concat:%s|%s|%s|%s'%(fn_1,fn_2,fn_3,fn_4),
          '-an', # Tells FFMPEG not to expect any audio
          '-c', 'copy',
          #'-s', 'WxH',
          #'-r', '30',
          #'-pix_fmt', 'yuv420p'
          ]
      command.append('labber.mp4')
      print(command)
      sp.call(command)

if (__name__ == '__main__'):
  parser = argparse.ArgumentParser(description='Make snapshot of blood vessel network creation')  
  parser.add_argument('filename', type=argparse.FileType('r'), default=sys.stdin)
  parser.add_argument('debug_filename', type=argparse.FileType('r'), default=sys.stdin) 
  parser.add_argument('-g', '--grp_pattern', default='*/vessels')
  parser.add_argument('--noPov', default=False, action='store_true')
  parser.add_argument('-v', '--movie', default=False, action='store_true')

  goodArguments, otherArguments = parser.parse_known_args()
  #qsub.parse_args(otherArguments)
  dbg_fn = goodArguments.debug_filename.name
  #pattern = sys.argv[2]
  options = getattr(krebs.povrayRenderSettings,'dbg_vessels')
  f = h5py.File(dbg_fn,'r')
  dirs_list=[]
  myDebug=False
  dirs_list.append(myutils.walkh5(f['/'], 'growth_initial_*'))
  dirs_list.append(myutils.walkh5(f['/'], 'after_initial_growth')) 
  dirs_list.append(myutils.walkh5(f['/'], 'after_initial_capillaries')) 
  dirs_list.append(myutils.walkh5(f['/'], 'after_initial_calcflow')) 
  if myDebug:
    dirs_list.append(myutils.walkh5(f['/'], 'after_remodel_00_0000*'))
    dirs_list.append(myutils.walkh5(f['/'], 'after_remodel_01_0000*'))
  else:
    dirs_list.append(myutils.walkh5(f['/'], 'after_remodel_*'))
  dirs_list.append(myutils.walkh5(f['/'], 'with_capillaries_hit_*'))
  dirs_list.append(myutils.walkh5(f['/'], 'without_capillaries_hit_*'))
  
  dirs_list.append(myutils.walkh5(f['/'], 'growth_hit_*'))
  dirs_list.append(myutils.walkh5(f['/'], 'after_growth_hit_*'))
  dirs_list.append(myutils.walkh5(f['/'], 'after_capillaries_hit_*'))
  dirs_list.append(myutils.walkh5(f['/'], 'after_calcflow_hit_*'))
  
  for dirs in dirs_list:
    for grpname in dirs:
      vesselgrp = f[grpname]
      output_fn = splitext(basename(dbg_fn))[0]+posixpath.basename(vesselgrp.name)
      if not goodArguments.noPov:
        renderSliceWithDistribution(vesselgrp, output_fn ,options )
  
#  for grpname in dirs2:
#    vesselgrp = f[grpname]
#    output_fn = splitext(basename(fn))[0]+posixpath.basename(vesselgrp.name)
#    if debugWithPov and 1:
#      renderSliceWithDistribution(vesselgrp, output_fn ,options )
#  for grpname in dirs3:
#    vesselgrp = f[grpname]
#    output_fn = splitext(basename(fn))[0]+posixpath.basename(vesselgrp.name)
#    if debugWithPov and 1:
#      renderSliceWithDistribution(vesselgrp, output_fn ,options )
#  for grpname in dirs4:
#    vesselgrp = f[grpname]
#    output_fn = splitext(basename(fn))[0]+posixpath.basename(vesselgrp.name)
#    if debugWithPov and 1:
#      renderSliceWithDistribution(vesselgrp, output_fn ,options )
#  for grpname in dirs5:
#    vesselgrp = f[grpname]
#    output_fn = splitext(basename(fn))[0]+posixpath.basename(vesselgrp.name)
#    if debugWithPov and 0:
#      renderSliceWithDistribution(vesselgrp, output_fn ,options )
    
  if goodArguments.movie:
    create_movie(goodArguments)
    
  
