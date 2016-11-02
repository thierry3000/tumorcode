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

"""
  This is the main script in order to generate raytraced images of vessel networks.
  It submitts stuff via qsub.
  Can be configured for tumor runs as well as single vessel network files.
  It can make videos of the former.
"""
if __name__ == '__main__':
  import os.path, sys
  sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'..'))
  
from os.path import basename, splitext, join
import sys, os
import posixpath
import h5py
import h5files
import numpy
import krebsutils
from mystruct import Struct
from copy import deepcopy
import myutils
import identifycluster

def clientfunc(job, cwd):
  os.chdir(cwd)
  job.render()

def estimateRuntimeAndMemory(g):
  ''' time in hours. memory in mB '''
  if 'po2vessels' in g:
    N = g['po2vessels'].shape[1]
  elif 'conc' in g:
    N = len(g.parent['iff/vessels/edges/node_a_index'])
  else:
    if 'vessels' in g: g = g['vessels']
    N = len(g['edges/node_a_index'])
  t = (N*5. / 60.)/(100000)
  #m = (1200. * N) / (300000) #this seams to be to low for tumors!
  m = (2400. * N) / (300000)
  if m < 1000:
    m=1000
  return t, m


class RenderJob(object):
  def __init__(self, f, group_name, postfix, params = dict()):
    self.fn = f.filename
    self.postfix = postfix
    self.group_name = group_name
    if identifycluster.getname() == 'snowden':
      my_threads = 16
    else:
      my_threads = 4
    self.params = dict(
      res=(1400,1400),
      aa=4,
      colored_slice=True,
      out_alpha=True,
      num_threads=my_threads, 
      temp_file_dir = '/tmp',
      cam = 'topdown'
    )
    self.params.update(params)
    self.runtime_and_mem = estimateRuntimeAndMemory(f[self.group_name])

  @property
  def imageFilename(self):
    return '%s-%s%s.png' % (splitext(basename(self.fn))[0], myutils.sanitize_posixpath(self.group_name).replace('/','_'), ('_'+self.postfix) if self.postfix else '')

  def render(self):
    """ run povray. This should be called on the computing node. """
    with h5files.open(self.fn,'a') as f:
      if 'po2vessels' in f[self.group_name]:
        from krebs.povrayRenderOxygenDetailed import renderScene
        renderScene(f[self.group_name],
                    self.imageFilename,
                    self.params)

      elif 'tumor' in f[self.group_name]:
        from krebs.povrayRenderTumor import renderScene
        renderScene(f[self.group_name]['vessels'],
                    f[self.group_name]['tumor'],
                    self.imageFilename,
                    **self.params)
#      else:
#        from krebs.povrayRenderVessels import renderScene
#        renderScene(f[self.group_name],
#                    self.imageFilename,
#                    **self.params)
      #renderScene(drug_grp, imagefn, parameters)
      elif 'conc' in f[self.group_name]:
        from krebs.povrayRenderIff import renderScene
        renderScene(f[self.group_name],
                    self.imageFilename,
                    self.params)
      else:
        from krebs.povrayRenderVessels import render_different_data_types
        render_different_data_types(f[self.group_name],
                    **self.params )



if __name__ == '__main__':
  import qsub
  import argparse
  parser = argparse.ArgumentParser(description='Povray wrapper')
  
  parser.add_argument('vesselFileNames', nargs='*', type=argparse.FileType('r'), default=sys.stdin)
  parser.add_argument('grp_pattern')
  parser.add_argument("-d","--data", dest="datalist", help="which data (pressure, flow, shearforce, hematocrit) as comma separated list", default='pressure', action="store")
  parser.add_argument("-f","--filter-uncirculated", dest="filteruncirculated", help="filter uncirculated vessels", default=False, action="store_true")
  parser.add_argument("--filter-radius-high-pass", dest="filterradiushighpass", action="store", default = -1., type=float)
  parser.add_argument("--filter-radius-low-pass", dest="filterradiuslowpass", action="store",  default = -1., type=float)     
  parser.add_argument("--no-overlay", dest="overlay", default = True, action="store_false")
  parser.add_argument("--dpi", dest="dpi", default=None, action="store")
  parser.add_argument("--format", dest="format", default=None, action="store")
  parser.add_argument("-c","--cam", dest="cam", help="camera mode: topdown, pie, topdown_slice", default='topdown_slice', action="store", type=str)
  parser.add_argument("-u","--auc", dest="plot_auc", help="for area under curve, we have only a single timepoint", default=False, action="store_true")
  parser.add_argument("-a","--auto_color", dest="auto_colorscale", help="", default=False, action="store_true")  
  
  goodArguments, otherArguments = parser.parse_known_args()
  qsub.parse_args(otherArguments)
  
  #options, args = parser.parse_args()
  parameters = dict(
    cam = goodArguments.cam,
    datalist = map(lambda s: s, map(str.strip, goodArguments.datalist.split(','))),
    filteruncirculated = goodArguments.filteruncirculated,
    filterradiushighpass = goodArguments.filterradiushighpass,
    filterradiuslowpass = goodArguments.filterradiuslowpass,
    overlay = goodArguments.overlay,
    plot_auc = goodArguments.plot_auc,
    auto_colorscale = goodArguments.auto_colorscale,
  )

  #create filename due to former standards
  filenames=[]
  for fn in goodArguments.vesselFileNames:
    filenames.append(fn.name)
  pattern = goodArguments.grp_pattern
  postfix = ''

  jobs = []
  a = jobs.append

  for fn in filenames:
    with h5files.open(fn,'r') as f:
      paths = myutils.walkh5(f, pattern)
      if goodArguments.plot_auc:
        path = myutils.walkh5(f, 'out0000')
        j = RenderJob(f, path[0], postfix, parameters)
      for path in paths:
        j = RenderJob(f, path, postfix, parameters)
        a(j)

  for job in jobs:
    t, m = job.runtime_and_mem
    print 'submit %s, %i mb, %f h' % (job.imageFilename, m, t)
    qsub.submit(
      qsub.func(clientfunc, job, os.getcwd()),
      name='job_render_'+basename(job.imageFilename),
      num_cpus=job.params['num_threads'],
      mem=('%iMB' % m),
      days=t/24.)
