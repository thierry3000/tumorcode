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
#! /usr/bin/env python
# -*- coding: utf-8 -*-
if __name__ == '__main__':
  import os.path, sys
  sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../..'))
  
import os, sys
import numpy as np
import h5py
from copy import deepcopy

import myutils
import dicttoinfo
import krebsutils
import extensions

def centerslice(a):
  return a[:,:,a.shape[2]/2]

class Measure(object):
  @staticmethod
  def obtain_dataset(f, name, shape):
    try:
      return f[name]
    except KeyError:
      return f.create_dataset(name, shape, dtype = np.float32)

  def __init__(self, outfilename, params):
    self.lastt = 0.
    self.outfilename = outfilename
    exposure_thresholds_reference_conc_in = params.pop('exposure_thresholds_reference_conc_in', 1.)
    exposure_thresholds_reference_conc_ex = params.pop('exposure_thresholds_reference_conc_ex', 1.)
    baselevels = np.asarray([ 0.01, 0.1, 0.25, 0.5, 0.75, 1., 10. ])
    self.exposure_thresholds = (baselevels * exposure_thresholds_reference_conc_ex,
                                baselevels * exposure_thresholds_reference_conc_in)
    self.moviefilename = params.pop('moviefilename', None)
    self.movieintervall = params.pop('movieintervall', None)
    self.nextmovie_t = 0
    self.movie_out_number = 0


  def __call__(self, *args, **kwargs):
    (t, _), conc = args
    """should measure locally at very short time steps
     *AUC,
     *max conc
     *exposure time for c_max greater than
      + c_1 \approx 0
      + c_2 > 0
      + c_i ...
      + c_n < c_inject
     conc[0] is extracellular intrinsic
     conc[1] is intracellular intrinsic
    """
    # this code is for recording integral and other metrics at short intervalls
    dt = t - self.lastt
    self.lastt = t
    with h5py.File(self.outfilename, 'r+') as f:
      shape = conc[0].shape
      g = f.require_group('measurements').require_group('drug_local_integral')
      for i, in_ex in enumerate(['ex','in']):
        ds = Measure.obtain_dataset(g, 'auc_'+in_ex, shape)
        ds[:] = np.asarray(ds) + dt * conc[i] # numerical integration as step function -> AUC metric
        ds = Measure.obtain_dataset(g, 'c_max_'+in_ex, shape)
        ds[:] = np.maximum(np.asarray(ds), conc[i]) # maximal conc. metric
      for i, in_ex in enumerate(['ex','in']):
        for j, exposure_conc in enumerate(self.exposure_thresholds[i]):
          ds = Measure.obtain_dataset(g, 'exposure_time_%02i_%s' % (j,in_ex), shape)
          ds[:] = np.asarray(ds) + dt * (conc[i] > exposure_conc)
          ds.attrs['exposure_conc'] = exposure_conc
          ds.attrs['exposure_conc_num'] = j
      g.attrs['tend'] = t
    # this code is for recording slices through the conc. distributions at short intervalls
    if self.moviefilename and self.movieintervall and (t - self.nextmovie_t >=  -0.1*dt):
      self.nextmovie_t += self.movieintervall
      with h5py.File(self.moviefilename, 'w' if self.movie_out_number==0 else 'r+') as f:
        g = f.require_group('frames')
        idx = g.attrs.get('num_frames', 0)
        g.attrs.modify('num_frames', idx+1)
        print 'write t = %f, idx = %i' % (t, idx)
        gt = g.create_group('%04i' % idx)
        gt.attrs['time'] = t
        gt.attrs['number'] = idx
        gt.create_dataset('conc_ex', data = centerslice(conc[0]), compression = 9, dtype = np.float32)
        gt.create_dataset('conc_in', data = centerslice(conc[1]), compression = 9, dtype = np.float32)
      self.movie_out_number += 1


def run_iffsim(params):
  params = deepcopy(params)
  outfilename = params.pop('fn_out')
  measure_params = params.pop('ift_measure', dict())

  krebsutils.set_num_threads(params.pop('num_threads', 1))
  krebsutils.run_iffsim(dicttoinfo.dicttoinfo(params), str(outfilename.encode('utf-8')), Measure(outfilename, measure_params))


if __name__ == '__main__':
  param_fn = sys.argv[1]
  if param_fn.endswith('.py'):
    import imp
    param_module = imp.load_source('param_module', param_fn)
    params = deepcopy(param_module.params)
  elif param_fn.endswith('.json'):
    import json
    with open(param_fn, 'r') as f:
      params = json.load(f)
  run_iffsim(params)

