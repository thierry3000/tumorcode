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
import os,sys
from os.path import basename, dirname, join, splitext, commonprefix
if __name__=='__main__': sys.path.append(join(dirname(__file__),'..'))
import h5py
import numpy as np
import extensions # for hdf5 support in np.asarray
import krebsutils
import myutils
import posixpath
import math
import collections


import matplotlib
import matplotlib.cm
import matplotlib.pyplot as pyplot
import mpl_utils


from krebs.analyzeGeneral import DataVesselGlobal, DataTumorTissueSingle, DataDistanceFromCenter, DataBasicVessel, DataVesselSamples, DataVesselRadial, BinsSpecRange, BinsSpecArray, obtain_distmap_, generate_samples, combineSamples, HdfCacheRadialDistribution, CalcPhiVessels, calc_distmap


filename, pattern = sys.argv[1], sys.argv[2]


dataman = myutils.DataManager(20, [ DataTumorTissueSingle(), DataDistanceFromCenter(), DataBasicVessel(), DataVesselSamples(), DataVesselRadial(), DataVesselGlobal()])
f = h5files.open(filename)
group = f[pattern]
gvessels, gtumor = group['vessels'], group['tumor']
tumor_ld = dataman.obtain_data('ld', gtumor.file)

position    = dataman.obtain_data('basic_vessel_samples', 'position', gvessels, 30.)
weight_smpl = dataman.obtain_data('basic_vessel_samples', 'weight', gvessels, 30.)
flags       = dataman.obtain_data('basic_vessel_samples', 'flags', gvessels, 30.)
  
dist_smpl, distmap, mask, ld   = dataman.obtain_data('distancemap_samples', gvessels, gtumor, 30., 'radial', None)

if 0:
  fig = pyplot.figure()
  ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
  plt = ax.scatter(position[:,0], position[:,1], c = weight_smpl, s = 5, edgecolor=None, cmap = matplotlib.cm.jet, marker = 'o')
  ax.set_axis_bgcolor('k')
  fig.colorbar(plt)
  pyplot.show()


mask = mask & myutils.bbitwise_and(flags, krebsutils.CIRCULATED)
dist_smpl   = dist_smpl[mask]
weight_smpl = weight_smpl[mask]

bins = BinsSpecRange(-10000., 10000., 100.).arange()
a = myutils.MeanValueArray.fromHistogram1d(bins, dist_smpl, weight_smpl, weight_smpl) #np.ones_like(dist_smpl)
b = myutils.MeanValueArray.fromHistogram1d(bins, distmap.ravel(), np.ones_like(distmap.ravel()))
a.cnt = b.cnt.copy()
a.sum *= 1./(tumor_ld.scale**3)
a.sqr *= a.sum**2

c = a.sum / (bins[1:]**2-bins[:-1]**2)

fig, axes = pyplot.subplots(2,1)
axes[0].plot(bins[1:], a.avg*1.e6)
axes[1].plot(bins[1:], c*1.e6)
pyplot.show()