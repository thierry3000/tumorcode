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
#!/usr/bin/env python2
# -*- coding: utf-8 -*-
if __name__ == '__main__':
  import os.path, sys
  sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'..'))
import os,sys
from os.path import basename, splitext
import h5py
import numpy as np
import extensions # for asarray with h5py support
import krebsutils
import posixpath
import math
import glob
from collections import defaultdict
from pprint import pprint

import matplotlib
import matplotlib.pyplot as pyplot
import mpl_utils


import myutils

from analyzeGeneral import DataTumorTissueSingle, BinsSpecArray, DataVesselRadial, DataVesselSamples, DataBasicVessel, DataDistanceFromCenter



### ---------------- ------- -----------------------
### cached mvd/volume density measurements
### ---------------- ------- -----------------------
#class DataBloodVolume(object):
#  keywords = [
#    'phi_vessels_precise',  # map of the local volume fraction (with high accuracy)
#    'phi_vessel_per_region', # vd in bins
#    'mvd_per_regions', # mvd in bins
#  ]
#
#  def __init__(self, fm_filename):
#    self.fm_filename = fm_filename
#
#  def obtain_data(self, dataman, dataname, *args):
#    f, group = args[0], args[1]
#    fileid = 'file'+myutils.checksum(basename(f.filename))
#    obtain_data = lambda *args: dataman.obtain_data(args[0], f, *args[1:])
#    ld = obtain_data('ld')
#    #####
#    if dataname == 'phi_vessels_precise':
#      def read(gmeasure, name):
#        return np.asarray(gmeasure[name])
#      def write(gmeasure, name):
#        phi_vessels = calc_phi_vessels(f[group]['vessels'], ld, scaling = 1., samples_per_cell = 10)
#        gmeasure.create_dataset(name, data = phi_vessels, compression = 9)
#      with h5py.File(self.fm_filename, 'a') as fm:
#        return myutils.hdf_data_caching(read, write, fm, (fileid, group, 'phi_vessels_precise'), (0,0,1))
#    #####
#    if dataname == 'phi_vessel_per_region':
#      bin_spec = args[2]
#      distances = bin_spec.arange()
#      volel = (ld.scale**3)
#      dist = obtain_data('fieldvariable','dist_tumor', group)
#      phi = obtain_data('phi_vessels_precise', group)
#      res = []
#      for i in xrange(len(distances)-1):
#        a, b = distances[i], distances[i+1]
#        mask = np.logical_and(dist > a, dist <= b)
##        pyplot.imshow(imslice(mask * phi))
##        pyplot.show()
#        q = np.sum(phi[mask])*volel
#        masked_dist = dist[mask]
#        res.append((q, np.sum(mask)*volel, np.amin(masked_dist), np.amax(masked_dist)))
#      return np.asarray(res)
##    ####
##    if dataname == 'mvd_per_regions':
##      bins = args[2]
##      paramid = 'PARAMID'+myutils.checksum(bins)
##      def read(gmeasure, name):
##        return myutils.MeanValueArray.read(gmeasure, name)
##      def write(gmeasure, name):
##        dist = obtain_data('fieldvariable','dist_tumor', group)
##        ld = obtain_data('ld')
##        a = calcMVDAsLineDensity(f[group]['vessels'], dist, ld, bins)
##        a.write(gmeasure, name)
##        dataman('basic_vessel_radial', 'mvd', vesselgroup, tumorgroup, 10., bins_spec, distance_distribution_name, None, cachelocation)
##      with h5py.File(self.fm_filename, 'a') as fm:
##        return myutils.hdf_data_caching(read, write, fm, (fileid, group, 'mvd_per_regions', paramid), (0,0,2,0))

## ---------------- ------- -----------------------
## averaging utility
## ---------------- ------- -----------------------
def obtain_averaged_phi(dataman, fmeasure, groups, bins_spec, distancemap_spec):
  l = []
  for group in groups:
    cachelocation = (fmeasure, '%s_PARAMID%s_FILEID%s' % (group.name, str(hash(bins_spec)), myutils.checksum(group.file.filename)))
    data = dataman.obtain_data('basic_vessel_radial', 'phi_vessels', group['vessels'], group['tumor'], 10., bins_spec, distancemap_spec, None, cachelocation).copy()
    l.append(data)
  return myutils.MeanValueArray.fromSummation(x.avg for x in l)

def obtain_averaged_mvd(dataman, fmeasure, groups, bins_spec, distancemap_spec):
  l = []
  for group in groups:
    cachelocation = (fmeasure, '%s_PARAMID%s_FILEID%s' % (group.name, str(hash(bins_spec)), myutils.checksum(group.file.filename)))
    data = dataman.obtain_data('basic_vessel_radial', 'mvd', group['vessels'], group['tumor'], 10., bins_spec, distancemap_spec, None, cachelocation)
    l.append(data)
  return myutils.MeanValueArray.fromSummation(x.avg for x in l)

def obtain_average_radial_distance(dataman, fmeasure, groups, bins_spec, distancemap_spec):
  l = []
  for group in groups:
    cachelocation = (fmeasure, '%s_PARAMID%s_FILEID%s' % (group.name, str(hash(bins_spec)), myutils.checksum(group.file.filename)))
    data = dataman.obtain_data('basic_vessel_radial', 'radial_distance', group['vessels'], group['tumor'], 30., bins_spec, distancemap_spec, None, cachelocation)
    l.append(data)
  return myutils.MeanValueArray.fromSummation(x.avg for x in l)


def plot_data(ax, bins, data, std, color, label):
  data = np.ma.filled(data, 0.)
  std  = np.ma.filled(std, 0.)
  xy, clist = [], []
  for i in xrange(len(data)):
    clist.append(0.5*(bins[i]+bins[i+1]))
    xy.append((bins[i], data[i]))
    xy.append((bins[i+1], data[i]))
  xy = np.asarray(xy).transpose()
  ax.plot(xy[0], xy[1], color = color, lw = 1.)
  ax.errorbar(clist, data, yerr=std, lw = 0, color = color, mew=2., elinewidth=1., capsize = 10., label = label)
  for c, d in zip(clist, data):
    ax.text(c, d, '$'+myutils.f2s(d, prec=2, latex=True)+'$', color = color)


def plot_mvd_in_large_bins1(pdfpages, groupslist):
  distancemap_spec = 'levelset'
  for bin_spec in [
      BinsSpecArray([-1.e5, -200, 0., 200., 1.e6]),
      BinsSpecArray([-1.e5, 0., 1.e6])
    ]:
    fig, axes = pyplot.subplots(2,1, figsize = (8,6))
    realbins = {}

    ax = axes[0]
    for i, (groupkey, groups) in enumerate(groupslist):
      data, std = obtain_averaged_phi(dataman, fmeasure, groups, bin_spec, distancemap_spec)
      realbins[groupkey] = ranges = determine_real_bins(data, bin_spec)
      time  = dataman.obtain_data('time', groups[0])
      plot_data(ax, realbins[groupkey], data[:,0], std[:,0], 'rgk'[i], 't = %i h' % time)
    lim = ax.get_ylim()
    ax.vlines(ranges[1:-1], lim[0], lim[1])
    ax.set(ylabel = r'rBV')
    ax.legend()

    ax = axes[1]
    for i, (groupkey, groups) in enumerate(groupslist):
      mvd, std = obtain_averaged_mvd(dataman, fmeasure, groups, bin_spec, distancemap_spec)
      time  = dataman.obtain_data('time', groups[0])
      plot_data(ax, realbins[groupkey], mvd, std, 'rgk'[i], 't = %i h' % time)
    lim = ax.get_ylim()
    ax.vlines(ranges[1:-1], lim[0], lim[1])
    #ax.set(ylim = lim, xlim = (ax.get_xlim()[0], 1000))
    ax.set(xlabel = r'distance from rim [$\mu m$]', ylabel = r'mvd [1/mm$^2$]')
    ax.legend()

    pdfpages.savefig(fig)


def plot_mvd_in_large_bins2(pdfpages, groupslist):
  tumorradi = {}
  for key, groups in groupslist:
    rlist = map(lambda g: dataman.obtain_data('approximate_tumor_radius', g['tumor']), groups)
    tumorradi[key] = np.average(rlist)
  distancemap_spec = 'radial'

  for i, (groupkey, groups) in enumerate(groupslist):
    fig, axes = pyplot.subplots(2,1, figsize = (8,6))

    time  = dataman.obtain_data('time', groups[0])
    rtumor = tumorradi[groupkey]
    for j, bin_spec in enumerate([
        BinsSpecArray([-1.e5, rtumor-200, rtumor, rtumor+200., 1.e6]),
        BinsSpecArray([-1.e5, 200+rtumor, 1.e6])
      ]):
      distances = obtain_average_radial_distance(dataman, fmeasure, groups, bin_spec, distancemap_spec)

      ranges = bin_spec.arange()
      ranges[0] = 0
      ranges[-1] = distances.avg[-1]*0.5

      ax = axes[0]
      phi = obtain_averaged_phi(dataman, fmeasure, groups, bin_spec, distancemap_spec)
      plot_data(ax, ranges, phi.avg, phi.std, 'rg'[j], 't = %i h' % time)
      lim = ax.get_ylim()
      ax.vlines(ranges[1:-1], lim[0], lim[1])
      ax.set(ylabel = r'rBV')
      ax.legend()

      ax = axes[1]
      mvd = obtain_averaged_mvd(dataman, fmeasure, groups, bin_spec, distancemap_spec)
      plot_data(ax, ranges, mvd.avg*1e6, mvd.std*1e6, 'rg'[j], 't = %i h' % time)
      lim = ax.get_ylim()
      ax.vlines(ranges[1:-1], lim[0], lim[1])
      #ax.set(ylim = lim, xlim = (ax.get_xlim()[0], 1000))
      ax.set(xlabel = r'distance from rim [$\mu m$]', ylabel = r'mvd [1/mm$^2$]')
      ax.legend()

    pdfpages.savefig(fig)





if __name__ == '__main__':
  filenames = sys.argv[1:]
  files = [ h5files.open(fn, 'r') for fn in filenames ]
  fmeasure = h5files.open('analyzeBloodVolumeMeasurements.h5', 'a')
  dataman = myutils.DataManager(10, [DataTumorTissueSingle(),
                                      DataVesselRadial(),
                                      DataVesselSamples(),
                                      DataBasicVessel(),
                                      DataDistanceFromCenter()])

  #allgroups = sorted(filter(lambda k: k.startswith('out'), files[0].keys()))
  allgroups = defaultdict(list)
  for f in files:
    keys = filter(lambda k: k.startswith('out'), f.keys())
    for k in keys:
      allgroups[k].append(f[k])
  allgroups = [ (k, allgroups[k]) for k in sorted(allgroups.keys()) ]
  groupslist = [allgroups[-1], allgroups[len(allgroups)/2]]
  outfn = 'bloodVolume-%s.pdf' % splitext(basename(filenames[0]))[0]
  pprint(groupslist)
  print '-> %s' % (outfn)
  with mpl_utils.PdfWriter(outfn) as pdfpages:
    if 0:
      plot_mvd_in_large_bins1(pdfpages, groupslist)
    else:
      plot_mvd_in_large_bins2(pdfpages, groupslist)
