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
import h5files
import numpy as np
import extensions # for hdf5 support in np.asarray
import krebsutils
import myutils
import posixpath
import math
import collections
import itertools

import matplotlib
import matplotlib.cm
import matplotlib.pyplot as pyplot
import mpl_utils

from quantities import Prettyfier

from krebs.analyzeGeneral import DataVesselGlobal, DataTumorTissueSingle, DataDistanceFromCenter, DataBasicVessel, DataVesselSamples, DataVesselRadial, BinsSpecRange, BinsSpecArray, obtain_distmap_, generate_samples, combineSamples, HdfCacheRadialDistribution, CalcPhiVessels, calc_distmap
from krebs.analyzeBloodFlow import DataTumorBloodFlow


from myutils import f2l


class EnsembleItem(object):
  def __init__(self, **kwargs):
    self.path = ''
    self.time = None
    self.gvessels = None
    self.gtumor = None
    self.group = None
    for k, v in kwargs.iteritems():
      setattr(self, k, v)


class EnsembleFiles(object):
  def __init__(self, dataman, filenames, pattern):
    files     = [h5files.open(fn, 'r+') for fn in filenames]
    items     = []
    has_tumor = True
    for f in files:
      paths = myutils.walkh5(f['.'], pattern)
      for path in paths:
        g = f[path]
        if g.attrs.get('CLASS', None) == 'GRAPH':
          gvessels = g
          gtumor   = None
          try:
            source = h5files.openLink(g, 'SOURCE')
            gtumor = source.parent['tumor']
            g = source.parent
          except Exception, e:
            raise RuntimeError('tried to get tumor data but failed:' + str(e))
        else:
          gvessels, gtumor = g['vessels'], (g['tumor'] if 'tumor' in g else None)
        e = EnsembleItem(path = path, gvessels = gvessels, gtumor = gtumor, group = g)
        e.time = g.attrs['time']
        has_tumor = has_tumor and gtumor is not None
        e.vessel_system_length = dataman.obtain_data('vessel_system_length', gvessels)
        items.append(e)
    if has_tumor:
      d = collections.defaultdict(list) # path -> list of EnsembleItem
      for e in items:
        d[e.path].append(e)
      tumor_snapshot_times = dict((k,np.average(map(lambda e: e.time, v))) for k,v in d.items())
      tumor_snapshot_order = sorted(tumor_snapshot_times.keys(), key = (lambda path: tumor_snapshot_times[path]))
      tumor_snapshots      = [(d[path], path, tumor_snapshot_times[path]) for path in tumor_snapshot_order]
    self.files = files
    self.items = items
    self.tumor_snapshots = tumor_snapshots # list of tuple(items, path, time)
    self.has_tumor = has_tumor
    ld = krebsutils.read_lattice_data_from_hdf(items[0].gvessels['lattice'])
    self.world_size = ld.GetWorldSize()



def GetAverageApproximateTumorRadius(dataman, ensemble, ensembleitems):
  l = map(lambda item: dataman.obtain_data('approximate_tumor_radius', item.gtumor),
          ensembleitems)
  return np.average(l)


def CollectAllRadialData(dataman, ensembleitems, measurementinfo):
  bin_spec = BinsSpecRange(-100000., 100000., 100.) #MakeBinSpec(distancemap_spec, item.gtumor, only_tumor_vs_normal)
  curves = collections.defaultdict(list)
  for item in ensembleitems:
    for name in ['mvd','velocity','phi_vessels', 'shearforce', 'radius', 'hematocrit']:
      data = dataman.obtain_data('basic_vessel_radial', name, item.gvessels, item.gtumor, measurementinfo['sample_length'], bin_spec, measurementinfo['distancemap_spec'], None, measurementinfo['cachelocation_callback'](item.group))
      curves[name].append(data)
  print ' ... finished getting radial curves'
  return bin_spec, curves



def PlotRadialCurves(pdfwriter, bins_spec, snapshotlist, measurementinfo, world_size):
  charsize = 12/90.
  figscale = mpl_utils.a4size[0]*0.33
  fig, axes = pyplot.subplots(2, 2, dpi = 90, figsize=(figscale*2,figscale*1.5))
  mpl_utils.subplots_adjust_abs(fig, left=4*charsize, right=-2*charsize, top=-3.*charsize, bottom=5.*charsize, hspace=8.*charsize, wspace=12*charsize)
  axes = axes.ravel()

  bins = bins_spec.arange()*1.e-3

  FmtTime = lambda time: 't='+f2l(round(time))

  default_colors = 'rgbcmyk'

  def plot(ax, name, scalefactor = 1., label = None, colors = default_colors, errorbars = True, zero_ylim = True):
    for color, i, (time, tumor_radius, curves) in itertools.izip(itertools.cycle(colors), itertools.count(), snapshotlist):
      curve = myutils.MeanValueArray.fromSummation(map(lambda x: x.avg, curves[name]))
      label = FmtTime(time)
      mask = ~curve.avg.mask
      if errorbars:
        mpl_utils.errorbar(ax, bins[mask], scalefactor*curve.avg[mask], yerr=scalefactor*curve.std_mean[mask], label = label,
                           marker = None, color = color, every = 2)
      else:
        ax.plot(bins[mask], scalefactor*curve.avg[mask], label = label, color=colors[i])
    if zero_ylim:
      _, ymax = ax.get_ylim()
      ax.set_ylim((0., ymax))
    if measurementinfo['distancemap_spec']=='levelset':
      ax.set(xlim=(-2.,2.))
      mpl_utils.add_crosshair(ax, (0.,None), ls=':')
    else:
      mpl_utils.add_crosshair(ax, (0.5e-3*world_size[0], None))
      ax.set(xlim=(0, 0.5e-3*world_size[0]))
      for color, i, (time, tumor_radius, curves) in itertools.izip(itertools.cycle(colors), itertools.count(), snapshotlist):
        mpl_utils.add_crosshair(ax, (tumor_radius*1.e-3, None), ls=':', color = color)
      

  def fmt_axis_labels(ax, name, scalefactor = None, hidex = True):
    ax.set(ylabel = (r'$%s$ %s' % (Prettyfier.get_sym(name), Prettyfier.get_bunit(name))) + ((r' $\times\, %s$' % f2l(scalefactor, exponential=True)) if scalefactor else ''),
           title  = Prettyfier.get_label(name))
    if hidex: ax.set(xticklabels = [])
    else:
      xsym = r'\phi' if measurementinfo['distancemap_spec']=='levelset' else r'|x|'
      ax.set(xlabel = r'$%s$ [$mm$]' % xsym)

#  if measurementinfo.distancemap_spec == 'radial':
#    axtwin = axes[0].twinx()
#    plot(axtwin, 'theta_tumor', colors = [(0.5,0.5,0.5)]*7, errorbars=False)
  plot(axes[0], 'mvd', label = True)
  fmt_axis_labels(axes[0], 'mvd')
  theLegend = axes[0].legend()
  theLegend.get_title().set_fontsize('2')
  pyplot.setp(axes[0].get_legend().get_texts(),fontsize='2')
  theLegend.__sizes = [1]

  plot(axes[1], 'phi_vessels')
  fmt_axis_labels(axes[1], 'phi_vessels')

  plot(axes[2], 'radius')
  fmt_axis_labels(axes[2], 'radius', hidex = False)

  plot(axes[3], 'velocity', scalefactor = 1.0e-3)
  fmt_axis_labels(axes[3], 'velocity_mm', hidex = False)

  pdfwriter.savefig(fig, postfix='_o2radial')




def doit(filenames, pattern):
  dataman = myutils.DataManager(20, [ DataTumorTissueSingle(), DataDistanceFromCenter(), DataBasicVessel(), DataVesselSamples(), DataVesselRadial(), DataVesselGlobal(), DataTumorBloodFlow()])

  ensemble = EnsembleFiles(dataman, filenames, pattern)
  if ensemble.has_tumor:
    print 'paths: ', map(lambda (_t0, path, _t1): path, ensemble.tumor_snapshots)
  else:
    print 'paths: ', set(map(lambda e: e.path, ensemble.items))

  prefix, suffix = myutils.splitcommonpresuffix(map(lambda s: basename(s), filenames))
  outputbasename, _ = splitext(prefix+suffix)
  fn_measure = outputbasename+'-radial-cache.h5'
  f_measure = h5files.open(fn_measure, 'a')
  def cachelocation(g):
    path = posixpath.join('FileCS_'+myutils.checksum(basename(g.file.filename)), g.name.strip(posixpath.sep))
    return (f_measure, path)

  #name = ensemble.items[0].path.split('/')[-1]
  measurementinfo = dict(sample_length = 30.,
                         cachelocation_callback = cachelocation,
                         distancemap_spec = 'radial')

  print 'getting radial curves'
  stuff = []
  stuff2 = []
  for items, path, time in ensemble.tumor_snapshots:
    bins_spec, curves = CollectAllRadialData(dataman, items, measurementinfo)
    tumorradius = GetAverageApproximateTumorRadius(dataman, ensemble, items)
    stuff.append((time, tumorradius, curves))
    stuff2.append((time, tumorradius, curves, path))

  #output_filename+= '_'+name
  with mpl_utils.PdfWriter(outputbasename+'-radial.pdf') as pdfwriter:
    PlotRadialCurves(pdfwriter, bins_spec, stuff, measurementinfo, ensemble.world_size)
  
  with h5py.File(outputbasename+'-radial.h5', 'w') as f:
    f.attrs['COMMONPREFIX'] = os.path.commonprefix(map(lambda s: basename(s), filenames))
    f.create_dataset('files', data = filenames)
    f.create_dataset('groups', data = np.asarray(set(map(lambda item: item.path, ensemble.items)), dtype = np.str))
    for i, ( time, tumorradius, curves, path) in enumerate(stuff2):
      g = f.create_group(path.strip('/').replace('/','-'))
      g.attrs['TUMORRADIUS'] = tumorradius
      g.attrs['TIME'] = time
      for name, curve in curves.iteritems():
        curve = myutils.MeanValueArray.fromSummation(map(lambda x: x.avg, curve))
        g.create_dataset(name+'/avg', data = np.asarray(curve.avg))
        g.create_dataset(name+'/std', data = np.asarray(curve.std))
        g.create_dataset(name+'/mask', data = ~curve.avg.mask)
        g.create_dataset(name+'/std_mean', data = np.asarray(curve.std_mean))
      g.create_dataset('bins', data = bins_spec.arange())

if __name__ == '__main__':
  krebsutils.set_num_threads(2)
  filenames, pattern = sys.argv[1:-1], sys.argv[-1]
  doit(filenames, pattern)