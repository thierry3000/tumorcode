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
  
  
import os,sys
from os.path import join, basename, dirname
import h5py
import numpy as np
import math
import posixpath
from pprint import pprint
import itertools
from scipy.optimize import leastsq
import collections
import glob
import md5
import cPickle
import pprint
import qsub

import extensions # for asarray with h5py support
import krebsutils
#import time as timer

#import vtkcommon
from mystruct import Struct
import myutils

from plotBulkTissue import commonOutputName, contour, imslice, imshow, colorbar
from mpl_utils import PageWriter
from plotIff import LabelFactory, ColorMaps, fig_numbering
import analyzeGeneral

import plotIff


import matplotlib
import matplotlib.pyplot as pyplot
import mpl_utils
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.gridspec as gridspec

LF = LabelFactory
CM = ColorMaps

bins_dist = np.arange(-10000., 10000., 30.)
bin_centers = np.average((bins_dist[1:],bins_dist[:-1]), axis=0)

mastersize = mpl_utils.a4size
f2s = myutils.f2s
gridcolor = (0.7,0.7,0.7)

def text2(ax, txt):
  ax.text(0.01, 0.9, txt, ha = "left", transform = ax.transAxes)



class DataDrugAverages(object):
  keywords = [
    'drug_radial_average', 'drug_vs_vessel_distance_average',
    'drug_vs_vessel_distance_no_tumor_sep_average',
    'drug_global_series_average',
    'radial_average',
    'exposure_histograms_average',
    'global_local_integral_average',
  ]

  def obtain_data(self, dataman, dataname, *args):
    files = args[0]
    args = args[1:]
    assert isinstance(files, list)
    #####
    if dataname in ('drug_radial_average', 'drug_vs_vessel_distance_average', 'drug_vs_vessel_distance_no_tumor_sep_average'):
      group, name = args # time group
      avg_conc = myutils.MeanValueArray.empty()
      for f in files:
        dx, conc = dataman.obtain_data(dataname[:-len('_average')], f, group, name)
        avg_conc += conc.avg
      return dx, avg_conc

    ####
    if dataname == 'drug_global_series_average':
      cons = lambda : myutils.MeanValueArray.empty()
      d = collections.defaultdict(cons)
      for f in files:
        df = dataman.obtain_data('drug_global_series', f)
        for k, v in df.iteritems():
          #d[k] += v.avg
          if d[k]:
            l = min(len(d[k]), len(v.sum))
            d[k] = v[:l]+d[k][:l]
          else:
            d[k] = v
      return d

    if dataname == 'radial_average':
      avg_dat = myutils.MeanValueArray.empty()
      for f in files:
        dx, dat = dataman.obtain_data(('radial', args[0], args[1]), f, *args[2:])
        avg_dat += dat.avg
      return dx, avg_dat

    ###
    if dataname == 'exposure_histograms_average':
      datalist = []
      for f in files:
        cthres, times, histo = dataman.obtain_data(dataname[:-len('_average')], f)
        datalist.append(histo)
      data = np.asarray(datalist)
      avg = np.average(datalist, axis = 0)
      std = np.std(datalist, axis = 0)
      return cthres, times, avg, std

    ####
    if dataname == 'global_local_integral_average':
      # example:
      # data = obtain_data('global_local_integral_average')
      # avg, std = data['auc_in']['max']
      # data[whatever]['avg'] contains the local spatial fluctuations in ``std'', not the std of the average
      names = ['auc_in', 'auc_ex', 'c_max_in', 'c_max_ex']
      tmp = collections.defaultdict(lambda : collections.defaultdict(list))
      for f in files:
        data = dataman.obtain_data('local_integral', f)
        #tc = self.obtain_data(f,'tissue_composition')
        #mask = (tc.phi_viabletumor > 0.01).ravel()
        #mask = dataman.obtain_data('mask_viabletumor',f).ravel()
        t_c = dataman('tissue_composition',f)
        mask = np.asarray(t_c['phi_viabletumor']) > 0.01
        #mask = dataman('mask_viabletumor')
        for name in names:
          #d = np.asarray(data[name]).ravel()[mask]
          d = np.asarray(data[name])[mask]
          tmp[name]['max'].append(d.max())
          tmp[name]['min'].append(d.min())
          tmp[name]['avg'].append(myutils.MeanValueArray.inOneBin(d))
      for name in names:
        for k in ['max', 'min']:
          tmp[name][k] = np.average(tmp[name][k]), np.std(tmp[name][k])
        avg = myutils.MeanValueArray.empty()
        for v in tmp[name]['avg']:
          avg += v
        tmp[name]['avg'] = float(avg.avg), float(avg.std)
      return tmp
    raise KeyError('request of unknown data %s' % ((dataname,)+args))






class DataDrugSingle(object):
  keywords = [
  'peclet_number',
  'conc', 'conc_ex', 'conc_cell', 'src_flow_vess', 'src_diff_vess', 'src_lymph', 'conc_ratio',
  'drug_vessel_correlation', 'drug_vessel_correlation_minus_avg',
  'drug_vs_vessel_distance', 'drug_vs_vessel_distance_no_tumor_sep', 'drug_global', 'drug_global_series', 'time', 'local_integral', 'exposure_histograms', 'drug_radial',
  ]
  for w in ['tumor_distance', 'vessel_distance']:
    for v in ['vessels', 'tissue', 'auc_in', 'c_max_in']:
      keywords.append(('radial', w, v))

  def obtain_data(self, dataman, dataname, *args):
    f, args = args[0], args[1:]
    obtain_data = lambda *args: dataman.obtain_data(args[0], f, *args[1:])
    ld = obtain_data('ld')

    #####
    if dataname == 'peclet_number':
      lengthscale = args[0]
      params = f['parameters/ift'].attrs
      kdiff = float(params['kdiff'])
      vel = obtain_data('iffvelocity_magnitude')
      return vel*(lengthscale/kdiff)

    #####
    if dataname in ['conc', 'conc_ex', 'conc_cell', 'src_flow_vess', 'src_diff_vess', 'src_lymph']:
      group = args[0]
      if len(args)>1 and args[1] == 'imslice':
        return imslice(f[group][dataname])
      else:
        return np.asarray(f[group][dataname])

    #####
    if dataname == 'conc_ratio':
      acx = obtain_data('conc_ex', *args)
      acc = obtain_data('conc_cell', *args)
      return np.ma.masked_invalid(acc / acx)

    #####
    if dataname in ('drug_vessel_correlation', 'drug_vessel_correlation_minus_avg'):
      def read(gmeasure, group):
        corr = myutils.MeanValueArray.read(gmeasure, group)
        r = np.asarray(gmeasure[group]['r'])
        return r, corr

      def write(gmeasure, group):
        phi_vessels = obtain_data('phi_vessels') # v
        conc        = obtain_data('conc', group)  # c
        mask        = obtain_data('mask_tumor') # only consider pairs of v and c  within tumor
        # compute <v*c(r)> vs r
        r, n, c, c2 = krebsutils.radial_correlation(phi_vessels, conc, 10,
                                                    bins_per_unit=2,
                                                    return_raw_data = True,
                                                    subtract_avg = dataname.endswith('minus_avg'),
                                                    mask = mask)
        r *= ld.scale
        corr = myutils.MeanValueArray(n, c, c2)
        corr.write(gmeasure, group)
        gmeasure[group].create_dataset('r', data = r)

      group, = args # time group

      fm = myutils.MeasurementFile(f, h5files)
      return myutils.hdf_data_caching(read, write, fm, ('measurements', dataname, group), (0,4,0))

    ######
    if dataname == 'drug_vs_vessel_distance_no_tumor_sep':
      def read(gmeasure, dsname):
        corr = myutils.MeanValueArray.read(gmeasure, dsname)
        r = np.asarray(gmeasure['r'])
        return r, corr

      def write(gmeasure, dsname):
        tc     = obtain_data('tissue_composition')
        conc   = obtain_data(name, group)
        #do not sort inside an outsite tumor kyle 2014
        #mask   = (tc['dist_tumor'] < 0.).ravel()
        dist   = tc['dist_vessels'].ravel()
        conc   = conc.ravel()
        bins = np.arange(0., 200., 10.)
        
        corr = myutils.MeanValueArray.fromHistogram1d(bins, dist, conc)
        corr.write(gmeasure, dsname)
        r = np.average((bins[1:], bins[:-1]), axis=0)
        if not 'r' in gmeasure:
          gmeasure.create_dataset('r', data = r)

      group, name, = args # time group

      fm = myutils.MeasurementFile(f, h5files)
      return myutils.hdf_data_caching(read, write, fm, ('measurements', dataname, group, name), (0,9,0,0))
    if dataname == 'drug_vs_vessel_distance':
      def read(gmeasure, dsname):
        corr = myutils.MeanValueArray.read(gmeasure, dsname)
        r = np.asarray(gmeasure['r'])
        return r, corr

      def write(gmeasure, dsname):
        tc     = obtain_data('tissue_composition')
        conc   = obtain_data(name, group)
        mask   = (tc['dist_tumor'] < 0.).ravel()
        dist   = tc['dist_vessels'].ravel()[mask]
        conc   = conc.ravel()[mask]

        bins = np.arange(-1000., 1000., 20.)

        corr = myutils.MeanValueArray.fromHistogram1d(bins, dist, conc)
        corr.write(gmeasure, dsname)
        r = np.average((bins[1:], bins[:-1]), axis=0)
        if not 'r' in gmeasure:
          gmeasure.create_dataset('r', data = r)

      group, name, = args # time group

      fm = myutils.MeasurementFile(f, h5files)
      return myutils.hdf_data_caching(read, write, fm, ('measurements', dataname, group, name), (0,9,0,0))

    if dataname[0] == 'radial':
      dataname, name_xdata, name_ydata = dataname
      name_combined = '%s_vs_%s' % (name_ydata, name_xdata)

      def read(gmeasure, group):
        corr = myutils.MeanValueArray.read(gmeasure, group)
        r = np.asarray(gmeasure[group]['r'])
        return r, corr

      def write(gmeasure, group):
        tc = obtain_data('tissue_composition')
        name_to_data = {
          'tumor_distance'  :
            lambda : tc['dist_tumor'].ravel(),
          'vessel_distance' :
            lambda : tc['dist_vessels'].ravel(),
          'vessels' :
            lambda : tc['phi_vessels'].ravel(),
          'tissue' :
            lambda : tc['phi_viabletumor'].ravel(),
          'auc_in' :
            lambda : np.asarray(obtain_data('local_integral')['auc_in']).ravel(),
          'c_max_in' :
            lambda : np.asarray(obtain_data('local_integral')['c_max_in']).ravel(),
        }

        x, y = name_to_data[name_xdata](), name_to_data[name_ydata]()

        if name_xdata == 'vessel_distance':
          mask   = (tc['dist_tumor'] < 0.).ravel()
          x = x[mask]
          y = y[mask]

        bins = np.arange(-1000., 1000., 20.)
        corr = myutils.MeanValueArray.fromHistogram1d(bins, x, y)
        corr.write(gmeasure, group)
        r = np.average((bins[1:], bins[:-1]), axis=0)
        gmeasure[group].create_dataset('r', data = r)

      fm = myutils.MeasurementFile(f, h5files)
      return myutils.hdf_data_caching(read, write, fm, ('measurements', name_combined), (0,8))

    ######
    if dataname == 'drug_radial':
      def read(gmeasure, dsname):
        conc = myutils.MeanValueArray.read(gmeasure, dsname)
        dx = np.asarray(gmeasure['dx'])
        return dx, conc

      def write(gmeasure, dsname):
        dist_tumor  = obtain_data('tissue_composition')['dist_tumor']
        conc        = obtain_data(name, group)
        # compute
        radial_conc = myutils.MeanValueArray.fromHistogram1d(bins_dist, dist_tumor.ravel(), conc.ravel())
        radial_conc.write(gmeasure, dsname)
        if not 'dx' in gmeasure:
          gmeasure.create_dataset('dx', data = bin_centers)

      group, name = args # time group
      fm = myutils.MeasurementFile(f, h5files)
      return myutils.hdf_data_caching(read, write, fm, ('measurements', 'drug_radial', group, name), (0, 3, 0, 0))

    #####
    if dataname == 'drug_global':
      def read(gmeasure, group):
        s = dict( # this shit goes over a group and converts 0-dim datasets to floats
          (k, float(np.asarray(v))) for k, v in gmeasure[group].iteritems()
        )
        return s

      def write(gmeasure, group):
        def doit(g, name, c, **opts):
          r = dict()
          #r = g.create_group(name)
          r['min'] = c.min()
          r['max'] = c.max()
          r['avg'] = np.ma.average(c)
          r['sum'] = c.sum()
          f = opts.get('mask')
          ctum = c[f]
          r['tum_min'] = ctum.min()
          r['tum_max'] = ctum.max()
          r['tum_avg'] = np.ma.average(ctum)
          r['tum_sum'] = ctum.sum()
          # add this stuff to the group
          for k, v in r.iteritems():
            g[name+'_'+k] = v # this creates a 0-dimensional hdf dataset
        # make a flat group with globally averaged ata
        g = gmeasure.create_group(group)
        g['time'] = f[group].attrs['time'] / 3600.
        g['c_in'] = f[group].attrs['DRUG_INFLOW_CONC']
        mask       = obtain_data('mask_tumor')
        for name in ['conc', 'conc_ex', 'conc_cell', 'src_flow_vess', 'conc_ratio' ]:
          c = obtain_data(name, group)
          doit(g, name, c, mask = mask)

      group, = args # time group

      fm = myutils.MeasurementFile(f, h5files)
      return myutils.hdf_data_caching(read, write, fm, ('measurements', 'drug_global', group), (0,4,0))

    ####
    if dataname == 'drug_global_series':
      def read(gmeasure, group):
        return dict(
          (k, np.asarray(v)) for k, v in gmeasure[group].iteritems()
        )

      def write(gmeasure, group):
        groups = myutils.getTimeSortedGroups(f['/'], 'out', key='time')
        l = [ obtain_data('drug_global', posixpath.basename(g.name))
              for g in groups ]
        l = myutils.zipListOfDicts(l)
        g = gmeasure.create_group(group)
        for k, v in l.iteritems():
          g.create_dataset(k, data = v)

      fm = myutils.MeasurementFile(f, h5files)
      l = myutils.hdf_data_caching(read, write, fm, ('measurements', 'drug_global_series'), (0, 1))
      for k, v in l.iteritems():
        l[k] = myutils.MeanValueArray(np.ones_like(v), v, v*v)
      return l

    ####
    if dataname == 'time':
      group, = args
      return f[group].attrs['time']/3600.

    ####
    if dataname == 'local_integral':
      g = f['measurements/drug_local_integral']
      return dict((k,v) for k, v in g.iteritems())

    ####
    if dataname == 'exposure_histograms':
      data = obtain_data('local_integral')
      mask = obtain_data('mask_viabletumor').ravel()
      exposure_data = []
      try:
        for i in range(1000):
          exposure_data.append(data['exposure_time_%02i_in' % i])
      except KeyError:
        pass
      max_time = f['measurements/drug_local_integral'].attrs['tend']
      binsize = 60.*30.
      bins = np.arange(-binsize*0.5, max_time+binsize, 60.*30.)
      t = np.average((bins[1:], bins[:-1]), axis=0)
      cthres = np.asarray([d.attrs['exposure_conc'] for d in exposure_data])
      ret = np.empty((len(cthres), len(t)))
      for i, d in enumerate(exposure_data):
        exptimes = np.asarray(d).ravel()[mask]
        normalization = 1./len(exptimes)
        h, bins = np.histogram(exptimes, bins = bins)
        h = np.asarray(h, dtype = np.float) * normalization
        ret[i][:] = h
      return cthres, t, ret

    raise KeyError('request of unknown data %s' % ((dataname,)+args))






def fit_func_exp(x_in, y_in, initial_params):
  """
    fit f_p(x) = p[0]*exp(-x/p[2]) + p[1]
    returns tuple of fit parameters p and fit function using the parameters:
      (p', f_p'(x))
  """
  def func(p, x): # the fit function
    return p[0]*np.exp(-x/p[2])+p[1]

  def exponential_fit(x, y, initial_params):
    objective_func = lambda p, x, y: (func(p,x)-y)
    p, success = leastsq(objective_func, initial_params, args=(x, y))
    return p, success

  p, success = exponential_fit(x_in, y_in, initial_params)
  if not success:
    raise RuntimeError("fitting failed")

  ret_func = lambda x: func(p, x)
  return p, ret_func


#def gieve_delivery_average(dataman, files, mask_factory, name):
#    mean, width = [], []
#    for f in files:
#      mask = mask_factory(f)
#      dat = dataman('local_integral',f)
#      mask = mask.ravel()
#      c = np.asarray(dat[name]).ravel()
#      c = c[mask]
#      #datalist.append(c)
#      mean.append(np.average(c))
#      width.append(np.std(c))
#    return np.average(mean), np.std(mean), np.average(width), np.std(width)
    #datalist = np.concatenate(datalist)
    #return np.average(datalist), np.std(datalist)



class DataDrugAverages2(object):
  keywords = [
    'drug_delivery_histogram',
    'drug_delivery_avg',
  ]

  def __init__(self, path):
    self.path = path

  def obtain_data(self, dataman, dataname, *args):
    files = args[0]
    args = args[1:]
    assert isinstance(files, list)

    averaged_data_filename = join(self.path,'measuredglobal.h5')
    with h5py.File(averaged_data_filename, 'a') as favg:
      #####
      if dataname in ('drug_delivery_histogram', 'drug_delivery_avg'):
        (mask_factory, regionid), quantity = args

        def write(gmeasure, groupname):
          datalist = []
          means, widths = [], []
          for f in files:
            mask = mask_factory(f)
            dat = dataman('local_integral',f)
            mask = mask.ravel()
            c = np.asarray(dat[quantity]).ravel()
            c = c[mask]
            datalist.append(c)
            means.append(np.average(c))
            widths.append(np.std(c))
          datalist = np.concatenate(datalist)

          mean = np.average(datalist)
          std = np.std(datalist)
          print 'mean %s = %f' % (str(args), mean)

          bins = np.linspace(-10, 10, 400)
          datalist /= mean #make the histogram relative to the mean
          h, bins = np.histogram(datalist, bins)
          h = h/((bins[1:]-bins[:-1]) * len(datalist))
          #print 'propability distribution check:', np.sum(h*(bins[1:]-bins[:-1]))
          g = gmeasure.create_group(groupname)
          g.attrs['files'] = str(','.join(f.filename for f in files))
          g.attrs['regionid'] = regionid
          g.attrs['quantity'] = quantity
          for name, val in zip('mean std avg_mean std_mean avg_width std_width'.split(), [mean,std,np.average(means), np.std(means), np.average(widths), np.std(widths)]):
            g.create_dataset(name, data = val)
          g.create_dataset('bins', data = bins)
          g.create_dataset('h', data = h)

        if dataname == 'drug_delivery_histogram':
          def read(gmeasure, groupname):
            g = gmeasure[groupname]
            bins,h = np.asarray(g['bins']), np.asarray(g['h'])
            x = np.average((bins[1:], bins[:-1]), axis=0)
            c = np.ma.array(h, mask = (h <= 1.e-3))
            return x, c, float(g['mean'][...])
        else:
          def read(gmeasure, groupname):
            g = gmeasure[groupname]
            r = tuple(float(g[name][...]) for name in 'avg_mean std_mean avg_width std_width'.split())
            return r


        identifier = md5.new(cPickle.dumps(tuple(f.filename for f in files)+(regionid, quantity))).hexdigest()
        return myutils.hdf_data_caching(read, write, favg, ('drug_delivery_histogram', identifier), (0,2))



class DataDrugMovie(object):
  keywords = [
    'movieframes', 'movieinfo',
  ]
  def obtain_data(self, dataman, dataname, *args):
    f, args = args[0], args[1:]
    #obtain_data = lambda *args: dataman.obtain_data(args[0], f, *args[1:])
    movie_fn = join(dirname(f.filename), 'moviedata', 'iffdrugmv'+basename(f.filename)[len('iffdrug'):])
    print movie_fn
    with h5py.File(movie_fn, 'r') as movie_f:
      getgroup = lambda i: movie_f['frames']['%04i' % i]
      if dataname == 'movieinfo':
        n = movie_f['frames'].attrs['num_frames']
        times = [ getgroup(i).attrs['time'] for i in xrange(n) ]
        return times
      elif dataname == 'movieframes':
        idx_list, quantity = args[0], args[1]
        return [ np.asarray(getgroup(idx)[quantity]) for idx in idx_list ]



def gieve_histogram(ax, files, dataman, mask_factory, name, (xlabel, ylabel, meanlabel), pdfpages):
    x, c, mean = dataman('drug_delivery_histogram', files, mask_factory, name)
    ax.bar(x, c, width = x[1]-x[0], lw = 0., color = '0.5')
    ax.set(xlabel = xlabel,
           ylabel = ylabel,
           xlim = (0., 3.),
           ylim = (0., 2.5))
    mpl_utils.add_crosshair(ax, (0, None), color = gridcolor)
    ax.grid(linestyle=':', linewidth=0.5, color=gridcolor)
    ax.text(1., min(ax.get_ylim()[1]*0.8, c.max()*1.1), r'%s $=%s$' % (meanlabel, f2s(mean, latex=True)), ha='center') #fontsize=10
    ax.vlines([1.], 0., 2.5, color = 'k')


def labelfactory(quali, region):
  return ('',#LF.math(quali+'/'+LF.avgOver(quali,region)),
          LF.math(ur'%s(%s)' % (LF.probForIn(region),quali)),
          LF.math(LF.avgOver(quali,region)),)

def mask_tumorregion_factory(dataman, min_micron, max_micron):
  def fun(f):
    tc = dataman('tumor_composition', f)
    d = tc['dist_tumor']
    return np.logical_and(d > min_micron, d < max_micron)
  return fun, __name__+cPickle.dumps((min_micron, max_micron))

mask_tumorboundary_factory = lambda dataman: mask_tumorregion_factory(dataman, -200., 30.)
mask_tumorcenter_factory = lambda dataman: mask_tumorregion_factory(dataman, -100000., -500.)
mask_viabletumor_factory = lambda dataman: (lambda f: dataman('mask_viabletumor', f), 'mask_viabletumor')



def plot_averaged_data(files, dataman, pdfpages):
  is_averaging = len(files)>1
  # now time series data
  gg = []
  for f in files:
    groups = myutils.getTimeSortedGroups(f['/'], 'out', key='time')
    groupbasenames = [ posixpath.basename(group.name) for group in groups ]
    gg.append(set(groupbasenames))
  groupbasenames = set.intersection(*gg)
  groups = [ files[0][g] for g in groupbasenames ]
  time = lambda q: dataman.obtain_data('time',files[0], q)

  timepoints = [ 0.1, 1., 2., 4., 8.,24., 72.]
  mygroups = myutils.closest_items(groupbasenames, lambda q: time(q), timepoints)

  if 1 and 'measurements/drug_local_integral' in f: # auc & max histogram

    fig, axes = mpl_utils.subplots_abs_mm((mastersize[0]/mpl_utils.mm_to_inch*0.8,
                                 mastersize[0]/mpl_utils.mm_to_inch*0.55),
                                          3, 2,
                                          10, 15, 65, 28, 15, 8)
    axes = axes.ravel()
    gieve_histogram(axes[0], files, dataman, mask_viabletumor_factory(dataman), 'auc_in', labelfactory(LF.dqualiauc,'V'), pdfpages)
    gieve_histogram(axes[1], files, dataman, mask_viabletumor_factory(dataman), 'c_max_in', labelfactory(LF.dqualimax,'V'), pdfpages)
    gieve_histogram(axes[2], files, dataman, mask_tumorboundary_factory(dataman), 'auc_in', labelfactory(LF.dqualiauc,'TB'), pdfpages)
    gieve_histogram(axes[3], files, dataman, mask_tumorboundary_factory(dataman), 'c_max_in', labelfactory(LF.dqualimax,'TB'), pdfpages)
    gieve_histogram(axes[4], files, dataman, mask_tumorcenter_factory(dataman), 'auc_in', labelfactory(LF.dqualiauc,'TI'), pdfpages)
    gieve_histogram(axes[5], files, dataman, mask_tumorcenter_factory(dataman), 'c_max_in', labelfactory(LF.dqualimax,'TI'), pdfpages)
    for i,ax in enumerate(axes):
      text2(ax, fig_numbering[i])
    pdfpages.savefig(fig, 'qualihisto')

  getdata = lambda g, name: dataman('drug_radial_average', files, g, name) if is_averaging else dataman('drug_radial', files[0], g, name)
  label_radial_conc = LF.math(LF.avgOver('s', r'\theta'))

  def mkfig(nrows, ncols, h = 0.16):
    fig, axes = mpl_utils.subplots_abs_mm((mastersize[0]/mpl_utils.mm_to_inch*0.5*ncols,
                                 mastersize[0]/mpl_utils.mm_to_inch*0.22/0.16*h*nrows),
                                          nrows, ncols,
                                          10, 20, mpl_utils.a4size[0]/mpl_utils.mm_to_inch*0.38, mpl_utils.a4size[0]/mpl_utils.mm_to_inch*h, 15, 5)
    return fig, axes

  def plot(ax, style_num, x, y, **kwargs):
    markers = 'sDo<>d2h'
    #markers=markers+r'$\bowtie$'+r'$\circlearrowleft$',
    colors  = 'krgbmyck'
    #colors = colors+'#eeefff'
    xscale = kwargs.pop('xscale', 1.e-3)
    pp = mpl_utils.errorbar(ax, xscale*x, y,
                marker = markers[style_num], color  = colors[style_num],
                **kwargs)
    return pp


  def plt_drug_radial(ax):
    # drug radial distributions vs tumor distance
    x, y = dataman('radial_average', files, 'tumor_distance', 'vessels') if is_averaging else dataman(('radial', 'tumor_distance', 'vessels'), files[0])
    x = x * 1.e-3
    #ax2 = ax.twinx()
    #ax2.set(yticks = [], ylim = (0., 1.))

    plts = []
    plts += [ mpl_utils.errorbar(ax, x, y.avg, linestyle = '-', label = '$\phi_v$', color = (0.5,0.5,0.5), marker = None) ]

    for style_num, group in enumerate(mygroups):
      x, c = getdata(group, 'conc')
      plts += [ plot(ax, style_num, x, c.avg, yerr = c.std, label = r't = %s' % f2s(time(group)), every = 5) ]

    ax.set(xlim = (lambda (x0, x1): (max(x0, -1.5), min(x1, 1.)))(ax.get_xlim()),
           xlabel = r'$\theta$ [mm]', ylabel = label_radial_conc)
    ax.legend(plts, [p.get_label() for p in plts])


  def plt_src_radial(ax):
    factors = {
      'src_flow_vess' : 1 * 1.e7,
      'src_diff_vess' : 1 * 1.e5,
      'src_lymph'     : 1 * 1.e5,
    }
    group = mygroups[2]
    for style_num, name in enumerate(['src_flow_vess', 'src_diff_vess', 'src_lymph']):
      x, c = getdata(group, name)
      c = c * factors[name]
      p = plot(ax, style_num, x, c.avg, yerr = c.std, label = r'%s $\times %s$' % (name, f2s(factors[name], exponential=True)), every = 5)

    ax.set(xlim = (lambda (x0, x1): (max(x0, -1.5), min(x1, 1.)))(ax.get_xlim()),
           xlabel = r'$\theta$ [mm]', ylabel = label_radial_conc,
           title = 't = %s' % f2s(time(group)))
    ax.legend()


  def plt_drug_vs_vessels(ax): # drug vs distance from vessels
    getdata = lambda g: dataman('drug_vs_vessel_distance_average', files, g, 'conc') if is_averaging else dataman('drug_vs_vessel_distance', files[0], g, 'conc')

    label_radial_conc = LF.math(LF.avgOver('s', ur'\rho'))

#    fig, ax = pyplot.subplots(1,1, figsize = 0.5*np.asarray((mastersize[0], mastersize[0])))
#    fig.subplots_adjust(left = 0.15)
    ax.set(xlabel = ur'$\rho$ [\u03BCm]', ylabel = label_radial_conc)

    x, tum = dataman('radial_average', files, 'vessel_distance', 'tissue') if is_averaging else dataman(('radial', 'vessel_distance', 'tissue'), files[0])
    #ax2 = ax.twinx()
    p = mpl_utils.errorbar(ax, x, tum.avg, linestyle = '-', label = '$\phi_T$', color = (0.5,0.5,0.5), marker = None)
    #ax2.set(ylim=(0.,1.), yticklabels = [])

    for style_num, group in enumerate(mygroups):
      x, c = getdata(group)
      title = r't = %s' % f2s(time(group))

      p = plot(ax, style_num, x, c.avg, yerr = c.std, linestyle = 'none', label = title, xscale = 1.)

      if time(group) < 5.:
        #fit
        mask = np.logical_and(x>=10, x<1000)
        x_in, y_in = x[mask], c.avg[mask]
        (a, b, l), func = fit_func_exp(x_in, y_in, (1., 0., 100.))
        #print 'fit to correlation t = %s params = %s' % (f2s(time(group)), str((a,b,l)))
        #plot
        title = ur'$\propto\,e^{-\rho / %s}$' % myutils.f2s(l, prec=2, latex=True)
        x0, x1 = ax.get_xlim()
        xx = np.arange(0., x1, 10)
        ax.plot(xx, func(xx), color = p[0].get_color(), scalex = False, scaley = False, label = title)
    ax.set(xlim = (-100., 300.))
    ax.legend()
    
  def plt_drug_vs_vessels_timeline(axes): # drug vs distance from vessels
    getdata = lambda g: dataman('drug_vs_vessel_distance_average', files, g, 'conc') if is_averaging else dataman('drug_vs_vessel_distance', files[0], g, 'conc')
    for style_num, group in enumerate(mygroups):
      ax = axes[style_num]
      title = r't = %s' % f2s(time(group))
      #tissue thing
      label_radial_conc = LF.math(LF.avgOver('s', ur'\rho'))
      ax.set(xlabel = ur'$\rho$ [\u03BCm]', ylabel = label_radial_conc)
      x, data_from_manager = dataman('radial_average', files, 'vessel_distance', 'tissue') if is_averaging else dataman(('radial', 'vessel_distance', 'tissue'), files[0])
      p = mpl_utils.errorbar(ax, x, data_from_manager.avg, linestyle = '-', label = '$\phi_T$', color = (0.5,0.5,0.5), marker = None)
      #solutes       
      x, data_from_manager = getdata(group)
      title = r't = %s' % f2s(time(group))
      p = plot(ax, style_num, x, data_from_manager.avg, yerr = data_from_manager.std, linestyle = 'none', label = title, xscale = 1.)
      ax.legend()
      ax.legend(fontsize='x-small')
  def plt_drug_vs_vessels_timeline_no_tumor_sep(axes): # drug vs distance from vessels
    getdata = lambda g: dataman('drug_vs_vessel_distance_no_tumor_sep_average', files, g, 'conc') if is_averaging else dataman('drug_vs_vessel_distance_no_tumor_sep', files[0], g, 'conc')
    for style_num, group in enumerate(mygroups):
      ax = axes[style_num-1]
      title = r't = %s' % f2s(time(group))
      #tissue thing
#      label_radial_conc = LF.math(LF.avgOver('s', ur'\rho'))
#      ax.set(xlabel = ur'$\rho$ [\u03BCm]', ylabel = label_radial_conc)
#      x, data_from_manager = dataman('radial_average', files, 'vessel_distance', 'tissue') if is_averaging else dataman(('radial', 'vessel_distance', 'tissue'), files[0])
#      p = mpl_utils.errorbar(ax, x, data_from_manager.avg, linestyle = '-', label = '$\phi_T$', color = (0.5,0.5,0.5), marker = None)
      #solutes       
      x, data_from_manager = getdata(group)
      title = r't = %s' % f2s(time(group))
      p = plot(ax, style_num, x, data_from_manager.avg, yerr = data_from_manager.std, linestyle = 'none', label = title, xscale = 1.)
      ax.legend()
      ax.legend(fontsize='x-small')


  def plt_quali(ax, quali): # auc, max_conc vs distance from vessels
    label = { 'auc_in' : LF.math(LF.avgOver(LF.dqualiauc, LF.levelsetfunc)),
              'c_max_in' : LF.math(LF.avgOver(LF.dqualimax, LF.levelsetfunc)) }

    x, y = dataman('radial_average', files, 'tumor_distance', 'vessels') if is_averaging else dataman(('radial', 'tumor_distance', 'vessels'), files[0])
    ax2 = ax.twinx()
    mpl_utils.errorbar(ax2, x*1.e-3, y.avg, linestyle = '-', label = '$\phi_v$', color = (0.7,0.7,0.7), marker = None)
    #ax2.legend(loc = mpl_utils.loc.upper_right)
    ax2.set(yticks = [])



    x, y = dataman('radial_average', files, 'tumor_distance', quali) if is_averaging else dataman(('radial', 'tumor_distance', quali), files[0])
    if quali.startswith('auc'): y *= (1./3600.)
    mpl_utils.errorbar(ax, x*1.e-3, y.avg, yerr = y.std, every = 5, color = 'k', label = label[quali])
    ax.set(ylabel = label[quali])
    ax.set(xlim = (-1., 1.))
    ax.set(xlabel = ur'$\theta$ [mm]')

#  def plt_quali(ax): # auc, max_conc vs distance from vessels
#    label = { 'auc_in' : LF.math(LF.avgOver(LF.dqualiauc, LF.levelsetfunc)),
#              'c_max_in' : LF.math(LF.avgOver(LF.dqualimax, LF.levelsetfunc)) }
#
#    x, y = dataman('radial_average', files, 'tumor_distance', 'vessels') if is_averaging else dataman(('radial', 'tumor_distance', 'vessels'), files[0])
#    ax2 = ax.twinx()
#    mpl_utils.errorbar(ax2, x*1.e-3, y.avg, linestyle = '-', label = '$v$', color = (0.7,0.7,0.7), marker = None)
#    ax2.legend(loc = mpl_utils.loc.upper_right)
#    ax2.set(yticks = [])
#
#    x, y = dataman('radial_average', files, 'tumor_distance', 'c_max_in') if is_averaging else dataman(('radial', 'tumor_distance', 'c_max_in'), files[0])
#    mpl_utils.errorbar(ax, x*1.e-3, y.avg, yerr = y.std, every = 5, color = 'k', label = label['c_max_in'])
#    ax.set(ylabel = label['c_max_in'])
#    ax.set(xlim = (-1., 1.))
#    ax.set(xlabel = ur'$\theta$ [mm]')
#
#    ax = ax.twinx()
#    x, y = dataman('radial_average', files, 'tumor_distance', 'auc_in') if is_averaging else dataman(('radial', 'tumor_distance', 'auc_in'), files[0])
#    mpl_utils.errorbar(ax, x*1.e-3, y.avg, yerr = y.std, every = 5, color = 'k', label = label['auc_in'])
#    ax.set(ylabel = label['auc_in'])


  fig, axes = mkfig(1, 2, 0.22)
  axes = axes.ravel()

  for i, ax, func in zip(xrange(len(axes)), axes, [plt_drug_radial, plt_drug_vs_vessels]):#, lambda ax: plt_quali(ax, 'auc_in'), lambda ax: plt_quali(ax, 'c_max_in') ]):
    func(ax)
    ax.grid(linestyle=':', linewidth=0.5, color=gridcolor)
    mpl_utils.add_crosshair(ax, (0, 0), color = gridcolor)
    text2(ax, fig_numbering[i])
    #ax.legend(loc=2)
    ax.legend(fontsize='x-small')
  pdfpages.savefig(fig, 'drugradial')
  
  fig, axes = mkfig(7, 1, 0.22)
  axes = axes.ravel()
  plt_drug_vs_vessels_timeline(axes)
  pdfpages.savefig(fig, 'drug_from_vessel_radial')
  
  fig, axes = mkfig(7, 1, 0.22)
  axes = axes.ravel()
  plt_drug_vs_vessels_timeline_no_tumor_sep(axes)
  pdfpages.savefig(fig, 'drug_from_vessel_radial_no_tumor_sep')

  fig, axes = mkfig(2,1, 0.16)
  axes = axes.ravel()

  plt_quali(axes[0], 'c_max_in')
  plt_quali(axes[1], 'auc_in')
  axes[0].set(xticklabels = [], xlabel = '')
  for i, ax in enumerate(axes.ravel()):
    ax.grid(linestyle=':', linewidth=0.5, color=gridcolor)
    mpl_utils.add_crosshair(ax, (0, 0), color = gridcolor)
    text2(ax, fig_numbering[i])
  pdfpages.savefig(fig, 'qualiradial')


#  mpl_utils.add_crosshair(ax, (0, 0), color = gridcolor)
#  ax.grid(linestyle=':', linewidth=0.5, color=gridcolor)


  if 'measurements/drug_local_integral' in f and 0:
    if is_averaging:
      cthres, times, histo, yerr = dataman('exposure_histograms_average', files)
    else:
      cthres, times, histo = dataman('exposure_histograms', files[0])
      yerr = np.zeros_like(histo)

    histo = np.cumsum(histo[:,::-1], axis=1)[:,::-1]
    times = times/3600.

    fig, ax = pyplot.subplots(1, 1, figsize = (mastersize[0]*0.25, mastersize[0]*0.25))
    for i in range(6):
      mpl_utils.errorbar(ax, times, histo[i], yerr = yerr[i], label = '$c_{exp}$ = $%s$' % f2s(cthres[i]))
      #ax.plot(times, histo[i], label = '$c_{exp}$ = $%s$' % f2s(cthres[i]))
    ax.legend()
    pdfpages.savefig(fig)

    # 3d plot
    fig = pyplot.figure(figsize = (mastersize[0]*0.5, mastersize[0]*0.5))
#        cthres, times, histo = dataman('exposure_histograms')
#        histo = np.cumsum(histo[:,::-1], axis=1)[:,::-1]
    cthres = cthres[:-1]
    histo = histo[:-1,:]

    times = times/3600.

    ax = fig.add_subplot(111,
                           projection='3d',
                           xlabel=r'$c_{exp}$', ylabel=r'$t_{exp}$',
                           title=r'$P$',
                           elev = 30., azim = 45.)

    X, Y = np.meshgrid(cthres, times)
    Z = histo.transpose()

    ax.plot_wireframe(X, Y, Z, color = 'k', rstride = 20)
    #ax.plot_surface(X, Y, Z, cmap = CM.grey, vmin = 0., vmax = 1., linewidth=0)
    #ax.contour(X, Y, Z, colors=['r', 'r', 'r', 'r'], levels = [ 0.2, 0.4, 0.6, 0.8 ])
    #ax.set(xticks = cthres[1:-1:1], xlim = (0., 1.), ylim = (0., times[-1]), zlim = (0., 1.))
    pdfpages.savefig(fig)






def plot_single_values(files, dataman, pdfpages):
    ## ----------   single values -------
    datalists = collections.defaultdict(list)
    for f in files:
      #tc = dataman(f, 'tissue_composition')
      mask = dataman('mask_tumor', f)
      peclet_number = dataman('peclet_number', f, 1.)
      peclet_number = 1./peclet_number
      datalists['tumor_peclet_length'].append(np.average((peclet_number.ravel()[mask.ravel()])))
#      for name in ['conc_cell_tum_min', 'conc_cell_tum_max', 'conc_cell_tum_avg']:
#        gx = dataman(f, 'drug_global_series')
#        v = np.amax(gx[name].avg) # max value over all times
#        datalists[name].append(v)
      gx = dataman('drug_global_series', f)
      mask = gx['time'].avg > 3.
      datalists['avg_ratio'].append(np.average(gx['conc_ratio_tum_avg'].avg[mask]))
    gd = {}
    for k, v in datalists.iteritems():
      gd[k] = (np.average(v), np.std(v))

    intdata = dataman('global_local_integral_average', files)
    for k,v in intdata.iteritems():
      for k2, v2 in v.iteritems():
        v2 = np.asarray(v2)
        if k.startswith('auc_'):
          v2 *= 1./3600.
        gd['%s_%s' % (k, k2)] = v2

    lf2s = lambda q: myutils.f2s(q, latex = True)
    ld = dataman('ld', files[0])
    s  = ld.GetWorldSize()
    s0 = r'system size: %s $\times$ %s $\times$ %s [$mm$]' % (lf2s(s[0]), lf2s(s[1]), lf2s(s[2]))

    namemap = {
      'tumor_peclet_length' : r'$\langle k_d/|v| \rangle$',
#      'conc_cell_tum_min'   : r'$\langle min_\Omega(s_2) \rangle$',
#      'conc_cell_tum_max'   : r'$\langle max_\Omega(s_2) \rangle$',
#      'conc_cell_tum_avg'   : r'$\langle \langle s_2 \rangle_\Omega \rangle$',
      'c_max_in_min'   : r'$\langle min_\Omega max_t (s_2) \rangle$',
      'c_max_in_max'   : r'$\langle max_\Omega max_t (s_2) \rangle$',
      'c_max_in_avg'   : r'$\langle \langle s_2 \rangle_\Omega \rangle$',
      'auc_in_max'          : r'$\langle max_\Omega(auc(s_2)) \rangle$',
      'auc_in_min'          : r'$\langle min_\Omega(auc(s_2)) \rangle$',
      'auc_in_avg'          : r'$\langle \langle (auc(s_2) \rangle_\Omega \rangle$',
      'avg_ratio'           : r'$\langle \langle s_2 / s_1 \rangle_{t>3h} \rangle$',
    }

    lines = [ s0 ]
    for k, v in sorted(namemap.items(), key = lambda (k,v): k):
      x, y = gd[k]
      lines.append(r'%s = $%s \pm %s$ (%s%%)' % (v, lf2s(x), lf2s(y), lf2s(y/x) if abs(x)>1.e-13 else 'inf'))

    fig = pyplot.figure(figsize = (mastersize[0], mastersize[0]))
    fig.text(0.05, 0.5, '\n'.join(lines), fontsize = 10)
    pdfpages.savefig(fig)




def plot_global_average_over_time(files, dataman, pdfpages):
  data_ = dataman('drug_global_series_average', files)
  data = Struct((k,v.avg) for k,v in data_.iteritems())
  data_std = Struct((k,v.std) for k,v in data_.iteritems())

  def plotvstime_err(ax, y, yerr, fmt, **kwargs):
    return mpl_utils.errorbar(ax, data.time, y, yerr = yerr, fmt = fmt,every = 5, **kwargs)

  def plotvstime(ax, y, fmt, **kwargs):
    return plotvstime_err(ax, y, None, fmt, **kwargs)

  ### ---- this is a plot made of two parts with different time scales ----
  if 1:
#    fig, axes = pyplot.subplots(3, 1, figsize = (mastersize[0]*0.5, mastersize[0]), sharex = True)

    gs = gridspec.GridSpec(3,1, height_ratios=[4,1,1])
    gs.update(left=0.15, right=0.9, hspace = 0.05, bottom = 0.1, top = .9)
    fig = pyplot.figure(figsize = (mastersize[0]*0.5, mastersize[0]*0.5))
    ax = fig.add_subplot(gs[0])
    axes = [ax] + [fig.add_subplot(gs[i], sharex = ax) for i in range(1,3)]
    for ax in axes[0:-1]:
      for t in ax.get_xticklabels():
        t.set_visible(False)

    label_smax = LF.math(LF.maxOverTumor('s'))
    label_savg = LF.math(LF.avgOverTumor('s'))
    label_smin = LF.math(LF.minOverTumor('s'))
    label_sin  = ur'$s^v$'
    label_savgcell = LF.math(LF.avgOverTumor('s_2'))
    label_sratio = LF.math(LF.avgOverTumor('s_2/s_1'))


    ax = axes[0]
    #plotvstime(data.c_in, '-o', label=r'$c_{in}$', color = 'r')
    plotvstime_err(ax, data.conc_tum_avg     , data_std.conc_tum_avg     , '-s', label=label_savg       , color = 'k')
    plotvstime_err(ax, data.conc_tum_min     , data_std.conc_tum_min     , '-<', label=label_smin       , color = 'b')
    plotvstime_err(ax, data.conc_tum_max, data_std.conc_tum_max, '->', label=label_smax, color = 'm')
    plotvstime(ax, data.c_in, '-o', label=label_sin, color = 'r')
    ax.set(ylim = (0., 1.))
    ax.legend(loc = mpl_utils.loc.center_left)

#    axins = mpl_utils.inset_axes(ax,
#                       width="30%", # width = 30% of parent_bbox
#                       height="50%", # height : 1 inch
#                       loc=mpl_utils.loc.upper_right)
#    #axins.yaxis.set(visible = False)
#    mpl_utils.errorbar(axins, data.time, data.c_in, fmt = '-o', label=r'$c_{in}$', color = 'r')
#    axins.set(xlim = (0., 6.), ylim = (0., 1.), yticks = [ 0., 1.])

    axins = mpl_utils.inset_axes(ax,
                       width="100%", # width = 30% of parent_bbox
                       height="100%", # height : 1 inch
                       loc=mpl_utils.loc.upper_left,
                       bbox_to_anchor = (0.5, 0.5, 0.5, .5),
                       bbox_transform = ax.transAxes)
    #axins.yaxis.set(visible = False)
    plotvstime_err(axins, data.conc_tum_avg     , data_std.conc_tum_avg     , '-s', label=label_savg       , color = 'k')
    plotvstime_err(axins, data.conc_tum_min     , data_std.conc_tum_min     , '-<', label=label_smin       , color = 'b')
    plotvstime_err(axins, data.conc_tum_max, data_std.conc_tum_max, '->', label=label_smax, color = 'm')
    plotvstime(axins, data.c_in, '-o', label=label_sin, color = 'r')
    axins.set(xlim = (0., 8.))
    l = MaxNLocator(nbins = 4)
    axins.yaxis.set_major_locator(l)
    #axins.grid(linestyle=':', linewidth=0.5, color=gridcolor)

    ax = axes[1]
    plotvstime_err(ax, data.conc_cell_tum_avg, data_std.conc_cell_tum_avg, '-s', label=label_savgcell, color = 'k')
    ax.legend(loc = mpl_utils.loc.center_right)

    ax = axes[2]
    plotvstime_err(ax, data.conc_ratio_avg, data_std.conc_ratio_avg, '-', label=label_sratio, color = 'k')
    ax.legend(loc = mpl_utils.loc.center_right)

    axes[-1].set(xlabel = r't [h]')
    for ax in axes:
      #ax.grid(linestyle=':', linewidth=0.5, color=gridcolor)
      l = MaxNLocator(nbins = 4, prune = 'lower')
      ax.yaxis.set_major_locator(l)

    pdfpages.savefig(fig,'drugvstime')
  ### ---- ----

  else:
    ### ---- plot time dependent global averaged data ----
    fig, axes = pyplot.subplots(5, 1, figsize = 0.5*mastersize, sharex = True)
    axes = axes.ravel()

    ax = axes[0]
    ax.set(title = 'conc in tumor')
    plotvstime(data.c_in, '-o', label=r'$c_{in}$', color = 'r')
    plotvstime_err(data.conc_tum_avg     , data_std.conc_tum_avg     , '-s', label=r'$c_{avg}$'       , color = 'k')
    plotvstime_err(data.conc_tum_min     , data_std.conc_tum_min     , '-<', label=r'$c_{min}$'       , color = 'k')

    ax = axes[1]
    plotvstime_err(data.conc_ex_tum_avg  , data_std.conc_ex_tum_avg  , '-s', label=r'$c_{avg}$ (fluid)', color = 'm')

    ax = axes[2]
    plotvstime_err(data.conc_cell_tum_avg, data_std.conc_cell_tum_avg, '-s', label=r'$c_{avg}$ (cell)', color = 'b')

    ax = axes[3]
    plotvstime_err(data.conc_tum_max, data_std.conc_tum_max, '->', label=r'$c_{max}$', color = 'k')

    ax = axes[4]
    plotvstime_err(data.conc_ratio_avg, data_std.conc_ratio_avg, '-', label=r'$c_{avg}$ (ratio)', color = 'k')

    ax.set(xlabel = r't [h]')  #xscale = 'log', xscale = 'log'
    for ax in axes:
      ax.legend()

    pdfpages.savefig(fig, 'drugvstime')

def plot_drug_with_vessels(f, dataman_, pdfpages):
  groups = myutils.getTimeSortedGroups(f['/'], 'out', key='time')
  groupbasenames = [ posixpath.basename(group.name) for group in groups ]
  dataman = lambda *args: dataman_(args[0], f, *args[1:])
  time = lambda q: dataman('time', q)
  
  timepoints = [ 0.5, 1., 3., 6., 12.,  24., 48., 96. ]
  timepoints = [ 1. , 2. , 3. ]
  timepoints = np.arange(0.5,18.6,2.)
  mygroups = myutils.closest_items(groupbasenames, lambda q: (f[q].attrs['time']/3600.), timepoints)
  mygroups = np.unique(mygroups)

  ld = dataman('ld')
  tc = dataman('tissue_composition')

  def mkfig(nrows = 1, ncols = 1): # sizes set up for image plots
    fig, axes = mpl_utils.subplots_abs_mm((mastersize[0]/mpl_utils.mm_to_inch*0.51*ncols+10,
                                 mastersize[0]/mpl_utils.mm_to_inch*0.43*nrows+10),
                                          nrows, ncols,
                                          10, 10, mpl_utils.a4size[0]/mpl_utils.mm_to_inch*0.4, mpl_utils.a4size[0]/mpl_utils.mm_to_inch*0.4, 15, 4)
    mpl_utils.remove_frame(*axes.ravel())
    if nrows == 1 and ncols == 1: axes = axes[0,0]
    return fig, axes

  def textleft(ax, txt):
    ax.text(0., 1.03, txt, ha = "left", transform = ax.transAxes)
  def textright(ax, txt):
      ax.text(0.5, 1.03, txt, ha = "center", transform = ax.transAxes)

#  imgs = dict(
#    conc = dataman('conc', mygroups[0], 'imslice'),
#    vesselfraction = tc['phi_vessels'],
#  )
#  for k,v in imgs.iteritems():
#    imgs[k] = imslice(v)


  def imshow_(a, cmap, crange, vscale):
    a = imgs[a]
    if crange is None:
      crange = (a.min(), a.max())
    elif crange == 'zero-centered':
      q = np.abs(a).max()
      crange = (-q,q)
    if vscale <> 1.:
      a = a*vscale
    return imshow(ax, a, ld, vmin=vscale*crange[0], vmax=vscale*crange[1], cmap = cmap)

  def contour_(a):
    a = imgs[a]
    return contour(ax, a, ld, levels = [0])
  
  
  #http://stackoverflow.com/questions/10127284/overlay-imshow-plots-in-matplotlib
  if 0:
    # colors
    color1 = matplotlib.colors.colorConverter.to_rgba('white')
    color2 = matplotlib.colors.colorConverter.to_rgba('black')
    #colormap
    cmap1 = matplotlib.colors.LinearSegmentedColormap.from_list('my_cmap',['red','blue'],256)
    cmap2 = matplotlib.colors.LinearSegmentedColormap.from_list('my_cmap2',[color1,color2],256)
  if 1:
    # colors
    color1 = matplotlib.colors.colorConverter.to_rgba('white')
    color1 = (0,0,0,0)
    color2 = matplotlib.colors.colorConverter.to_rgba('black')
    #colormap
    cmap1 = matplotlib.colors.LinearSegmentedColormap.from_list('my_cmap',[color1,'red'],256)
    cmap2 = matplotlib.colors.LinearSegmentedColormap.from_list('my_cmap2',[color1,'blue'],256)
    
  nrows=0
  ncols=0
  if len(timepoints) > 15:
    print('only 15 timepoints are supported!!!!')
  elif len(timepoints) > 10:
    nrows=3
    ncols=5
  elif len(timepoints) > 5:
    nrows=2
    ncols=5
  elif len(timepoints) > 1:
    nrows=1
    ncols=len(timepoints)
  fig, axes = mkfig(nrows,ncols)
  axes = axes.ravel()
  #plt_conc_with_vessels(ax)
  first_min=0
  first_max=0
  for num, group in enumerate(mygroups):
    drug_data = dataman('conc_ex', group, 'imslice')
    ax=axes[num]
    vessel_data = imslice(tc['phi_vessels'])
    #for log plot
    vessel_data=vessel_data+0.05
    #p1 = imshow(ax, drug_data, ld, vmin=0., vmax=None, cmap=cmap1)
    #p2 = imshow(ax, vessel_data, ld, vmin = 0., vmax = 1., cmap=cmap2)
    if num==0:
      first_min=drug_data.min()
      first_max=drug_data.max()
    #p1 = imshow(ax, drug_data, ld, cmap=cmap1, norm=matplotlib.colors.LogNorm(vmin=first_min,vmax=first_max))
    #p2 = imshow(ax, vessel_data, ld, cmap=cmap2, norm=matplotlib.colors.LogNorm(vmin=0.01,vmax=1))
    p1 = imshow(ax, drug_data, ld, cmap=cmap1, vmin=first_min,vmax=first_max)
    p2 = imshow(ax, vessel_data, ld, cmap=cmap2 )    
    contour(ax, imslice(tc['dist_viabletumor']), ld, levels=[0], colors='w')
    textleft(ax, 'time: %s h'%timepoints[num])
    if num==0:
      save_this_handle_for_colorbar = p1
  colorbar(fig, axes[len(mygroups)-1], save_this_handle_for_colorbar)
  pdfpages.savefig(fig, 'vessel_drug_overlay')
  pyplot.close()


def plot_snapshots(f, dataman_, pdfpages):
  groups = myutils.getTimeSortedGroups(f['/'], 'out', key='time')
  groupbasenames = [ posixpath.basename(group.name) for group in groups ]
  dataman = lambda *args: dataman_(args[0], f, *args[1:])
  time = lambda q: dataman('time', q)

  ld = dataman('ld')
  tc = dataman('tissue_composition')

  def mkfig(nrows = 1, ncols = 1): # sizes set up for image plots
    fig, axes = mpl_utils.subplots_abs_mm((mastersize[0]/mpl_utils.mm_to_inch*0.51*ncols+10,
                                 mastersize[0]/mpl_utils.mm_to_inch*0.43*nrows+10),
                                          nrows, ncols,
                                          10, 10, mpl_utils.a4size[0]/mpl_utils.mm_to_inch*0.4, mpl_utils.a4size[0]/mpl_utils.mm_to_inch*0.4, 15, 4)
    mpl_utils.remove_frame(*axes.ravel())
    if nrows == 1 and ncols == 1: axes = axes[0,0]
    return fig, axes

  def textleft(ax, txt):
    ax.text(0., 1.03, txt, ha = "left", transform = ax.transAxes)

  def plt_phitum(ax):
    ax.set(title=ur'$\phi_T$')
    imshow(ax, imslice(tc['phi_viabletumor']), ld, vmin=0., vmax=1., cmap=CM.grey)
    #contour(ax, imslice(tc['dist_tumor']), ld, levels=[0.], colors='r')
    #contour(ax, imslice(tc['dist_necro']), ld, levels=[0.], colors='g')
    contour(ax, imslice(tc['dist_viabletumor']), ld, levels=[0], colors='r')

  def plt_distmap(ax):
    ax.set(title = ur'$\theta \,\, [mm]$')
    # dist_necro is <0 in necrotic regions, dist_tumor is <0 in tumor+necro
    p1 = imshow(ax, 1.e-3*imslice(tc['dist_tumor']), ld, cmap=CM.spectral) #vmin=-50*ld.scale*1.e-3, vmax=50*ld.Scale()*1.e-3
    contour(ax, imslice(tc['dist_tumor']), ld, levels=[0], colors='w')
    colorbar(fig, ax, p1)

  def plt_distmap_vess(ax):
    ax.set(title = ur'$\theta \,\, [mm]$')
    # dist_necro is <0 in necrotic regions, dist_tumor is <0 in tumor+necro
    p1 = imshow(ax, 1.e-3*imslice(tc['dist_vessels']), ld, cmap=CM.spectral) #vmin=-6*ld.scale, vmax=6*ld.Scale()
    contour(ax, imslice(tc['dist_vessels']), ld, levels=[0], colors='w')
    colorbar(fig, ax, p1)

  def plt_peclet(ax):
    peclet_number = dataman('peclet_number', 1.)
    peclet_number = 1./peclet_number
    peclet_number = imslice(peclet_number)
    peclet_number = np.log10(peclet_number)
    p1 = imshow(ax, peclet_number, ld, cmap=CM.spectral)
    ax.set(title = '$log_{10}(k_d/|v|)$')
    contour(ax, imslice(tc['dist_viabletumor']), ld, levels=[0], colors='w')
    colorbar(fig, ax, p1)

  fig, axes = mkfig(3, 1)
  axes = axes.ravel()
  for i, ax in enumerate(axes):
    ax.set_title('(a) (b) (c)'.split()[i]+' '+ax.get_title())
  plt_phitum(axes[0])
  plt_distmap(axes[1])
  plt_distmap_vess(axes[2])
  pdfpages.savefig(fig, 'distmaps')

  fig, ax = mkfig()
  plt_peclet(ax)
  pdfpages.savefig(fig, 'peclet')

  if 'measurements/drug_local_integral' in f:  # this is just to check, probably not for publishing
    # plot the maximal drug conc and the exposure times
    data = dataman('local_integral')

    def plt_cmax_ex(ax):
      ax.set(title = '$ECMAX$')
      p1 = imshow(ax, imslice(data['c_max_ex']), ld, vmin=0., vmax=None, cmap=CM.spectral)
      contour(ax, imslice(tc['dist_viabletumor']), ld, levels=[0], colors='w')
      colorbar(fig, ax, p1)
      textleft(ax, fig_numbering[0])

    def plt_auc_ex(ax):
      ax.set(title = '$ECAUC$')
      p1 = imshow(ax, imslice(data['auc_ex']), ld, vmin=0., vmax=None, cmap=CM.spectral)
      contour(ax, imslice(tc['dist_viabletumor']), ld, levels=[0], colors='w')
      colorbar(fig, ax, p1)
      textleft(ax, fig_numbering[1])

    fig, axes = mkfig(1, 2)
    axes = axes.ravel()
    plt_cmax_ex(axes[0])
    plt_auc_ex(axes[1])
    pdfpages.savefig(fig, 'qualiexmaps')

#    fig, ax = mkfig()
#    ax.set(title = '$max_t(s_2)/max_t(s_1)$')
#    a = imslice(data['c_max_in'])
#    b = imslice(data['c_max_ex'])
#    p1 = imshow(ax, (a/b), ld, vmin=0., vmax=None, cmap=CM.grey)
#    contour(ax, imslice(tc['dist_viabletumor']), ld, levels=[0], colors='w')
#    colorbar(fig, ax, p1)
#    pdfpages.savefig(fig)

    def plt_cmaxin(ax):
      ax.set(title = LF.math(LF.dqualimax))
      p1 = imshow(ax, imslice(data['c_max_in']), ld, vmin=0., vmax=None, cmap=CM.spectral)
      contour(ax, imslice(tc['dist_viabletumor']), ld, levels=[0], colors='w')
      colorbar(fig, ax, p1)
      textleft(ax, fig_numbering[0])

    def plt_aucin(ax):
      ax.set(title = LF.math(LF.dqualiauc))
      p1 = imshow(ax, imslice(data['auc_in'])*(1./3600.), ld, vmin=0., vmax=None, cmap=CM.spectral)
      contour(ax, imslice(tc['dist_viabletumor']), ld, levels=[0], colors='w')
      colorbar(fig, ax, p1)
      #from basic_units import mum
      mpl_utils.add_sizebar(ax, 0.5)
      textleft(ax, fig_numbering[1])

    fig, axes = mkfig(1, 2)
    axes = axes.ravel()
    #axes[1].xaxis.set_units(mum)
    #axes[2].yaxis.set_units(mum)
    plt_cmaxin(axes[0])
    plt_aucin(axes[1])
    pdfpages.savefig(fig, 'qualiinmaps')
    pyplot.close()

  if 'measurements/drug_local_integral' in f and 0:  # this is just to check, probably not for publishing
    fig, axes = pyplot.subplots(2, 3, subplot_kw = dict(aspect='equal', adjustable='box-forced'), figsize = (mastersize[0], mastersize[0]))
    axes = axes.ravel()
    mpl_utils.remove_frame(*axes)

    data = dataman('local_integral')

    for i in range(6):
      name = 'exposure_time_%02i_in' % i
      ax = axes[i]
      ax.set(title = '$c_{exp}$ = %s' % f2s(data[name].attrs['exposure_conc']))
      p1 = imshow(ax, imslice(data[name]), ld, vmin=0., vmax=None, cmap=CM.grey)
      contour(ax, imslice(tc['dist_viabletumor']), ld, levels=[0], colors='w')
      colorbar(fig, ax, p1)

    pdfpages.savefig(fig)
    pyplot.close()

  if 1:
    # snapshot panel
    timepoints = [ 0.5, 1., 3., 6., 12.,  24., 48., 96. ]
    mygroups = myutils.closest_items(groupbasenames, lambda q: (f[q].attrs['time']/3600.), timepoints)
    mygroups = np.unique(mygroups)

    fig, allaxes = pyplot.subplots(2, 4, subplot_kw = dict(aspect='equal', adjustable='box-forced'), figsize = (mastersize[0], mastersize[0]*0.52))
    mpl_utils.subplots_adjust_abs(fig,
                                  wspace = mpl_utils.mm_to_inch*2.,
                                  hspace = mpl_utils.mm_to_inch*2.,
                                  left   = mpl_utils.mm_to_inch*10.,
                                  right  = -mpl_utils.mm_to_inch*20.,
                                  top    = -mpl_utils.mm_to_inch*10.,
                                  bottom = mpl_utils.mm_to_inch*10.
                                  )
    axes = allaxes.ravel()
    mpl_utils.remove_frame(*axes)

    maxval = max(dataman('conc', group, 'imslice').max() for group in mygroups)

    for num, group in enumerate(mygroups):
      ax = axes[num]
      conc = dataman('conc', group, 'imslice')
      p1 = imshow(ax, conc, ld, vmin = 0., vmax = maxval, cmap = CM.spectral) #
      contour(ax, imslice(tc['dist_viabletumor']), ld, levels=[0], colors='w')
      ax.set(xticklabels = [], yticklabels = [])

      at = AnchoredText('%s h' % f2s(dataman('time', group)), loc=2 if num<=3 else 3, frameon=True, borderpad = 0.1)
      at.patch.set_boxstyle("square,pad=0.")
      at.patch.set_linewidth(0)
      ax.add_artist(at)

      if num == 3:
        cax = mpl_utils.inset_axes(ax, width="5%", height="100%", loc=3,
                                   bbox_to_anchor = (1.02, 0, 1, 1), bbox_transform = ax.transAxes, borderpad = 0)
        fig.colorbar(p1, cax = cax)
      if num == 7:
        mpl_utils.add_sizebar(ax)
    pdfpages.savefig(fig, 'drugsnapshots')
    

def plot_snapshots_linear(f, dataman_, pdfpages):
  groups = myutils.getTimeSortedGroups(f['/'], 'out', key='time')
  groupbasenames = [ posixpath.basename(group.name) for group in groups ]
  dataman = lambda *args: dataman_(args[0], f, *args[1:])
  time = lambda q: dataman('time', q)

  ld = dataman('ld')
  tc = dataman('tissue_composition')

  if 1:
    # snapshot panel
    timepoints = [ 0.5, 1., 3., 6., 12.,  24., 48., 96. ]
    timepoints = [ 0.5, 1., 1.5, 2., 2.5,  3., 3.5, 4. ]
    timepoints = np.arange(0.5,6.5,0.5)
    timepoints = np.arange(0.2,3,0.2)
    mygroups = myutils.closest_items(groupbasenames, lambda q: (f[q].attrs['time']/3600.), timepoints)
    mygroups = np.unique(mygroups)

    fig, allaxes = pyplot.subplots(3, 4, subplot_kw = dict(aspect='equal', adjustable='box-forced'), figsize = (mastersize[0], mastersize[0]*0.52))
    mpl_utils.subplots_adjust_abs(fig,
                                  wspace = mpl_utils.mm_to_inch*2.,
                                  hspace = mpl_utils.mm_to_inch*2.,
                                  left   = mpl_utils.mm_to_inch*10.,
                                  right  = -mpl_utils.mm_to_inch*20.,
                                  top    = -mpl_utils.mm_to_inch*10.,
                                  bottom = mpl_utils.mm_to_inch*10.
                                  )
    axes = allaxes.ravel()
    mpl_utils.remove_frame(*axes)

    maxval = max(dataman('conc', group, 'imslice').max() for group in mygroups)

    for num, group in enumerate(mygroups):
      ax = axes[num]
      conc = dataman('conc', group, 'imslice')
      p1 = imshow(ax, conc, ld, vmin = 0., vmax = maxval, cmap = CM.spectral) #
      contour(ax, imslice(tc['dist_viabletumor']), ld, levels=[0], colors='w')
      ax.set(xticklabels = [], yticklabels = [])

      at = mpl_utils.AnchoredText('%s h' % f2s(dataman('time', group)), loc=2 if num<=3 else 3, frameon=True, borderpad = 0.1)
      at.patch.set_boxstyle("square,pad=0.")
      at.patch.set_linewidth(0)
      ax.add_artist(at)

      if num == 3:
        cax = mpl_utils.inset_axes(ax, width="5%", height="100%", loc=3,
                                   bbox_to_anchor = (1.02, 0, 1, 1), bbox_transform = ax.transAxes, borderpad = 0)
        fig.colorbar(p1, cax = cax)
      if num == (len(timepoints)-1):
        mpl_utils.add_sizebar(ax)
    pdfpages.savefig(fig, 'drugsnapshots')
    
    timepoints = np.arange(1,13,1)
    mygroups = myutils.closest_items(groupbasenames, lambda q: (f[q].attrs['time']/3600.), timepoints)
    mygroups = np.unique(mygroups)

    fig, allaxes = pyplot.subplots(3, 4, subplot_kw = dict(aspect='equal', adjustable='box-forced'), figsize = (mastersize[0], mastersize[0]*0.52))
    mpl_utils.subplots_adjust_abs(fig,
                                  wspace = mpl_utils.mm_to_inch*2.,
                                  hspace = mpl_utils.mm_to_inch*2.,
                                  left   = mpl_utils.mm_to_inch*10.,
                                  right  = -mpl_utils.mm_to_inch*20.,
                                  top    = -mpl_utils.mm_to_inch*10.,
                                  bottom = mpl_utils.mm_to_inch*10.
                                  )
    axes = allaxes.ravel()
    mpl_utils.remove_frame(*axes)

    maxval = max(dataman('conc', group, 'imslice').max() for group in mygroups)

    for num, group in enumerate(mygroups):
      ax = axes[num]
      conc = dataman('conc', group, 'imslice')
      p1 = imshow(ax, conc, ld, vmin = 0., vmax = maxval, cmap = CM.spectral) #
      contour(ax, imslice(tc['dist_viabletumor']), ld, levels=[0], colors='w')
      ax.set(xticklabels = [], yticklabels = [])

      at = mpl_utils.AnchoredText('%s h' % f2s(dataman('time', group)), loc=2 if num<=3 else 3, frameon=True, borderpad = 0.1)
      at.patch.set_boxstyle("square,pad=0.")
      at.patch.set_linewidth(0)
      ax.add_artist(at)

      if num == 3:
        cax = mpl_utils.inset_axes(ax, width="5%", height="100%", loc=3,
                                   bbox_to_anchor = (1.02, 0, 1, 1), bbox_transform = ax.transAxes, borderpad = 0)
        fig.colorbar(p1, cax = cax)
      if num == (len(timepoints)-1):
        mpl_utils.add_sizebar(ax)
    pdfpages.savefig(fig, 'drugsnapshots2')
    
    timepoints = [ 0.1, 1., 8, 26 ]
    mygroups = myutils.closest_items(groupbasenames, lambda q: (f[q].attrs['time']/3600.), timepoints)
    mygroups = np.unique(mygroups)

    fig, allaxes = pyplot.subplots(1, 4, subplot_kw = dict(aspect='equal', adjustable='box-forced'), figsize = (mastersize[0], mastersize[0]*0.52))
    mpl_utils.subplots_adjust_abs(fig,
                                  wspace = mpl_utils.mm_to_inch*2.,
                                  hspace = mpl_utils.mm_to_inch*2.,
                                  left   = mpl_utils.mm_to_inch*10.,
                                  right  = -mpl_utils.mm_to_inch*20.,
                                  top    = -mpl_utils.mm_to_inch*10.,
                                  bottom = mpl_utils.mm_to_inch*10.
                                  )
    axes = allaxes.ravel()
    mpl_utils.remove_frame(*axes)

    maxval = max(dataman('conc', group, 'imslice').max() for group in mygroups)

    for num, group in enumerate(mygroups):
      ax = axes[num]
      conc = dataman('conc', group, 'imslice')
      p1 = imshow(ax, conc, ld, vmin = 0., vmax = maxval, cmap = CM.spectral) #
      contour(ax, imslice(tc['dist_viabletumor']), ld, levels=[0], colors='w')
      ax.set(xticklabels = [], yticklabels = [])

      at = mpl_utils.AnchoredText('%s h' % f2s(dataman('time', group)), loc=2 if num<=3 else 3, frameon=True, borderpad = 0.1)
      at.patch.set_boxstyle("square,pad=0.")
      at.patch.set_linewidth(0)
      ax.add_artist(at)

      if num == 3:
        cax = mpl_utils.inset_axes(ax, width="5%", height="100%", loc=3,
                                   bbox_to_anchor = (1.02, 0, 1, 1), bbox_transform = ax.transAxes, borderpad = 0)
        fig.colorbar(p1, cax = cax)
      if num == (len(timepoints)-1):
        mpl_utils.add_sizebar(ax)
    pdfpages.savefig(fig, 'drugsnapshots3')
    
    timepoints = [ 0.1, 8., 24, 72 ]
    mygroups = myutils.closest_items(groupbasenames, lambda q: (f[q].attrs['time']/3600.), timepoints)
    mygroups = np.unique(mygroups)

    fig, allaxes = pyplot.subplots(1, 4, subplot_kw = dict(aspect='equal', adjustable='box-forced'), figsize = (mastersize[0], mastersize[0]*0.52))
    mpl_utils.subplots_adjust_abs(fig,
                                  wspace = mpl_utils.mm_to_inch*2.,
                                  hspace = mpl_utils.mm_to_inch*2.,
                                  left   = mpl_utils.mm_to_inch*10.,
                                  right  = -mpl_utils.mm_to_inch*20.,
                                  top    = -mpl_utils.mm_to_inch*10.,
                                  bottom = mpl_utils.mm_to_inch*10.
                                  )
    axes = allaxes.ravel()
    mpl_utils.remove_frame(*axes)

    maxval = max(dataman('conc', group, 'imslice').max() for group in mygroups)

    for num, group in enumerate(mygroups):
      ax = axes[num]
      conc = dataman('conc', group, 'imslice')
      p1 = imshow(ax, conc, ld, vmin = 0., vmax = maxval, cmap = CM.spectral) #
      contour(ax, imslice(tc['dist_viabletumor']), ld, levels=[0], colors='w')
      ax.set(xticklabels = [], yticklabels = [])

      at = mpl_utils.AnchoredText('%s h' % f2s(dataman('time', group)), loc=2 if num<=3 else 3, frameon=True, borderpad = 0.1)
      at.patch.set_boxstyle("square,pad=0.")
      at.patch.set_linewidth(0)
      ax.add_artist(at)

      if num == 3:
        cax = mpl_utils.inset_axes(ax, width="5%", height="100%", loc=3,
                                   bbox_to_anchor = (1.02, 0, 1, 1), bbox_transform = ax.transAxes, borderpad = 0)
        fig.colorbar(p1, cax = cax)
      if num == (len(timepoints)-1):
        mpl_utils.add_sizebar(ax)
    pdfpages.savefig(fig, 'drugsnapshots4')



def measure_and_plot(filenames):
  #krebsutils.set_num_threads(3)
  rc = matplotlib.rc
  rc('figure', **{'subplot.left' : 0.15,
                  'subplot.right' : 1.-0.15,
                  'subplot.bottom' : 0.2,
                  'subplot.top' : 1.-0.1,
                  'subplot.wspace' : 0.2,
                  'subplot.hspace' : 0.2})

  of = commonOutputName(filenames)
  path = dirname(filenames[0])
  print 'processing: ', filenames
  print 'out:', of
  with PageWriter(of) as pdfpages:
    files = [ h5py.File(fn, 'r+') for fn in filenames ]
    dataman = myutils.DataManager(100, [ analyzeGeneral.DataBasicVessel(),
                                         DataDrugSingle(), 
                                         DataDrugAverages(), 
                                         plotIff.DataTissue(), 
                                         plotIff.DataGlobalIff(), 
                                         DataDrugAverages2(path) ])
    #plot_single_values(files, dataman, pdfpages)
    plot_global_average_over_time(files, dataman, pdfpages)
    plot_averaged_data(files, dataman, pdfpages)
    #plot_snapshots(files[0], dataman, pdfpages)
    plot_snapshots_linear(files[0], dataman, pdfpages)
    plot_drug_with_vessels(files[0], dataman, pdfpages)



def plot_comparison(path):
  rc = matplotlib.rc
  rc('figure', **{'subplot.left' : 0.1,
                  'subplot.right' : 1.-0.1,
                  'subplot.bottom' : 0.2,
                  'subplot.top' : 1.-0.1,
                  'subplot.wspace' : 0.2,
                  'subplot.hspace' : 0.2})

  allfilenames = {
  'default' : 'iffdrug_variant42-inj_*.h5',
  'vX' : 'iffdrug_variant31-inj_*.h5',
  'vY' :'iffdrug_variant42-inf_*.h5',
  'vY2' : 'iffdrug_variant42-linf_*.h5',
  'vA01' : 'iffdrug_variant51-inj_*.h5',
  'vA04' : 'iffdrug_variant52-inj_*.h5',
  #'vA12' : 'iffdrug_variant53-inj_*.h5',
  'vA08'  : 'iffdrug_variantA08-inj_*.h5',
  'vB01' : 'iffdrug_variant61-inj_*.h5',
  #('iffdrug_variant62-inj_*.h5'),
  'vC04' : 'iffdrug_variant71-inj_*.h5',
  'vC07' : 'iffdrug_variant72-inj_*.h5',
  'vZ' : 'iffdrug_variant91-inj_*.h5',
#  'vE02' : 'iffdrug_variantE02-inj_*.h5',
  'vE02' : 'iffdrug_variantE02-inj_hierar_test1-fcc.25x25.80-sample00.h5',

  }
  allfiles = {}
  for k,v in allfilenames.iteritems():
    allfilenames[k] = glob.glob(join(path, v))
    allfiles[k] = list(h5py.File(fn,'r') for fn in allfilenames[k])
    assert allfiles[k] # not empty
  readablename = [
    ('default' , 'b.c.'),
    ('vX' , '(i)'), # 'heavy',
    ('vY' , '(ii-a)'), #'24h inf.',
    ('vY2' , '(ii-b)'), #'cont. inf.',
    ('vZ' , '(iii)'), #'no conv.',
    ('vA01' , '(iv-a)'), #'cap. perm. x10',
    ('vA04' , '(iv-b)'), #'cap. perm x0.1',
    ('vA08' , '(iv-c)'), #'cap. perm x0.01',
    ('vA12' , '(iv-d)'), #'cap. perm x0.001',
    ('vE02' , '(v)'), #'K x10',
    ('vB01' , '(vi)'), #'lymph. x10',
    ('vC04' , '(vii-a)'), #'tum. lymph. x0.1',
    ('vC07' , '(vii-b)'), #'tum. lymph x1',
  ]
  readablename = dict(readablename)

#\item \label{lbl:var_heavier_drug} Heavier drug particles.
#\item \label{lbl:var_infusion} Prolonged infusions.
#\item \label{lbl:var_no_convection} Neglected Convection.
#\item \label{lbl:var_vessel_permeability} Vascular permeability.
#\item \label{lbl:var_conductivity} Hydraudic conductivity of interstitium.
#\item \label{lbl:var_lymphatics} Amount of normal lymphatics.
#\item \label{lbl:var_tumor_lymphatics} Tumor lymphatics.

  dataman = myutils.DataManager(100, [DataDrugSingle(), DataDrugAverages(), plotIff.DataTissue(), plotIff.DataGlobalIff(), plotIff.DataGlobalAverageIff(), DataDrugAverages2(path)])
  averaged_data_filename = join(path,'measuredglobal.h5')
  with h5py.File(averaged_data_filename, 'a') as favg:
    if 1: # generate data
      print 'generating global data cache'
      correlate_data = []
      for idx in 'default vZ vX vE02 vB01 vA01 vA04 vA08 vC04 vC07 vY vY2'.split():
        files = allfiles[idx]
        print 'obtaining data for global correlation from %s' % os.path.commonprefix([f.filename for f in files])
        tmp = collections.defaultdict(list)
        for f in files:
          p = dataman('peclet_number', f, 500.)
          tmp['peclet'].append(np.average(p))
        for k,v in tmp.iteritems():
          tmp[k] = np.asarray((np.average(v), np.std(v)))
        global_iff_data = dataman('iff_global_average', files)
        for k,(avg, std) in global_iff_data.iteritems():
          tmp[k] = np.asarray((avg, std))
        for name in ['auc_in', 'c_max_in']:
          ma, ms, wa, ws = dataman('drug_delivery_avg', files, mask_tumorcenter_factory(dataman), name)
          tmp[name+'_c_avg'] = np.asarray((ma, ms))
          tmp[name+'_c_std'] = np.asarray((wa/ma, ws/ma))
          ma, ms, wa, ws = dataman('drug_delivery_avg', files, mask_tumorboundary_factory(dataman), name)
          tmp[name+'_b_avg'] = np.asarray((ma, ms))
          tmp[name+'_b_std'] = np.asarray((wa/ma, ws/ma))
          ma, ms, wa, ws = dataman('drug_delivery_avg', files, mask_viabletumor_factory(dataman), name)
          tmp[name+'_avg'] = np.asarray((ma, ms))
          tmp[name+'_std'] = np.asarray((wa/ma, ws/ma))
        tmp.update(kdiff = np.asarray((files[0]['parameters/ift'].attrs['kdiff'], 0)))
        tmp.update(sys_name = idx)
        correlate_data.append(tmp)
      correlate_data = myutils.zipListOfDicts(correlate_data)
      # write
      g = favg.recreate_group('global_comparison')
      for k, v in correlate_data.iteritems():
        g.create_dataset(k, data = v)
    else:
      correlate_data = {}
      for k, v in favg['global_comparison'].iteritems():
        correlate_data[k] = np.asarray(v)

  #pprint.pprint(correlate_data)

  d = dict((correlate_data['sys_name'][i], i) for i in range(len(correlate_data['sys_name'])))
  #indices_to_plot = 'default vX vB01 vE02 vA01 vA04 vA08 vC04 vC07 vZ'.split() # everything but the infusions
  indices_to_plot = 'default vX vZ vA01 vA04 vA08 vE02 vB01 vC04 vC07'.split() # everything but the infusions
  indices_to_plot = [ d[k] for k in indices_to_plot ]

  correlate_data['sys_name'] = np.vectorize(lambda x: readablename[x], otypes = [np.object])(correlate_data['sys_name'])
  correlate_data['auc_in_avg'] *= 1./3600.
  correlate_data['auc_in_b_avg'] *= 1./3600.
  correlate_data['auc_in_c_avg'] *= 1./3600.
  correlate_data['tumor_out_per_vol'] *= 1.e3
  correlate_data['tumor_src_plus_per_vol'] *= 1.e3

  labels = dict(
    tumor_src_plus_per_vol = LF.math(LF.avgOverTumor(LF.source+'_l^+')),
    tumor_out_per_vol = LF.math(LF.avgOverTumor(LF.source+'_l')),

    auc_in_avg = LF.math(LF.avgOverViableTumor(LF.dqualiauc)),
    auc_in_std = LF.math(LF.stdOverViableTumor(LF.dqualiauc)),
    c_max_in_avg = LF.math(LF.avgOverViableTumor(LF.dqualimax)),
    c_max_in_std = LF.math(LF.stdOverViableTumor(LF.dqualimax)),

    auc_in_b_avg = LF.math(LF.avgOver(LF.dqualiauc, 'TB')),
    auc_in_b_std = LF.math(LF.stdOver(LF.dqualiauc, 'TB')),
    c_max_in_b_avg = LF.math(LF.avgOver(LF.dqualimax, 'TB')),
    c_max_in_b_std = LF.math(LF.stdOver(LF.dqualimax, 'TB')),

    auc_in_c_avg = LF.math(LF.avgOver(LF.dqualiauc, 'TC')),
    auc_in_c_std = LF.math(LF.stdOver(LF.dqualiauc, 'TC')),
    c_max_in_c_avg = LF.math(LF.avgOver(LF.dqualimax, 'TC')),
    c_max_in_c_std = LF.math(LF.stdOver(LF.dqualimax, 'TC')),

    q_auc_in_avg = LF.math(LF.avgOver(LF.dqualiauc, r'\xi')),
    q_auc_in_std = LF.math(LF.stdOver(LF.dqualiauc, r'\xi')),
    q_c_max_in_avg = LF.math(LF.avgOver(LF.dqualimax, r'\xi')),
    q_c_max_in_std = LF.math(LF.stdOver(LF.dqualimax, r'\xi')),
  )

  print 'plotting ...'

  def plot(ax, idx, idy):
    x, xerr = correlate_data[idx][indices_to_plot,...].transpose()
    y, yerr = correlate_data[idy][indices_to_plot,...].transpose()
    mpl_utils.errorbar(ax, x, y, xerr=xerr, yerr=yerr, lw = 0, marker = 'x', markersize = 10)
    ax.set(xlim = (-0.05 * x.max(), 1.05*x.max()), ylim = (-0.05*y.max(), 1.05*y.max()))
    ax.set(xlabel = labels[idx], ylabel = labels[idy])
    ax.grid(linestyle=':', linewidth=0.5, color=gridcolor)
    #ax.set(title = '%s vs. %s' % (idy, idx))

  def plot_by_name(ax, idy, moment):
    names = correlate_data['sys_name'][indices_to_plot,...]
    index_of_bc = 0
    x = np.arange(len(indices_to_plot))
    x = x.max() - x
    lim = (x.min()-0.6, x.max()+0.6)

    ax.vlines(correlate_data[idy+'_'+moment][index_of_bc][0], lim[0], lim[1], colors = 'k', zorder = 1)
    ax.vlines(correlate_data[idy+'_b_'+moment][index_of_bc][0], lim[0], lim[1], colors = 'm', zorder = 1)
    ax.vlines(correlate_data[idy+'_c_'+moment][index_of_bc][0], lim[0], lim[1], colors = 'b', zorder = 1)

    y, yerr = correlate_data[idy+'_'+moment][indices_to_plot,...].transpose()
    ax.barh(x, y, xerr=yerr, height = 0.3, align='center', color = 'k', label = r'$\xi = V$', ecolor = 'k')
    y, yerr = correlate_data[idy+'_b_'+moment][indices_to_plot,...].transpose()
    ax.barh(x-0.3, y, xerr=yerr, height = 0.3, align='center', color = 'm', label = r'$\xi = TB$', ecolor = 'm')
    y, yerr = correlate_data[idy+'_c_'+moment][indices_to_plot,...].transpose()
    ax.barh(x+0.3, y, xerr=yerr, height = 0.3, align='center', color = 'b', label = r'$\xi = TI$', ecolor = 'b')

    ax.set(yticklabels = names, yticks = x,
           xlabel = labels['q_'+idy+'_'+moment],
           ylim = lim)
    ax.grid(linestyle=':', linewidth=0.5, color=gridcolor)


  def layoutedfigure():
    size = mpl_utils.a4size[0]*np.asarray((0.5, 0.43))
    relmm = lambda s, x: (x * mpl_utils.mm_to_inch)/s
    width = mpl_utils.a4size[0] / size[0] * 0.15
    height = mpl_utils.a4size[0] / size[1] * 0.15
    left = relmm(size[0], 15)
    bottom = relmm(size[1], 10)
    fig = pyplot.figure(figsize = size)
    axes = [
      fig.add_axes([left, 0.5+bottom, width, height]),
      fig.add_axes([0.5+left, 0.5+bottom, width, height]),
      fig.add_axes([left, bottom, width, height]),
      fig.add_axes([0.5+left, bottom, width, height]),
    ]
    return fig, axes

  with PageWriter('comparison') as pdfpages:

#    fig, axes = pyplot.subplots(2, 2, figsize = (a4size[0], a4size[0]))
#    axes = axes.ravel()

    print 'correlation with flow ...'

    for srcterm in ['tumor_src_plus_per_vol', 'tumor_out_per_vol']:
      for region in ['','b_', 'c_']:
        fig, axes = layoutedfigure()
        plot(axes[0], srcterm, 'auc_in_%savg' % region)
        plot(axes[1], srcterm, 'auc_in_%sstd' % region)
        plot(axes[2], srcterm, 'c_max_in_%savg' % region)
        plot(axes[3], srcterm, 'c_max_in_%sstd' % region)
        for i, ax in enumerate(axes):
          ax.text(0.95, 0.9, fig_numbering[i], transform = ax.transAxes, ha = 'right')
        pdfpages.savefig(fig)

    print 'barcharts ...'

    if 1:
      fig, axes = pyplot.subplots(2, 2, figsize = (mastersize[0]*0.8, mastersize[0]*0.7))
      mpl_utils.subplots_adjust_abs(fig, top = -5 * mpl_utils.mm_to_inch, bottom = 10 * mpl_utils.mm_to_inch,
                                         left = 15 * mpl_utils.mm_to_inch, right = -5 * mpl_utils.mm_to_inch)
      axes = axes.ravel()
      plot_by_name(axes[0], 'auc_in','avg')
      plot_by_name(axes[1], 'auc_in','std')
      plot_by_name(axes[2], 'c_max_in','avg')
      plot_by_name(axes[3], 'c_max_in','std')
      axes[0].legend(loc = mpl_utils.loc.upper_right)
      for i, ax in enumerate(axes):
        ax.text(0.00, 1.02, fig_numbering[i], transform = ax.transAxes, ha = 'left')
      pdfpages.savefig(fig)

    print 'histograms ...'

    if 1:
      order = 'default vX vY vY2 vB01 vE02 vA01 vA04 vA08 vC04 vC07 vZ'.split()
      def histogramtable(figtitle, mask, quantity, labels):
        fig, axes = pyplot.subplots(7, 2, figsize = (mastersize[0]*0.8, mastersize[0]*0.75))
        mpl_utils.subplots_adjust_abs(fig, top = -5 * mpl_utils.mm_to_inch, bottom = 10 * mpl_utils.mm_to_inch,
                                           left = 15 * mpl_utils.mm_to_inch, right = -5 * mpl_utils.mm_to_inch)

        axes = axes.ravel()
        for ax in axes:
          ax.set(visible = False)

        items = [ (k,allfiles[k]) for k in order ]
        for confignum, (name, files) in enumerate(items):
          ax = axes[confignum]
          ax.set(visible = True)
          gieve_histogram(ax, files, dataman, mask, quantity, labels, pdfpages)
          ax.text(0.95, 0.8, readablename[name], transform = axes[confignum].transAxes, ha = 'right')
          if confignum <> len(items)-1:
            ax.set(xticklabels = [], yticklabels = [], xlabel = '', ylabel = '')
        #fig.subplots_adjust(bottom = 0.1, left = 0.15, right = 0.9, wspace = 0.25, hspace = 0.3)
        pdfpages.savefig(fig)

      histogramtable('region: viable tumor, quantity: auc',  mask_viabletumor_factory(dataman), 'auc_in', labelfactory(LF.dqualiauc,'V') )
      histogramtable('region: viable tumor, quantity: cmax', mask_viabletumor_factory(dataman), 'c_max_in', labelfactory(LF.dqualimax,'V') )
      histogramtable('region: tumor rim, quantity: auc', mask_tumorboundary_factory(dataman), 'auc_in', labelfactory(LF.dqualiauc,'TB'))
      histogramtable('region: tumor rim, quantity: cmax', mask_tumorboundary_factory(dataman), 'c_max_in', labelfactory(LF.dqualimax,'TB'))
      histogramtable('region: tumor interior, quantity: auc', mask_tumorcenter_factory(dataman), 'auc_in', labelfactory(LF.dqualiauc,'TI'))
      histogramtable('region: tumor interior, quantity: cmax', mask_tumorcenter_factory(dataman), 'c_max_in', labelfactory(LF.dqualimax,'TI'))


if __name__ == "__main__":
  import argparse
  parser = argparse.ArgumentParser(description='Analyze Drug distributions.')  
  parser.add_argument('iffFileNames', nargs='+', type=argparse.FileType('r'), default=sys.stdin, help='iff files to calculate')   
  goodArguments, otherArguments = parser.parse_known_args()
  qsub.parse_args(otherArguments)
  
  #create filename due to former standards
  filenames=[]
  for fn in goodArguments.iffFileNames:
    filenames.append(fn.name)
  try:
    for fn in filenames:
      if not os.path.isfile(fn):
        raise AssertionError('The file %s is not present!'%fn)
  except Exception, e:
    print e.message
    sys.exit(-1)
  measure_and_plot(filenames)
