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

''' yet this a copy of plotsForPaper.py of the 
    detailedo2Analysis
    maybe one could turn this into comparison of
    yes Adaption vs. no Adaption ???
'''
if __name__ == '__main__':
  import os.path, sys
  sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../..'))
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
from copy import copy

import matplotlib
import matplotlib.cm
import matplotlib.pyplot as pyplot
import mpl_utils

import scipy.optimize

from krebs.quantities import Prettyfier
from krebs import quantities
#DataVesselGlobal, DataTumorTissueSingle, DataDistanceFromCenter, DataBasicVessel, DataVesselSamples, DataVesselRadial, BinsSpecRange, BinsSpecArray, obtain_distmap_, generate_samples, combineSamples, HdfCacheRadialDistribution, CalcPhiVessels, calc_distmap
from krebs.analyzeGeneral import BinsSpecRange, BinsSpecArray
from krebs import analyzeGeneral
from krebs import analyzeBloodFlow
from krebs import detailedo2Analysis
from krebs import detailedo2

from myutils import f2l, f2s

#vesselTypeColors = [
#   '#e50000', #red
#   '#9a0eea', #violet
#   '#0343df', #blue
#   '#f97306', #orange
#   '#677a04', #olive green
#   '#ceb301', #mustard
#   '#04d8b2', #aquamarine
#   '#06470c', #forest greeen
#   '#840000', #dark red
#   '#607c8e', #blue grey
#]
#alternateVesselTypeColors = 'kkkkkkkkkkkkkkkkkkkkk'
vesselTypeColors = 'kkkkkkkkkkkkkkkkkkk' 
alternateVesselTypeColors = ['0.5',]*10
vesselTypeMarkers = '<>*osd^+x'


RewriteVesselLabel = lambda s: 'RC%i' % (ord(s)-ord('A')+1)


def PlotRadialCurves(pdfwriter, bins_spec, snapshotlist, measurementinfo):
  charsize = 12/90.
  figscale = 1.8
  fig, axes = pyplot.subplots(4, 2, dpi = 90, figsize=(5.*figscale,5.*figscale))
  mpl_utils.subplots_adjust_abs(fig, left=6*charsize, right=-2*charsize, top=-3.*charsize, bottom=8.*charsize, hspace=8.*charsize, wspace=12*charsize)
  axes = axes.ravel()

  bins = bins_spec.arange()
  bins_x = 0.5*(bins[:-1]+bins[1:])

  FmtTime = lambda time: 't='+f2l(round(time))

  default_colors = 'rgbcmyk'

  def plot(ax, name, scalefactor = 1., label = None, colors = default_colors, errorbars = True, zero_ylim = True):
    for i, (time, tumor_radius, curves) in enumerate(snapshotlist):
      curve = myutils.MeanValueArray.fromSummation(map(lambda x: x.avg, curves[name]))
      label = FmtTime(time)
      mask = ~curve.avg.mask #>0
      #print name, curve.avg[mask]
      if errorbars:
        mpl_utils.errorbar(ax, bins[mask], scalefactor*curve.avg[mask], yerr=scalefactor*curve.std_mean[mask], label = label,
                           marker = None, color = colors[i], every = 2)
      else:
        ax.plot(bins[mask], scalefactor*curve.avg[mask], label = label, color=colors[i])
    if zero_ylim:
      _, ymax = ax.get_ylim()
      ax.set_ylim((0., ymax))
    if measurementinfo.distancemap_spec=='levelset':
      ax.set(xlim=(-2.,2.))
      mpl_utils.add_crosshair(ax, (0.,None), ls=':')
    else:
      #ax.set(xlim=(0.,4.))
      for i, (time, tumor_radius, curves) in enumerate(snapshotlist):
        mpl_utils.add_crosshair(ax, (tumor_radius*1.e-3, None), ls=':')

  def fmt_axis_labels(ax, name, scalefactor = None, hidex = True):
    ax.set(ylabel = (r'$%s$ %s' % (Prettyfier.get_sym(name), Prettyfier.get_bunit(name))) + ((r' $\times\, %s$' % f2l(scalefactor, exponential=True)) if scalefactor else ''),
           title  = Prettyfier.get_label(name))
    if hidex: ax.set(xticklabels = [])
    else:
      xsym = r'\phi' if measurementinfo.distancemap_spec=='levelset' else r'|x|'
      ax.set(xlabel = r'$%s$ [$mm$]' % xsym)

#  if measurementinfo.distancemap_spec == 'radial':
#    axtwin = axes[0].twinx()
#    plot(axtwin, 'theta_tumor', colors = [(0.5,0.5,0.5)]*7, errorbars=False)
  plot(axes[0], 'mvd', label = True)
  fmt_axis_labels(axes[0], 'mvd')
  axes[0].legend()

#  plot(axes[1], 'velocity', scalefactor = 1.e-3)
#  fmt_axis_labels(axes[1], 'velocity', scalefactor = 1.e3,)
  plot(axes[1], 'phi_vessels')
  fmt_axis_labels(axes[1], 'phi_vessels')

  plot(axes[2], 'po2')
  fmt_axis_labels(axes[2], 'po2',)

  plot(axes[3], 'po2_tissue')
  fmt_axis_labels(axes[3], 'po2_tissue',)

  plot(axes[4], 'sat_via_hb_ratio')
  fmt_axis_labels(axes[4], 'sat_via_hb_ratio')

  plot(axes[5], 'dS_dx', zero_ylim = False, scalefactor = -1.)
  fmt_axis_labels(axes[5], 'dS_dx')

  plot(axes[6], 'vfhb_oxy')
  fmt_axis_labels(axes[6], 'vfhb_oxy', hidex = False)

  plot(axes[7], 'vfhb_deoxy')
  fmt_axis_labels(axes[7], 'vfhb_deoxy', hidex = False)

  pdfwriter.savefig(fig, postfix='_o2radial')


def FitExpFunction(x_in, y_in, initial_params):
  from scipy.optimize import leastsq
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


#def PlotOxygenVsVessels(pdfwriter, (curves_by_bins, phi_vessels_by_bins), binspec):
#  charsize = 12/90.
#  figscale = 1.8
#  fig, ax = pyplot.subplots(1, 1, dpi = 90, figsize=(2.5*figscale,3.*figscale))
#  mpl_utils.subplots_adjust_abs(fig, left=6*charsize, right=-2*charsize, top=-3.*charsize, bottom=8.*charsize, hspace=8.*charsize, wspace=12*charsize)
#  ax2 = ax.twinx()
#
#  bins = binspec.arange()
#  bins = np.average((bins[1:],bins[:-1]), axis=0)
#  colors = 'rgbcmyk'
#  labels = ['tumor', 'boundary', 'normal']
#  for i in curves_by_bins.keys():
#    curve = curves_by_bins[i]
#    phi_vessels = phi_vessels_by_bins[i]
#    label = labels[i]
#    scalefactor = 1.
#    mpl_utils.errorbar(ax, bins, curve.avg, yerr=curve.std_mean, label = label,
#                           marker = None, color = colors[i], every = 2)
#    mpl_utils.errorbar(ax2, bins, phi_vessels.avg, yerr=phi_vessels.std_mean, label = None,
#                           marker = None, color = colors[i], every = 2, ls = ':')
#    mask = (bins>0) & ~curve.avg.mask
#    params, func = FitExpFunction(bins[mask], curve.avg[mask], (30., 0., 100.))
#    x = np.linspace(bins[mask].min(), bins[mask].max(), 200)
#    ax.plot(x, func(x), color = colors[i], linewidth=2., label = '~exp(-x/%s)'%f2l(params[2]))
#  ax.legend()
#  ax2.set(ylim = (0., 1.))
#  mpl_utils.add_crosshair(ax, (0.,None), ls=':')
#  pdfwriter.savefig(fig, postfix='_o2vsvessels')



def PlotLocalScatterData(pdfwriter, smpl, (items, path, label)):
  mask = myutils.bbitwise_and(smpl['flags'], krebsutils.CIRCULATED)
  for k, v in smpl.items():
    smpl[k] = v[mask]

  mask_arterial = myutils.bbitwise_and(smpl['flags'], krebsutils.ARTERY)
  mask_capillaries = myutils.bbitwise_and(smpl['flags'], krebsutils.CAPILLARY)

  if True: # has tumor
    mask_tumor = myutils.bbitwise_and(smpl['flags'], krebsutils.WITHIN_TUMOR)
    mask_arterial = mask_arterial & ~mask_tumor
  else:
    q = smpl['radius'][mask_capillaries]
    smpl['radius'][mask_capillaries] = q*myutils.random_sign(q.shape, q.dtype)
  smpl['radius'][mask_arterial] *= -1.

  def plot(ax, smpl, namex, namey, mask, color, **kwargs):
    x = smpl[namex]
    y = smpl[namey]
    s = smpl['weight']
    if mask is not None:
      x, y ,s = x[mask], y[mask], s[mask]
    if len(s)<=0 or np.amax(s)<=0.: return
    s = np.clip(s*1./100, 0.2, 10.)
    #if 'label' in kwargs:
    #  ax.scatter(x[:1], y[:1], s = 10., alpha=1., marker='o', color=color, label = kwargs.pop('label'), **kwargs)
    ax.scatter(x, y, s=s, alpha=1., rasterized=True, marker='o', color=color, **kwargs)

  charsize = 12./90.
  fig, axes = pyplot.subplots(2, 1, dpi = 300, figsize = mpl_utils.a4size*np.asarray([0.4, 0.35]))
  if True:
    plot(axes[0], smpl, 'radius', 'po2', ~mask_tumor, 'k', label = 'normal')
    plot(axes[1], smpl, 'radius', 'sat', ~mask_tumor, 'k')
    plot(axes[0], smpl, 'radius', 'po2',  mask_tumor, 'r', label = 'tumor')
    plot(axes[1], smpl, 'radius', 'sat',  mask_tumor, 'r')
  else:
    plot(axes[0], smpl, 'radius', 'po2', None, 'k')
    plot(axes[1], smpl, 'radius', 'sat', None, 'k')

  for ax in axes:
    ax.set(xlabel = r'r [$\mu m$]')
    ax.yaxis.grid(True, which='major', color = '#f0f0f0', ls = '-')
    ax.xaxis.grid(True, which='major', color = '#f0f0f0', ls = '-')
  axes[0].set(ylabel = r'P [mmHg]')
  axes[1].set(ylabel = r'S')
  axes[0].legend(loc = mpl_utils.loc.lower_right)
  pyplot.tight_layout()
  pdfwriter.savefig(fig, postfix='_o2scatter1')
  pyplot.close("all")

  fig, axes = pyplot.subplots(2,1, dpi = 300, figsize = mpl_utils.a4size*np.asarray([0.4, 0.35]))
  if True:
    plot(axes[0], smpl, 'velocity', 'po2', ~mask_tumor, 'k', label = 'normal')
    plot(axes[0], smpl, 'velocity', 'po2', mask_tumor, 'r', label = 'tumor')
    plot(axes[1], smpl, 'velocity', 'sat', ~mask_tumor, 'k')
    plot(axes[1], smpl, 'velocity', 'sat', mask_tumor, 'r')
  else:
    plot(axes[0], smpl, 'velocity', 'po2', None, 'r')
    plot(axes[1], smpl, 'velocity', 'sat', None, 'r')
  for ax in axes:
    ax.set(xscale = 'log', xlabel = r'Velocity [$\mu m / s$]')
    ax.yaxis.grid(True, which='major', color = '#f0f0f0', ls = '-')
    ax.xaxis.grid(True, which='major', color = '#f0f0f0', ls = '-')
  axes[0].set(ylabel = r'P [mmHg]')
  axes[1].set(ylabel = r'S')
  axes[0].legend(loc = mpl_utils.loc.lower_right)
  pyplot.tight_layout()
  pdfwriter.savefig(fig, postfix='_o2scatter2')
  pyplot.close("all")


  fig, axes = pyplot.subplots(2,1, dpi = 300, figsize = mpl_utils.a4size*np.asarray([0.4, 0.35]))
  mpl_utils.subplots_adjust_abs(fig, left=6*charsize, right=-2*charsize, top=-5.*charsize, bottom=5.*charsize, hspace=20.*charsize)
  #if has_tumor:
  if True:
    plot(axes[0], smpl, 'flow', 'po2', ~mask_tumor, 'k', label = 'normal')
    plot(axes[0], smpl, 'flow', 'po2', mask_tumor, 'r', label = 'tumor')
    plot(axes[1], smpl, 'flow', 'sat', ~mask_tumor, 'k')
    plot(axes[1], smpl, 'flow', 'sat', mask_tumor, 'r')
  else:
    plot(axes[0], smpl, 'flow', 'po2', None, 'r')
    plot(axes[1], smpl, 'flow', 'sat', None, 'r')
  for ax in axes:
    ax.set(xscale = 'log', xlabel = r'Flow [$\mu m^3 / s$]')
    ax.yaxis.grid(True, which='major', color = '#f0f0f0', ls = '-')
    ax.xaxis.grid(True, which='major', color = '#f0f0f0', ls = '-')
  axes[0].set(ylabel = r'P [mmHg]')
  axes[1].set(ylabel = r'S')
  axes[0].legend(loc = mpl_utils.loc.lower_right)
  pyplot.tight_layout()
  pdfwriter.savefig(fig, postfix='_o2scatter3')
  pyplot.close("all")


  fig, ax = pyplot.subplots(1,1, dpi = 300, figsize = mpl_utils.a4size*np.asarray([0.4, 0.17]))
  plot(ax, smpl, 'radius', 'velocity', ~mask_tumor, 'k', label = 'normal')
  plot(ax, smpl, 'radius', 'velocity', mask_tumor, 'r', label = 'tumor')

  ax.set(yscale = 'log', xlabel = r'Radius [$\mu m$]', ylabel = r'Velocity [$\mu m/s$]')
  ax.yaxis.grid(True, which='major', color = '#f0f0f0', ls = '-')
  ax.xaxis.grid(True, which='major', color = '#f0f0f0', ls = '-')
  
  ax.set()
  ax.legend(loc = mpl_utils.loc.lower_right)
  mpl_utils.tight_layout(fig)
  pdfwriter.savefig(fig, postfix='_o2scatter4')
  pyplot.close("all")


def FormatParameters_(parameters):
  #mbr = lambda s: Prettyfier.mm(Prettyfier.br(s))
  usingMM = parameters.get('michaelis_menten_uptake', False)
  if usingMM:
    parameters = filter(lambda (k,v): not k.startswith('rd_'), parameters.iteritems())
  else:
    parameters = filter(lambda (k,v): not k.startswith('mmcons_'), parameters.iteritems())
  parameters = sorted(parameters, key = lambda (k, v): k)
  #result = ['++++ Oxygen Parameters ++++']
  result = []
  for k, v in parameters:
    k = 'o2p_'+k
    try:
      v, unit = Prettyfier.get_value_unit(k, v, mathMode = True, brackets = True)
      sym = Prettyfier.get_msym(k)
    except KeyError:
      continue
    if isinstance(v, float):
      v = f2l(v)
    else:
      v = str(v)
    #result.append(name + ':')
    result.append('%s = $%s$ %s' % (sym, v, unit))
  return result

def FormatParameters(po2group):
  parameters = detailedo2.readParameters(po2group)
  return FormatParameters_(parameters)


def PlotGlobalData(pdfwriter, items0, t0, data0glob, items1, t1, data1glob, data1tumor):
  multipliers = collections.defaultdict(lambda : float(1))

  FmtTime = lambda time: 't='+f2l(round(time))

  def ToStr(name, avg, std):
    avg = avg*multipliers[name]
    std = std*multipliers[name]
    return f2l(avg), f2l(std)

  def FormatStuff(data):
    #mbr = lambda s: Prettyfier.mm(Prettyfier.br(s))
    result_string = []
    sorted_data = sorted(data.items(), key = lambda (k,v): k)
    for name, values in sorted_data:
      #if name in 'Jout_tv Jin_root Jout_root Jout_cons': continue
      values, unit = Prettyfier.get_value_unit(name, values, mathMode = True, brackets = True)
      avg, std = np.average(values), np.std(values)
      avg, std = ToStr(name, avg, std)
      result_string.append(r'$<%s>$ = $%s \pm %s$ %s' %
        (Prettyfier.get_sym(name), avg, std, unit))
    return result_string

  result_string = ['', '  Parameters  '.center(20, '+') ]
  result_string += FormatParameters(items1[0].po2group)
  if data0glob:
    result_string += ['', ('Initial Network '+FmtTime(t0)).center(20, '+') ]
    result_string += FormatStuff(data0glob)

  fig, _ = mpl_utils.MakeTextPage(result_string, figsize = (mpl_utils.a4size[0], mpl_utils.a4size[1]))
  mpl_utils.subplots_adjust(fig, abs=True, unit='mm', top=-5, bottom=5)
  pdfwriter.savefig(fig, postfix='_o2ginitlobal')

  result_string = [ FmtTime(t1), '  Tumor  '.center(20, '+') ]  
  result_string += FormatStuff(data1tumor)
  
  fig, _ = mpl_utils.MakeTextPage(result_string, figsize = (mpl_utils.a4size[0], mpl_utils.a4size[1]))
  mpl_utils.subplots_adjust(fig, abs=True, unit='mm', top=-5, bottom=5)
  pdfwriter.savefig(fig, postfix='_o2tum')

  result_string = [ FmtTime(t1), '  Global  '.center(20, '+') ]
  result_string += FormatStuff(data1glob)

  fig, _ = mpl_utils.MakeTextPage(result_string, figsize = (mpl_utils.a4size[0], mpl_utils.a4size[1]))
  mpl_utils.subplots_adjust(fig, abs=True, unit='mm', top=-5, bottom=5)
  pdfwriter.savefig(fig, postfix='_o2tum')




class PlotCompareInitialNetworksScatter(object):
    def __init__(self, t0, data0glob, t1, data1glob, data1tumor, items0, items1, xdata, ydata):
        self.stuff0 = t0, data0glob, items0
        self.stuff1 = t1, data1glob, data1tumor, items1
        self.xdataName  = xdata
        self.ydataName  = ydata
        self.lastDataDisplayedName = '', ''
        self.fig, self.ax = pyplot.subplots(1,1, figsize = mpl_utils.a4size*np.asarray([0.4,0.22]), dpi = 90.)
        mpl_utils.subplots_adjust(self.fig, abs=True, unit='mm', top=-5, bottom=15, left=15, right=-5)
        self.numDatasetsDisplayed = 0
        self.bounds = None
        if data0glob:
          dataN = dict(data0glob)
          labelDataN = 'Initial'
        else:
          dataN = dict(data1glob)
          labelDataN = 'Global'
        dataT = dict(data1glob)
        dataT.update(data1tumor)
        self.data_ = {
          't' : dataT, 'n' : dataN
        }
        self.dataLabel = { 't' : 'Tumor', 'n' :  labelDataN, 'd' : 'Diff' }
        self.initialVesselTypeData = np.asarray(list(i.initialVesselType for i in items1))
        initialVesselTypes = list(set(self.initialVesselTypeData)) # eliminate duplicates
        self.initialVesselTypes = sorted(initialVesselTypes, key = lambda c: ord(c)) # just sort it
        self.addLegend = False
        self.xunit_multi, self.xunit = Prettyfier.get_value_unit(self.xdataName, 1)
        self.yunit_multi, self.yunit = Prettyfier.get_value_unit(self.ydataName, 1)
        self.initialVesselTypeLabels = dict((s, RewriteVesselLabel(s)) for s in  initialVesselTypes) # A -> RC1, B -> RC2 ... etc


    def data(self, origin, name):
      if origin <> 'd':
        return self.data_[origin][name]
      else:
        return self.data_['t'][name] - self.data_['n'][name]


    def addPlot(self, xnetwork, ynetwork, colored = True, legend = False):
      self.addLegend = self.addLegend or legend
      #http://xkcd.com/color/rgb/
      if colored:
        colors = vesselTypeColors 
      else:
        colors = alternateVesselTypeColors
      #t0, data0glob, data0normal, items0 = self.stuff0
      #t1, data1glob, data1tumor, data1normal, items1 =self.stuff1
      ax = self.ax

      self.lastDataDisplayedName = xnetwork, ynetwork

      dataX = self.xunit_multi*self.data(xnetwork, self.xdataName)
      dataY = self.yunit_multi*self.data(ynetwork, self.ydataName)
      
      self.bounds = mpl_utils.AccumulateXYBounds(np.amin(dataX), np.amax(dataX), np.amin(dataY), np.amax(dataY), previous = self.bounds)
      
      for i, c in enumerate(self.initialVesselTypes):
        mask = c == self.initialVesselTypeData
        label = self.initialVesselTypeLabels[c] if legend else None
        ax.scatter(dataX[mask], dataY[mask], marker = vesselTypeMarkers[i], c = colors[i], facecolor='none', edgecolor = colors[i], linewidth=1., s = 12, label = label)

      self.numDatasetsDisplayed += 1
      return self
    
    def addLinearFit(self, legend = False):
      self.addLegend = self.addLegend or legend
      xnetwork, ynetwork = self.lastDataDisplayedName
      dataX = self.xunit_multi*self.data(xnetwork, self.xdataName)
      dataY = self.yunit_multi*self.data(ynetwork, self.ydataName)
      def fun(x, p):
        return x*p[0]
      def objective(p):
        return fun(dataX, p) - dataY
      (params, _) = scipy.optimize.leastsq(objective, (1,))
      xb, yb = self.bounds
      x = np.linspace(xb[0], xb[1], 2)
      y = fun(x, params)
      label = ('$%s \cdot %s$' % (f2l(params[0]), Prettyfier.get_sym(self.xdataName))) if legend else None
      self.ax.plot(x, y, label = label, color = 'r', lw = 1.5)
      return self
    
    def addFinalStuff(self):
      ax = self.ax
      #print self.xdata, self.ydata, self.bounds, xb, yb
      xb, yb = self.bounds
      #if self.xdata == self.ydata:
      #  xb = min(xb[0], yb[0]), xb[1]
      #  yb = min(xb[0], yb[0]), yb[1]
      mpl_utils.SetSensibleAxLimits(ax, xb, yb)
      
      if self.xdataName == self.ydataName and self.numDatasetsDisplayed==1: # add diagonal
        ax.plot(xb, xb, color = 'k', lw = 0.5)

      if self.addLegend:
        ax.legend(loc=mpl_utils.loc.upper_right, prop={'size':6})

      ax.grid(linestyle=':', linewidth=0.5, color = '0.5')

      br = lambda s: ('$[%s]$' % s) if s else ''
      xnetwork, ynetwork = self.lastDataDisplayedName
      if self.numDatasetsDisplayed == 1:
        ax.set(xlabel = '$%s$ %s (%s)' % (Prettyfier.get_sym(self.xdataName), br(self.xunit), self.dataLabel[xnetwork]),
               ylabel = '$%s$ %s (%s)' % (Prettyfier.get_sym(self.ydataName), br(self.yunit), self.dataLabel[ynetwork]))
      else:
        ax.set(xlabel = '$%s$ %s' % (Prettyfier.get_sym(self.xdataName), br(self.xunit)),
               ylabel = '$%s$ %s' % (Prettyfier.get_sym(self.ydataName), br(self.yunit)))
      return self
    
    def ticks(self, axis, count):
      loc = matplotlib.ticker.MaxNLocator(nbins=count)
      ax = self.ax.xaxis if axis=='x' else self.ax.yaxis
      ax.set_major_locator(copy(loc))
      return self
    
    def write(self, pdfwriter):
      self.addFinalStuff()
      xnetwork, ynetwork = self.lastDataDisplayedName
      pdfwriter.savefig(self.fig, dpi = 90., postfix='_%s_%s_%s_%s' % (xnetwork, self.xdataName, ynetwork, self.ydataName))
      pyplot.close(self.fig)
  




def ExportGlobalDataAsText(filenamebase, snapshotlist, measurementinfo): #, ensemble_items, t0, data0glob, data0normal, t1, data1glob, data1tumor, data1normal)
  FmtTime = lambda time: f2s(round(time))
  s = []
  out = s.append
  for ensembleitems, globaldata, normaldata, tumordata, time in snapshotlist:
    out('t='+FmtTime(time))
    out('\n')
    data = globaldata.values() + normaldata.values() + tumordata.values()
    data = zip(*data)
    labels = map(lambda s: 'global_'+s, globaldata.keys()) + map(lambda s: 'normal_'+s, normaldata.keys()) + map(lambda s: 'tumor_'+s, tumordata.keys())
    out(' '.join(labels))
    out('\n')
    for l in data:
      out(' '.join(map(str,l)))
      out('\n')
    out('\n')
  with open(filenamebase+'_data.txt', 'w') as f:
    f.writelines(s)
  del s[:]
  labels = set(globaldata.keys() + normaldata.keys() + tumordata.keys())
  unitlabels = map(lambda s: 'unit "%s"="%s"\n' % (s, Prettyfier.get_unit(s)), labels)
  for x in unitlabels:
    out(x)
  with open(filenamebase+'_unit.txt', 'w') as f:
    f.writelines(s)


def ExportGlobalDataAsH5(filenamebase, snapshotlist, measurementinfo): #, ensemble_items, t0, data0glob, data0normal, t1, data1glob, data1tumor, data1normal)  
  try:
    os.remove(filenamebase+'_globalh5export.h5')
  except OSError:
    pass
  with h5py.File(filenamebase+'_globalh5export.h5', 'w-') as f:
    for i, (ensembleitems, globaldata, tumordata, time) in enumerate(snapshotlist):
      gbase = f.create_group(['initial', 'final'][i])
      for name, data in zip(['global', 'tumor'], [globaldata, tumordata]):
        g = gbase.create_group(name)
        for k, v in data.iteritems():
          g.create_dataset(k, data = v)
        
    filenames = map(lambda item: item.po2group.file.filename.encode('ascii'), ensembleitems)
    f.create_dataset('filenames', data = filenames)
      

#def make_radial_check_plots(pdfwriter, po2groups):
#  bin_spec = BinsSpecRange(-10000., 10000., 100.)
#
#  curves = collections.defaultdict(list)
#  for g in po2groups[-1:]:
#    gvessels, gtumor = detailedo2.OpenVesselAndTumorGroups(g)
#    pos_smpl  = dataman.basic_vessel_samples('position', gvessels, sample_length)
#    dist_smpl, distmap, mask   = dataman.obtain_data('distancemap_samples', gvessels, gtumor, sample_length, 'levelset')
#    tumor_ld = dataman.obtain_data('ld', gtumor.file)
#    x = pos_smpl[:,0]
#    y = pos_smpl[:,1]
#    z = pos_smpl[:,2]
#    dist_smpl[~mask] = 10000.
#    mask = (z>-10.) & (z<10.)
#    x = x[mask]
#    y = y[mask]
#    dist_smpl = dist_smpl[mask]
#    print distmap.max()
#    print tumor_ld
#    print x.min(), x.max()
#    plotBulkTissue.imshow(pyplot.gca(), plotBulkTissue.imslice(distmap), tumor_ld, vmin=np.min(distmap), vmax=np.max(distmap))
#    pyplot.scatter(x,y, s=1., alpha=1., rasterized=True, marker='x', c=dist_smpl, cmap=matplotlib.cm.jet, vmin=np.min(distmap), vmax=np.max(distmap))
#    pyplot.show()
##    data = dataman.obtain_data('detailedPO2_radial', g, sample_length, bin_spec, distancemap_spec, cachelocation(g.po2group))
##    for k,v in data.iteritems():
##      curves[k].append(v.cnt)
#    for name in ['mvd','velocity']:
#      data = dataman.obtain_data('basic_vessel_radial', name, gvessels, gtumor, sample_length, bin_spec, 'levelset', None, cachelocation(g))
#      curves[name].append(data.cnt)
#
#  num_vessel_smpl = myutils.MeanValueArray.fromSummation(curves['velocity'])
#  num_vol_smpl    = myutils.MeanValueArray.fromSummation(curves['mvd'])
#
#  bins = bin_spec.arange()
#  bins = 1.e-3*np.average((bins[1:],bins[:-1]), axis=0)
#
#  pyplot.close()
#  figsize = mpl_utils.a4size[0], mpl_utils.a4size[1]*0.5
#  charsize = 12./90.
#  fig, ax = pyplot.subplots(3,1, figsize = figsize, dpi = 90.)
#  mpl_utils.subplots_adjust_abs(fig, left=5*charsize, right=-2*charsize, top=-3.*charsize, bottom=5.*charsize, wspace=18.*charsize)
#  ax[0].plot(bins, np.power(num_vessel_smpl.avg, 1./2.))
#  #ax[0].set(yscale = 'log')
#  ax[0].plot(bins, np.power(num_vol_smpl.avg, 1./2.))
#  #ax[1].set(yscale = 'log')
#  ax[1].plot(bins, num_vessel_smpl.avg/num_vol_smpl.avg)
#  pyplot.show()


def PlotHistograms(pdfwriter, mainCombinedHistogramGroup, region, title):
  names = ['sat', 'po2_tissue', 'radius', 'velocity']
  # regions are 'all', 'tum', 'norm'
  colors = [
     '#e50000', #red
     '#9a0eea', #violet
     '#0343df', #blue
     '#f97306', #orange
     '#677a04', #olive green
     '#ceb301', #mustard
     '#04d8b2', #aquamarine
     '#06470c', #forest greeen
     '#840000', #dark red
     '#607c8e', #blue grey
  ]
  def scatterInStyle(ax, x, y, i, label):
    ax.scatter(x, y, marker = '<>*osd^+x'[i], c = colors[i], facecolor='none', edgecolor = colors[i], linewidth=0.5, s = 5., label = label, zorder = 2)
  
  networkTypes = 'ABCDEFGHI'
  
  fig, axes = pyplot.subplots(4,1, figsize = mpl_utils.a4size*np.asarray([0.33,0.5]), dpi = 90.)
  fig.suptitle(title)
  
  for i, name in enumerate(names):
    histo = myutils.MeanValueArray.read(mainCombinedHistogramGroup['ensemble'][name][region])
    histo *= 1.0/np.sum(histo.sum) # do normalization here
    bins  = np.asarray(mainCombinedHistogramGroup['ensemble'][name][region]['bins'])
    ax = axes[i]
    ax.bar(bins[:-1], histo.sum, width = (bins[1:]-bins[:-1]), color = 'white', zorder = 1)
    for j, network in enumerate(networkTypes):
      histo = myutils.MeanValueArray.read(mainCombinedHistogramGroup['networkTypes'][network][name][region])
      histo *= 1.0/np.sum(histo.sum) # do normalization here
      bins  = np.asarray(mainCombinedHistogramGroup['networkTypes'][network][name][region]['bins'])
      x = 0.5*(bins[:-1]+bins[1:])
      scatterInStyle(ax, x, histo.sum, j, RewriteVesselLabel(network))
    xlabel = '$%s$ %s' % (Prettyfier.get_sym(name), Prettyfier.get_bunit(name))
    ax.set(xlabel = xlabel, ylabel = 'freq.', xlim = (bins[0], bins[-1]))
  
  pyplot.tight_layout()
  pdfwriter.savefig(fig)



##############################################################################
##############################################################################
sample_length = 30.


# note: use insertVesselConfigInO2File.py script (in scripts repository) to copy vesselfile message over. It is the type name in it.
def GetVesselTypeLabel(po2group):
  import re
  if ('VESSELFILE_MESSAGE' in po2group.file.attrs.keys()):
    msg = po2group.file.attrs['VESSELFILE_MESSAGE']
    m = re.search('type(\w)', msg)
  else:
    m = re.search('type(\w)', po2group.file.filename)
  if not m:
    return 'unkown'
  return m.group(1)


class EnsembleItem(object):
  def __init__(self, **kwargs):
    self.path = ''
    self.time = None
    self.po2group = None
    self.gvessels = None
    self.gtumor = None
    self.vessel_system_length = 0.
    self.initialVesselType = ''
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
        po2group = f[path]
        gvessels, gtumor = detailedo2.OpenVesselAndTumorGroups(po2group)
        e = EnsembleItem(path = path, po2group = f[path], gvessels = gvessels, gtumor = gtumor)
        if 'SOURCE' in po2group:
          #print("attention thierry hack!!!")
          if( 'time' in po2group['SOURCE'].attrs.keys()):
            source = h5files.openLink(po2group, 'SOURCE')
            t = source.attrs['time']
            e.time = t
        has_tumor = has_tumor and gtumor is not None
        e.vessel_system_length = dataman.obtain_data('vessel_system_length', gvessels)
        e.initialVesselType = GetVesselTypeLabel(po2group)
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
    self.has_tumor = has_tumor
    if has_tumor:
      self.tumor_snapshots = tumor_snapshots # list of tuple(items, path, time)
    
    self.o2ConfigName = set(item.po2group.attrs.get('O2_CONFIG_NAME',None) for item in items)
    if len(self.o2ConfigName) != 1:
      raise RuntimeError("Detected different O2_CONFIG_NAMES %s. You don't want to mix configurations, do you?" % self.o2ConfigName)
    self.o2ConfigName = self.o2ConfigName.pop()


class MeasurementInfo(object):
  def __init__(self, **kwargs):
    self.sample_length = 30.
    self.cachelocation_callback = None
    self.distancemap_spec = 'radial'
    for k,v in kwargs.iteritems():
      setattr(self, k, v)

#
#def MakeBinSpec(distancemap_spec, gtumor, only_tumor_vs_normal):
#  if gtumor:
#    if only_tumor_vs_normal:   # two bins: tumor - normal
#      if distancemap_spec == 'radial':
#        r = dataman.obtain_data('approximate_tumor_radius', item.gtumor)
#        bin_spec = BinsSpecArray([-100000., r, 40000.])
#      else: # levelset as distance -> centered around 0
#        bin_spec = BinsSpecArray([-100000.,0., 100000])
#    else: # many small bins
#      bin_spec = BinsSpecRange(-100000., 100000., 100.)
#  else:
#    assert only_tumor_vs_normal is None
#    assert distancemap_spec == 'radial'
#    bin_spec = BinsSpecArray([-100000.,0., 100000])
#  return bin_spec


def GetAverageApproximateTumorRadius(dataman, ensembleitems):
  l = map(lambda item: dataman.obtain_data('approximate_tumor_radius', item.gtumor),
          ensembleitems)
  return np.average(l)


def CollectAllRadialData(dataman, ensembleitems, measurementinfo):
  bin_spec = BinsSpecRange(-100000., 100000., 100.) #MakeBinSpec(distancemap_spec, item.gtumor, only_tumor_vs_normal)
  curves = collections.defaultdict(list)
  for item in ensembleitems:
    # get o2 data first
    print 'generating radial curves for',item.po2group.file.filename,':',item.po2group.name
    data = dataman.obtain_data('detailedPO2_radial', item.po2group, measurementinfo.sample_length, bin_spec, measurementinfo.distancemap_spec, measurementinfo.cachelocation_callback(item.po2group))
    try:
      del data['extpo2'] # delete if present, else ignore
    except KeyError:
      pass
    for k,v in data.iteritems():
      curves[k].append(v)
    for name in ['mvd','velocity','phi_vessels', 'shearforce', 'radius']:
      data = dataman.obtain_data('basic_vessel_radial', name, item.gvessels, item.gtumor, measurementinfo.sample_length, bin_spec, measurementinfo.distancemap_spec, None, measurementinfo.cachelocation_callback(item.po2group))
      curves[name].append(data)
  print ' ... finished getting radial curves'
  return bin_spec, curves


def ConvertDictValuesToNumpyArrays_(x):
  for k, v in x.iteritems():
    x[k] = np.asarray(v)
  return x


def CollectAllGlobalData(dataman, ensembleitems, measurementinfo):
  # these are maps from 'dataname' -> list of data points; one point per tumor
  curves_tumor = collections.defaultdict(list)
  curves_total = collections.defaultdict(list)
  # go and collect all data
  for item in ensembleitems:
    print 'generating global data for',item.po2group.file.filename,':',item.po2group.name
    # first real global data
    for prop in ['e1','e2','e3','Jin_root', 'Jout_root', 'Jout_tv', 'tv_cons', 'Jout_cons']:  #'po2','sat','gtv', 'jtv','mro2', 'po2_tissue', 'oef', 'chb', 'chb_oxy', 'chb_deoxy', 'sat_via_hb_ratio'
      data = dataman.obtain_data('detailedPO2_global', prop, item.po2group, measurementinfo.sample_length, measurementinfo.cachelocation_callback(item.po2group))
      curves_total[prop].append(data)
    # then radial curves, with two bins, one for tumor and one for normal tissue
    r = dataman.obtain_data('approximate_tumor_radius', item.gtumor)
    curves_total['approximate_tumor_radius'].append(r)
    if measurementinfo.distancemap_spec == 'radial':
      bin_spec_tum = BinsSpecArray([-100000., r, 40000.])
    else: # levelset as distance -> centered around 0
      bin_spec_tum = BinsSpecArray([-100000.,0., 100000])
    bin_spec_total = BinsSpecArray([-1000000.,1000000]) # everything
    # get the data and put it in the right list
    for curves, bin_spec in zip([curves_total, curves_tumor], [bin_spec_total, bin_spec_tum]):
      datadict = dataman.obtain_data('detailedPO2_radial', item.po2group, measurementinfo.sample_length, bin_spec, measurementinfo.distancemap_spec, measurementinfo.cachelocation_callback(item.po2group))
      for name, data in datadict.iteritems():
        curves[name].append(data.avg[0])
      for name in ['mvd','velocity','phi_vessels', 'shearforce', 'radius', 'S_rho', 'hematocrit']:
        data = dataman.obtain_data('basic_vessel_radial', name, item.gvessels, item.gtumor, measurementinfo.sample_length, bin_spec, measurementinfo.distancemap_spec, None, measurementinfo.cachelocation_callback(item.po2group))
        curves[name].append(data.avg[0])
      data = dataman.obtain_data('detailedPO2_peff_radial', item.po2group, measurementinfo.sample_length, bin_spec, measurementinfo.distancemap_spec, measurementinfo.cachelocation_callback(item.po2group))
      curves['Peff'].append(data.avg[0])
      peffSrho = dataman.obtain_data('detailedPO2_peffSrho_radial', item.po2group, measurementinfo.sample_length, bin_spec, measurementinfo.distancemap_spec, measurementinfo.cachelocation_callback(item.po2group))
      curves['PeffSrho'].append(peffSrho.avg[0])    
    del curves, bin_spec
    oefdata = detailedo2Analysis.ObtainOxygenExtractionFraction(dataman, item.po2group, measurementinfo.cachelocation_callback(item.po2group))
    curves_tumor['oef'].append(oefdata['oef_tumor'])
    curves_total['oef'].append(oefdata['oef_total'])
    rBFdata = dataman.obtain_data('blood_flow_rbf', item.gvessels, item.gtumor, measurementinfo.cachelocation_callback(item.po2group))
    curves_tumor['rBF'].append(rBFdata['rBF_tumor'])
    curves_total['rBF'].append(rBFdata['rBF_total'])
    curves_total['Sin'].append(oefdata['total_o2_in']/rBFdata['total_flow_in']*60.*1.e-12)
    curves_tumor['Sin'].append(oefdata['tumor_o2_in']/rBFdata['flow_in']*60.*1.e-12)
    curves_tumor['rJin'].append(oefdata['tumor_o2_in']/rBFdata['tumor_volume']*60.*1.e-9*1.e3)
    curves_total['rJin'].append(oefdata['total_o2_in']/rBFdata['total_volume']*60.*1.e-9*1.e3)

  curves_total = ConvertDictValuesToNumpyArrays_(curves_total)
  curves_tumor = ConvertDictValuesToNumpyArrays_(curves_tumor)
  
  parameters = dataman.obtain_data('detailedPO2Parameters', ensembleitems[0].po2group)
  hema = parameters['calcflow']['inletHematocrit']*np.ones_like(curves_total['Sin'])
  #print detailedo2.ConcentrationToPO2(curves_total['Sin'], hema, parameters)
  curves_total['Sin'] = detailedo2.PO2ToSaturation(detailedo2.ConcentrationToPO2(curves_total['Sin'], hema, parameters), parameters)
  curves_tumor['Sin'] = detailedo2.PO2ToSaturation(detailedo2.ConcentrationToPO2(curves_tumor['Sin'], hema, parameters), parameters)

  for c in [curves_tumor, curves_total]: # this depends on rBF so we have to use the bulk of tissue
    c['mtt']  = c['phi_vessels']/c['rBF']
    c['S_rho_over_rBV'] = c['S_rho']/c['phi_vessels']
    c['PeffSrhoMTT'] = c['PeffSrho']*c['mtt']
    c['kExpFun'] = c['Sin']*(1.0 - np.exp(-c['PeffSrhoMTT']))/c['PeffSrhoMTT']
    c['Y_plus_oef'] = c['sat_via_hb_ratio'] + c['oef']
  
  print ' ... finished getting global data'
  return curves_total, curves_tumor


def CollectInterdependentDependentData(dataman, items0, items1, curves_tumor, measurementinfo):
  from krebs.analyzeBloodFlow import ComputeIsoTumorSpherePerfusionScaleFactor    

  newdata = collections.defaultdict(list)
  for i, item0, item1 in zip(itertools.count(), items0, items1):
    # item 0 is t=0, item1 is final tumor
#    rescaledPerfusion = ComputeIsoTumorSphereRescaledPerfusion(dataman, 
#      item0.gvessels, item1.gvessels, item1.gtumor, 
#      measurementinfo.cachelocation_callback(item0.po2group),
#      measurementinfo.cachelocation_callback(item1.po2group))
#    rescaledPerfusion = float(rescaledPerfusion[...]) * 60. # 1/s -> 1/min

    scaleFactor = ComputeIsoTumorSpherePerfusionScaleFactor(dataman, 
      item0.gvessels, item1.gvessels, item1.gtumor, 
      measurementinfo.cachelocation_callback(item0.po2group),
      measurementinfo.cachelocation_callback(item1.po2group))
    rBF  = curves_tumor['rBF'][i]
    rJin = curves_tumor['rJin'][i]
    OEF  = curves_tumor['oef'][i]
    newdata['scaled_rBF'].append(scaleFactor * rBF)
    newdata['scaled_rJin'].append(scaleFactor * rJin)
    newdata['scaled_oef'].append(OEF / scaleFactor)
  
  newdata = ConvertDictValuesToNumpyArrays_(newdata)
  curves_tumor.update(newdata)
  

#def CollectLocalScatterData(dataman, ensembleitems, measurementinfo, every):
#  data = collections.defaultdict(list)
#  for item in ensembleitems:
#    print 'generating scatter data for',item.po2group.file.filename,':',item.po2group.name
#    for name in ['po2','sat']:
#      d = dataman.obtain_data('detailedPO2_samples', name, item.po2group, measurementinfo.sample_length, every, detailedo2Analysis.MakeSampleLocation(item.po2group))
#      data[name].append(d.copy())
#    for name in [ 'weight', 'flags', 'flow', 'radius', 'velocity', 'shearforce' ]:
#      d = dataman.obtain_data('basic_vessel_samples', name, item.gvessels, sample_length)
#      data[name].append(d[::every].copy()) 
#  for k, v in data.items():
#    data[k] = np.concatenate(v)
#  return data


class ComputesRegionalHistograms(object):
  def __init__(self, dataman, vesselgroup, tumorgroup, sample_length, ld):
    print 'ComputesRegionalHistograms for %s' % (str(vesselgroup))
    self.vesselgroup = vesselgroup
    self.tumorgroup  = tumorgroup
    vesselDistanceSamples, distanceMap, vesselOutsideMask, ld = \
      dataman.obtain_data('distancemap_samples', vesselgroup, tumorgroup, sample_length, 'levelset', ld)
    # note: mask is true for samples inside the distancemap box
    self.maskInTumorVesselSamples = vesselOutsideMask & (vesselDistanceSamples<0.)
    self.maskInTumorTissue        = distanceMap.ravel() < 0.
    self.weight = dataman.obtain_data('basic_vessel_samples', 'weight', vesselgroup, sample_length)
    #self.maskCirculated = analyzeGeneral.GetMaskUncirculated(dataman, vesselgroup)
    flags = dataman.obtain_data('basic_vessel_samples', 'flags', vesselgroup, sample_length)
    self.maskCirculated = myutils.bbitwise_and(flags, krebsutils.CIRCULATED)
    print 'cnt in tum %s of total %s', (np.count_nonzero(self.maskInTumorVesselSamples), len(self.maskInTumorVesselSamples))
  
  def MakeHistogramsForVessels(self, dataSamples, bins):
    mc         = self.maskCirculated
    histoAll   = myutils.MeanValueArray.fromHistogram1d(bins, 
                                                        dataSamples[mc], 
                                                        self.weight[mc])
    histoTum   = myutils.MeanValueArray.fromHistogram1d(bins, 
                                                        dataSamples[mc & self.maskInTumorVesselSamples], 
                                                        self.weight[mc & self.maskInTumorVesselSamples])
    mask       = mc & ~self.maskInTumorVesselSamples
    histoNorm  = myutils.MeanValueArray.fromHistogram1d(bins, 
                                                        dataSamples[mask], 
                                                        self.weight[mask])
    return histoAll, histoTum, histoNorm
  
  def MakeHistogramForTissue(self, dataSamples, bins):
    histoAll  = myutils.MeanValueArray.fromHistogram1d(bins,
                                                       dataSamples.ravel(),
                                                       1.)
    histoTum  = myutils.MeanValueArray.fromHistogram1d(bins,
                                                       dataSamples.ravel()[self.maskInTumorTissue],
                                                       1.)
    histoNorm = myutils.MeanValueArray.fromHistogram1d(bins,
                                                       dataSamples.ravel()[~self.maskInTumorTissue],
                                                       1.)
    return histoAll, histoTum, histoNorm


@myutils.UsesDataManager
def GetComputesRegionalHistograms(dataman, po2group):
  vesselgroup, tumorgroup = detailedo2.OpenVesselAndTumorGroups(po2group)
  ld = krebsutils.read_lattice_data_from_hdf(po2group['field_ld'])
  computesRegionalHistograms = ComputesRegionalHistograms(dataman, vesselgroup, tumorgroup, sample_length, ld)
  return computesRegionalHistograms



def GetDataSamplesForHistograms(dataman, name, po2group, measurementinfo):    
    gvessels, gtumor = detailedo2.OpenVesselAndTumorGroups(po2group)   
    if name == 'po2_tissue':
      return np.asarray(po2group['po2field']).ravel()
    elif name == 'sat':
      return dataman.obtain_data('detailedPO2_samples', name, po2group, measurementinfo.sample_length, 1, detailedo2Analysis.MakeSampleLocation(po2group))
    elif name in ['radius', 'velocity' ]:
      return dataman.obtain_data('basic_vessel_samples', name, gvessels, measurementinfo.sample_length)
    else:
      assert False


@myutils.UsesDataManager
def ComputeRegionalHistogramsOfPo2Group(dataman, name, region, po2group, binspec, measurementinfo):    
  def read(gmeasure, groupname):
    return myutils.MeanValueArray.read(gmeasure[groupname], region)  
  
  def write(gmeasure, groupname):
    print 'computing histogram %s of %s' % (name, str(po2group))
    computesRegionalHistograms = GetComputesRegionalHistograms(dataman, po2group)
    bins = np.linspace(*binspec)
    smpl = GetDataSamplesForHistograms(dataman, name, po2group, measurementinfo)
    if name == 'po2_tissue':
      func = computesRegionalHistograms.MakeHistogramForTissue  
    else:
      func = computesRegionalHistograms.MakeHistogramsForVessels  
    histoAll, histoTum, histoNorm = func(smpl, bins)
    
    gmeasure = gmeasure.create_group(groupname)
    histoAll.write(gmeasure, 'all')
    histoTum.write(gmeasure, 'tum')
    histoNorm.write(gmeasure, 'norm')
  
  cachelocation = measurementinfo.cachelocation_callback(po2group)
  version       = myutils.checksum(binspec, 3)
  return myutils.hdf_data_caching(read, write, cachelocation[0], ('histograms', cachelocation[1], name), (None, None, version))
  

def ComputeHistogramsOfPo2Items(dataman, items, measurementinfo, outputGroup):
    rangesAndBinCount = {
      'sat' : (0., 1.,20),
      'po2_tissue' : (0., 50., 25),
      'radius' : (0., 20., 20),
      'velocity' : (0., 1000., 20),
    }
    combinedHistograms = collections.defaultdict(lambda: myutils.MeanValueArray.empty())
    names = rangesAndBinCount.keys()
    combinedHistogramsByNetworkType = collections.defaultdict(lambda: myutils.MeanValueArray.empty())
    for item in items:
      for name in names:
        for region in ['all','tum', 'norm']:
          h = ComputeRegionalHistogramsOfPo2Group(dataman, name, region, item.po2group, rangesAndBinCount[name], measurementinfo)
          combinedHistograms[name, region] += h
          combinedHistogramsByNetworkType[name, region, item.initialVesselType] += h
    allInitialVesselTypes = set(t for (_,_,t) in combinedHistogramsByNetworkType.iterkeys())

    try:
      del outputGroup['ensemble']
    except KeyError:
      pass
    try:
      del outputGroup['networkTypes']
    except KeyError:
      pass    
    
    outputGroup.file.flush()
    for name in names:
      bins = np.linspace(*rangesAndBinCount[name])
      for region in ['all','tum','norm']:
        combinedHistograms[name,region].write(outputGroup,'ensemble/%s/%s' % (name, region))
        outputGroup.create_dataset('ensemble/%s/%s/bins' % (name, region), data = bins)
        for vt in allInitialVesselTypes:
          combinedHistogramsByNetworkType[name, region, vt].write(outputGroup, 'networkTypes/%s/%s/%s' % (vt, name, region))
          outputGroup.create_dataset('networkTypes/%s/%s/%s/bins' % (vt, name, region), data = bins)
        
    outputGroup.file.flush()
  



def GetConvergenceData(po2group):
  if not 'iterations' in po2group: return None
  def read(gmeasure, name):
    return dict(map(lambda (k,v): (k,np.asarray(v)), gmeasure[name].iteritems()))
  def write(gmeasure, name):
    def read_iter_snapshot(g):
      return dict(map(lambda (k,v): (k, float(v[()])), g.iteritems()))
    data = map(read_iter_snapshot, po2group['iterations'].itervalues())
    data = sorted(data, key = lambda d: d['iteration'])
    data = myutils.zipListOfDicts(data, numpy_output=True)
    gmeasure = gmeasure.create_group(name)
    for k, v in data.iteritems():
      gmeasure.create_dataset(k, data = v, compression = 9)
  return myutils.hdf_data_caching(read, write, po2group, ('iterationsTranspose',), (1,))


def PlotConvergenceData(pdfwriter, dataman, (items, path, time)):
  fig = pyplot.figure(figsize = mpl_utils.a4size*np.asarray((0.5,0.25)))
  ax  = fig.add_axes([0.2, 0.2, 0.7, 0.6])
  for item in items:
    iterations = GetConvergenceData(item.po2group)
    ycurve = iterations['delta_vessM'] + iterations['delta_fieldM']
    ax.plot(iterations['iteration'], ycurve, color = 'k')
  title = 'Convergence t = $%s$' % f2l(time)
  ax.set(xlabel = 'Iteration $n$', ylabel = r'$|P^{(n+1)} - P^{(n)}|_\infty + |P_t^{(n+1)} - P_t^{(n)}|_\infty$', title = title, yscale = 'log')
  pdfwriter.savefig(fig, postfix='_convergence_'+path)


def doit(filenames, pattern, normalTissueEnsembleInput = None):
  dataman = myutils.DataManager(50, map(lambda x: x(), detailedo2Analysis.O2DataHandlers) + [ analyzeGeneral.DataTumorTissueSingle(), analyzeGeneral.DataDistanceFromCenter(), analyzeGeneral.DataBasicVessel(), analyzeGeneral.DataVesselSamples(), analyzeGeneral.DataVesselRadial(), analyzeGeneral.DataVesselGlobal(), analyzeBloodFlow.DataTumorBloodFlow()])
  ensemble = EnsembleFiles(dataman, filenames, pattern)
  
  if ensemble.has_tumor:
    if len(ensemble.tumor_snapshots) == 1:
      print 'paths: ', map(lambda (_t0, path, _t1): path, ensemble.tumor_snapshots)
      items0, _, t0 = None, None, None
      items1, _, t1 = ensemble.tumor_snapshots[0] 
    else:
      items0, _, t0 = ensemble.tumor_snapshots[0]  # initial
      items1, _, t1 = ensemble.tumor_snapshots[-1] # final
  else:
    normalTissueEnsembleInput = filenames, pattern
  if normalTissueEnsembleInput:
    filenamesNormal, patternNormal = normalTissueEnsembleInput
    normalTissueEnsemble = EnsembleFiles(dataman, filenamesNormal, patternNormal)
    items0, _, t0 = normalTissueEnsemble.tumor_snapshots[0]

  out_prefix, out_suffix = myutils.splitcommonpresuffix(map(lambda s: basename(s), filenames))
  output_base_filename = splitext(out_prefix+out_suffix)[0]
  if ensemble.o2ConfigName:
    fn_measure = 'detailedo2_%s_common.h5' % ensemble.o2ConfigName
  else:
    fn_measure = 'detailedo2_common.h5'
    
  f_measure = h5files.open(fn_measure, 'a')
  def cachelocation(g):
    path = posixpath.join('FileCS_'+myutils.checksum(basename(g.file.filename)), g.name.strip(posixpath.sep))
    return (f_measure, path)

  measurementinfo = MeasurementInfo(sample_length = 30.,
                                    cachelocation_callback = cachelocation,
                                    distancemap_spec = 'radial')

  with mpl_utils.PageWriter(output_base_filename+'.pdf', fileformats = ['pdf']) as pdfwriter:
  
#    if 0:
#      groups = groups_by_path[pathorder[-1]]
#      po2groups = [v.po2group for v in groups]
#      binspec_tumor = BinsSpecArray([-10000., -200., 200., 10000])
#      binspec_vessels = BinsSpecRange(-10000., 10000., 30.)
#      data = CalcOxygenVsVesselsCached(dataman, po2groups, binspec_tumor, binspec_vessels, (f_measure, 'files'+myutils.checksum([g.name for g in po2groups])))
#      PlotOxygenVsVessels(pdfwriter, data, binspec_vessels)

    if 0:
        print 'getting radial curves'
        bins_spec, curves0 = CollectAllRadialData(dataman, items0, measurementinfo)
        bins_spec, curves1 = CollectAllRadialData(dataman, items1, measurementinfo)
        tumor_radius0 = GetAverageApproximateTumorRadius(dataman, items0)
        tumor_radius1 = GetAverageApproximateTumorRadius(dataman, items1)
        PlotRadialCurves(pdfwriter, bins_spec, [(t0, tumor_radius0, curves0), (t1, tumor_radius1, curves1)], measurementinfo)

    if 0:
      every = int(np.sum(item.vessel_system_length for item in ensemble.tumor_snapshots[-1][0])/(30.*20000.))
      print 'getting samples, every = %i' % every
      smpl = CollectLocalScatterData(dataman, ensemble.tumor_snapshots[-1][0], measurementinfo, every)
      print 'making scatterplots'
      PlotLocalScatterData(pdfwriter, smpl, ensemble.tumor_snapshots[-1])
      del smpl


    if 1:
      print 'getting global data'
      if items0:
        data0glob, data0tumor = CollectAllGlobalData(dataman, items0, measurementinfo)
      else:
        data0glob, data0tumor = None, None
      data1glob, data1tumor = CollectAllGlobalData(dataman, items1, measurementinfo)
      CollectInterdependentDependentData(dataman, items0, items1, data1tumor, measurementinfo)

      if 1:
        PlotGlobalData(pdfwriter, items0, t0, data0glob, items1, t1, data1glob, data1tumor)
        #to do:
        # k=Peff*S*rho*MTT 
        # S*rho die Gefäßoberfläche pro Gewebevolumen
        # P= mass transfer coefficient /oxygen solubility in blood
        # P=gamma/alpha=Sauerstoff-Diffusionskoeffizient* Nusselt number/Radius P=D_p*nu/r
        #  P_eff=P*beta/(1+beta)
        #  beta = sol. O2 (Plasma) / bound O2 (Hemoglobin)
        #  rho = 1/gewebevolumen
        # plot S*rho/rBV
        def plotit(x,y):
          return PlotCompareInitialNetworksScatter(t0, data0glob, t1, data1glob, data1tumor, items0, items1, x, y)
        plotit('chb', 'chb').addPlot('n','t', legend = True).write(pdfwriter)
        plotit('sat_via_hb_ratio', 'sat_via_hb_ratio').addPlot('n','t').write(pdfwriter)
        plotit('chb', 'sat_via_hb_ratio').addPlot('t','t').addPlot('n','n', colored = False).write(pdfwriter)
        plotit('rBF', 'sat_via_hb_ratio').addPlot('n','n').ticks('x',6).write(pdfwriter)
        plotit('rBF', 'sat_via_hb_ratio').addPlot('t','t').ticks('x',6).write(pdfwriter)
        plotit('rBF', 'chb').addPlot('n','t').ticks('x',6).write(pdfwriter)
        plotit('rBF', 'rBF').addPlot('n','t').ticks('x',6).write(pdfwriter)
        plotit('rBF', 'oef').addPlot('n','n').ticks('x',6).write(pdfwriter)
        plotit('rBF', 'oef').addPlot('t','t').ticks('x',6).write(pdfwriter)
        plotit('rBF', 'Y_plus_oef').addPlot('n','n').write(pdfwriter)
        plotit('rBF', 'Y_plus_oef').addPlot('t','t').write(pdfwriter)

        plotit('sat', 'sat').addPlot('n','t').write(pdfwriter)        
        
        plotit('chb', 'rBF').addPlot('t','t').write(pdfwriter)
        plotit('chb', 'rBF').addPlot('n','n').write(pdfwriter)
        plotit('mvd', 'mvd').addPlot('n','t').write(pdfwriter)
        plotit('mvd', 'sat_via_hb_ratio').addPlot('t','t').addPlot('n','n', colored = False).write(pdfwriter)
        plotit('oef', 'oef').addPlot('n','t').write(pdfwriter)
        #plotit('rBF', 'phi_vessels').addPlot('n','n').addLinearFit(legend = True).ticks('x',6).write(pdfwriter)
        #plotit('rBF', 'phi_vessels').addPlot('t','t').addLinearFit(legend = True).ticks('x',6).write(pdfwriter)
        plotit('chb', 'mtt').addPlot('t','t').write(pdfwriter)
        plotit('chb', 'mtt').addPlot('n','n').write(pdfwriter)
        plotit('mtt', 'mtt').addPlot('n','t').write(pdfwriter)
        
        plotit('Peff', 'Peff').addPlot('n','t').write(pdfwriter)
        plotit('S_rho', 'S_rho').addPlot('n','t').write(pdfwriter)
        plotit('S_rho_over_rBV', 'S_rho_over_rBV').addPlot('n','t').write(pdfwriter)
        plotit('chb',         'PeffSrhoMTT').addPlot('n','n').write(pdfwriter)
        plotit('chb',         'PeffSrhoMTT').addPlot('t','t').write(pdfwriter)
        plotit('PeffSrhoMTT', 'PeffSrhoMTT').addPlot('n','t').write(pdfwriter)
        plotit('chb',         'kExpFun').addPlot('n','n').write(pdfwriter)
        plotit('chb',         'kExpFun').addPlot('t','t').write(pdfwriter)
        plotit('phi_vessels', 'Sin').addPlot('n','n').ticks('x',4).write(pdfwriter)
        plotit('phi_vessels', 'Sin').addPlot('t','t').ticks('x',4).write(pdfwriter)
        plotit('phi_vessels', 'S_rho_over_rBV').ticks('x',4).addPlot('t','t').write(pdfwriter)
        plotit('phi_vessels', 'S_rho_over_rBV').ticks('x',4).addPlot('n','n').write(pdfwriter)
        plotit('phi_vessels', 'S_rho').addPlot('n','n').ticks('x',4).write(pdfwriter)
        plotit('phi_vessels', 'S_rho').addPlot('t','t').ticks('x',4).write(pdfwriter)
        plotit('phi_vessels', 'oef').addPlot('n','n').ticks('x',4).write(pdfwriter)
        plotit('phi_vessels', 'oef').addPlot('t','t').ticks('x',4).write(pdfwriter)
        plotit('phi_vessels', 'rBF').addPlot('t','t').addPlot('n','n', colored = False).ticks('x',6).write(pdfwriter)
        plotit('phi_vessels', 'rBF').addPlot('n','n').addLinearFit(legend = True).ticks('x',6).write(pdfwriter)
        plotit('phi_vessels', 'rBF').addPlot('t','t').addLinearFit(legend = True).ticks('x',6).write(pdfwriter)
        plotit('phi_vessels', 'phi_vessels').addPlot('n','t').ticks('x',6).write(pdfwriter)
        plotit('phi_vessels', 'Peff').addPlot('n','n').ticks('x',4).write(pdfwriter)
        plotit('phi_vessels', 'Peff').addPlot('t','t').ticks('x',4).write(pdfwriter)
        plotit('phi_vessels', 'sat_via_hb_ratio').addPlot('t','t').addPlot('n','n', colored = False).write(pdfwriter)
        plotit('phi_vessels', 'scaled_rBF').addPlot('t','t').ticks('x',4).write(pdfwriter)

        # difference plots        
#        plotit('mvd', 'sat_via_hb_ratio').addPlot('n','d').write(pdfwriter)
#        plotit('chb', 'sat_via_hb_ratio').addPlot('n','d').write(pdfwriter)
#        plotit('rBF', 'sat_via_hb_ratio').addPlot('n','d').write(pdfwriter)
#        plotit('mvd', 'sat_via_hb_ratio').addPlot('t','d').write(pdfwriter)
#        plotit('chb', 'sat_via_hb_ratio').addPlot('t','d').write(pdfwriter)
#        plotit('rBF', 'sat_via_hb_ratio').addPlot('t','d').write(pdfwriter)

        
        curves = collections.defaultdict(list)
        for item in items1:
          for mangled_name, param_name in [('o2p_mmcons_m0_tum','mmcons_m0_tum')]:
            curves[mangled_name].append(item.po2group['parameters'][param_name][()])
          curves['o2p_rmax'].append(float(item.gtumor.file['parameters/vessels'].attrs['radMax']))
          curves['o2p_vesselCompressionFactor'].append(float(item.gtumor.file['parameters/vessels'].attrs['vesselCompressionFactor']))
        for k, v in curves.items():
          curves[k] = np.asarray(v)
        data1tumor.update(curves)
        
        if np.std(curves['o2p_mmcons_m0_tum']) > 0.01*np.average(curves['o2p_mmcons_m0_tum']):
          PlotCompareInitialNetworksScatter(t0, data0glob, t1, data1glob, data1tumor, items0, items1, 'o2p_mmcons_m0_tum', 'sat_via_hb_ratio').addPlot('t','t').ticks('x',5).write(pdfwriter)
        if np.std(curves['o2p_vesselCompressionFactor']) > 0.01*np.average(curves['o2p_vesselCompressionFactor']):
          PlotCompareInitialNetworksScatter(t0, data0glob, t1, data1glob, data1tumor, items0, items1, 'radius', 'sat_via_hb_ratio').addPlot('t','t').write(pdfwriter)
          PlotCompareInitialNetworksScatter(t0, data0glob, t1, data1glob, data1tumor, items0, items1, 'o2p_vesselCompressionFactor', 'sat_via_hb_ratio').addPlot('t','t').write(pdfwriter)
          PlotCompareInitialNetworksScatter(t0, data0glob, t1, data1glob, data1tumor, items0, items1, 'o2p_vesselCompressionFactor', 'chb').addPlot('t','t').write(pdfwriter)
          PlotCompareInitialNetworksScatter(t0, data0glob, t1, data1glob, data1tumor, items0, items1, 'o2p_vesselCompressionFactor', 'mvd').addPlot('t','t').write(pdfwriter)
          #PlotCompareInitialNetworksScatter(t0, data0glob, t1, data1glob, data1tumor, items0, items1, 'o2p_vesselCompressionFactor', 'sat_via_hb_ratio').addPlot('t','d').write(pdfwriter)
        if np.std(curves['o2p_rmax']) > 0.01*np.average(curves['o2p_rmax']):
          PlotCompareInitialNetworksScatter(t0, data0glob, t1, data1glob, data1tumor, items0, items1, 'o2p_rmax', 'sat_via_hb_ratio').addPlot('t','t').write(pdfwriter)
          PlotCompareInitialNetworksScatter(t0, data0glob, t1, data1glob, data1tumor, items0, items1, 'o2p_rmax', 'chb').addPlot('t','t').write(pdfwriter)
          PlotCompareInitialNetworksScatter(t0, data0glob, t1, data1glob, data1tumor, items0, items1, 'o2p_rmax', 'mvd').addPlot('t','t').write(pdfwriter)

      if items0:
        print 'exporting ascii'
        snapshotlist = [
          (items0, data0glob, data0tumor, t0),
          (items1, data1glob, data1tumor, t1),
        ]
        ExportGlobalDataAsH5(output_base_filename, snapshotlist, measurementinfo)
        
    if 0:
      #try:
      #  histogramGroupFinal   = f_measure['combinedHistogramsFinal']
      #  histogramGroupInitial = f_measure['combinedHistogramsInitial']
      #except KeyError:
      
      histogramGroupFinal   = f_measure.recreate_group('combinedHistogramsFinal')
      histogramGroupInitial = f_measure.recreate_group('combinedHistogramsInitial')        
      ComputeHistogramsOfPo2Items(dataman, items1, measurementinfo, histogramGroupFinal)
      ComputeHistogramsOfPo2Items(dataman, items0, measurementinfo, histogramGroupInitial)
      PlotHistograms(pdfwriter, histogramGroupFinal, 'tum', 'Tumor')
      PlotHistograms(pdfwriter, histogramGroupInitial, 'all', 'Initial')
        
    if 1:
      print 'making convergence data plot'
      for items in ensemble.tumor_snapshots:
        PlotConvergenceData(pdfwriter, dataman, items)

if __name__ == '__main__':
  krebsutils.set_num_threads(2)
  filenames, pattern = sys.argv[1:-1], sys.argv[-1]
  doit(filenames, pattern)
