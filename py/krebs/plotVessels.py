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
from os.path import basename, commonprefix
import time
import krebsutils
import h5py
import h5files
import numpy as np
import itertools
import extensions
import collections
import posixpath
import math

import mpl_utils
from mystruct import Struct
import myutils
from myutils import f2l

import matplotlib

from quantities import Prettyfier
import analyzeGeneral
import analyzeBloodFlow

import krebsjobs.submitVesseltreeCalibration

##----------------------------------------------------------------------
## to plot the measuremnt data within the original vessel tree files
##----------------------------------------------------------------------
class VesselData(object):
  def __init__(self):
    self.data = Struct(
      mvd_by_iter = [],
      radii_prob = [],
      lengths_by_rad = [],
      lengths_prob = [],
      num_branches_by_rad = []
    )
  def add(self, group):
    for name, datas in self.data.iteritems():
      ds = group[name]
      datas.append(np.asarray(ds))
  def getAvgCurve(self, name):
    return np.average(self.data[name], axis=0)
  def __getitem__(self, name):
    return self.data[name]


def plot_topological_stats_avg(data, output_writer):
  fig = matplotlib.figure.Figure()
  axes = [fig.add_subplot(x) for x in range(321,327)]
  #fig, axes = pyplot.subplots(3, 2)
  #axes = axes.ravel()

  def plot_curve(ax, name, kwargs_subplot):
    ax.set(**kwargs_subplot)
    d = data.getAvgCurve(name)
    ax.plot(d[0,:,0], d[1,:,0])

  def plot_many(ax, name, kwargs_subplot):
    ax.set(**kwargs_subplot)
    dd = data[name]
    for d in dd:
      ax.plot(d[0,:,0], d[1,:,0])

  plot_many(axes[0],'mvd_by_iter', dict(ylabel='mvd', xlabel='iteration', title='MVD evolution over time'))
  plot_curve(axes[1],'radii_prob', dict(ylabel='prob', xlabel='rad', yscale='log', title='radii distribution'))
  plot_curve(axes[2],'lengths_by_rad', dict(ylabel='length', xlabel='rad', title='branch length vs radius'))
  plot_curve(axes[3],'lengths_prob', dict(ylabel='prob', xlabel='length', yscale='log', title='branch length distribution'))
  plot_curve(axes[4],'num_branches_by_rad', dict(ylabel='count', xlabel='rad', title='#branches vs radius'))
  del plot_many
  del plot_curve
  output_writer.savefig(fig)



##----------------------------------------------------------------------
## plot flow and other data of any vessel network
##----------------------------------------------------------------------


def negated_arterial_radii(radii, flags):
  indices = np.nonzero(np.bitwise_and(flags, krebsutils.ARTERY))
  radii = radii.copy()
  radii[indices] *= -1
  indices = np.nonzero(np.bitwise_and(flags, krebsutils.CAPILLARY))
  radii[indices] *= np.sign(np.random.random_sample(indices[0].shape)-0.5)
  return radii


def plot_flow_data(data, pdfpages):
  radii = negated_arterial_radii(data.edges['radius'], data.edges['flags'])
  weights = data.edges['length']
  weights /= np.sum(weights)
  radii_bins = np.linspace(np.amin(radii)-1, np.amax(radii)+1, 50)

  def plot_vs_r(name, kwargs_subplot):
    val_avg, val_std = scatter_histogram(radii, data.edges[name], radii_bins, weights)
    fig = pyplot.figure()
    p = fig.add_subplot(111, **kwargs_subplot)
    p.errorbar(np.average((radii_bins[:-1],radii_bins[1:]),axis=0), val_avg, yerr=val_std, fmt='x', elinewidth=1., capsize=1.)
    #p.bar(radii_bins[:-1], val_avg, width=radii_bins[1:]-radii_bins[:-1] )
    pdfpages.savefig(fig)

  plot_vs_r('pressure', dict(ylabel='p', xlabel='r', title='pressure vs radius'))
  plot_vs_r('flow', dict(ylabel='q', xlabel='r', title='flow rate vs radius'))
  plot_vs_r('shearforce', dict(ylabel='f', xlabel='r', title='shear force vs radius'))
  del plot_vs_r

  press = data.edges['pressure']
  press_bins = np.linspace(np.amin(press)-0.1, np.amax(press)+0.1, 50)
  val_avg, val_std = scatter_histogram(press, data.edges['shearforce'], press_bins, weights)

  fig = matplotlib.figure()
  p = fig.add_subplot(111, ylabel='f', xlabel='p', title='shear force vs pressure')
  p.errorbar(np.average((press_bins[:-1],press_bins[1:]),axis=0), val_avg, yerr=val_std, fmt='x', elinewidth=1., capsize=1.)
  #p.bar(radii_bins[:-1], val_avg, width=radii_bins[1:]-radii_bins[:-1] )
  pdfpages.savefig(fig, postfix='_vesselflowdata')


def generate_samples(graph, name, association, scale):
  DATA_LINEAR = krebsutils.VesselSamplingFlags.DATA_LINEAR
  DATA_CONST = krebsutils.VesselSamplingFlags.DATA_CONST
  DATA_PER_NODE = krebsutils.VesselSamplingFlags.DATA_PER_NODE
  if not len(graph.edgelist):
    return np.asarray([], dtype=float)
  if association == 'edges':
    data = krebsutils.edge_to_node_property(int(np.amax(graph.edgelist)+1), graph.edgelist, graph.edges[name], 'avg')
  else:
    data = graph.nodes[name]
#  return data
  return krebsutils.sample_edges(graph.nodes['position'], graph.edgelist, data, scale, DATA_LINEAR | DATA_PER_NODE)


def getMultiScatter(scale, vesselgroups):
  from collections import defaultdict as ddict
  import math
  datas = ddict(list)
  for f in vesselgroups:
    #ld = krebsutils.read_lattice_data_from_hdf(krebsutils.find_lattice_group_(f['vessels']))

    graph = analyzeGeneral.read_vessels_data(f, ['position', 'flags', 'radius', 'pressure', 'flow', 'shearforce'])

#    dd = krebsutils.calc_vessel_hydrodynamics(f['vessels'])
#    graph.nodes['pressure'] = dd[0]
#    graph.edges['flow'] = dd[1]
#    graph.edges['shearforce'] = dd[2]

    graph = graph.get_filtered(
      edge_indices = np.nonzero(graph.edges['flags'] & krebsutils.CIRCULATED)[0]
    )
    flags = graph.edges['flags']
    indices = {}
    for n, f in zip(['art','vein','capi'], [krebsutils.ARTERY,krebsutils.VEIN,krebsutils.CAPILLARY]):
      indices[n] = np.where(np.logical_and(flags & f, 0==(flags & krebsutils.WITHIN_TUMOR)))[0]
    indices['tum'] = np.where(flags & krebsutils.WITHIN_TUMOR)[0]

    print '%i edges, sample length %f' % (len(graph.edgelist), scale)

    #print graph.nodes['position']

    mydata = {}
    for subgroup in ['art', 'vein', 'capi']:
      subgraph = graph.get_filtered(edge_indices=indices[subgroup])
      mydata[('pressure', subgroup)] = generate_samples(subgraph, 'pressure', 'nodes', scale)
      mydata[('flow', subgroup)] = q = generate_samples(subgraph, 'flow', 'edges', scale)
      mydata[('shearforce', subgroup)] = generate_samples(subgraph, 'shearforce', 'edges', scale)
      mydata[('radius', subgroup)] = r = generate_samples(subgraph, 'radius', 'edges', scale)
      mydata[('velocity', subgroup)] = q/(r*r*math.pi)
    r = mydata[('radius', 'capi')]
    r *= np.sign(np.random.random_sample(r.shape)-0.5)
#    r = mydata[('radius', 'tum')]
#    r *= np.sign(np.random.random_sample(r.shape)-0.5)
    r = mydata[('radius', 'art')]
    r *= -1.
    for k, v in mydata.iteritems():
      datas[k].append(v)
  res = dict( (k,np.hstack(v)) for k,v in datas.iteritems() )
  return res


def plotMultiScatterBeauty(measures, pdfpages, withRadiusPressureRelation=True):
    # data is stored in dict as measureres[(dataname,side)]=dataarray
    avtypes = ['art','capi','vein','tum']
    if not any( side=='tum' for name,side in measures.iterkeys() ):
        avtypes.remove('tum')
    getindx = lambda dat: np.where(np.logical_and(dat>-80.0,dat<120.0))[0]
    r = dict((t, getindx(measures[('radius',t)])) for t in avtypes)

    n2title = {
        'flow':(r'flow rate [$\mu$m$^3$/s]'),
        'pressure':(r'pressure [kPa]'),
        'radius':(r'vessel radius [$\mu$m]'),
        'logradius':(r'$log_{10}$(vessel radius) [$\mu$m]'),
        'shearforce':(r'shear force [kPa]'),
        'velocity':(r'blood velocity [$\mu$m/s]'),
    }

    def plot(ax, xid, yid,lx=False,ly=False):
        ax.set(xlabel=(n2title['log'+xid] if lx else n2title[xid]), ylabel=n2title[yid])
        ax.set(yscale = 'log' if ly else 'linear') #xscale = 'log' if lx else 'linear',
        ax.grid(linestyle=':', linewidth=0.5)
        for side in avtypes:
            valx = measures[(xid,side)][r[side]]
            valy = measures[(yid,side)][r[side]]
            if lx:
              valx = np.sign(valx) * np.log10(np.abs(valx))
            q = ax.scatter(valx, valy, s=0.5, alpha=0.2, label=side, color='k', edgecolors='k', rasterized=True)
    fig = matplotlib.figure.Figure(figsize = (mpl_utils.a4size[0]*0.8, mpl_utils.a4size[0]*0.8))
    axes = [fig.add_subplot(x) for x in range(221,225)]
    #fig, axes = pyplot.subplots(2,2,figsize = (mpl_utils.a4size[0]*0.8, mpl_utils.a4size[0]*0.8))
    #axes = axes.ravel()
    fig.subplots_adjust(wspace=0.3, hspace=0.25)
    plot(axes[0],'radius','shearforce',ly=True)
    plot(axes[1],'radius','pressure')
    plot(axes[2],'radius','velocity', ly=True)
    plot(axes[3], 'radius', 'flow', ly=True, lx=True)
    r = np.linspace(-100., 100., 200)
    pressure = map(lambda r: krebsutils.PressureRadiusRelation(abs(r), bool(r<0)), r)
    if(withRadiusPressureRelation):
      axes[1].plot(r, pressure, label = 'Pressure Radius Relation')
    pdfpages.savefig(fig, dpi=320, postfix='_flowsccatter')
  

@myutils.UsesDataManager
def generateRadiusHistogram(dataman, vesselgroups, destination_group,destination_name, filterflags=None):
  def process(vesselgroups):
    bins = np.logspace(-1., 1.1, 50, base=10.)
    result = []
    for g in vesselgroups:
      r = dataman.obtain_data('basic_vessel_samples', 'radius', g, 30.)
      w = dataman.obtain_data('basic_vessel_samples', 'weight', g, 30.)
      f = dataman.obtain_data('basic_vessel_samples', 'flags', g, 30.)
      i = myutils.bbitwise_and(f, krebsutils.CIRCULATED)
      if filterflags is not None:
        i &= myutils.bbitwise_and(f, filterflags)
      h    = myutils.MeanValueArray.fromHistogram1d(bins, r[i], w[i])
      result.append(h)
    result = myutils.MeanValueArray.fromSummation(result)
    #ax.bar(bins[:-1], result.sum, width=(bins[1]-bins[0]))
    y = result.sum
    y /= np.sum(y)
    y /= (bins[1:]-bins[:-1])
    return bins[:-1],y
  def write(gmeasure, groupname):
    gmeasure = gmeasure.create_group(groupname)
    h,bin_edges = process(vesselgroups)
    gmeasure.create_dataset('h', data=h)
    gmeasure.create_dataset('bin_edges', data=bin_edges)
    
  def read(gmeasure,groupname):
    gmeasure=gmeasure[groupname]
    return gmeasure['h'],gmeasure['bin_edges']
  ret = myutils.hdf_data_caching(read,write,destination_group, (destination_name,),(1,))
  return ret


def PlotRadiusHistogram(ax, dataman, vesselgroups, filterflags=None):
  #bins = np.linspace(0., 200., 20)
  #bins = np.asarray([ 0., 4, 8, 16, 32, 64, 128, 256])
  bins = np.logspace(-1., 2., 50, base=10.)
  result = []
  for g in vesselgroups:
    r = dataman.obtain_data('basic_vessel_samples', 'radius', g, 30.)
    w = dataman.obtain_data('basic_vessel_samples', 'weight', g, 30.)
    f = dataman.obtain_data('basic_vessel_samples', 'flags', g, 30.)
    i = myutils.bbitwise_and(f, krebsutils.CIRCULATED)
    if filterflags is not None:
      i &= myutils.bbitwise_and(f, filterflags)
    h    = myutils.MeanValueArray.fromHistogram1d(bins, r[i], w[i])
    result.append(h)
  result = myutils.MeanValueArray.fromSummation(result)
  #ax.bar(bins[:-1], result.sum, width=(bins[1]-bins[0]))
  y = result.sum
  y /= np.sum(y)
  y /= (bins[1:]-bins[:-1])
  ax.step(bins[:-1], y, where='post')


def PlotRadiusHistogram2(dataman, vesselgroups, pdfpages):
  fig = matplotlib.figure.Figure(figsize = (mpl_utils.a4size[0]*0.8, mpl_utils.a4size[0]*0.4))
  ax = fig.add_subplot(111)
  PlotRadiusHistogram(ax, dataman, vesselgroups)
  ax.set(xlabel='r [$\mu m$]', ylabel='p', title = 'Radius Histogram')
  fig.subplots_adjust(bottom=0.2)
  pdfpages.savefig(fig, postfix='_radiushisto')


def FormatGeometricAndPerfusionData(geometric_data, perfusion_data):
  def printstuff_(sym, a, mult, unit=''):
     text.append('$%s$ = $%s \pm %s$ %s' % (sym, f2l(np.average(a)*mult, exponential=False), f2l(np.std(a)*mult, exponential=False), unit))
  text = []
  rbv, a, v, c = geometric_data[:4]
  printstuff_('phi_vessels'  , rbv, 1.)
  printstuff_('rBV_a', a, 1.)
  printstuff_('rBV_v', v, 1.)
  printstuff_('rBV_c', c, 1.)
  fv = v/(a+v)
  printstuff_('f_v', fv, 1.)
  printstuff_('mvd'  , geometric_data[4], 1e6, '$1/mm^2$')
  printstuff_('mvd_a', geometric_data[5], 1e6, '$1/mm^2$')
  printstuff_('mvd_v', geometric_data[6], 1e6, '$1/mm^2$')
  printstuff_('mvd_c', geometric_data[7], 1e6, '$1/mm^2$')
  printstuff_('cap_distance', geometric_data[8], 1, '$\mu m$')
  printstuff_('<r>', geometric_data[9], 1, '$\mu m$')
  printstuff_('fBF', perfusion_data, 60., '$ml\,Blood\,ml^{-1} min^{-1}$')
  return text


class Format_(object):
  exponents = collections.defaultdict(lambda : None)
  exponents.update({
    #'chb' : -6,  'chb_oxy' : -6, 'chb_deoxy' : -6,
    'jtv' : -3, 'mro2' : -3
  })
  multi = collections.defaultdict(lambda : 1)
  multi.update({
    'shearforce' : 1e3
  })
  
  def __call__(self, name, v):
    avg, std = np.average(v, axis=0), np.std(v, axis=0)
    #if exponent is not None:
    exponent = Format.exponents[name]
    multi    = Format.multi[name]
    if exponent:
      exponent_str = ('\,10^{%i}' % exponent) if exponent<>0 else ''
      f = math.pow(10., exponent)
      s = r'%s \pm %s%s' % (f2l(avg/f*multi, exponential=False), f2l(std/f*multi, exponential=False), exponent_str)
    else:
      s = r'%s \pm %s\,' % (f2l(avg*multi), f2l(std*multi))
    return s

Format = Format_()

def fmt_(b):
    for i,c in enumerate('xyz'):
      yield '%s = %s .. %s mm' % (c, f2l(b[0][i]*1.e-3, exponential=False), f2l(b[1][i]*1.e-3, exponential=False))

def PrintGlobalDataWithOxygen(pdfpages, po2groups, vesselgroups, f_measure, dataman):
  from analyzeBloodVolumeSimple import cylinderCollectionVolumeDensity
  import detailedo2Analysis.plotsForPaper
  import detailedo2
  
  sample_length = detailedo2Analysis.plotsForPaper.sample_length
  def cachelocation(g):
    path = posixpath.join('FileCS_'+myutils.checksum(basename(g.file.filename)), g.name.strip(posixpath.sep))
    return (f_measure, path)

  f2l = myutils.f2l

  bbox_vessels = list()
  bbox_field   = list()
  data_by_name = collections.defaultdict(list)
  O2_prop_list = ['po2','sat', 'sat_via_hb_ratio','gtv', 'jtv','mro2', 'po2_tissue', 'chb', 'chb_oxy', 'chb_deoxy', 'oef', 'e1', 'e3', 'Jin_root', 'Jout_root', 'Jout_tv', 'Jout_cons', 'sat_vein', 'sat_art', 'sat_capi']
  prop_list = ['mvd', 'mvd_a', 'mvd_v', 'mvd_c', 'rBV', 'rBV_a', 'rBV_v', 'rBV_c', 'venous_rBV_fraction', 'rBF', 'meanCapillaryDistance', 'mean_r']  
  #prop_list2 = ['shearforce', 'velocity']
  #prop_list2 = ['velocity']

  for po2group in po2groups:
    for prop in O2_prop_list:
      data = dataman.obtain_data('detailedPO2_global', prop, po2group, sample_length, cachelocation(po2group))
      data_by_name[prop].append(data)
    ld = krebsutils.read_lattice_data_from_hdf(po2group['field_ld'])
    bbox_field.append(ld.worldBox)
    gvessels, _ = detailedo2.OpenVesselAndTumorGroups(po2group)
#    for prop in prop_list2:
#      data = dataman.obtain_data('basic_vessel_global', prop, gvessels, cachelocation(gvessels))
#      data_by_name[prop].append(data)
    ld = krebsutils.read_lattice_data_from_hdf(gvessels['lattice'])
    bbox_vessels.append(ld.worldBox)
    rbv, a, v, c = cylinderCollectionVolumeDensity(gvessels)
    sa, sv, sc = data_by_name['sat_art'][-1], data_by_name['sat_vein'][-1], data_by_name['sat_capi'][-1]
    data_by_name['sat_estimated_by_acv'].append((sa*a+sv*v+sc*c)/rbv)
  O2_prop_list.append('sat_estimated_by_acv')  
 
  try:
    os.remove('initialvesseldata.h5')
  except OSError:
    pass
  with h5py.File('initialvesseldata.h5', 'w-') as f:
    for k, v in data_by_name.iteritems():
      f.create_dataset(k, data = v)

  

#  def fmt_(v, exponent=None, multi=1.):
#    avg, std = np.average(v, axis=0), np.std(v, axis=0)
#    if exponent is not None:
#      exponent_str = ('\,10^{%i}' % exponent) if exponent<>0 else ''
#      f = math.pow(10., exponent)
#      s = r'%s \pm %s%s' % (f2l(avg/f*multi, exponential=False), f2l(std/f*multi, exponential=False), exponent_str)
#    else:
#      s = r'%s \pm %s\,' % (f2l(avg*multi), f2l(std*multi))
#    return s
  
  result_string = []
  for name,v in myutils.iterate_items(data_by_name, O2_prop_list):
    result_string.append(r'$<%s>$ = $%s$%s' %
      (Prettyfier.get_sym(name), Format(name, v), Prettyfier.get_munit(name)))

#  text = FormatGeometricAndPerfusionData()

  prop_list2 = ['avg_cap_dist','mvd']
  #bbox_vessels = list()
  #bbox_field   = list()
  
  for name in prop_list2:
    data = []
    for gvessels in vesselgroups:
      data.append(dataman.obtain_data('basic_vessel_global', name, gvessels, cachelocation(gvessels)))
      #ld = krebsutils.read_lattice_data_from_hdf(vesselgroups[0]['lattice'])
      #bbox_vessels.append(ld.worldBox)
    result_string.append(r'$<%s>$ = $%s$%s' %
      (Prettyfier.get_sym(name), Format(name, data), Prettyfier.get_munit(name)))
        
  

  bbox_vessels = np.average(bbox_vessels, axis=0).reshape(3,2).transpose()
  bbox_field = np.average(bbox_field, axis=0).reshape(3,2).transpose()
  result_string += ['Vessel System Bounding Box'] + list(fmt_(bbox_vessels)) + ['FD-Grid Bounding Box'] + list(fmt_(bbox_field))

  fig, _ = mpl_utils.MakeTextPage(result_string, figsize = (mpl_utils.a4size[0], mpl_utils.a4size[0]))
  pdfpages.savefig(fig, postfix='_vesselsglobal')

  fig = matplotlib.figure.Figure(figsize = (mpl_utils.a4size[0]*0.5, mpl_utils.a4size[0]*0.5))
  ax = fig.add_axes([0.1, 0.2, 0.8, 0.75])
  for po2group in po2groups:
    iterations = detailedo2Analysis.plotsForPaper.GetConvergenceData(po2group)
    x = iterations['iteration']
    y = iterations['delta_vessM'] + iterations['delta_fieldM']
    ax.plot(x, y)
  ax.set(yscale = 'log', xlabel = 'iterations', ylabel = r'$|p(new) - p(old)|_\infty$')
  pdfpages.savefig(fig, postfix='_convergence')
  
#def GenerateHistogramData(dataman, vesselgroups, po2groups, f_measure):
#  sample_length = 30.
#  def cachelocation(g):
#    path = posixpath.join('FileCS_'+myutils.checksum(basename(g.file.filename)), g.name.strip(posixpath.sep))
#    return (f_measure, path)
#
#  result = {}
#  for name in ['phi_vessels', 'mvd']:
#    result[name] = map(lambda gvessels: dataman.obtain_data('basic_vessel_global', name, gvessels, cachelocation),
#                       vesselgroups)
#  result['rBF'] = map(getTotalPerfusion, vesselgroups)
#
#  for k, l in result.items():
#    result[k] = myutils.MeanValueArray.fromHistogram1d(binspecs[k], l)
#
#  for name in ['radius', 'shearforce', 'velocity']:
#    def make_(gvessels):
#      smpl = dataman.obtain_data('basic_vessel_samples', name, gvessels, sample_length)
#      mask = analyzeGeneral.GetMaskUncirculated(dataman, gvessels)
#      return analyzeGeneral.ComputeSampleHistogram(
#        dataman,
#        gvessels,
#        binspecs[name],
#        smpl,
#        mask,
#        sample_length)
#    result[name] = myutils.MeanValueArray.fromSummation(
#      map(make_, vesselgroups))
#  return result
def PrintGlobalData(pdfpages, vesselgroups, f_measure, dataman):
  from analyzeBloodVolumeSimple import cylinderCollectionVolumeDensity
  import detailedo2Analysis.plotsForPaper
  import detailedo2
  
  sample_length = detailedo2Analysis.plotsForPaper.sample_length
  def cachelocation(g):
    path = posixpath.join('FileCS_'+myutils.checksum(basename(g.file.filename)), g.name.strip(posixpath.sep))
    return (f_measure, path)

  f2l = myutils.f2l

  bbox_vessels = list()
  data_by_name = collections.defaultdict(list)
  #prop_list = ['mvd', 'mvd_a', 'mvd_v', 'mvd_c', 'rBV', 'rBV_a', 'rBV_v', 'rBV_c', 'venous_rBV_fraction', 'rBF', 'meanCapillaryDistance', 'mean_r']  
  #prop_list = ['rBV',  'rBF', ]  
  #prop_list = 'phi_a phi_v phi_c mvd_a mvd_v mvd_c mean_r'.split()
 
  prop_list2 = ['radius','shearforce','velocity','flow','avg_cap_dist','mvd_linedensity','phi_vessels','total_perfusion']
 
  try:
    os.remove('initialvesseldata.h5')
  except OSError:
    pass
    with h5py.File('initialvesseldata.h5', 'w-') as f:
      for k, v in data_by_name.iteritems():
        f.create_dataset(k, data = v)
  
  result_string = []

  #bins_spec   = analyzeGeneral.BinsSpecRange(100., 600., 100.)
  def suggest_bins_from_world(ld):
    avg_size=np.average(ld.GetWorldSize())
    center = avg_size/2.
    half_center = center/2.
    #lower=center-0.3*avg_size
    #upper=center+0.3*avg_size
    #stepsize = 0.1*center
    #bins_spec   = analyzeGeneral.BinsSpecRange(center-stepsize, center, stepsize)
    bins_spec   = analyzeGeneral.BinsSpecRange(half_center, center, 200.)
    return bins_spec
  for name in prop_list2:
    data = []
    for gvessels in vesselgroups:
      data.append(dataman.obtain_data('basic_vessel_global', name, gvessels, cachelocation(gvessels)))
      ld_vessels = krebsutils.read_lattice_data_from_hdf(gvessels['lattice'])
      bbox_vessels.append(ld_vessels.worldBox)
    result_string.append(r'$<%s>$ = $%s$%s' %
      (Prettyfier.get_sym(name), Format(name, data), Prettyfier.get_munit(name)))
  mvd_exp=[]
  for gvessels in vesselgroups:
    ld = krebsutils.read_lattice_data_from_hdf(gvessels.parent['field_ld'])
    mvd_sampling_results, mvd_bins = dataman.obtain_data('sphere_vessel_density',  gvessels, None, suggest_bins_from_world(ld), 'radial', ld, cachelocation(gvessels))
    #print(mvd_sampling_results)    
    mvd_exp.append(np.mean(np.asarray(mvd_sampling_results)*1e6))
  #mvd_exp=np.asarray(mvd_exp)
  #print(mvd_exp)
  result_string.append(r'$<%s>$ = $%s$%s' %
    (Prettyfier.get_sym('mvd_exp'), Format('mvd_exp', mvd_exp), Prettyfier.get_munit('mvd_exp')))
      
  ld = krebsutils.read_lattice_data_from_hdf(vesselgroups[0]['lattice'])
  bbox_vessels.append(ld.worldBox)

  bbox_vessels = np.average(bbox_vessels, axis=0).reshape(3,2).transpose()
  result_string += ['Vessel System Bounding Box'] + list(fmt_(bbox_vessels))


  fig, _ = mpl_utils.MakeTextPage(result_string, figsize = (mpl_utils.a4size[0], mpl_utils.a4size[0]))
  pdfpages.savefig(fig, postfix='_vesselsglobal')

#  fig = matplotlib.figure.Figure(figsize = (mpl_utils.a4size[0]*0.5, mpl_utils.a4size[0]*0.5))
#  ax = fig.add_axes([0.1, 0.2, 0.8, 0.75])
#  for po2group in po2groups:
#    iterations = detailedo2Analysis.plotsForPaper.GetConvergenceData(po2group)
#    x = iterations['iteration']
#    y = iterations['delta_vessM'] + iterations['delta_fieldM']
#    ax.plot(x, y)
#  ax.set(yscale = 'log', xlabel = 'iterations', ylabel = r'$|p(new) - p(old)|_\infty$')
#  pdfpages.savefig(fig, postfix='_convergence')

def plot_geometric_stuff_on_RC(dataman, f_measure, filenames, options, pdfpages):

  destination_group = f_measure.require_group('adaption/geometry/rBV')  
  
  counter = 0

  if(options.two):
    typelist = 'typeD- typeE- typeG- typeH-'
  else:
    typelist = 'typeA- typeB- typeC- typeF- typeI-'
  if(options.all_types):
    typelist = 'typeA- typeB- typeC- typeD- typeE- typeF- typeG- typeH- typeI-'
  
    
  if(options.single):
    typelist = 'typeF- '
  if 0:#for developing
    typelist = 'typeA- typeB-'
    
  ### make it right, filter
  typelist = 'typeA- typeB- typeC- typeD- typeE- typeF- typeG- typeH- typeI-'
  #filteredFiles = filter( lambda fn: t in fn,filenames)
  reduceTypeList = []
  for fn in filenames:
    for t in typelist.split():
      if t in fn and not t in reduceTypeList:
        reduceTypeList.append(t)
  rBVs_by_type = dict()
  avg_cap_dists_by_type = dict()
  mean_r_by_type = dict()
  perfusions_by_type = dict()
  big_dict = {}
  def find_enlargement_factors(filenames):
    factors = []
    files = [h5files.open(fn, 'r+') for fn in filenames]
    for f in files:
      if 'enlargeFactor' in f.attrs.keys():
        if not f.attrs.get('enlargeFactor') in factors:
          factors.append(f.attrs.get('enlargeFactor'))
    return factors
  theFactors = find_enlargement_factors(filenames)
  theFactors[:0]=[0]
  endities = 'rBVs avg_cap_dist mean_rs perfusion_s'
  if len(theFactors)>0:
    for fac in theFactors: #,reduceTypeList, endities.split()):
      #big_dict[str(fac)][t]=dict('rBVs'=[],'avg_cap_dist'=[],'mean_rs'=[],'perfusion_s'=[])
      big_dict[str(fac)]= {}
      for t in reduceTypeList:      
        big_dict[str(fac)][t]={}
        for endity in endities.split():
          big_dict[str(fac)][t][endity] = {}
          big_dict[str(fac)][t][endity]=[]#,'avg_cap_dist'=[],'mean_rs'=[],'perfusion_s'=[])
      
  for t in reduceTypeList:
    print('data for type: %s' % t)
    filteredFiles = filter( lambda fn: t in fn,filenames) 
    files = [h5files.open(fn, 'r+') for fn in filteredFiles]

    if(len(files)>0):#that means no file of dedicated type is left after filter
      if len(theFactors)==0:
        rBVs = []
        avg_cap_dists = []
        mean_rs = []
        perfusion_s = []
        for f in files:
          #rBV_of_file = analyzeGeneral.generate_rBV_of_group(dataman, destination_group, f)
          g_vessels = f['vessels']
          aCache = cachelocation(g_vessels,f_measure)
          rBV_of_file = dataman.obtain_data('basic_vessel_global','phi_vessels',g_vessels,aCache)        
          avg_cap_dist_of_file= dataman.obtain_data('basic_vessel_global','avg_cap_dist',g_vessels,aCache)         
          mean_r_of_file= dataman.obtain_data('geometric_data','geometric_detailed',g_vessels,aCache) 
          perfusion_of_file = dataman.obtain_data('basic_vessel_global','total_perfusion',g_vessels,aCache)
          perfusion_of_file =  np.asarray(perfusion_of_file)*60 #to minutes       
          perfusion_of_file = perfusion_of_file *100 #to per 100cm^3        
          rBVs.append(rBV_of_file)
          avg_cap_dists.append(avg_cap_dist_of_file)
          mean_rs.append(mean_r_of_file[6])
          perfusion_s.append(perfusion_of_file)
        rBVs_by_type[t] = rBVs
        avg_cap_dists_by_type[t] = avg_cap_dists
        mean_r_by_type[t] = mean_rs
        perfusions_by_type[t] = perfusion_s
        #data = generate_murray_data_plot3(dataman, files, destination_group, t[:-1])      
      else:
#        rBVs = []
#        avg_cap_dists = []
#        mean_rs = []
#        perfusion_s = []
        for f in files:
          #rBV_of_file = analyzeGeneral.generate_rBV_of_group(dataman, destination_group, f)
          g_vessels = f['vessels']
          if 'enlargeFactor' in f.attrs.keys():
            current_fac = f.attrs.get('enlargeFactor')
          else:
            current_fac = 0
          aCache = cachelocation(g_vessels,f_measure)
          rBV_of_file = dataman.obtain_data('basic_vessel_global','phi_vessels',g_vessels,aCache)        
          avg_cap_dist_of_file= dataman.obtain_data('basic_vessel_global','avg_cap_dist',g_vessels,aCache)         
          mean_r_of_file= dataman.obtain_data('geometric_data','geometric_detailed',g_vessels,aCache) 
          perfusion_of_file = dataman.obtain_data('basic_vessel_global','total_perfusion',g_vessels,aCache)
          perfusion_of_file =  np.asarray(perfusion_of_file)*60 #to minutes       
          perfusion_of_file = perfusion_of_file *100 #to per 100cm^3
          datas=[rBV_of_file,avg_cap_dist_of_file,mean_r_of_file[6],perfusion_of_file]
#          rBVs.append(rBV_of_file)
#          avg_cap_dists.append(avg_cap_dist_of_file)
#          mean_rs.append(mean_r_of_file[6])
#          perfusion_s.append(perfusion_of_file)
          for (endity,data) in zip(endities.split(),datas):
            a=big_dict[str(current_fac)][t][endity]
            a=list(a)
            a.append(data)
            big_dict[str(current_fac)][t][endity]=a
        #rBVs_by_type[t] = rBVs
        #avg_cap_dists_by_type[t] = avg_cap_dists
        #mean_r_by_type[t] = mean_rs
        #perfusions_by_type[t] = perfusion_s
    else:
      continue
    
    counter= counter +1    

  print("coutner: %i"%counter)

  def type_to_label(typelist):
    type_to_label_dict={"typeA-": "RC1","typeB-": "RC2","typeC-": "RC3","typeD-": "RC4","typeE-": "RC5","typeF-": "RC6","typeG-": "RC7","typeH-": "RC8","typeI-": "RC9"}
    label = []
    for t in typelist.split():
      label.append(type_to_label_dict[t])
    return label
  
  '''matplotlib rBV'''
  for endity in endities.split():
    fig = matplotlib.figure.Figure()
    ax = fig.add_subplot(111)
    
    for (k,current_fac) in enumerate(theFactors):
      data_vector_rBV = []
      for t in reduceTypeList:
        data_vector_rBV.append(np.asarray(big_dict[str(current_fac)][t][endity]))
  #print(data_vector_rBV)
  #    'rBVs avg_cap_dist mean_rs perfusion_s'
      #ax.boxplot(data_vector_rBV,positions=np.array(xrange(len(data_vector_rBV)))*2.0-0.4+k, sym='', widths=0.6)
      ax.boxplot(data_vector_rBV,positions=np.array(xrange(len(data_vector_rBV)))+k*0.25, sym='', widths=0.2)
      ax.set_xlim([-2,9])
    #ax.axhspan(0,.10,facecolor='b',alpha=0.5)
    ax.set_xticklabels(type_to_label(typelist))
    ax.set_xticks(np.array(xrange(len(data_vector_rBV)))+0.35)
    ax.set_xlabel('Configuration')
    ax.set_ylabel(endity)
    #ax.set_title('Vessel Volume Fractions by RC')
  
    #fig.tight_layout()
    #pdfpages.savefig(fig, bbox_extra_artists=(lgd,), bbox_inches='tight')
    #pdfpages.savefig(fig, bbox_inches='tight')
    pdfpages.savefig(fig)
  
  if 0:
    '''matplotlib avg cap dist'''
    fig2 = matplotlib.figure.Figure()
    ax2 = fig2.add_subplot(111)
    data_vector_cap = []
    for t in typelist.split():
      data_vector_cap.append(np.asarray(avg_cap_dists_by_type[t]))
    #print(data_vector_cap)
    ax2.boxplot(data_vector_cap)
    
    ax2.set_xticklabels(type_to_label(typelist))
    ax2.set_xlabel('Configuration')
    ax2.set_ylabel('avg cap dist (by MVD)')
    ax2.set_title('Average capillary distance by RC')
  
    fig2.tight_layout()
    pdfpages.savefig(fig2, bbox_inches='tight')
    
    '''matplotlib mean r '''
    fig3 = matplotlib.figure.Figure()
    ax3 = fig3.add_subplot(111)
    data_vector_meanr = []
    for t in typelist.split():
      data_vector_meanr.append(np.asarray(mean_r_by_type[t]))
    #print(data_vector_meanr)
    ax3.boxplot(data_vector_meanr)
    
    ax3.set_xticklabels(type_to_label(typelist))
    ax3.set_xlabel('Configuration')
    ax3.set_ylabel('<r>')
    ax3.set_title('Mean radius by RC')
  
    fig3.tight_layout()
    pdfpages.savefig(fig3, bbox_inches='tight')
    
    '''matplotlib perfusion '''
    fig4 = matplotlib.figure.Figure()
    ax4 = fig4.add_subplot(111)
    data_perfusion = []
    for t in typelist.split():
      data_perfusion.append(np.asarray(perfusions_by_type[t]))
    #print(data_vector_meanr)
    ax4.boxplot(data_perfusion)
    
    ax4.set_xticklabels(type_to_label(typelist))
    ax4.set_xlabel('Configuration')
    ax4.set_ylabel(r'rBF in $\frac{ml}{min 100cm^3}$')
    ax4.set_title('Perfusion by RC')
  
    fig4.tight_layout()
    pdfpages.savefig(fig4, bbox_inches='tight')
def FormatParameters(root):
  from quantities import Prettyfier
  parametergroup = root['parameters']
  parameters = dict(parametergroup.attrs)
  parameters = sorted(parameters.items(), key = lambda (k, v): k)
  result = ['++++ Vessel Network Creation Parameters ++++']
  for k, v in parameters:
    k = 'vgp_'+k
    try:
      name, sym, unit = Prettyfier.get_label(k), Prettyfier.get_msym(k), Prettyfier.get_munit(k)
    except KeyError:
      continue
    if isinstance(v, float):
      v = f2l(v)
    else:
      v = str(v)
    #name = name.ljust(60)
    #print len(name)
    result.append(
      '%s:\n%s = $%s$ %s' % (name, sym, v, unit))
  return result

''' create a cache '''
def cachelocation(g,f_measure):
  path = posixpath.join('FileCS_'+myutils.checksum(basename(g.file.filename)), g.name.strip(posixpath.sep))
  return (f_measure, path)
def DoIt(filenames, pattern, with_o2):
  fn_measure = basename(commonprefix(filenames))
  fn_measure = myutils.strip_from_end(fn_measure, '.h5')
  fn_measure = myutils.strip_from_end(fn_measure, '-type')
  if with_o2:
    fn_measure = myutils.strip_from_end(fn_measure, '_detailedpo2')

  
  files = [h5files.open(fn, 'a') for fn in filenames]
  f_measure = h5files.open('plotVessels_chache.h5', 'a', search = False)
  groups = list(itertools.chain.from_iterable(myutils.walkh5(f, pattern, return_h5objects=True) for f in files))
  if len(groups)<=0:
    print 'no matching groups in hdf file(s)'
    sys.exit(0)

  if with_o2:
    name = posixpath.commonprefix(map(lambda g: g.name, groups))
    name = myutils.strip_from_start(name, '/po2/vessels').replace('/', '-')
    fn_measure += name

  with mpl_utils.PdfWriter(fn_measure+'.pdf') as pdfpages:
    rc = matplotlib.rc
    rc('font', size = 8.)
    rc('axes', titlesize = 10., labelsize = 8.)

    if with_o2:
      import detailedo2Analysis as o2analysis
      import detailedo2Analysis.plotsForPaper
      import detailedo2
      dataman = myutils.DataManager(20, [ o2analysis.DataDetailedPO2(), analyzeGeneral.DataTumorTissueSingle(), analyzeGeneral.DataDistanceFromCenter(), analyzeGeneral.DataBasicVessel(), analyzeGeneral.DataVesselSamples(), analyzeGeneral.DataVesselRadial(), analyzeGeneral.DataVesselGlobal()])

      vesselgroups = list(detailedo2.OpenVesselAndTumorGroups(g)[0] for g in groups)
      #original_vesselgroups = list(h5files.openLink(g, 'SOURCE') for g in vesselgroups)
      if 1:
        PrintGlobalDataWithOxygen(pdfpages, groups, vesselgroups, f_measure, dataman)

        '''FormatParameters makes the network creation parameters
            that does not work, if we have an o2 file'''
        #text = FormatParameters(original_vesselgroups[0].file)
        text = [' ']
        text += detailedo2Analysis.plotsForPaper.FormatParameters(groups[0])
        fig, _ = mpl_utils.MakeTextPage(text,figsize = (mpl_utils.a4size[0], mpl_utils.a4size[0]))
        pdfpages.savefig(fig, postfix='_vesselsparams')
      if 1:
        res = getMultiScatter(300. * len(filenames), vesselgroups)
        plotMultiScatterBeauty(res, pdfpages)

    else:
      dataman = myutils.DataManager(20, [analyzeGeneral.DataTumorTissueSingle(), 
                                          analyzeGeneral.DataVesselRadial(), 
                                          analyzeGeneral.DataDistanceFromCenter(),
                                          analyzeBloodFlow.DataTumorBloodFlow(),
                                          analyzeGeneral.DataBasicVessel(),
                                          analyzeGeneral.DataVesselSamples(),
                                          analyzeGeneral.DataVesselGlobal()
                                          ])
      #dataman = myutils.DataManager(20, [ analyzeGeneral.DataBasicVessel(), analyzeGeneral.DataVesselSamples(), analyzeGeneral.DataVesselGlobal()])
      vesselgroups = groups
    
      if 1:
        res = getMultiScatter(300. * len(filenames), vesselgroups)
        plotMultiScatterBeauty(res, pdfpages)
      if 1:
        PlotRadiusHistogram2(dataman, vesselgroups, pdfpages)
            
      if 0:
        text = FormatGeometricAndPerfusionData()

        prop_list2 = ['shearforce', 'velocity', 'avg_cap_dist',]
        bbox_vessels = list()
        bbox_field   = list()
        
        for name in prop_list2:
          data = []
          for gvessels in vesselgroups:
            data.append(dataman.obtain_data('basic_vessel_global', name, gvessels, cachelocation(gvessels)))
            ld = krebsutils.read_lattice_data_from_hdf(vesselgroups[0]['lattice'])
            bbox_vessels.append(ld.worldBox)
          text.append(r'$<%s>$ = $%s$%s' %
            (Prettyfier.get_sym(name), Format(name, data), Prettyfier.get_munit(name)))
        
        
        
        #bbox_field = np.average(bbox_field, axis=0).reshape(3,2).transpose()
        bbox_vessels = np.average(bbox_vessels, axis=0).reshape(3,2).transpose()        
        text += ['Vessel System Bounding Box'] + list(fmt_(bbox_vessels)) + ['FD-Grid Bounding Box'] # + list(fmt_(bbox_field))
        fig, _ = mpl_utils.MakeTextPage(text,figsize = (mpl_utils.a4size[0]*0.8, mpl_utils.a4size[0]*0.8))
        pdfpages.savefig(fig, postfix='_vesselsglobal')

        text = FormatParameters(vesselgroups[0].file)
        fig, _ = mpl_utils.MakeTextPage(text,figsize = (mpl_utils.a4size[0]*0.8, mpl_utils.a4size[0]*0.8))
        pdfpages.savefig(fig, postfix='_vesselsparams')
      if 0 and all(map(lambda g: 'data' in g.parent, vesselgroups)):
        data = VesselData()
        for g in vesselgroups:
          data.add(g.parent['data'])
        plot_topological_stats_avg(data, pdfpages)
      if 0: #reproduce swine
        plot_geometric_stuff_on_RC(dataman, f_measure, filenames, options, pdfpages)
      if 1:
        PrintGlobalData(pdfpages, vesselgroups, f_measure, dataman)


if __name__ == "__main__":
  import argparse
  parser = argparse.ArgumentParser(description='Plot/ Analyze infos about vessel network.')
  parser.add_argument('vesselFileNames', nargs='*', type=argparse.FileType('r'), default=sys.stdin, help='Vessel file to calculate')  
  parser.add_argument('grp_pattern',help='Where to find the vessel group in the file')    
  parser.add_argument("-O","--with-o2", dest="with_o2", help="look at detailed o2 data", default=False, action="store_true")
  parser.add_argument("-T","--only_two_root", dest="two", help="flag to change the considered types", default=False, action="store_true")  
  parser.add_argument("-a","--with_all_types", dest="all_types", help="take all types",default=False, action="store_true")  
  parser.add_argument("-s","--singel_type", dest="single", help="", default=False, action="store_true")    
  goodArguments, otherArguments = parser.parse_known_args()

  try:
    dirs = set()
    for fn in goodArguments.vesselFileNames:
      if not os.path.isfile(fn.name):
        raise AssertionError('The file %s is not present!'%fn)
      with h5py.File(fn.name, 'r') as f:
        d = myutils.walkh5(f, goodArguments.grp_pattern)
        if not len(d)>0:
          raise AssertionError('pattern "%s" not found in "%s"!' % (grp_pattern, fn))
        else:
          dirs = set.union(dirs,d)
  except Exception, e:
    print e.message
    sys.exit(-1)
  
  print('Resolved groups: %s' % ','.join(dirs))
  #create filename due to former standards
  filenames=[]
  for fn in goodArguments.vesselFileNames:
    filenames.append(fn.name)
    ''' hack to center'''
    with h5py.File(fn.name, 'r+') as f:
      # centering is needed because quantities are analyzed in dependence on the distance from the system origin!!
      krebsjobs.submitVesseltreeCalibration.CenterTheLattice(f, 'field_ld')
      krebsjobs.submitVesseltreeCalibration.CenterTheLattice(f, 'vessels/lattice')
      f.flush()
      krebsjobs.submitVesseltreeCalibration.ObtainDataOfVesselFile(f)
  DoIt(filenames, goodArguments.grp_pattern, goodArguments.with_o2)
