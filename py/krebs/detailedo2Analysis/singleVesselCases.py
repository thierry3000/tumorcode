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
  sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','..'))

import os,sys
import h5py
import h5files
import numpy as np
import extensions # for hdf5 support in np.asarray
import krebsutils
import myutils
import posixpath
import math
from pprint import pprint
from copy import copy, deepcopy
from collections import namedtuple

import matplotlib
import matplotlib.cm
import matplotlib.pyplot as pyplot
import matplotlib.ticker
import mpl_utils

import scipy.optimize

from myutils import sanitize_posixpath, f2l

from krebs.analyzeGeneral import  generate_samples, DataBasicVessel, DataVesselSamples, DataVesselGlobal
from krebs.quantities import Prettyfier
from krebs import plotBulkTissue
from krebs import detailedo2
from krebs import vesselgenerator
from krebs import detailedo2Analysis
from krebs.detailedo2Analysis import DataDetailedPO2
from krebs.detailedo2Analysis import singleVesselParameterSets
from krebs.detailedo2Analysis import plotsForPaper

sample_length = 2.

def niceColors(n):
  return map(matplotlib.cm.jet, np.linspace(0., 1., n))

def niceMarker(i, color = 'k'):
  return mpl_utils.styleMarker('x+^dso><*'[i], color = color, linewidth = 0.5)


def CenterLattice(g):
  '''needs h5 goup g, which is deleted and rewritten where the new lattice is centered
  '''
  print '------',g.name+'--------'
  ld = krebsutils.read_lattice_data_from_hdf(g)
  #meanwhile worldBox is 3D
  thisBox=ld.worldBox
  bmin, bmax = thisBox[0:2]
  offset = -0.5*(bmin+bmax) # which is the negative center
  #offset[0] = 0 # don't move in x
  orig = ld.GetOriginPosition()
  orig += offset
  ld.SetOriginPosition(orig)
  parent = g.parent
  name   = str(g.name)
  del g
  del parent[name]
  #krebsutils.write_lattice_data_to_hdf(parent, name, ld)
  fn=str(parent.file.filename)
  #path=str(vesselgroup.name)
  krebsutils.write_lattice_data_to_hdf_by_filename(fn, name, ld)
  return offset


def GenerateSingleCapillaryWPo2(dataman, f, config_name, config_version, (bfparams, po2params, params)):
  version = myutils.checksum(1,config_version, params, bfparams, po2params)

  def write(gmeasure, name):
    gmeasure = gmeasure.create_group(name)
    # --- create single capillary vessel system ---
    p = params['pgrad']*params['size'][0]*params['scale']
    krebsutils.vesselgen_generate_single(gmeasure, params['size'], params['direction_mode'], params['scale'], params['ld_type'], params['r'], 0.1, p, params['size'][0]/10, bfparams)
    vesselgroup = gmeasure['vessels']
    if params['direction_mode'] == 0:
      CenterLattice(vesselgroup['lattice'])
    vesselgroup.attrs['VERSION'] = version
    # --- setup po2 stuff ---
    po2group = gmeasure.create_group('po2')
    po2group.attrs['SOURCE_FILE'] = f.filename
    po2group.attrs['SOURCE_PATH'] = str(vesselgroup.name)
    po2group['SOURCE']            = h5py.SoftLink(vesselgroup)
    # ---- write all parameters ---
    pgroup = gmeasure['parameters']
    for k,v in params.iteritems():
      pgroup.attrs[k] = v
    pgroup2 = pgroup.create_group('calcflow')
    for k,v in bfparams.iteritems():
      pgroup2.attrs[k] = v
    # ----
    detailedo2.computePO2_(po2group, vesselgroup, None, po2params)
    # generated samples here as well for completeness and ease of debugging. put the sample data under gmeasure/name/samples.
    fluxes = dataman.obtain_data('detailedPO2_total_fluxes', '', po2group, sample_length, 1, (gmeasure.parent, name))

  def read(gmeasure, name):
    return gmeasure[name]


  return myutils.hdf_data_caching(read, write, f,
                                  (config_name,),
                                  (version,))



def GenerateSingleCapillarySamples(dataman, po2group, cachelocation, properties = ['po2','sat', 'extpo2']):
  vesselgroup, _ = detailedo2.OpenVesselAndTumorGroups(po2group)
  smpl = dict(
    weight = dataman.obtain_data('basic_vessel_samples', 'weight', vesselgroup, sample_length),
    pos = dataman.obtain_data('basic_vessel_samples', 'position', vesselgroup, sample_length),
  )
  for prop in properties:
    smpl[prop] = dataman.obtain_data('detailedPO2_samples', prop, po2group, sample_length, 1, cachelocation)
  # sort by x-position
  x = smpl['pos'][:,0]
  ordering = np.argsort(x)
  pos = smpl['pos'][ordering,...]
  x0 = 0
  dv = (pos[-1] - x0)
  dv /= np.linalg.norm(dv)
  x = np.asarray(map(lambda x_: np.dot(dv, x_-x0), pos))
  for k, v in smpl.items():
    smpl[k] = v[ordering]
  smpl['x'] = x
  return smpl



def GenerateGlobalMeasurementOutputStrings(dataman, po2group, cachelocation):
  prop_list = ['gtv', 'jtv','mro2', 'po2_tissue', 'oef', 'e1', 'e3', 'Jin_root', 'Jout_root', 'Jout_tv', 'Jout_cons']

  Get = lambda prop: dataman.obtain_data('detailedPO2_global', prop, po2group, sample_length, cachelocation)
  data = map(Get, prop_list)
  Fmt = lambda (name, val): r'$%s = %s$ %s' % (Prettyfier.get_sym(name), f2l(val), Prettyfier.get_munit_legacy(name))
  result_string = map(Fmt, zip(prop_list, data))
  return result_string

def mm_(s):
  '''math mode'''
  return '$'+s+'$' if s else s

def fmt(name, sym, value, unit):
  value = f2l(value)
  if unit:
    value = '%s %s' % (mm_(value), mm_(unit))
  return '%s %s = %s' % (name, mm_(sym), value)



def GenerateParameterOutputStrings(dataman, group, cachelocation):
  params = myutils.hdf_read_dict_hierarchy(group['po2']['parameters'])
  inletPo2 = params['po2init_r0']

  for name in 'po2init_r0 po2init_dr po2init_cutoff mmcons_k_tum mmcons_k_necro mmcons_m0_tum mmcons_m0_necro'.split():
    del params[name]
  text = plotsForPaper.FormatParameters_(params)

  samples = GenerateSingleCapillarySamples(dataman, group['po2'], cachelocation)
  inletPressure = np.amax(dataman.obtain_data('vessel_graph_property', group['vessels'], 'nodes', 'pressure')[0])
  capillaryLength = np.linalg.norm((samples['pos'][0] - samples['pos'][-1]))
  velocity = dataman.obtain_data('basic_vessel_global', 'velocity', group['vessels'], cachelocation)
  flow = dataman.obtain_data('basic_vessel_global', 'flow', group['vessels'], cachelocation)

  params = myutils.hdf_read_dict_hierarchy_attr(group['parameters'])

  parameterlist = [
    ('Blood Pressure @ Inlet', '', inletPressure, 'kPa'),
    ('Blood Pressure Gradient', '', params['pgrad'], 'kPa\,/\,\mu m'),
    ('Capillary Length', '', capillaryLength, '\mu m'),
    (Prettyfier.get_label('velocity'), '', velocity, Prettyfier.get_unit_legacy('velocity')),
    (Prettyfier.get_label('flow'), '', flow, Prettyfier.get_unit_legacy('flow')),
    ('Discharge Hematocrit', '', params['calcflow']['inletHematocrit'], ''),
    ('Capillary Radius', '', params['r'], '\mu m'),
    ('Oxygen PO2 @ Inlet', '', inletPo2, 'mmHg'),
  ]
  text += ['------------------------']
  text += map(lambda x: fmt(*x), parameterlist)
  return text


class Plotty(object):
  def __init__(self, dataman, group, label):
    self.label   = label
    self.dataman = dataman
    self.group   = group

    self.cachelocation = (group.file, group.name)
    self.po2group, self.vesselgroup = group['po2'], group['vessels']

    #self.po2vessels, self.po2fieldld, self.po2field, _ = dataman.obtain_data('detailedPO2', self.po2group)
    self.po2fieldld = krebsutils.read_lattice_data_from_hdf(self.po2group['field_ld'])
    self.po2vessels = self.po2group['po2vessels']
    self.po2field   = self.po2group['po2field']
    self.worldbb = self.po2fieldld.worldBox
    self.fluxes = dataman.obtain_data('detailedPO2_total_fluxes', '', self.po2group, sample_length, 1, self.cachelocation)
    self.samples = GenerateSingleCapillarySamples(dataman, self.po2group, self.cachelocation)
    self.samplesx = 1.e-3*self.samples['x']
    self.numSamples = len(self.samples['x'])

    ld = self.po2fieldld
    x0, x1, y0, y1, z0, z1 = ld.box
    self.x_field = np.linspace(ld.LatticeToWorld((x0, 0, 0))[0], ld.LatticeToWorld((x1,0,0))[0], x1-x0+1)
    self.y_field = np.linspace(ld.LatticeToWorld((0, y0, 0))[1], ld.LatticeToWorld((0, y1,0))[1], y1-y0+1)
    
    if not 'LIMITS' in self.po2field.attrs:
      self.po2field.attrs['LIMITS'] = (np.amin(self.po2field), np.amax(self.po2field))
    if not 'LIMITS' in self.po2vessels.attrs:
      self.po2vessels.attrs['LIMITS'] = (np.amin(self.po2vessels), np.amax(self.po2vessels))
    pf0, pf1 = self.po2field.attrs['LIMITS']
    pv0, pv1 = self.po2vessels.attrs['LIMITS']
    self.po2range = (min(pf0, pv0), max(pf1, pv1))
    self.params_ = None

  def GetFieldCenterLine(self):
    ldidx = self.po2fieldld.WorldToLattice((0., 0., 0.,))
    arridx = ldidx - self.po2fieldld.box.reshape(3,2)[:,0]
    p_field = np.asarray(self.po2field[:,arridx[1],arridx[2]])
    return self.x_field, p_field

  def PlotField(self, ax, sampleIndex = None, xcoord = None, color = 'k', label = '', addXPositionLabel = False, addXMark = False, addYMark = False, marker = None):
    if sampleIndex is None:
      sampleIndex = np.searchsorted(self.samples['x'], xcoord)
    if marker is not None:
      marker = niceMarker(marker, color)
    else:
      marker = dict(color = color)
    p = self.samples['pos'][sampleIndex]
    idx = self.po2fieldld.WorldToLattice(p)
    idx -= self.po2fieldld.box.reshape(3,2)[:,0]
    yv_field = np.asarray(self.po2field[idx[0],:,idx[2]])
    if addXPositionLabel:
      label += ' x = '+myutils.f2s(self.samples['x'][sampleIndex])+' $\mu m$'
    if label:
      marker['label'] = label
    ax.plot(self.y_field, yv_field, **marker)
    if addXMark:
      ax.axvline(p[1], color = 'k', ls = ':')
    if addYMark:
      po2 = self.samples['po2'][sampleIndex]
      ax.axhline(po2, ls = '-', color = color)

  def PlotFieldLongitudinal(self, ax, **kwargs): # through center?!
    x_field, p_field = self.GetFieldCenterLine()
    ax.plot(x_field, p_field, **kwargs)

  def PlotImage(self, ax):
    p = self.samples['pos'][self.numSamples/2]
    idx = self.po2fieldld.WorldToLattice(p)
    idx -= self.po2fieldld.box.reshape(3,2)[:,0]
    tmp = np.asarray(self.po2field[idx[0],:,:])
    i = np.argmax(tmp)
    _, z = np.unravel_index(i, tmp.shape)
    plotBulkTissue.imshow(ax, np.asarray(self.po2field[:,:,z]).transpose(), self.po2fieldld, interpolation = 'bilinear', crange=self.po2range)

  def AddStatsPage(self, pdfwriter):
    text = GenerateGlobalMeasurementOutputStrings(self.dataman, self.po2group, self.cachelocation)
    text += ['-----------------------']
    text += GenerateParameterOutputStrings(self.dataman, self.group, self.cachelocation)
    fig, _ = mpl_utils.MakeTextPage(text,figsize = (mpl_utils.a4size[0]*0.5, mpl_utils.a4size[1]*0.75))
    pdfwriter.savefig(fig, postfix='_data')

  def PlotSamples(self, ax, name, **kwargs):
    ax.plot(self.samplesx, self.samples[name], markevery=50, **kwargs)

  @property
  def params(self):
    if not self.params_:
      self.params_ = myutils.hdf_read_dict_hierarchy(self.po2group['parameters'])
    return self.params_
    
    
def plotAnalyzeConvergence(dataman, pdfwriter, plotties): # convergence w.r.t. lattice spacing
  def powerfunc(x, p):
    return (np.power(x, p[0])-1.)*p[1]
  
  def powerfit(x, y):
    def objectivefunction1(p):
      return (y - powerfunc(x, p))
    (resultcoeff, _) = scipy.optimize.leastsq(objectivefunction1, (0.5,1.))
    return resultcoeff
    
  fig, ax = pyplot.subplots(1,1, figsize = (mpl_utils.a4size[0]*0.25, mpl_utils.a4size[0]*0.2))
  conv_pfield = []
  for i, plt in enumerate(plotties):
    x, p = plt.GetFieldCenterLine()
    #params = myutils.hdf_read_dict_hierarchy(plt.po2group['parameters'])
    conv_pfield.append((plt.params['grid_lattice_const'],p[len(p)/2]))
  conv_pfield = np.asarray(conv_pfield).transpose()
  #lattice_constants = conv_pfield[0].copy()
  conv_pfield[1] = np.abs(conv_pfield[1] - conv_pfield[1,0])
  conv_pfield[0] /= conv_pfield[0,0]
  conv_coeff = powerfit(conv_pfield[0], conv_pfield[1])    
  ax.plot(conv_pfield[0], conv_pfield[1], color = 'k', marker = 'x', lw = 0)
  x = np.linspace(0, 6., 60)
  y = powerfunc(x, conv_coeff)
  ax.plot(x, y, label = '$(x^{%s}-1)\cdot%s$' % tuple(map(f2l, conv_coeff)), color = 'k')
  ax.set(xlabel = '$x = h / (10 \mu m)$', ylabel = '$|P_t(h) / P_t(h = 10\mu m - 1|$', ylim = (-2, 10))
  ax.legend()
  pyplot.tight_layout()
  pdfwriter.savefig(fig, postfix='_pfield_convergence')


def plotAnalyzeIterativeConvergence(dataman, pdfwriter, plotties):
  def plotIter(ax, group, **kwargs):
    iterations = plotsForPaper.GetConvergenceData(group)
    x = iterations['iteration']
    y = iterations['delta_vessM'] + iterations['delta_fieldM']
    ax.plot(x, y, **kwargs)
    
  fig, ax = pyplot.subplots(1,1, figsize = (mpl_utils.a4size[0]*0.4, mpl_utils.a4size[0]*0.4))
  mpl_utils.subplots_adjust_abs(fig, unit = 'mm', left = 15., right=-5., bottom = 10., top = -5.)
  
  ax.set(xlabel = 'iteration', ylabel = r'$|p(new) - p(old)|_\infty$', yscale = 'log')
  for i, plt in enumerate(plotties):
    plotIter(ax, plt.po2group, label = plt.label, **niceMarker(i, 'k'))
  ax.set(yscale = 'log', title = 'Convergence')
  ax.legend(loc = mpl_utils.loc.lower_left)
  axins = mpl_utils.inset_axes(ax,
                     width="30%", # width = 30% of parent_bbox
                     height="30%", # height : 1 inch
                     loc=mpl_utils.loc.upper_right)
  for i, plt in enumerate(plotties):
    marker = niceMarker(i, 'k')
    axins.plot([i], [plt.fluxes['e3']], **marker)
  axins.set(yscale = 'log', title = r'Error $e3$ [%]') #, ylim = (1.e-6, 100.0))
  #axins.xaxis.set(visible = False)    
  pdfwriter.savefig(fig, postfix='_iterations')


def plot_single_capillary(dataman, group, useInsets = False):
    # used to plot reproduction cases of nair 1989
    plt = Plotty(dataman, group, '')
    with mpl_utils.PdfWriter('sv-'+os.path.basename(group.name)+'.pdf') as pdfwriter:
      fig0, (ax0, ax2) = pyplot.subplots(2,1, figsize = (mpl_utils.a4size[0]*0.35, mpl_utils.a4size[1]*0.30))
      fig1, ax1 = pyplot.subplots(1,1, figsize = (mpl_utils.a4size[0]*0.35, mpl_utils.a4size[1]*0.15))

      plt.PlotSamples(ax0,'po2', label = 'P', **niceMarker(0, 'k'))
      plt.PlotSamples(ax0,'extpo2', label = 'P$_t$', **niceMarker(1, 'k'))
      plt.PlotSamples(ax2,'sat', label = None, **niceMarker(2,'k'))


      c = niceColors(3)
      plt.PlotField(ax1, xcoord =   200, color = c[0], addXMark = True, addYMark = True, marker = 1, label = 'x = 0.2 mm')
      plt.PlotField(ax1, xcoord =  2000, color = c[1], addXMark = True, addYMark = True, marker = 1, label = 'x = 1   mm')
      plt.PlotField(ax1, xcoord =  3800, color = c[2], addXMark = True, addYMark = True, marker = 1, label = 'x = 3.8   mm')

      ax0.legend()
      ax1.legend()
      
      ax2.set(xlabel = 'x [mm]', ylabel = 'S')
      ax0.set(xlabel = 'x [mm]', ylabel = 'PO$_2$ [mmHg]')
      ax1.set(xlabel = 'y$_1$ [$\mu$m]', ylabel = 'PO$_2$ [mmHg]', xlim = (-200., 200.)) #title = 'Tissue $PO_2$ - Transversal' ylim = (40., 100.)

      for ax in [ax0, ax1, ax2]:
        ax.yaxis.grid(True, which='major', color = '#f0f0f0', ls = '-')

      mpl_utils.tight_layout(fig0)
      mpl_utils.tight_layout(fig1)
      
      pdfwriter.savefig(fig0, postfix='_long')
      pdfwriter.savefig(fig1, postfix='_trans')

      plt.AddStatsPage(pdfwriter)
    pyplot.close('all')


def plot_single_capillary_long(dataman, group, useInsets = False):
    plt = Plotty(dataman, group, '')
    with mpl_utils.PdfWriter('sv-'+os.path.basename(group.name)+'.pdf') as pdfwriter:
      fig0, ax0 = pyplot.subplots(1,1, figsize = (mpl_utils.a4size[0]*0.35, mpl_utils.a4size[1]*0.15))
      fig1, ax1 = pyplot.subplots(1,1, figsize = (mpl_utils.a4size[0]*0.35, mpl_utils.a4size[1]*0.15))

      plt.PlotSamples(ax0,'po2', label = 'P', **niceMarker(0, 'k'))
      plt.PlotSamples(ax0,'extpo2', label = 'P$_t$', **niceMarker(1, 'k'))

      c = niceColors(4)
      plt.PlotField(ax1, xcoord =   100, color = c[0], addXMark = True, addYMark = True, marker = 1, label = 'x = 0.1 mm')
      plt.PlotField(ax1, xcoord =  1000, color = c[1], addXMark = True, addYMark = True, marker = 1, label = 'x = 1   mm')
      plt.PlotField(ax1, xcoord =  4000, color = c[2], addXMark = True, addYMark = True, marker = 1, label = 'x = 4   mm')
      plt.PlotField(ax1, xcoord = 48000, color = c[3], addXMark = True, addYMark = True, marker = 1, label = 'x = 48 mm')

      ax0.legend()
      ax1.legend()
      
      #plt.StyleSamplePlot(ax0, title = '$PO_2$ Longitudinal', yquantity = 'po2')
      #plt.StyleTransversalPlot(ax1, title = 'Tissue $PO_2$ - Transversal')

      ax0.set(xlabel = 'x [mm]', ylabel = 'PO$_2$ [mmHg]')
      ax1.set(xlabel = 'y$_1$ [$\mu$m]', ylabel = 'PO$_2$ [mmHg]', ylim = (40., 100.), xlim = (-200., 200.)) #title = 'Tissue $PO_2$ - Transversal'

      for ax in [ax0, ax1]:
        ax.yaxis.grid(True, which='major', color = '#f0f0f0', ls = '-')

      mpl_utils.tight_layout(fig0)
      mpl_utils.tight_layout(fig1)
      
      pdfwriter.savefig(fig0, postfix='_long')
      pdfwriter.savefig(fig1, postfix='_trans')

      plt.AddStatsPage(pdfwriter)
    pyplot.close('all')



def compare_single_capillary2(dataman, outfilename, grouplist, axisticks):
  plotties = []
  for group, label in grouplist:
    plotties.append(Plotty(dataman, group, label))
  colors = niceColors(len(grouplist))
  with mpl_utils.PdfWriter('comp-'+outfilename+'.pdf') as pdfwriter:
    fig, axes = pyplot.subplots(3,1, figsize = (mpl_utils.a4size[0]*0.35, mpl_utils.a4size[1]*0.45))
    axes = axes.ravel()

    for i, plt in enumerate(plotties):
      plt.PlotSamples(axes[0], 'po2', label = plt.label, **niceMarker(0, colors[i]))
      plt.PlotSamples(axes[1], 'extpo2', **niceMarker(1, colors[i]))
      plt.PlotSamples(axes[2], 'sat', label = plt.label, **niceMarker(2, colors[i]))
    axes[0].legend(loc = mpl_utils.loc.lower_left)
    axes[0].set(ylabel = r'P [mmHg]')
    axes[1].set(ylabel = r'P$_t$ [mmHg]')
    axes[2].set(ylabel = r'S')
    axes[2].set(ylim = (0.80, 1.), xlabel = '$x$ [mm]')
    for ax in axes[0:2]:
      ax.set(ylim = (50., 100.))
      #ax.xaxis.set(visible=False)
      ax.set(xticklabels = [])
    for ax in axes:
      ax.yaxis.grid(True, which='major', color = '#f0f0f0', ls = '-')
    pyplot.tight_layout()
    pdfwriter.savefig(fig, postfix='_po2')

    fig, ax = pyplot.subplots(1,1, figsize = (mpl_utils.a4size[0]*0.35, mpl_utils.a4size[1]*0.15))
    for i, plt in enumerate(plotties):
      plt.PlotField(ax, sampleIndex = plt.numSamples/2, color = colors[i], addXMark = False, addYMark = True, marker = True) #xcoord = 500
    #ax.set(title = r'$P_t$ Transv.')
    #ax.tick_params(axis='both', which='major', labelsize=8)
    ax.set(xlabel = 'y$_1$ [$\mu$ m]', ylabel = 'PO$_2$ [mmHg]', ylim = (40., 100.), xlim = (-200., 200.)) #title = 'Tissue $PO_2$ - Transversal'
    ax.yaxis.grid(True, which='major', color = '#f0f0f0', ls = '-')
    pyplot.tight_layout()
    pdfwriter.savefig(fig, postfix='_transv')
#    axins = mpl_utils.inset_axes(axes[0],
#                       width="60%", # width = 30% of parent_bbox
#                       height="30%", # height : 1 inch
#                       loc=mpl_utils.loc.upper_right)
#    axins.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(nbins = 3))
#    axins.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(nbins = 3))
#    axins.set(title = r'$P_t$ Transv.')
#    axins.tick_params(axis='both', which='major', labelsize=8)
    plotties[0].AddStatsPage(pdfwriter)
    plotAnalyzeIterativeConvergence(dataman, pdfwriter, plotties)

  pyplot.close('all')






def compare_single_capillary3(dataman, outfilename, grouplist, axisticks):
  plotties = []
  for group, label in grouplist:
    plotties.append(Plotty(dataman, group, label))
  colors = niceColors(len(grouplist))

  with mpl_utils.PdfWriter('comp-'+outfilename+'.pdf') as pdfwriter:
    fig, axes = pyplot.subplots(3,1, figsize = (mpl_utils.a4size[0]*0.4, mpl_utils.a4size[1]*0.5))
    axes = axes.ravel()
    
    for i, plt in enumerate(plotties):
      plt.PlotSamples(axes[0], 'extpo2', colors[i], 0, None)
      plt.PlotSamples(axes[0], 'po2', colors[i], 1, plt.label)
    axes[0].legend(loc = mpl_utils.loc.lower_left)
    axes[0].set(xlabel = 'x [$\mu$m]', ylabel = 'mmHg', title =  'PO$_2$ - Longitudinal')
    axes[0].set(ylim = (50., 100.))
    #axes[1].set(xlabel = 'x [$\mu m$]', ylabel = Prettyfier.get_bunit('sat'), title = Prettyfier.get_msym('sat')+' - Longitudinal')

    for i, plt in enumerate(plotties):
      plt.PlotFieldLongitudinal(axes[1], color = colors[i])
      plt.PlotSamples(axes[1], 'po2', colors[i], 1, None)
    axes[1].legend(loc = mpl_utils.loc.lower_left)
    axes[1].set(xlabel = 'x [$\mu$m]', ylabel = 'mmHg', title =  'PO$_2$ - Longitudinal')

    for i, plt in enumerate(plotties):
      plt.PlotField(axes[2], sampleIndex = plt.numSamples/2, color = colors[i], addXMark = True, addYMark = True, marker = True) #xcoord = 500
    #axins.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(nbins = 3))
    #axins.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(nbins = 3))
    axes[2].set(title = r'P$_t$ Transv.')
    axes[2].tick_params(axis='both', which='major', labelsize=8)
    axes[2].set(xlabel = 'y$_1$ [$\mu$m]', ylabel = 'mmHg', title = 'Tissue PO$_2$ - Transversal', ylim = (0., 100.))
    for ax in axes:
      ax.yaxis.grid(True, which='major', color = '#f0f0f0', ls = '-')
    pyplot.tight_layout()
    pdfwriter.savefig(fig, postfix='_curves')

    plotAnalyzeConvergence(dataman, pdfwriter, plotties)

    plotties[0].AddStatsPage(pdfwriter)

    plotAnalyzeIterativeConvergence(dataman, pdfwriter, plotties)
    
  pyplot.close('all')





if __name__ == '__main__':
  krebsutils.set_num_threads(2)
  dataman = myutils.DataManager(20, [DataDetailedPO2(),DataBasicVessel(), DataVesselSamples(), DataVesselGlobal()])

  fn = 'vessel-single-all.h5'
  #os.unlink(fn)
  f = h5files.open(fn,'a')
  
  GenerateSingleCapillaryWPo2(dataman, f, 'nair_uptake', 14, singleVesselParameterSets.nair_uptake)
  plot_single_capillary(dataman, f['nair_uptake'], useInsets = True)

  GenerateSingleCapillaryWPo2(dataman, f, 'nair_release', 14, singleVesselParameterSets.nair_release)
  plot_single_capillary(dataman, f['nair_release'], useInsets = True)

  grouplist = []
  for name in [ 'moschandreou_case%02i' % i for i in xrange(6) ]:
    params = getattr(singleVesselParameterSets, name)
    r = params.paramsTube['r']
    label = 'r = %s' % f2l(r)
    GenerateSingleCapillaryWPo2(dataman, f, name, 14, params)
    grouplist.append((f[name], label))
  compare_single_capillary2(dataman, 'moschandreou_cases', grouplist, None)
  plot_single_capillary(dataman, f['moschandreou_case02'], useInsets = True)

  GenerateSingleCapillaryWPo2(dataman, f, 'moschandreou_extra_long', 14, singleVesselParameterSets.moschandreou_extra_long)
  plot_single_capillary_long(dataman, f['moschandreou_extra_long'], useInsets = True)

#  GenerateSingleCapillaryWPo2(dataman, f, 'moschandreou_diag_base', 14, singleVesselParameterSets.moschandreou_diag_base)
#  GenerateSingleCapillaryWPo2(dataman, f, 'moschandreou_diag'    , 14, singleVesselParameterSets.moschandreou_diag)
#  compare_single_capillary2(dataman, 'moschandreou_diag', [(f['moschandreou_diag_base'], 'parallel'),
#                                                          (f['moschandreou_diag'], 'diag')], None)

#  GenerateSingleCapillaryWPo2(dataman, f, 'test_thin', 7, singleVesselParameterSets.basecase)
#  GenerateSingleCapillaryWPo2(dataman, f, 'test_straight', 7, singleVesselParameterSets.test_straight)
#  GenerateSingleCapillaryWPo2(dataman, f, 'test_diag', 7, singleVesselParameterSets.test_diag)
#  plot_single_capillary(dataman, f['test_diag'])
#  plot_single_capillary(dataman, f['test_straight'])
#  plot_single_capillary(dataman, f['test_thin'], useInsets = True)
#  compare_single_capillary2(dataman, 'test_diag', [(f['test_straight'], 'parallel'), (f['test_diag'], 'diag')], None)

#  def CompareGridSizeSeries(basename):
#    grouplist = []
#    for i in xrange(0,6):
#      name = '%s-num%02i' % (basename, i)
#      params = getattr(singleVesselParameterSets, name)
#      GenerateSingleCapillaryWPo2(dataman, f, name, 9, params)
#      label = 'h = $%s$' % f2l(params.paramspo2['grid_lattice_const'])
#      grouplist.append((f[name], label))
#      plot_single_capillary(dataman, f[name], useInsets = True)
#    compare_single_capillary3(dataman, basename, grouplist, axisticks = None)
#
#  CompareGridSizeSeries('hseries1')
#  CompareGridSizeSeries('hseries1diag')

#  GenerateSingleCapillaryWPo2(dataman, f, 'hseries1-num00', 9, getattr(singleVesselParameterSets, 'hseries1-num00'))
#  plot_single_capillary(dataman, f['hseries1-num00'], useInsets = True)

#  GenerateSingleCapillaryWPo2(dataman, f, 'hseries1-num05', 9, getattr(singleVesselParameterSets, 'hseries1-num05'))
#  plot_single_capillary(dataman, f['hseries1-num05'], useInsets = True)

#  GenerateSingleCapillaryWPo2(dataman, f, 'thierry_case', 17, singleVesselParameterSets.thierry_case)
#  plot_single_capillary(dataman, f['thierry_case'], useInsets = True)