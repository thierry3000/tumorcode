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
import numpy as np
import extensions # for hdf5 support in np.asarray
import krebsutils
import myutils
import posixpath
import math
import collections
import itertools

import detailedo2
import detailedo2Analysis
import analyzeGeneral
import analyzeBloodFlow
import qsub

import matplotlib.pyplot
import mpl_utils

from detailedo2Analysis.plotsForPaper import ComputeRegionalHistogramsOfPo2Group

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

class MeasurementInfo(object):
  def __init__(self, **kwargs):
    self.sample_length = 30.
    self.cachelocation_callback = None
    self.distancemap_spec = 'radial'
    for k,v in kwargs.iteritems():
      setattr(self, k, v)

class EnsembleItem(object):
  def __init__(self, **kwargs):
    self.path = ''
    self.time = None
    self.po2group_w_a = None
    self.po2group_no_a = None
    self.gvessels_w_a = None
    self.gvessels_no_a = None
    self.gtumor = None
    self.vessel_system_length = 0.
    self.initialVesselType = ''
    for k, v in kwargs.iteritems():
      setattr(self, k, v)
      
def GetVesselTypeLabel(po2group):
  import re
  if 'VESSELFILE_MESSAGE' in po2group.file.attrs.keys():
    msg = po2group.file.attrs['VESSELFILE_MESSAGE']
  else:
    #generic try
    msg = po2group.file.filename
  m = re.search('type(\w)', msg)
  if not m:
    return 'unkown'
  return m.group(1)
  
class EnsembleFiles(object):
  def __init__(self, dataman, filenames, pattern):
    files     = [h5files.open(fn, 'r+') for fn in filenames]
    items     = []
    has_tumor = True
    for f in files:
      paths = myutils.walkh5(f['.'], pattern)
      for path in paths:
        po2group_w_a = f[path+'/vessels_after_adaption']
        gvessels_w_a, gtumor = detailedo2.OpenVesselAndTumorGroups(po2group_w_a)
        po2group_no_a = f[path+'/recomputed']
        gvessels_no_a, gtumor = detailedo2.OpenVesselAndTumorGroups(po2group_no_a)
        e = EnsembleItem(path = path, po2group_w_a = po2group_w_a, gvessels_w_a = gvessels_w_a, po2group_no_a = po2group_no_a, gvessels_no_a = gvessels_no_a,gtumor = gtumor)
#        if 'SOURCE' in po2group:
#          source = h5files.openLink(po2group, 'SOURCE')
#          if 'time' in source.attrs.keys():
#            t = source.attrs['time']
#            e.time = t
        has_tumor = has_tumor and gtumor is not None
        e.vessel_system_length = dataman.obtain_data('vessel_system_length', gvessels_no_a)
        e.initialVesselType = GetVesselTypeLabel(po2group_no_a)
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
    if has_tumor:
      self.tumor_snapshots = tumor_snapshots # list of tuple(items, path, time)
      self.has_tumor = has_tumor
    self.o2ConfigName = set(item.po2group_no_a.attrs.get('O2_CONFIG_NAME',None) for item in items)
    if len(self.o2ConfigName) != 1:
      raise RuntimeError("Detected different O2_CONFIG_NAMES %s. You don't want to mix configurations, do you?" % self.o2ConfigName)
    self.o2ConfigName = self.o2ConfigName.pop()

def compare_tissue_saturation(dataman, ensemble, pdfwriter):
  for item in ensemble.items:
    _, po2ld, po2field_no_a, parameters = dataman.obtain_data('detailedPO2', item.po2group_no_a)
    _, po2ld, po2field_w_a, parameters = dataman.obtain_data('detailedPO2', item.po2group_w_a)  
  
  fig = matplotlib.pyplot.figure()
  fig.suptitle('tissue O2')
  ax1 = fig.add_subplot(111)
  no_bins=100
  ax1.hist(np.asarray(po2field_no_a).ravel()-np.asarray(po2field_w_a).ravel(),no_bins)
  ax1.grid(True)
  
  pdfwriter.savefig(fig, postfix='_tissue_saturation')
  
  
def compare_vessel_saturation(dataman, ensemble, pdfwriter):
  for item in ensemble.items:
    _, po2ld, po2field_no_a, parameters = dataman.obtain_data('detailedPO2', item.po2group_no_a)
    _, po2ld, po2field_w_a, parameters = dataman.obtain_data('detailedPO2', item.po2group_w_a)  
  
  fig = matplotlib.pyplot.figure()
  fig.suptitle('tissue O2')
  ax1 = fig.add_subplot(111)
  no_bins=100
  ax1.hist(np.asarray(po2field_no_a).ravel()-np.asarray(po2field_w_a).ravel(),no_bins)
  ax1.grid(True)
  
  pdfwriter.savefig(fig, postfix='_tissue_saturation')


 
def doit(filenames):
  dataman = myutils.DataManager(50, map(lambda x: x(), detailedo2Analysis.O2DataHandlers) + [ analyzeGeneral.DataTumorTissueSingle(), analyzeGeneral.DataDistanceFromCenter(), analyzeGeneral.DataBasicVessel(), analyzeGeneral.DataVesselSamples(), analyzeGeneral.DataVesselRadial(), analyzeGeneral.DataVesselGlobal(), analyzeBloodFlow.DataTumorBloodFlow()])
  ensemble = EnsembleFiles(dataman, filenames,'po2/adaption/' )
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
    if 0:    
      compare_tissue_saturation(dataman, ensemble, pdfwriter)
    
    if 0:
      #try:
      #  histogramGroupFinal   = f_measure['combinedHistogramsFinal']
      #  histogramGroupInitial = f_measure['combinedHistogramsInitial']
      #except KeyError:
      
      #histogramGroupFinal   = f_measure.recreate_group('combinedHistogramsFinal')
      histogramGroupInitial = f_measure.recreate_group('combinedHistogramsInitial')        
      #ComputeHistogramsOfPo2Items(dataman, ensemble.items, measurementinfo, histogramGroupFinal)
      ComputeHistogramsOfPo2Items(dataman, ensemble.items, measurementinfo, histogramGroupInitial)
      #PlotHistograms(pdfwriter, histogramGroupFinal, 'tum', 'Tumor')
      PlotHistograms(pdfwriter, histogramGroupInitial, 'all', 'Initial')

if __name__ == '__main__':
  krebsutils.set_num_threads(2)
  import argparse
  parser = argparse.ArgumentParser(description='Compare O2 from adation and no adaption')  
  parser.add_argument('detailedO2FileNames', nargs='+', type=argparse.FileType('r'), default=sys.stdin, help='detailedO2 files to calculate')   
  #parser.add_argument('grp_pattern',help='Where to find the oxygendata. Usually this is somthing with po2/out*')        
  goodArguments, otherArguments = parser.parse_known_args()
  qsub.parse_args(otherArguments)
  
  #create filename due to former standards
  filenames=[]
  for fn in goodArguments.detailedO2FileNames:
    filenames.append(fn.name)
  try:
    for fn in filenames:
      if not os.path.isfile(fn):
        raise AssertionError('The file %s is not present!'%fn)
  except Exception, e:
    print e.message
    sys.exit(-1)
  doit(filenames)