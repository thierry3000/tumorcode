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
# -*- coding: utf-8 -*-
if __name__ == '__main__':
  import os.path, sys
  sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','..'))

import os,sys
from os.path import join, basename, dirname, splitext
import h5py
import h5files
import numpy as np
import extensions # for asarray with h5py support
import krebsutils
import posixpath
import math
import glob
import itertools
from pprint import pprint
import collections

import scipy.ndimage
import scipy.signal
import scipy.stats

from krebs import analyzeGeneral

import myutils
import mpl_utils

import matplotlib
import matplotlib.pyplot as pyplot

from analyzeBloodPressureMvdCorrelation import DataPressureMvdCorrelation
from detailedo2Analysis.plotsForPaper import GetVesselTypeLabel, RewriteVesselLabel, vesselTypeColors, vesselTypeMarkers


if __name__ == '__main__':
  import optparse  #Note: Deprecated since version 2.7. Use argparse instead
  parser = optparse.OptionParser()
  parser.add_option("-n",dest="num_samples", help="number of samples", default=10000, type = "int")
  options, args = parser.parse_args()
  filenames = args[:-2]
  pattern_before = args[-2]
  pattern_after  = args[-1]  
  print '-----------  looking for files  ----------'
  thegroups = []
  thefiles   = {}
  for filename in filenames:
    f = h5files.open(filename)
    thefiles[filename] = f
    print 'opened -- ', filename,'/',
    paths_before = myutils.walkh5(f['.'], pattern_before)
    paths_after  = myutils.walkh5(f['.'], pattern_after)
    print paths_before, paths_after
    for path_before, path_after in zip(paths_before, paths_after):
      vesselgroup_before = f[path_before]
      vesselgroup_after  = f[path_after]
      tumorgroup  = analyzeGeneral.try_find_tumor_group_from_vesselgroup(vesselgroup_after)
      assert tumorgroup
      thegroups.append([vesselgroup_before, vesselgroup_after, tumorgroup])  
  
  prefix, suffix = myutils.splitcommonpresuffix(map(lambda s: basename(s), filenames))
  outputbasename, _ = splitext(prefix+suffix)  
  
  fn_measure = join(dirname(outputbasename), 'common-mvd-grad-map-cache.h5')
  f_measure = h5files.open(fn_measure, 'a')
  
  def cachelocation(dataname, vesselgroup, tumorgroup, version):
    path = myutils.splitPath(posixpath.join(splitext(basename(vesselgroup.file.filename))[0], vesselgroup.name.strip(posixpath.sep)))+(dataname,)
    return (f_measure, path)

  def cachelocationEnsemble(dataname, groups):
    groupchecksum = myutils.checksum(*map(lambda g: str(g.file.filename+'/'+g.name), sum(groups, [])))
    path = ('%s_%s' % (dataname, groupchecksum),)
    return (f_measure, path)
  
  dataman = myutils.DataManager(20, [ analyzeGeneral.DataTumorTissueSingle(), 
                                      analyzeGeneral.DataDistanceFromCenter(), 
                                      analyzeGeneral.DataBasicVessel(), 
                                      analyzeGeneral.DataVesselSamples(),
                                      DataPressureMvdCorrelation(200., 5, 30., cachelocation, cachelocationEnsemble)])

  print '-----------computing sammples------------'
  localSamples = dataman.obtain_data('intervascular_map_correlations', thegroups)
  globalSamples = dataman.obtain_data('intervascular_global_correlations', thegroups)
  print '--------------plotting ---------------'

  with mpl_utils.PageWriter(outputbasename+'_mvd-grad.pdf', fileformats=['svg']) as pdfwriter:
    fig, axes = pyplot.subplots(1,2, figsize = mpl_utils.a4size*np.asarray((0.8,0.2)))
    ax = axes[0]

    dataPage = []
  
    samples = localSamples
    mask = np.asarray(samples['mask'])        
    maski = ~mask
    maski[np.random.rand(len(mask)) < 0.95] = False  
    r = np.random.rand(len(mask))
    dilution = float(options.num_samples) / len(mask)
    mask[r > dilution] = False
    maski[r > dilution] = False  
    normalizationConst = 1./np.amax(samples['grad'])
    
    filemap = np.asarray(samples['filemapping']).reshape(-1,3)
    numberToLabel = {}
    indexToNumber = []
    for fn, _, _ in filemap:
      f = thefiles[fn]
      vesseltype = GetVesselTypeLabel(f['/'])
      number = ord(vesseltype)-ord('A')
      label = RewriteVesselLabel(vesseltype)
      numberToLabel[number] = label
      indexToNumber.append(number)
    indexToNumber = np.asarray(indexToNumber)
    sampleNumbers = indexToNumber[np.asarray(samples['fileindex'])]  # which root config (in integers as ord(vesseltype)) does a sample belong to
        

    normalCorrelation, normalP = scipy.stats.pearsonr(samples['mvd'][maski], samples['grad'][maski])
    tumorCorrelation, tumorP = scipy.stats.pearsonr(samples['mvd'][mask], samples['grad'][mask])

    stylestuff = dict(markerfacecolor='none', linewidth=0, ms = 3.)
    for number, label in numberToLabel.items():
#      combimask = np.logical_and(sampleNumbers == number, maski)
#      ax.plot(samples['mvd'][combimask], samples['grad'][combimask], color = 'k', markeredgecolor = 'k', marker = vesselTypeMarkers[number], zorder = 2, **stylestuff)
      combimask = np.logical_and(sampleNumbers == number, mask)
      ax.plot(samples['mvd'][combimask], normalizationConst*samples['grad'][combimask], color = vesselTypeColors[number], markeredgecolor = vesselTypeColors[number], label = label, zorder = 4, marker = vesselTypeMarkers[number], **stylestuff)
    
    dataPage += ['Local:',
                 '<(mvd - <mvd>)(grad - <grad>)> = %f, p = %f (normal)' % (normalCorrelation, normalP),
                 '<(mvd - <mvd>)(grad - <grad>)> = %f, p = %f (tumor)' % (tumorCorrelation, tumorP),
                ]
    ax.legend()

    del mask, maski, dilution    
    
    ax = axes[1]
    
    samples = globalSamples
    mask = samples['mask']
    maski = ~mask
    normalizationConst = 1./np.amax(samples['grad'])

    sampleNumbers = indexToNumber[np.asarray(samples['fileindex'])]

    normalCorrelation, normalP = scipy.stats.pearsonr(samples['mvd'][maski], samples['grad'][maski])
    tumorCorrelation, tumorP = scipy.stats.pearsonr(samples['mvd'][mask], samples['grad'][mask])
    
    stylestuff['ms'] = 5.
    for number, label in numberToLabel.items():
#      combimask = np.logical_and(sampleNumbers == number, maski)
#      ax.plot(samples['mvd'][combimask], samples['grad'][combimask], color = 'k', markeredgecolor = 'k', marker = vesselTypeMarkers[number], zorder = 2, **stylestuff)
      combimask = np.logical_and(sampleNumbers == number, mask)
      ax.plot(samples['mvd'][combimask], normalizationConst*samples['grad'][combimask], color = vesselTypeColors[number], markeredgecolor = vesselTypeColors[number], label = label, zorder = 4, marker = vesselTypeMarkers[number], **stylestuff)
    ax.legend()

    dataPage += ['Global:',
                 '<(mvd - <mvd>)(grad - <grad>)> = %f, p = %f (normal)' % (normalCorrelation, normalP),
                 '<(mvd - <mvd>)(grad - <grad>)> = %f, p = %f (tumor)' % (tumorCorrelation, tumorP),
                ]
    ax.text(0.01, 0.95, '\n'.join(dataPage), verticalalignment = 'top')
    pdfwriter.savefig(fig)
