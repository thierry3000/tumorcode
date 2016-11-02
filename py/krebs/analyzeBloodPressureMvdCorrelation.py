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
  sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'..'))
  
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

import myutils
import mpl_utils
import analyzeGeneral

import matplotlib
import matplotlib.pyplot as pyplot



class DataPressureMvdCorrelation(object):  
    keywords = [
      'local_mvd_map', 'intervascular_pressure_map', 'intervascular_map_common_ld', 'intervascular_map_tumor_mask', 'intervascular_map_correlations', 'intervascular_global_correlations'
    ]

    def makeLD(self, vesselgroup):
      ldvessels = krebsutils.read_lattice_data_from_hdf(vesselgroup['lattice'])
      fieldld = krebsutils.SetupFieldLattice(ldvessels.worldBox, 3, self.bin_size, 0.)
      fieldldFine = krebsutils.SetupFieldLattice(fieldld.worldBox, 3, self.bin_size / self.fine_bin_subdivision, 0.)
      return fieldld, fieldldFine
    
    def makeCacheLocation(self, dataname, args, name, ver = 1):
      vesselgroup, tumorgroup, fieldLd, fieldLdFine = args
      myfile, mycachelocation = self.cachelocationFactory(dataname, vesselgroup, tumorgroup, myutils.checksum(fieldLd.box, fieldLd.scale, self.fine_bin_subdivision))
      version         = (None,)*len(mycachelocation) + (ver,)
      mycachelocation = mycachelocation + (name,)
      return myfile, mycachelocation, version   
    
    def __init__(self, bin_size, fine_bin_subdivision, sample_length, cachelocationFactory, cachelocationEnsembleFactory):
      self.bin_size = bin_size
      self.fine_bin_subdivision = fine_bin_subdivision
      self.sample_length = sample_length
      self.cachelocationFactory = cachelocationFactory
      self.cachelocationEnsembleFactory = cachelocationEnsembleFactory

    def obtain_data(self, dataman, dataname, *args):
      if dataname == 'intervascular_map_common_ld':
        vesselgroup, tumorgroup = args
        fieldLd, fieldLdFine = self.makeLD(vesselgroup)
        return fieldLd, fieldLdFine

      if dataname == 'intervascular_map_tumor_mask':
        vesselgroup, tumorgroup, fieldLd, fieldLdFine = args
        print 'intervascular_map_tumor_mask', str(vesselgroup)
        def read(gmeasure, groupname):
          return gmeasure[groupname]          
        def write(gmeasure, groupname):
          distmap, _ = analyzeGeneral.obtain_distmap_(dataman, tumorgroup, 'levelset', fieldLd)
          mask = distmap < -fieldLd.scale-100.
          gmeasure.create_dataset(groupname, data = mask, compression = 9)          
        myfile, mycachelocation, version = self.makeCacheLocation(dataname, args, 'mask', 3)
        return myutils.hdf_data_caching(read, write, myfile, mycachelocation, version)
              
      if dataname == 'local_mvd_map':
        vesselgroup, tumorgroup, fieldLd, fieldLdFine = args
        print 'local_mvd_map', str(vesselgroup)
        def read(gmeasure, groupname):
          return gmeasure[groupname]
        def write(gmeasure, groupname):
          weight   = dataman.obtain_data('basic_vessel_samples', 'weight', vesselgroup, self.sample_length)
          flags    = dataman.obtain_data('basic_vessel_samples', 'flags', vesselgroup, self.sample_length)
          position = dataman.obtain_data('basic_vessel_samples', 'position', vesselgroup, self.sample_length)            
          mask = myutils.bbitwise_and(flags, krebsutils.CIRCULATED)
          flags    = flags[mask]
          position = position[mask,...]
          weight   = weight[mask]
          ####### put into bins
          eps = 1.0-1.e-15
          x0, x1, y0, y1, z0, z1 = fieldLd.worldBox
          ranges = [np.arange(x0, x1, fieldLd.scale*eps),
                    np.arange(y0, y1, fieldLd.scale*eps),
                    np.arange(z0, z1, fieldLd.scale*eps),
                    ]
          mvd, _ = np.histogramdd(position, bins = ranges, weights = weight)
          mvd *= 1.e6/(fieldLd.scale**3)
          ####### save
          gmeasure.create_dataset(groupname, data = mvd, compression = 9, dtype = np.float32)
        myfile, mycachelocation, version = self.makeCacheLocation(dataname, args, 'mvd')
        return myutils.hdf_data_caching(read, write, myfile, mycachelocation, version)
      
      if dataname == 'intervascular_pressure_map':
        vesselgroup, tumorgroup, fieldLd, fieldLdFine = args
        print 'intervascular_pressure_map', str(vesselgroup)
        def read(gmeasure, groupname):
          return gmeasure[groupname]      
        def write(gmeasure, groupname):
          graph = dataman.obtain_data('vessel_graph', vesselgroup, ['position', 'radius', 'flags', 'pressure'])
          graph = graph.get_filtered(myutils.bbitwise_and(graph['flags'], krebsutils.CIRCULATED))
          edgevalues = graph['pressure']
          edgevalues = edgevalues[graph.edgelist]
          # main calculation
          thefield = krebsutils.CalcIntervascularInterpolationField(graph.edgelist, graph['radius'], graph['position'], edgevalues, fieldLdFine, 1.)
          del edgevalues
          # gradient of the interpolated blood pressure
          thegrad     = scipy.ndimage.gaussian_filter(thefield, 1.0, 0, mode='nearest')
          gradfield_  = krebsutils.field_gradient(thegrad)
          thegrad = np.sum([np.square(g) for g in gradfield_], axis = 0)
          thegrad = np.sqrt(thegrad)
          del gradfield_, thefield          
          # now we scale down the highres version, TODO: make it work with other than 3 dimensions
          # first local average
          m = self.fine_bin_subdivision
          kernel = np.ones((m,m,m), dtype = np.float32)
          thegrad = scipy.signal.fftconvolve(thegrad, kernel, mode = 'valid')
          # then pick every m'th which contains the average of m finer boxes combined
          thegrad = np.ascontiguousarray(thegrad[::m,::m,::m])
          assert all(thegrad.shape == fieldLd.shape)
          gmeasure.create_dataset(groupname, data = thegrad, compression = 9, dtype = np.float32)          
        myfile, mycachelocation, version = self.makeCacheLocation(dataname, args, 'grad')
        return myutils.hdf_data_caching(read, write, myfile, mycachelocation, version)
      
      if dataname == 'intervascular_map_correlations':
        listofgroups, = args
        def read(gmeasure, groupname):
          gmeasure = gmeasure[groupname]
          return gmeasure
        def write(gmeasure, groupname):
          allSamples = []
          filemapping = []
          for i,(vesselgroup_before, vesselgroup_after, tumorgroup_after) in enumerate(listofgroups):
            fieldLd, fieldLdFine = dataman.obtain_data('intervascular_map_common_ld', vesselgroup_before, None)
            mask = dataman.obtain_data('intervascular_map_tumor_mask', vesselgroup_after, tumorgroup_after, fieldLd, fieldLdFine)
            mvd = dataman.obtain_data('local_mvd_map', vesselgroup_after, None, fieldLd, fieldLdFine)
            grad = dataman.obtain_data('intervascular_pressure_map', vesselgroup_before, None, fieldLd, fieldLdFine)
            mvd = np.asarray(mvd).ravel()
            grad = np.asarray(grad).ravel()
            mask = np.asarray(mask).ravel()
            data = (i * np.ones(mask.shape, dtype = np.int), mask, mvd, grad)
            allSamples += zip(*data)
            stuff = map(unicode.encode, (vesselgroup_before.file.filename, vesselgroup_before.name, vesselgroup_after.name))
            filemapping += stuff            
          allSamples = zip(*allSamples)
          gmeasure = gmeasure.create_group(groupname)
          gmeasure.create_dataset('fileindex', data = allSamples[0]) # index into filemapping array
          gmeasure.create_dataset('mask'     , data = allSamples[1]) # equal true for samples within tumor (?)
          gmeasure.create_dataset('mvd'      , data = allSamples[2])
          gmeasure.create_dataset('grad'     , data = allSamples[3])
          gmeasure.create_dataset('filemapping', data = np.asarray(filemapping))
        myfile, mycachelocation, = self.cachelocationEnsembleFactory(dataname, listofgroups)
        version = (None,)*(len(mycachelocation)-1) + (4,)
        return myutils.hdf_data_caching(read, write, myfile, mycachelocation, version)

      if dataname == 'intervascular_global_correlations':
        listofgroups, = args
        group = dataman.obtain_data('intervascular_map_correlations', listofgroups)
        sorteddata = collections.defaultdict(list)
        stuff = map(lambda s: np.asarray(group[s]), 'fileindex mask mvd grad'.split())
        stuff = zip(*stuff) # transposed
        for fileindex, mask, mvd, grad in stuff:
          sorteddata[fileindex, mask].append((mvd, grad))
        result = []
        for (fileindex,mask), mvd_grad in sorteddata.iteritems():
          mvd_grad = np.average(mvd_grad, axis=0) # result: 2 elements: (mvd, grad)
          result.append((fileindex, mask, mvd_grad[0], mvd_grad[1]))
        result = sorted(result, key = lambda t: t[0])
        result = zip(*result)
        result = dict(zip('fileindex mask mvd grad'.split(), map(np.asarray, result)))
        return result # returns dict of fileindex, mask, mvd, grad


if __name__ == '__main__':
  import optparse  #Note: Deprecated since version 2.7. Use argparse instead
  parser = optparse.OptionParser()
  parser.add_option("-n",dest="num_samples", help="number of samples", default=10000, type = "int")
  parser.add_option("--picz", dest="picz", help="make picz", default = False, action = 'store_true')
  options, args = parser.parse_args()
  filenames = args[:-2]
  pattern_before = args[-2]
  pattern_after  = args[-1]  
  print '-----------  looking for files  ----------'
  thegroups = []
  for filename in filenames:
    f = h5files.open(filename)
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

  if not options.picz:
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
    
        normalCorrelation, normalP = scipy.stats.pearsonr(samples['mvd'][maski], samples['grad'][maski])
        tumorCorrelation, tumorP = scipy.stats.pearsonr(samples['mvd'][mask], samples['grad'][mask])
        
        ax.plot(samples['mvd'][maski], samples['grad'][maski], color = 'k', label = 'normal', lw = 0, marker = 'x', ms = 2.)
        ax.plot(samples['mvd'][mask], samples['grad'][mask], color = 'r', label = 'tumor', lw = 0, marker = '+', ms = 4.)
        
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
    
        normalCorrelation, normalP = scipy.stats.pearsonr(samples['mvd'][maski], samples['grad'][maski])
        tumorCorrelation, tumorP = scipy.stats.pearsonr(samples['mvd'][mask], samples['grad'][mask])
        
        ax.plot(samples['mvd'][maski], samples['grad'][maski], color = 'k', label = 'normal', lw = 0, marker = 'x', ms = 2.)
        ax.plot(samples['mvd'][mask], samples['grad'][mask], color = 'r', label = 'tumor', lw = 0, marker = '+', ms = 4.)
        ax.legend()
    
        dataPage += ['Global:',
                     '<(mvd - <mvd>)(grad - <grad>)> = %f, p = %f (normal)' % (normalCorrelation, normalP),
                     '<(mvd - <mvd>)(grad - <grad>)> = %f, p = %f (tumor)' % (tumorCorrelation, tumorP),
                    ]
        ax.text(0.01, 0.95, dataPage, verticalalignment = 'top')
        pdfwriter.savefig(fig)
  else:
    maxGrad = 0.
    maxMVD  = 0.
    picz = []
    for i,(vesselgroup_before, vesselgroup_after, tumorgroup_after) in enumerate(thegroups):
      fieldLd, fieldLdFine = dataman.obtain_data('intervascular_map_common_ld', vesselgroup_before, None)
      mask = dataman.obtain_data('intervascular_map_tumor_mask', vesselgroup_after, tumorgroup_after, fieldLd, fieldLdFine)
      mvd = dataman.obtain_data('local_mvd_map', vesselgroup_after, None, fieldLd, fieldLdFine)
      grad = dataman.obtain_data('intervascular_pressure_map', vesselgroup_before, None, fieldLd, fieldLdFine)
      z = mask.shape[2]//2
      mask = np.transpose(np.asarray(mask[:,:,z], dtype = np.float32))
      mvd  = np.transpose(mvd[:,:,z])
      grad = np.transpose(grad[:,:,z])
      maxGrad = max(maxGrad, np.amax(grad))
      maxMVD = max(maxMVD, np.amax(mvd))
      extents = fieldLd.worldBox[:4]*1.0e-3
      picz.append((vesselgroup_before.file.filename, mask, mvd, grad, extents))
    h5files.closeall()

    picz = list(reversed(picz))
    with mpl_utils.PageWriter(outputbasename+'_mvd-grad-picz.pdf', fileformats = ['svg']) as pdfwriter:
      while picz:
        fig, axes = pyplot.subplots(5, 4, figsize = mpl_utils.a4size)
        axes = np.reshape(axes, (-1,2))
        for ax0, ax1 in axes:
          if not picz:
            break
          filename, mask, mvd, grad, extents = picz.pop()
          #now plot mvd, and gradient
          plt0 = ax0.imshow(
            mvd,
            extent = extents,
            origin = 'lower',
            interpolation = 'nearest',
            vmin = 0.,
            vmax = maxMVD)
          plt1 = ax1.imshow(
            grad,
            extent = extents,
            origin = 'lower',
            interpolation = 'nearest',
            vmin = 0.,
            vmax = maxGrad)
          ax0.contour(
            mask,
            levels = [0.],
            extent = extents,
            linewidths = 0.25,
            colors = 'w')
          ax1.contour(
            mask,
            levels = [0.],
            extent = extents,
            linewidths = 0.25,
            colors = 'w')
          mpl_utils.remove_frame(ax0, ax1)
          ax0.set(title = 'MVD    ('+filename+')')
          ax1.set(title = 'Grad')
          divider = mpl_utils.make_axes_locatable(ax0)
          cax = divider.append_axes("right", size = "5%", pad = 0.05)
          fig.colorbar(plt0, cax = cax)
          divider = mpl_utils.make_axes_locatable(ax1)
          cax = divider.append_axes("right", size = "5%", pad = 0.05)
          fig.colorbar(plt1, cax = cax)   
        #pyplot.tight_layout()
        pdfwriter.savefig(fig)