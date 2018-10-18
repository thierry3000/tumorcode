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
  
import os, sys
from os.path import join, basename, dirname
import h5py
import numpy as np
import extensions # for asarray with h5py support
import krebsutils
import posixpath
import math
import glob
import time
import md5
import re
import cPickle
from collections import defaultdict
from pprint import pprint
import qsub

#import vtkcommon
from mystruct import Struct
import myutils

from plotBulkTissue import commonOutputName, colorbar, contour, imslice, imshow
from analyzeGeneral import calc_distmap, CalcPhiVessels, DataDistanceFromCenter, DataBasicVessel

import matplotlib
import matplotlib.pyplot as pyplot
import mpl_utils

import plotVessels


##############################################################################
### labels and color setup
##############################################################################

mastersize = mpl_utils.a4size
a4size = mpl_utils.a4size
gridcolor = (0.7,0.7,0.7)

def mk_CM_():
  grey = matplotlib.cm.gray_r
  grey_r = matplotlib.cm.gray
  spectral = matplotlib.cm.Spectral
  return Struct(locals())
ColorMaps = mk_CM_()
CM = ColorMaps


#see http://assorted-experience.blogspot.com/2007/07/custom-colormaps.html
#xa = np.linspace(0., 1., num = 256)
#ya = 1.-xa #np.power(1.-xa, 1./2.4)
#l = tuple((x, y, y) for x, y in zip(xa, ya))
#cm_grey = matplotlib.colors.LinearSegmentedColormap('greys_gamma', {
#    'red' :   l,
#    'green' : l,
#    'blue' :  l
#  })

#colors = [(1,1,1), (.3,.3,.9), (.3,.9,0.3), (.9,.3,.3), (.9,.9,.5)]
colorconverter = matplotlib.colors.ColorConverter()
colors = '#ffffff #446bdb #7a7a7a #cd3737 #fffebe'.split()
#colors = '#ffffff #446bdb #bb37c7 #ea7243 #fffebe'.split()
#colors = '#ffffff #8b96ff #868686 #740000 #000000'.split()
colors = [ colorconverter.to_rgb(c) for c in colors ]
#vals = [0, 1.e-3, 1.e-2, 1.e-1, 1.e-0, 1]
vals = [0., 0.25, 0.5, 0.75, 1.]
#cm_spectral = matplotlib.colors.LinearSegmentedColormap('special log', {
#  'red'   : [(x, r, r) for (x,(r,g,b)) in zip(vals,colors)],
#  'green' : [(x, g, g) for (x,(r,g,b)) in zip(vals,colors)],
#  'blue'  : [(x, b, b) for (x,(r,g,b)) in zip(vals,colors)],
#}, N = 256, gamma = 1.)
fig_numbering = [ x for x in 'ABCDEFGHIJ' ]


def mk_LF_():
  levelsetfunc = ur'\theta'
  dqualimax = ur'ICMAX'
  dqualiauc = ur'ICAUC'
  avgOver = lambda x,y: ur'\langle %s \rangle_{%s}' % (x,y)
  stdOver = lambda x,y: ur'\mathrm{std}_{%s}(%s)' % (y,x)
  avgOverTumor = lambda x: ur'\langle %s \rangle_T' % x
  minOverTumor = lambda x: ur'\min_T(%s)' % x
  maxOverTumor = lambda x: ur'\max_T(%s)' % x
  avgOverViableTumor = lambda x: ur'\langle %s \rangle_V' % x
  stdOverViableTumor = lambda x: ur'\mathrm{std}_V(%s)' % x
  probForIn = lambda x: ur'\tilde{p}_{%s}' % x
  probForInViableTumorAuc = ur'\tilde{p}(%s)' % dqualiauc
  probForInViableTumorMax = ur'\tilde{p}(%s)' % dqualimax
  math = lambda x: ur'$%s$' % x
  source = 'Q'
  transmural_flow_coeff = ur'L^v_l 2 \pi r'
  transmural_flow_per_length = ur'L^v_l 2 \pi r_v (p_v - p_i)'
  relative_extravasation_flow = ur'\lambda'
  vessel_flow_rate = ur'q_v'

  return Struct(locals())

LabelFactory = mk_LF_()

LF = LabelFactory

lf2s = lambda x: myutils.f2s(x, latex=True)

def text1(ax, txt):
  ax.text(0.1, 0.9, txt, ha = "left", transform = ax.transAxes)

def text2(ax, txt):
  ax.text(0.01, 0.9, txt, ha = "left", transform = ax.transAxes)


##############################################################################
### measurement and cache management
##############################################################################


bins_rad = np.arange(0, 20000, 30.)
bins_dist = np.arange(-10000., 10000., 30.)
bins_r = dict(vs_r = bins_rad, vs_dr = bins_dist)
sample_length = 50.
far = 1.e10
stored_sample_fraction_normal = 0.5 #0.1
stored_sample_fraction_tumor  = 0.5 #1.


def downcast(a):
  if a.dtype == np.float64: a = np.asanyarray(a, dtype = np.float32)
  return a


def getVesselsTotalLength(vessels):
    pos = vessels['position']
    e = vessels.edgelist
    dp = pos[e[:,0],:] - pos[e[:,1],:]
    l = krebsutils.vector_lengths(dp)
    return np.sum(l)


class DataTissue(object):
  keywords = [
    'ld', 'phi_vessels', 'tumor_composition', 'dist_vessels', 'tissue_composition',
    'iffvelocity', 'iffvelocity_magnitude', 'iffvelocity_outward', 'iff_pressure', 'iff_sources',
    'mask_tumor', 'mask_viabletumor',
  ]

  def obtain_data(self, dataman, dataname, *args):
    f, args = args[0], args[1:]
    obtain_data = lambda *args: dataman.obtain_data(args[0], f, args[1:])

    if dataname == 'ld':
      ld = krebsutils.read_lattice_data_from_hdf(f['field_ld'])
      return ld

    ld = obtain_data('ld')

    #####
    if dataname == 'phi_vessels':
      def read(gmeasure, name):
        return np.asarray(gmeasure[name])

      def write(gmeasure, name):
        phi_vessels = CalcPhiVessels(dataman, f['iff/vessels'], ld, scaling = 1.)
        gmeasure.create_dataset(name, data = phi_vessels, compression = 9)

      return myutils.hdf_data_caching(read, write, f, ('measurements', 'phi_vessels'), (None,2))

    #####
    if dataname == 'dist_vessels':
      phi_vessels = obtain_data('phi_vessels')
      return krebsutils.distancemap(np.asarray(phi_vessels > 0.01, dtype=np.float32))*ld.scale

    #####
    if dataname == 'tumor_composition':
      if 'iff/phi_cells' in f:
        phi_cells = np.asarray(f['iff/phi_cells'])
      else:
        phi_cells = 0.5 * np.ones(ld.shape)
      tumor_contour_level = 0.5 * np.average(phi_cells)
      theta_tumor = np.asarray(f['iff/theta_tumor'])
      theta_necro = np.asarray(f['iff/theta_necro'])
      phi_viabletumor = theta_tumor * (phi_cells - theta_necro)

      def read(gmeasure, groupname):
        return np.asarray(gmeasure[groupname])

      def write(gmeasure, groupname):
        dist_tumor = calc_distmap(theta_tumor, ld, 0.5)
        gmeasure.create_dataset(groupname, data = dist_tumor, compression = 9)

      dist_tumor = myutils.hdf_data_caching(read, write, f, ('measurements', 'dist_tumor'), (None,1))

      def write(gmeasure, groupname):
        dist_necro = calc_distmap(theta_necro, ld, tumor_contour_level)
        gmeasure.create_dataset(groupname, data = dist_necro, compression = 9)
        
      dist_necro = myutils.hdf_data_caching(read, write, f, ('measurements', 'dist_necro'), (None,1))

      dist_viabletumor = np.maximum(dist_tumor, -dist_necro)

      return dict(
        phi_cells = phi_cells,
        phi_viabletumor = phi_viabletumor,
        phi_necro = theta_necro,
        dist_viabletumor = dist_viabletumor,
        dist_tumor = dist_tumor,
        dist_necro = dist_necro)

    #####
    if dataname == 'tissue_composition':
      tc = dict(obtain_data('tumor_composition'))
      tc.update(phi_vessels = obtain_data('phi_vessels'),
                dist_vessels = obtain_data('dist_vessels'))
      return tc

    #####
    if dataname == 'iffvelocity':
      v = np.asarray(f['iff/iff_velocity'])
      v = np.rollaxis(v, 3, 0)
      return v # first dimension is the vector component of the velocity

    #####
    if dataname == 'iffvelocity_magnitude':
      v = obtain_data('iffvelocity')
      vv = np.sqrt(v[0,...]**2 + v[1,...]**2 + v[2,...]**2)
      return vv

    #####
    if dataname == 'iffvelocity_outward':
      v = obtain_data('iffvelocity')
      distmap = obtain_data('tumor_composition')['dist_tumor']
      vnorm = np.zeros(distmap.shape, dtype=np.float32)
      dgrad = krebsutils.field_gradient(distmap, spacing=ld.scale)
      for vv, g in zip(v, dgrad):
        vnorm += vv*g
      return vnorm

    #####
    if dataname in ('iff_pressure', 'iff_sources'):
      return f['iff'][dataname]

    #####
    if dataname == 'mask_tumor':
      return np.asarray(f['iff/theta_tumor']) > 0.01

    #####
    if dataname == 'mask_viabletumor':
      def read(gmeasure, dsname):
        return np.asarray(gmeasure[dsname])

      def write(gmeasure, dsname):
        tc = obtain_data('tumor_composition')
        mask = np.asarray(tc['phi_viabletumor']) > 0.01
        gmeasure.create_dataset(dsname, data = mask, compression = 4)

      return myutils.hdf_data_caching(read, write, f, ('measurements', 'mask_viabletumor'), (None,1))


class DataGlobalIff(object):
  keywords = [
    'iff_global'
  ]

  def obtain_data(self, dataman, dataname, *args):
    f, args = args[0], args[1:]
    obtain_data = lambda *args: dataman.obtain_data(args[0], f, args[1:])

    def read(gmeasure, groupname):
      return dict((k, float(v[...])) for (k,v) in gmeasure[groupname].iteritems())

    def write(gmeasure, groupname):
      ld = obtain_data('ld')
      distmap = obtain_data('tumor_composition')['dist_tumor']

      res = {}
    #  e = extractVtkFields.Extractor(f, ['iff_flux_out', 'distmap'])
    #  (surface, _), = e.asVtkDataSets(return_ld_path = True)
    #  surface.GetCellData().SetActiveScalars('distmap')
    #  surface = vtkcommon.vtkCellDataToPointData(surface)
    #  surface = vtkcommon.vtkContour(surface, [0.])
    #  res_cd, res_pd, res_surf = vtkcommon.vtkIntegrateData(surface)
    #  res = dict((surface.GetPointData().GetArray(i).GetName(), dat) for i,dat in enumerate(res_pd))
    #  res['area'] = res_surf
    #  res['iff_flux_out_per_area'] = res['iff_flux_out'] / res['area']
    #  surface = None

      volelem = (ld.scale**3)

      #src = np.asarray(f['iff_sources'])
      #tum = np.asarray(f['theta_tumor'])
      #plus = np.asarray(src > 0, dtype=float)
    #  res['tumor_src_plus'] = np.sum(src*tum*plus)*volelem
    #  res['tumor_src_minus'] = np.sum(src*tum*(1-plus))*volelem
    #  res['ext_src_plus'] = np.sum(src*(1-tum)*plus)*volelem
    #  res['ext_src_minus'] = np.sum(src*(1-tum)*(1-plus))*volelem

      print 'glob, distmap: %f %f, ld.scale %f' % (distmap.min(), distmap.max(), ld.scale)
      src_plus = np.asarray(f['iff/iff_sources_vess_in'])
      src_lymph = np.asarray(f['iff/iff_sources_lymph_out'])
      src_minus_vess = np.asarray(f['iff/iff_sources_vess_out'])
      #theta_tum = np.asarray(f['iff/theta_tumor'])
      src_minus = src_minus_vess+src_lymph
      tum = np.asarray(distmap<-1.*ld.scale, dtype=float)
      res['tumor_volume'] = np.sum(tum)*volelem
      res['sys_volume'] = np.prod(ld.shape)*volelem
      res['tumor_src_plus'] = np.sum(src_plus*tum)*volelem
      res['tumor_src_minus'] = np.sum(src_minus*tum)*volelem
      res['ext_src_plus'] = np.sum(src_plus*(1-tum))*volelem
      res['ext_src_minus'] = np.sum(src_minus*(1-tum))*volelem

      res['src_plus'] = np.sum(src_plus)*volelem
      res['lymph_src_minus'] = np.sum(src_lymph)*volelem
      res['vess_src_minus'] = np.sum(src_minus_vess)*volelem

      res['tumor_lymph_src_minus'] = np.sum(src_lymph*tum)*volelem
      res['ext_lymph_src_minus'] = res['lymph_src_minus']
      res['tumor_vess_src_minus'] = np.sum(src_minus_vess*tum)*volelem
      res['ext_vess_src_minus'] = np.sum(src_minus_vess*(1.-tum))*volelem

      for n in ['tumor_src_plus', 'tumor_src_minus', 'tumor_lymph_src_minus', 'tumor_vess_src_minus']:
        res[n+'_per_vol'] = res[n]/res['tumor_volume']
      for n in ['ext_src_plus', 'ext_src_minus', 'ext_lymph_src_minus', 'ext_vess_src_minus']:
        res[n+'_per_vol'] = res[n]/(res['sys_volume'] - res['tumor_volume']+1.e-13)
      for n in ['src_plus', 'lymph_src_minus', 'vess_src_minus' ]:
        res[n+'_per_vol'] = res[n]/res['sys_volume']
      g = gmeasure.create_group(groupname)
      for k,v in res.iteritems():
        g.create_dataset(k, data = v)

    ret = myutils.hdf_data_caching(read, write, f, ('measurements','global',), (None, 2,))
    a = ret['tumor_src_plus_per_vol']
    b = ret['tumor_src_minus_per_vol']
    ret['tumor_out_per_vol'] = a+b
    ret['tumor_src'] = ret['tumor_src_plus']+ret['tumor_src_minus']
    ret['ext_src'] = ret['ext_src_plus']+ret['ext_src_minus']
    return ret

# this needs some work ... (see analyzeBloodFlow)
#class DataTumorBloodFlow(object):
#  keywords = [
#    'tumor_blood_flow', 'extravasation_vs_blood_flow'
#  ]
#
#  def obtain_data(self, dataman, dataname, *args):
#    f, args = args[0], args[1:]
#    obtain_data = lambda *args: dataman.obtain_data(args[0], f, args[1:])
#
#    if dataname == 'tumor_blood_flow':
#      def read(gmeasure, groupname):
#        return dict((k, float(v[...])) for (k,v) in gmeasure[groupname].iteritems())
#
#      def write(gmeasure, groupname):
#        ld = obtain_data('ld')
#        distmap = obtain_data('tumor_composition')['dist_tumor']
#
#        vessels = krebsutils.read_vesselgraph(f['iff/vessels'], ['flow', 'pressure', 'position', 'flags'])
#        pos = vessels['position']
#        dist = krebsutils.sample_field(pos, distmap, ld, linear_interpolation=True)
#        press = vessels['pressure']
#        flow = vessels['flow']
#        flags = vessels['flags']
#
##        loc_in = []
##        loc_out = []
#
#        total_flow = 0.
#        flow_in, flow_out = 0., 0.
#        for i, (a,b) in enumerate(vessels.edgelist):
#          if not (flags[i] & krebsutils.CIRCULATED): continue
#          if dist[a]<0 and dist[b]>0: # a is in the tumor
#            if press[a]<press[b]:
#              flow_in += flow[i]
##              loc_in.append(np.average((pos[a], pos[b]), axis=0))
#            else:
#              flow_out += flow[i]
##              loc_out.append(np.average((pos[a], pos[b]), axis=0))
#          elif dist[a]>0 and dist[b]<0: # also in tumor
#            if press[a]>press[b]:
#              flow_in += flow[i]
##              loc_in.append(np.average((pos[a], pos[b]), axis=0))
#            else:
#              flow_out += flow[i]
#              #loc_out.append(np.average((pos[a], pos[b]), axis=0))
#          if (flags[i] & krebsutils.BOUNDARY) and (flags[i] & krebsutils.ARTERY):
#            total_flow += flow[i]
#            #loc_out.append(np.average((pos[a], pos[b]), axis=0))
#
#        res = dict(flow_in = flow_in, flow_out = flow_out, total_flow = total_flow)
#
#        g = gmeasure.create_group(groupname)
#        for k,v in res.iteritems():
#          g.create_dataset(k, data = v)
#
##              loc_out.append(np.average((pos[a], pos[b]), axis=0))
#
##        loc_in = np.asarray(loc_in)
##        loc_out = np.asarray(loc_out)
#
##        fig, ax = pyplot.subplots(1,1)
##        ax.plot(loc_in[:,0], loc_in[:,1], lw = 0, ms=3, marker = 'x', color = 'r')
##        ax.plot(loc_out[:,0], loc_out[:,1], lw = 0, ms=3, marker = 'x', color = 'b')
##        pyplot.show()
#
#      fm = myutils.MeasurementFile(f, h5files)
#      ret = myutils.hdf_data_caching(read, write, fm, ('tumor_blood_flow',), (2,))
#      return ret



def averaged(getter, sources, names):
  '''Computes mean of different curves or values. It uses MeanValueArray.avg, i.e. returns the mean of averages'''
  all_res = defaultdict(myutils.MeanValueArray.empty)
  for src in sources:
    for name in names:
      dat = getter(src, name)
      if isinstance(dat, myutils.MeanValueArray):
        all_res[name] += dat.avg
      else:
        all_res[name] += dat
  for k, v in all_res.iteritems():
    all_res[k] = (v.avg, v.std)
  return all_res



class DataGlobalAverageIff(object):
  keywords = [
    'iff_global_average'
  ]

  def obtain_data(self, dataman, dataname, *args):
    files, args = args[0], args[1:]
    datlist = defaultdict(list)
    for f in files:
      d = dataman('iff_global', f)
      for k,v in d.iteritems():
        datlist[k].append(v)
    for k,v in datlist.iteritems():
      datlist[k] = np.average(datlist[k]), np.std(datlist[k])
    return datlist



class DataIffAverages(object):
  keywords = [
    'iff_global_cached_average', 'iff_radial_cached_average'
  ]

  def __init__(self, path):
    self.path = path

  def obtain_data(self, dataman, dataname, *args):
    files = args[0]
    args = args[1:]
    assert isinstance(files, list)

    identifier = myutils.checksum(*list(basename(f.filename) for f in files))

    averaged_data_filename = join(self.path,'measuredglobal.h5')
    with h5files.open(averaged_data_filename, 'a') as favg:
      #####
      if dataname in ('iff_radial_cached_average'):
        def write(gmeasure, groupname):
          # NOTE: averaged returns a dict which maps to tuples of (avg, std), both entries containing arrays
          bins = dataman('iff_radial',files[0],'vs_dr','bins')
          radial = averaged(lambda f,name: dataman('iff_radial', f, 'vs_dr', name), files,
                            ['iff_pressure', 'ivp', 'ivp_minus_ifp', 'iff_velocity_out', 'iff_velocity_mag', 'iff_sources'])

          g = gmeasure.create_group(groupname)
          g.create_dataset('bins', data = bins, compression = 4)
          for k, (a, b) in radial.iteritems():
            h = g.create_group(k)
            h.create_dataset('avg', data = a, compression = 4)
            h.create_dataset('std', data = b, compression = 4)

        def read(gmeasure, groupname):
          g = gmeasure[groupname]
          items = dict(g.iteritems())
          bins = np.asarray(items.pop('bins'))
          data = dict((k, (np.asarray(h['avg']), np.asarray(h['std']))) for k, h in items.iteritems())
          return bins, data

        return myutils.hdf_data_caching(read, write, favg, (dataname, identifier), (0,1))

      if dataname in ('iff_global_cached_average'):
        def write(gmeasure, groupname):
          keys = dataman('iff_global',files[0]).keys()
          data = averaged(lambda f, name: dataman('iff_global',f)[name], files, keys)
          tbl = [ (k, float(v[0][...]), float(v[1][...])) for k,v in data.iteritems() ]
          tbl = np.array(tbl, dtype = [('name', np.str_, 128), ('avg', float), ('std', float)])
          gmeasure.create_dataset(groupname, data = tbl)

        def read(gmeasure, groupname):
          tbl = gmeasure[groupname]
          data = dict((r[0], (r[1], r[2])) for r in tbl)
          return data

        return myutils.hdf_data_caching(read, write, favg, (dataname, identifier), (0,1))



class DataRadialIff(object):
  keywords = [
    'iff_radial',
    'iff_samples',
  ]

  def obtain_data(self, dataman, dataname, *args):
    f, args = args[0], args[1:]
    obtain_data = lambda *args: dataman.obtain_data(args[0], f, args[1:])

    def write(gmeasure, groupname):
      ld = obtain_data('ld')
      gmeasure.try_del('radial')
      gmeasure.try_del('samples')

      vesselgroup = f['iff/vessels']
      if not 'flow' in vesselgroup['edges']:
        def replace_ds(name, data):
          if name in vesselgroup:
            vesselgroup[name][...] = data
          else:
            vesselgroup.create_dataset(name, data = data, compression = 9)
        (pressure, flow, hematocrit, flags) = krebsutils.calc_vessel_hydrodynamics(vesselgroup, return_flags=True)
        vesselgroup.create_dataset('edges/flow', data = flow, compression = 9)
        replace_ds('nodes/pressure', pressure)
        replace_ds('edges/flags', flags)

      # read circulated vessels
      vessels = krebsutils.read_vesselgraph(f['iff/vessels'], ['flow','pressure','conductivity', 'position', 'flags', 'wall_conductivity', 'radius'])
      vessels = vessels.get_filtered(edge_indices = np.nonzero(vessels.edges['flags'] & krebsutils.CIRCULATED)[0])
      q = vessels['pressure']
      print q.min(), q.max(), np.average(q)

      # distance maps
      radialmap = krebsutils.make_radial_field(ld)
      distmap = obtain_data('tumor_composition')['dist_tumor']
      map_r = dict(vs_r = radialmap, vs_dr = distmap)
      #
      smpl_pos  = plotVessels.generate_samples(vessels, 'position', 'nodes', sample_length)
      smpl_dist = krebsutils.sample_field(smpl_pos, distmap, ld, linear_interpolation=True, extrapolation_value = far)
      smpl_rad  = krebsutils.sample_field(smpl_pos, radialmap, ld, linear_interpolation=True, extrapolation_value = far)
      smpl_r = dict(vs_r = smpl_rad, vs_dr = smpl_dist)

      gmeasure.require_group('radial').require_group('vs_r').create_dataset('bins', data = np.average((bins_rad[:-1],bins_rad[1:]), axis=0))
      gmeasure.require_group('radial').require_group('vs_dr').create_dataset('bins', data = np.average((bins_dist[:-1],bins_dist[1:]), axis=0))

      samplesdata = dict()
      sample_mask0 = np.logical_and(smpl_dist<0, np.random.random_sample(smpl_dist.shape)<stored_sample_fraction_tumor)
      sample_mask1 = np.logical_and(smpl_dist>=0, np.random.random_sample(smpl_dist.shape)<stored_sample_fraction_normal)
      sample_mask = np.logical_or(sample_mask0, sample_mask1)
      del sample_mask0, sample_mask1

      def store_samples(name, smpl):
        samplesdata[name] = downcast(smpl[sample_mask])

      store_samples('r', smpl_rad)
      store_samples('dr', smpl_dist)

    #  def store_histogram_data(flavor, name, (s_avg, s_std, s_sqr)):
    #    dstdata[flavor][name] = dict(avg = s_avg, std = s_std, sqr = s_sqr)
      def store_histogram_data(flavor, name, d):
        d.write(gmeasure['radial'][flavor], name)

      def store_vesseldata(name, samples):
        for flavor in ['vs_r', 'vs_dr']:
          #d = myutils.scatter_histogram(smpl_r[flavor], samples, bins_r[flavor], 1.)
          #store_histogram_data(flavor, name, d)
          d = myutils.MeanValueArray.fromHistogram1d(bins_r[flavor], smpl_r[flavor], samples)
          store_histogram_data(flavor, name, d)
        store_samples(name, samples)

      def store_fielddata(name, samples, mask = None):
        for flavor in ['vs_r', 'vs_dr']:
          x, y = map_r[flavor], samples
          if mask is not None:
            x, y = x[mask], y[mask]
          #d = myutils.scatter_histogram(x.ravel(), y.ravel(), bins_r[flavor], 1.)
          d = myutils.MeanValueArray.fromHistogram1d(bins_r[flavor], x.ravel(), y.ravel())
          store_histogram_data(flavor, name, d)

      g = f['iff']
      ifpressure = np.asarray(g['iff_pressure'])
      store_fielddata('iff_pressure', ifpressure)

      smpl_ifp = krebsutils.sample_field(smpl_pos.astype(np.float32), ifpressure.astype(np.float32), ld, linear_interpolation=True)
      del ifpressure

      store_fielddata('iff_sources', np.asarray(g['iff_sources']))
      store_fielddata('iff_sources_vess_out', np.asarray(g['iff_sources_vess_out']))
      store_fielddata('iff_sources_vess_in', np.asarray(g['iff_sources_vess_in']))
      store_fielddata('iff_sources_lymph_out', np.asarray(g['iff_sources_lymph_out']))

      smpl_pressure = plotVessels.generate_samples(vessels, 'pressure', 'nodes', sample_length)
      store_vesseldata('ivp', smpl_pressure)

      store_vesseldata('ivp_minus_ifp', smpl_pressure - smpl_ifp)
      #del smpl_pressure, smpl_ifp, smpl_pos

      smpl_vess_wcond = plotVessels.generate_samples(vessels, 'wall_conductivity', 'edges', sample_length)
      smpl_vess_radii = plotVessels.generate_samples(vessels, 'radius', 'edges', sample_length)
      smpl_vess_wcond *= math.pi*2.*smpl_vess_radii
      smpl_vess_cond = plotVessels.generate_samples(vessels, 'conductivity', 'edges', sample_length)
      smpl_vess_flow = plotVessels.generate_samples(vessels, 'flow', 'edges', sample_length)
      smpl_vess_outflux = smpl_vess_wcond*(smpl_pressure-smpl_ifp)
      smpl_outflux_per_flow = smpl_vess_outflux / smpl_vess_flow
      smpl_vess_lambda = np.sqrt(smpl_vess_cond / smpl_vess_wcond)

      store_vesseldata('v_radii', smpl_vess_radii)
      store_vesseldata('v_cond', smpl_vess_cond)
      store_vesseldata('v_wcond', smpl_vess_wcond)
      store_vesseldata('v_flow', smpl_vess_flow)
      store_vesseldata('v_lambda', smpl_vess_lambda)
      store_vesseldata('v_outflux', smpl_vess_outflux)
      store_vesseldata('v_outflux_per_flow', smpl_outflux_per_flow)

      del smpl_ifp, smpl_pressure
      del smpl_vess_cond, smpl_vess_flow, smpl_vess_outflux, smpl_outflux_per_flow

      iff_velocity_out = dataman('iffvelocity_outward', f)
      iff_velocity_mag = dataman('iffvelocity_magnitude', f)
      store_fielddata('iff_velocity_out', iff_velocity_out)
      store_fielddata('iff_velocity_mag', iff_velocity_mag)

      del iff_velocity_out, iff_velocity_mag
      myutils.hdf_write_dict_hierarchy(gmeasure, 'samples', samplesdata)

    def read(gmeasure, groupname):
      q = dataname[len('iff_'):]
      if q == 'radial':
        flavor, subname = args[0], args[1]
        if subname != 'bins':
          return myutils.MeanValueArray.read(gmeasure[q][flavor], subname)
        else:
          return np.asarray(gmeasure[q][flavor][subname])
      else:
        r = args[1]
        rnd = np.random.mtrand.RandomState(abs(hash(f.filename)) & (0xffffffff))
        arr = np.asarray(gmeasure[q][args[0]])
        #arr = arr[rnd.randint(0, high = len(arr), size = len(arr)*r)]
        arr = arr[rnd.randint(0, high = len(arr))]
        return arr

    return myutils.hdf_data_caching(read, write, f, ('measurements', 'radial',), (None, 1,))


@myutils.UsesDataManager
def ComputeIfpVsIffCorrelationDataLocal(dataman, iff_file):
  boundary_region_distances = -100., 100.  
  tumor_center_distance = 500.
  iff_bins = np.linspace(0., 5., 100)
  ifp_bins = np.linspace(0., 20., 20)
  
  def read(gmeasure, gname):
    return gmeasure[gname]

  def write(gmeasure, gname):
    print('compute Ifp vs Iff correlation for %s' % iff_file)
    ld = krebsutils.read_lattice_data_from_hdf(iff_file['field_ld'])

    iff_radial_field = dataman.obtain_data('iffvelocity_outward', iff_file)
    iff_radial_field = np.asarray(iff_radial_field).ravel()  # because masking works only (?) with 1d arrays
    distmap          = dataman.obtain_data('tumor_composition', iff_file)['dist_tumor']
    mask             = np.logical_and(boundary_region_distances[0] < distmap, distmap < boundary_region_distances[1])
    mask             = mask.ravel()
    iff_radial_field = iff_radial_field[mask]
    iff_histo        = myutils.MeanValueArray.fromHistogram1d(iff_bins, iff_radial_field, np.ones_like(iff_radial_field))
    
    ifp_field        = dataman.obtain_data('iff_pressure', iff_file)
    ifp_field        = np.asarray(ifp_field).ravel()
    distmap          = dataman.obtain_data('distance_from_center_distribution', ld)
    distmap          = distmap.ravel()
    mask             = distmap < tumor_center_distance
    ifp_field        = ifp_field[mask]
    ifp_histo        = myutils.MeanValueArray.fromHistogram1d(ifp_bins, ifp_field, np.ones_like(ifp_field))
    
    gmeasure = gmeasure.create_group(gname)
    g = iff_histo.write(gmeasure, 'iff_histo')
    g.attrs['mean'] = np.average(iff_radial_field)
    g.attrs['std']  = np.std(iff_radial_field)
    g = ifp_histo.write(gmeasure, 'ifp_histo')
    g.attrs['mean'] = np.average(ifp_field)
    g.attrs['std']  = np.std(ifp_field)  

  ret = myutils.hdf_data_caching(read, write, iff_file, ('measurements', 'ifp_vs_iff_correlation',), (None, 1,))  
  return ret



def ComputeIfpVsIffCorrelationData(dataman, global_data_filename, filenamelist):
  '''note: function closes all files it touches'''
  version_identifier = myutils.checksum(filenamelist, 1)

  def read(gmeasure, gname):
    # return those two datasets in a dict with their names, as numpy arrays
    return dict(map(lambda s: (s,np.asarray(gmeasure[gname][s])), ['iff_means', 'ifp_means']))

  def write(gmeasure, gname):
    result = []  
    for filename in filenamelist:
      with h5files.open(filename, 'a') as f:
        resultgroup = ComputeIfpVsIffCorrelationDataLocal(dataman, f)
        iff_mean = resultgroup['iff_histo'].attrs['mean']
        #iff_std  = resultgroup['iff_histo'].attrs['std']
        ifp_mean = resultgroup['ifp_histo'].attrs['mean']
        #ifp_std  = resultgroup['ifp_histo'].attrs['std']
        del resultgroup
        result.append((iff_mean, ifp_mean))
    
    iff_means, ifp_means = zip(*result)
    gmeasure = gmeasure.create_group(gname)
    gmeasure.create_dataset('filenames', data = filenames)
    gmeasure.create_dataset('iff_means', data = iff_means)
    gmeasure.create_dataset('ifp_means', data = ifp_means)
    
  with h5files.open(global_data_filename, 'a') as global_data_file:
    ret = myutils.hdf_data_caching(read, write, global_data_file, ('ifp_vs_iff_correlation',), (version_identifier,))    
    return ret
      

##############################################################################
### plotting from down here
##############################################################################


def plot_global_header(files, dataman, pdfpages):
  results = []
  for f in files:
    res = dataman('iff_global', f)
    # tumor blood flow disable because blood flow must be computed for it
    #res.update(dataman('tumor_blood_flow', f))
    #res['extravasated_fraction'] = res['tumor_src_plus']/res['flow_in']
    results.append(res)
  results = myutils.zipListOfDicts(results)
  for k, v in results.iteritems():
    results[k] = np.average(v), np.std(v)

  as_str = lambda name, scale: r'$%s \pm %s$' % (lf2s(results[name][0]*scale), lf2s(results[name][1]*scale))
  vol_p_sec = r'[$\mu m^3 / s$]'

  ld = dataman('ld', files[0])
  s = krebsutils.LatticeDataGetWorldSize(ld)*1.e-3
  ss = [ r'system size: $%s \times %s\times %s$ [$mm$]' % (lf2s(s[0]), lf2s(s[1]), lf2s(s[2])),
         r'tumor volume: %s [$\mu m^3$]' % (as_str('tumor_volume',1.)),
         '','total fluid fluxes',
         r'tumor iff sources (total) = %s %s' % (as_str('tumor_src', 1), vol_p_sec) ,
         r'tumor iff sources (in)    = %s %s' % (as_str('tumor_src_plus',1), vol_p_sec),
         r'tumor iff sources (out)   = %s %s' % (as_str('tumor_src_minus',1), vol_p_sec),
         r'tissue iff sources (total) = %s %s' % (as_str('ext_src', 1), vol_p_sec) ,
         r'tissue iff sources (in)    = %s %s' % (as_str('ext_src_plus',1), vol_p_sec),
         r'tissue iff sources (out)   = %s %s' % (as_str('ext_src_minus',1), vol_p_sec) ]
  if 'flow_in' in results:
    ss += [
         r'blood flow, in = %s %s' % (as_str('flow_in',1), vol_p_sec),
         r'blood flow, out = %s %s' % (as_str('flow_out',1), vol_p_sec),
         r'extravasated fluid fraction: %s %s' % (as_str('extravasated_fraction',1), vol_p_sec) ]

  if 0: # enable again when tumor_blood_flow works
    tumor_blood_flow = dataman('tumor_blood_flow', f)
    blood_in, blood_out = tumor_blood_flow['flow_in'], tumor_blood_flow['flow_out']
    s5 = r'blood flow, in = $%s$ [$\mu m^3 /s$], out = $%s$ [$\mu m^3 /s$]' % (lf2s(blood_in), lf2s(blood_out))
    extravasation_vs_blood_ratio = res['tumor_src_plus'] / blood_in
    s6 = r'extravasated fluid fraction: $%s$' % lf2s(extravasation_vs_blood_ratio)
  fig = pyplot.figure(figsize=(5,5))
  fig.text(0.05, 0.5, '\n'.join(ss))
  pdfpages.savefig(fig)


def plot_radial(files, dataman, pdfpages):

  radialgetter = lambda f,name: dataman('iff_radial', f, 'vs_dr', name)

  radial = averaged(
    radialgetter,
    files,
    ['iff_pressure', 'ivp', 'ivp_minus_ifp', 'iff_velocity_out', 'iff_velocity_mag', 'v_wcond', 'iff_sources', 'iff_sources_vess_in', 'iff_sources_vess_out', 'iff_sources_lymph_out'])
  bins = (1./1000.) * dataman('iff_radial', files[0], 'vs_dr', 'bins')
  xlim = -1.5, 1.0
  mask = np.logical_and(bins < xlim[1], bins >= xlim[0])

  def plot(name, **kwargs):
    f = kwargs.pop('value_prefactor', 1.)
    avg, std = radial[name]
    avg, std = avg[mask], std[mask]
    kwargs.update(every = 5, marker = kwargs.pop('marker', 's'), color = kwargs.pop('color', 'k'))
    ret = mpl_utils.errorbar(ax, bins[mask], f*avg, yerr = f*std, **kwargs)
    return ret

  fig, axes = pyplot.subplots(3, 1, figsize = (mastersize[0]*0.5, mastersize[0]*3*0.25))
  axes = axes.ravel()
  mpl_utils.subplots_adjust_abs(fig, left = mpl_utils.mm_to_inch*20,
                                right = -mpl_utils.mm_to_inch*10,
                                top  = -mpl_utils.mm_to_inch*5,
                                bottom = mpl_utils.mm_to_inch*10,
                                hspace = mpl_utils.mm_to_inch*30,)

  ax = axes[0]
  ax.set(xticklabels=[])
  ax.set(ylabel = r'[kPa]') #, xlabel = r'$\theta$ [mm]', title = 'pressure')
  plot('iff_pressure', label = r'$p_i$', marker='o', color = 'r')
  plot('ivp', label = '$p_v$', marker = 's', color = 'k')
  plot('ivp_minus_ifp', label = '$p_v - p_i$', marker='>', color = 'b')
  text2(ax, fig_numbering[0])
  ax.legend()

  ax = axes[1]
  ax.set(ylabel = ur'[\u03BCm/s]') #, xlabel = r'$\theta$ [mm]', title = 'velocity')
  ax.set(xticklabels=[])
  plot('iff_velocity_out', label = r'$v_{||}$', marker = 'o', color = 'r')
  plot('iff_velocity_mag', label = r'$|v|$', marker = 's', color = 'b')
  text2(ax, fig_numbering[1])
  ax.legend()

  ax = axes[2]
  ax.set(ylabel = r' [s$^{-1}$]', xlabel = r'$\theta$ [mm]') #, title = 'sources')
  plot('iff_sources', label = r'$%s_l\,\times 10^3$' % LF.source, value_prefactor = 1.e3, marker='o', color = 'r')
  plot('iff_sources_vess_in', label = r'$%s_{lv}^+\,\times 10^3$' % LF.source, value_prefactor = 1.e3, marker='s', color = 'g')
  plot('iff_sources_vess_out', label = r'$%s_{lv}^-\,\times 10^4$' % LF.source, value_prefactor = 1.e4, marker='<', color = 'b')
  plot('iff_sources_lymph_out', label = r'$%s_{ll}^-\,\times 10^3$' % LF.source, value_prefactor = 1.e3, marker='d', color = 'm')
  text2(ax, fig_numbering[2])
  ax.legend(loc = 4)

  for ax in axes:
    ax.set_xlim(*xlim)
    ax.grid(linestyle=':', linewidth=0.5, color=gridcolor)
    mpl_utils.add_crosshair(ax, (0,0), color = gridcolor)

  pdfpages.savefig(fig)

  if 1:
    yscale = 'log'
    samplegetter = lambda f,name: dataman('iff_samples', f, name, 0.05)
    samples = defaultdict(list)
    for f in files:
      dr = samplegetter(f, 'dr')
      mask = np.logical_and(dr < 1.e6, dr > -1.e6) # mask entries which are outside of the numerical grid
      samples['dr'].append(dr[mask])
      for name in ['v_wcond', 'v_flow', 'v_outflux', 'v_outflux_per_flow', 'v_lambda']:
        smpl = samplegetter(f, name)
        samples[name].append(smpl[mask])
    for k,v in samples.iteritems():
      samples[k] = np.concatenate(v)
    #radial = averaged(radialgetter, files, ['v_wcond', 'v_flow', 'v_outflux', 'v_outflux_per_flow', 'v_lambda'])

    def mkfig(nrows, ncols):
      fig, axes = mpl_utils.subplots_abs_mm((mastersize[0]/mpl_utils.mm_to_inch*0.5*ncols,
                                   mastersize[0]/mpl_utils.mm_to_inch*0.2*nrows),
                                            nrows, ncols,
                                            10, 20, a4size[0]/mpl_utils.mm_to_inch*0.38, a4size[0]/mpl_utils.mm_to_inch*0.16, 15, 5)
      return fig, axes

    def plot(ax, name, **scatter_kwargs):
      f = scatter_kwargs.pop('value_prefactor', 1.)
      y = samples[name] * f
      x = samples['dr'] * 1.e-3
      ret = ax.scatter(x, y, s=0.5, alpha = 0.6, color = scatter_kwargs.pop('color', 'k'), rasterized=True, edgecolors='none', **scatter_kwargs)
      return ret

    fig, axes = mkfig(4, 1)
    axes = axes.ravel()

    ax = axes[0]
    ax.set(ylabel=ur'[$\u03BCm^3 / \u03BCm\,kPa\,s$]', yscale = yscale) #title=r'transmural flow coefficient'
    ax.set(xticklabels=[])
    plot(ax, 'v_wcond')
    text1(ax, LF.math(LF.transmural_flow_coeff))
    text2(ax, fig_numbering[0])
    #ax.legend()

    ax = axes[1]
    ax.set(ylabel=ur'[$\u03BCm^3 / \u03BCm\,s$]', yscale = yscale) #title=r'transmural flux per length'
    ax.set(xticklabels=[])
    plot(ax, 'v_outflux')
    text1(ax, LF.math(LF.transmural_flow_per_length))
    text2(ax, fig_numbering[1])
    ax.set(ylim = (1.e-2, 1.e3))
    #ax.legend()

    ax = axes[2]
    ax.set(xticklabels=[])
    ax.set(ylabel=ur'[$\u03BCm^3/s$]', yscale = yscale)  #xlabel=r'$\theta$ [mm]' title=r'vessel flow rate'
    plot(ax, 'v_flow')
    text1(ax, LF.math(LF.vessel_flow_rate))
    ax.set(ylim = (1.e0, 1.e10))
    text2(ax, fig_numbering[2])
    #ax.legend()

    ax = axes[3]
    ax.set(xlabel=r'$\theta$ [mm]', ylabel=ur'[\u03BCm]', yscale = yscale) #title=r'$\lambda$'
    plot(ax, 'v_lambda') #, value_prefactor = 1.e-3)
    text1(ax, LF.math(LF.relative_extravasation_flow))
    text2(ax, fig_numbering[3])
    #ax.legend()

  #  ax.set(title=r'relative transmural loss', xlabel=r'$\Delta r$ [mm]', ylabel=r'[$\mu m^{-1}$]', yscale = 'log')
  #  plot('v_outflux_per_flow')

    for ax in axes:
      ax.set_xlim(*xlim)
      ax.grid(linestyle=':', linewidth=0.5, color=gridcolor)
      mpl_utils.add_crosshair(ax, (0,None), color = gridcolor)

    pdfpages.savefig(fig)



def plot_single_images2(f, dataman, pdfpages):
    ## first get lattice data
    ld = dataman('ld', f)

    #### vessel volume fraction #####
#    vessels = krebsutils.read_vesselgraph(f['vessels'], ['radius','position', 'flags'])
#    mask = np.nonzero(vessels.edges['flags'] & krebsutils.CIRCULATED)[0]
#    vessels = vessels.get_filtered(edge_indices = mask)
#    vessels.nodes['position']
#    vessels.edges['radius']
#    vessel_fraction = krebsutils.make_vessel_volume_fraction_field(
#      vessels.nodes['position'],
#      vessels.edgelist,
#      vessels.edges['radius'],
#      ld)
#    del vessels
#    ld.Rescale(1.e-3)
    #vessel_fraction = dataman('phi_vessels', f)

    tc = dataman('tissue_composition', f)
    imgs = dict(
      iff_pressure = dataman('iff_pressure', f),
      iff_sources = dataman('iff_sources', f),
      #theta_tumor = tc['phi_cells'],
      #distmap = tc['dist_tumor'],
      #theta_necro = tc['phi_necro'],
      dist_viabletumor = tc['dist_viabletumor'],
      iff_velocity_mag = dataman('iffvelocity_magnitude', f),
      iff_velocity_out = dataman('iffvelocity_outward', f),
      vesselfraction = tc['phi_vessels'],
      iff_velocity = dataman('iffvelocity', f)[0]
    )
    for k,v in imgs.iteritems():
      imgs[k] = imslice(v)

    def mkfig(nrows, ncols): # sizes set up for image plots
      fig, axes = mpl_utils.subplots_abs_mm((mastersize[0]/mpl_utils.mm_to_inch*0.33*ncols,
                                   mastersize[0]/mpl_utils.mm_to_inch*0.28*nrows),
                                            nrows, ncols,
                                            2, 5, a4size[0]/mpl_utils.mm_to_inch*0.25, a4size[0]/mpl_utils.mm_to_inch*0.25, 15, 4)
      mpl_utils.remove_frame(*axes.ravel())
      return fig, axes


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

    def textright(ax, txt):
      ax.text(0.5, 1.03, txt, ha = "center", transform = ax.transAxes)

    def textleft(ax, txt):
      ax.text(0., 1.03, txt, ha = "left", transform = ax.transAxes)

    def plt_phi_vessels(ax):
      ax.set(title = '$\phi_v$')
      p1 = imshow(ax, imslice(tc['phi_vessels']), ld, vmin = 0., vmax = 1., cmap=CM.grey)
      contour(ax, imslice(tc['dist_viabletumor']), ld, levels=[0], colors='w')
      colorbar(fig, ax, p1)
      mpl_utils.add_sizebar(ax)
      textleft(ax, fig_numbering[0])

#r'v_x' r'$v_{||}$' r'$|v|$' $p_i$ $\Gamma_l$

    def plt_ifp(ax):
      ax.set(title=ur'$p_i$')
      p = imshow_('iff_pressure', matplotlib.cm.jet, None, 1.)
      _, cax = colorbar(fig, ax, p)
      contour_('dist_viabletumor')
      textright(cax, ur'[kPa]')
      textleft(ax, fig_numbering[1])

    def plt_ifv(ax):
      ax.set(title=ur'$v_x$')
      p = imshow_('iff_velocity', matplotlib.cm.seismic, 'zero-centered', 1.)
      contour_('dist_viabletumor')
      _, cax = colorbar(fig, ax, p)
      textright(cax, ur'[\u03BCm/s]')
      textleft(ax, fig_numbering[3])

    def plt_ifsrc(ax):
      ax.set(title=ur'$%s_l$' % LF.source)
      p = imshow_('iff_sources', matplotlib.cm.seismic, 'zero-centered', 1.e3)
      contour_('dist_viabletumor')
      _, cax = colorbar(fig, ax, p)
      textright(cax, ur'[$10^3$ \u03BCm$^3$/\u03BCm$^3$s]')
      textleft(ax, fig_numbering[2])

    def plt_ifvmag(ax):
      ax.set(title=ur'$|v|$')
      p = imshow_('iff_velocity_mag', matplotlib.cm.gray, None, 1.)
      contour_('dist_viabletumor')
      _, cax = colorbar(fig, ax, p)
      textright(cax, ur'[\u03BCm/s]')
      textleft(ax, fig_numbering[4])

    def plt_ifvout(ax):
      ax.set(title=ur'$v_{||}$')
      p = imshow_('iff_velocity_out', matplotlib.cm.gray, None, 1.)
      contour_('dist_viabletumor')
      _, cax = colorbar(fig, ax, p)
      textright(cax, ur'[\u03BCm/s]')
      textleft(ax, fig_numbering[5])

    fig, axes = mkfig(2,3)
    for func, ax in zip([plt_phi_vessels, plt_ifp, plt_ifsrc, plt_ifv, plt_ifvmag, plt_ifvout],
                        axes.ravel()):
      func(ax)
    pdfpages.savefig(fig)





def fit_resist(x_in, y_in, initial_params, weights):
  """ fit resistor chain model to total fluxes
  """
  from scipy.optimize import leastsq

  def func(p, x): # the fit function
    return p[0]/(np.power(x, -1.)+1./p[1]) + p[2]

  def do_the_fit(x, y, initial_params):
    objective_func = lambda p, x, y: (func(p,x)-y)*weights
    p, success = leastsq(objective_func, initial_params, args=(x, y))
    return p, success

  p, success = do_the_fit(x_in, y_in, initial_params)
  if not success:
    raise RuntimeError("fitting failed")

  ret_func = lambda x: func(p, x)
  return p, ret_func





def plotComparisonsForIffPaper(pdfpages, path):
  g = lambda fn: glob.glob(os.path.join(path, fn))
  allfiles = dict(
    default = g('iff_default_*.h5'),
    vA01 = g('iff_variantA01_*.h5'),
    vA02 = g('iff_variantA02_*.h5'),
    vA03 = g('iff_variantA03_*.h5'),
    vA04 = g('iff_variantA04_*.h5'),
    vA05 = g('iff_variantA05_*.h5'),
    #vA06 = g('measure_iff_variantA06_*.h5'),
    vA07 = g('iff_variantA07_*.h5'),
    vA08 = g('iff_variantA08_*.h5'),
    vA09 = g('iff_variantA09_*.h5'),
    vA10 = g('iff_variantA10_*.h5'),
    vA11 = g('iff_variantA11_*.h5'),
    vA12 = g('iff_variantA12_*.h5'),

    vB01 = g('iff_variantB01_*.h5'),
    vB02 = g('iff_variantB02_*.h5'),
    vB03 = g('iff_variantB03_*.h5'),
    vB04 = g('iff_variantB04_*.h5'),
    vB05 = g('iff_variantB05_*.h5'),
    vB06 = g('iff_variantB06_*.h5'),
    vB07 = g('iff_variantB07_*.h5'),
    vB08 = g('iff_variantB08_*.h5'),
    vB09 = g('iff_variantB09_*.h5'),

    vC01 = g('iff_variantC01_*.h5'),
    vC02 = g('iff_variantC02_*.h5'),
    vC03 = g('iff_variantC03_*.h5'),
    vC04 = g('iff_variantC04_*.h5'),
    vC05 = g('iff_variantC05_*.h5'),
    vC06 = g('iff_variantC06_*.h5'),

    vD01 = g('iff_variantD01_*.h5'),
    vD02 = g('iff_variantD02_*.h5'),
    vD03 = g('iff_variantD03_*.h5'),
    vD04 = g('iff_variantD04_*.h5'),
    vD05 = g('iff_variantD05_*.h5'),
    vD06 = g('iff_variantD06_*.h5'),
  )
  for i in range(1, 7):
    allfiles['vE%02i' %i] = g('iff_variantE%02i_*.h5' % i)

  for k, v in allfiles.iteritems():
    assert v, ("missing files for case %s" % k)

  dataman = myutils.DataManager(100, [DataTissue(), DataGlobalIff(), DataRadialIff(), DataIffAverages(path), DataTumorBloodFlow()])

  allavg = {}
  for k, v in allfiles.iteritems():
    files = [ h5py.File(fn, 'r') for fn in v ]
    bins, radial = dataman('iff_radial_cached_average', files)
    data = dataman('iff_global_cached_average', files)
    for f in files: f.close()
    a, da = data['tumor_src_plus_per_vol']
    b, db = data['tumor_src_minus_per_vol']
    data['tumor_out_per_vol'] = (a+b, da+db)
    allavg[k] = (data, bins, radial)


  run_labels = [
    #('v16', r'$\times 32$', 32),
    ('vA01', r'$\times 10$', 10),
    ('vA05', r'$\times 4$',  4),
    ('vA02', r'$\times 2$',  2),
    ('default', r'$\times 1$', 1),
    ('vA03', r'$\times \frac{1}{2}$', 0.5),
    ('vA04', r'$\times \frac{1}{10}$', 0.1),
    ('vA07', r'$\times \frac{1}{20}$', 0.05),
    ('vA08', r'$\times \frac{1}{100}$', 0.01),
    ('vA09', r'$\times \frac{1}{4}$', 0.25),
    ('vA10', r'$\times \frac{1}{200}$', 1./200.),
    ('vA11', r'$\times \frac{1}{500}$', 1./500.),
    ('vA12', r'$\times \frac{1}{1000}$', 1./1000.),

    ('vB07', r'$\times 100$', 100),
    ('vB06', r'$\times 50$', 50),
    ('vB01', r'$\times 10$', 10),
    ('vB05', r'$\times 5$', 5),
    ('vB02', r'$\times 2$', 2),
    ('default', r'$\times 1$', 1),
    ('vB03', r'$\times \frac{1}{2}$', .5),
    ('vB04', r'$\times \frac{1}{10}$', .1),
    ('vB08', r'$\times \frac{1}{20}$', .05),
    ('vB09', r'$\times \frac{1}{100}$', .01),

    ('vC01', r'$\times \frac{1}{100}$', 1./100.),
    ('vC02', r'$\times \frac{1}{50}$', 1./50.),
    ('vC03', r'$\times \frac{1}{20}$', 1./20.),
    ('vC04', r'$\times \frac{1}{10}$', 1./10.),
    ('vC05', r'$\times \frac{1}{5}$', 1./5.),
    ('vC06', r'$\times \frac{1}{2}$', 1./2.),

    ('vD01', r'$\times 100$', 100.),
    ('vD02', r'$\times 50$', 50.),
    ('vD03', r'$\times 20$', 20.),
    ('vD04', r'$\times 10$', 10.),
    ('vD05', r'$\times 5$', 5.),
    ('vD06', r'$\times 2$', 2.),

    ('vE01', r'$\times 100$', 100.),
    ('vE02', r'$\times 10$', 10.),
    ('vE03', r'$\times 50$', 50.),
    ('vE04', r'$\times 20$', 20.),
    ('vE05', r'$\times 5$', 5.),
    ('vE06', r'$\times 2$', 2.),
  ]
  run_labels = dict((x[0], x) for x in run_labels)

  def plotradial(ax, runname, dataname, **kwargs):
    xlim = -1/5, 1.
    f = kwargs.pop('value_prefactor', 1.)
    _, bins, radial = allavg[runname]
    avg, std = radial[dataname]
    bins = bins * 1./1000.
    mask = np.logical_and(bins < xlim[1], bins >= -1.5)
    avg, std = avg[mask], std[mask]
    return ax.plot(bins[mask], f*avg, **kwargs)

  #####################################################
  #####################################################
  def plot_pressure_velocity_etc_radial(title, runs):
    def textleft(ax, txt):
      ax.text(0.01, 0.9, txt, ha = "left", transform = ax.transAxes)

    if 0:
      fig, axes = pyplot.subplots(3, 1, figsize = (mastersize[0]*0.5, mastersize[0]*2*0.25))
      axes = axes.ravel()
      mpl_utils.subplots_adjust_abs(fig, left = mpl_utils.mm_to_inch*20,
                                    right = -mpl_utils.mm_to_inch*10,
                                    top  = -mpl_utils.mm_to_inch*5,
                                    bottom = mpl_utils.mm_to_inch*10,
                                    hspace = mpl_utils.mm_to_inch*30,)
    else:
      fig, axes = pyplot.subplots(2, 1, figsize = (6.83/3., 6.83/2.))
      axes = axes.ravel()
      mpl_utils.subplots_adjust_abs(fig,
                                    left = mpl_utils.mm_to_inch*11,
                                    right = -mpl_utils.mm_to_inch*4,
                                    top  = -mpl_utils.mm_to_inch*5,
                                    bottom = mpl_utils.mm_to_inch*10,
                                    hspace = mpl_utils.mm_to_inch*30,)

    fig.suptitle(title)

    ax = axes[0]
    ax.set(ylabel = ur' [kPa]') # xlabel = r'$\theta$ [mm]'
    for name, label, _ in runs:
      plotradial(ax, name, 'iff_pressure', label = label)
    plotradial(ax, 'default', 'ivp', label = 'vessels', color = 'k')

    #text2(ax, fig_numbering[0])
    #text1(ax, ur'$p_i$')


    ax = axes[1]
    ax.set(ylabel = ur'[\u03BCm/s]', xlabel = r'$\theta$ [mm]')
    for name, label, _ in runs:
      p = plotradial(ax, name, 'iff_velocity_out', label = label)
      plotradial(ax, name, 'iff_velocity_mag', color = p[0].get_color(), linestyle=':')
    #text2(ax, fig_numbering[1])
    #text1(ax, ur'$v_{||}$ and $|v|$')
    ax.legend()

    if 0:
      ax = axes[2]
      ax.set(ylabel = ur'$\times 10^3$ [s$^{-1}$]', xlabel = r'$\theta$ [mm]')
      for name, label, _ in runs:
        plotradial(ax, name, 'iff_sources', label = label, value_prefactor = 1.e3)
      #text2(ax, fig_numbering[2])
      text1(ax, LF.math(LF.source+'_l'))

    for ax in axes[0:1]: ax.set(xticklabels = [])
    for ax in axes:
      mpl_utils.add_crosshair(ax, (0,0), color = gridcolor)
      ax.grid(linestyle=':', linewidth=0.5, color= gridcolor)

    pdfpages.savefig(fig)

  #####################################################

  runs = [ run_labels[x] for x in ['vA01', 'default', 'vA04', 'vA08', 'vA12']] # 'vA07', 'vA08']]
  plot_pressure_velocity_etc_radial('leakynesses', runs)

  #############################################

  runs = [ run_labels[x] for x in ['vB07', 'vB01', 'default', 'vB04', 'vB09' ]]
  plot_pressure_velocity_etc_radial('lymphatic wall permeabilities', runs)

  #############################################

  runs = [ run_labels[x] for x in ['default', 'vC06', 'vC05', 'vC03', 'vC01' ]]
  plot_pressure_velocity_etc_radial('tumor lymphatics', runs)

  #############################################

  runs = [ run_labels[x] for x in ['default', 'vE02', 'vE01']]
  plot_pressure_velocity_etc_radial('tumor conductivity', runs)

  #############################################
  ##### global average curves  ##########3

  def plot_global_single_curve(ax, dataname, x, y, **kwargs):
    x = np.asarray(x)
    factor = kwargs.pop('value_scale', 1.)
    dofit = kwargs.pop('dofit', False)
    #y0, y0std = np.asarray(allavg['default'][0][dataname])
    y = np.asarray([allavg[name][0][dataname] for name in y])
    y, ystd = y.transpose()
    #y /= y0
    #ystd /= y0
    #y -= 1
    y *= factor
    ystd *= factor
    y = np.abs(y)
    label = kwargs.pop('label', '')
    label += ur'$\,\times %s$' % myutils.f2s(factor, latex=True, exponential=True)
    kwargs.update(label = label)
    mpl_utils.errorbar(ax, x, y, yerr = ystd, **kwargs)
#    if dofit:
#      p, func = fit_resist(x,y, (1.,1., 0.), weights = 1./x)
#      x_fitted = np.linspace(np.amin(x), np.amax(x), 100)
#      y_fitted = func(x_fitted)
#      ax.plot(x_fitted, y_fitted)

  props = {
    #'iff_flux_out_per_area' :  dict(color = 'r', label = r'(i)     boundary', marker = 'o'), #, value_scale = 0.001),
    'tumor_out_per_vol'             :  dict(color = 'r', label = r'$%s_{l}$' % LF.source, marker = 'o'),
    'tumor_src_plus_per_vol'        :  dict(color = 'k', label = r'$%s_{l}^+$' % LF.source, marker = '<'),
    'tumor_src_minus_per_vol'       :  dict(color = 'g', label = r'$-%s_{l}^-$' % LF.source, marker = '>'),
    'tumor_vess_src_minus_per_vol'  :  dict(color = 'b', label = r'$-%s_{lv}^-$' % LF.source, marker = 'd'),
    'ext_src_plus_per_vol'   :  dict(color = 'c', label = r'normal, in', marker = 's'),
    'ext_src_minus_per_vol'  :  dict(color = 'b', label = r'normal, out', marker = 'd'),
    'lymph_src_minus_per_vol': dict(color = 'r', label = r'$-%s_{ll}^-$' % LF.source, marker = 's'),
    'vess_src_minus_per_vol' : dict(color = 'g', label = r'$-%s_{lv}^-$' % LF.source, marker = 'o'),
    'src_plus_per_vol' : dict(color = 'b', label = r'$-%s_{lv}^+$' % LF.source, marker = '<'),
  }
  for d in props.itervalues():
    d.update(value_scale = 1.e3)
  props['vess_src_minus_per_vol']['value_scale'] = 1.e4
  props['tumor_vess_src_minus_per_vol']['value_scale'] = 1.e4

  #############################################
  #############################################

  fig, axes = pyplot.subplots(1, 4, figsize = (mastersize[0], mastersize[0]*0.33))
  axes = axes.ravel()
  for ax in axes:
    ax.set(xscale = 'log') #, yscale = 'log')
    axes = axes.ravel()
  mpl_utils.subplots_adjust_abs(fig, left = mpl_utils.mm_to_inch*18,
                                      right = -mpl_utils.mm_to_inch*5,
                                      top  = -mpl_utils.mm_to_inch*5,
                                      bottom = mpl_utils.mm_to_inch*13,
                                      wspace = mpl_utils.mm_to_inch*50,)

  def plot_global_measure(ax, title, xlabel, runs, datasets):
    x = zip(*runs)[2]
    y = zip(*runs)[0]
    #limits = dict(ylim = (-1, 2.))
    limits = dict()
    ax.set(xlabel = xlabel, **limits) #title = title
    #ax.set(yscale = 'log')
    for name in datasets:
        plot_global_single_curve(ax, name, x, y, **props[name])

  #############################################

  datasets = [ 'tumor_src_plus_per_vol', 'tumor_vess_src_minus_per_vol', 'tumor_src_minus_per_vol', 'tumor_out_per_vol' ]
  runs = [ run_labels[x] for x in [ 'vA01', 'vA05', 'vA02', 'default', 'vA03', 'vA09', 'vA04','vA07', 'vA08', 'vA10', 'vA11', 'vA12' ] ]

  plot_global_measure(axes[0], 'leakyness', ur'relative $\lambda_{l,T}$', runs, datasets)

  #############################################

  runs = [ run_labels[x] for x in [ 'vB07', 'vB06', 'vB01', 'vB05', 'vB02', 'default', 'vB03', 'vB04', 'vB08', 'vB09' ] ]
  plot_global_measure(axes[2], 'lymph permeability', ur'relative $L^{(L)}_l$', runs, datasets)

  #############################################

  runs = [ run_labels[x] for x in [ 'vC01', 'vC02', 'vC03', 'vC04', 'vC05', 'vC06'] ]
  plot_global_measure(axes[3], 'tumor lymphatics', ur'$S^{(L)}_T / S^{(L)}_N$', runs, datasets)

  #############################################

  runs = [ run_labels[x] for x in [  'default', 'vE06', 'vE05', 'vE02', 'vE04', 'vE03', 'vE01'] ]
  plot_global_measure(axes[1], 'tumor conductivity', ur'relative $K_l$', runs, datasets)

  def text(ax, txt):
    ax.text(0.05, 0.9, txt, ha = "left", transform = ax.transAxes)

  axes[0].legend(loc = mpl_utils.loc.upper_right)
  axes[0].set(ylabel = ur'[s$^{-1}$]')
  for i, ax in enumerate(axes):
    mpl_utils.add_crosshair(ax, (1,0), color = gridcolor)
    ax.grid(linestyle=':', linewidth=0.5, color= gridcolor)
    text(ax, fig_numbering[i])
  pdfpages.savefig(fig)


def plotIfpVsIffCorrelationData(data, filenames, pdfpages):
  iff_means, ifp_means = data['iff_means'], data['ifp_means']
  
  class VariantFromFilename(object):
    regex1 = re.compile(r'_variant(\w)(\d\d)_')
    regex2 = re.compile(r'_(default)_')
    
    @staticmethod
    def find(filename):
      filename = basename(filename)
      match = VariantFromFilename.regex1.search(filename)
      if match:
        # returns a tuple: (primary variant, i.e. what variable was changed | some text that distiguishes different parameters within the variant)
        return match.group(1), match.group(2)
      match = VariantFromFilename.regex2.search(filename)
      if match:
        return match.group(1), ''
      return 'unkown', ''
  
  variants_to_indices = defaultdict(list)
  for i, fn in enumerate(filenames):
    variant, variant2 = VariantFromFilename.find(fn)
    variants_to_indices[(variant,variant2)].append(i)
  for k, v in variants_to_indices.items():
    variants_to_indices[k] = np.asarray(v)
  
  color_map = dict(default = 'k')
  color_map.update(zip('ABCDE', 'rgbcmy'))
  # in C i apparently varied the lymphatic permeability as well
  label_map = dict(zip('ABCDE', ['(A) $\lambda_{l,T}$', '(C) $L_l^{(L)}$', '(D) $S_T^{(L)} / S_N^{(L)}$', 'just tissue permeability', '(C) $K_l$ & $L_l^{(L)}$']))  
  label_map.update(default = 'Base Case')
  
  markers = '<>do1234sphH^vD*x'
  #markers = r'$\triangleleft$ $\triangleright$ $\diamond$ $\boxtimes$ $\circledcirc$ $\otimes$ $\dagger$ d o s p h H ^ v D'.split(' ')
  #print markers

  marker_map2 = {
    'default' : [
      ('', r'', 1, markers[15]),
    ],
    'A': [
      ('01', r'$\times 10$', 10, markers[0]),
      ('05', r'$\times 4$',  4, markers[1]),
      ('02', r'$\times 2$',  2, markers[2]),
      ('03', r'$\times \frac{1}{2}$', 0.5, markers[3]),
      ('09', r'$\times \frac{1}{4}$', 0.25, markers[7]),
      ('04', r'$\times \frac{1}{10}$', 0.1, markers[4]),
      ('07', r'$\times \frac{1}{20}$', 0.05, markers[5]),
      ('08', r'$\times \frac{1}{100}$', 0.01, markers[6]),
      ('10', r'$\times \frac{1}{200}$', 1./200., markers[8]),
      ('11', r'$\times \frac{1}{500}$', 1./500., markers[9]),
      ('12', r'$\times \frac{1}{1000}$', 1./1000., markers[10]),
    ],
    'B': [
    ('07', r'$\times 100$', 100, markers[11]),
    ('06', r'$\times 50$', 50, markers[12]),
    ('01', r'$\times 10$', 10, markers[0]),
    ('05', r'$\times 5$', 5, markers[1]),
    ('02', r'$\times 2$', 2, markers[2]),
    ('03', r'$\times \frac{1}{2}$', .5, markers[3]),
    ('04', r'$\times \frac{1}{10}$', .1, markers[4]),
    ('08', r'$\times \frac{1}{20}$', .05, markers[5]),
    ('09', r'$\times \frac{1}{100}$', .01, markers[6]),
    ],
    'C': [
    ('01', r'$\times \frac{1}{100}$', 1./100., markers[6]),
    ('02', r'$\times \frac{1}{50}$', 1./50., markers[13]),
    ('03', r'$\times \frac{1}{20}$', 1./20., markers[5]),
    ('04', r'$\times \frac{1}{10}$', 1./10., markers[4]),
    ('05', r'$\times \frac{1}{5}$', 1./5., markers[7]),
    ('06', r'$\times \frac{1}{2}$', 1./2., markers[3]),
    ('07', r'$\times \frac{1}{2}$', 1., markers[16]),
    ],
    'D': [
    ('01', r'$\times 100$', 100., markers[11]),
    ('02', r'$\times 50$', 50., markers[12]),
    ('03', r'$\times 20$', 20., markers[14]),
    ('04', r'$\times 10$', 10., markers[0]),
    ('05', r'$\times 5$', 5., markers[1]),
    ('06', r'$\times 2$', 2., markers[2]),
    ],
    'E': [
    ('01', r'$\times 100$', 100., markers[11]),
    ('02', r'$\times 10$', 10., markers[0]),
    ('03', r'$\times 50$', 50., markers[12]),
    ('04', r'$\times 20$', 20., markers[14]),
    ('05', r'$\times 5$', 5., markers[1]),
    ('06', r'$\times 2$', 2., markers[2]),
    ]
  }
  marker_map = defaultdict(list)
  for k, l in marker_map2.items():
    for v, a, b, c in l:
      marker_map[c].append(round(b, 5))
    l = dict([(v, (a,b,c)) for (v,a,b,c) in l])
    marker_map2[k] = l

  for k, v in marker_map.items():
    v = set(v)
    marker_map[k] = v
    #print('%s -> %s' % (k, v))

  has_legend = {}
  
  W, H = mpl_utils.a4size[0]/2, mpl_utils.a4size[0]/2
  order = ['E','D','A','B','C','default']
  order = dict(zip(order, xrange(len(order))))

  variants_to_indices = sorted(variants_to_indices.items(), key = lambda ((v1,v2),lst): order[v1])
  variants_to_indices = filter(lambda ((k1,k2),v): k1 not in 'D', variants_to_indices)
  # remove cases D and E
  variants_to_indices2 = filter(lambda ((k1,k2),v): k1 not in 'DE', variants_to_indices)
  
  fig, ax = pyplot.subplots(1,1, figsize = (W/2, H/2))
  for (v1,v2), indices in variants_to_indices:
    x = iff_means[indices]
    y = ifp_means[indices]
    label = None if has_legend.get(v1, False) else label_map[v1]; has_legend[v1] = True # label only once per variation set
    color = color_map[v1]
    ax.scatter(x, y, s=10, label = label, **mpl_utils.styleMarkerSP(marker_map2[v1][v2][2], color = color, linewidth = 0.5)) #, edgecolors='none', **scatter_kwargs)
  #ax.set(xlabel = r'IFF [$\mu m /s$]', ylabel = r'IFP [kPa]')
  ax.legend()
  ax.grid(linestyle=':', linewidth=0.5, color=gridcolor)
  #loc = matplotlib.ticker.MaxNLocator(nbins=4)
  #loc = matplotlib.ticker.LinearLocator(6)
  loc = matplotlib.ticker.FixedLocator([0.,1.,2.,3.,4.])
  ax.xaxis.set_major_locator(loc)
  pdfpages.savefig(fig)

  # and plot ... copy pasted ..
  fig, ax = pyplot.subplots(1,1, figsize = (W,H))
  for (v1,v2), indices in variants_to_indices2:
    x = iff_means[indices]
    y = ifp_means[indices]
    label = None if has_legend.get(v1, False) else label_map[v1]; has_legend[v1] = True # label only once per variation set
    color = color_map[v1]
    ax.scatter(x, y, s=10, label = label, **mpl_utils.styleMarkerSP(marker_map2[v1][v2][2], color = color, linewidth = 0.5)) #, edgecolors='none', **scatter_kwargs)
  ax.set(xlabel = r'IFF [$\mu m /s$]', ylabel = r'IFP [kPa]')
  ax.grid(linestyle=':', linewidth=0.5, color=gridcolor)
  #ax.legend()
  pdfpages.savefig(fig)
  
  # plot the same stuff again but now showing the average and standard deviations

  def errorbars(ax, v, x, y, **args):
    avgx = np.average(x)
    #minx = -(np.amin(x)-avgx)
    #maxx = np.amax(x)-avgx
    minx = maxx = np.std(x)
    avgy = np.average(y)
    #miny = -(np.amin(y)-avgy)
    #maxy = np.amax(y)-avgy
    miny = maxy = np.std(y)
    ax.errorbar([avgx], [avgy], xerr = [[minx], [maxx]], yerr = [[miny],[maxy]], zorder = order[v], color = args['color'], marker = '') #marker = args['marker'])
    ax.scatter([avgx], [avgy], **args)
  

  fig, ax = pyplot.subplots(1,1, figsize = (W/2, H/2))
  for (v1,v2), indices in variants_to_indices:
    x = iff_means[indices]
    y = ifp_means[indices]
    label = None if has_legend.get(v1, False) else label_map[v1]; has_legend[v1] = True # label only once per variation set
    color = color_map[v1]
    errorbars(ax, v1, x, y, s=20, label = label, **mpl_utils.styleMarkerSP(marker_map2[v1][v2][2], color = color, linewidth = 0.5)) #, edgecolors='none', **scatter_kwargs)
  ax.legend()
  ax.grid(linestyle=':', linewidth=0.5, color=gridcolor)
  loc = matplotlib.ticker.FixedLocator([0.,1.,2.,3.,4.])
  ax.xaxis.set_major_locator(loc)
  pdfpages.savefig(fig)

  fig, ax = pyplot.subplots(1,1, figsize = (W,H))
  for (v1,v2), indices in variants_to_indices2:
    x = iff_means[indices]
    y = ifp_means[indices]
    label = None if has_legend.get(v1, False) else label_map[v1]; has_legend[v1] = True # label only once per variation set
    color = color_map[v1]
    errorbars(ax, v1, x, y, s=20, label = label, **mpl_utils.styleMarkerSP(marker_map2[v1][v2][2], color = color, linewidth = 0.5)) #, edgecolors='none', **scatter_kwargs)
  ax.set(xlabel = r'IFF [$\mu m /s$]', ylabel = r'IFP [kPa]')
  ax.grid(linestyle=':', linewidth=0.5, color=gridcolor)
  pdfpages.savefig(fig)

  
  fig, ax = pyplot.subplots(1,1, figsize = (W,H))
  for i, (m, vals) in enumerate(marker_map.items()):
    x = np.log10(list(vals))
    y = np.atleast_1d(i*np.ones_like(x)) 
    ax.scatter(x, y, label = (','.join(('%0.4f' % x) for x in vals)), **mpl_utils.styleMarkerSP(m, color = 'k', linewidth = 0.5))
  ax.legend()
  pdfpages.savefig(fig)

##########################################################################
##########################################################################

def GetDataManDataInstances():
  return [DataTissue(), DataGlobalIff(), DataRadialIff(), DataDistanceFromCenter(), DataBasicVessel()]

def do_plotting(filenames):
  rc = matplotlib.rc
  rc('figure', **{'subplot.left' : 0.15,
                  'subplot.right' : 1.-0.15,
                  'subplot.bottom' : 0.2,
                  'subplot.top' : 1.-0.1,
                  'subplot.wspace' : 0.2,
                  'subplot.hspace' : 0.2})
  print 'plotting:', str(filenames)
  fnmeasure = commonOutputName(filenames)
  dataman = myutils.DataManager(100, GetDataManDataInstances())  

  with mpl_utils.PageWriter(fnmeasure+'global.pdf', fileformats=['svg', 'pdf']) as pdfpages:  
    if 0:
      global_ifp_vs_iff_correlation_filename = join(dirname(fnmeasure),'measuredglobal.h5')
      print 'store Ifp vs Iff Correlation in %s' % global_ifp_vs_iff_correlation_filename
      data = ComputeIfpVsIffCorrelationData(dataman, global_ifp_vs_iff_correlation_filename, filenames)
      plotIfpVsIffCorrelationData(data, filenames, pdfpages)
      
    if 1:  
      files = [ h5py.File(fn, 'r+') for fn in filenames ]
      plot_global_header(files, dataman, pdfpages)
      plot_radial(files, dataman, pdfpages)
      plot_single_images2(files[0], dataman, pdfpages)


if __name__ == "__main__":
  import argparse
  parser = argparse.ArgumentParser(description='Analyze IFF distributions.')  
  #parser.add_argument('Iffparams', help = 'choose the parameter for the simulation')  
  parser.add_argument('iffFileNames', nargs='+', type=argparse.FileType('r'), default=sys.stdin, help='iff files to calculate')   
  #parser.add_argument('grp_pattern',help='Where to find the tumor. Usually this is somthing with out*')      
  #this enables access to the default values  
  #atest = parser.add_argument('-a', '--analyze', help = 'loop through all files analyze data and make plot', default=False, action='store_true')  
  #parser.add_argument('-m', '--memory', help= 'Memory assigned by the queing system', type=str, default = '2GB')
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
  
  do_plotting(filenames)

else:
  def plotIffForPaper():
    rc = matplotlib.rc
    rc('figure', **{'subplot.left' : 0.15,
                    'subplot.right' : 1.-0.15,
                    'subplot.bottom' : 0.2,
                    'subplot.top' : 1.-0.1,
                    'subplot.wspace' : 0.2,
                    'subplot.hspace' : 0.2})

    path = sys.argv[1]
    with mpl_utils.PageWriter('forpaper', fileformats=['svg', 'pdf']) as pdfpages:
      plotComparisonsForIffPaper(pdfpages, path)

#    g = lambda fn: glob.glob(os.path.join(path, fn))
#    do_plotting(g('iff_default_*.h5'))
#    do_plotting(g('iff_variantA08_*.h5')) # vess permeability x 0.01
#    do_plotting(g('iff_variantA12_*.h5')) # vess permeability x 0.001
#    do_plotting(g('iff_variantA01_*.h5')) # vess permeability x 10
#    do_plotting(g('iff_variantB01_*.h5')) # lymph permeability x 10
#    do_plotting(g('iff_variantB04_*.h5')) # lymph permeability x 0.1
    #do_plotting(g('iff_variantC01_*.h5')) # tumor lymphatic density = 1/100 * tissue lymphatic density
    #do_plotting(g('iff_variantC06_*.h5')) # tumor lymphatic density = 0.5 * tissue lymphatic density
#    do_plotting(g('iff_variantC04_*.h5')) # tumor lymphatic density =
#    do_plotting(g('iff_variantC07_*.h5')) # tumor lymphatic density =
#    do_plotting(g('iff_variantE02_*.h5')) # tissue & lymphatic permeability x10
#    do_plotting(g('iff_variantE01_*.h5')) # tissue & lymphatic permeability x100

