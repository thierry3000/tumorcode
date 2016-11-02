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
#!/usr/bin/env python
# -*- coding: utf-8 -*-
if __name__ == '__main__':
  import os.path, sys
  sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'..'))
  
import os,sys
from os.path import basename, splitext, commonprefix
import h5py
import h5files
import os,sys
import numpy as np
import collections
from mystruct import Struct
import extensions # for hdf5 support in np.asarray
import myutils
import krebsutils
import posixpath
import math
import pprint
from myutils import f2s, f2l
import random


from plotBulkTissue import commonOutputName, colorbar, contour, imslice, imshow, with_contour_factory, with_cb
from analyzeGeneral import DataBasicVessel, DataVesselSamples, DataVesselRadial, DataTumorTissueSingle, DataDistanceFromCenter, BinsSpecRange
from plotVessels import PlotRadiusHistogram
from quantities import Prettyfier

import matplotlib
import matplotlib.pyplot as pyplot
import mpl_utils

import analyzeGeneral


class DataTumor3dRendering(object):
  keywords = ['3drendering']

  def obtain_data(self, dataman, dataname, *args):
    if dataname == '3drendering':
      import cStringIO
      import povrayRenderTumor
      import PIL.Image as Image

      f, group, showVessels, showTumor = args[0], args[1], args[2], args[3]

      ld = dataman.obtain_data('ld', f)
      is_cubic = ld.shape[2]*2 > ld.shape[0]

      def read(gmeasure, groupname):
        ds = np.asarray(gmeasure[groupname])
        memfile = cStringIO.StringIO(ds)
        img = Image.open(memfile)
        arr = np.array(img)
        return arr

      def write(gmeasure, groupname):
        tempfn = '%s-%x.png' % (splitext(f.filename)[0], random.getrandbits(32))
        vesselgroup = f[group]['vessels'] if showVessels else None
        tumorgroup = analyzeGeneral.tumor_group(f, group) if showTumor else None
        povrayRenderTumor.renderScene(vesselgroup, tumorgroup,
              tempfn,
              res=(800, 800),
              aa=1,
              cam= 'corner' if is_cubic else 'topdown',
              colored_slice=True,
              out_alpha=False,
              num_threads=3,
              temp_file_dir = '/tmp',
              camera_distance = 2. if is_cubic else 1.45
              )
        #fn = 'prez3d-stf-test-st0-MOtcx10ncx1-73877c0b.png'
        stringfile = open(tempfn).read()
        bytestream = np.fromstring(stringfile, dtype=np.byte)
        gmeasure.create_dataset(groupname, (len(bytestream),), data = bytestream)

      dataname = '3drendering' + ('_w_vessels' if showVessels else '') + ('' if showTumor else '_notum')
      fm = myutils.MeasurementFile(f, h5files)
      ret = myutils.hdf_data_caching(read, write, fm, (group,dataname), (0, 11))
      return ret


class DataTumorMeasureCurvature(object):
  keywords = ['field_curvature', 'field_curvatureradius', 'field_interfacedelta', 'shape_metrics', 'curvature_samples']

  def obtain_data(self, dataman, dataname, *args):
    f, args = args[0], args[1:]
    obtain_data = lambda *args: dataman.obtain_data(args[0], f, *args[1:])
    ld = obtain_data('ld')
    if dataname == 'field_interfacedelta':
      group, = args
      distmap = obtain_data('fieldvariable', 'dist_tumor', group)
      return krebsutils.smooth_delta(distmap, 0.5*ld.scale)
    ###
    if dataname == 'field_curvature':
      group, = args
      distmap = obtain_data('fieldvariable', 'dist_tumor', group)
      curvature = krebsutils.curvature(ld, distmap, False, False)
      return curvature
    ###
    if dataname == 'field_curvatureradius':
      group, = args
      distmap = obtain_data('fieldvariable', 'dist_tumor', group)
      curvature = krebsutils.curvature(ld, distmap, False, False)
      curvature = np.ma.array(curvature, mask = curvature<>0.)
      # kappa = (1/r0 + 1/r1), r0=r1 -> r = 2/kappa
      curvature = 2.* np.ma.power(curvature, -1.)
      return curvature
    ###
    if dataname == 'shape_metrics':
      group, = args
      def read(gmeasure, groupname):
        return dict((k, float(v[...])) for (k,v) in gmeasure[groupname].iteritems())
      def write(gmeasure, groupname):
        dim = 2 if ld.shape[2]==1 else 3
        cellvol = ld.scale**dim
        if 0:
          vtkds, = extractVtkFields.Extractor(statedata['tumor'], ['ls']).asVtkDataSets()
          area = integrate_surface(vtkds)
          del vtkds
        theta_tumor = obtain_data('fieldvariable','theta_tumor', group)
        vol = np.sum(theta_tumor)*cellvol
        delta = obtain_data('field_interfacedelta', group)
        # grid cells contributing to the interfacial area
        interface_indices = np.nonzero(delta)
        interface_weights = delta[interface_indices]
        # area estimate
        area = np.sum(delta)*cellvol
        # radius estimate
        radialmap = krebsutils.make_radial_field(ld)
        radius = np.average(radialmap[interface_indices], weights = interface_weights)
        # comparison with the area of a sphere/cylinder with the same volume as the shape
        if dim == 3:
          sphere_equiv_radius = math.pow(vol*3./(4.*math.pi), 1./3.)
          sphere_equiv_area = 4.*math.pi*(sphere_equiv_radius ** 2)
          cylinder_equiv_radius = math.sqrt(vol / ld.GetWorldSize()[2]  / math.pi)
          cylinder_equiv_area = 2.*math.pi*cylinder_equiv_radius*ld.GetWorldSize()[2]
        else:
          # pi r^2 = A -> r = (A/pi)^(1/2)
          sphere_equiv_radius = math.pow(vol/math.pi, 1./2.)
          sphere_equiv_area   = 2.*math.pi * sphere_equiv_radius
          cylinder_equiv_radius = sphere_equiv_radius
          cylinder_equiv_area   = sphere_equiv_area
        # how much of a perfect sphere/cylinder the shape is
        sphericity = sphere_equiv_area / area
        cylindericity = cylinder_equiv_area / area

        res = dict(area = area, volume = vol,
                 radius=radius, sphericity = sphericity, cylindericity = cylindericity,
                 sphere_equiv_radius = sphere_equiv_radius, sphere_equiv_area = sphere_equiv_area,
                 cylinder_equiv_radius = cylinder_equiv_radius, cylinder_equiv_area = cylinder_equiv_area
                 )
        g = gmeasure.create_group(groupname)
        for k,v in res.iteritems():
          g.create_dataset(k, data = v)
        gmeasure.file.flush()
      fm = myutils.MeasurementFile(f, h5files)
      ret = myutils.hdf_data_caching(read, write, fm, (group,'shape_metrics'), (0,1))
      return ret
    ###
    if dataname == 'curvature_samples':
      group, = args
      delta = obtain_data('field_interfacedelta', group)
      interface_indices = np.nonzero(delta)
      curv = obtain_data('field_curvature', group)
      return curv[interface_indices], delta[interface_indices]


#class DataTissueRadial(object):
#  keywords = [ 'tissue_radial']
#
#  def obtain_data(self, dataman, dataname, *args):
#    f, group = args[0], args[1]
#    fileid = 'file'+myutils.checksum(basename(f.filename))
#    obtain_data = lambda *args: dataman.obtain_data(args[0], f, *args[1:])
#    ld = obtain_data('ld')
#
#    if dataname == 'tissue_radial':
#      param_sample_length = 50.
#      cellvol = ld.scale**3
#      bins_rad = np.arange(0, 20000, 30.)
#      bins_dist = np.arange(-10000., 10000., 30.)
#
#      def write(gmeasure, groupname):
#        import plotVessels
#        grp_radial = gmeasure.create_group(groupname)
#        #grp_radial = gmeasure.create_group('radial')
#        #grp_samples = gmeasure.create_group('samples')
#
#        distmap = obtain_data('fieldvariable', 'dist_tumor', group)
#        radialmap = krebsutils.make_radial_field(ld)
#
#        far = 1.e10
#        # load vessels
#        vessels = krebsutils.read_vesselgraph(f[group]['vessels'], ['position',  'flags', 'shearforce', 'radius', 'flow', 'maturation'])
#        vessels = vessels.get_filtered(edge_indices = np.bitwise_and(vessels.edges['flags'],krebsutils.CIRCULATED) != 0)
#        # generate sample positions
#        s = plotVessels.generate_samples(vessels, 'position', 'nodes', param_sample_length)
#        # obtain distance field value for each position sample
#        dist = krebsutils.sample_field(s, distmap, ld, linear_interpolation=True, extrapolation_value = far)
#        rad = krebsutils.sample_field(s, radialmap, ld, linear_interpolation=True, extrapolation_value = far)
#        del s
#
#        grp_radial.create_group('vs_r').create_dataset('bins', data = np.average((bins_rad[:-1],bins_rad[1:]), axis=0))
#        grp_radial.create_group('vs_dr').create_dataset('bins', data = np.average((bins_dist[:-1],bins_dist[1:]), axis=0))
#
#        for name in ['shearforce', 'radius', 'flow', 'maturation' ]:
#          s = plotVessels.generate_samples(vessels, name, 'edges', param_sample_length)
#          myutils.MeanValueArray.fromHistogram1d(bins_rad, rad, s).write(grp_radial['vs_r'], name)
#          myutils.MeanValueArray.fromHistogram1d(bins_dist, dist, s).write(grp_radial['vs_dr'], name)
#
#        def make_mvd(flavor, s, bins, smap, mask_bound):
#          a = myutils.MeanValueArray.fromHistogram1d(bins, s, np.ones_like(s))
#          b = myutils.MeanValueArray.fromHistogram1d(bins, smap.ravel(), np.ones_like(smap.ravel()))
#          a.cnt = b.cnt.copy()
#          a.sum *= param_sample_length/cellvol
#          a.sqr *= a.sum**2
#          a.write(grp_radial[flavor], 'mvd')
#
#        make_mvd('vs_r', rad, bins_rad, radialmap, 6000.)
#        make_mvd('vs_dr', dist, bins_dist, distmap, 1000.)
#    #    pyplot.errorbar(dstdata['vs_dr']['bins'], dstdata['vs_dr']['mvd_avg'], yerr=dstdata['vs_dr']['mvd_std'])
#    #    pyplot.show()
#
#        del rad, dist
#
#        for name in ['phi_cells', 'sources', 'phi_tumor', 'press', 'oxy', 'gf', 'phi_vessels']:
#          a = obtain_data('fieldvariable', name, group)
#          myutils.MeanValueArray.fromHistogram1d(bins_rad, np.ravel(radialmap), np.ravel(a)).write(grp_radial['vs_r'], name)
#          myutils.MeanValueArray.fromHistogram1d(bins_dist, np.ravel(distmap), np.ravel(a)).write(grp_radial['vs_dr'], name)
#
#        def make_velocity_projection(flavor, distancemap, bins, v):
#          vnorm = np.zeros(distancemap.shape, dtype=np.float32)
#          dgrad = krebsutils.field_gradient(distancemap, spacing=ld.scale)
#          for vv, g in zip(v, dgrad):
#            vnorm += vv*g
#          del dgrad
#          myutils.MeanValueArray.fromHistogram1d(bins, np.ravel(distancemap), np.ravel(vnorm)).write(grp_radial[flavor], 'vel')
#
#        # velocities, projected radially outward
#        vel_field = obtain_data('fieldvariable', 'vel', group)
#        vel_field = tuple(vel_field[:,:,:,i] for i in range(3))
#        make_velocity_projection('vs_r', radialmap, bins_rad, vel_field)
#        make_velocity_projection('vs_dr', distmap, bins_dist, vel_field)
#        #del vel_field
#
#      def read(gmeasure, groupname):
#        grp = gmeasure[groupname]
#        res = {}
#        for flavor in ['vs_r', 'vs_dr']:
#          for key in grp[flavor].iterkeys():
#            if key == 'bins': continue
#            res[flavor,key] = myutils.MeanValueArray.read(grp[flavor], key)
#          bins = np.asarray(grp[flavor]['bins'])
#          bins = np.average((bins[1:], bins[:-1]), axis=0)
#          res[flavor,'bins'] = myutils.MeanValueArray.fromArray(bins)
#        return res
#
#      fm = myutils.MeasurementFile(f)
#      return myutils.hdf_data_caching(read, write, fm, (group, 'radial'), (0,2))


def PlotRadial(pdfwriter, dataman, resultfiles, distance_distribution_name):
  binspec = BinsSpecRange(-10000., 10000., 100.)
  curves_by_path_and_name = myutils.MultiDict(list)
  for f in files:
    fm = myutils.MeasurementFile(f.f, h5files)
    for outgroup_name in f.groupnames:
      gvessels = f.f[outgroup_name]['vessels']
      gtumor   = f.f[outgroup_name]['tumor']
      for name in [ 'mvd', 'phi_vessels', 'shearforce', 'velocity', 'radius', 'maturation', 'flow' ]:
        d = dataman.obtain_data('basic_vessel_radial', name, gvessels, gtumor, 30., binspec, distance_distribution_name, None, (fm,outgroup_name))
        curves_by_path_and_name[outgroup_name, name].append(d)
      #tumor_radius = f.f[outgroup_name]['tumor'].attrs['TUMOR_RADIUS'] # this needs work ...
      tumor_radius = dataman.obtain_data('approximate_tumor_radius', gtumor)
      curves_by_path_and_name[outgroup_name, 'approximate_tumor_radius'].append(myutils.MeanValueArray(1,tumor_radius,tumor_radius**2))
  for k, v in curves_by_path_and_name.items():
    curves_by_path_and_name[k] = myutils.MeanValueArray.fromSummation(v_.avg for v_ in v)
  # ---- here goes the plotting -----#
  charsize = 12/90.
  figscale = 1.8
  fig, axes = pyplot.subplots(3, 2, dpi = 90, figsize=(7.*figscale,5.*figscale))
  mpl_utils.subplots_adjust_abs(fig, left=6*charsize, right=-2*charsize, top=-3.*charsize, bottom=8.*charsize, hspace=8.*charsize, wspace=12*charsize)
  axes = axes.ravel()

  default_colors = 'rgbcmykrgbcmykrgbcmyk'
  pathorder = resultfiles[0].groupnames
  labels_by_path = dict((g,'t=%s' % f2s(resultfiles[0].get_group_time(g))) for g in resultfiles[0].groupnames )
  bins = binspec.arange()
  bins = 1.e-3*np.average((bins[1:], bins[:-1]), axis=0)
  ld = dataman.obtain_data('ld', resultfiles[0].f)
  rmax = ld.GetWorldSize()[0]*0.25*1.e-3

  def plot(ax, name, scalefactor = 1., label = None, colors = default_colors, errorbars = True, zero_ylim = True):
    for i, path in enumerate(pathorder):
      curve = curves_by_path_and_name[path,name]
      label = labels_by_path[path] if label else None
      mask = ~curve.avg.mask
      if errorbars:
        mpl_utils.errorbar(ax, bins[mask], scalefactor*curve.avg[mask], yerr=scalefactor*curve.std_mean[mask], label = label,
                           marker = None, color = colors[i], every = 2)
      else:
        ax.plot(bins[mask], scalefactor*curve.avg[mask], label = label, color=colors[i])
    if zero_ylim:
      _, ymax = ax.get_ylim()
      ax.set_ylim((0., ymax))
    if distance_distribution_name=='levelset':
      ax.set(xlim=(-2.,2.))
      mpl_utils.add_crosshair(ax, (0.,None), ls=':')
    else:
      ax.set(xlim=(0.,rmax))
      for i, path in enumerate(pathorder):
        mpl_utils.add_crosshair(ax, (curves_by_path_and_name[path,'approximate_tumor_radius'].avg*1.e-3,None), ls=':')

  def fmt_axis_labels(ax, name, scalefactor = None, hidex = True):
    ax.set(ylabel = (r'$%s$ %s' % (Prettyfier.get_sym(name), Prettyfier.get_bunit(name))) + ((r' $\times\, %s$' % f2l(scalefactor, exponential=True)) if scalefactor else ''),
           title  = Prettyfier.get_label(name))
    if hidex: ax.set(xticklabels = [])
    else:
      xsym = r'\phi' if distance_distribution_name=='levelset' else r'|x|'
      ax.set(xlabel = r'$%s$ [$mm$]' % xsym)

  plot(axes[0], 'mvd', label = True, scalefactor = 1.e6)
  fmt_axis_labels(axes[0], 'mvd')
  axes[0].legend()

  plot(axes[1], 'phi_vessels', scalefactor = 1.)
  fmt_axis_labels(axes[1], 'phi_vessels', scalefactor = 1)

  plot(axes[2], 'shearforce')
  fmt_axis_labels(axes[2], 'shearforce',)

  plot(axes[3], 'velocity')
  fmt_axis_labels(axes[3], 'velocity',)

  plot(axes[4], 'radius')
  fmt_axis_labels(axes[4], 'radius', hidex=False)

  plot(axes[5], 'flow', scalefactor = 1)
  fmt_axis_labels(axes[5], 'flow', scalefactor = 1, hidex=False)

  #plot(axes[5], 'maturation', scalefactor = 1)
  #fmt_axis_labels(axes[5], 'maturation', scalefactor = 1, hidex=False)

  pdfwriter.savefig(fig, postfix='_radial')
  # ---- done plotting -----#



class ResultFile(object):
  def __init__(self, f, dataman, pattern=None):
    self.f = f
    self.dataman = dataman
    self.obtain_data = lambda *args: self.dataman.obtain_data(args[0], f, *args[1:])
    self.get_group_time = lambda q: self.obtain_data('time', q)
    groups = myutils.getTimeSortedGroups(f['/'], 'out', key='time')
    groups = [ posixpath.basename(group.name) for group in groups ]
    if pattern:
      s = set(myutils.walkh5(f, pattern))
      groups = [ g for g in groups if g in s ]
    if len(groups)==0:
      raise ValueError('There is no tumor group found in this file!')
    self.groupnames = groups


def obtain_volumes(f, group, dataman):
  ld = dataman('ld',f)
  dist_viabletumor = dataman('fieldvariable',f, 'dist_viabletumor', group)
  dist_tumor = dataman('fieldvariable',f, 'dist_tumor', group)
  pyplot.imshow(imslice(dist_viabletumor))
  pyplot.colorbar()
  pyplot.show()
  nsites_total = np.prod(ld.shape);
  nsites_tumor = np.sum(np.asarray(dist_tumor < 0., dtype=np.uint))
  nsites_viabletumor = np.sum(np.asarray(dist_viabletumor < 0., dtype=np.uint))
  return np.prod(ld.GetWorldSize()), float(nsites_tumor)/nsites_total, float(nsites_viabletumor)/nsites_total


def printVolumina(filenames, group):
  data = []
  dataman = myutils.DataManager(100, [DataTumorTissueSingle()])
  for fn in filenames:
    with h5py.File(fn, 'r') as f:
      d = obtain_volumes(f, group, dataman)
      data.append(d)
  data = np.asarray(data).transpose()
  v_tot = data[0]
  vr_tumor = data[1]
  vr_viable = data[2]
  def printavg(name, data):
    avg, std = np.average(data), np.std(data)
    print '%s = %s +/- %s' % (name, f2s(avg), f2s(std))
  printavg('V_Sys', v_tot)
  printavg('V_trum', vr_tumor*v_tot)
  printavg('V_viable',vr_viable*v_tot)
  printavg('rV_trum', vr_tumor)
  printavg('rV_viable',vr_viable)


def plotSnapshots(resultfile, pdfpages):
  parameters = resultfile.f['parameters/tumor']
  param_oxy_death = float(parameters.attrs['o2_necro_threshold'])
  param_oxy_prol = float(parameters.attrs['o2_prol_threshold'])
  ld = resultfile.obtain_data('ld')
  is3d = ld.shape[2]>1
  has_vessels = 'vessels' in resultfile.f[resultfile.groupnames[-1]]
  for groupname in resultfile.groupnames[:]:
    #dpi = 90
    fig, axes = pyplot.subplots(2,2, sharex=True, sharey=True, figsize=(8,7))
    fig.subplots_adjust(left=0.10, right=0.95, top=0.95, bottom=0.10, wspace=0.2, hspace=0.2)
    #mpl_utils.remove_frame(*axes.ravel())

    dist_viabletumor = resultfile.obtain_data('fieldvariable', 'dist_viabletumor', groupname, 'imslice')
    oxy = resultfile.obtain_data('fieldvariable', 'oxy', groupname, 'imslice')

    with_contour = with_contour_factory(dist_viabletumor, ld, worldscaling=1.e-3)

    @with_cb
    def imshow_cb(ax, data, ld, cmap, crange, vscale):
      return imshow(ax, data, ld, cmap=cmap, vscale=vscale, crange=crange, worldscaling=1.e-3)

    @with_contour
    def plot_phi_cells(ax):
      data = resultfile.obtain_data('fieldvariable', 'phi_viabletumor', groupname, 'imslice')
      ax.set(title = ur'$\phi$')
      imshow_cb(ax, data, ld, matplotlib.cm.coolwarm, (0.4, 0.5), 1.)
      contour(ax, oxy, ld, levels = [param_oxy_death, param_oxy_prol], colors = ('r', 'b'),  worldscaling=1.e-3)

    def plot_oxy(ax):
      ax.set(title = '$O_2$')
      imshow_cb(ax, oxy, ld, matplotlib.cm.Reds, (0.,oxy.max()), 1.)
      contour(ax, oxy, ld, levels = [param_oxy_death, param_oxy_prol], colors = ('r', 'b'),  worldscaling=1.e-3)

    def plot_sources(ax):
      data = resultfile.obtain_data('fieldvariable', 'sources', groupname, 'imslice')
      ax.set(title = ur'$\Gamma_\phi$')
      imshow_cb(ax, data, ld, matplotlib.cm.coolwarm, 'zero-centered', 1.e3)
      contour(ax, oxy, ld, levels = [param_oxy_death, param_oxy_prol], colors = ('r', 'b'),  worldscaling=1.e-3)

    def plot_distmap(ax):
      data = resultfile.obtain_data('fieldvariable', 'dist_tumor', groupname, 'imslice')
      ax.set(title = 'distance map')
      imshow_cb(ax, data, ld, matplotlib.cm.coolwarm, (-5.*ld.scale, 5.*ld.scale), 1.)

    def plot_curv(ax):
      data = imslice(resultfile.obtain_data('field_curvatureradius', groupname))
      delta = imslice(resultfile.obtain_data('field_interfacedelta', groupname))*ld.scale
      ax.set(title = '$\kappa$')
      imshow_cb(ax, data*delta, ld, matplotlib.cm.coolwarm, 'zero-centered', 1.)

    @with_contour
    def plot_vesselvolume(ax):
      data = resultfile.obtain_data('fieldvariable','phi_vessels', groupname, 'imslice')
      ax.set(title = ur'$\phi_v$')
      imshow_cb(ax, data, ld, matplotlib.cm.coolwarm, (0., 1.), 1.)
      contour(ax, oxy, ld, levels = [param_oxy_death], colors = 'r',  worldscaling=1.e-3)

    plot_phi_cells(axes[0,0])
    plot_oxy(axes[0,1])
    plot_sources(axes[1,0])
    #plot_distmap(axes[1,0])
    if has_vessels:
      plot_vesselvolume(axes[1,1])
    else:
      plot_curv(axes[1,1])
    pdfpages.savefig(fig)


def linearFit(x, y):
  x = np.asarray(x)
  y = np.asarray(y)
  from scipy.optimize import leastsq
  func = lambda p, x, y: (x*p[0]+p[1]-y)
  p, success = leastsq(func, (1, 0), args=(x, y))
  return p


def plotShapeMetrics(resultfile, pdfpages):
  ld = resultfile.obtain_data('ld')
  # get data as time series
  data = collections.defaultdict(list)
  for groupname in resultfile.groupnames:
    metrics = resultfile.obtain_data('shape_metrics', groupname)
    for k,v in metrics.iteritems():
      data[k].append(v)
    data['time'].append(resultfile.obtain_data('time',groupname))

  # plot it
  for groupname in resultfile.groupnames[-1:]:
    metrics = resultfile.obtain_data('shape_metrics', groupname)
    curv, curv_weights = resultfile.obtain_data('curvature_samples', groupname)
    #pprint.pprint(metrics)
    fig, axes = pyplot.subplots(2,2, figsize = (8, 6))
    text = []
    text += [ 'time : {}'.format(f2s(resultfile.obtain_data('time', groupname))) ]
    text += [ 'surface tension : {}'.format(f2s(float(resultfile.f['parameters/tumor'].attrs['surface_tension']))) ]
    text += [ 'sphericity : {}'.format(f2s(metrics['sphericity'])) ]
    text += [ 'approx. radius : {}'.format(f2s(metrics['radius']))]
    text = '\n'.join(text)
    axes[0,0].text(0.1,0.1,text, fontsize = 14)
    mpl_utils.remove_frame(axes[0,0])

    h, bins = np.histogram(curv, bins=10, weights = curv_weights, normed = True)
    bins = 0.5*(bins[1:]+bins[:-1])
    axes[0,1].bar(bins, h, width=bins[1]-bins[0])
    axes[0,1].set(title = ur'$\kappa$', xlabel = ur'[1/\u03BCm]', ylabel = 'N')

    axes[1,0].plot(data['time'], data['radius'], marker = 'x')
    axes[1,0].set(xlabel = '$t$', ylabel = r'$r_{tum}$')
    v, y0 = linearFit(data['time'], data['radius'])
    x = np.asarray([0., data['time'][-1]*1.1])
    y = y0 + v*x
    axes[1,0].plot(x, y, label = ur'v = %s [\u03BCm/h]' % f2s(v))
    axes[1,0].legend()

    axes[1,1].set_visible(False)

    fig.subplots_adjust(left=0.10, right=0.95, top=0.95, bottom=0.10, wspace=0.2, hspace=0.2)
    pdfpages.savefig(fig)






def plot3dRendering(resultfile, pdfpages, showVessels = False, showTumor = True):
    arr = resultfile.obtain_data('3drendering', resultfile.groupnames[-1], showVessels, showTumor)
    fig, ax = pyplot.subplots(1, 1, figsize = (8,7))
    fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1, wspace=0., hspace=0.)
    ax.imshow(arr, interpolation = 'none')
    mpl_utils.remove_frame(ax)
    pdfpages.savefig(fig)



def plotSingleRun(fn):
  out = splitext(basename(fn))[0]
  f = h5files.open(fn, 'r')
  rc = matplotlib.rc
  rc('figure', figsize=(7,7), dpi=100)
  rc('font', size = 8.)
  rc('axes', titlesize = 10., labelsize = 10.)
  rc('pdf', compression = 6, fonttype = 42)
  rc('figure', **{'subplot.left' : 0.02,
                  'subplot.right' : 1.-0.05,
                  'subplot.bottom' : 0.01,
                  'subplot.top' : 1.-0.05,
                  'subplot.wspace' : 0.1,
                  'subplot.hspace' : 0.1})
  #rc('savefig', facecolor='none', edgecolor='none', dpi=100)
  rc('savefig', dpi=100)
  rc('font', **{'family':'sans-serif'}) #,'sans-serif':['Helvetica']})
  rc('path', simplify_threshold=0.01)
  #rc('text', usetex=True, **{ 'latex.unicode' : True })

  with mpl_utils.PdfWriter(out+".pdf") as pdfpages:
    dataman = myutils.DataManager(100, [DataTumorTissueSingle(), DataTumorMeasureCurvature(), DataTumor3dRendering(), DataTissueRadial(), DataTissueRadialAveraged()])
    resultfile = ResultFile(f, dataman)
    is3d = resultfile.obtain_data('ld').shape[2]>1
    if is3d:
      plot3dRendering(resultfile, pdfpages, showVessels = True, showTumor = False)
      plot3dRendering(resultfile, pdfpages)
      plot3dRendering(resultfile, pdfpages, showVessels = True)
    plotSnapshots(resultfile, pdfpages)
    #plotShapeMetrics(resultfile, pdfpages)
#    plotRadial(resultfile, pdfpages, 'vs_dr')
#    plotRadial(resultfile, pdfpages, 'vs_r')


def plotSnapshotMatrix(resultfile_matrix, pdfpages, xlabels =  [], ylabels = []):
  nrows, ncols = resultfile_matrix.shape
  fig, axes = pyplot.subplots(nrows, ncols, sharex=True, sharey=True, figsize=(ncols * 3.+ 2., nrows * 3 + 2.))
  for ax, resultfile in zip(axes.ravel(), resultfile_matrix.ravel()):
    ld   = resultfile.obtain_data('ld')
    data = resultfile.obtain_data('fieldvariable', 'dist_tumor', resultfile.groupnames[-1], 'imslice')
    #imshow_(ax, data, ld, matplotlib.cm.coolwarm, None, 1.)
    contour(ax, data, ld, levels = [0.], colors = 'k')
    ax.set(xticklabels = [], yticklabels = [])
  for xlabel, ax in zip(xlabels, axes[-1,:]):
    ax.set(xlabel = xlabel)
  for ylabel, ax in zip(ylabels, axes[:,0]):
    ax.set(ylabel = ylabel)
  pdfpages.savefig(fig)


def fix_time(g):
  if 'TIME' in g.attrs and posixpath.basename(g.name).startswith('out'):
    t = g.attrs['TIME']
    del g.attrs['TIME']
  else:
    t = None
  if not 'time' in g.attrs and t is not None:
    g.attrs['time'] = t
  g.file.flush()
  if isinstance(g, (h5py.Group,h5py.File)):
    for child in g.itervalues():
      fix_time(child)


def plotComparison():
  postfixes1 = [
    '-st0',
    '-st100',
    '-st1000',
    '-st10000',
  ]
  labels1 = [ '$0$', '$10^{0}$', '$10^{2}$', '$10^{3}$', '$10^{4}$']
  postfixes2 = [
    '-MOtcx1ncx1',
    '-MOtcx10ncx1',
  ]
  labels2 = [ '$5k\,/\,5k$', '$50k\,/\,5k$' ]
  filename_matrix = np.asarray([
    [ 'prez2d-stf-test'+p1+p2+'.h5' for p2 in postfixes2 ] for p1 in postfixes1
  ])
  outfn = 'prez2d-stf-test'

  dataman = myutils.DataManager(100, [DataTumorTissueSingle(), DataTumorMeasureCurvature(), DataTumor3dRendering()])

  file_matrix = np.asarray([ ResultFile(h5py.File(fn, 'r+'),dataman) for fn in filename_matrix.ravel()], dtype=np.object).reshape(filename_matrix.shape)

  with mpl_utils.SinglePageWriter(outfn+".pdf") as pdfpages:
    plotSnapshotMatrix(file_matrix, pdfpages, xlabels = labels2, ylabels = labels1)


if __name__ == "__main__":
  filenames = sys.argv[1:-1]
  pattern = sys.argv[-1]
  if 0:
    plotSingleRun(sys.argv[1])
  else:
    dataman = myutils.DataManager(100, [DataTumorTissueSingle(), DataVesselRadial(), DataVesselSamples(), DataBasicVessel(), DataDistanceFromCenter()])
    files = [ ResultFile(h5files.open(fn, 'r+'), dataman, pattern) for fn in filenames ]
    for f in files:
      fix_time(f.f)
    common_filename = myutils.sanitize_posixpath(splitext(commonprefix(filenames))[0])
    with mpl_utils.PageWriter(common_filename+'.pdf', formats = ['pdf']) as pdfwriter:
      PlotRadial(pdfwriter, dataman, files, 'radial')

      #vesselgroups = [ f.f[f.groupnames[-1]]['vessels'] for f in files ]
      #fig, ax = pyplot.subplots(1,1, figsize = (mpl_utils.a4size[0]*0.8, mpl_utils.a4size[0]*0.4))
      #PlotRadiusHistogram(ax, dataman, vesselgroups, krebsutils.WITHIN_TUMOR)
      #pdfwriter.savefig(fig, postfix='_radiushisto')
  if 0:
    printVolumina(sys.argv[2:], sys.argv[1])