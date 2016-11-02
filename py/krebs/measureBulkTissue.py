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
import os,sys
from os.path import join, basename, dirname, splitext
#import time
if __name__ == '__main__':
  sys.path.append(join(dirname(__file__),'..'))
import h5py
import h5files
import numpy as np
import vtkcommon
import extensions # for asarray with h5py support
import krebsutils
import math
from mystruct import Struct
import myutils
import posixpath
from copy import deepcopy
from collections import defaultdict
from pprint import pprint

#from plotBulkTissue import commonOutputName, LabelFactory, colorbar, contour, imslice, imshow, ColorMaps
from plotBulkTissue import commonOutputName

from plotIff import mk_LF_, mk_CM_

from matplotlib.ticker import MaxNLocator
import matplotlib
import matplotlib.pyplot as pyplot
import mpl_utils

LabelFactory=mk_LF_()
LF = LabelFactory
ColorMaps = mk_CM_()
CM = ColorMaps

mastersize = mpl_utils.a4size
a4size = mpl_utils.a4size
gridcolor = (0.7,0.7,0.7)
f2s = lambda x: myutils.f2s(x, latex=True)

import extractVtkFields, plotVessels

bins_rad = np.arange(0, 20000, 30.)
bins_dist = np.arange(-10000., 10000., 30.)



def get_tumld(tumorgroup):
  p = tumorgroup['ptc'].attrs['LATTICE_PATH']
  tum_ld = krebsutils.read_lattice_data_from_hdf(tumorgroup.file[p])
  #field_ld is in root
  #tum_ld = krebsutils.read_lattice_data_from_hdf(tumorgroup.parent.parent['field_ld'])
  #tum_ld = krebsutils.read_lattice_data_from_hdf(tumorgroup.file[tumorgroup['conc'].attrs['LATTICE_PATH']])  
  return tum_ld

def calc_distmap(tumorgroup):
  ld = get_tumld(tumorgroup)
  distmap = np.asarray(tumorgroup['ls']) > 0
  distmap = krebsutils.flood_fill(distmap, (0,0,0))
  distmap = np.logical_not(distmap)
  distmap = krebsutils.distancemap(distmap)*ld.scale
  return distmap


def integrate_surface(vtkds):
  surface = vtkcommon.vtkScaleDataSet(vtkds, 1)
  surface.GetCellData().SetActiveScalars('ls')
  surface = vtkcommon.vtkCellDataToPointData(surface)
  surface = vtkcommon.vtkContour(surface, [0.])
  res_cd, res_pd, res_surf = vtkcommon.vtkIntegrateData(surface)
  return res_surf




def generate_data(statedata, dstgroup):
  """
    main routine that does all the measurement. Enable Individual parts. Computed date overwrites old data.
  """
  print statedata.file.filename, statedata.name

  time = statedata.attrs['time']
  tum_grp = statedata['tumor']
  tum_ld  = get_tumld(tum_grp)
  cellvol = tum_ld.scale**3
  distmap = calc_distmap(tum_grp)
  radialmap = krebsutils.make_radial_field(tum_ld)
  axial_max_rad = np.max(np.abs(tum_ld.worldBox.ravel()))
  ptc = np.asarray(tum_grp['ptc'])

  if 1:
    ## shape data: rim surface area, volume
    vtkds, = extractVtkFields.Extractor(statedata['tumor'], ['ls']).asVtkDataSets()
    area = integrate_surface(vtkds)
    del vtkds
    vol = np.sum(ptc)*cellvol
    radius = np.average(radialmap[np.nonzero(np.logical_and(distmap>-2*tum_ld.scale, distmap<2*tum_ld.Scale()))])
    sphere_equiv_radius = math.pow(vol*3./(4.*math.pi), 1./3.)
    sphere_equiv_area = 4.*math.pi*(sphere_equiv_radius ** 2)
    sphericity = sphere_equiv_area / area
    #cylinder vol = pi r^2 h
    cylinder_equiv_radius = math.sqrt(vol / tum_ld.GetWorldSize()[2]  / math.pi)
    cylinder_equiv_area = 2.*math.pi*cylinder_equiv_radius*tum_ld.GetWorldSize()[2]
    cylindericity = cylinder_equiv_area / area

    d = dict(time = time, area = area, volume = vol,
             radius=radius, sphericity = sphericity, cylindericity = cylindericity,
             sphere_equiv_radius = sphere_equiv_radius, sphere_equiv_area = sphere_equiv_area,
             cylinder_equiv_radius = cylinder_equiv_radius, cylinder_equiv_area = cylinder_equiv_area
             )
    pprint(d)

    print 'geometry data for %s/%s' % (statedata.file.filename, statedata.name)
    g = dstgroup.recreate_group('geometry')
    myutils.hierarchy_to_hdf(g,'.',d)


  if 1:
    print 'regenerating radial for %s/%s' % (statedata.file.filename, statedata.name)
    dstdata = dstgroup.recreate_group('radial')

    ## radial data; vessels
    vessels = krebsutils.read_vesselgraph(statedata['vessels'], ['position',  'flags', 'shearforce', 'radius', 'flow', 'maturation'])
    vessels = vessels.get_filtered(edge_indices = (vessels.edges['flags'] & krebsutils.CIRCULATED != 0))
    sample_length = 50.

    far = 1.e10

    s = plotVessels.generate_samples(vessels, 'position', 'nodes', sample_length)
    dist = krebsutils.sample_field(s, distmap, tum_ld, linear_interpolation=True, extrapolation_value = far)
    rad = krebsutils.sample_field(s, radialmap, tum_ld, linear_interpolation=True, extrapolation_value = far)
    del s

    dstdata.recreate_group('vs_r').create_dataset('bins', data = np.average((bins_rad[:-1],bins_rad[1:]), axis=0))
    dstdata.recreate_group('vs_dr').create_dataset('bins', data = np.average((bins_dist[:-1],bins_dist[1:]), axis=0))
#    dstdata['vs_r'] = dict(bins = np.average((bins_rad[:-1],bins_rad[1:]), axis=0))
#    dstdata['vs_dr'] = dict(bins = np.average((bins_dist[:-1],bins_dist[1:]), axis=0))

#    def add_histogram_data(dst, name, (s_avg, s_std, s_sqr)):
#      dst[name] = dict(avg = s_avg, std = s_std, sqr = s_sqr)

    for name in ['shearforce', 'radius', 'flow', 'maturation' ]:
      s = plotVessels.generate_samples(vessels, name, 'edges', sample_length)
      myutils.MeanValueArray.fromHistogram1d(bins_rad, rad, s).write(dstdata['vs_r'], name)
      myutils.MeanValueArray.fromHistogram1d(bins_dist, dist, s).write(dstdata['vs_dr'], name)
#      d = myutils.scatter_histogram(rad, s, bins_rad, 1.)
#      add_histogram_data(dstdata['vs_r'], name, d)
#      d = myutils.scatter_histogram(dist, s, bins_dist, 1.)
#      add_histogram_data(dstdata['vs_dr'], name, d)

    def make_mvd(flavor, s, bins, smap, mask_bound):
      a = myutils.MeanValueArray.fromHistogram1d(bins, s, np.ones_like(s))
      b = myutils.MeanValueArray.fromHistogram1d(bins, smap.ravel(), np.ones_like(smap.ravel()))
      a.cnt = b.cnt.copy()
      a.sum *= sample_length/cellvol
      a.sqr *= a.sum**2
      a.write(dstdata[flavor], 'mvd')
#      d = myutils.scatter_histogram(xdata, np.ones_like(xdata), bins=bins, xdata2 = np.ravel(xdata2))
#      m = dstdata[flavor]['bins'] >  mask_bound
#      for x in d[:2]:
#        x *= sample_length/cellvol
#        x.mask |= m
#      for x in d[2:]:
#        x *= (sample_length/cellvol)**2
#        x.mask |= m
#      add_histogram_data(dstdata[flavor], 'mvd', d)

    make_mvd('vs_r', rad, bins_rad, radialmap, axial_max_rad)
    make_mvd('vs_dr', dist, bins_dist, distmap, 1000.)
#    pyplot.errorbar(dstdata['vs_dr']['bins'], dstdata['vs_dr']['mvd_avg'], yerr=dstdata['vs_dr']['mvd_std'])
#    pyplot.show()

    del rad, dist

    for name in ['conc', 'sources', 'ptc', 'press']:
      a = np.ravel(np.asarray(tum_grp[name]))
      myutils.MeanValueArray.fromHistogram1d(bins_rad, np.ravel(radialmap), a).write(dstdata['vs_r'], name)
      myutils.MeanValueArray.fromHistogram1d(bins_dist, np.ravel(distmap), a).write(dstdata['vs_dr'], name)
#      d = myutils.scatter_histogram(np.ravel(distmap), np.ravel(np.asarray(tum_grp[name])), bins_dist, 1.)
#      add_histogram_data(dstdata['vs_dr'], name, d)
#      d = myutils.scatter_histogram(np.ravel(radialmap), np.ravel(np.asarray(tum_grp[name])), bins_rad, 1.)
#      add_histogram_data(dstdata['vs_r'], name, d)

    # velocities, projected radially outward
    vel_field = tuple(np.asarray(tum_grp['vel'][:,:,:,i]) for i in range(3))
    vnorm = np.zeros(distmap.shape, dtype=np.float32)
    dgrad = krebsutils.field_gradient(radialmap, spacing=tum_ld.scale)
    for vv, g in zip(vel_field, dgrad):
      vnorm += vv*g
    del vel_field, dgrad

    phi_tumor = ptc * np.asarray(tum_grp['conc'])

    oxy = np.asarray(statedata['fieldOxy'])

    gf = np.array(statedata['fieldGf'])

    for name, field in [('phi_tumor', phi_tumor),
                        ('vel', vnorm),
                        ('oxy', oxy),
                        ('gf', gf)]:
      a = np.ravel(field)
      myutils.MeanValueArray.fromHistogram1d(bins_rad, np.ravel(radialmap), a).write(dstdata['vs_r'], name)
      myutils.MeanValueArray.fromHistogram1d(bins_dist, np.ravel(distmap), a).write(dstdata['vs_dr'], name)

#      d = myutils.scatter_histogram(np.ravel(distmap), np.ravel(field), bins_dist, 1.)
#      add_histogram_data(dstdata['vs_dr'], name, d)
#      d = myutils.scatter_histogram(np.ravel(radialmap), np.ravel(field), bins_rad, 1.)
#      add_histogram_data(dstdata['vs_r'], name, d)


def averaged(data):
  groups, path = data
  avg = [ np.asarray(g[path]) for g in groups ]
  std = np.std(avg, axis=0)
  avg = np.average(avg, axis=0)
  return avg, std

def averaged_global(data, name):
  groups_by_time, times = data
  res = [ averaged((groups_by_time[t], name)) for t in times ]
  return np.asarray(res).transpose()




def averaged_radial(data, flavor, names):
  groups_by_time, times = data
  all_res = {}
  for name in names:
    for t in times:
      avg, std, sqr = [], [], []
      for group in groups_by_time[t]:
        g = group['radial/%s/%s' % (flavor, name)]
        avg.append(myutils.MeanValueArray.read(g).avg)
        #std.append(np.asarray(g['std']))
        #sqr.append(np.asarray(g['sqr']))
      std = np.ma.std(avg, axis=0)
      avg = np.ma.average(avg, axis=0)
      all_res[name, t] = avg, std
  return all_res



def plot_many(filenames, pdfpages):
  groups_by_time = defaultdict(list)
  files = [ h5py.File(fn, 'r') for fn in filenames ]
  for f in files:
    groups = myutils.getTimeSortedGroups(f['.'])
    for g in groups:
      groups_by_time[round(g.attrs['time'])].append(g)

  times = groups_by_time.keys()
  times.sort()
  times = np.asarray(times)

  radius, radius_std = averaged_global((groups_by_time, times), 'geometry/radius')
  volume, volume_std = averaged_global((groups_by_time, times), 'geometry/volume')
  sphere_equiv_radius, sphere_equiv_radius_std = averaged_global((groups_by_time, times), 'geometry/sphere_equiv_radius')
  #sphere_equiv_radius, sphere_equiv_radius_std = averaged_global((groups_by_time, times), 'geometry/cylinder_equiv_radius')
  sphericity, sphericity_std = averaged_global((groups_by_time, times), 'geometry/sphericity')

  # estimate velocity
  data = np.asarray((times, radius))
  if len(data[0])>1:
    from scipy.optimize import leastsq
    func = lambda p, x, y: (x*p[0]+p[1]-y)
    p, success = leastsq(func, (1, 0), args=(data[0], data[1]))
    velocity = p[0]
    print 'estimated velocity: %f'  % velocity
    print 'fit params: %s' % str(p)
  else:
    velocity = 0.

  if 1:
    #plot by time
    fig, axes = pyplot.subplots(2,2, figsize = (mastersize[0]*0.5, mastersize[0]*0.5))
    mpl_utils.subplots_adjust_abs(fig, left = mpl_utils.mm_to_inch*13,
                                  right = -mpl_utils.mm_to_inch*5,
                                  top  = -mpl_utils.mm_to_inch*5,
                                  bottom = mpl_utils.mm_to_inch*10,
                                  hspace = mpl_utils.mm_to_inch*30,
                                  wspace = mpl_utils.mm_to_inch*40,)
    axes = axes.ravel()

    def plt_rad(ax):
      ax.set(ylabel='[mm]', xlabel='t [h]')
      mpl_utils.errorbar(ax, times, 1.e-3*radius, yerr=1.e-3*radius_std, label = 'r', lw = 0., marker = 'x', color = 'k', markersize = 5.)
      label = u'$r_0 + v_{fit} t$\n$v_{fit} = %s$ [\u03BCm/h]' % f2s(velocity)
      ax.plot(times, 1.e-3*(p[1] + p[0] * times), label = label, color = 'r')
      ax.legend()
      #ax.text(0.6, 0.2, r'$v_{fit} = %s$' % f2s(velocity), transform = ax.transAxes)

    def plt_vol(ax):
      ax.set(ylabel='volume', xlabel='t [h]')
      ax.errorbar(times, 1.e-9*volume, yerr=1.e-9*volume_std)

    def plt_sprad(ax):
      ax.set(ylabel='sphere equiv. radius', xlabel='t [h]')
      ax.errorbar(times, 1.e-3*sphere_equiv_radius, yerr=1.e-3*sphere_equiv_radius_std)

    def plt_sphereicity(ax):
      ax.set(ylabel='sphericity', xlabel='t [h]')
      ax.errorbar(times, sphericity, yerr=sphericity_std)

    for ax, func in zip(axes, [plt_rad, plt_vol, plt_sprad, plt_sphereicity]):
      l = MaxNLocator(nbins = 4)
      ax.xaxis.set_major_locator(l)
      func(ax)

    pdfpages.savefig(fig)


  times = times[0::2]
  radial = averaged_radial((groups_by_time, times), 'vs_dr', ['phi_tumor', 'mvd', 'radius', 'shearforce', 'flow', 'sources', 'vel', 'oxy', 'maturation'])
  bins = np.asarray(groups_by_time[times[0]][0]['radial/vs_dr/bins'])
  bins *= 1./1000.
  xlim = -1.5, 1.0
  mask = np.logical_and(bins < xlim[1], bins >= xlim[0])


  def plot_times(ax, name, **kwargs):
    f = kwargs.pop('value_prefactor', 1.)
    colors = 'rgbmk'
    markers = 'os<>d'
    for i, t in enumerate(times):
      avg, std = radial[name, t]
      avg, std = avg[mask], std[mask]
      #ax.errorbar(bins[mask], f*avg, yerr=f*std, label = 't = %s' % f2s(t), **kwargs)
      mpl_utils.errorbar(ax, bins[mask], f*avg, yerr=f*std, label = 't = $%s$' % f2s(t),
                         marker = markers[i], color = colors[i], every = 5,
                         **kwargs)

  def text1(ax, txt):
    ax.text(0.95, 0.9, txt, ha = "right", transform = ax.transAxes)

  def text2(ax, txt):
    ax.text(0.01, 0.9, txt, ha = "left", transform = ax.transAxes)

  def mkfig(nrows, ncols):
    fig, axes = mpl_utils.subplots_abs_mm((mastersize[0]/mpl_utils.mm_to_inch*0.5*ncols,
                                 mastersize[0]/mpl_utils.mm_to_inch*0.2*nrows),
                                          nrows, ncols,
                                          10, 20, a4size[0]/mpl_utils.mm_to_inch*0.38, a4size[0]/mpl_utils.mm_to_inch*0.16, 15, 5)
    return fig, axes

#  fig, axes = pyplot.subplots(4, 2, figsize = (mastersize[0], mastersize[0]*0.25*4.))
#  mpl_utils.subplots_adjust_abs(fig, left = mpl_utils.mm_to_inch*20,
#                                right = -mpl_utils.mm_to_inch*10,
#                                top  = -mpl_utils.mm_to_inch*5,
#                                bottom = mpl_utils.mm_to_inch*10,
#                                hspace = mpl_utils.mm_to_inch*30,)

  def plt_mvd(ax):
    # mvd can be written as N^3 * L0/N * 3 / V = (V=L0^3) ... = N^2 / L0^2
    ax.set(ylabel = ur'$\times 10^3$ [\u03BCm$^{-2}$]', xlim = xlim)
    plot_times(ax,'mvd', value_prefactor = 1e3)
    text1(ax, r'$L/V$')
    ax.legend(loc = mpl_utils.loc.lower_left, frameon = True)

  def plt_rad(ax):
    ax.set(ylabel = ur'[\u03BCm]', xlim = xlim)
    plot_times(ax,'radius')
    text1(ax, r'$r_v$')

  def plt_vel(ax):
    ax.set(ylabel = ur'[\u03BCm/h]', xlim = xlim)
    text1(ax, r'$v_\phi$')
    plot_times(ax,'vel')

  def plt_tum(ax):
    ax.set(ylabel = ur'', xlim = xlim, ylim = (-0.1, 0.7))
    plot_times(ax,'phi_tumor')
    text1(ax, r'$\phi_t$')

  def plt_oxy(ax):
    ax.set(ylabel = ur'', xlim = xlim)
    plot_times(ax,'oxy')
    text1(ax,r'$c_o$')

  def plt_sf(ax):
    ax.set(ylabel = ur'[Pa]', xlim = xlim)
    plot_times(ax,'shearforce', value_prefactor = 1e3)
    text1(ax,r'$f_v$')

  def plt_wall(ax):
    ax.set(ylabel = ur'[\u03BCm]', xlim = xlim)
    plot_times(ax,'maturation')
    text1(ax,r'$w_v$')

  def plt_qv(ax):
    ax.set(ylabel = ur'[\u03BCm]', xlim = xlim)
    ax.set(xlabel = ur'$\theta$ [mm]')
    plot_times(ax,'flow')
    text1(ax,r'$q_v$')

  fig, axes = mkfig(3, 2)
  plt_mvd(axes[0,0])
  plt_rad(axes[0,1])
  plt_vel(axes[1,0])
  plt_oxy(axes[1,1])
  plt_wall(axes[2,0])
  plt_sf(axes[2,1])

  for ax in axes[2,:]: ax.set(xlabel = ur'$\theta$ [mm]')

  axes = axes.ravel()
  for i, ax in enumerate(axes):
    ax.grid(linestyle=':', linewidth=0.5, color=gridcolor)
    mpl_utils.add_crosshair(ax, (0,0), color = gridcolor)
    if not ax.get_xlabel():
      ax.set(xticklabels = [])
    text2(ax, '(%s)' % 'abcdefghij'[i])

  pdfpages.savefig(fig)






if __name__ == '__main__':
  filenames = sys.argv[1:]
  if not myutils.is_measurement_file(filenames[0]):
    for fn in filenames:
      with h5py.File(fn, 'r') as f:
        with myutils.MeasurementFile(f, h5files) as fdst:
          grps = myutils.getTimeSortedGroups(f['.'], "out")
          for grp in grps:
            dstgroup = myutils.require_snapshot_group_(fdst, grp)
            generate_data(grp, dstgroup)
  else:
    rc = matplotlib.rc
    rc('figure', **{'subplot.left' : 0.15,
                    'subplot.right' : 1.-0.15,
                    'subplot.bottom' : 0.2,
                    'subplot.top' : 1.-0.1,
                    'subplot.wspace' : 0.2,
                    'subplot.hspace' : 0.2})
    fnmeasure = commonOutputName(filenames)
    with mpl_utils.PdfWriter(fnmeasure+'.pdf') as pdfpages:
      plot_many(filenames, pdfpages)
