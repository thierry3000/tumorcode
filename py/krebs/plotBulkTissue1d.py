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
import h5py
import os,sys
import posixpath
import numpy as np
import collections
import matplotlib
import matplotlib.pyplot as plt
import mpl_utils
import extensions

from mystruct import Struct
import myutils


dataconfigs = [
    Struct(name='kdiff_cells', rng='auto', title=r'kdiff_cells'),
    Struct(name='kdiff_obst', rng='auto', title=r'kdiff_obst'),
    Struct(name='necro', rng=(-0.1,1.1), title=r'necro'),
    Struct(name='totalconc', rng=None, title=r'vol. ratio $\phi + m$'),
    Struct(name='conc', rng=(-0.1,1.1), title=r'vol. ratio $\phi$'),
    Struct(name='conc_necro', rng=(-0.1,1.1), title=r'vol. ratio $d$'),
    Struct(name='obstacle', rng=(-0.1,1.1), title=r'obst. ratio $m$'),
    Struct(name='oxy', rng=(-0.1,1.1), title=r'O$_2\, c$'),
    Struct(name='oxy_sources', rng=(-0.1,1.1), title=r'$\Gamma_c$'),
    Struct(name='fieldOxy', rng=(-0.1,1.1), title=r'O$_2$'),
    Struct(name='vel_0', rng='auto', coords='xface', title=r'$v_x$'),
    Struct(name='ls', rng='zero-centered', title=r'$\theta$'),
    Struct(name='vel', rng='zero-centered', title=r'$v_x$'),
    Struct(name='press', rng='auto', title=r'pressure $\Sigma$'),
    Struct(name='ptc', rng=(0., 1.), title=r'$\Theta_\epsilon(\theta)$'),
    Struct(name='sources', rng='zero-centered', title=r'$\Gamma_\phi \times 10^3 $', scale= 1000.),
    Struct(name='vessel_volume_fraction', rng=(0.,1.), title=r'vol. ratio $v$')
  ]
dataconfigs = dict((d.name, d) for d in dataconfigs)



class Filedata1d:
  def __init__(self, fn):
    self.fn = fn
    self.f  = f = h5py.File(fn, 'r')
    self.ld = ld = f['field_ld'].attrs
    self.groupnames = getSortedTimeGroupNames(f)
    self.xcoords = (np.arange(ld['SIZE'][0]) * ld['SCALE'] + ld['WORLD_OFFSET'][0])*0.001
    self.xfacecoords = ((np.arange(ld['SIZE'][0]+1)-0.5) * ld['SCALE'] + ld['WORLD_OFFSET'][0])*0.001
    self.datanames = f[self.groupnames[-1]].keys()
    self.bounds_cache = {}

  def findBoundsOverTime(self, name):
    try:
      return self.bounds_cache[name]
    except KeyError:
      dc = dataconfigs[name]
      bounds = [], []
      for gn in self.groupnames[1:]:
        ds = np.asarray(self.f[gn][name])
        bounds[0].append(ds.min())
        bounds[1].append(ds.max())
      bounds = min(*bounds[0])*dc.get('scale', 1.), max(*bounds[1])*dc.get('scale', 1.)
      self.bounds_cache[name] = bounds
      return bounds

  def addToPlot(self, plt, name, groupname, *args, **kwargs):
    dc = dataconfigs[name]
    group = self.f[groupname]
    if name == 'totalconc':
      arr = np.asarray(group['conc'][:,0,0])
      if 'obstacle' in self.datanames:
        arr += np.asarray(group['obstacle'][:,0,0])
      if 'conc_necro' in self.datanames:
        arr += np.asarray(group['conc_necro'][:,0,0])
    else:
      arr = np.asarray(group[dc.name][:,0,0]) * dc.get('scale', 1.)
    p, = plt.plot(
        self.xcoords if dc.get('coords','x')=='x' else self.xfacecoords,
        arr,
        *args, **kwargs
      )
    if dc.rng:
      if dc.rng == 'zero-centered':
        a, b = self.findBoundsOverTime(name)
        v = max(abs(a), abs(b))
        #v = np.max(np.abs(arr))
        plt.set_ylim(-v, v)
      elif dc.rng == 'auto':
        a, b = self.findBoundsOverTime(name)
        a -= 0.1*(b-a)
        b += 0.1*(b-a)
        plt.set_ylim((a-0.1*(b-a)), b+0.1*(b-a))
      else:
        plt.set_ylim(*dc.rng)
    return p



def add_twinplot(fig, subplotparams, fd, groupname, datanames):
  ax1 = fig.add_subplot(*subplotparams)
  ax2 = ax1.twinx()
  p2 = fd.addToPlot(ax2, 'ptc', groupname, 'k--')
  ax2.yaxis.set_visible(False)
  p1 = []
  for n in datanames:
    p1 += [fd.addToPlot(ax1, n, groupname)]
  return ax1, ax2, p1, p2


def plotForPresent(fd, writer, data):
  dpi = 80
  rc = matplotlib.rc
  rc('figure', figsize = list(0.75/dpi*np.asarray([1024, 768])), dpi=dpi)
  rc('font', size = 12.)
  rc('axes') #, titlesize = 12., labelsize = 12.)
  rc('pdf', compression = 6, fonttype = 42)
  rc('figure', **{'subplot.left' : 0.04,
                  'subplot.right' : 1.-0.04,
                  'subplot.bottom' : 0.06,
                  'subplot.top' : 1.-0.06,
                  'subplot.wspace' : 0.1,
                  'subplot.hspace' : 0.1})
  rc('savefig', facecolor='none', edgecolor='none', dpi=dpi)
  rc('font', **{'family':'sans-serif'}) #,'sans-serif':['Helvetica']})
  rc('path', simplify_threshold=0.01)
  #rc('lines', linewidth=2.)

  for groupname in fd.groupnames:
    print 'plt', groupname
    group = fd.f[groupname]

    fig = matplotlib.figure.Figure()
    fig.suptitle('$t$=%s h' % (myutils.f2s(group.attrs['time'])))



    #datanames = set(['conc', 'obstacle', 'vel_0', 'press', 'sources', 'ls'])
    if 1:
      ax1, ax2, (p1,p2), p3 = add_twinplot(fig, (2,2,1), fd, groupname, ['conc','obstacle'])
      ax1.legend()
      ax1.xaxis.set_visible(False)
    else:
      ax1, ax2, (p1,), p3 = add_twinplot(fig, (2,2,1), fd, groupname, ['conc',])
      ax1.legend()
      ax1.xaxis.set_visible(False)

    ax1, ax2, (p1,), _ = add_twinplot(fig, (2,2,2), fd, groupname, ['press'])
    ax1.legend([p1], [dataconfigs['press'].title])
    ax1.xaxis.set_visible(False)

    ax1, ax2, (p1,), _ = add_twinplot(fig, (2,2,3), fd, groupname, ['vel_0'])
    ax1.legend([p1], [dataconfigs['vel_0'].title])

    ax1, ax2, (p1,), _ = add_twinplot(fig, (2,2,4), fd, groupname, ['sources', 'oxy_sources'])
    ax1.legend([p1], [dataconfigs['sources'].title, dataconfigs['oxy_sources'].title])


    pdf.savefig(fig)




class Subplotter(object):
  def __init__(self, fd, groupname, layout):
    self.layout = layout
    self.plot_num = 1
    self.fd = fd
    self.groupname = groupname
    self.fig = plt.figure()

  def subplot(self, (x, y), twinplots=False, **kwargs):
    lx, ly = self.layout

    ax1 = self.fig.add_subplot(ly, lx, x+lx*y+1, **kwargs)

    if twinplots: ax2 = ax1.twinx()

    if y < ly-1:
      ax1.set(xticklabels=[])
    else:
      ax1.set_xlabel(r'$x [mm]$')

    self.plot_num += 1

    ax1.set(xlim=(self.fd.xcoords[0], self.fd.xcoords[-1]))
    ax1.grid(linestyle=':', linewidth=0.5, color=(0.7,0.7,0.7))

    if not twinplots: return ax1

    ax2.yaxis.set_visible(False)

    return ax1, ax2


  def plt(self, ax, name, *args, **kwargs):
    return self.fd.addToPlot(ax, name, self.groupname, *args, **kwargs)



def simplePlot(fd, groupname, title):
  group = fd.f[groupname]
  datanames = ['conc',
               #'necro',
               #'ncells',
               'vel_0',
               #'press',
               'sources',
               #'ls',
               #'obstacle',
               'oxy',
               ]
  datanames = [ n for n in datanames if n in fd.datanames ]
  #layout = (1, len(datanames))
  layout = 2, 2
  layout_pos = [ (0, 0), (1, 0), (0, 1), (1, 1) ]

  subplotter = Subplotter(fd, groupname, layout)
  subplotter.fig.suptitle('%s $t$=%0.f h' % (title, group.attrs['time']))

  for i, name in enumerate(datanames):
    ax1, ax2 = subplotter.subplot(layout_pos[i], True) #ylabel=dataconfigs[name].title)

    p1 = subplotter.plt(ax1, name, 'g-', label=dataconfigs[name].title)
    p2 = subplotter.plt(ax2, 'ptc', 'b--')
    # fucking n00b programming with if statements, checking what to plot and setting special options!
    if name == 'conc':
      #ax1.set(ylabel = r'$vol. ratios$')
      pp = [p1]
      p1.set(label=r'$t+n+d$')
      p2.set(label=r'$\theta_t$')
      if 'necro' in fd.datanames:
        p3 = subplotter.plt(ax1, 'necro', 'm-', label=r'$d$')
        pp.append(p3)
      if 'obstacle' in fd.datanames:
        p4 = subplotter.plt(ax1, 'obstacle', 'r-', label=r'ecm/vessels')
        pp.append(p4)
      if 'necro' in fd.datanames or 'obstacle' in fd.datanames:
        p5 = subplotter.plt(ax1, 'totalconc', 'k-', label=r'not water')
        pp.append(p5)
      pp.append(p2)
      ax1.legend(*zip(*[(p, p.get_label()) for p in pp]), loc = 1)

    elif name == 'vel_0':
      ax3 = ax1.twinx()
      ax3.yaxis.set_visible(False)
      ax1.set(xlim=(subplotter.fd.xcoords[0], subplotter.fd.xcoords[-1]), label=dataconfigs['press'].title)
      p3 = subplotter.plt(ax3, 'press', 'm-', label=dataconfigs['press'].title)
      ax1.legend([p1,p3],[p1.get_label(), p3.get_label()], loc =1)

    elif name == 'oxy' and 'oxy_sources' in fd.datanames:
      p3 = subplotter.plt(ax1, 'oxy_sources', 'm-', label=dataconfigs['oxy_sources'].title)
      ax1.legend(loc = 4)

    else:
      ax1.legend(loc = 4)

  fig = subplotter.fig
  del subplotter
  return fig




def plotit(fn):
  out_fn = os.path.splitext(fn)[0]
  fd = Filedata1d(fn)
  if 1:
    data = [
      (g.attrs['time'], g.attrs['tumor_radius']) for g in [
        fd.f[gn] for gn in fd.groupnames
      ]
    ]
    if len(fd.groupnames) >= 2:
      maxv = max(np.amax(np.abs(fd.f[gn]['vel_0'])) for gn in fd.groupnames[2:])

      data = np.asarray(data).transpose()
      from scipy.optimize import leastsq
      func = lambda p, x, y: (x*p[0]+p[1]-y)
      p, success = leastsq(func, (1, 0), args=(data[0], data[1]))
      estv = p[0]
      print 'estimated velocity: %f, max vel.: %f'  % (estv, maxv)
    else:
      estv, maxv = 0, 0

  if 0:
    rc = matplotlib.rc
    rc('figure', figsize=(8.3, 11.7))
    rc('font', size = 13)
    #rc('axes', titlesize = 12., labelsize = 12.)
    rc('legend', fontsize = 'small')
  else:
    dpi = 80
    rc = matplotlib.rc
    rc('figure', figsize = list(0.75/dpi*np.asarray([1024, 768])), dpi=dpi)
    rc('font', size = 12.)
    rc('axes') #, titlesize = 12., labelsize = 12.)
    rc('pdf', compression = 6, fonttype = 42)
    rc('figure', **{'subplot.left' : 0.04,
                    'subplot.right' : 1.-0.04,
                    'subplot.bottom' : 0.06,
                    'subplot.top' : 1.-0.06,
                    'subplot.wspace' : 0.1,
                    'subplot.hspace' : 0.1})
    rc('savefig', facecolor='none', edgecolor='none', dpi=dpi)
    rc('font', **{'family':'sans-serif'}) #,'sans-serif':['Helvetica']})
    rc('path', simplify_threshold=0.01)

  with mpl_utils.PdfWriter(out_fn) as pdf:
    for g in fd.groupnames:
      print fn, g
      #title = r'%s ($v_e = %s, v_{max} = %s$)' % (splitext(basename(fn))[0], myutils.f2s(estv), myutils.f2s(maxv))
      title = ''
      fig = simplePlot(fd, g, title)
      pdf.savefig(fig)
      plt.close()





#def comparePlot(fds, fdtitles, groupname, datanames, title):
#  firstgroup = fds[0].f[groupname]
#  fig = matplotlib.figure.Figure(figsize=(8.3, 11.7))
#  fig.suptitle('%s $t$=%0.f h' % (title, firstgroup.attrs['time']))
#  for i, is_first, is_last, name in myutils.enumerate_seq(datanames):
#    ax1 = fig.add_subplot(len(datanames),1,i+1)
#    for fd, fdtitle in zip(fds, fdtitles):
#      p1 = fd.addToPlot(ax1, name, groupname)
#      p1.set_label(fdtitle)
#    ax1.set_ylabel(dataconfigs[name].title)
#    if not is_last:
#      ax1.xaxis.set_visible(False)
#    if is_last:
#      ax1.set_xlabel(r'$x [mm]$')
#    if is_first:
#      ax1.legend()
#  return fig



if __name__ == "__main__":
  fn = sys.argv[1]
  plotit(fn)