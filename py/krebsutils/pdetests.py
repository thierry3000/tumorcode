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
if __name__ == '__main__':
  sys.path.append(join(dirname(__file__),'..'))
import h5py
import numpy as np
import math
import itertools
import collections
import myutils
import extensions
import pprint

import krebsutils
#imports_ = 'convectionTest convectionLevelset levelsetRedistancing zalesakDisk rotationalVelocityField calcCurvature surfaceTensionForceTest EST_testOperatorApplication EST_testOperatorSolve steppersConvergenceTest'.split()
imports_ = 'convectionTest convectionLevelset levelsetRedistancing zalesakDisk rotationalVelocityField calcCurvature surfaceTensionForceTest'.split()
# load functions from libkrebs_ into local namespace
# i do this because libkrebs can be the debug or release version of the library dependign on the debug switch in the python interpreter
# so do not directly use 'import libkrebs_.blabla'!
locals().update( (f,getattr(krebsutils.libkrebs, f)) for f in imports_)

import myutils
from dicttoinfo import dicttoinfo, Vec

try:
  import mpl_utils
  import matplotlib
  matplotlib.use('Qt4Agg')
  import matplotlib.pyplot as pyplot
  has_matplotlib = True
except ImportError:
  has_matplotlib = False


#############################################
## utils
##############################################

def make_x(ld):
  box = ld.box(ravel = True)
  x = np.empty((box[1]-box[0]+1,))
  for i in xrange(box[0], box[1]+1):
    x[i] = ld.LatticeToWorld((i,0,0))[0]
  return x

def make_ld(params):
  n = params['ncells']
  h = (params['x1']-params['x0'])/n
  ld = krebsutils.LatticeDataQuad3d((0, n-1, 0, 0, 0, 0), h)
  ld.SetCellCentering((True, False, False))
  ld.SetOriginPosition  ((params['x0'],0.,0.))
  return ld


class Results(object):
  def __init__(self, name):
    f = h5py.File(name+'.h5','r')
    self.ld = ld = krebsutils.read_lattice_data_from_hdf(f['field_ld'])
    self.groups = myutils.getTimeSortedGroups(f['.'], 'out')
    self.f = f
  def __len__(self):
    return len(self.groups)
  def __getitem__(self, key):
    return self.groups[key]
  def __iter__(self):
    return self.groups.__iter__()


def plot1dresult(axes, res, index, style, **kwargs):
  x = make_x(res.ld)
  g = res[index]
  t, u = g.attrs['time'], np.asarray(g['conc'])
  axes.plot(x, u[:,0,0], style, **kwargs)


#def sampleFxt1d(t, xr, func):
#  a = np.empty_like(xr)
#  for i, x in enumerate(xr):
#    a[i] = func(x, t)
#  return a


def _ncycle(iterable):
    """
    Method to create a generator from an iterable.  It keeps the
    current position of the iterable in memory.  Each time the
    next() method for the iterable is called, it will return the
    next item.  If there are no more items, it will cycle to the
    first item.
    """

    for item in itertools.cycle(iterable):
        yield item

markers = _ncycle("osv<>^hdx*")

##############################################
## 0d stepper scheme test
##############################################

def runStepperConvergenceTests():
  out_fn = "stepper_test"

  steppersConvergenceTest(out_fn + '.h5')

  pdf = mpl_utils.PdfWriter(out_fn + ".pdf")

  f = h5py.File(out_fn + '.h5', 'r')

#  model_lambda = f.attrs['lambda']
  bytime = collections.defaultdict(list)
  bymethod = collections.defaultdict(list)
  for g in f['/'].itervalues():
    if g.name == '/exact':
      exact = g
    else:
      bytime[g.attrs['dt']].append(g)
      bymethod[g.attrs['method']].append(g)

  for dt, groups in sorted(bytime.iteritems(), key = lambda (k, v): k):
    if dt < exact['x'][-1]*0.01:
      continue
    fig = pyplot.figure(figsize=(11.7, 8.3))
    fig.suptitle('dt = %f' % dt)
    plot = fig.add_subplot(1,1,1)
    plot.set_yscale("log")
    plot.set_xlabel('x')
    plot.set_ylabel('y')
    for marker, g in zip(markers, groups):
      plot.plot(np.asarray(g['x']), np.asarray(g['y']), label = g.attrs['method'], marker=marker)
    plot.plot(np.asarray(exact['x']), np.asarray(exact['y']), label = 'exact', linewidth=2)
    plot.legend(loc = 'upper left')
    pdf.savefig(fig)

  err = collections.defaultdict(collections.defaultdict)
  for g in f['/'].itervalues():
      if g.name == '/exact': continue
      val = abs(g['y'][-1] - exact['y'][-1])
      m, t = g.attrs['method'], g.attrs['dt']
      err[m][t] = val
      print m, t, val
  for k, d in err.iteritems():
    err[k] = np.asarray(sorted(d.iteritems(), key = lambda (k, v): k)).transpose()

  fig = pyplot.figure(figsize=(11.7, 8.3))
  plot = fig.add_subplot(1,1,1)
  plot.set_yscale("log")
  plot.set_xscale("log")
  plot.set_xlabel('dt')
  plot.set_ylabel('e')
  for marker, (method, (t, y)) in zip(markers, err.iteritems()):
    plot.plot(t, y, label = method, marker=marker)
  plot.legend(loc = 'lower right')
  pdf.savefig(fig)

  pdf.close()



##############################################
## numerical advection scheme analysis
##############################################

def compute1dAdvectionModel(name, params, conc_func, velocity_func):
  params['fn_out'] = name
  f_conc = lambda x,y,z: conc_func(x)
  f_vel = lambda x,y,z: (velocity_func(x), 0, 0)
  convection(dicttoinfo(params),
                               make_ld(params),
                               dict(conc_profile = f_conc,
                                    velocity_profile = f_vel))


def runConvergenceCheck1d(conc_func, vel_func, analytic_res, params):
  """
    plot error from solution of a 1d linear advection equation
  """
  ncells_list = np.asarray([ 32, 64, 128, 256, 512, 1024 ])
  res = []
  # check convergence of an advection scheme
  for num, ncells in enumerate(ncells_list):
    params['ncells'] = ncells
    compute1dAdvectionModel('convconvergence%i' % num, params, conc_func, vel_func)
    results = Results('convconvergence%i' % num)
    u = results[-1]['conc'][:,0,0]
    t = results[-1].attrs['time']
    x = make_x(make_ld(params))
    u_exact = np.frompyfunc(analytic_res, 2, 1)(x, t)
    i0, i1 = int(ncells*(0.5-0.10)),int(ncells*(0.5+0.05))
    #print x[i0], x[i1]
    # the error is the L2 norm of the difference from teh exact solution in the center part
    # where it is smooth. At the discontinuities it degenerates to 1st order.
    # one can check that, too by removing [i0:i1]
    err = np.linalg.norm((u-u_exact)[i0:i1])*((x[-1]-x[0])/ncells)
    res.append(err)
    pyplot.plot(x, u-u_exact, label='%i' % num)
  pyplot.legend()
  pyplot.show()
  # plot the error together with some power laws
  pyplot.gca().set(yscale="log", xscale='log')
  pyplot.plot(ncells_list, res)
  pyplot.plot(ncells_list, np.power(ncells_list, -3.), label='3')
  pyplot.plot(ncells_list, np.power(ncells_list, -2.), label='2')
  pyplot.plot(ncells_list, np.power(ncells_list, -1.), label='1')
  pyplot.legend()
  pyplot.show()


def runMethodsComparison1d(conc_func, vel_func, analytic_res, params):
  x = make_x(make_ld(params))
  u0 = np.frompyfunc(analytic_res, 2, 1)(x, 0.)
  u1 = np.frompyfunc(analytic_res, 2, 1)(x, 100.)
  pyplot.plot(x, u0, '-k')
  pyplot.plot(x, u1, '-k')

  styles = zip(['r', 'g', 'b', 'y' ], ['*','x','+','<'])
  for i, method in enumerate([ 'eno2', 'weno3', 'eno3', 'weno5' ]):
    params['advection_method'] = method
    compute1dAdvectionModel('convtest_'+method, params, conc_func, vel_func)
    results = Results('convtest_'+method)
    c, m = styles[i]
    plot1dresult(pyplot.gca(), results, -1, '-'+c+m, markeredgecolor=c, label = method)

  mass0 = np.sum(results[0]['conc'])
  mass1 = np.sum(results[-1]['conc'])
  print 'mass start->end: %f -> %f' % (mass0, mass1)

  pyplot.legend()
  pyplot.show()


def runAdvectionSchemeAnalysis():
  def conc_func(x):
    c, w = 80., 15.
    if not -w < x-c < w: return 0.
    return 1. + 0.25 * math.cos((x-c+w/2)/w*math.pi)
  if 0:
    vel_func = lambda x: 0.0075*x
    #analytic_res = lambda x, t: conc_func(x-vel_func(x)*t)
    def analytic_res(x,t):
      # x(t,x0) = x0 * exp(a * t)
      # ln(x) = ln(x0 * exp(a*t)) = ln(x0) + a*t
      # x0 = x * exp(-a*t)
      x0 = x*math.exp(-0.0075*t)
      return conc_func(x0)
  if 1:
    vel_func = lambda x: -0.75
    analytic_res = lambda x,t : conc_func(x+0.75*t)

  params = dict(
    ncells = 100,
    x0 = -50.,
    x1 = 150.,
    tend = 100,
    out_intervall = 20.,
    stepper = 'rk3',
    conservative = False
  )

  runMethodsComparison1d(conc_func,  vel_func, analytic_res, params)

  # convergence test
  params['advection_method'] = "weno5"
  runConvergenceCheck1d(conc_func, vel_func, analytic_res, params)

##############################################
## levelset scheme analysis
##############################################
#see http://assorted-experience.blogspot.com/2007/07/custom-colormaps.html
cm_ls = matplotlib.colors.LinearSegmentedColormap('levelset', {
    'red' :   ((0,1,1), (0.5,0  ,1  ), (1,0,0)),
    'green' : ((0,1,1), (0.5,0  ,1  ), (1,0,0)),
    'blue' :  ((0,1,1), (0.5,0  ,1. ), (1,0,0))
  })


def runLevelsetRedistancingScheme():
    name = 'zalesak_redistanceing'
    params = dict(fn_out = 'zalesak_redistanceing')

    ld = krebsutils.LatticeDataQuad3d((0, 99, 0, 99, 0, 0), 1./100.)
    ld.SetCellCentering((True, True, False))

    f_distort = lambda w: (-1 if w<0 else 1)*(abs(w)**(1./3.))
    f_conc = lambda x,y,z: f_distort(zalesakDisk((x,y,z), 0.33, (0.5, 0.5, 0.)))

#    if not os.path.exists(name+'.h5'):
    krebsutils.set_num_threads(2)
    levelsetRedistancing(dicttoinfo(params), ld, dict(conc_profile=f_conc))

#    with mpl_utils.SinglePageWriter(name+'.png') as pdfpages:
#      results = Results(name)
#
#      fig, ax = pyplot.subplots(1,1, figsize = (800./90., 600/90.), dpi = 90.)
#      a = results.f['active']
#      a = a[:,:,a.shape[2]/2]
#      img = ax.imshow(np.asarray(a), origin='bottom', cmap = matplotlib.cm.coolwarm, extent=ld.worldBox[:4]) #vmin=-0.1, vmax=0.1
#      fig.colorbar(img)
#      pdfpages.savefig(fig, dpi = 90)
#      pyplot.close()
#
#      for grp in results[::10]:
#        a = grp['ls']
#        a = a[:,:,a.shape[2]/2]
#        amax = np.amax(np.abs(a))
#        fig, ax = pyplot.subplots(1,1, figsize = (800./90., 600/90.), dpi = 90.)
#        img = ax.imshow(np.asarray(a), origin='bottom', cmap = cm_ls, extent=ld.worldBox[:4], vmin=-amax, vmax=amax)
#        fig.colorbar(img)
#        pdfpages.savefig(fig, dpi = 90)
#        pyplot.close()


def runLevelsetAdvectionScheme():
    gradient_approx_method = 'upwind'
    stepper = 'rk3'
    advection_method = "weno5"
    velocity_type = 'shear' #'rotation'

    name = 'zalesak_%s_%s_%s_%s' % (velocity_type, gradient_approx_method, stepper, advection_method)
    params = dict(
      tend = math.pi*2.,
      out_intervall = math.pi*2./12,
      stepper = stepper,
      fn_out = name,
      advection_method = advection_method,
      gradient_approx_method = gradient_approx_method,
      num_threads = 3,
    )
    ld = krebsutils.LatticeDataQuad3d((0, 99, 0, 99, 0, 0), 1./100.)
    ld.SetCellCentering((True, True, False))

    f_conc = lambda x,y,z: zalesakDisk((x,y,z), 1.)
    if velocity_type=='shear':
      f_vel  = lambda x,y,z: rotationalVelocityField((x,y,z), (0,0,.1), (0.5, 0.5, 0.), True)
    else:
      f_vel  = lambda x,y,z: rotationalVelocityField((x,y,z), (0,0,1), (0.5, 0.5, 0.), False)

    convectionLevelset(dicttoinfo(params), ld, dict(conc_profile = f_conc, velocity_profile = f_vel))

    def levelsetshow(a):
      #t = np.max(np.abs(a))
      pyplot.imshow(np.asarray(a), origin='bottom', cmap = cm_ls, vmin=-0.1, vmax=0.1, extent=ld.worldBox[:4])

    rc = matplotlib.rc
    rc('font', size = 8.)
    rc('axes', titlesize = 10., labelsize = 10.)
    rc('pdf', compression = 6, fonttype = 42)
    rc('figure', **{'subplot.wspace' : 0.2,
                    'subplot.hspace' : 0.25})

    results = Results(name)
    pyplot.suptitle(name)
    pyplot.subplot(221, xlabel='x', ylabel='y', title='t=0')
    levelsetshow(results[0]['ls'][:,:,0])
    pyplot.colorbar()

    pyplot.subplot(222, title=r't=2$\pi$')
    levelsetshow(results[-1]['ls'][:,:,0])
    pyplot.colorbar()

#    t, m = zip(*[ (g.attrs['time'], g.attrs['mass']) for g in results ])
#    pyplot.subplot(223, title='mass loss', xlabel='t', ylabel='m', ylim=(0, max(m)*1.1))
#    pyplot.plot(t, m, 'x-k')
    t, gnorm = results.f['gradientmeasure/time'], results.f['gradientmeasure/gradnorm']
    pyplot.subplot(223, title=r'$||\nabla d| -1|_\inf$', xlabel='t', ylabel='', xlim=(-0.1*params['tend'],params['tend']*1.1))
    pyplot.plot(t, gnorm, '-k')

    pyplot.subplot(224, title=r'contours t=0 .. 2$\pi$', xlabel='x', ylabel='y')
    for i, g in enumerate(results[::2]):
      col = 'k' if i==0 else 'r'
      pyplot.contour(np.asarray(g['ls'][:,:,0]), [0.], extent=ld.worldBox[:4], colors = col)

    pyplot.savefig(name+'.pdf', dpi=180., papertype='a4')
    pyplot.show()



##############################################
## advection diffusion reaction scheme analysis
##############################################

def circle(x,y,z, center, radius, width):
  r = math.sqrt((x-center[0])**2 + (y-center[1])**2 + (z-center[2])**2)
  return min(1, max(0, (radius-r)/width))


def convectionDiffusionReaction(name, params, ld, **kwargs):
  params['fn_out'] = name
  libkrebs.convectionDiffusionReaction(dicttoinfo(params),ld, kwargs)


def runAdvectionDiffusionReactionAnalysis():
  params = dict(stepper = 'vsimexbdf2',
                tend = 1.,
                out_intervall = 0.1)

  ld = krebsutils.LatticeDataQuad3d((0, 49, 0, 49, 0, 0), 1./50.)
  ld.SetCellCentering((True, True, False))

  #funcs have the arguments (x,y,z,t,conc)
  velocity = lambda x,y,z,t,c: (0.5, 0, 0)
  src_impl_clinear = lambda x,y,z,t,c: -0.01
  src_impl_cconst = lambda x,y,z,t,c: 0.
  src_expl = lambda x,y,z,t,c: 0.
  kdiff = lambda x,y,z,t,c: 0.0
  conc_profile = lambda x,y,z,t,c,: circle(x,y,z, (0.2, 0.5, 0), 0.2, 0.01)

  name = 'cdrtest'
  convectionDiffusionReaction(name, params, ld,
                              conc_profile = conc_profile,
                              velocity = velocity,
                              src_impl_clinear = src_impl_clinear,
                              src_impl_cconst = src_impl_cconst,
                              src_expl = src_expl,
                              kdiff = kdiff)

  results = Results(name)

  mm = (1000000, -100000)
  for i, g in enumerate(results):
    a = np.asarray(g['conc'])
    mm = min(mm[0], a.min()), max(mm[1], a.max())
  print mm

  for i, g in enumerate(results):
    a = g['conc'][:,:,0]
    pyplot.imshow(np.asarray(a), vmin=mm[0], vmax=mm[1], extent=worldBox[:4])
    pyplot.show()


##############################################
##############################################

def testCurvature():
  ld = krebsutils.LatticeDataQuad3d((-10, 10, -10, 10, -10, 10), 1./20)
  ld.SetCellCentering((True, True, True))
  shape = ld.box
  print shape
  X, Y, Z = ld.scale*np.mgrid[shape[0,0]:shape[0,1]+1, shape[1,0]:shape[1,1]+1, shape[2,0]:shape[2,1]+1]
  radius = 0.25
  phi = np.sqrt(Y*Y+Z*Z+X*X) - radius
  curv, force = calcCurvature(ld, np.asarray(phi, dtype=np.float32), False, True)
  mask = np.logical_or(-ld.scale*2. > phi, phi > ld.scale*2)
  curv[mask] = 0.
  curv[0,...] = 0.
  curv[-1,...] = 0.
  curv[:,0,:] = 0
  curv[:,-1,:] = 0
  curv[:,:,0] = 0
  curv[:,:,-1] = 0

  def proj_(a, idx):
    if idx == 2:
      z = phi.shape[2]/2
      ext = ld.worldBox[[0,1,2,3]]
      return a[:,:,z], ext
      #---------------
    elif idx == 0:
      x = phi.shape[0]/2
      ext = ld.worldBox[[2,3,4,5]]
      return a[x,:,:], ext
      #---------------
    elif idx == 1:
      y = a.shape[1]/2
      ext = ld.worldBox[[1,2,4,5]]
      return a[:,y,:].transpose(), ext
  proj = lambda a: proj_(a,1)

  if 0:
    phi_x, ext = proj(phi)
    pyplot.subplot(211)
    pyplot.imshow(phi_x, extent=ext)
    pyplot.colorbar()
    pyplot.contour(phi_x, levels=[0.], extent=ext)
    curv_x, ext = proj(curv)
    pyplot.subplot(212)
    pyplot.imshow(curv_x, extent=ext)
    pyplot.colorbar()
    pyplot.contour(curv_x, levels=[2./radius], extent=ext)
    pyplot.show()

  f1, ext = proj_(force[0], 2)
  f2, ext = proj_(force[1], 2)
  pyplot.quiver(f1.transpose(), f2.transpose())
  pyplot.show()

##############################################
##############################################

def circle_distfunc(x,y,z, radius, center):
  r = math.sqrt((x-center[0])**2 + (y-center[1])**2 + (z-center[2])**2)
  return radius-r


def runSTFTest():
    params = dict(
      tend = 0.01,
      out_intervall = 0.0001,
      stepper = 'impeuler',
      fn_out = 'sfttest',
      num_threads = 2,
      kdiff = 10.,
      kstf = 0.5,
    )
    ld = krebsutils.LatticeDataQuad3d((0, 49, 0, 49, 0, 0), 1./50.)
    ld.SetCellCentering((True, True, False))
    extent = ld.worldBox[:4]
    #f_conc = lambda x,y,z: zalesakDisk((x,y,z), 0.3, (0.5,0.5,0.))
    f_conc = lambda x,y,z: circle_distfunc(x,y,z, 0.3, (0.5, 0.5, 0.))
    if not os.path.isfile('sfttest.h5'):
      surfaceTensionForceTest(dicttoinfo(params), ld, dict(conc_profile = f_conc))

    rc = matplotlib.rc
    rc('font', size = 8.)
    rc('axes', titlesize = 10., labelsize = 10.)
    rc('figure', **{'subplot.wspace' : 0.15,
                    'subplot.hspace' : 0.15})
    rc('savefig', facecolor = 'white')
    dpi = 100.
    results = Results('sfttest')
    ldface = (krebsutils.read_lattice_data_from_hdf(results.f['face_ld_0']),
              krebsutils.read_lattice_data_from_hdf(results.f['face_ld_1']))
    extface = tuple(q.worldBox[:4] for q in ldface)
    data = collections.defaultdict(list)
    for idx, group in enumerate(results):
      def show(a, ext):
        pyplot.imshow(np.asarray(a[...,0]).transpose(), extent = ext, origin='bottom', interpolation='nearest')
      def contour():
        return pyplot.contour(np.asarray(group['ls'][...,0]).transpose(), levels=[0.], extent = extent, origin='lower')
      print "plotting %i" % idx
      if 1:
        pyplot.figure(figsize = (1200./dpi,1000./dpi))

#        pyplot.subplot(221)
#        show(group['ls'], extent)
#        pyplot.colorbar()
#        contour()
#        pyplot.subplot(222)
#        show(group['rho'], extent)
#        pyplot.colorbar()
#        contour()

        pyplot.subplot(221)
        pyplot.plot(np.linspace(extent[0],extent[1],50), group['ls'][:,25,0])
        pyplot.gca().set(xlabel = 'x', ylabel = r'$\phi$')
        pyplot.gca().grid(color=(0.5,)*3)

        pyplot.subplot(222)
        pyplot.plot(np.linspace(extent[0],extent[1],50), group['rho'][:,25,0])
        pyplot.gca().set(xlabel = 'x', ylabel = r'$\rho$')
        pyplot.gca().grid(color=(0.5,)*3)

        pyplot.subplot(223)
        show(group['force_0'],extface[0])
        pyplot.colorbar()
        contour()
        pyplot.subplot(224)
        show(group['force_1'],extface[1])
        pyplot.colorbar()
        contour()
        pyplot.savefig("sfttest%04i.png" % idx, dpi = dpi)

      data['time'].append(group.attrs['time'])
      data['mass'].append(group.attrs['mass'])
      data['area'].append(group.attrs['area'])

    pyplot.figure(figsize = (1200./dpi,600./dpi))
    pyplot.subplot(121)
    pyplot.plot(data['time'], data['mass'], marker = 'x')
    pyplot.gca().set(xlabel = 't', ylabel = 'mass')
    pyplot.subplot(122)
    pyplot.plot(data['time'], data['area'], marker = '+')
    pyplot.gca().set(xlabel = 't', ylabel = 'area')
    pyplot.savefig('vstime.png')


def runSTFStabilityTest():
    '''
      this function runs a series of pde solutions with different time steps in order to find the proper constant
      for stability. A scaling with a power of deltax is implied, but the prefactor must be determined.
      When the method is unstable the solution explodes. So the final total mass in the system is plotted in order to view this.
    '''
    def run(deltax, dtscaling, filename):
      dim = 2
      kdiff = deltax**2  # cancel out deltax**2 in the diffusion CFL condition -> max_dt_diffusion = 1/dim
      kstf = 3000. # make it large so the timestep restriction for the stf dominates
      dt_max_diff = 0.8*0.5/kdiff*(deltax**2)/(dim);
      dt = dtscaling/(kdiff*dim*kstf)*(deltax**3) # the assumed scaling with various factors. dtscaling is determined experimentally.
      params = dict(
        dt = dt,
        tend = dt*20,
        out_intervall = dt * 10,
        stepper = 'euler',
        fn_out = filename,
        num_threads = 2,
        kdiff = kdiff,
        kstf = kstf,
        levelset_reinit = False,
      )
      print 'running with', dt, dt_max_diff
      assert dt < dt_max_diff
      # initial condition is a spherical region.
      pprint.pprint(params)
      size = (50,50,20)
      ld = krebsutils.LatticeDataQuad3d((0, size[0]-1, 0, size[1]-1, 0, size[2]-1), deltax)
      ld.SetCellCentering((True, True, False))
      f_conc = lambda x,y,z: circle_distfunc(x,y,z, size[0]*0.3*deltax, (size[0]*0.5*deltax, size[1]*0.5*deltax, size[2]*0.5*deltax))
      surfaceTensionForceTest(dicttoinfo(params), ld, dict(conc_profile = f_conc))

    def plot(filename):
      rc = matplotlib.rc
      rc('font', size = 8.)
      rc('axes', titlesize = 10., labelsize = 10.)
      rc('figure', **{'subplot.wspace' : 0.15,
                      'subplot.hspace' : 0.15})
      rc('savefig', facecolor = 'white')
      dpi = 100.
      results = Results(filename)
      ld = krebsutils.read_lattice_data_from_hdf(results.f['field_ld'])
      size = ld.shape
      extent = ld.worldBox[:4]
      ldface = (krebsutils.read_lattice_data_from_hdf(results.f['face_ld_0']),
                krebsutils.read_lattice_data_from_hdf(results.f['face_ld_1']))
      extface = tuple(q.worldBox[:4] for q in ldface)
      data = collections.defaultdict(list)
      for idx, group in enumerate(results):
        def show(a, ext):
          pyplot.imshow(np.asarray(a[...,size[2]/2]).transpose(), extent = ext, origin='bottom', interpolation='nearest')
        def contour():
          return pyplot.contour(np.asarray(group['ls'][...,size[2]/2]).transpose(), levels=[0.], extent = extent, origin='lower')
        print "plotting %i" % idx
        if 1:
          pyplot.figure(figsize = (1200./dpi,1000./dpi))

          pyplot.subplot(221)
          show(group['ls'], extent)
          pyplot.colorbar()
          contour()
          pyplot.subplot(222)
          show(group['rho'], extent)
          pyplot.colorbar()
          contour()

          pyplot.subplot(223)
          show(group['force_0'],extface[0])
          pyplot.colorbar()
          contour()
          pyplot.subplot(224)
          show(group['force_1'],extface[1])
          pyplot.colorbar()
          contour()
          pyplot.savefig("%s_%i.png" % (filename, idx), dpi = dpi)

        data['time'].append(group.attrs['time'])
        data['mass'].append(group.attrs['mass'])
        data['area'].append(group.attrs['area'])

    def determine_final_total_mass(filename):
      results = Results(filename)
      if results.f.attrs['hard_fail']: # it diverged over a certain threshold
        total_mass = 1.e6
      else:
        total_mass = np.sum(np.abs(np.asarray(results[-1]['rho'])))
      print 'analyze result:', total_mass
      return total_mass


    data = []
    deltax_list = [ 100., 10., 1., 0.1, 0.01 ]
    dtscaling_list = np.arange(1., 13., 1.)
    for deltax in deltax_list:
      for dtscaling in dtscaling_list:
        filename = "stfstability_%f_%f" % (deltax, dtscaling)
        if not os.path.isfile(filename+'.h5'):
          run(deltax, dtscaling, filename)
          #plot(filename)
        data.append((deltax, dtscaling, determine_final_total_mass(filename)))
    print data
    data = np.asarray(data)

    # scatter plot of total mass vs. parameters
    pyplot.close()
    pyplot.scatter(data[:,0], data[:,1], s = 200, c = np.log10(data[:,2]), cmap = matplotlib.cm.jet)
    pyplot.gca().set(xscale = 'log')
    pyplot.gca().set(xlabel = 'deltax', ylabel = 'dtscaling', title = '$log_{10}$(total mass)')
    pyplot.colorbar()
    pyplot.savefig("stability_analysis.png")
    # results indicate that the stability threshold is
    # dtscaling = 4, hence
    # dt <= 4./(kdiff*dim*kstf)*(deltax**3)


def check_surface_metrics():
  deltax = 1./20.
  size = 20.
  radius = 0.3
  if 0:
    field_func = lambda x,y,z: -circle_distfunc(x,y,z, radius, (0.5, 0.5, 0))
    x0,x1,y0,y1 = 0, size, 0, size
    X, Y = deltax*np.mgrid[x0:x1+1,y0:y1+1]
    data = np.vectorize(field_func)(X, Y, 0)

    delta = krebsutils.smooth_delta(data, deltax)

    line_estimate = np.sum(delta)*deltax**2
    actual_circumference = radius*math.pi*2

    pyplot.imshow(delta)
    pyplot.colorbar()
    pyplot.show()
  else:
    field_func = lambda x,y,z: -circle_distfunc(x,y,z, radius, (0.5, 0.5, 0.5))
    x0,x1,y0,y1,z0,z1 = 0, size, 0, size, 0, size
    X, Y, Z = deltax*np.mgrid[x0:x1+1,y0:y1+1,z0:z1+1]
    data = np.vectorize(field_func)(X, Y, Z)

    delta = krebsutils.smooth_delta(data, deltax)

    line_estimate = np.sum(delta)*deltax**3
    actual_circumference = (radius**2)*math.pi*4

    print 'estimate: ', line_estimate
    print 'actual: ', actual_circumference


def runTestEllipticSolver():
  try:
    make_image = sys.argv[1] == '--img'
  except:
    make_image = False
  if make_image:
    size = (128, 128, 1)
    krebsutils.set_num_threads(4)
    time, mem = EST_testOperatorApplication(size, True)
    time, mem = EST_testOperatorSolve(size, True)
  else:
    if sys.flags.debug:
      size = (128, 128, 32)
      krebsutils.set_num_threads(4)
      time, mem = EST_testOperatorSolve(size, False)
    else:
      size = (128, 128, 128)
      times = collections.defaultdict(list)
      num_repeat = 5
      for num_threads in [1,2,3,4]:
        krebsutils.set_num_threads(num_threads)
        print 'np %i' % num_threads,
        for i in range(num_repeat):
          print '.',
          time, mem = EST_testOperatorSolve(size, False)
          times[num_threads].append(time)
        print ''
      for k, v in times.iteritems():
        times[k] = np.average(v)
      base = times[1]
      for k in sorted(times.keys()):
        v = times[k]
        print 'np = %i, time: %f, speedup: %f' % (k, v, base/v)



if __name__ == "__main__":
  #runTestEllipticSolver()
  check_surface_metrics()
#  runSTFTest()
  #runSTFStabilityTest()
#  testCurvature()
#  runStepperConvergenceTests()
#  runAdvectionSchemeAnalysis()
  #runLevelsetRedistancingScheme()
  #runLevelsetAdvectionScheme()
#  runAdvectionDiffusionReactionAnalysis()
