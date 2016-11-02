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

import os
import sys
from os.path import join
import qsub
from copy import copy, deepcopy

qsub.parse_args(sys.argv)



hostname = os.uname()[1]
if hostname == 'durga':
  home = '/localdisk/mwelter/'
  resultdir = '/localdisk/mwelter/phasefield_instability_jobs/results'
  cmd = join(home,'phasefield_instability_jobs','build2drel','phasefield_instability')
  cmd3d = join(home,'phasefield_instability_jobs','build3drel','phasefield_instability')
elif hostname.startswith('sleipnir'):
  home = '/home/mwelter/'
  resultdir  = '/home/mwelter/storage1/instability'
  resultdir2 = '/home/mwelter/storage1/instability2'
  cmd = join(home,'tumorcode','buildopt','phasefield_instability')
else:
  home = '/home/mwelter/'
  resultdir  = '/home/mwelter/results/instability'
  cmd = join(home,'tumorcode','buildopt','phasefield_instability')
  
  #cmd3d = join(home,'tumorcode','deal.ii-tests','phasefield_instability','build3dopt','phasefield_instability')

relparams = dict(
  num_threads = 1,
  domain_size = 500,
  coarse_lattice_const = 100,
  finest_lattice_const = 0.5,
  out_intervall = 10,
  max_time = 1000,
  triangulation_smoothing_level = 2, # 1 by default
  interface_width = 1.5, #1,
  initial_radius = 10,
  oxygen_diffusion_range = 5,
  oxygen_homogeneous_solution = 0.9,
  proliferation_coefficient = 0.1,
  source_term_type = 0,
  flow_basisfunction_degree = 1
)


def run(id, **additional_params):
  p = deepcopy(relparams)
  p.update(additional_params)
  p['out_fn'] = join(p.pop('dir',resultdir), id)
  mycmd = p.pop('cmd', cmd)
  days = p.pop('days', None)
  mem = p.pop('mem', None)
  qsub.qsub('job_'+id, cmd, p, num_cpus = p.get('num_threads',1), days = days, mem = mem)


def runbessa(name, **args):
  a = dict(
    dir=resultdir2,
    source_term_type=1,
    flow_basisfunction_degree=0,
    out_intervall = 100,
    max_time = 20000)
  run(name, **a)


if True:
  # must check gridding parameters and consequential solutions
  def thisrun(name, **kwargs):
    c = dict(
      domain_size = 200.,
      coarse_lattice_const = 100,
      initial_radius = 8,
      max_time = 1000,
      out_intervall = 10,
      tumor_mobility = 20.,
      surface_tension = 0.1,
      finest_lattice_const = 0.5,
      triangulation_smoothing_level= 1,
      num_threads = 2)
    c.update(kwargs)
    run('phasetum-test-'+name+"-st0.1", **c)
  thisrun('lc1', finest_lattice_const=1)
  thisrun('gsmooth0', triangulation_smoothing_level=0)
  thisrun('base', triangulation_smoothing_level=1)
  thisrun('gsmooth2', triangulation_smoothing_level=2)
  thisrun('lc2', finest_lattice_const=2)
  thisrun('iw0.33', interface_width=0.33)
  thisrun('iw1', interface_width=1)
  thisrun('iw3-lc1', interface_width=3, finest_lattice_const=1)
  thisrun('iw3', interface_width=3)
  thisrun('iw3-lc2', interface_width=3, finest_lattice_const=2)
  thisrun('iw9', interface_width=9.)
  


if False:
  # very long sim time, should give nice fingering
  for visc in [ 1., 10. ]:
    for sigma in [ 0., 0.01, 0.1, 1 ]:
        name = 'phasetum-bessa-vc%02i-sig%03i' % (visc, 100*sigma)
        runbessa(name, tumor_mobility = visc, surface_tension = sigma)

if False:
  for visc in [ 1, 2., 10, 20, 100 ]:
    for sigma, sigmastr in [ (-1., "0"), (0.01, "001") ]:
      name = 'phasetum-hele-shaw-%04ix-sigma-%s' % (visc * 10, sigmastr)
      print name
      run(name, tumor_mobility = visc, surface_tension = sigma)

if False:
  for visc in [ 1, 1.1, 1.5, 2., 10, 20, 100 ]:
    for ro2 in [ 1., 3., 5., 10 ]:
      name = 'phasetum-hele-shaw-%04ix-ro%02i' % (visc * 10, ro2)
      run(name, relparams, tumor_mobility = visc, oxygen_diffusion_range = ro2, surface_tension = -1, num_threads = 2, interface_width = 1.5, triangulation_smoothing_level = 2)

if False:
  for size in [ 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600 ]:
    run('phasetum-scaling-%i' % size,
	relparams, initial_radius = 5, out_intervall=20, max_time=100, domain_size = size)

if False:
  run('phasetum-hele-shaw-10x-sigma-01', relparams, interface_width = 1.5, triangulation_smoothing_level = 2, tumor_mobility = 10, surface_tension = 0.1, num_threads = 1)
  run('phasetum-hele-shaw-10x-sigma-001', relparams, interface_width = 1.5, triangulation_smoothing_level = 2, tumor_mobility = 10, surface_tension = 0.01, num_threads = 1)
  run('phasetum-hele-shaw-10x-sigma-0001', relparams, interface_width = 1.5, triangulation_smoothing_level = 2, tumor_mobility = 10, surface_tension = 0.001, num_threads = 2)
  run('phasetum-hele-shaw-10x-sigma-1', relparams, interface_width = 1.5, triangulation_smoothing_level = 2, tumor_mobility = 10, surface_tension = 1)
  run('phasetum-hele-shaw-10x-sigma-10', relparams, interface_width = 1.5, triangulation_smoothing_level = 2, tumor_mobility = 10, surface_tension = 10)

if False:
  run('phasetum-hale-shaw-100x_medres'      , relparams, finest_lattice_const = 1, triangulation_smoothing_level = 2, tumor_mobility = 100)
  run('phasetum-hale-shaw-100x_oxr20_medres', relparams, finest_lattice_const = 1, triangulation_smoothing_level = 2, tumor_mobility = 100, oxygen_diffusion_range = 20)
  run('phasetum-hale-shaw-10x_oxr20_medres' , relparams, finest_lattice_const = 1, triangulation_smoothing_level = 1, tumor_mobility = 10, oxygen_diffusion_range = 20)

if False:
  run('phasetum-hale-shaw-10x_medres'       , relparams, finest_lattice_const = 1, triangulation_smoothing_level = 1, tumor_mobility = 10)
  run('phasetum-medres'                     , relparams, finest_lattice_const = 1)

if False:
  run('phasetum-anisotropic099lowres', relparams, isotropy = 0.99, finest_lattice_const = 1)
  run('phasetum-anisotropicsphere09lowres', relparams, isotropy = 0.9, sphere_domain = 1, finest_lattice_const = 1)
  run('phasetum-anisotropicsphere099lowres', relparams, isotropy = 0.99, sphere_domain = 1, finest_lattice_const = 1)

if False:
  run('phasetum-sinks10lowres', relparams, pressure_diffusion_range = 10, finest_lattice_const = 2)
  run('phasetum-sinks20lowres', relparams, pressure_diffusion_range = 20, finest_lattice_const = 2)

if False:
  run('phasetum-anisotropic09lowres', relparams, isotropy = 0.9, finest_lattice_const = 2)
  run('phasetum-anisotropic08lowres', relparams, isotropy = 0.8, finest_lattice_const = 2)
  run('phasetum-rt1-oxr5', relparams, interface_relaxation_time = 1.)

if False:
  run('phasetum-spheredomlowres', relparams, sphere_domain = 1, finest_lattice_const = 2)

if False:
  run('phasetum-phasemob', relparams, phase_dependent_permeability = 1)
  run('phasetum-nosmooth', relparams, triangulation_smoothing_level = 0)
  run('phasetum-spheredom', relparams, sphere_domain = 1)

if False:
  # base case runs
  parset = [
    (0.2, 1),
    (0.2, 3),
    (0.2, 5),
    (0.2, 10),
    (0.2, 20),
    (0.02, 5),
  ]
  for (rt, oxr) in parset:
    run("phasetum-rt%0.2f-oxr%i" % (rt,oxr), relparams, 
      interface_relaxation_time = rt,
      oxygen_diffusion_range = oxr
    )

if 0:
  # parameter group injection in existing file
  def run(fn, baseparams, **override_params):
    fn = fn[len('phasetum-'):]+'.h5'
    p = baseparams.copy()
    p.update(override_params)
    f = h5py.File(fn, 'r+') # read/write file must exist
    if 'parameters' in f:
      print 'skipped', fn
      return
    print 'processing', fn
    g = f.require_group('parameters')
    for k, v in p.iteritems():
      g.attrs[k] = v
    f.close()

########################################################################
###  3d mode                                                          ##
########################################################################

params3d = dict(
  num_threads = 8,
  coarse_lattice_const = 125,
  finest_lattice_const = 0.5,
  domain_size = 250,
  interface_width = 1.5,
  triangulation_smoothing_level = 1,
  
  out_intervall = 0.5,
  max_time = 1000,
  
  oxygen_diffusion_range = 5,
  oxygen_homogeneous_solution = 0.9,
  initial_radius = 10,
  proliferation_coefficient = 0.1,
  source_term_type = 0
)

if False:
  run('basic-3d', params3d, cmd = cmd3d, days = 32, mem='10gb')
