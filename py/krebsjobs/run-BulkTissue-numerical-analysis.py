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

import os, sys
from os.path import join, basename, dirname, splitext
import qsub
from dicttoinfo import dicttoinfo, Vec
from copy import deepcopy, copy
import glob
#from pprint import pprint
import myutils
import mpl_utils
import krebs.plotBulkTissue1d
import krebs.plotBulkTissue2d

exe = '/home/mwelter/tumorcode/buildopt/tum-bulktissue-no-vessels'
#exe = '/home/mwelter/tumorcode/build/src/programs/tum-bulktissue-no-vessels'
dstdir = ''

qsub.parse_args(sys.argv)

def mkdir(outfn):
  p = dirname(outfn)
  if p and not os.path.isdir(p):
    print 'creating', p
    os.makedirs(p)


def run(c):
  c = deepcopy(c) # so it can be changed for the simulation without disturbing custom data
  name = c.pop('name')
  fn = join(dstdir, name)
  c['fn_out'] = fn
  if not os.path.isfile(name+'.h5'):
    mkdir(c["fn_out"])
    res = os.system("ulimit -c unlimited; %s -pa '%s'" % (exe, dicttoinfo(c)))
    if res > 0: sys.exit(res)
  #  qsub.submit(qsub.exe([exe, '-pa', "'%s'" % dicttoinfo(c)]),
  #              name = name,
  #              days = 2,
  #              num_cpus = 2,
  #              change_cwd = True)

def plot(name, dim=1):
  import pickle
  checks = {}
  try:
    f = open('filelist.pickle', 'r')
    checks = pickle.load(f)
    f.close()
  except:
    pass
  fn = join(dstdir,name)
  try:
    t = str(os.path.getmtime(fn+'.h5'))
  except OSError:
    return
  if not fn in checks or t <> checks[fn] or not os.path.isfile(fn+'.pdf'):
    if dim == 1:
      krebs.plotBulkTissue1d.plotit(fn+'.h5')
    elif dim == 2:
      krebs.plotBulkTissue2d.doit(fn+'.h5')
    checks[fn] = t
    f = open('filelist.pickle', 'w')
    pickle.dump(checks, f)
    f.close()


basic = dict(
  fn_out = 'unset_filename',
  lattice_size = Vec((200, 1, 1)),
  lattice_scale = 30.,
  tend = 1000.,
  out_intervall = 100,
  save_hdf = True,
  save_images = False,
  tumor = dict(
    tumor_size = 1000.,
    time_prol = 24.,
    time_death = 240.,
    time_necrosis = 4800000.,
    use_necrotic_regions = True,
    ncells_norm = 0.4,
    ncells_tumor = 0.5,
    ncells_sigma = 0.1,
    ncells_ecm = 0.2,
    cell_mobility_ecmstar = 0.2,
    cell_mobility_ecmstarstar = 0.4,
    cell_mobility = 3600.,
    cell_mobility_tumor = 3600.,
    write_levelset_function = 1,
    write_face_velocities = 1,
    timestep_factor = 1.,
    wrong_model = 0,
    source_model_version = 3,
    mixture_stepper = 'impeuler'
  ),
  #o2_range = 100.,
  o2_rel_tumor_source_density = 1.,
  o2_level_normal = 0.8,
  use_o2_source_decay = False,
  num_threads = 3,
)
  
def sp(c, params=dict(), tumorparams=dict()):
  c.update(params)
  c['tumor'].update(tumorparams)

def multi_sp(runs, postfix, params=dict(), tumorparams=dict(), return_copy=True):
  d = dict()
  for k, c in runs.iteritems():
    c = deepcopy(c)
    c['name'] = c['name']+postfix
    sp(c, params, tumorparams)
    d[c['name']] = c
  if return_copy:
    return d
  else:
    runs.update(d)
  
  
simple = deepcopy(basic)
simple['name'] = 'simple'
simple['tumor']['cell_mobility'] = 5000.
simple['tumor']['cell_mobility_tumor'] = 5000
simple['tumor']['ncells_tumor'] = 0.6

obst = deepcopy(basic)
obst['name'] = 'obst'
obst['tumor']['sigma_model_version'] = 2
obst['test_obstacle'] = 1
obst['tend'] = 1600

obst2 = deepcopy(obst)
obst2['name'] = 'obst2'
obst2['tumor']['cell_mobility_ecmstar'] = -1
obst2['tumor']['cell_mobility_ecmstarstar'] = -1


steady = deepcopy(basic)
steady['name'] = 'steady'
steady['tend'] = 500.
steady['tumor']['ncells_norm'] = 0.4
steady['tumor']['ncells_tumor'] = 0.4
steady['tumor']['wrong_model'] = True
steady['test_obstacle'] = 1

lowoxy = deepcopy(basic)
lowoxy['name'] = 'lowoxy'
lowoxy['o2_level_normal'] = 0.9
lowoxy['use_o2_source_decay'] = True
lowoxy['o2_rel_tumor_source_density'] = 0.
lowoxy['tumor']['o2_necro_threshold'] = 0.03
lowoxy['tumor']['o2_prol_threshold'] = 0.3  

necro = deepcopy(basic)
necro['name'] = 'necro'
necro['tumor']['time_necrosis'] = 48
necro['tumor']['o2_necro_threshold'] = 0.03
necro['tumor']['o2_prol_threshold'] = 0.3
necro['use_o2_source_decay'] = True
necro['o2_rel_tumor_source_density'] = 0.

necro2 = deepcopy(necro)
necro2['name'] = 'necro2'
necro2['tumor']['time_death_tumor'] = 24000.
necro2['tumor']['cell_mobility'] = 5000
necro2['tumor']['cell_mobility_tumor'] = 5000

necro3 = deepcopy(necro2)
necro3['name'] = 'necro3'
necro3['tumor']['ncells_norm'] = 0.4
necro3['tumor']['ncells_tumor'] = 0.6
necro3['o2_range_tumor'] = 50  
necro3['o2_source_decay_time'] = 48
 

noisy  = deepcopy(necro3)
noisy['name'] = 'noisy'
noisy['tumor']['sigma_model_version'] = 2
noisy['tumor']['ecm_noise_density'] = 0.01
noisy['tumor']['ecm_noise_std'] = 1.
 
#necro4 = deepcopy(necro3)
#necro4['name'] = 'necro4'
#necro4['o2_source_decay_time'] = 48*2 

necro5 = deepcopy(necro2)
necro5['name'] = 'necro5'
necro5['tumor']['time_death_tumor'] = 240.


growtharrest = deepcopy(basic)
growtharrest['name'] = 'growtharrest'
growtharrest['tend'] = 1000.
growtharrest['o2_level_normal'] = 0.7
growtharrest['o2_range_normal'] = 200.
growtharrest['o2_range_tumor'] = 50.
growtharrest['o2_range_necro'] = 200.
growtharrest['tumor']['tumor_size'] = 100.
growtharrest['tumor']['o2_necro_threshold'] = 0.04
growtharrest['tumor']['o2_prol_threshold'] = 0.4
growtharrest['tumor']['time_death_tumor'] = 24000000.
growtharrest['o2_source_decay_time'] = 48
growtharrest['use_o2_source_decay'] = False
growtharrest['o2_rel_tumor_source_density'] = 0.

runs = [lowoxy, simple, steady, necro, necro2, necro3, obst, noisy, necro5, obst2, growtharrest ]
runs = dict((r['name'], r) for r in runs)

runs_imex = multi_sp(runs, 'imex', tumorparams=dict(mixture_stepper='vsimexbdf2', timestep_factor=0.5))

runs_nonecroregions = multi_sp(runs, '_noregions', tumorparams=dict(use_necrotic_regions=False))

runs.update(runs_imex)
runs.update(runs_nonecroregions)

for r in runs.itervalues():
  run(r)
  plot(r['name'])

necro3_2d_rect = deepcopy(necro3)
necro3_2d_rect['name'] += '_rect'
necro3_2d_rect['tumor']['init_shape'] = 'rectangle'

runs_2d = [ necro2, necro3, necro3_2d_rect, runs_nonecroregions['necro3_noregions'], runs_nonecroregions['necro2_noregions'] ]
runs_2d = dict((r['name'], r) for r in runs_2d)
runs_2d = multi_sp(runs_2d, '_2d', params=dict(lattice_size=Vec((200, 200, 1)), num_threads=4),tumorparams=dict(mixture_stepper='vsimexbdf2', timestep_factor=0.5))


for r in runs_2d.itervalues():
  run(r)
  plot(r['name'], dim=2)




def run_numerical_tests():
  
  convergence_step_widths = [ 0.01, 0.02, 0.04, 0.08, 0.16, 0.32, 0.64, 1. ]
  stepper_methods = [ 'euler', 'impeuler', 'vsimexbdf2', 'rk3', 'rk4' ]
  gridsizes = [ 30, 5, 1 ]
  
  def run_different_steppers(name, cc):
    for stepper in stepper_methods:
      c = deepcopy(cc)
      c['stepper'] = stepper
      if stepper == 'vsimexbdf2':
        c['tumor']['timestep_factor'] = 0.2
      run(name+'-'+stepper, c)
  
  
  def run_3d(name, cc, size, stepper):
    c = deepcopy(cc)
    c['tumor']['timestep_factor'] = 0.1 if (stepper == 'vsimexbdf2') else 1.
    c['stepper'] = stepper
    c['lattice_size'] = Vec((size,)*3)
    c['tumor']['interface_width'] = 30.
    c['tend'] = 20
    c['out_intervall'] = 20
    c['tumor']['write_face_velocities'] = False
    run(name+"-s%03i-%s" % (size, stepper), c)
  
  
  def run_3d_performance_tests():
    for l in [ 40, 70, 100, 120, 150 ]:
      run_3d('prez3d', cfg, l, 'vsimexbdf2')  
      run_3d('prez3d', cfg, l, 'impeuler')  
  
  
  def run_lattice_and_stepwidth(name, cc, s, h, stepper):
    c = deepcopy(cc)
    c['tend'] = 1000.
    c['lattice_size'] = Vec((int(100*30/s), 1, 1))
    c['lattice_scale'] = s
    c['stepper'] = stepper
    c['tumor']['timestep_factor'] = h
    c['tumor']['elliptic_solver_params'] = { 'use_multigrid' : True }
    run(name+"-conv-s%02i-h%0.3f-%s" % (s,h,stepper), c)
  
  
  def run_convergence_test(name, cc):
    cc = deepcopy(cc)
    cc['o2_rel_tumor_source_density'] = 1.
    name = name+'-maxo2'
  
    for lw in [ 5, 10, 15, 30, 60 ]:
      run_lattice_and_stepwidth(name, cc, lw, 1, 'impeuler')
    run_lattice_and_stepwidth(name, cc, 5, 1, 'rk4')
    
  #  cc['tumor']['interface_width'] = 60.
  #  for lw in [ 5, 10, 15, 30, 60 ]:
  #    run_lattice_and_stepwidth(name+'-reg30', cc, lw, 1, 'impeuler')
  #  run_lattice_and_stepwidth(name+'-reg30', cc, 5, 1, 'rk4')
  #  
  #  cc['tumor']['interface_width'] = 120.
  #  for lw in [ 5, 10, 15, 30, 60 ]:
  #    run_lattice_and_stepwidth(name+'-reg60', cc, lw, 1, 'impeuler')
  #  run_lattice_and_stepwidth(name+'-reg60', cc, 5, 1, 'rk4')
  
  
  
  if 0:
    # test steppers
    #run_different_steppers('prez1d', cfg)
    #run_3d_performance_tests()
    # check how the solution depends on the lattice spacing and time step length
    run_convergence_test('prez1d', cfg)
  
  
  
  if 1:
  
    rc = matplotlib.rc
    rc('font', size = 12)
    rc('axes', titlesize = 12., labelsize = 12.) 
    rc('legend', fontsize = 9)
    plt.subplots_adjust(wspace=0.5)
    
    if 0: #  plot individuals
      for fn in ['prez1d-%s.h5' % n for n in stepper_methods]:
        print 'plotting', fn
        pdf = mpl_utils.PdfWriter(splitext(fn)[0] + ".pdf")
        fd = Filedata1d(fn)
        for g in fd.groupnames:
          simplePlot(pdf, fd, g, fd.f['parameters'].attrs['stepper'])
        pdf.close()
      
    if 0: # plot comparison of different methods
      pdf = mpl_utils.PdfWriter("compare1.pdf")
      fdnames = ['euler', 'vsimexbdf2', 'rk3']
      fds = [ Filedata1d('prez1d-%s.h5' % n) for n in fdnames ]
      for g in fds[0].groupnames[-1:]:
        print 'plotting', g
        comparePlot(pdf, fds, fdnames, g, ['ptc', 'ls', 'vel_0', 'press'], 'Methods')
      pdf.close()    
  
    if 1: # plot comparison of different time steps
      pdf = mpl_utils.PdfWriter("compare3.pdf")    
      
      fd_ref = Filedata1d('prez1d-maxo2-conv-s05-h1.000-rk4.h5')
      grp_name = fd_ref.groupnames[-1]
      
      def compare(name, title, configs):
        title_format = r'$l$=%i, $h=%.1f h_{euler}$, %s'
        fds, fdtitles = [], []
        for s, h, m in configs:
          fn ="%s-conv-s%02i-h%0.3f-%s.h5" % (name, s,h,m)
          print 'plotting', fn
          fds.append(Filedata1d(fn))
          fdtitles.append(title_format % (s,h,m))
        iwidth = float(fds[0].f['parameters']['tumor'].attrs.get('interface_width',-1))
        iwidth = r"$4 \cdot l$" if iwidth < 0 else (r"$%.1f \mu m$" % iwidth)
        fig = comparePlot(fds, fdtitles, grp_name, ['ptc', 'ls', 'vel_0', 'press', 'conc'], 'Stepwidths'+title+(' $\delta=$%s  ' % iwidth))
        pdf.savefig(fig)
      
  #    compare('prez1d', '', [
  #      (1,1,'vsimexbdf2'),
  #      (5,1,'vsimexbdf2'),
  #      (30, 1, 'vsimexbdf2')
  #      ])
  #
  #    compare('prez1d', '', [
  #      (1,1.,'vsimexbdf2'),
  #      (10,1.,'rk4'),
  #      (30,1., 'rk4')])
  
      stepper = 'impeuler'
      stepw = 1.
      compare('prez1d-maxo2', '  $max_{o2}$  ', [
        (5,stepw,stepper),
        (10,stepw,stepper),
        (15,stepw,stepper),
        (30,stepw,stepper),
        (60,stepw,stepper),
        (5,1,'rk4'),
        ])
    
      compare('prez1d-maxo2-reg30', r'  $max_{o2}$  ', [
        (5,stepw,stepper),
        (10,stepw,stepper),
        (15,stepw,stepper),
        (30,stepw,stepper),
        (60,stepw,stepper),
        (5,1,'rk4'),
      ])
  
      compare('prez1d-maxo2-reg60', r'  $max_{o2}$  ', [
        (5,stepw,stepper),
        (10,stepw,stepper),
        (15,stepw,stepper),
        (30,stepw,stepper),
        (60,stepw,stepper),
        (5,1,'rk4'),
      ])
  
      pdf.close()
  
    if 0: # memory and runtime
      pdf = mpl_utils.PdfWriter("performance1.pdf")
      data = []
      for fn in  glob.glob('prez3d-*.h5'):
        print fn
        fd = Filedata1d(fn)
        g = fd.f[fd.groupnames[-1]]
        s = Struct(dict(
          (k, g.attrs[k]) for k in ['mem_vsize', 'mem_rss', 'real_time', 'time'],
        ))
        s['stepper'] = fd.f.attrs['stepper']
        s['size'] = fd.f['field_ld'].attrs['SIZEX']
        data.append(s)
      data2 = collections.defaultdict(list)
      for d in data:
        data2[(d['stepper'])].append(d)
      for k, l in data2.iteritems():
        l.sort(key = lambda d: d.size)
        l = myutils.zipListOfDicts(l, numpy_output=True)
        l['mem_rss'] /= (1024.**2)
        l['mem_vsize'] /= (1024.**2)
        l['real_time'] /= 1000.
        l['dofs'] = l['size']**3
        data2[k] = l
        
      fig = plt.figure()
      plt.xlabel('lateral size')
      plt.ylabel('mem [Mb]')
      for stepper, l in data2.iteritems():
        plt.plot(l['size'], l['mem_vsize'], label='%s, rss' % (stepper))
        plt.plot(l['size'], l['mem_rss'], label='%s, vsize' % (stepper))
      plt.legend()
      pdf.savefig(fig)
  
      fig = plt.figure()
      plt.xlabel('lateral size')
      plt.ylabel('mem / sites [bytes]')
      for stepper, l in data2.iteritems():
        plt.plot(l['size'], l['mem_vsize']*(1024**2)/l['dofs'], label='%s, rss' % (stepper))
        plt.plot(l['size'], l['mem_rss']*(1024**2)/l['dofs'], label='%s, vsize' % (stepper))
      plt.legend()
      pdf.savefig(fig)
  
      fig = plt.figure()
      plt.xlabel('lateral size')
      plt.ylabel('run time / sim time [s]')
      for stepper, l in data2.iteritems():
        plt.plot(l['size'], l['real_time']/l['time'], label='%s' % (stepper))
      plt.legend()
      pdf.savefig(fig)
  
      fig = plt.figure()
      plt.xlabel('lateral size')
      plt.ylabel('run time / (sim time $\cdot$ sites) [$\mu$s]')
      for stepper, l in data2.iteritems():
        plt.plot(l['size'], l['real_time']*(1000**2)/l['dofs']/l['time'], label='%s' % (stepper))
      plt.legend()
      pdf.savefig(fig)
      
      pdf.close()
      
  #    ind = np.arange(N)
  #    width = 1.
  #    data = zip_dicts(data, numpy_output=True)
  #    data['mem_rss'] /= 1024**2
  #    data['mem_vsize'] /= 1024**2
  #    fig = plt.figure(figsize = (4,4))
  #    plt.bar(ind, data['mem_rss'], width/2., label='rss')
  #    plt.bar(ind+width/2., data['mem_vsize'], width/2., label='vsize')
  #    plt.title('Memory Consumption')
  #    plt.xticks(ind+width/2., stepper_methods)
  #    plt.legend()
  #    fig.autofmt_xdate()
  #    pdf.savefig(fig)
  #    plt.close()
  #    fig = plt.figure(figsize = (4,4))
  #    plt.bar(ind, data['real_time'], width)
  #    plt.title('Runtime [ms]')
  #    plt.xticks(ind+width/2., stepper_methods)
  #    fig.autofmt_xdate()
  #    pdf.savefig(fig)
  #    pdf.close()
