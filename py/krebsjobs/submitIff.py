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
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'..'))
import identifycluster
if identifycluster.name == 'snowden':
  import matplotlib
  matplotlib.use('Agg') 

import qsub
from copy import deepcopy
import dicttoinfo
import numpy as np
import h5py

import krebsjobs.parameters.parameterSetsIff as paramIff

import myutils

def centerslice(a):
  return a[:,:,a.shape[2]/2]

class Measure(object):
  @staticmethod
  def obtain_dataset(f, name, shape):
    try:
      return f[name]
    except KeyError:
      return f.create_dataset(name, shape, dtype = np.float32)

  def __init__(self, outfilename, params):
    self.lastt = 0.
    self.outfilename = outfilename
    exposure_thresholds_reference_conc_in = params.pop('exposure_thresholds_reference_conc_in', 1.)
    exposure_thresholds_reference_conc_ex = params.pop('exposure_thresholds_reference_conc_ex', 1.)
    baselevels = np.asarray([ 0.01, 0.1, 0.25, 0.5, 0.75, 1., 10. ])
    self.exposure_thresholds = (baselevels * exposure_thresholds_reference_conc_ex,
                                baselevels * exposure_thresholds_reference_conc_in)
    self.moviefilename = params.pop('moviefilename', None)
    self.movieintervall = params.pop('movieintervall', None)
    self.nextmovie_t = 0
    self.movie_out_number = 0


  def __call__(self, *args, **kwargs):
    (t, _), conc = args
    """should measure locally at very short time steps
     *AUC,
     *max conc
     *exposure time for c_max greater than
      + c_1 \approx 0
      + c_2 > 0
      + c_i ...
      + c_n < c_inject
     conc[0] is extracellular intrinsic
     conc[1] is intracellular intrinsic
    """
    # this code is for recording integral and other metrics at short intervalls
    dt = t - self.lastt
    self.lastt = t
    print("try to open: %s" % self.outfilename)
    with h5py.File(self.outfilename, 'r+') as f:
      shape = conc[0].shape
      g = f.require_group('measurements').require_group('drug_local_integral')
      for i, in_ex in enumerate(['ex','in']):
        ds = Measure.obtain_dataset(g, 'auc_'+in_ex, shape)
        ds[:] = np.asarray(ds) + dt * conc[i] # numerical integration as step function -> AUC metric
        ds = Measure.obtain_dataset(g, 'c_max_'+in_ex, shape)
        ds[:] = np.maximum(np.asarray(ds), conc[i]) # maximal conc. metric
      for i, in_ex in enumerate(['ex','in']):
        for j, exposure_conc in enumerate(self.exposure_thresholds[i]):
          ds = Measure.obtain_dataset(g, 'exposure_time_%02i_%s' % (j,in_ex), shape)
          ds[:] = np.asarray(ds) + dt * (conc[i] > exposure_conc)
          ds.attrs['exposure_conc'] = exposure_conc
          ds.attrs['exposure_conc_num'] = j
      g.attrs['tend'] = t
    # this code is for recording slices through the conc. distributions at short intervalls
    if self.moviefilename and self.movieintervall and (t - self.nextmovie_t >=  -0.1*dt):
      self.nextmovie_t += self.movieintervall
      with h5py.File(self.moviefilename, 'w' if self.movie_out_number==0 else 'r+') as f:
        g = f.require_group('frames')
        idx = g.attrs.get('num_frames', 0)
        g.attrs.modify('num_frames', idx+1)
        print 'write t = %f, idx = %i' % (t, idx)
        gt = g.create_group('%04i' % idx)
        gt.attrs['time'] = t
        gt.attrs['number'] = idx
        gt.create_dataset('conc_ex', data = centerslice(conc[0]), compression = 9, dtype = np.float32)
        gt.create_dataset('conc_in', data = centerslice(conc[1]), compression = 9, dtype = np.float32)
      self.movie_out_number += 1

def runs_on_client(name,config):
#  import krebs.iffsim
#  krebs.iffsim.run_iffsim(config)
  from krebs.iff import iff_cpp
  import krebsutils
  ''' now config is asumed to be imported from a py dict '''
  params = deepcopy(config)
  outfilename = params.pop('fn_out')
  measure_params = params.pop('ift_measure', dict())

  #set_num_treads is DEPRECATED we use the environment setting now
  #krebsutils.set_num_threads(params.pop('num_threads', 1))
  #krebsutils.run_iffsim(dicttoinfo.dicttoinfo(params), str(outfilename.encode('utf-8')), Measure(outfilename, measure_params))
  iff_cpp.run_iffsim(dicttoinfo.dicttoinfo(params), str(outfilename.encode('utf-8')), Measure(outfilename, measure_params))
  if __debug__:
    print('iff_cpp returned')


def run_simple(name, config):
  if 1:
    days = 10. if config['ift'] else 0.5
    qsub.submit(qsub.func(runs_on_client, name, config),
                name = 'job_iff'+name,
                num_cpus = 2,#not used on c++ side anymore
                days = days,
                change_cwd = True,)
  else:
    fn = splitext(config['fn_out'])[0]
    import json
    with open(fn+'.json','w') as f: f.write(json.dumps(config, indent=4))
#      open(fn+'.info','w').write(dicttoinfo(config))

def run(outdir, name, tumorfn, config, p=dict(), piff=dict(), pdrug=dict(), grp_pattern = ''):
  group = grp_pattern  
  dirs = set()
  for fn in filenames:
    with h5py.File(fn, 'r') as f:
      d = myutils.walkh5(f, grp_pattern)
      assert len(d), ' pattern "%s" not found in "%s"!' % (grp_pattern, fn)
      dirs =set.union(dirs, d)
  print 'and resolved groups therein: %s' % ','.join(dirs)  
  
  c = deepcopy(config)
  #c.update(p)
  #c['iff'].update(piff)
  #c['ift'].update(pdrug)
  if 'adaption' in group:
    setvesselfile(c, tumorfn, group, '', None)
  else:
    setvesselfile(c, tumorfn, group, 'vessels', 'tumor')

  stripped_name = splitext(basename(tumorfn))[0]
  stripped_name = stripped_name[len('tum-'):]

  c['fn_out'] = join(outdir, name+'_'+stripped_name+'.h5')
  try:
    moviename = c['ift_measure']['moviefilename']
  except:
    pass
  else:
    c['ift_measure']['moviefilename'] = join(outdir, name+'_'+moviename+'_'+stripped_name+'.h5')

  print '----------------------------------------'
  print 'submitting %s with params' % c['fn_out']
  myutils.pprint(c)
  print '----------------------------------------'

  run_simple(c['fn_out'], c)
  print('---- run_simple returned ----')



def setvesselfile(c, fn, commonpath, vesselpath, tumorpath):
  import posixpath
  c['fn_tumor'] = fn
  c['h5_path_vessel'] = posixpath.join(commonpath, vesselpath)
  c['h5_path_lattice'] = posixpath.join(commonpath, vesselpath, 'lattice')
  if tumorpath is not None:  
    c['h5_path_tumor'] = posixpath.join(commonpath, tumorpath)

def updated_config(c, p=dict(), piff=dict(), pdrug=dict()):
  c = deepcopy(c)
  c.update(p)
  c['iff'].update(piff)
  c['ift'].update(pdrug)
  return c

if 0:
  paramset_name = args[1]
  files = args[2:-1]
  group = args[-1:]
  group = group[0]
  for fn in files:
    if 0: # parameter variations for paper
      def run2(name, p=dict(), piff=dict()):
        run('', 'iff_'+name, fn, defaultconfig, p, piff)
      cfgs = []
      f = defaultconfig['iff']['iff_capillary_permeability_tumor']
      cfgs += [
        ['default', dict(), dict()],
        ['variantA01', dict(), dict(iff_capillary_permeability_tumor = f * 10.)],
        ['variantA05', dict(), dict(iff_capillary_permeability_tumor = f * 4.)],
        ['variantA02', dict(), dict(iff_capillary_permeability_tumor = f * 2.)],
        ['variantA03', dict(), dict(iff_capillary_permeability_tumor = f * .5)],
        ['variantA09', dict(), dict(iff_capillary_permeability_tumor = f * .25)],
        ['variantA04', dict(), dict(iff_capillary_permeability_tumor = f * .1)],
        ['variantA07', dict(), dict(iff_capillary_permeability_tumor = f * .05)],
        ['variantA08', dict(), dict(iff_capillary_permeability_tumor = f * .01)],
        ['variantA10', dict(), dict(iff_capillary_permeability_tumor = f / 200.)],
        ['variantA11', dict(), dict(iff_capillary_permeability_tumor = f / 500.)],
        ['variantA12', dict(), dict(iff_capillary_permeability_tumor = f /1000.)],
      ]
      f = defaultconfig['iff']['lymph_permeability']
      cfgs += [
        ['variantB07', dict(), dict(lymph_permeability = f * 100.)],
        ['variantB06', dict(), dict(lymph_permeability = f * 50.)],
        ['variantB01', dict(), dict(lymph_permeability = f * 10.)],
        ['variantB05', dict(), dict(lymph_permeability = f * 5.)],
        ['variantB02', dict(), dict(lymph_permeability = f * 2.)],
        ['variantB03', dict(), dict(lymph_permeability = f * .5)],
        ['variantB04', dict(), dict(lymph_permeability = f * .1)],
        ['variantB08', dict(), dict(lymph_permeability = f * .05)],
        ['variantB09', dict(), dict(lymph_permeability = f * .01)],
      ]
      f = defaultconfig['iff']['lymph_surf_per_vol']
      cfgs += [
        ['variantC01', dict(), dict(lymph_surf_per_vol_tumor = f / 100)],
        ['variantC02', dict(), dict(lymph_surf_per_vol_tumor = f / 50)],
        ['variantC03', dict(), dict(lymph_surf_per_vol_tumor = f / 20)],
        ['variantC04', dict(), dict(lymph_surf_per_vol_tumor = f / 10)],
        ['variantC05', dict(), dict(lymph_surf_per_vol_tumor = f / 5)],
        ['variantC06', dict(), dict(lymph_surf_per_vol_tumor = f / 2)],
        ['variantC07', dict(), dict(lymph_surf_per_vol_tumor = f)],
      ]
      f = defaultconfig['iff']['iff_cond_tissue']
      for i, x in enumerate([100, 50, 20, 10, 5, 2]):
        cfgs.append(
          ['variantD%02i' % (i+1), dict(), dict(iff_cond_tissue = f * x, iff_cond_tumor = f * x, iff_cond_necro = f*10 * x)]
        )
      f = defaultconfig['iff']['iff_cond_tissue']
      g = defaultconfig['iff']['lymph_permeability']
      q = list(enumerate([100, 10, 50, 20, 5, 2]))
      for i, x in q[2:]:
        cfgs.append(
          ['variantE%02i' % (i+1), dict(), dict(iff_cond_tissue = f * x, iff_cond_tumor = f * x, iff_cond_necro = f*10 * x, lymph_permeability = g * x)]
        )
      for name, p, piff in cfgs:
        run2(name, p, piff)
    else: # single run with defaults
      ''' here a singel run with
      PARAMETERS is starte'''
      run('', 'iff_defaultconfig2', fn, iff_defaultconfig2, dict(), dict())


if 0: # drug stuff
  group = args[1]
  for fn in args[2:]:

    def run3(name, cfg):
      run('', 'iffdrug_'+name, fn, cfg)

    cfgs = [
       ['variant42', dict(), dict(), dict() ], # timescale 10s and 100min (very close to sinek2009 parameters)

       ['variant31', dict(), dict(), dict(kdiff = 0.16, capillary_permeability_normal = 0.00017*0.01, capillary_permeability_tumor = 0.017*0.01, comprates_k12 = 0.001, comprates_k21 = 0.00002, ) ], # 100x slower diffusion
       ['variant32', dict(), dict(), dict(kdiff = 0.16, capillary_permeability_normal = 0.00017*0.01, capillary_permeability_tumor = 0.017*0.01, comprates_k12 = 0, comprates_k21 = 0, ) ], # 100x slower diffusion
    ]
    f = defaultconfig2['iff']['iff_capillary_permeability_tumor']
    g = defaultconfig2['ift']['capillary_permeability_tumor']
    cfgs += [
             ['variant51', dict(), dict(iff_capillary_permeability_tumor= f * 10), dict(capillary_permeability_tumor = g * 10)],
             ['variant52', dict(), dict(iff_capillary_permeability_tumor= f * 0.1), dict(capillary_permeability_tumor = g * 0.1)],
             ['variant53', dict(), dict(iff_capillary_permeability_tumor= f * 0.001), dict(capillary_permeability_tumor = g * 0.001)],
             ['variantA08', dict(), dict(iff_capillary_permeability_tumor= f * 0.01), dict(capillary_permeability_tumor = g * 0.01)],
              ]
    f = defaultconfig['iff']['lymph_permeability']
    cfgs += [['variant61', dict(), dict(lymph_permeability= f * 10), dict()],
             ['variant62', dict(), dict(lymph_permeability= f * 0.1), dict()]]
    f = defaultconfig['iff']['lymph_surf_per_vol']
    cfgs += [['variant71', dict(), dict(lymph_surf_per_vol_tumor = f * 0.1), dict()],
             ['variant72', dict(), dict(lymph_surf_per_vol_tumor = f * 1), dict()]]
    f = defaultconfig2['iff']['iff_cond_tissue']
    g = defaultconfig2['iff']['lymph_permeability']
    h = defaultconfig2['ift']['kdiff']
    cfgs += [
      ['variantE01', dict(), dict(iff_cond_tissue = f * 10, iff_cond_tumor = f*10, iff_cond_necro = f * 10), dict()],
      ['variantE02', dict(), dict(iff_cond_tissue = f*10, iff_cond_tumor = f*10, iff_cond_necro = f*10, lymph_permeability = g*10), dict(kdiff = h*10)]
    ]
    cfgs.append(
      ['variant91', dict(), dict(), dict(convective_transport_enabled = False) ]
    )
    cfgs = dict((c[0], (c[0], updated_config(defaultconfig2, c[1], c[2], c[3])))  for c in cfgs)

    infusion = dict(inject_t = 24 * 3600, inject_max = 1, inject_mode = 'DF_INJECT_MODE_JUMP')
    longinfusion = dict(inject_t = 1000 * 3600, inject_max = 1, inject_mode = 'DF_INJECT_MODE_JUMP')
    decay = dict(inject_t = 1. * 3600, inject_max = 1, inject_mode = 'DF_INJECT_MODE_EXP')

    def make_infusion((name, cc), is_long):
      c_inf = deepcopy(cc)
      c_inf['ift'].update(longinfusion if is_long else infusion)
      f = c_inf['ift']['comprates_k12']/(c_inf['ift']['comprates_k21']+1.e-13)*c_inf['ift']['inject_max']
      c_inf['ift_measure'].update(exposure_thresholds_reference_conc_in = f)
      name = name+('-linf' if is_long else '-inf')
      return name, c_inf

    def make_bolus((name, cc)):
      c_inj = deepcopy(cc)
      c_inj['ift'].update(decay)
      name = name+'-inj'
      return name, c_inj

    def make_movie((name, cc)):
      cc = deepcopy(cc)
      cc['ift_measure'].update(movieintervall = 60. * 10., moviefilename = 'movie')
      return name, cc

    cfgs2 = [
       make_movie(make_bolus(cfgs['variant42'])),
       make_movie(make_infusion(cfgs['variant42'], True)),
       make_movie(make_bolus(cfgs['variant91'])),
       make_movie(make_bolus(cfgs['variant31'])),
#        make_infusion(cfgs['variant42'], True),
#        make_bolus(cfgs['variantA08']),
#        make_bolus(cfgs['variantE02'])
    ]
    cfgs = cfgs2
    for c in cfgs:
      run3(*c)
        
if not qsub.is_client and __name__=='__main__':
  import argparse
  parser = argparse.ArgumentParser(description='Compute IFF distributions and drugs analyze them.')  
  parser.add_argument('Iffparams', help = 'choose the parameter for the simulation,\n possible configs are found in /krebsjobs/parameters/parameterSetsIff.py')  
  parser.add_argument('tumorFileNames', nargs='*', type=argparse.FileType('r'), default=sys.stdin, help='tumor files to calculate')   
  parser.add_argument('grp_pattern',help='Where to find the tumor. Usually this is somthing with out*')      
  #this enables access to the default values  
  #atest = parser.add_argument('-a', '--analyze', help = 'loop through all files analyze data and make plot', default=False, action='store_true')  
  #parser.add_argument('-m', '--memory', help= 'Memory assigned by the queing system', type=str, default = '2GB')
  goodArguments, otherArguments = parser.parse_known_args()
  qsub.parse_args(otherArguments)
  
  #create filename due to former standards
  filenames=[]
  for fn in goodArguments.tumorFileNames:
    filenames.append(fn.name)
  
  def find_tum(name, obj):
    if 'tumor' in name:
      return True
        
  dirs = set()
  try:
    if not goodArguments.Iffparams in dir(paramIff):
      raise AssertionError('Unknown parameter set %s!' % goodArguments.Iffparams)
    for fn in filenames:
      if not os.path.isfile(fn):
        raise AssertionError('The file %s is not present!'%fn)
        #experiment without tumor
      with h5py.File(fn, 'r') as f:
        if not f.visititems(find_tum):
          raise AssertionError('Are you sure this is a tumor file???')
  except Exception, e:
    print e.message
    sys.exit(-1)
 
  theParams = getattr(paramIff, goodArguments.Iffparams)
  theParams['parameterset_name'] = goodArguments.Iffparams
  for fn in filenames:
    #keep the syntax for compatibility with old stuff
    run('', goodArguments.Iffparams, fn, theParams, dict(), dict(), grp_pattern=goodArguments.grp_pattern)
