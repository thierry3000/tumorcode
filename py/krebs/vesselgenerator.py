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
#! /usr/bin/env python2
# -*- coding: utf-8 -*-
import os, sys
from os.path import join, basename, dirname, splitext
sys.path.append(join(dirname(__file__),'..'))

import math
import random
import numpy as np
import cStringIO
from mystruct import Struct
# for checksum generation
import zlib
import cPickle
# for tower sampling
import bisect
# to run vesselgen
import krebsutils

# well defined random seed
random.seed(123456)
#
#num_threads = 2

class VD(Struct):
  '''represents parameters for the vessel generator'''
  def __init__(self, **kwargs):
    Struct.__init__(self)
    for k, v, in dict(name=None, shape=None, scale=80., roots=[], seed=0, latticetype="quad",
                      num_hierarchical_iterations=0,
                      tip_radius_arterial = 3.5,
                      tip_radius_capi = 3.5,
                      tip_radius_vein = 4.5,
                      max_sprout_radius_artery = 8.,
                      max_sprout_radius_vein = 8., # these are the max radi to which capillaries can connect
                      murray_alpha_vein = 3.,
                      murray_alpha_artery = 3.,
                      capillariesUntilLevel = 0,
                      o2range = 300.,
                      outfilename = 'unnamed',
                      full_debug_output = False,
                      ensemble_index = 0,
                      generate_more_capillaries = False,
                      calcflow = None,
                      num_threads = 1,
                      changeRateThreshold = 1.e-3,  # this may be too low for 2d configurations!!!
                      ).iteritems():
      self[k] = kwargs.pop(k, v)
    assert not kwargs
    self.__setattr__ = VD.checked_setattr

  @staticmethod
  def checked_setattr(self, k, v):
    assert k in self
    Struct.__setattr__(self, k, v)

  @property
  def num_iter(vd):
    """
      estimates maximum number of iterations
      proportional to largest lateral grid size
    """
    return max(*vd.shape)*100

  def generate_info_string(vesseldata):
      '''
        write a config as boost .info file,
        to be read by the boost::property_tree library
      '''
      sio = cStringIO.StringIO()
      def a(s, *args):
        sio.write(s % args)
        sio.write('\n')
      size, scale, roots = vesseldata.shape, vesseldata.scale, vesseldata.roots
      a('lattice_size "<%i, %i, %i>"', size[0], size[1], size[2])
      a('lattice_spacing %f', scale)
      a('lattice_type "%s"', vesseldata.latticetype)
      a('seed %u', vesseldata.seed)
      a('message "%s"', vesseldata.name)
      a('out_fn "%s"', vesseldata.outfilename)
      a('full_debug_output %s', 'true' if vesseldata.full_debug_output else 'false')
      a('max_num_iter %i', vesseldata.num_iter)
      a('num_threads %i', vesseldata.num_threads)
      a('num_hierarchical_iterations %i', vesseldata.num_hierarchical_iterations)
      a('max_sprout_radius_artery %s', vesseldata.max_sprout_radius_artery)
      a('max_sprout_radius_vein %s', vesseldata.max_sprout_radius_vein)
      a('radius_vein %s', vesseldata.tip_radius_vein)
      a('radius_capi %s', vesseldata.tip_radius_capi)
      a('radius_artery %s', vesseldata.tip_radius_arterial)
      a('capillariesUntilLevel %i', vesseldata.capillariesUntilLevel)
      a('changeRateThreshold %f', vesseldata.changeRateThreshold)
      #a('radius_capi 4.')
      #a('radius_artery 4.5')
      a('murray_alpha_vein %s', vesseldata.murray_alpha_vein)
      a('murray_alpha_artery %s', vesseldata.murray_alpha_artery)
      a('ensemble_index %i' % vesseldata.ensemble_index)
      a('generate_more_capillaries %i' % vesseldata.generate_more_capillaries)
      a('o2 {')
      a('  diffusion_range %s', vesseldata.o2range)
      a('}')
      if vesseldata.calcflow:
        a('calcflow {')
        a('  viscosityPlasma %s', vesseldata.calcflow['viscosityPlasma'])
        a('  rheology %s', vesseldata.calcflow['rheology'])
        a('  inletHematocrit %s', vesseldata.calcflow['inletHematocrit'])
        a('  includePhaseSeparationEffect %s', vesseldata.calcflow['includePhaseSeparationEffect'])
        a('}')
      a('roots {')
      for q in roots:
        (x,y,z),t = q[:2]
        x,y,z = x%size[0],y%size[1],z%size[2]
        assert t in ['a','v'], 't must be a or v'
        a('"" {')
        a('    p "<%i, %i, %i>"',x ,y, z)
        a('    t %s', t)
        if len(q)>2:
          a('    l %i', q[2])
        else:
          a('    l 1')
        a('  }')
      a('}')
      return sio.getvalue()


##==============================================================================
##              helper functions
##==============================================================================

#def fix_dim(vd):
#  """
#    make the shape 3d
#  """
#  s = vd.shape
#  if len(s) == 1:  # a 2d system given by a single lateral size
#    vd.shape = (s[0],s[0])
#  elif len(s) == 2:
#    vd.shape = (s[0],s[0], s[1])  # a 3d system given by a lateral size in the ground plane and z-height
#  return vd


def fix_worldsize(vd):
  f = 1./0.8 if vd.latticetype == 'fcc' else 1
  vd.orgshape = vd.shape
  x,y,z = vd.shape
  y = int(f*y+0.5)
  if z*2 > x:
    z = int(f*z+0.5) # cubic case
  vd.shape = (x,y,z)
  return vd


def make_root_length_random(vd, l_min_frac, l_max_frac):
  def get_root_axis(root, shape):
    pos = root[0]
    p = np.abs(np.asarray(pos,dtype=float)-np.asarray(shape,dtype=float)/2.)
    return np.argmax(p);

  for i, root in enumerate(vd.roots):
    axis = get_root_axis(root, vd.shape)
    sys_len = vd.shape[axis]
    l_min = l_min_frac * sys_len
    l_max = l_max_frac * sys_len
    l = random.randint(int(l_min), int(l_max))
    vd.roots[i] = root[:2]+(l,)
  return vd


def make_fn(vd):
  def subd_(s, n):
    return subd_((s-1)*2+1, n-1) if n>0 else s
  subd = lambda s: subd_(s, vd.num_hierarchical_iterations)
  sx = subd(vd.shape[0])
  sz = vd.shape[2]
  if sz > 1:
    sz = subd(sz)
  return '%s-%ix%iL%i' % (vd.name, sx, sz, vd.scale)



def weighted_choice(items, weights):
  '''
    tower sampling -- returns the selected item
    tested with

    bins = np.asarray([0,0,0], dtype=np.float64)
    N = 1000
    for i in xrange(N):
      q = weighted_choice([0,1,2], [1,2,3])
      bins[q] += 1
    print bins / N * 6.
  '''
  w = np.asarray(weights, dtype=np.float64)
  cw = np.cumsum(w)
  cw = cw[:-1]/cw[-1]
  r = random.random()
  i = bisect.bisect(cw, r)
  return items[i]


def checksum(obj):
  ''' python object to checksum, used to generate random seeds from some configuration '''
  return zlib.crc32(cPickle.dumps(obj)) & 0xffffffff


def getAbsolutePosition(pos, shape):
  ''' from range [0, 1] to lattice index'''
  x,y,z = pos
  eps = 1.0e-5
  return int((x-eps)*shape[0])%shape[0], int((y-eps)*shape[1])%shape[1], int((z-eps)*shape[2])%shape[2]


def getBorderMask(shape, num_axes):
  '''returns a integer nd-array where border sites get 1, sites in the bulk get 0'''
  realdim = len(shape)
  assert realdim in (2,3), "must be 2 or 3 dimensional"
  assert 0 <= num_axes <= realdim, "for dim it must be 2 <= dim <= realdim"
  assert shape[1] >= 3 or num_axes < 2
  assert shape[2] >= 3 or num_axes < 3
  mask = np.ones(shape, dtype=np.int8)
  if num_axes >= 1:
      innermask = mask[1:-1,...]
  else:
      innermask = mask[1:,...]
  if num_axes >= 2:
      innermask = innermask[:,1:-1,...]
  if num_axes >= 3:
      innermask = innermask[...,1:-1]
  innermask[...] = 0
  return mask


def iterBorderSites(shape, num_axes):
  mask = getBorderMask(shape, num_axes)
  for i, x in np.ndenumerate(mask):
    if not x: continue
    yield i


def generate_random_distributed_roots(shape, count, axes, sides, min_spacing):
  """
    random distributed points on the boundary of the rectangular region

    test:
    l = generate_random_distributed_roots((100, 100, 100), count=10, axes=3, sides=2, min_spacing=2)
    for x in l:
      print x

    axes and sides --  integers determining at which faces samples are generated
                       i.e. axes=1, sides=2, means at left and right face
  """
  l = []
  al = np.asarray(shape, dtype=np.float32)
  face_weights = (al[1]*al[2], al[0]*al[2], al[0]*al[1])[:axes]
  axes = np.asarray([0,1,2])[:axes]
  for i in range(count):
    for trial_num in range(100):
      axis = weighted_choice(axes, face_weights)
      side = random.randint(0,sides-1)
      axu = (axis+1)%3
      axv = (axis+2)%3
      u, v = random.randint(1, shape[axu]-2), random.randint(1, shape[axv]-2)
      p = np.zeros((3,), dtype=np.int)
      p[axis] = 0 if side==0 else shape[axis]-1
      p[axu] = u
      p[axv] = v
      close = False
      for x in l:
        if np.linalg.norm(x-p) < min_spacing:
            close = True
            break
      if not close:
        break
    if close:
        raise RuntimeError('failed to find place for source %i after 100 trials' % i)
    else:
        l.append(p)
  return l

##==============================================================================
##              top level tree generation functions
##==============================================================================

def generate_onesided_alternating_config(vd):
  roots = []
  for p in iterBorderSites(vd.shape, 1):
    if p[0] > vd.shape[0]/2: continue
    q = reduce((lambda a,b: a^(b&1)), p, 1)
    t = 'a' if q else 'v'
    roots.append((p, t))
  vd.roots = roots
  vd.seed  = checksum((vd.shape, vd.scale, vd.ensemble_index, vd.seed))
  return vd


def generate_onesided_alternating_config_w_stems(vd, l_min, l_max):
  vd = generate_onesided_alternating_config(vd)
  vd = make_root_length_random(vd, l_min, l_max)
  vd.seed = checksum((vd.seed, vd.roots))
  return vd


def generate_alternating_config(vd, axes):
  roots = []
  for p in iterBorderSites(vd.shape, axes):
    q = reduce((lambda a,b: a^(b&1)), p, 1)
    t = 'a' if q else 'v'
    roots.append((p, t))
  vd.roots = roots
  vd.seed  = checksum((vd.shape, vd.scale, vd.ensemble_index, vd.seed))
  return vd


def generate_alternating_config_w_stems(vd, axes, l_min, l_max):
  """
    roots on all side faces specified by axes = 1,2 or 3
    with fractional length l_min to l_max
  """
  vd = generate_alternating_config(vd, axes)
  vd = make_root_length_random(vd, l_min, l_max)
  vd.seed = checksum((vd.seed, vd.roots))
  return vd


def generate_alternating_diluted_config(vd):
  vd = generate_alternating_config(vd, 2)
  random.seed(vd.seed)
  np.random.seed(vd.seed)
  roots = np.asarray(vd.roots, np.object)
  while True:
    r = random.randint(2, len(roots))
    #indices = np.random.choice(len(roots), size=r, replace=False)
    #new_roots = roots[indices]
    new_roots = random.sample(roots, r)
    #print 'select %i of %i roots'  % (r, len(roots))
    has_artery = any((t == 'a') for (p,t) in new_roots)
    has_vein   = any((t == 'v') for (p,t) in new_roots)
    if not has_artery or not has_vein:
      continue
    else:
      break
  vd.roots = list(new_roots)
  return vd
  
  #filter(lambda (p,t): t=='a', vd.roots), np.object)
  #filter(lambda (p,t): t=='v', vd.roots), np.object)
  

def generate_random_config(vd, count, axes, sides):
  vd.seed = checksum((vd.shape, vd.scale, count, axes, sides, vd.ensemble_index, vd.seed))
  random.seed(vd.seed)
  fl = np.asarray(vd.shape, dtype=np.float32)
  fl = np.asarray((min(fl[1],fl[2]), min(fl[0],fl[2]), min(fl[0],fl[1])))
  min_spacing = 0.25*np.amin(fl[:axes])/count
  pos = generate_random_distributed_roots(vd.shape, count, axes, sides, min_spacing)
  roots = [ (p, 'a') for p in pos[:len(pos)/2] ] + [ (p, 'v') for p in pos[len(pos)/2:] ]
  vd.roots = roots
  return vd


configs1 = [
    ('baumlecfg01',[ (0.,0.33,0.33,'a'),(0.,0.66,0.66,'v') ] ), # a/v pair on one side
    ('baumlecfg02',[ (0.,0.25,0.33,'a'),(1.,0.75,0.66,'v') ] ),  # a/v pair on opposing sides
    ('baumlecfg03',[ (0.,0.25,0.33,'a'),(1.,0.5,0.66,'v'),(0.,0.75,0.33,'a') ] ), # 2a/1v on opposing sides
    ('baumlecfg04',[ (0.,0.25,0.33,'v'),(1.,0.5,0.66,'a'),(0.,0.75,0.33,'v') ] ), # 1a/2v on opposing sides
    ('baumlecfg05',[ (0.,0.25,0.33,'a'),(0.,0.5,0.66,'v'),(0.,0.75,0.33,'a') ] ), # 2a/1v on same sides
    ('baumlecfg06',[ (0.,0.25,0.33,'v'),(0.,0.5,0.66,'a'),(0.,0.75,0.33,'v') ] ), # 1a/2v on same sides
    ('baumlecfg07',[ (0.,0.33,0.33,'a'),(0.66,0.,0.66,'v'),(1.,0.66,0.33,'a'),(0.66,1.,0.66,'v') ] ), # a/v/a/v in circular arrangement
    ('baumlecfg08',[ (0.,1./5,0.33,'a'),(1.,2./5,0.66,'v'),(0.,3./5,0.66,'a'),(1.,4./5,0.33,'v') ] ), # 2a/2v on  opposing sides
    ('baumlecfg09',[ (0.,1./5,0.33,'a'),(1.,2./5,0.66,'a'),(0.,3./5,0.66,'v'),(1.,4./5,0.33,'v') ] ), # 2 x a/v pairs on  opposing sides
    ('baumlecfg10',[ (0.,1./5,0.33,'a'),(1.,2./5,0.66,'v'),(0.,3./5,0.66,'v'),(1.,4./5,0.33,'a') ] ), # 2 x a/v pairs on  opposing sides - different permutation
    ('baumlecfg11',[ (0.,0.1,0.33,'a'),(0.,0.9,0.66,'v') ] ), # a/v pair on one side at the edges
    ('baumlecfg12',[ (0.,0.1,0.33,'a'),(1.,0.9,0.66,'v') ] ),  # a/v pair on opposing sides at the edges
]
configs1dict = dict((t[0], t) for t in configs1)

def generate_handmade1_config(vd, cfg):
  cfg_name, roots = cfg
  name = vd.name
  vd.roots = [ (getAbsolutePosition(q[:3], vd.shape), q[3]) for q in roots ]
  vd.seed = checksum((vd.shape, vd.scale, vd.roots, vd.ensemble_index, vd.seed))
  return vd


def generate_handmade1_config_w_stems(vd, cfg, l_min, l_max):
  vd = generate_handmade1_config(vd, cfg)
  vd = make_root_length_random(vd, l_min, l_max)
  vd.seed = checksum((vd.seed, vd.roots))
  return vd

def generate_symmetric_config():
  krebsutils.vesselgen_generate_symmetric()
  
if __name__ == '__main__':
  paramset11 = dict(
        tip_radius_arterial = 2.5,
        tip_radius_capi = 2.5,
        tip_radius_vein = 3.8,
        murray_alpha_vein = 3.,
        murray_alpha_artery = 3.,
        scale = 90.,
        max_sprout_radius_artery = 4.,
        max_sprout_radius_vein = 4.,
        calcflow = dict(
          viscosityPlasma = 1.2e-6,
          rheology = 'RheologyForHuman',
          inletHematocrit = 0.45,
          includePhaseSeparationEffect = False,
        ),
        changeRateThreshold = 1.e-3
      )  
  if 0:
    paramset11.update(
      shape=(50,50,3),
      latticetype='fcc',
      num_hierarchical_iterations=0,
      name='vessels-gf', 
      ensemble_index=0,
      seed = 1298762,
      o2range = 200.,
      max_sprout_radius_artery = 20.,
      max_sprout_radius_vein = 20.,
    )
    paramset11.update(
      name='vessels-q2d-gf-hc', 
      shape=(32,32,2),
      num_hierarchical_iterations = 2,
    )
    paramset11.update(
      name='vessels-2d-gf-hc', 
      shape=(32,32,1),
      num_hierarchical_iterations = 2,
      # 2d sims take longer to converge since branches have to 
      # squeeze through narrow "corridors".     
      changeRateThreshold = 0.1e-3, # lower than for 3d networks.
    )
    vd = VD(**paramset11)
    vd = fix_worldsize(vd)
    vd = generate_handmade1_config(vd, configs1dict['baumlecfg12'])
    
    #vd = generate_alternating_config(vd, 1)
    #vd = generate_handmade1_config_w_stems(vd, configs1dict['baumlecfg01'], vd.shape[0]*0.5, vd.shape[0]*0.8)
    vd.outfilename = vd.name
    #vd.full_debug_output = True
    configstring = vd.generate_info_string()
    print 'running: ', configstring
    krebsutils.run_vesselgen(configstring)
  if 0:
    size = (10, 10, 10)
    scale = 150.
    pgrad = 0.2
    p = pgrad*size[0]*scale
    r = 6.
    krebsutils.vesselgen_generate_single('vess-single-capillary.h5', size, scale, 'quad', r, 0., p)

  if 1:
    scale = 150.
    exponent_of_two = 1
    import h5py
    myFile = h5py.File('vess-cub.h5','w')
    krebsutils.vesselgen_generate_symmetric(myFile, exponent_of_two, scale, paramset11['calcflow'], False)
