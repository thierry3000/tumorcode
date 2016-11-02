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

import numpy as np
import krebsutils as krebs
import extensions
from PIL import Image

import matplotlib
matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt


def make_test_array(c):
  a = np.zeros(2*c, dtype=bool)
  for i in np.ndindex(a.shape):
    r = np.linalg.norm(i-c) / np.linalg.norm(c)
    a[i] = r>0.2 and r<0.5
  return a


def testfill():
  
  c = np.asarray((25, 25, 25))
  a = make_test_array(c)
  x = krebs.flood_fill(a, (0,0,0))
  img = Image.fromArray(255*np.asarray(x[:,:,c[2]],dtype=np.uint8))
  img.save('filled2.png')  
  x = krebs.flood_fill(a, c)
  img = Image.fromArray(255*np.asarray(x[:,:,c[2]],dtype=np.uint8))
  img.save('filled1.png')  
  img = Image.fromArray(255*np.asarray(a[:,:,c[2]],dtype=np.uint8))
  img.save('tobefilled.png')

def testdistancemap():
  c = np.asarray((25, 25, 25))
  a = make_test_array(c)
  x = krebs.distancemap(a)
  print "range: %f %f" % (x.min(), x.max())

  grad = krebs.field_gradient(x, spacing=1.)
  for i, g in enumerate(grad):
    print "grad range %i: %f %f" % (i, g.min(), g.max())
  gnorm = np.sum(g*g for g in grad)
  
  for i in [c[2]]: #xrange(x.shape[2]):
    plt.subplot(221)
    p = plt.imshow(x[:,:,i], vmin=x.min(), vmax=6, interpolation='nearest')
    plt.colorbar(p)
    plt.subplot(222)
    p = plt.imshow(a[:,:,i], interpolation='nearest')
    plt.colorbar(p)
    plt.subplot(223)
    p = plt.imshow(grad[0][:,:,i], interpolation='nearest')
    plt.colorbar(p)
    plt.subplot(224)
    p = plt.imshow(gnorm[:,:,i], interpolation='nearest')
    plt.colorbar(p)
    plt.show()
  


def testfieldsampling():
  c = np.asarray((40, 20, 1))
  bb = np.vstack((-c,c-1)).transpose().reshape((6,))
  ld = krebs.LatticeDataQuad3d(bb, 1.)
  ld.SetCellCentering((True, True, True))
  wbb = ld.worldBox
  print ld  , wbb
  
  num = (bb[1]-bb[0])*(bb[3]-bb[2])
  pos = np.random.uniform(wbb[0], wbb[1], num), np.random.uniform(wbb[2], wbb[3], num), 0.5*np.ones(num)
  pos = np.asarray(pos,dtype=np.float32).transpose()

  a = np.asarray(make_test_array(c), dtype=np.float32)

  smpl = krebs.sample_field(pos, a, ld, linear_interpolation=True)

  print 'arange = %f %f' % (a.min(), a.max())
  print 'smpl range = %f %f' % (smpl.min(), smpl.max())

  img = a[:,::-1,1]
  img = img.transpose()
  plt.imshow(img, extent = ld.worldBox, interpolation='nearest', vmin=0, vmax=1)

  plt.scatter(pos[:,0], pos[:,1], c=smpl, vmin=0., vmax=1.)
  
  plt.show()


def testcorrelation():
  a1 = np.zeros((50, 50))
  a1[::3,::3] = 1
  #a1[:,:] = 1

  a2 = a1
  #a2 = np.zeros((50, 50))  
  #a2[10:20, 10:20] = 1
  
  r, c, cerr = krebs.radial_correlation(a1, a1, 10, bins_per_unit=100)
  r4, c4, c4err = krebs.radial_correlation(a1, a1, 10, bins_per_unit=2)

  from scipy.signal  import gaussian
  from scipy.signal  import convolve
  gfilter = gaussian(100, 25.)
  gfilter /= np.sum(gfilter)
  c2 = convolve(c, gfilter, mode='same')
  
  ax = plt.subplot2grid((2,2), (0, 0))
  ax.imshow(a1, interpolation='none')
  
  ax = plt.subplot2grid((2,2), (0, 1))  
  ax.imshow(a2, interpolation='none')  
  
  ax = plt.subplot2grid((2,2), (1, 0), colspan=2)
  
  ax.errorbar(r, c, cerr, color='0.8')
  
  ax.plot(r, c2, 'g-', linewidth = 0.5)
  
  ax.errorbar(r4, c4, c4err, linewidth = 2.)
  
  ax.set(ylim = (-0.1, 0.1))
  
  plt.show()


def run():
  #testfill()
  #testdistancemap()
  #testfieldsampling()
  testcorrelation()
  
if __name__ == '__main__':
  run()