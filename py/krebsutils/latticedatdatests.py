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

"""
  this is stuff i used to develop and test the
  FCC lattice representation (LatticeDataFCC)
"""


import numpy as np
import krebsutils as krebs
import math
import collections


def TestBasics():
  x0, x1, y0, y1, z0, z1 = bb = (0,4, 0,4 ,-3,6)
  ld = krebs.LatticeDataQuad3d(bb, 100.)
  print ld.GetWorldBox()
  print ld.worldBox
  print ld.GetBox()
  print ld.box
  print ld.scale
  print ld.GetScale()
  print ld.shape
  print ld.GetWorldSize()
  print ld.worldSize

  ld2 = krebs.LatticeDataQuad3d(bb, 100.)  
  print ld == ld2
  
  ld3 = krebs.LatticeDataQuad3d((0,3,0,4,-3,6), 100.)
  print ld == ld3
  
  

def generate_unitcell(type):
  if type == 'cubic':
    l = list(np.ndindex(2,2,2))
  elif type == 'fcc':
    l = list(np.ndindex(2,2,2))
    l += [(1,0.5,0.5),(0,0.5,0.5),
	  (0.5,1,0.5),(0.5,0,0.5),
	  (0.5,.5,1),(0.5,0.5,0)]
  elif type == 'bcc':
    l = list(np.ndindex(2,2,2))
    l += [(0.5, 0.5, 0.5)]
  return np.asarray(l, dtype=float)


def plotunitcells(cell, mat, scale):
  import matplotlib
  matplotlib.use('Qt4Agg')
  from mpl_toolkits.mplot3d import Axes3D
  import matplotlib.pyplot as plt
  pts = []
  for p0 in np.ndindex(5,5,5):
    for pu in cell:
      p = np.asarray(p0) + np.asarray(pu)
      delta = np.asarray(tuple(sum(pi*mij for (pi,mij) in zip(p, mj)) for mj in mat))
      p += delta
      p = np.asarray(tuple(pi*si for (pi,si) in zip(p, scale)))
      pts.append(p)
  pts = np.asarray(pts)

  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')
  ax.plot(pts[:,0],pts[:,1],pts[:,2],'*')
  bd = (0, np.amax(pts))
  ax.set(xlim=bd, ylim=bd, zlim=bd)
  plt.show()


def testunicells():
  hexa_c0 = math.sqrt(6.)/3.
  hexa_c1 = math.sqrt(3.)/2.
  hexa_c2 = hexa_c1 / 3.
  hexa_mat = ((0, 0.5, 0.5), (0, 0, 1./3.), (0,0,0))
  hexa_scale = (1, hexa_c1, hexa_c0)
  print hexa_scale
  plotunitcells(generate_unitcell('bcc'), hexa_mat, hexa_scale)


def boxindex(bb):
  x0, x1, y0, y1, z0, z1 = bb
  for (i,j,k) in np.ndindex(x1-x0+1, y1-y0+1, z1-z0+1):
    yield (i+x0, j+y0, k+z0)

fcc_points_of_interest =  [ (0,0,0), (0,0,1), (0,0,2),
                            (0,1,0), (0,1,1), (0,1,2) ]

def generatefccneighbors():
  x0, x1, y0, y1, z0, z1 = bb = (-3,3, -3,3 ,-3,3)
  ld = krebs.LatticeDataFCC(bb, 1.)
  for p0 in fcc_points_of_interest:
    nbs = []
    for p1 in boxindex(bb):
      wp1 = ld.LatticeToWorld(p1) - ld.LatticeToWorld(p0)
      if (np.linalg.norm(wp1)-1.)**2 < 0.1:
        nbs.append(np.asarray(p1) - np.asarray(p0))
    nbs.sort(key = lambda (i,j,k): k*100+10*j+i)
    if 0:
      print '%i nbs of %s = {' % (len(nbs), p0,)
      print ','.join(str(p) for p in nbs)
      print '}'
    else:
      x,y,z = p0
      print 'k=0;'
      for num, (i,j,k) in enumerate(nbs):
        print 'vnb[%i][%i][k++] = Int3(%2i,%2i,%2i);' % (y,z,i,j,k)

def testhcplattice():
  x0, x1, y0, y1, z0, z1 = bb = (0,4, 0,4 ,-1,3)
  ld = krebs.LatticeDataFCC(bb, 1.)
  print type(ld), ld

  if 1:
    print 'neighbors'
    nbcnt = ld.NbCount()
    for i in xrange(nbcnt):
      p = ld.NbLattice((0,0,0), i)
      wp = ld.LatticeToWorld(p)
      print '%i -> %s = <%.1f, %.1f, %.1f>' % ((i, p)+tuple(wp))

  import matplotlib
  matplotlib.use('Qt4Agg')
  from mpl_toolkits.mplot3d import Axes3D
  import matplotlib.pyplot as plt

  def plotsiteneighbors(ax, ld, p, axes, **kwargs):
    wp = ld.LatticeToWorld(p)
    for i in xrange(ld.NbCount()):
      wq = ld.LatticeToWorld(ld.NbLattice(p, i))
      a = np.asarray([wp,0.5*(wq-wp)+wp])
      if len(axes)==3:
        ax.plot(a[:,axes[0]], a[:,axes[1]], a[:,axes[2]], **kwargs)
      else:
        ax.plot(a[:,axes[0]], a[:,axes[1]], **kwargs)


  def plotneighbors(fig, subplotparams, axes):
    colors = 'rgbcmyk'
    if len(axes) == 3:
      ax = fig.add_subplot(subplotparams, projection='3d')
    else:
      ax = fig.add_subplot(subplotparams)
    points = collections.defaultdict(list)
    for num, (i,j,k) in enumerate(np.ndindex(x1-x0+1, y1-y0+1, z1-z0+1)):
      i,j,k = (i+x0, j+y0, k+z0)
      plotsiteneighbors(ax, ld, (i,j,k), axes=axes, color=colors[k % len(colors)])
      points[k].append(ld.LatticeToWorld((i,j,k)))
    for k, pts in points.iteritems():
      pts = np.asarray(pts)
      if len(axes)==3:
        ax.plot(pts[:,axes[0]], pts[:,axes[1]], pts[:,axes[2]], '.', color=colors[k % len(colors)], markersize=20)
      else:
        ax.plot(pts[:,axes[0]], pts[:,axes[1]], '.', color=colors[k % len(colors)], markersize=20)

    wbb = ld.worldBox

    wbb = wbb.reshape((3,2))
    ax.set_xlim(*wbb[axes[0]])
    ax.set_ylim(*wbb[axes[1]])
    if len(axes) == 3:
      ax.set_zlim(*wbb[axes[2]])
    ax.set_aspect(1.)
    axlabels = ['x','y','z']
    ax.set_xlabel(axlabels[axes[0]])
    ax.set_ylabel(axlabels[axes[1]])

  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')
  plotsiteneighbors(ax, ld, (0,0,0), (0,1,2), color = 'k')
  ax.set(xlabel = 'x', ylabel = 'y', zlabel = 'z')
  plt.show()

#  fig = plt.figure()
#  plotneighbors(fig, 221, (0,1))
#  plotneighbors(fig, 222, (0,2))
#  plotneighbors(fig, 223, (1,2))
#  plotneighbors(fig, 224, (0,1,2))
#  plt.show()


def printFCCLatticeDirections():
  LD = krebs.LatticeDataFCC
  sz = (4,4,4)
  ld = LD((0, sz[0]-1, 0, sz[1]-1, 0, sz[2]-1), 1.)
  for d in xrange(ld.NbCount()):
    indexv = ld.NbLattice((0,0,0), d)
    wp     = ld.LatticeToWorld(indexv)
    print '%i -> %s' % (d, wp)


def march(ld, p, d, l):
  t1 = p
  for q in xrange(l):
    #print t1,
    t1 = ld.NbLattice(t1, d)
  #print t1
  return t1


def testhcplattice2():
  LD = krebs.LatticeDataFCC
  sz = (4,4,4)
  ld = LD((0, sz[0]-1, 0, sz[1]-1, 0, sz[2]-1), 1.)

  print ld.GetAxisDirLen((2,15,1),(3,14,0))

  if 0:
    for p0 in fcc_points_of_interest:
      if 1:
        for d_test in xrange(ld.NbCount()):
          p1 = ld.NbLattice(p0, d_test)
          d, l = ld.GetAxisDirLen(p0, p1)
          if d is not None:
            t1 = march(ld, p0, d, l)
            ok = all(t1 == p1)
            print '%s -> %s has dir %s and length %s, marching result: %s, %s' % (str(p0),str(p1), d, l, str(t1), 'ok' if ok else 'FAIL')
    if 0:
      inner_region = set([tuple(ld.NbLattice(p0, d1)) for d1 in xrange(ld.NbCount())] + [p0])
      for d1 in xrange(ld.NbCount()):
        p1 = ld.NbLattice(p0, d1)
        for d2 in xrange(ld.NbCount()):
          p2 = ld.NbLattice(p1, d2)
          if tuple(p2) in inner_region: continue
          d, l = ld.GetAxisDirLen(p0, p2)
          if d is not None:
            t1 = march(ld, p0, d, l)
            if d<>d1 or d<>d2 or any(t1<>p2):
              print 'FAIL','%s -> %s has dir (%s,%s)->%s and length %s, marching result: %s' % (str(p0),str(p2), d1,d2 , d, l, str(t1))
            else:
              print 'OK','%s -> %s has dir %s and length %s'  % (str(p0),str(p2), d, l)


if __name__ == '__main__':
  TestBasics()
  #generatefccneighbors()
  #testhcplattice()
  #testhcplattice2()
  #printFCCLatticeDirections()
