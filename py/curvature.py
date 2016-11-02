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
# -*- coding: utf-8 -*-
#!/usr/bin/env python

import sys
import numpy as np
import numpy.linalg as linalg
import math

_do_dbg_out = False
def _dbg_out(s, *args):
  if _do_dbg_out:
    sys.stdout.write(s)
    sys.stdout.write('\n')


def estimateFrameAL(points):
  N, dim = points.shape
  assert N&1 == 1, "number of points must be uneven"
  assert N >= 3, "need more than 3 points"
  ic = N/2
  p0 = points[ic].copy()
  for i in xrange(N):
    points[i] -= p0
  #np.subtract(points, p0, points)
  l = np.empty(N, dtype=np.double)
  w = np.empty(N, dtype=np.double)
  l[ic] = 0.
  w[ic] = 1.
  for dir in (-1, 1):
    for i in xrange(N/2):
      i1 = ic + dir * i
      i2 = ic + dir * (i+1)
      li = dir * linalg.norm(points[i1]-points[i2])
      assert abs(li)>1.e-18, "degenerate points"
      l[i2] = l[i1] + li
      #w[i2] = 1.
      #w[i2] = 1./(i+2)
      w[i2] = math.exp(-float(i+1)*4./N)
  _dbg_out(str(points))
  _dbg_out('l = %s' % l)
  _dbg_out('w = %s' % w)
  dr = np.zeros((3,))
  ddr = np.zeros_like(dr)
  A = np.zeros((2,2), dtype = np.double)
  t = np.zeros((2,), dtype = np.double)
  A[0,0] = np.sum(w * np.power(l, 2.))
  A[0,1] = 0.5 * np.sum(w * np.power(l, 3.))
  A[1,1] = 0.25 * np.sum(w * np.power(l, 4.))
  A[1,0] = A[0,1]
  _dbg_out('A = %s' % A)
  for d in xrange(dim):  
    t[0] = np.sum(w * l * points[:,d])
    t[1] = 0.5 * np.sum(w * np.power(l, 2.) * points[:,d])
    _dbg_out('t[%i] = %s' % (d, t))
    dr[d], ddr[d] = linalg.solve(A, t)
#    det = A[0,0]*A[1,1] - A[0,1]*A[1,0]
#    dr[d] = (A[1,1]*t[0] - A[0,1]*t[1])/det
#    ddr[d]= (A[0,0]*t[1] - A[1,0]*t[0])/det
  _dbg_out('dr  = %s,\nddr = %s' % (dr, ddr))
  return dr, ddr


def estimateCurvatureAL(points):
  dr, ddr = estimateFrameAL(points)
  curv = linalg.norm(ddr)
  #curv = linalg.norm(np.cross(dr, ddr)/math.pow(linalg.norm(dr), 3))
  #_dbg_out('curv = %f' % curv)
  return curv  


#def estimateCurvature(points):
#  _dbg_out("first estimate")
#  dr, ddr = estimateFrame(points)
#  ndr = linalg.norm(dr)
#  nddr = linalg.norm(ddr)
#  if nddr <= 1.e-19:
#    return 0.
#  dr *= 1./ndr
#  ddr *= 1./nddr
#  _dbg_out("lengths= %s, %s" % (ndr, nddr))
#  z = np.cross(dr, ddr)
#  N, dim = points.shape
#  M = dr[:dim], ddr[:dim], z
#  _dbg_out("frame\nx=%s, y=%s, z=%s" % M)
#  tpoints = np.empty_like(points)
#  for d in xrange(dim):
#    for i in xrange(N):
#      tpoints[i][d] = np.dot(M[d],points[i])
#  _dbg_out("real estimate")
#  return estimateCurvatureSimple(tpoints)      


#def estimateFrameAndProject1(points):
#  N, dim = points.shape
#  assert N&1 == 1, "number of points must be uneven"
#  assert N >= 3, "need more than 3 points"
#  ic = N/2
#  p0 = points[ic].copy()
#  for i in xrange(N):
#    points[i] -= p0
#  if dim == 3:
#    normal, _, _, _ = linalg.lstsq(points, np.zeros(N))
#    i = np.argmin(np.abs(normal))
#    q = np.zeros(3)
#    q[i] = 1
#    X = np.cross(normal, q)
#    Y = np.cross(normal, X)
#    M = (X,Y)
#    pp = np.empty(N, 2)
#    for d in xrange(2):
#      for i in xrange(N):
#        pp[i][d] = np.dot(M[d],points[i])    
#  return pp


class CurvEst:
  def __init__(self, points):
    self.org_points = points
    N, dim = self.N, self.dim = points.shape
    assert N&1 == 1, "number of points must be uneven"
    assert N >= 3, "need more than 3 points"
    ic = self.ic = N/2
    self.p0 = points[ic].copy()
  def estimateFrameByCovarianceTensor(self):
    N, dim = self.N, self.dim
    M = np.zeros((dim, dim))
    w_sum = 0.
    for i in xrange(N):
      w = math.exp(-float(abs(i-self.ic))*2./N)
      p = self.org_points[i] - self.p0
      M += w * np.outer(p,p)
      w_sum += w
#      print 'w%i=%f' % (i, w)
    M *= 1./w_sum
    ew, v = linalg.eigh(M)
    f = zip(ew, v)
    f.sort(key = lambda (ew,v): -ew)
    ew, v = zip(*f)
    self.frame_scale_inv = 1./math.sqrt(ew[0])
    self.frame = np.asarray(tuple(self.frame_scale_inv * v[i] for i in xrange(dim)))
    return self.p0, ew, v
  def projectPointsToFrame(self):
    N, dim = self.N, self.dim
    pp = np.empty((N, 2))
    for d in xrange(2):
      for i in xrange(N):
        pp[i][d] = np.dot(self.frame[d], self.org_points[i] - self.p0)    
    self.points = pp
    ycenter = np.sum(self.points[:,1])
    if ycenter < 0.:
      self.points[:,1] *= -1.
  def computeCurvatureByParabolicFit(self):
    #r = np.sum(self.points[:,1] * np.power(self.points[:,0], 2.))
    #a = 0.5 * np.sum(np.power(self.points[:,0], 4.))
    #k = abs(r/a)
    A = np.asarray((self.points[:,0], np.power(self.points[:,0], 2.))).transpose()
    r = self.points[:,1]
    (a, b), _, _, _ = linalg.lstsq(A, r)
    self.a, self.b = (a,b)
#    return abs(2.*b)*self.frame_scale_inv
    return 2*math.sqrt(b*b / math.pow(1+a*a, 3))*self.frame_scale_inv
  def computeCurvatureByCircleFit(self):
    #WARING: fails for saddle points
    A = 2. * self.points[:,1:]
    b = np.power(self.points[:,0], 2.) + np.power(self.points[:,1], 2.)
    (r,) , _, _, _ = linalg.lstsq(A, b)
    self.circle_fit_r = r
    return 1./abs(r) * self.frame_scale_inv
  def computeCurvatureByCircleFit2(self):
    #WARING: fails for saddle points
    p0 = self.points[self.ic].copy()
    A = 2. * (self.points - p0)
    b = np.sum(np.power(self.points, 2.), axis=1)
    p , _, _, _ = linalg.lstsq(A, b)
    self.circle_fit_p = p
    return 1./linalg.norm(p) * self.frame_scale_inv


def estimateCurvature(points):
    curv = CurvEst(points)
    curv.estimateFrameByCovarianceTensor()
    curv.projectPointsToFrame()
    #return curv.computeCurvatureByCircleFit2()
    #return curv.computeCurvatureByCircleFit()
    return curv.computeCurvatureByParabolicFit()


def estimateCurvatures(components, points, filter_width):
  curvatures = np.zeros(len(points))
  for component in components:
    L = len(component)
    N = min(filter_width, (L-1) | 1)
    if N < 3: continue
    if component[0] == component[-1]:
      c = np.empty(L+filter_width-1, np.int)
      B = filter_width/2
      c[:B] = component[-B-1:-1]
      c[B:L+B] = component
      c[L+B:] = component[1:B+1]
      rng = B,L+B
    else:
      c = component
      rng = 1,L-1
#    print "ext index range = %s" % (c,)
#    print "component: %s" % component
    for i in xrange(*rng):
      W = min(filter_width/2, i, len(c)-i-1)
      indx = c[i-W:i+W+1]
#      print 'indx = %s' % indx
      pts = points[indx,:]
      #print "pts around index %i = %s" % (c[i],pts)
      curvatures[c[i]] = estimateCurvature(pts)
    if component[0] != component[-1]:
      curvatures[c[0]] = curvatures[c[1]]
      curvatures[c[-1]] = curvatures[c[-2]]
  return curvatures


def _test1():
  N = 5
  straightline = np.zeros((N, 3))
  straightline[:,0] = np.arange(N)
  estimateCurvature(straightline)

def _test2():
  N = 5
  angle = 90.
  radius = 0.1
  circle = np.zeros((N, 2))
  angles = (angle / 180. * math.pi /(N-1)) * np.arange(N)
  circle[:,0] = np.cos(angles)
  circle[:,1] = np.sin(angles)
  circle *= radius
  c = estimateCurvature(circle)
  print 'N=%i, c=%f, error=%f' % (N,c,abs(c - 1./radius)*radius)


def _plotbox(points):
  curv = CurvEst(points)
  
  p0, ew, v = curv.estimateFrameByCovarianceTensor()
  pyplot.plot(points[:,0], points[:,1], 'o-')
  px = p0 + v[0] * math.sqrt(ew[0])
  py = p0 + v[1] * math.sqrt(ew[1])
  pyplot.plot([p0[0], px[0]], [p0[1], px[1]])  
  pyplot.plot([p0[0], py[0]], [p0[1], py[1]])
  pyplot.show()

  curv.projectPointsToFrame()
  k = curv.computeCurvatureByParabolicFit()
  k2 = curv.computeCurvatureByCircleFit()
  k3 = curv.computeCurvatureByCircleFit2()
  print 'kpara=%f, kcirc=%f, kcirc2=%f' % (k, k2, k3)
  
#  print curv.__dict__  
  checkx = np.arange(-2, 2, 0.05)
  check = curv.b * np.power(checkx, 2.) + curv.a * checkx
  pts = curv.points
  pyplot.plot(pts[:,0], pts[:,1], 'x-')
  pyplot.plot(checkx, check)
  pyplot.plot(checkx, curv.a * checkx)
  r = curv.circle_fit_r
  pyplot.plot(checkx, r-np.sqrt(r*r - np.power(checkx, 2.)))
  angles = np.arange(0., math.pi*2., 0.1)
  r = linalg.norm(curv.circle_fit_p)
  cx = np.cos(angles) * r + curv.circle_fit_p[0]
  cy = np.sin(angles) * r + curv.circle_fit_p[1]
  pyplot.plot(cx, cy)
  pyplot.show()
  return k3

  
if __name__=='__main__':
  import matplotlib.pyplot as pyplot  

#  N = 41
#  a0 = -0.
#  r0 = 5.
#  a1 = 90.
#  r1 = 1.
#  angle_jitter = 0.3
#  radi_jitter = 0.01
#  angles = (((a1-a0)/(N-1)) * (np.arange(N) + np.random.randn(N)*angle_jitter) + a0)*(math.pi/180.)
#  radi = (np.arange(N)) * (r1-r0)/(N-1) + r0 + (np.random.randn(N)*radi_jitter*(r0+r1)*0.5)
#  circle = np.zeros((N, 2))
#  circle[:,0] = radi*np.cos(angles)
#  circle[:,1] = radi*np.sin(angles)
#  _plotbox(circle)

#  pyplot.plot(circle[:,0], circle[:,1],'o-')
#  pyplot.show()
#  if False:
#    indices = range(N)
#  else:
#    indices = range(N+1)
#    indices[-1] = indices[0]
#  curv = estimateCurvatures([indices], circle, 31)  
#  print curv
#  pyplot.plot(angles, curv)
#  pyplot.plot(angles, circle[:,0])
#  pyplot.plot(angles, circle[:,1])
#  pyplot.show()
  
  f = open("fuck.txt",'r')
  l = f.readlines()
  points = []
  for i, q in enumerate(l):
    if not q: continue
    try:
      a, b = q.split(' ')
    except:
      continue
    points.append((float(a), float(b)))
  points = np.asarray(points)
  
  indices = np.arange(len(points))
  curv = estimateCurvatures([indices], points, 31)  
#  pyplot.plot(indices, curv)
#  pyplot.plot(indices, points[:,0])
#  pyplot.plot(indices, points[:,1])
#  pyplot.show()

  for i in xrange(len(points)):
    print "%03i: %f" % (i, curv[i])

  ic = 632
  _plotbox(points[ic-15:ic+15+1,:])