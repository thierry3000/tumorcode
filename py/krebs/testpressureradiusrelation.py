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
from os.path import join, basename, dirname, splitext, commonprefix
import time
if __name__ == '__main__':
  sys.path.append(join(dirname(__file__),'..'))
import krebsutils
import h5py
import h5files
import numpy as np
import itertools
import extensions
import mpl_utils
import myutils
import math
import quantities

import matplotlib
import matplotlib.pyplot as pyplot

#from Marieb Human Anatomy & Physiology,
#blood pressure [mmHg]
#arterioles: 80 - 35 
#capillaries: 35 - 15
#venules:     15 - 10
#veins:       10 - 0

# thickness table from Marieb Human Anatomy & Physiology
# diameter [um],  thickness [um]
# 15000           1000              elastic artery
# 6000            1000              muscular artery
# 37              6                 arteriole
# 9               0.5               capillary
# 20              1                 venule
# 5000            500               vein
 
# from secomb or whereever it is from
def PressureRadiusRelation(rad):
  A2 = 18.
  A1 = 89.
  x0 = -21.
  dx = 16
  p = A2 + (A1-A2)/(1.+math.exp((rad-x0)/dx))
  p *= 0.1*1.33
  return p
 

# improved accuracy for large vessels,
# calibrated by eye measure according to Marieb
# has a steeper slope than my original curve (dx=10)
def PressureRadiusRelationHuman(rad):
  A2 = 0.
  A1 = 89.
  x0 = -17.
  dx = 10.
  p = A2 + (A1-A2)/(1.+math.exp((rad-x0)/dx))
  x0 = -21.
  dx = 2000
  pp = A2 + (A1-A2)/(1.+math.exp((rad-x0)/dx))
  p = 0.2*pp+0.8*p
  p *= 0.1*1.33
  return p
 

mariebTable = np.asarray([
  (-6000./2., 80., 90.),
  (-37./2., 35., 80.),
  (-9./2., 15., 35.),
  (20./2., 10., 15.),
  (5000./2., 10., 0.),
])
mariebTable[:,1:] *= quantities.mmHg.asNumber(quantities.kPa)

def varmap(x):
  #return np.sign(x)*np.log(np.abs(x)+1.)
  return x

if 0:
  fig, ax = pyplot.subplots(1,1)
  
  x = np.linspace(-6000., 6000., 6000)
  y = map(lambda x: PressureRadiusRelationHuman(x), x)
  ax.plot(varmap(x), y, color = 'r')
  
  y = map(lambda x: PressureRadiusRelation(x), x)
  ax.plot(varmap(x), y, color = 'g')
  
  #y = map(krebsutils.GetHealthyVesselWallThickness, x)
  #pyplot.plot(x, y)
  
  ax.vlines(varmap(mariebTable[:,0]), mariebTable[:,1], mariebTable[:,2], color = 'k')
  
  ax.set(xlim = (-100., 100.))
  pyplot.show()

if 1:
  Human = "RheologyForHuman"
  Rat   = "RheologyForRats"
  New   = "RheologySecomb2005"
  
  fig, ax = pyplot.subplots(2, 1)
  
  r = np.linspace(2., 200., 200)
  h = 0.45
  nu_h = [ krebsutils.CalcRelativeViscosity(r_, h, Human) for r_ in r ]
  nu_r = [ krebsutils.CalcRelativeViscosity(r_, h, Rat) for r_ in r ]
#  nu_ht = [ krebsutils.CalcRelativeViscosityByTable(r_, h, Human) for r_ in r ]
#  nu_rt = [ krebsutils.CalcRelativeViscosityByTable(r_, h, Rat) for r_ in r ]
  ax[0].plot(r, nu_h, label = 'human')
  ax[0].plot(r, nu_r, label = 'rat')
#  ax[0].plot(r, nu_ht, label = 'human table')
#  ax[0].plot(r, nu_rt, label = 'rat table')
  ax[0].legend()
  
  r = 500
  h = np.linspace(0., 0.6, 200)
  nu_h = [ krebsutils.CalcRelativeViscosity(r, h_, Human) for h_ in h]
  nu_r = [ krebsutils.CalcRelativeViscosity(r, h_, Rat) for h_ in h]
#  nu_ht = [ krebsutils.CalcRelativeViscosityByTable(r, h_, Human) for h_ in h]
#  nu_rt = [ krebsutils.CalcRelativeViscosityByTable(r, h_, Rat) for h_ in h]
  ax[1].plot(h, nu_h, label = 'human')
  ax[1].plot(h, nu_r, label = 'rat')
#  ax[1].plot(h, nu_ht, label = 'human table')
#  ax[1].plot(h, nu_rt, label = 'rat table')
  ax[1].legend()
  pyplot.show()
  
  rlist = np.concatenate((np.linspace(2, 20, 100, endpoint=False), np.linspace(20, 100.,  50)))
  hlist = [0.2, 0.4, 0.45, 0.6, 0.8]

  fig, ax = pyplot.subplots(1,1)  
  for h in hlist:
    eta_rel = map(lambda r: krebsutils.CalcRelativeViscosity(r, h, New), rlist)
    ax.plot(rlist, eta_rel)
  pyplot.show()
  
  #hd = np.linspace(0.3, 0.5, 20.)
#  r = np.linspace(2., 100., 200)
#  ht_human = [ krebsutils.CalcFahraeusEffect(0.45, r_, Human) for r_ in r]
#  ht_rat   = [ krebsutils.CalcFahraeusEffect(0.45, r_, Rat) for r_ in r]
#  fig, ax = pyplot.subplots(1,1)
#  ax.plot(r, ht_human, label='human')
#  ax.plot(r, ht_rat, label='rat')
#  ax.legend()
#  pyplot.show()