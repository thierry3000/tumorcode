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
  sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','..'))
  
import os,sys
import numpy as np
#from .. import extensions # for hdf5 support in np.asarray
import extensions
import krebsutils
import myutils
import posixpath
import math
from copy import copy, deepcopy
from collections import namedtuple
f2l = myutils.f2l
#import krebs.detailedo2Analysis

import matplotlib
import matplotlib.cm
import matplotlib.pyplot as pyplot
import mpl_utils
import matplotlib.ticker

from krebs.detailedo2Analysis import singleVesselParameterSets

import scipy.optimize


f2l_ = lambda x: f2l(x, prec=2)


def MTC2Nusselt(mtc, r, params):
  return 2.0*r / (params['D_plasma']*params['solubility_plasma']) * mtc

def MTC2Sherwood(mtc, r, params):
  return 2.0*r / (params['D_tissue']*params['solubility_tissue']) * mtc

def SherwoodToNusselt(sh, r, params):
  return sh*((params['D_tissue']*params['solubility_tissue'])/(params['D_plasma']*params['solubility_plasma']))

def NuToMTC(nu, r, params):
  return (params['D_tissue']*params['solubility_tissue']) / (2.0*r) * nu;

#  paramspo2Override = dict(
#    massTransferCoefficientModelNumber = 1,
#    conductivity_coeff1 = 7.2,
#    conductivity_coeff2 = 4.7
#  )

#hellums_1996, safaeian thesis:
nusselt_number_vs_r = np.asarray([
  (3.5*0.5, 1.2),
  (5.4*0.5, 1.7),
  (8.1*0.5, 2.5),
  (10*0.5, 2.25),
  (15*0.5, 2.1),
  (20.0*0.5, 3.25),
  (29.0*0.5, 4.0),
  (50.0*0.5, 4.6),
  (100.0*0.5, 4.8),
])

#moschandreou_2011
sherwood_numbers_tissue = np.asarray([
  (6., 0.6),
  (13.5, 0.8),
  (27., 1.1),
  (50., 1.53),
])

sherwood_numbers = np.asarray([
  (6., 2.64),
  (10., 3.42),
  (13.5, 3.5),
  (27., 4.375),
  (40., 4.625),
  (50., 4.75),
])

sherwood_numbers_extra = np.asarray([
  (1500.0, 200.0),
  (3000.0, 50.0)
])

#D (um)  K (mmHg.um.s/um3 O2), i think it is MTC = 1/(K 2 pi r)
originalSecombDataOfK = '''
    4           2.45
    5           2.00
    6           1.74
    8           1.44
   10           1.27
   12           1.16
   14           1.08
   16           1.02
   18           0.974
   20           0.935
   22           0.895
   24           0.853
   28           0.815
   32           0.779
   40           0.750
   50           0.730
   54		    0.747
'''
originalSecombDataOfK = np.asarray(map(float, filter(lambda s: s, map(str.strip, originalSecombDataOfK.split(' ')))))
originalSecombDataOfK = originalSecombDataOfK.reshape(-1,2)

secombMTC = originalSecombDataOfK.copy()
d = secombMTC[:,0]
K = secombMTC[:,1]
secombMTC[:,1] = np.power(math.pi * K * d, -1.)
secombMTC[:,0] = 0.5 * d
del d, K

print 'secomb MTC bounds: (%f - %f um) = %f - %f' % (secombMTC[:,0].min(), secombMTC[:,0].max(), secombMTC[:,1].min() ,secombMTC[:,1].max())


def BadFitNusseltNumber(r):  # used in single vessel experiments
  #coeff = (21., 11.)  # result of the data fit to single vessels
  p1, p2 = (7.2, 8.0)
  return p2*(1.0 - np.exp(-r/p1))


paramspo2base = deepcopy(singleVesselParameterSets.basecase.paramspo2)

# examine fit curves for nusselt numbers and mtc

#r_arr = np.linspace(1.5, 3000., 1000)
r_arr = np.logspace(0., math.log10(4000.0), 200)

nuFromMoschandreou = sherwood_numbers.copy()
nuFromMoschandreou[:,1] = SherwoodToNusselt(sherwood_numbers[:,1], sherwood_numbers[:,0], paramspo2base)

nuFromSecomb = secombMTC.copy()
nuFromSecomb[:,1] = MTC2Nusselt(secombMTC[:,1], secombMTC[:,0], paramspo2base) 

nuExtraLarge = sherwood_numbers_extra.copy()
nuExtraLarge[:,1] = SherwoodToNusselt(sherwood_numbers_extra[:,1], sherwood_numbers_extra[:,0], paramspo2base)

print 'nusselt Extra Large = ', nuExtraLarge

actuallyUsedCoeffsForActuallyUsedMTC = (0.0052*1.1, 16, 0.0078/math.exp(-4./16)*1.6) # = (0.00572, 16, 0.016024637200263012)
print 'actuallyUsedCoeffsForActuallyUsedMTC =', actuallyUsedCoeffsForActuallyUsedMTC

def actuallyUsedMTC(r):
  p1, p2 ,p3  = actuallyUsedCoeffsForActuallyUsedMTC
  return np.exp(-r/p2)*p3 + p1

def actuallyUsedNu(r):
  return MTC2Nusselt(actuallyUsedMTC(r), r, paramspo2base)

def fitfunction(r, coeff_):
  p1, p2, p3 = coeff_
  return r*(np.exp(-r/p2)*p3 + p1)

#fitDataX = np.concatenate((sherwood_numbers[:,0], nusselt_number_vs_r[:,0], sherwood_numbers_extra[:,0]))
#fitDataY = np.concatenate((nuFromMoschandreou, nusselt_number_vs_r[:,1], nuFromSherwoodExtra))
#fitDataW = np.ones_like(fitDataX)
#fitDataW[len(sherwood_numbers)+4] = 0.0 # the outlier
#def objectivefunction1(p):
#  return fitDataW*(fitDataY - fitfunction(fitDataX, p))
#(fitcoeff1, _) = scipy.optimize.leastsq(objectivefunction1, (0.,20.,0.1))
#overRFit1 = fitfunction(r_arr, fitcoeff1)
#fitlabel1 = r'r\,(exp(-r/{1})*{2} + {0})'.format(*map(f2l_,fitcoeff1))

fitData = np.concatenate((nusselt_number_vs_r, nuFromMoschandreou), axis = 0)
#fitData = nuFromSecomb
fitDataX = fitData[:,0]
fitDataY = fitData[:,1]

fitDataW = np.ones_like(fitDataX)
fitDataW[len(sherwood_numbers)+4] = 0.0 # the outlier

if 0:  # this is a better version of the function below.  It takes into account that Nu should increase linearly with the radius for very large radii
    def fitfunction2(r, coeff_):
      p1, p2, p3 = coeff_
      return p2*(1.0 - np.exp(-r/p1)) + p3 * r
    initialCoeff = (10.0, 3.0, 0.)
    def objectivefunction2(coeff_):
      return fitDataW*(fitDataY - fitfunction2(fitDataX, coeff_))
    (nusselt_fitcoeff, stuff) = scipy.optimize.leastsq(objectivefunction2, initialCoeff)
    fitlabel2 = '{1}\,(1 - exp(-r/{0})) + {2} r'.format(*map(f2l_, nusselt_fitcoeff))
    print 'nusselt_fitcoeff for %s = %s' % (fitlabel2, nusselt_fitcoeff)
else:
    def fitfunction2(r, coeff_):
      p1, p2 = coeff_
      return p2*(1.0 - np.exp(-r/p1))
    initialCoeff = (10.0, 3.0)
    def objectivefunction2(coeff_):
      return fitDataW*(fitDataY - fitfunction2(fitDataX, coeff_))
    (nusselt_fitcoeff, stuff) = scipy.optimize.leastsq(objectivefunction2, initialCoeff)
    fitlabel2 = '{1}\,(1 - exp(-r/{0}))'.format(*map(f2l_, nusselt_fitcoeff))
    print 'nusselt_fitcoeff for %s = %s' % (fitlabel2, nusselt_fitcoeff)
    

print 'Nu @r = 3000 = ', fitfunction2(3000., nusselt_fitcoeff)

nuFromMtc_arr = fitfunction2(r_arr, nusselt_fitcoeff)

mask = (r_arr<60.)  & (r_arr>2.)
a = NuToMTC(nuFromMtc_arr[mask], r_arr[mask], paramspo2base)
r = r_arr[mask]
print 'used in Single Vessel Plots: MTC(r = %f) = %f,  MTC(r = %f) = %f' % (r[0], a[0], r[-1], a[-1])

actuallUsedNu_arr = actuallyUsedNu(r_arr)

if 1:
  fig, axes = pyplot.subplots(2, 1, figsize = (mpl_utils.a4size[0]*0.4, mpl_utils.a4size[0]*0.6))
  scale = lambda r, values:   NuToMTC(values, r, paramspo2base)
  ax = axes[0]
  #ax.plot(r_arr, scale(r_arr, actuallUsedNu_arr), label = '$MTC^{**}/r$', color='k')
  #ax.plot(r_arr, scale(r_arr, nuFromMtc_arr), label = '$MTC^*/r$', color = 'r')
  ax.plot(r_arr, scale(r_arr, fitfunction2(r_arr, nusselt_fitcoeff)), label = '$MTC(Nu = '+fitlabel2+')$', color = 'b')
  ax.plot(nusselt_number_vs_r[:,0], scale(nusselt_number_vs_r[:,0], nusselt_number_vs_r[:,1]), label = '$MTC$  (hellums)', lw = 0, ms = 5., marker = 'x', markeredgecolor = 'k', markerfacecolor = 'none', c = 'k')
  ax.plot(nuFromMoschandreou[:,0], scale(nuFromMoschandreou[:,0],nuFromMoschandreou[:,1]), label = '$MTC$ (moschandreou)', lw = 0, ms = 5., marker = 'o', markeredgecolor = 'k', markerfacecolor = 'none', c = 'k')
  #ax.plot(nuFromSecomb[:,0], scale(nuFromSecomb[:,0], nuFromSecomb[:,1]), label = '$MTC$ (secomb)', lw = 0, ms = 5., marker = 'd', markeredgecolor = 'k', markerfacecolor = 'none', c = 'k')
  #ax.plot(sherwood_numbers_extra[:,0], scale(sherwood_numbers_extra[:,0])*nuFromSherwoodExtra, label = '$nu/r$ (artery)', lw = 0, ms = 5., marker = '+', markeredgecolor = 'k', markerfacecolor = 'none', c = 'k')
  ax.set(xscale='log', 
         xlim = (1.,200))
  ax.legend()

  ax = axes[1]
  #ax.plot(r_arr, actuallUsedNu_arr,  label = '$nu^{**}$', color='k')
  #ax.plot(r_arr, nuFromMtc_arr, label = '$nu^*$', color = 'r')
  ax.plot(r_arr, fitfunction2(r_arr, nusselt_fitcoeff), label = '$'+fitlabel2 + '$', color = 'b')
  ax.plot(nusselt_number_vs_r[:,0], nusselt_number_vs_r[:,1], label = '$Nu$ (hellums)', lw = 0, ms = 5., marker = 'x', markeredgecolor = 'k', markerfacecolor = 'none', c = 'k')
  ax.plot(nuFromMoschandreou[:,0], nuFromMoschandreou[:,1], label = '$Nu$ (moschandreou)', lw = 0, ms = 5., marker = 'o', markeredgecolor = 'k', markerfacecolor = 'none', c = 'k')
  #ax.plot(nuFromSecomb[:,0], nuFromSecomb[:,1], label = '$Nu$ (secomb)', lw = 0, ms = 5., marker = 'd', markeredgecolor = 'k', markerfacecolor = 'none', c = 'k')
  #ax.plot(sherwood_numbers_extra[:,0], nuFromSherwoodExtra, label = '$nu$ (artery)', lw = 0, ms = 5., marker = '+', markeredgecolor = 'k', markerfacecolor = 'none', c = 'k')
  ax.set(xscale='log',
         xlim = (1.,200))
  ax.set(ylim=(0., 6.)) 
         #yscale = 'log')
  ax.set(xlabel = 'r [$\mu m$]')
  ax.legend(loc = mpl_utils.loc.upper_left)

  fig.savefig('transvascular_conductivity_constant.svg')
  #pyplot.show()
  
else:
  fig, ax = pyplot.subplots(1,1, figsize = (mpl_utils.a4size[0]*0.8, mpl_utils.a4size[0]*0.4))
  r = np.linspace(1., 200., 200)
  mtc = fitfunction(r, coeff)/r
  print np.amin(mtc), np.amax(mtc)
  ax.plot(r, mtc, label = 'mtc', color = 'k')
  pyplot.show()
