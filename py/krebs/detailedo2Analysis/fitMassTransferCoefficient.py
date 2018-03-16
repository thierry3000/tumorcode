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
import h5py
import numpy as np
import extensions # for hdf5 support in np.asarray
import krebsutils
import myutils
import posixpath
import math
from copy import copy, deepcopy
from collections import namedtuple
f2l = myutils.f2l

import matplotlib
import matplotlib.cm
import matplotlib.pyplot as pyplot
import mpl_utils
import matplotlib.ticker

from krebs import detailedo2
from krebs.detailedo2Analysis  import singleVesselParameterSets

import scipy.optimize

iteration = 0

fitCases = ['nair_uptake', 'nair_release'] + map(lambda i: 'moschandreou_case%02i' % i, range(6))
fitCases = map(lambda s: (s, singleVesselParameterSets.literature_cases[s]), fitCases)

moschandreou_s_vs_r = np.asarray([
  (6.,  0.848),
  (10., 0.891),
  (13.5,0.913),
  (27., 0.945),
  (40., 0.953),
  (50., 0.957),
]).transpose()

nair_release = np.asarray([
 (27.*0.5, 0.6),
]).transpose()

nair_uptake = np.asarray([
 (27.*0.5, 0.83),
]).transpose()
targetdata = np.concatenate((nair_uptake, nair_release, moschandreou_s_vs_r), axis=1)


def ProduceData(fitParameters, filename): 
  from krebs.analyzeGeneral import   DataBasicVessel, DataVesselSamples, DataVesselGlobal
  from krebs.detailedo2Analysis import DataDetailedPO2
  import krebs.detailedo2Analysis.singleVesselCases as singleVesselCases

  paramspo2Override = dict(
    massTransferCoefficientModelNumber = 1,
    conductivity_coeff1 = fitParameters[0],
    conductivity_coeff2 = fitParameters[1],
    conductivity_coeff3 = fitParameters[2],
  )

  f = h5files.open(filename,'a')

  krebsutils.set_num_threads(2)
  dataman = myutils.DataManager(20, [DataDetailedPO2(),DataBasicVessel(), DataVesselSamples(), DataVesselGlobal()])

  for k, params in fitCases:
    params = deepcopy(params)
    params.paramspo2.update(paramspo2Override)
    singleVesselCases.GenerateSingleCapillaryWPo2(dataman, f, k, 16, params)
  return f


def GetLastSatSample(group):
  dat = group['samples_and_fluxes/smpl_po2'][-1]
  parameters = myutils.hdf_read_dict_hierarchy(group['po2/parameters'])
  sat = detailedo2.PO2ToSaturation(dat, parameters)
  return group['parameters'].attrs['r'], sat[()]


def GetDataSeries(f):
  dat = []
  for i in xrange(6):
    g = f['moschandreou_case%02i' % i]
    r, s = GetLastSatSample(g)
    dat.append((r,s))
  dat_moschandreou = np.asarray(dat).transpose()
  g = f['nair_release']
  r, s = GetLastSatSample(g)
  nair_release = np.asarray([(r,s)]).transpose()
  g = f['nair_uptake']
  r, s = GetLastSatSample(g)
  nair_uptake = np.asarray([(r,s)]).transpose()
  return nair_uptake, nair_release, dat_moschandreou



def PlotDat((dat, coeffs), c):
  label = 'c1 = %f, c2 = %f, cg = %f' % (coeffs[0], coeffs[1], 0.)
  ax.plot(dat[0][0], dat[0][1], lw = 0, marker = '+', mew = 1, ms = 5, c = c)
  ax.plot(dat[1][0], dat[1][1], lw = 0, marker = 'x', mew = 1, ms = 5, c = c)
  ax.plot(dat[2][0], dat[2][1], lw = 0, marker = 'o', ms = 5., markeredgecolor = c, markerfacecolor = 'none', c = c, label = label)


if 0:
  nu_fit_coeffs = (6.07978304,  3.9202903,  0.01977186)

  #pmap = np.exp
  #pmapInv = np.log
  def pmap((p1, p2)):
    return math.exp(p1), math.exp(p2), nu_fit_coeffs[2]
  
  def pmapInv((p1, p2, p3))  :
    return math.log(p1), math.log(p2)
  
  def ObjectiveFunction(fitParameters):
    global iteration
    print '----- iteration %i ------' % iteration
    print 'ObjectiveFunction: mc=%s, c=%s' % (fitParameters, pmap(fitParameters))  
    fitParameters = pmap(fitParameters)  
    f = ProduceData(fitParameters, 'optimization_iteration%03i.h5' % iteration)
    d = GetDataSeries(f)
    d_ = np.concatenate(d, axis=1)[1]
    #pyplot.close('all')
    PlotDat((d, fitParameters), 'r')
    pyplot.draw()
    err = d_ - targetdata[1]
    #err[-2:] *= 20
    print 'currentdata = ', d_
    print 'currenterror = ', err
    print 'error = ',np.sum(err**2, axis=0)
    iteration += 1
    return err
  
  fig = pyplot.figure(figsize = (mpl_utils.a4size[0]*0.33, mpl_utils.a4size[0]*0.33))
  ax = fig.add_axes([0.2, 0.2, 0.75, 0.75])
  ax.plot(moschandreou_s_vs_r[0], moschandreou_s_vs_r[1], lw = 0, marker = 'o', ms = 5., c = 'k')
  ax.plot(nair_uptake[0], nair_uptake[1], lw = 2., marker = '+', ms = 5, c = 'k')
  ax.plot(nair_release[0], nair_release[1], lw = 2., marker = 'x', ms = 5, c = 'k')
  ax.set(xlabel = 'r [$\mu m$]', ylabel = 'S(x = 4 mm)')
  pyplot.show(block = False)
  
  initial_coeff, kwargs,  = pmapInv(nu_fit_coeffs), dict(epsfcn = 0.2, factor=5.)
  result = scipy.optimize.leastsq(ObjectiveFunction, initial_coeff, full_output=True, **kwargs)
  
  print result[3]
  coeff = pmap(result[0])
  result_f = ProduceData(coeff, 'optimization_iteration%03i.h5' % iteration)
  result_dat = GetDataSeries(result_f)
  PlotDat((result_dat, coeff), 'r')
  print 'RESULT = ',coeff
  pyplot.show()

else:
  #coeff = (0.0052*1.1, 16, 0.0078/math.exp(-4./16)*1.6)
  
  #coeff = (7.2, 4.0)
  
  #coeff = (21., 11.) # fit results based on initial guess (7.2, 4, 1.0)
  
  #coeff = (6.07978304,  3.9202903,  0.01977186)  # based on new fit to Nusselt Numbers from 19.11.2015
  #coeff = (5.240603622683324, 5.331639330237671, 0.01977186) # fitting the first two coefficients
  coeff = (7.99868309, 4.7419113, 0.)
  
  fig = pyplot.figure(figsize = (mpl_utils.a4size[0]*0.33, mpl_utils.a4size[0]*0.33))
  ax = fig.add_axes([0.2, 0.2, 0.75, 0.75])
  ax.plot(moschandreou_s_vs_r[0], moschandreou_s_vs_r[1], lw = 0, marker = 'o', ms = 5., markeredgecolor = 'k', markerfacecolor = 'none', c = 'k')
  ax.plot(nair_uptake[0], nair_uptake[1], lw = 2., marker = '+', ms = 5, c = 'k')
  ax.plot(nair_release[0], nair_release[1], lw = 2., marker = 'x', ms = 5, c = 'k')
  ax.set(xlabel = 'r [$\mu m$]', ylabel = 'S(x = 4 mm)')
  result_f = ProduceData(coeff, 'singelVesselCases.h5')
  result_dat = GetDataSeries(result_f)
  PlotDat((result_dat, coeff), 'r')
  pyplot.savefig('transascular_conductivity_datafit.svg')
  pyplot.show()


  