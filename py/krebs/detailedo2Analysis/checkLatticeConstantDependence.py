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
# -*- coding: utf-8 -*-
if __name__ == '__main__':
  import os.path, sys
  sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','..'))

import os,sys
import h5py
import h5files
import numpy as np
import extensions # for hdf5 support in np.asarray
import krebsutils
import myutils
import math
from pprint import pprint
from copy import copy, deepcopy

import matplotlib
import matplotlib.cm
import matplotlib.pyplot as pyplot
import mpl_utils

import scipy.optimize

from myutils import f2l

from krebs import detailedo2
from krebsjobs.parameters import parameterSetsO2

theTumorFile = 'tum-smallnetwork2.h5'
theTumorGroup = 'out0002'

krebsutils.set_num_threads(2)


def runDetailedO2(fn, pattern, o2params):
  print 'detailedo2 on %s / %s' % (fn, pattern)
  o2_refs = detailedo2.doit(fn, pattern, (o2params, o2params['name']))
  h5files.closeall() # just to be sure
  return o2_refs[0]

resultFileNames = []
latticeConstants = [ 20, 25, 30, 35, 40, 50 ]
baseParams = parameterSetsO2.breastv3
for lc in latticeConstants:
  params = deepcopy(baseParams)
  params['grid_lattice_const'] = float(lc)
  params['name'] = 'lattice-const-%02i' % lc
  params['calcflow']['includePhaseSeparationEffect'] = False
  tmp = runDetailedO2(theTumorFile, theTumorGroup, params)
  resultFileNames.append(tmp)

vesselpo2list = []
tissuepo2list = []
for fn, path in resultFileNames:
  f = h5files.open(fn)
  po2group = f[path]
  po2vessels = np.asarray(po2group['po2vessels'])
  vesselpo2list.append(po2vessels)
  po2field = np.asarray(po2group['po2field'])
  tmp = np.average(po2field)
  tissuepo2list.append(tmp)

h5files.closeall()

tmp = np.linalg.norm(vesselpo2list[0])
for i, po2vessels in enumerate(vesselpo2list[1:]):
  vesselpo2list[i+1] =  np.linalg.norm(po2vessels-vesselpo2list[0])/tmp
vesselpo2list[0] = 0

vesselpo2list = np.asarray(vesselpo2list, dtype = np.float)
tissuepo2list = np.asarray(tissuepo2list, dtype = np.float)
latticeConstants = np.asarray(latticeConstants, dtype = np.float)

tmp = tissuepo2list[0]
tissuepo2list = np.abs(tissuepo2list - tmp)/tmp

def powerfunc(x, p, q):
  return (np.power(x/q[0], p[0])-1.)*p[1]

def powerfit(x, y, pinit, q):
  def objectivefunction1(p):
    return (y - powerfunc(x, p, q))
  (resultcoeff, _) = scipy.optimize.leastsq(objectivefunction1, pinit)
  return resultcoeff, lambda x: powerfunc(x, resultcoeff, q)

po2vesselsCoeff, po2VesselFunc = powerfit(latticeConstants, vesselpo2list, (1., 2.), (latticeConstants[0],))
print po2vesselsCoeff

po2tissueCoeff, po2TissueFunc = powerfit(latticeConstants, tissuepo2list, (1.,2.), (latticeConstants[0],))
print po2tissueCoeff

x = np.linspace(0., 50., 100)

fig, axes = pyplot.subplots(2, 1, figsize = (mpl_utils.a4size[0]*0.4, mpl_utils.a4size[1]*0.3))
axes[0].plot(latticeConstants, vesselpo2list, label = r'$||P(h) - P^*||/||P^*||$', lw = 0, marker = 'x', color = 'k')
axes[0].plot(x, po2VesselFunc(x), label = r'$((h/20 \mu m)^{%s} - 1)*%s$' % tuple(map(f2l, po2vesselsCoeff)), color = 'k')
axes[0].legend(loc = mpl_utils.loc.lower_right)
axes[0].set(xlabel = r'$h$ [$\mu m$]', ylabel = r'mmHg')
axes[0].set(ylim = (-0.2, 0.2))
axes[1].plot(latticeConstants, tissuepo2list, label = r'$|<P_t(h)> - <P_t^*>|/|<P_t^*>|$', lw = 0, marker = 'x', color = 'k')
axes[1].plot(x, po2TissueFunc(x), label = r'$((h/20 \mu m)^{%s} - 1)*%s$' % tuple(map(f2l, po2tissueCoeff)), color = 'k')
axes[1].legend(loc = mpl_utils.loc.lower_right)
axes[1].set(xlabel = r'$h$ [$\mu m$]', ylabel = r'mmHg')
pyplot.tight_layout()
fig.savefig('fubar.pdf')
#pyplot.show()