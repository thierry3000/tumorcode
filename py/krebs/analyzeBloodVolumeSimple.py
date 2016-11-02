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
#!/usr/bin/env python2
# -*- coding: utf-8 -*-
if __name__ == '__main__':
  import os.path, sys
  sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'..'))
  
import os,sys
import h5py
import numpy as np
from scipy.spatial.distance import pdist
import extensions # for asarray with h5py support
import krebsutils
import posixpath
import math
import glob
from collections import defaultdict
from pprint import pprint

import matplotlib
import identifycluster
if (identifycluster.getname()=='snowden' or identifycluster.getname()=='durga'):
  matplotlib.use('agg')



#from plotBulkTissue import commonOutputName, OutputWriter, add_sizebar, colorbar, contour, imslice, imshow, ColorMaps, subplots_abs_mm, fig_numbering, calc_phi_vessels, calc_distmap
import myutils

from plotBulkTissue2d import DataTumorTissueSingle
from analyzeBloodFlow import DataTumorBloodFlow, obtain_averaged_blood_flow
from analyzeGeneral import RemoveArteriovenousFlagsFromCapillaries,totalLdVolume


## ---------------- ------- ----------------------------------------
## compute the mvd/volume density globally or locally in radial bins
## ---------------- ------- ----------------------------------------
def cylinderCollectionVolumeDensity(vesselgroup):
  vessels = krebsutils.read_vesselgraph(vesselgroup, ['flags', 'radius', 'length'])
  flags   = RemoveArteriovenousFlagsFromCapillaries(vessels['flags'])
  mask1 = myutils.bbitwise_and(flags, krebsutils.CIRCULATED)
  totalvol = totalLdVolume(vesselgroup)
  def compute(flagMask):
    if flagMask:
      mask = mask1 & myutils.bbitwise_and(flags, flagMask)
    else:
      mask = mask1
    length = np.asarray(vessels['length'][mask], dtype = np.float64)
    radius = np.asarray(vessels['radius'][mask], dtype = np.float64)
    vol = math.pi * np.power(radius, 2.) * length
    vol = np.sum(vol)
    return vol/totalvol
  return compute(0), compute(krebsutils.ARTERY), compute(krebsutils.VEIN), compute(krebsutils.CAPILLARY)

def cylinderCollectionLineDensity(vesselgroup):
  vessels = krebsutils.read_vesselgraph(vesselgroup, ['flags', 'length'])
  flags   = RemoveArteriovenousFlagsFromCapillaries(vessels['flags'])
  mask1 = myutils.bbitwise_and(flags, krebsutils.CIRCULATED)
  totalvol = totalLdVolume(vesselgroup)
  def compute(flagMask):
    if flagMask:
      mask = mask1 & myutils.bbitwise_and(flags, flagMask)
    else:
      mask = mask1
    length = np.asarray(vessels['length'][mask], dtype = np.float64)
    total = np.sum(length)
    return total/totalvol
  return compute(0), compute(krebsutils.ARTERY), compute(krebsutils.VEIN), compute(krebsutils.CAPILLARY)


## ---------------- ------- -----------------------
## plotting
## ---------------- ------- -----------------------
def plotShit(filenames):
  data = []
  for fn in filenames:
    with h5py.File(fn, 'r') as f:
      l = float(f['vessels/lattice'].attrs['SCALE'])
      rv = cylinderCollectionVolumeDensity(f['vessels'])
      mvd = cylinderCollectionLineDensity(f['vessels'])
      data.append((l, rv, mvd))
  data = np.asarray(sorted(data, key = lambda x: x[0])).transpose()

  fig, ax = pyplot.subplots(1,1)
  ax.plot(data[0], data[1])
  fig.savefig('shit.png')


def printBloodFlowSimple(filenames, group):
  files = [ h5py.File(fn, 'r') for fn in filenames ]
  fnmeasure = 'analyzeBloodVolumeMeasurements.h5'
  dataman = myutils.DataManager(100, [DataBloodVolume(fnmeasure), DataTumorTissueSingle(), DataTumorBloodFlow()])

  bv = obtain_averaged_blood_flow(dataman, files, group)
  # unit is 1/s, converted to 1/min -> x60
  def printit(id, label, unit):
    s = r'%s = %s +/- %s [%s]' % (
      label, myutils.f2s(bv[id][0]*60.), myutils.f2s(bv[id][1]*60.), unit)
    print s
  if 'tumor_flow' in bv:
    printit('tumor_flow', 'BF_{tumor}', 'um^3 Blood / min'),
    printit('tumor_flow_p_volume', 'rBF_{tumor}', 'ml Blood / (ml Tissue min)'),
    print 'tumor_volume %e' % bv['tumor_volume'][0]
  printit('total_flow', 'BF_{total}', 'um^3 Blood / min'),
  printit('total_flow_p_volume', 'rBF_{total}', 'ml Blood / (ml Tissue min)')
  print 'total_volume %e' % bv['total_volume'][0]



if __name__ == '__main__':
  group = sys.argv[1]
  dat = []
  for fn in sys.argv[2:]:
    with h5py.File(fn, 'r') as f:
      rv, rv_a, rv_v, rv_c = cylinderCollectionVolumeDensity(f[group])
      mvd, mvd_a, mvd_v, mvd_c = cylinderCollectionLineDensity(f[group])
      dat.append((rv, rv_a, rv_v, rv_c, mvd, mvd_a, mvd_v, mvd_c))
  dat = np.asarray(dat).transpose()
  def printstuff(name, a, mult):
    print '%s = %f +/- %f' % (name, np.average(a)*mult, np.std(a)*mult)
  rbv, a, v, c = dat[:4]
  printstuff('rBV'  , rbv, 1.)
  printstuff('rBV_a', a, 1.)
  printstuff('rBV_v', v, 1.)
  printstuff('rBV_c', c, 1.)
  fv = v/(a+v)
  printstuff('f_v', fv, 1.)
  printstuff('mvd', dat[4], 1e6)
  printstuff('mvd_a', dat[5], 1e6)
  printstuff('mvd_v', dat[6], 1e6)
  printstuff('mvd_c', dat[7], 1e6)

  #printBloodFlowSimple(sys.argv[2:], group)

# experimental results:
#  rBV(normal) = 0.015 (haut) ??????
#  rBV(normal) = 0.01 (burst)
# von döme
#  mvd(normal) = 80 [#/mm^2] (haut)
#  intercap-dist ~ 0.19
#  rBV = 0.006 (5 mm capillare)
#  rBV = 0.025 (10 mm capillare)

# angenommen intercap-dist = 0.08, rcap = 10 um -> rBV = 0.15
# intercap-dist = 0.08, rcap = 5 um -> 0.037

# mit neuesten ifp_paper bäumen (intercap-dist = 80 um)
# rBV = 0.08
# mvd = 470

# mit quad gitter und intercap-dist = 110 um, rcap = 8 um
#  rBV = 0.052
#  mvd = 258 [1/mm^2]
# was den erwartungen entspricht

#$rBF_{normal}$ & 0.06  \\
#$rBF_{tumor}$ & 0.32 \\
