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
  sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../krebs'))
  
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
import mpl_utils
import identifycluster
if (identifycluster.getname()=='snowden' or identifycluster.getname()=='durga'):
  matplotlib.use('agg')



#from plotBulkTissue import commonOutputName, OutputWriter, add_sizebar, colorbar, contour, imslice, imshow, ColorMaps, subplots_abs_mm, fig_numbering, calc_phi_vessels, calc_distmap
import myutils

#from plotBulkTissue2d import DataTumorTissueSingle
from analyzeBloodFlow import DataTumorBloodFlow, obtain_averaged_blood_flow
from analyzeGeneral import RemoveArteriovenousFlagsFromCapillaries,totalLdVolume


## ---------------- ------- ----------------------------------------
## compute the mvd/volume density globally or locally in radial bins
## ---------------- ------- ----------------------------------------
def cylinderCollectionVolumeDensity(vesselgroup):
  vessels = krebsutils.read_vesselgraph(vesselgroup, ['flags', 'radius', 'length'])
  flags   = RemoveArteriovenousFlagsFromCapillaries(vessels['flags'])
  flags = flags[:,0]
  mask1 = myutils.bbitwise_and(flags, krebsutils.CIRCULATED)
  totalvol = totalLdVolume(vesselgroup)
  def compute(flagMask):
    if flagMask:
      mask = mask1 & myutils.bbitwise_and(flags, flagMask)
    else:
      mask = mask1
    length = np.asarray(vessels['length'][mask], dtype = np.float64)
    radius = np.asarray(vessels['radius'][mask], dtype = np.float64)
    radius = radius[:,0]
    vol = math.pi * np.power(radius, 2.) * length
    vol = np.sum(vol)
    return vol/totalvol
  return compute(0), compute(krebsutils.ARTERY), compute(krebsutils.VEIN), compute(krebsutils.CAPILLARY)

def cylinderCollectionLineDensity(vesselgroup):
  vessels = krebsutils.read_vesselgraph(vesselgroup, ['flags', 'length'])
  flags   = RemoveArteriovenousFlagsFromCapillaries(vessels['flags'])
  flags = flags[:,0]
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

def surface2Volume(vesselgroup):
  vessels = krebsutils.read_vesselgraph(vesselgroup, ['flags', 'radius', 'length'])
  flags   = RemoveArteriovenousFlagsFromCapillaries(vessels['flags'])
  flags = flags[:,0]
  mask1 = myutils.bbitwise_and(flags, krebsutils.CIRCULATED)
  totalvol = totalLdVolume(vesselgroup)
  def compute(flagMask):
    if flagMask:
      mask = mask1 & myutils.bbitwise_and(flags, flagMask)
    else:
      mask = mask1
    length = np.asarray(vessels['length'][mask], dtype = np.float64)
    radius = np.asarray(vessels['radius'][mask], dtype = np.float64)
    radius = radius[:,0]
    vol = math.pi * np.power(radius, 2.) * length
    surface = 2* math.pi* radius* length
    #surface = np.sum(surface)
    #vol = np.sum(vol)
    s2v = surface/vol
    mystd = np.std(s2v)
    myavg = np.average(s2v)
    myrel = mystd/myavg
    print("spread: %f" % myrel)
    return np.average(surface/vol)
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

def boxplotFromData_rBV(data, pp, data2=None):
  print('do box')
  fig2, ax2 = matplotlib.pyplot.subplots(1,1)
  #ax2.boxplot([data[0,:], data[1,:],data[2,:], data[3,:]])
  violin1 = ax2.violinplot([data[0,:], data[1,:],data[2,:], data[3,:]], positions=[0,1,2,3],showmeans=False,showmedians=True)
  
  if not data2 is None:
    violin2 = ax2.violinplot([data2[0,:], data2[1,:],data2[2,:], data2[3,:]],positions=[0.5,1.5,2.5,3.5],showmeans=False,showmedians=True)
    ax2.set_xticks([0.25, 1.25, 2.25,3.25])
  else:
    ax2.set_xticks([0., 1., 2.,3.])
  ax2.set_xticklabels([r'$rBV$ in $\%$', r'artery', r'vein',r'capillary'])
  
  
  aColor = violin1['cmedians'].get_color()
  aColor = aColor[0,0:3]
  blue_patch = matplotlib.patches.Patch(color=aColor)
  if not data2 is None:
    aColor = violin2['cmedians'].get_color()
    aColor = aColor[0,0:3]
    red_patch = matplotlib.patches.Patch(color=aColor)
    label = ['tumorCode', 'adaption']
    fake_handels = [blue_patch, red_patch]
    ax2.legend(fake_handels, label)
#  else:
#    label = ['adaption']
#    fake_handels = [blue_patch]
#    ax2.legend(fake_handels, label)
  
  #ax2.legend(['bl','dlf'])
  pp.savefig(fig2)
  
  
def boxplotFromData_s2v(data, pp, data2=None):
  print('do box')
  fig2, ax2 = matplotlib.pyplot.subplots(1,1)
  #ax2.boxplot([data[0,:], data[1,:],data[2,:], data[3,:]])
  violin1 = ax2.violinplot([data[8,:], data[9,:],data[10,:], data[11,:]], positions=[0,1,2,3],showmeans=False,showmedians=True)
  
  if not data2 is None:
    violin2 = ax2.violinplot([data2[8,:], data2[9,:],data2[10,:], data2[11,:]],positions=[0.5,1.5,2.5,3.5],showmeans=False,showmedians=True)
    ax2.set_xticks([0.25, 1.25, 2.25,3.25])
    aColor = violin1['cmedians'].get_color()
    aColor = aColor[0,0:3]
    blue_patch = matplotlib.patches.Patch(color=aColor)
    aColor = violin2['cmedians'].get_color()
    aColor = aColor[0,0:3]
    red_patch = matplotlib.patches.Patch(color=aColor)
    label = ['tumorCode', 'adaption']
    fake_handels = [blue_patch, red_patch]
    ax2.legend(fake_handels, label, loc='upper center')
  else:
    ax2.set_xticks([0., 1., 2.,3.])
    aColor = violin1['cmedians'].get_color()
    aColor = aColor[0,0:3]
    blue_patch = matplotlib.patches.Patch(color=aColor)
    label = ['adaption']
#    fake_handels = [blue_patch]
#    ax2.legend(fake_handels, label, loc='upper center')
  ax2.set_xticklabels([r'all', r'artery', r'vein',r'capillary'])
  
  
  
  #ax2.legend(['bl','dlf'])
  pp.savefig(fig2)

def getDataFromFiles(filenames, groupname):
  dat = []
  for fn in filenames:
    with h5py.File(fn, 'r') as f:
      print('open file: %s at %s' %(fn,groupname))
      rv, rv_a, rv_v, rv_c = cylinderCollectionVolumeDensity(f[groupname])
      mvd, mvd_a, mvd_v, mvd_c = cylinderCollectionLineDensity(f[groupname])
      s2v, s2v_a, s2v_v, s2v_c = surface2Volume(f[groupname])
      dat.append((rv, rv_a, rv_v, rv_c, mvd, mvd_a, mvd_v, mvd_c,s2v,s2v_a,s2v_v,s2v_c))
  dat = np.asarray(dat).transpose()
  
  ### rBV to percentage
  dat[0:4,:] = 100*dat[0:4,:]
  return dat

if __name__ == '__main__':
  
  import argparse
  parser = argparse.ArgumentParser(description='Plot/ Analyze rBV surface to volume.')
  parser.add_argument('grp_pattern1',help='Where to find the vessel group in the file')    
  parser.add_argument('--grp_pattern2',default=None, help='Where to find the vessel group in the file') 
  parser.add_argument('vesselFileNames1', nargs='*', type=argparse.FileType('r'), default=sys.stdin, help='Vessel file to calculate')  
  
  #parser.add_argument('--vesselFileNames2', nargs='*', type=argparse.FileType('r'), default=sys.stdin, help='Vessel file to calculate')  
  #parser.add_argument('vesselFileNames2', nargs='+', type=argparse.FileType('r'), default=sys.stdin,help='Vessel file to calculate')  
  
    
  goodArguments, otherArguments = parser.parse_known_args()
  
  n_files = len(goodArguments.vesselFileNames1)
  print('found %i files' % n_files)
  
  if not goodArguments.grp_pattern2 is None:
    #this means we have 2 file per simulation 
    # adapted and non adapted
    n_files=n_files/2
    
  try:
    dirs = set()
    for fn in goodArguments.vesselFileNames1[0:n_files]:
      if not os.path.isfile(fn.name):
        raise AssertionError('The file %s is not present!'%fn)
      with h5py.File(fn.name, 'r') as f:
        d = myutils.walkh5(f, goodArguments.grp_pattern1)
        if not len(d)>0:
          raise AssertionError('pattern "%s" not found in "%s"!' % (goodArguments.grp_pattern1, fn))
        else:
          dirs = set.union(dirs,d)
  except Exception, e:
    print e.message
    sys.exit(-1)
  filenames1=[]
  for fn in goodArguments.vesselFileNames1[0:n_files]:
    filenames1.append(fn.name)
  print(filenames1)
  
  if not goodArguments.grp_pattern2 is None:
    try:
      dirs = set()
      for fn in goodArguments.vesselFileNames1[n_files:]:
        if not os.path.isfile(fn.name):
          raise AssertionError('The file %s is not present!'%fn)
        with h5py.File(fn.name, 'r') as f:
          d = myutils.walkh5(f, goodArguments.grp_pattern2)
          if not len(d)>0:
            raise AssertionError('pattern "%s" not found in "%s"!' % (goodArguments.grp_pattern2, fn))
          else:
            dirs = set.union(dirs,d)
    except Exception, e:
      print e.message
      sys.exit(-1)
    filenames2=[]
    for fn in goodArguments.vesselFileNames1[n_files:]:
      filenames2.append(fn.name)
    print(filenames2)
  #group = sys.argv[-1]
  #filenames = sys.argv[1:-1]
  #print(filenames)
  
  
  dat = getDataFromFiles(filenames1, goodArguments.grp_pattern1)
  if not goodArguments.grp_pattern2 is None:
    dat2 = getDataFromFiles(filenames2, goodArguments.grp_pattern2)
    outfilename='rBV_for_file_%s' % os.path.basename(filenames1[0])
    with mpl_utils.PdfWriter(outfilename + '.pdf') as pdfpages:  
      boxplotFromData_rBV(dat, pdfpages, data2 = dat2)
    outfilename='s2v_for_file_%s' % os.path.basename(filenames1[0])
    with mpl_utils.PdfWriter(outfilename + '.pdf') as pdfpages:  
      boxplotFromData_s2v(dat, pdfpages, data2 = dat2)
  else:
    #case of no adaption, or only adaption
    outfilename='rBV_for_files'
    with mpl_utils.PdfWriter(outfilename + '.pdf') as pdfpages:  
      boxplotFromData_rBV(dat, pdfpages)
    outfilename='s2v_for_files'
    with mpl_utils.PdfWriter(outfilename + '.pdf') as pdfpages:  
      boxplotFromData_s2v(dat, pdfpages)
    
  
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
  printstuff('s2v', dat[8], 1)
  printstuff('s2v_a', dat[9], 1)
  printstuff('s2v_v', dat[10], 1)
  printstuff('s2v_c', dat[11], 1)

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
