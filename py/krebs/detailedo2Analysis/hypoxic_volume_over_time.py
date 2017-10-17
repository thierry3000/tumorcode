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

"""
T.F. in coorperation with Gero Wedemann, Bertin Hoffmann, Uni Stralsund

this file estimates hypoxic volume for tumor calculations
and plots this volume over time,
for control purposes the size of the tumor is also plotted.

*** Movie
to create the corresponding growth movies I used the following steps:
  1) rename
  rename -n 's/-out0[0-1][0-9]0-/-/' *.png
  2) create movie
  -r rate             set frame rate (Hz value, fraction or abbreviation)
      20 mean 20 frames per second
  for the tumor)
  ffmpeg -r 20 -f image2 -s 1920x1080 -i "fakeTum-gero_3d_8mm-typeI-sample00-gero_3month_to_5mmd-out0%03d.png_slice.png" -vcodec libx264 -crf 25 -pix_fmt yuv420p growth.mp4
  for the oxygen, not every time point considered
  ffmpeg -r 2 -f image2 -s 1920x1080 -pattern_type glob -i "*po2vessels.png" -vcodec libx264 -crf 25 -pix_fmt yuv420p o2.mp4
"""
if __name__ == '__main__':
  import os.path, sys
  #sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'..'))
  sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../..'))
  
import sys
import numpy as np
#import matplotlib

#import identifycluster
#if (identifycluster.getname()=='snowden' or identifycluster.getname()=='durga'):
#  matplotlib.use('agg')
#import matplotlib.pyplot as plt
''' this imports matplotlib and takes care of the cluster settings'''
import mpl_utils 
import matplotlib.pyplot as plt

import qsub
import myutils
import h5files

from krebs.detailedo2 import PO2ToSaturation, OpenVesselAndTumorGroups, chb_of_rbcs
from krebs.detailedo2Analysis import DataDetailedPO2
from krebs.analyzeGeneral import DataBasicVessel

def estimate_ratio_hypoxic(oxy_grp,threshold):
  dataman = myutils.DataManager(2, [DataDetailedPO2(), DataBasicVessel()])

  gvessels, gtumor = OpenVesselAndTumorGroups(oxy_grp)

  po2vessels, po2field_ld, po2field, parameters = dataman('detailedPO2', oxy_grp)
  po2vessels = np.average(po2vessels, axis=0)
  print 'po2vessels:', po2vessels.min(), po2vessels.max()
  print 'po2field:', np.amin(po2field), np.amax(po2field)
  #tissueOxygen = np.asarray(oxy_grp['po2field'])
  #print(tissueOxygen)
  
  #I neglect 5 entries from the boarder
  oxy_np_field= np.asarray(po2field)
  border=15
  cropped_oxy = oxy_np_field[border:-border,border:-border,border:-border]
  hypoxic_Tissue = cropped_oxy<threshold
  hypoxic_counts = np.sum(hypoxic_Tissue[:])  
  number_of_boxes = cropped_oxy.shape[0]*cropped_oxy.shape[1]*cropped_oxy.shape[2]
  #times volume of each box
  cropped_volume = number_of_boxes*np.power(po2field_ld.scale,3)
  print('considerd volume of: %f mum^3' % cropped_volume)
  print('considerd volume of: %f mm^3' % (cropped_volume/1e9))
  hypoxic_fraction = float(hypoxic_counts)/float(number_of_boxes)
  print('hypoxic fraction: %s ' % hypoxic_fraction)
  hypoxic_volume = hypoxic_counts*np.power(po2field_ld.scale,3)
  #to mm
  hypoxic_volume = hypoxic_volume/1e9
  return hypoxic_fraction,hypoxic_volume
''' 
    The Clinical Importance of Assessing Tumor Hypoxia: Relationship of Tumor Hypoxia to Prognosis and Therapeutic Opportunities, doi:  10.1089/ars.2013.5378, 2014
    normoxisch: x > 2.5mmHg
    hypoxisch:  0.1 mmHg > x < 2.5 mmHg
    annoxisch: x < 0.1 mmHg
    '''
def estimate_annoxic_hypoxic_normoxic(oxy_grp, tumorradius,threshold1,threshold2):
  dataman = myutils.DataManager(2, [DataDetailedPO2(), DataBasicVessel()])

  gvessels, gtumor = OpenVesselAndTumorGroups(oxy_grp)

  po2vessels, po2field_ld, po2field, parameters = dataman('detailedPO2', oxy_grp)
  po2vessels = np.average(po2vessels, axis=0)
  print 'po2vessels:', po2vessels.min(), po2vessels.max()
  print 'po2field:', np.amin(po2field), np.amax(po2field)
  #tissueOxygen = np.asarray(oxy_grp['po2field'])
  #print(tissueOxygen)
  
  #I neglect 5 entries from the boarder
  oxy_np_field= np.asarray(po2field)
  #find minimal dimension
  min_dim = np.min(oxy_np_field.shape)
  print("tumorradius: %f, ld: %f" %(tumorradius,po2field_ld.scale))
  border0=int(np.floor(oxy_np_field.shape[0]/2)-np.ceil(tumorradius/po2field_ld.scale))
  border1=int(np.floor(oxy_np_field.shape[1]/2)-np.ceil(tumorradius/po2field_ld.scale))
  border2=int(np.floor(oxy_np_field.shape[2]/2)-np.ceil(tumorradius/po2field_ld.scale))
  
  print("border0: %i" %border0)
  print("border1: %i" %border1)
  print("border2: %i" %border2)
  cropped_oxy = oxy_np_field[border0:-border0,border1:-border1,border2:-border2]
  annoxic_Tissue = cropped_oxy<threshold1
  annoxic_counts = np.sum(annoxic_Tissue[:])
  hypoxic_Tissue = np.logical_and(cropped_oxy<threshold2, cropped_oxy>threshold1)
  #hypoxic_Tissue = hypoxic_Tissue>threshold1
  hypoxic_counts = np.sum(hypoxic_Tissue[:])
  normoxic_Tissue = cropped_oxy>threshold2
  normoxic_counts = np.sum(normoxic_Tissue[:])
  number_of_boxes = cropped_oxy.shape[0]*cropped_oxy.shape[1]*cropped_oxy.shape[2]
  ''' volume correction for not sampling a sphere, but a cube
  \frac{volume(sphere)}{volume(cube of 2 times radius)} = \pi/6
  '''
  volume_correction_factor = np.pi/6.
  #volume_correction_factor = 1.
  #times volume of each box
  cropped_volume = number_of_boxes*np.power(po2field_ld.scale,3)
  print('considerd volume of: %f mum^3' % cropped_volume)
  print('considerd volume of: %f mm^3' % (cropped_volume/1e9))
  hypoxic_fraction = float(hypoxic_counts)/float(number_of_boxes)
  print('hypoxic fraction: %s ' % hypoxic_fraction)
  annoxic_volume = annoxic_counts*np.power(po2field_ld.scale,3) * volume_correction_factor
  #to mm
  annoxic_volume = annoxic_volume/1e9
  normoxic_volume = normoxic_counts*np.power(po2field_ld.scale,3) * volume_correction_factor
  #to mm
  normoxic_volume = normoxic_volume/1e9
  hypoxic_volume = hypoxic_counts*np.power(po2field_ld.scale,3) * volume_correction_factor
  #to mm
  hypoxic_volume = hypoxic_volume/1e9
  #return hypoxic_fraction,hypoxic_volume
  return annoxic_volume,hypoxic_volume,normoxic_volume
def run(goodArguments):
  print('starting with arguments: %s' % goodArguments)
  no_files = len(goodArguments.oxygenFiles)
  
  hypoxicVolumes=[]
  tumorVolumes=[]
  threshold = 15
  test = myutils.MeanValueArray.empty()
  for aFile in goodArguments.oxygenFiles:
    with h5files.open(aFile.name) as f:
      try:
        if not 'po2' in f:
          raise AssertionError('no proper oxygen file: %s!'%f)
      except Exception, e:
        print e.message
        sys.exit(-1)
      paths = myutils.walkh5(f, 'po2/out*')
      print('found paths: %s' % paths)
      hypoxicVolumes_per_time=[]
      tumorVolumes_per_time=[]
      timepoints=[]
      for path in paths:
        hypoxicFraction,hypoxicTissueVolume = estimate_ratio_hypoxic(f[path], threshold)
        hypoxicVolumes_per_time.append(hypoxicTissueVolume)      
        t=f[path]['SOURCE'].attrs['time']
        r=f[path]['SOURCE/tumor'].attrs['TUMOR_RADIUS']
        volume=4/3.*3.1414*r*r*r/1e9
        tumorVolumes_per_time.append(volume)
        timepoints.append(t)
      hypoxicVolumes.append(hypoxicVolumes_per_time)
      tumorVolumes.append(tumorVolumes_per_time)
      
  fig1 = plt.figure()
  ax1 = fig1.add_subplot(111)
#  ax1.scatter(timepoints,np.mean(hypoxicVolumes,0),color='r',label=r"hypoxic($PO_2$<%.1f mmHg)" %threshold)
#  ax1.scatter(timepoints,np.mean(tumorVolumes,0),color='b',label=r"tumor radius")
  ax1.errorbar(timepoints,np.mean(hypoxicVolumes,0),np.std(hypoxicVolumes,0),color='r',label=r"hypoxic($PO_2$<%.1f mmHg)" %threshold)
  ax1.scatter(timepoints,np.mean(tumorVolumes,0),color='b',label=r"tumor radius")
  
  ax1.set_xlabel('time/ s')
  ax1.set_ylabel(r"volume/$mm^3$")
  #ax1.set_ylabel(r"hypoxic volume ($PO_2$<%.1f mmHg)/$mm^3$" %threshold)
  ax1.legend(loc=2)
  common = os.path.commonprefix([afile.name for afile in goodArguments.oxygenFiles])
  ax1.set_title('file: %s' % common)
  with mpl_utils.PdfWriter('hypoxic_%s.pdf' % common) as pdfpages:
    pdfpages.savefig(fig1, postfix='_vesselsglobal')
  
def run_out_in_single_file(goodArguments):
  print('starting with arguments: %s' % goodArguments)
  no_files = len(goodArguments.oxygenFiles)
  
  annoxicVolumes=[]
  hypoxicVolumes=[]
  normoxicVolumes=[]
  tumorVolumes=[]
  #threshold = 15
  threshold1 = 0.1
  threshold2 = 2.5
  test = myutils.MeanValueArray.empty()
  
  annoxicVolumes_per_time=[]
  hypoxicVolumes_per_time=[]
  normoxicVolumes_per_time=[]
  tumorVolumes_per_time=[]
  timepoints=[]
  
  for aFile in goodArguments.oxygenFiles:
    with h5files.open(aFile.name) as f:
      try:
        if not 'po2' in f:
          raise AssertionError('no proper oxygen file: %s!'%f)
      except Exception, e:
        print e.message
        sys.exit(-1)
      paths = myutils.walkh5(f, 'po2/out*')
      print('found paths: %s' % paths)
      
      for path in paths:
        #hypoxicFraction,hypoxicTissueVolume = estimate_ratio_hypoxic(f[path], threshold)
        t=f[path]['SOURCE'].attrs['time']
        r=f[path]['SOURCE/tumor'].attrs['TUMOR_RADIUS']
        volume=4/3.*3.1414*r*r*r/1e9
        #volume=(2*r)*(2*r)*(2*r)/1e9
        tumorVolumes_per_time.append(volume)
        t = t/(3600.0*24)#days
        timepoints.append(t)
        annoxicTissueVolume,hypoxicTissueVolume,normoxicTissueVolume = estimate_annoxic_hypoxic_normoxic(f[path],r,threshold1,threshold2)
        annoxicVolumes_per_time.append(annoxicTissueVolume)
        hypoxicVolumes_per_time.append(hypoxicTissueVolume)
        normoxicVolumes_per_time.append(normoxicTissueVolume)
      annoxicVolumes.append(annoxicVolumes_per_time)
      hypoxicVolumes.append(hypoxicVolumes_per_time)
      normoxicVolumes.append(normoxicVolumes_per_time)
      tumorVolumes.append(tumorVolumes_per_time)
      
  print("timepoints: %s" % timepoints)    
  fig1 = plt.figure()
  ax1 = fig1.add_subplot(111)
#  ax1.scatter(timepoints,np.mean(hypoxicVolumes,0),color='r',label=r"hypoxic($PO_2$<%.1f mmHg)" %threshold)
#  ax1.scatter(timepoints,np.mean(tumorVolumes,0),color='b',label=r"tumor radius")
  ax1.errorbar(timepoints,np.mean(annoxicVolumes,0),np.std(annoxicVolumes,0),color='k',label=r"annoxic($PO_2$<%.1f mmHg)" % float(threshold1))
  ax1.errorbar(timepoints,np.mean(hypoxicVolumes,0),np.std(hypoxicVolumes,0),color='r',label=r"hypoxic(%0.1f < $PO_2$ < %.1f mmHg)" % (float(threshold1),float(threshold2)))
  ax1.errorbar(timepoints,np.mean(normoxicVolumes,0),np.std(normoxicVolumes,0),color='b',label=r"normoxic($PO_2$>%.1f mmHg)" % float(threshold2))
  #ax1.errorbar(timepoints,np.mean(hypoxicVolumes,0),np.std(hypoxicVolumes,0),color='r',label=r"hypoxic($PO_2$<%.1f mmHg)" %threshold)
  ax1.scatter(timepoints,np.mean(tumorVolumes,0),color='b',label=r"tumor volume")
  
  ax1.set_xlabel('time/ d')
  ax1.set_ylabel(r"volume/$mm^3$")
  #ax1.set_ylabel(r"hypoxic volume ($PO_2$<%.1f mmHg)/$mm^3$" %threshold)
  ax1.legend(loc=2)
  common = os.path.commonprefix([afile.name for afile in goodArguments.oxygenFiles])
  ax1.set_title('file: %s' % common)
  with mpl_utils.PdfWriter('hypoxic_%s_threshold_%0.1f_%0.1f.pdf' % (common,threshold1,threshold2)) as pdfpages:
    pdfpages.savefig(fig1, postfix='_vesselsglobal')

if __name__ == '__main__':
  import argparse
  parser = argparse.ArgumentParser(description='Plot hypoxic volume over time')  
  #parser.add_argument('oxygenFile', help='Valid configuration are found in /py/krebsjobs/parameters/fakeTumorParams.py')
  parser.add_argument('oxygenFiles', nargs='+', type=argparse.FileType('r'), default=sys.stdin, help='oxygen files to calculate')

  goodArguments, otherArguments = parser.parse_known_args()
  qsub.parse_args(otherArguments)
  
  #tumorParameterName = goodArguments.tumParamSet
  #create filename due to former standards
  #filenames=[]
  #for fn in goodArguments.vesselFileNames:
  #  filenames.append(fn.name)
    
  try:
#    if not (tumorParameterName in dir(parameterSets)) and (not 'auto' in tumorParameterName):
#        raise AssertionError('Unknown parameter set %s!' % tumorParameterName)
    for f in goodArguments.oxygenFiles:
        if not os.path.isfile(f.name):
            raise AssertionError('The file %s is not present!'%fn)
  except Exception, e:
    print e.message
    sys.exit(-1)
    
  #run(goodArguments)
  run_out_in_single_file(goodArguments)
#  if not 'auto' in tumorParameterName:
#    factory = getattr(parameterSets, tumorParameterName)
#    if type(factory).__name__ == 'function':
#      configs = factory(len(filenames))
#      for fn, cfg in zip(filenames, configs):
#        run(fn, factory.name, cfg, '4GB', 2.)
#    else:
#      for fn in filenames:
#        #run(fn, tumorParameterName, factory, '4GB', 2.)
#        run(fn, tumorParameterName, factory, '2GB', 5.)
#  else:
#    for fn in filenames:
#      for t in typelist:
#        if t in fn:
#          print(tumorParameterName)
#          tumorParameterName = tumorParameterName[5:]#removes auto_
#          type_to_paramset = create_auto_dicts(tumorParameterName+'_')
#          tumorParameterName = type_to_paramset[t]
#          factory = getattr(parameterSets, tumorParameterName)
#          #use fixed timestep for tumor simulations
#          factory['adaption']['delta_t'] = 0.10
#        #run(fn, tumorParameterName, factory, '4GB', 2.)
#          run(fn, tumorParameterName, factory, '2GB', 5.)