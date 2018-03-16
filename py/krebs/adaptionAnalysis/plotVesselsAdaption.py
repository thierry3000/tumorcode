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
Created on Thu Feb 25 14:53:15 2016

@author: thierry
"""
if __name__ == '__main__':
  import os.path, sys
  sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../..'))
  sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'..'))
  
import os,sys
from os.path import join, basename, dirname, splitext, commonprefix
import time
import krebsutils
import h5py
import numpy as np
import itertools
import extensions
import collections
import posixpath
import math

import mpl_utils
from mystruct import Struct
import myutils
from myutils import f2l

import matplotlib
import matplotlib.pyplot as pyplot
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LogNorm 
from matplotlib.font_manager import FontProperties


from quantities import Prettyfier

from plotVessels import *
from analyzeMurray import *
from adaption import getVesselTypes
from analyzeGeneral import getGeometricData,getTotalPerfusion,generate_adaption_data_average_rBV

import adaption 
from scipy.optimize import fsolve
from scipy.optimize import minimize_scalar

#multiprocessing
#import multiprocessing
import threading

type_to_label_dict={"typeA-": "RC1","typeB-": "RC2","typeC-": "RC3","typeD-": "RC4","typeE-": "RC5","typeF-": "RC6","typeG-": "RC7","typeH-": "RC8","typeI-": "RC9"}  

def getVesselTypes(vessel_groups):
  veins = []
  arteries = []
  capillaries = []
  if ('edges' in vessel_groups.keys()):# in fact it is only one group!
    g = vessel_groups    
    flags = np.array(krebsutils.read_vessels_from_hdf(g,['flags'])[1])
    circulated = np.bitwise_and(flags,krebsutils.CIRCULATED)
    circulated_indeces = np.nonzero(circulated)
    veins = np.bitwise_and(flags,krebsutils.VEIN)
    veins = veins[circulated_indeces]
    arteries = np.bitwise_and(flags,krebsutils.ARTERY)
    arteries = arteries[circulated_indeces]
    capillaries = np.bitwise_and(flags,krebsutils.CAPILLARY)
    capillaries = capillaries[circulated_indeces]
    veins = np.nonzero(veins)
    arteries =np.nonzero(arteries)
    capillaries = np.nonzero(capillaries)
  else:
    for g in vessel_groups:
      flags = krebsutils.read_vessels_from_hdf(g,['flags'])
      goods = np.sum(np.bitwise_and(flags[1],krebsutils.VEIN))
      veins.append(goods)
      goods = np.sum(np.bitwise_and(flags[1],krebsutils.ARTERY))
      arteries.append(goods)
      goods = np.sum(np.bitwise_and(flags[1],krebsutils.CAPILLARY))
      capillaries.append(goods)
    
  return veins, arteries, capillaries
    
def printBarPlot_vesseltype_on_root_node_configuration(filenames, pdfpages):
  means_veins = []
  means_arteries = []
  means_capillaries = []
  
  std_veins = []
  std_arteries = []
  std_capillaries = []
  
  
  counter = 0
  for t in 'typeA- typeB- typeC- typeD- typeE- typeF- typeG- typeH- typeI-'.split():
    #print('rBF for type: %s' % t)
    filteredFiles = filter( lambda fn: t in fn,filenames)
  
    files = [h5files.open(fn, 'r+') for fn in filteredFiles]
    groups_of_type = list(itertools.chain.from_iterable(myutils.walkh5(f, 'adaption/recomputed', return_h5objects=True) for f in files))
    
    veins, arteries, capillaries = getVesselTypes(groups_of_type)
    
    #print("something: %i " % len(perfusion_data_without))
    #print("something: %i " % len(perfusion_data_with))
    if (len(veins)>0):
      means_veins.append( np.average(veins))
      std_veins.append( np.std(veins))
    
    if (len(arteries)>0):
      means_arteries.append( np.average(arteries))
      std_arteries.append( np.std(arteries))
      
    if (len(capillaries)>0):
      means_capillaries.append( np.average(capillaries))
      std_capillaries.append( np.std(capillaries))
    counter= counter +1
  
  print("means_arteries: %i " % len(means_arteries))
  print(means_arteries)
  print("means_veins: %i " % len(means_veins))
  print(means_veins)
  print("means_capillaries: %i " % len(means_capillaries))
  print(means_capillaries)
  if not counter == len(means_arteries):
    raise AssertionError("Not all configurations found!")
  plt = pyplot
  fig, ax = plt.subplots()
  index = np.arange(counter)
  bar_width = 0.25
  opacity = 0.4
  error_config = {'ecolor': '0.3'}
  rects1 = plt.bar(index, means_arteries, bar_width,
                 alpha=opacity,
                 color='r',
                 yerr=std_arteries,
                 error_kw=error_config,
                 label='Arteries')
  rects2 = plt.bar(index + bar_width, means_veins, bar_width,
                 alpha=opacity,
                 color='b',
                 yerr=std_veins,
                 error_kw=error_config,
                 label='Veins',
                 hatch="/")
  rects3 = plt.bar(index + 3*bar_width, means_capillaries, bar_width,
                 alpha=opacity,
                 color='y',
                 yerr=std_capillaries,
                 error_kw=error_config,
                 label='Capillaries',
                 hatch="/")

  plt.xlabel('Configuration')
  plt.ylabel(r'no. vessels')
  plt.title('Vessels per type and root node configuration')
  plt.xticks(index + bar_width, ('RC1', 'RC2', 'RC3', 'RC4', 'RC5', 'RC6', 'RC7', 'RC8', 'RC9'))
  #plt.legend(loc='best')
  lgd = ax.legend(loc='center left', bbox_to_anchor=(1,1))
  
  
  plt.tight_layout()
  pdfpages.savefig(fig, bbox_extra_artists=(lgd,), bbox_inches='tight')

  
def printBarPlot_rBF(dataman, fmeasure, filenames, options, pdfpages):
  means_without = []
  std_without = []
  
  means_with = []
  std_with = []

  

  destination_group = fmeasure.require_group('adaption/geometry/rBF')  
  
  counter = 0

  if(options.two):
    typelist = 'typeD- typeE- typeG- typeH-'
  else:
    typelist = 'typeA- typeB- typeC- typeF- typeI-'
  if(options.all_types):
    typelist = 'typeA- typeB- typeC- typeD- typeE- typeF- typeG- typeH- typeI-'
  if(options.single):
    typelist = 'typeF- '
    
  for t in typelist.split():
    print('rBF for type: %s' % t)
    filteredFiles = filter( lambda fn: t in fn,filenames) 
    files = [h5files.open(fn, 'r+') for fn in filteredFiles]

    if(len(files)>0):
      w_avg, w_std, wo_avg, wo_std = generate_adaption_data_average_rBF(dataman, files, destination_group, t[:-1] + '_rBF')
    else:
      continue
    means_without.append(wo_avg[0])
    std_without.append(wo_std[0])
    means_with.append(w_avg[0])
    std_with.append(w_std[0])
    
#    means_a_with.append(w_avg[1])
#    std_a_with.append(w_std[1])
#    means_a_without.append(wo_avg[1])
#    std_a_without.append(wo_std[1])
#    
#    means_v_with.append(w_avg[2])
#    std_v_with.append(w_std[2])
#    means_v_without.append(wo_avg[2])
#    std_v_without.append(wo_std[2])
#    
#    means_c_with.append(w_avg[3])
#    std_c_with.append(w_std[3])
#    means_c_without.append(wo_avg[3])
#    std_c_without.append(wo_std[3])

    counter= counter +1    
  
  print("coutner: %i"%counter)
  if counter >0:
    print("means_without_all: %i" % len(means_without))
    print("std_without_all: %i" % len(std_without))  
  
  if not counter == len(means_without):
    raise AssertionError("Not all configurations found!")
  plt = pyplot
  fig, ax = plt.subplots()
  index = np.arange(counter)
  bar_width = 0.3
  opacity = 0.4
  error_config = {'ecolor': '0.3'}
  rects1 = plt.bar(index, means_without, bar_width,
                 alpha=opacity,
                 color='w',
                 yerr=std_without,
                 error_kw=error_config,
                 label='MW Algo')
  rects2 = plt.bar(index + bar_width, means_with, bar_width,
                 alpha=opacity,
                 color='w',
                 yerr=std_with,
                 error_kw=error_config,
                 label='Adaption',
                 hatch="/")

  plt.xlabel('Configuration')
  plt.ylabel(r'rBF/ $ml\,Blood\,ml^{-1} min^{-1}$')
  plt.title('Perfusion by Root Node configuration and adaption')
  plt.xticks(index + 2*bar_width, ('RC1', 'RC2', 'RC3', 'RC4', 'RC5', 'RC6', 'RC7', 'RC8', 'RC9'))
  #plt.legend(loc='best')
  lgd = ax.legend(loc='center left', bbox_to_anchor=(1,1))
  
  
  plt.tight_layout()
  pdfpages.savefig(fig, bbox_extra_artists=(lgd,), bbox_inches='tight')

def PlotQdevs_unnormalized(filenames, pdfpages):
  #this is from
  # http://matplotlib.org/examples/axes_grid/inset_locator_demo.html
  from mpl_toolkits.axes_grid1.inset_locator import inset_axes, zoomed_inset_axes
  from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
  qdevVector=[]
  plt = pyplot
  
  if(options.two):
    typelist = 'typeD- typeE- typeG- typeH-'
  else:
    typelist = 'typeA- typeB- typeC- typeF- typeI-'
  if(options.all_types):
    typelist = 'typeA- typeB- typeC- typeD- typeE- typeF- typeG- typeH- typeI-'
  if(options.single):
    typelist = 'typeF- '
    
  for t in typelist.split():
    fig= plt.figure()
    ax1 = plt.axes()
    ax1.set_title(t)
    fig.add_subplot(ax1)
    print('Qdev for type: %s' % t)
     
    filteredFiles = filter( lambda fn: t in fn,filenames)
    files = [h5files.open(fn, 'r+') for fn in filteredFiles]
    groups_with_adaption = list(itertools.chain.from_iterable(myutils.walkh5(f, 'adaption/vessels_after_adaption', return_h5objects=True) for f in files))
    
    for agroup in groups_with_adaption:
      print("filename: %s" % agroup.file.filename)
      qdev_vec = np.zeros(2999)
      qdev_vec_read = np.asarray(agroup['qdev'])
      qdev_vec_read = np.sqrt(qdev_vec_read)
      qdev_vec[0:len(qdev_vec_read)]= qdev_vec_read
      
      if 'qdev_mat' not in locals():
        qdev_mat = qdev_vec
      else:
        qdev_mat = np.vstack((qdev_mat, qdev_vec))
    average_for_type = np.average(qdev_mat,0)
    std_for_type = np.std(qdev_mat,0)
    print("len avg: %i, len std: %i" % (len(average_for_type),len(std_for_type)))
    ax1.errorbar(np.arange(len(average_for_type)),average_for_type, yerr= std_for_type, errorevery=20)
    #ax1.set_yscale('log')    
    ax1.set_ylabel(r'$\sum_{i=1}^{no. vessels} (r_i-r_{i-1})^2$')    
    ax1.set_xlabel('Iteration')

    axins = inset_axes(ax1,
                   width="70%",  # width = 30% of parent_bbox
                   height="30%",  # height : 1 inch
                   loc=1)
    axins.errorbar(np.arange(len(average_for_type))[-1500:],average_for_type[-1500:], yerr= std_for_type[-1500:], errorevery=20) 
    #axins.set_yscale('log')    
    ax1.grid()    
        
    del qdev_mat
    del average_for_type
    del std_for_type
    plt.tight_layout()
    pdfpages.savefig(fig)
    
def PlotQdevs_normalized(filenames, options, pdfpages):
  #this is from
  # http://matplotlib.org/examples/axes_grid/inset_locator_demo.html
  from mpl_toolkits.axes_grid1.inset_locator import inset_axes, zoomed_inset_axes
  from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
  qdevVector=[]
  plt = pyplot
  
  if(options.two):
    typelist = 'typeD- typeE- typeG- typeH-'
  else:
    typelist = 'typeA- typeB- typeC- typeF- typeI-'
  if(options.all_types):
    typelist = 'typeA- typeB- typeC- typeD- typeE- typeF- typeG- typeH- typeI-'
  if(options.single):
    typelist = 'typeF- '
    
  for t in typelist.split():
    fig= plt.figure()
    ax1 = plt.axes()
    ax1.set_title(t)
    fig.add_subplot(ax1)
    print('Qdev for type: %s' % t)
     
    filteredFiles = filter( lambda fn: t in fn,filenames)
    files = [h5files.open(fn, 'r+') for fn in filteredFiles]
    groups_with_adaption = list(itertools.chain.from_iterable(myutils.walkh5(f, 'adaption/vessels_after_adaption', return_h5objects=True) for f in files))

    if(len(groups_with_adaption)==0):
      continue
    
    for agroup in groups_with_adaption:
      print("filename: %s" % agroup.file.filename)
      vessel_count = agroup['edges'].attrs['COUNT']
      nqdev_vec = np.zeros(300)
      nqdev_vec_read = np.asarray(agroup['nqdev'])
      if(len(nqdev_vec_read)<300):
        nqdev_vec[0:len(nqdev_vec_read)] = nqdev_vec_read
      else:
        continue

      if 'nqdev_mat' not in locals():
        nqdev_mat = nqdev_vec
      else:
        nqdev_mat = np.vstack((nqdev_mat, nqdev_vec))
  
    ax1.plot(nqdev_mat.transpose())
    ax1.set_ylabel(r'$\frac{1}{N}\sqrt{\sum_{i=1}^{no. vessels} (r_i-r_{i-1})^2}$')
    ax1.set_yscale('log')
    ax1.set_xlabel('Iteration')
  
    ax1.grid()    
    #plt.tight_layout()
    pdfpages.savefig(fig)

@myutils.UsesDataManager
def generate_adaption_data_of_group_rBF(datamanager, destination_group, f):
  datanames = 'with_adaption_rBF without_adaption_rBF'.split()
  # structure in HDF file:
  #            gmeasure/groupname/
  #                               rbv  <- dataset
  #                               a    <- dataset
  #                               v    <- dataset
  
  def read(gmeasure, groupname):
    gmeasure = gmeasure[groupname]
    return [gmeasure[name][()] for name in datanames ]
    
  def write(gmeasure, groupname):
    group_without_adaption = f['adaption/recomputed']
    group_with_adaption = f['adaption/vessels_after_adaption']
    gmeasure = gmeasure.create_group(groupname)
    for (group,name) in zip([group_without_adaption,group_with_adaption],datanames):
      perfusion_data = getTotalPerfusion([group])*60
      gmeasure.create_dataset(name, data=perfusion_data)
    #geometric_data = getGeometricData([vesselgroup])
#    rBF = perfusion_data
#    for name, data in zip(datanames, [rBF]):
#      gmeasure.create_dataset(name, data = data)
  
  ret = myutils.hdf_data_caching(read, write, destination_group, (f.filename), (1, ))  
  # so, gmeasure/groupname is "/measurements/adaption".(None, 1,) are version number specifications, 
  # one for each path component. None means no version number is checked. If version number is larger than
  # stored number, then data is recomputed instead of loaded.
  return ret
  
@myutils.UsesDataManager
def generate_adaption_data_average_rBF(datamanager, inputfiles, destination_group, destination_name):
  def process(vesselgroups):
    tmp = []
    for f in inputfiles:
      rbF_with,rBF_without = generate_adaption_data_of_group_rBF(datamanager, destination_group, f)
      tmp.append([rbF_with,rBF_without])
      #  generate_adaption_data_of_group_rBF(datamanager, destination_group, vesselgroup))
    avg_w,avg_wo = np.average(tmp, axis = 0)
    std_w,std_wo = np.std(tmp, axis = 0)    
    #avg_w = np.average(tmp[0], axis = 0)
    #std_w = np.std(tmp[0], axis = 0)
    return avg_w, std_w, avg_wo, std_wo
  
  def write(gmeasure, groupname):
    gmeasure = gmeasure.create_group(groupname)
    avg_with,std_with, avg_without, std_without = process(inputfiles)
    #for name in ['with_adaption', 'without_adaption']:
    #  avg_with, std_with, avg_without, std_without = process(groups)
    gmeasure.create_dataset('with_adaption_rBF_avg', data = avg_with)
    gmeasure.create_dataset('with_adaption_rBF_std', data = std_with)
    gmeasure.create_dataset('without_adaption_rBF_avg', data = avg_without)
    gmeasure.create_dataset('without_adaption_rBF_std', data = std_without)    
  
#    groups_without_adaption = [f['adaption/recomputed'] for f in inputfiles]
#    groups_with_adaption = [f['adaption/vessels_after_adaption'] for f in inputfiles]
#
#    for name, groups in zip(['with_adaption', 'without_adaption'], [groups_with_adaption, groups_without_adaption]):
#      avg, std = process(groups)
#      gmeasure.create_dataset(name+'_rBF_avg', data = avg)
#      gmeasure.create_dataset(name+'_rBF_std', data = std)
      # will give datasets:
      #     with_adaption_avg, with_adaption_std, without_adaption_avg, and without_adaption_std.
      # Each is tuple containing (rbv, a, v, c) returned by generate_adaption_data_of_group(...)
  
  def read(gmeasure, groupname):
    gmeasure = gmeasure[groupname]
    r1 = gmeasure['with_adaption_rBF_avg']
    r2 = gmeasure['with_adaption_rBF_std']
    r3 = gmeasure['without_adaption_rBF_avg']
    r4 = gmeasure['without_adaption_rBF_std']
    return (r1, r2, r3, r4)
  
  ret = myutils.hdf_data_caching(read, write, destination_group, (destination_name,), (1,))
  return ret


@myutils.UsesDataManager
def generate_adaption_data_of_group_rBV(datamanager, destination_group, f):
  datanames = 'with_adaption_rbv without_adaption_rbv'.split()
  # structure in HDF file:
  #            gmeasure/groupname/
  #                               rbv  <- dataset
  #                               a    <- dataset
  #                               v    <- dataset
  
  def read(gmeasure, groupname):
    gmeasure = gmeasure[groupname]
    return [gmeasure[name][()] for name in datanames ]
    
  def write(gmeasure, groupname):
    group_without_adaption = f['adaption/recomputed']
    group_with_adaption = f['adaption/vessels_after_adaption']
    gmeasure = gmeasure.create_group(groupname)
    for (group,name) in zip([group_with_adaption,group_without_adaption],datanames):
      geometric_data = getGeometricData([group])
      rbv, a, v, c = geometric_data[:4]
      gmeasure.create_dataset(name, data = rbv)
#    
#    for name, data in zip(datanames, [rbv, a, v, c]):
#      gmeasure.create_dataset(name, data = data)
  
  ret = myutils.hdf_data_caching(read, write, destination_group, (f.filename), (1, ))  
  # so, gmeasure/groupname is "/measurements/adaption".(None, 1,) are version number specifications, 
  # one for each path component. None means no version number is checked. If version number is larger than
  # stored number, then data is recomputed instead of loaded.
  return ret
    
  

  
  
@myutils.UsesDataManager
def generate_murray_data_hist(datamanager, inputfiles, destination_group, destination_name):
  def process(vesselgroups):
    list_of_stuff = 'tmp_daughter1 tmp_daughter2 tmp_mother'.split()
    k=0
    my_dict_stuff = dict()
    for vesselgroup in vesselgroups:
      for endity in list_of_stuff:
        if endity in my_dict_stuff:
          my_dict_stuff[endity] = np.concatenate((my_dict_stuff[endity],np.asarray(generate_murray_data_of_group(datamanager, vesselgroup,destination_group)[k])),axis=0)
        else:
          anarray = np.asarray(generate_murray_data_of_group(datamanager, vesselgroup,destination_group)[k])
          my_dict_stuff[endity] = anarray        
          #setattr(self,endity,anarray)
          k=k+1
#    for vesselgroup in vesselgroups:
#      if('tmp_daughter1' in locals()):
#        np.concatenate((tmp_daughter1,np.asarray(generate_murray_data_of_group(datamanager, vesselgroup,destination_group)[0])),axis=0)
#      else:
#        tmp_daughter1 = np.asarray(
#          generate_murray_data_of_group(datamanager, vesselgroup,destination_group)[0]) 
#      if('tmp_daughter2' in locals()):
#        np.concatenate((tmp_daughter2,np.asarray(
#          generate_murray_data_of_group(datamanager, vesselgroup,destination_group)[1])))
#      else:
#        tmp_daughter2 = np.asarray(
#          generate_murray_data_of_group(datamanager, vesselgroup,destination_group)[1])
#      if('tmp_mother' in locals()):
#        np.concatenate((tmp_mother,np.asarray(
#          generate_murray_data_of_group(datamanager, vesselgroup,destination_group)[2])))
#      else:
#        tmp_mother = np.asarray(
#          generate_murray_data_of_group(datamanager, vesselgroup,destination_group)[2])

    #xedges=np.logspace(1,1.8,70)
    xedges=np.linspace(2.,15,70)
    yedges=np.linspace(-0.5,0.5,25)
    def remove_zeros(x):
      good_indexes = np.where(x>0)[0]
      return x[good_indexes]
    theo_mother = map(lambda x,y: np.power(x**3+y**3,1/3.),my_dict_stuff['tmp_daughter1'], my_dict_stuff['tmp_daughter2'])
    theo_mother = np.asarray(theo_mother)
    my_dict_stuff['theo_mother']=theo_mother
    good_indexes = np.where(theo_mother>0)[0]
    
    H, x_edges, y_edges = np.histogram2d(my_dict_stuff['tmp_daughter1'][good_indexes],(my_dict_stuff['tmp_mother'][good_indexes]-my_dict_stuff['theo_mother'][good_indexes])/my_dict_stuff['tmp_daughter1'][good_indexes], bins=(xedges,yedges), normed=True)

    def calc_alphas_effective(daughter,mother):
      res = []
      threads = []
      def func(alpha,*data ):
        return np.power(data[0],alpha)-2*np.power(data[1], alpha)
      if 0:
        for (r_m,r_d) in zip (mother, daughter):
          params = (r_m,r_d)
          #res.append(fsolve(func,3,args = params))
          ares = minimize_scalar(func,bounds=(0.8,6),method='bounded',args = params)
          if (ares.status == 0):        
            res.append(ares.x)
      #try parallel version
      if 1:
        def worker(r_m,r_d):
          params = (r_m,r_d)
          ares = minimize_scalar(func,bounds=(0.1,10),method='bounded',args = params)
          if (ares.status == 0):        
            res.append(ares.x)
          
        for (r_m,r_d) in zip (mother, daughter):
          t = threading.Thread(target=worker,args=(r_m,r_d,))
          threads.append(t)
          t.start()
      return res
    def calc_betas(daughter,mother):
      return [r_d/r_m for (r_d,r_m) in zip(daughter,mother)]
    tmp_long_daughters = np.hstack((my_dict_stuff['tmp_daughter1'][good_indexes],my_dict_stuff['tmp_daughter2'][good_indexes]))
    tmp_long_mothers = np.hstack((my_dict_stuff['tmp_mother'][good_indexes],my_dict_stuff['tmp_mother'][good_indexes]))   
    betas = calc_betas(tmp_long_daughters,tmp_long_mothers)
    #alphas_effective = calc_alphas_effective(tmp_daughter1[good_indexes],tmp_mother[good_indexes])
    alphas_effective = calc_alphas_effective(tmp_long_daughters,tmp_long_mothers)    
    return H, x_edges,y_edges, alphas_effective,betas
  
  def write(gmeasure, groupname):
    gmeasure = gmeasure.create_group(groupname)    
  
    groups_with_adaption = [f['adaption/vessels_after_adaption'] for f in inputfiles]

    #for groups in groups_with_adaption:
    H, edges_x, edges_y, alphas_effective, betas = process(groups_with_adaption)
    gmeasure.create_dataset('H', data = H)
    gmeasure.create_dataset('edges_x', data = edges_x)
    gmeasure.create_dataset('edges_y', data = edges_y)
    gmeasure.create_dataset('alphas_effective', data = alphas_effective)
    gmeasure.create_dataset('betas', data = betas)    
    #gmeasure.create_dataset('edges_z', data = edges[2])
    #gmeasure.create_dataset('yedges', data = yedges)
      # will give datasets:
      #     with_adaption_avg, with_adaption_std, without_adaption_avg, and without_adaption_std.
      # Each is tuple containing (rbv, a, v, c) returned by generate_adaption_data_of_group(...)
  
  def read(gmeasure, groupname):
    gmeasure = gmeasure[groupname]
    r1 = gmeasure['H']
    r2 = gmeasure['edges_x']
    r3 = gmeasure['edges_y']
    r4 = gmeasure['alphas_effective']
    r5 = gmeasure['betas']
    #r4 = gmeasure['edges_z']
    #r3 = gmeasure['yedges']
    #r4 = gmeasure['without_adaption_daughter']
    return (r1, r2, r3, r4, r5)
  
  ret = myutils.hdf_data_caching(read, write, destination_group, (destination_name,), (1,))
  return ret
  

  
@myutils.UsesDataManager
def generate_murray_data_of_group(datamanager, vesselgroup, destination_group):
  datanames = 'daughter1 daughter2 mother daughter1_scale daughter2_scale mother_scale'.split() 
  # structure in HDF file:
  #            gmeasure/groupname/
  #                               rbv  <- dataset
  #                               a    <- dataset
  #                               v    <- dataset
  
  def read(gmeasure, groupname):
    gmeasure = gmeasure[groupname]
    return [gmeasure[name][()] for name in datanames ]
    
  def write(gmeasure, groupname):
    res = ku.get_Murray2_p(vesselgroup)
    daughter1 = res[0]
    daughter2 = res[1]
    mother = res[2]
    res2 = ku.get_Murray_scale(vesselgroup)
    daughter1_scale = res2[0]
    daughter2_scale = res2[1]
    mother_scale = res2[2]
    gmeasure = gmeasure.create_group(groupname)
    for name, data in zip(datanames, [daughter1, daughter2, mother,daughter1_scale, daughter2_scale, mother_scale]):
      gmeasure.create_dataset(name, data = data)
    
  
  ret = myutils.hdf_data_caching(read, write, destination_group,(vesselgroup.file.filename),(1,))  
  # so, gmeasure/groupname is "/measurements/adaption".(None, 1,) are version number specifications, 
  # one for each path component. None means no version number is checked. If version number is larger than
  # stored number, then data is recomputed instead of loaded.
  return ret

def printMurray(dataman, f_measure, filenames, options, pdfpages):
    filenames = adaption.get_files_with_successful_adaption(filenames)
    files = [h5files.open(fn, 'r') for fn in filenames]
    destination_group = f_measure.require_group('adaption/murray')
    if(options.two):
      typelist = 'typeD- typeE- typeG- typeH-'
    else:
      typelist = 'typeA- typeB- typeC- typeF- typeI-'
    if(options.all_types):
      typelist = 'typeA- typeB- typeC- typeD- typeE- typeF- typeG- typeH- typeI-'
      
    if(options.single):
      typelist = 'typeF- '
      
    for t in typelist.split():
      print('Murray for type: %s' % t)
      filteredFiles = filter( lambda fn: t in fn,filenames) 
      files = [h5files.open(fn, 'r+') for fn in filteredFiles]
  
      if(len(files)>0):#that means no file of dedicated type is left after filter
        hist_data, x_edges, y_edges, alphas_effective, betas = generate_murray_data_hist(dataman, files, destination_group, t[:-1])
      #data = generate_murray_data_plot3(dataman, files, destination_group, t[:-1])      
      else:
        continue
      
      if 1:
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.set_title(t)
        #fig.add_subplot(ax1)
        xedges = np.asarray(x_edges)
        yedges = np.asarray(y_edges)
        Z = np.asarray(hist_data)
        xcenters = xedges[:-1] + 0.5 * (xedges[1:] - xedges[:-1])
        ycenters = yedges[:-1] + 0.5 * (yedges[1:] - yedges[:-1])
        dx = xedges[1:]-xedges[:-1]
        dy = yedges[1:] - yedges[:-1]
        dV = np.outer(dx,dy)
        Z = dV * Z
        im = matplotlib.image.NonUniformImage(ax1, interpolation='bilinear')
        im.set_data(xcenters, ycenters, Z.transpose())
        ax1.images.append(im)
        ax1.set_xlim(xedges[0], xedges[-1])
        ax1.set_ylim(yedges[0], yedges[-1])
        #ax1.set_aspect('equal')
        ax1.set_xlabel(r'$r_{daughter}$')
        ax1.set_ylabel(r'$\frac{r_{mother,Secomb}-r_{mother,Murray}}{r_{mother,Murray}}$')
        fig.colorbar(im)

      pdfpages.savefig(fig)
@myutils.UsesDataManager
def generate_murray_alphas_histogram(datamanager, alphas, destination_group, destination_name):

  def write(gmeasure,groupname):
    gmeasure = gmeasure.create_group(groupname)
    x_edges = np.linspace(0.0,1.3,50)
    h, edges = np.histogram(np.reciprocal(alphas),bins=x_edges)
    gmeasure.create_dataset('h_alphas', data=h)
    gmeasure.create_dataset('x_edges_alphas', data= edges)
  def read(gmeasure, groupname):
    gmeasure = gmeasure[groupname]
    r1 = gmeasure['h_alphas']
    r2 = gmeasure['x_edges_alphas']
    return (r1,r2)
  '''
  cache number needs to be higher than other murray stuff
  since this will open the cache file for the second time
  '''
  ret = myutils.hdf_data_caching(read, write, destination_group, (destination_name,), (1,))
  return ret
def generate_murray_beta_histogram(datamanager, betas, destination_group, destination_name):

  def write(gmeasure,groupname):
    gmeasure = gmeasure.create_group(groupname)
    x_edges = np.linspace(0.0,2.5,50)
    h, edges = np.histogram(betas,bins=x_edges)
    gmeasure.create_dataset('h_betas', data=h)
    gmeasure.create_dataset('x_edges_betas', data= edges)
  def read(gmeasure, groupname):
    gmeasure = gmeasure[groupname]
    r1 = gmeasure['h_betas']
    r2 = gmeasure['x_edges_betas']
    return (r1,r2)
  '''
  cache number needs to be higher than other murray stuff
  since this will open the cache file for the second time
  '''
  ret = myutils.hdf_data_caching(read, write, destination_group, (destination_name,), (1,))
  return ret

def printMurray_alphas_effective(dataman, f_measure, filenames, options, pdfpages):
    filenames = adaption.get_files_with_successful_adaption(filenames)
    files = [h5files.open(fn, 'r') for fn in filenames]
    destination_group = f_measure.require_group('adaption/murray')
    if(options.two):
      typelist = 'typeD- typeE- typeG- typeH-'
    else:
      typelist = 'typeA- typeB- typeC- typeF- typeI-'
    if(options.all_types):
      typelist = 'typeA- typeB- typeC- typeD- typeE- typeF- typeG- typeH- typeI-'
      
    if(options.single):
      typelist = 'typeF- '
    
    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)
        
    for t in typelist.split():
      print('Murray for type: %s' % t)
      filteredFiles = filter( lambda fn: t in fn,filenames) 
      files = [h5files.open(fn, 'r+') for fn in filteredFiles]
  
      #dataman.obtain_data('basic_vessel_global', name, gvessels, cachelocation(gvessels))
      if(len(files)>0):#that means no file of dedicated type is left after filter
        hist_data, x_edges, y_edges, alphas_effective, betas = generate_murray_data_hist(dataman, files, destination_group, t[:-1])
      #data = generate_murray_data_plot3(dataman, files, destination_group, t[:-1])      
      else:
        continue
      
      if 1:
        #destination_group = destination_group.create_group('alpha_eff')
        #f_measure.require_group('adaption/alpha_eff')
        #alphas_effective = np.asarray(alphas_effective)
        h,x_edges = generate_murray_alphas_histogram(dataman, alphas_effective, destination_group, t[:-1]+'_alpha_hist')
        #rects = ax1.bar(x_edges[:-1], h, 0.5)
        x_edges = np.asarray(x_edges)
        h =np.asarray(h)
        ax1.plot((x_edges[1:]+x_edges[:-1])/2, h, '.' , label = type_to_label_dict[t])    
        
        #legend1 = ax1.legend( shadow=True,prop=fontP)
        #ax2.plot((x_edges_2[1:]+x_edges_2[:-1])/2, h2, '.' , label = type_to_label_dict[t])
        h2,x_edges2 = generate_murray_beta_histogram(dataman, betas, destination_group, t[:-1]+'_beta_hist')
        #rects = ax1.bar(x_edges[:-1], h, 0.5)
        x_edges2 = np.asarray(x_edges2)
        h2 =np.asarray(h2)
        
        #legend1 = ax1.legend( shadow=True,prop=fontP)
        ax2.plot((x_edges2[1:]+x_edges2[:-1])/2, h2, '.' , label = type_to_label_dict[t])    
    
    ax1.set_title(r'Effective murray $\frac{1}{\alpha}$ by RC')
    ax1.set_xlabel(r'$\alpha$')
    ax1.set_ylabel(r'$\# count$')
    fontP = FontProperties()
    fontP.set_size('x-small')
    ax1.legend(prop=fontP)
    
    ax2.set_title(r'Effective murray $\beta$ by RC')
    ax2.set_xlabel(r'$\beta$')
    ax2.set_ylabel(r'$\# count$')
    fontP = FontProperties()
    fontP.set_size('x-small')
    ax2.legend(prop=fontP)
        
    pdfpages.savefig(fig)

def printBarPlot_rBV(dataman, fmeasure, filenames, options, pdfpages):
  means_without_all = []
  std_without_all = []
  means_with_all = []
  std_with_all = []  
  
#  means_a_with = []
#  std_a_with = []
#  means_a_without = []
#  std_a_without = []  
#  
#  means_v_with =[]
#  std_v_with =[]
#  means_v_without = []
#  std_v_without = []
#  
#  means_c_with = []
#  std_c_with = []
#  means_c_without = []
#  std_c_without = []

  destination_group = fmeasure.require_group('adaption/geometry/rBV')  
  
  counter = 0

  if(options.two):
    typelist = 'typeD- typeE- typeG- typeH-'
  else:
    typelist = 'typeA- typeB- typeC- typeF- typeI-'
  if(options.all_types):
    typelist = 'typeA- typeB- typeC- typeD- typeE- typeF- typeG- typeH- typeI-'
    
  if(options.single):
    typelist = 'typeI- '
    
  for t in typelist.split():
    print('rBV for type: %s' % t)
    filteredFiles = filter( lambda fn: t in fn,filenames) 
    files = [h5files.open(fn, 'r+') for fn in filteredFiles]

    if(len(files)>0):#that means no file of dedicated type is left after filter
      w_avg, w_std, wo_avg, wo_std = generate_adaption_data_average_rBV(dataman, files, destination_group, t[:-1]+'_rBV_avg')
      #data = generate_murray_data_plot3(dataman, files, destination_group, t[:-1])      
    else:
      continue
    
    w_avg = np.asarray(w_avg)
    w_std = np.asarray(w_std)
    wo_avg = np.asarray(wo_avg)
    wo_std = np.asarray(wo_std)
    means_without_all.append(wo_avg[0])
    std_without_all.append(wo_std[0])
    means_with_all.append(w_avg[0])
    std_with_all.append(w_std[0])
    
#    means_a_with.append(w_avg[1])
#    std_a_with.append(w_std[1])
#    means_a_without.append(wo_avg[1])
#    std_a_without.append(wo_std[1])
#    
#    means_v_with.append(w_avg[2])
#    std_v_with.append(w_std[2])
#    means_v_without.append(wo_avg[2])
#    std_v_without.append(wo_std[2])
#    
#    means_c_with.append(w_avg[3])
#    std_c_with.append(w_std[3])
#    means_c_without.append(wo_avg[3])
#    std_c_without.append(wo_std[3])

    counter= counter +1    
  
  print("coutner: %i"%counter)
  print("means_without_all: %i" % len(means_without_all))
  print("std_without_all: %i" % len(std_without_all))
  plt = pyplot
  fig, ax = plt.subplots()
  index = np.arange(counter)
  bar_width = 0.4
  fig, ax = plt.subplots()
  rects1 = ax.bar(index,means_without_all, bar_width, yerr=std_without_all)
  rects2 = ax.bar(index+bar_width,means_with_all, bar_width, color='w', yerr=std_with_all, ecolor='b', hatch='/')
  ax.set_xticks(index + bar_width)
  def type_to_label(typelist):
    type_to_label_dict={"typeA-": "RC1","typeB-": "RC2","typeC-": "RC3","typeD-": "RC4","typeE-": "RC5","typeF-": "RC6","typeG-": "RC7","typeH-": "RC8","typeI-": "RC9"}
    label = []
    for t in typelist.split():
      label.append(type_to_label_dict[t])
    return label
  ax.set_xticklabels(type_to_label(typelist))
  ax.set_xlabel('Configuration')
  ax.legend((rects1[0],rects2[0]),('MW','Adaption'))
  ax.set_ylabel('rBV')
  ax.set_title('Volume Fractions of different vessel types by RC')

  plt.tight_layout()
  #pdfpages.savefig(fig, bbox_extra_artists=(lgd,), bbox_inches='tight')
  pdfpages.savefig(fig, bbox_inches='tight')
        

def PlotRadiusHistogram_with_cache_by_RC(dataman, f_measure, filenames, options, pdfpages):
  print("Ploting radius histogram!")
  

  counter = 0

  if(options.two):
    typelist = 'typeD- typeE- typeG- typeH-'
  else:
    typelist = 'typeA- typeB- typeC- typeF- typeI-'
  if(options.all_types):
    typelist = 'typeA- typeB- typeC- typeD- typeE- typeF- typeG- typeH- typeI-'
    
  if(options.single):
    typelist = 'typeF- '
    
  fig1, ax1 = pyplot.subplots(1,1, figsize = (mpl_utils.a4size[0]*0.8, mpl_utils.a4size[0]*0.4))
  fig2, ax2 = pyplot.subplots(1,1, figsize = (mpl_utils.a4size[0]*0.8, mpl_utils.a4size[0]*0.4))
  figs = [fig1,fig2]
  axes = [ax1,ax2]
  Hmat =[]
  for t in typelist.split():
    
    print('Radius hist for type: %s' % t)
    filteredFiles = filter( lambda fn: t in fn,filenames) 
    files = [h5files.open(fn, 'r+') for fn in filteredFiles]
    if(len(files)==0):
      continue
    
    groups_with_adaption_by_type = [f['adaption/vessels_after_adaption'] for f in files]
    groups_without_adaption_by_type = [f['adaption/recomputed'] for f in files]
    
    for (groups,name,ax) in zip([groups_with_adaption_by_type,groups_without_adaption_by_type],['Adaption','no_Adaption'],axes):
      destination_group = f_measure.require_group('histograms_by_RC/%s/%s'%(name,t))      
      h,bin_edges = generateRadiusHistogram(dataman,groups,destination_group,'radius_%s' % t)
      #Hmat.append(h)
  
    #ax.step(h,bin_edges, where='post')
      ax.step(h,bin_edges, where='post', label = type_to_label_dict[t])
      #ax.plot((bin_edges[1:]+bin_edges[:-1])/2, h, label = type_to_label_dict[t])
      ax.set(xlabel='r [$\mu m$]', ylabel='p', title = 'Radius Histogram -- %s ' % (name))
  fontP = FontProperties()
  fontP.set_size('x-small')
  legend1 = ax1.legend( shadow=True,prop=fontP)
  legend2 = ax2.legend( shadow=True,prop=fontP)
  fig1.subplots_adjust(bottom=0.2)
  fig2.subplots_adjust(bottom=0.2)
  
  pdfpages.savefig(fig1, postfix='_radiushisto')
  pdfpages.savefig(fig2, postfix='_radiushisto')
    
  
def PlotRadiusHistogram_with_cache(dataman, f_measure, filenames, options, pdfpages):
  print("Ploting radius histogram!")
  filenames = adaption.get_files_with_successful_adaption(filenames)
  files = [h5files.open(fn, 'r') for fn in filenames]
  groups_with_adaption = [f['adaption/vessels_after_adaption'] for f in files]
  groups_without_adaption = [f['adaption/recomputed'] for f in files]
  destination_group = f_measure.require_group('histograms')
  
  
  
  for (groups,name) in zip([groups_without_adaption,groups_with_adaption],['no_Adaption','Adaption']):
    h,bin_edges = generateRadiusHistogram(dataman,groups,destination_group,'radius_%s' % name)
  
    fig, ax = pyplot.subplots(1,1, figsize = (mpl_utils.a4size[0]*0.8, mpl_utils.a4size[0]*0.4))
    #ax.step(h,bin_edges, where='post')
    ax.step(h,bin_edges)
    ax.set(xlabel='r [$\mu m$]', ylabel='p', title = 'Radius Histogram -- %s' % name)
    fig.subplots_adjust(bottom=0.2)
    pdfpages.savefig(fig, postfix='_radiushisto')

def DoIt(filenames, options):
  fn_measure = basename(commonprefix(filenames))
  fn_measure = myutils.strip_from_end(fn_measure, '.h5')
  fn_measure = myutils.strip_from_end(fn_measure, '-type')

  f_measure = h5files.open('adaption_common.h5', 'a', search = False)
  
  files = [h5files.open(fn, 'r') for fn in filenames]
  
  groups_with_adaption = [f['/vessels_after_adaption'] for f in files]
  #this data is stored with the adaption
  files_without_adaption = []
  for afile in files:
    completeFilename=str(afile['/vessels_after_adaption/parameters'].attrs['cwd'])+'/'+str(afile['/vessels_after_adaption/parameters'].attrs['vesselFileName'])
    files_without_adaption.append(
        h5files.open(completeFilename))
    
  groups_without_adaption = [f['vessels'] for f in files_without_adaption]

  with mpl_utils.PdfWriter('adaption_' + fn_measure+'.pdf') as pdfpages:
    import analyzeGeneral
    dataman = myutils.DataManager(20, [ analyzeGeneral.DataBasicVessel(), analyzeGeneral.DataVesselSamples(), analyzeGeneral.DataVesselGlobal()])
#    vesselgroups_without = groups_without_adaption
#    vesselgroups_with = groups_with_adaption
    
    geometric_data_before = getGeometricData(groups_without_adaption)
    perfusion_data_before = getTotalPerfusion(groups_without_adaption)*60
    
    geometric_data_after = getGeometricData(groups_with_adaption)
    perfusion_data_after = getTotalPerfusion(groups_with_adaption)*60
    
    if 1:
      res_without = getMultiScatter(300. * len(filenames), groups_without_adaption)
      plotMultiScatterBeauty(res_without, pdfpages)
      res_with = getMultiScatter(300. * len(filenames), groups_with_adaption)
      plotMultiScatterBeauty(res_with, pdfpages, withRadiusPressureRelation=False)
    
#    if 0:
#      PlotRadiusHistogram2(dataman, groups_without_adaption, pdfpages)
#      PlotRadiusHistogram2(dataman, groups_with_adaption, pdfpages)
#  
    if 0:
      #printBarPlot_rBV(filenames, pdfpages) -->alt
      printBarPlot_rBV(dataman, f_measure, filenames, options, pdfpages) #-->mw enginered
      
    if 0:  
      printBarPlot_rBF(dataman, f_measure, filenames, options, pdfpages)
    
    if 0:
      PlotRadiusHistogram_with_cache_by_RC(dataman, f_measure, filenames, options, pdfpages)
      #PlotRadiusHistogram_with_cache(dataman, f_measure, filenames, options, pdfpages)
    
    if 0:
      printBarPlot_vesseltype_on_root_node_configuration(filenames, pdfpages)
    
    ## is it worth to update this ???
    if 1:
      #importing the murray thing
      #DoGetMurray(filenames, pdfpages)
      printMurray(dataman, f_measure, filenames, options, pdfpages)
    if 1:      
      printMurray_alphas_effective(dataman, f_measure, filenames, options, pdfpages)
    if 1:
      text = ["Before adaption"]
      text2 = FormatGeometricAndPerfusionData(geometric_data_before, perfusion_data_before)
      text = text+text2
      def cachelocation(g):
        path = posixpath.join('FileCS_'+myutils.checksum(basename(g.file.filename)), g.name.strip(posixpath.sep))
        return (f_measure, path)

      prop_list2 = ['shearforce', 'velocity']        
      for name in prop_list2:
        data = []
        for gvessels in groups_without_adaption:
          data.append(dataman.obtain_data('basic_vessel_global', name, gvessels, cachelocation(gvessels)))
        text.append(r'$<%s>$ = $%s$%s' %
          (Prettyfier.get_sym(name), Format(name, data), Prettyfier.get_munit(name)))
    
      fig, _ = mpl_utils.MakeTextPage(text,figsize = (mpl_utils.a4size[0]*0.8, mpl_utils.a4size[0]*0.8))
      pdfpages.savefig(fig, postfix='_vesselsglobal')

      text = FormatParameters(groups_without_adaption[0].file)
      fig, _ = mpl_utils.MakeTextPage(text,figsize = (mpl_utils.a4size[0]*0.8, mpl_utils.a4size[0]*0.8))
      pdfpages.savefig(fig, postfix='_vesselsparams')
    

      text = ["After adaption"]
      text2 = FormatGeometricAndPerfusionData(geometric_data_after, perfusion_data_after)
      text = text+text2
      def cachelocation(g):
        path = posixpath.join('FileCS_'+myutils.checksum(basename(g.file.filename)), g.name.strip(posixpath.sep))
        return (f_measure, path)

      prop_list2 = ['shearforce', 'velocity']        
      for name in prop_list2:
        data = []
        for gvessels in groups_with_adaption:
          data.append(dataman.obtain_data('basic_vessel_global', name, gvessels, cachelocation(gvessels)))
        text.append(r'$<%s>$ = $%s$%s' %
          (Prettyfier.get_sym(name), Format(name, data), Prettyfier.get_munit(name)))
    
      fig, _ = mpl_utils.MakeTextPage(text,figsize = (mpl_utils.a4size[0]*0.8, mpl_utils.a4size[0]*0.8))
      pdfpages.savefig(fig, postfix='_vesselsglobal')

      text = FormatParameters(groups_without_adaption[0].file)
      fig, _ = mpl_utils.MakeTextPage(text,figsize = (mpl_utils.a4size[0]*0.8, mpl_utils.a4size[0]*0.8))
      pdfpages.savefig(fig, postfix='_vesselsparams')
    
        
    if 0 and all(map(lambda g: 'data' in g.parent, groups_without_adaption)):
      data = VesselData()
      for g in groups_without_adaption:
        data.add(g.parent['data'])
      plot_topological_stats_avg(data, pdfpages)
      
    if 1:
      #PlotQdevs_unnormalized(filenames, pdfpages)
      PlotQdevs_normalized(filenames, options, pdfpages)

if __name__ == "__main__":
  import optparse  #Note: Deprecated since version 2.7. Use argparse instead
  parser = optparse.OptionParser()
  parser.add_option("-O","--with-o2", dest="with_o2", help="look at detailed o2 data", default=False, action="store_true")
  parser.add_option("-T","--only_two_root", dest="two", help="flag to change the considered types", default=False, action="store_true")  
  parser.add_option("-a","--with_all_types", dest="all_types", help="take all types",default=False, action="store_true")  
  parser.add_option("-s","--singel_type", dest="single", help="", default=False, action="store_true")  
  options, args = parser.parse_args()

  #filenames, pattern = args[:-1], args[-1]
  filenames = args[:]
  DoIt(filenames, options)

