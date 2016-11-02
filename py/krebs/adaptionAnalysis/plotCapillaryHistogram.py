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
if __name__ == '__main__':
  import os.path, sys
  sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../..'))
  sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'..'))
from os.path import basename, commonprefix
import time
import krebsutils
import h5py
import h5files
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
#pyplot.style.use('ggplot')

from quantities import Prettyfier

from plotVessels import *
from plotVesselsAdaption import getVesselTypes
from analyzeMurray import *

from scipy.optimize import fsolve
from scipy.optimize import minimize_scalar


type_to_label_dict={"typeA-": "RC1","typeB-": "RC2","typeC-": "RC3","typeD-": "RC4","typeE-": "RC5","typeF-": "RC6","typeG-": "RC7","typeH-": "RC8","typeI-": "RC9"}
vesselTypeMarkers = '<>*osd^+x'
colors = [
     '#e50000', #red
     '#9a0eea', #violet
     '#0343df', #blue
     '#f97306', #orange
     '#677a04', #olive green
     '#ceb301', #mustard
     '#04d8b2', #aquamarine
     '#06470c', #forest greeen
     '#840000', #dark red
     '#607c8e', #blue grey
  ]

double_log = True

##### RADIUS
@myutils.UsesDataManager
def get_capillary_radius(datamanager, destination_group, f):
  datanames = 'nA_capillary_radius yA_capillary_radius'.split()
  def process(grp):
    veins, arteries, capillaries = getVesselTypes(grp)
    edges, radii = krebsutils.read_vessels_from_hdf(grp,['radius'])
    return radii[capillaries]
    
  def read(gmeasure, groupname):
    gmeasure = gmeasure[groupname]
    return [gmeasure[name][()] for name in datanames ]
    
  def write(gmeasure, groupname):
    gmeasure=gmeasure.create_group(groupname)
    group_without_adaption = f['adaption/recomputed']
    group_with_adaption = f['adaption/vessels_after_adaption']
    for (group,name) in zip([group_without_adaption,group_with_adaption],datanames):
      radii = process(group)
      gmeasure.create_dataset(name, data=radii)
  ret = myutils.hdf_data_caching(read, write, destination_group, (f.filename), (1, ))
  return ret
@myutils.UsesDataManager
def generate_capillary_hist_radius(dataman, inputfiles, destination_group, destination_name):
  def process(inputfiles):
    #groups_with_adaption = [f['adaption/vessels_after_adaption'] for f in files]
    #groups_without_adaption = [f['adaption/recomputed'] for f in files]
    
    all_capillary_radii_nA = []
    all_capillary_radii_yA = []
    for f in inputfiles:
      nA_radii, yA_radii = get_capillary_radius(dataman,destination_group,f)
      all_capillary_radii_nA = np.hstack((nA_radii,all_capillary_radii_nA))
      all_capillary_radii_yA = np.hstack((yA_radii,all_capillary_radii_yA))
    return (all_capillary_radii_nA,all_capillary_radii_yA)
  def write(gmeasure,groupname):
    gmeasure = gmeasure.create_group(groupname)
    x_edges_nA = np.linspace(2,8,50)
    x_edges_yA = np.linspace(2,8,50)
    nA_radii,yA_radii, = process(inputfiles)
    
    h1, x_edges_nA = np.histogram(nA_radii,bins=x_edges_nA,density=True)
    h2, x_edges_yA = np.histogram(yA_radii,bins=x_edges_yA,density=True)
    gmeasure.create_dataset('h_nA_radii', data=h1)
    gmeasure.create_dataset('h_yA_radii', data=h2)
    gmeasure.create_dataset('x_edges_nA_radii', data= x_edges_nA)
    gmeasure.create_dataset('x_edges_yA_radii', data= x_edges_yA)
  def read(gmeasure, groupname):
    gmeasure = gmeasure[groupname]
    r1 = gmeasure['h_nA_radii']
    r2 = gmeasure['h_yA_radii']
    r3 = gmeasure['x_edges_nA_radii']
    r4 = gmeasure['x_edges_yA_radii']
    return (r1,r2,r3,r4)
  '''
  cache number needs to be higher than other murray stuff
  since this will open the cache file for the second time
  '''
  ret = myutils.hdf_data_caching(read, write, destination_group, (destination_name,), (2,))
  return ret

def plotCapillaryRadiusHistogram(dataman, f_measure, filenames, options, pdfpages):
    filenames = adaption.get_files_with_successful_adaption(filenames)
    files = [h5files.open(fn, 'r') for fn in filenames]
    destination_group = f_measure.require_group('cap_radius')
    if(options.two):
      typelist = 'typeD- typeE- typeG- typeH-'
    else:
      typelist = 'typeA- typeB- typeC- typeF- typeI-'
    if(options.all_types):
      typelist = 'typeA- typeB- typeC- typeD- typeE- typeF- typeG- typeH- typeI-'
      
    if(options.single):
      typelist = 'typeF- '
      
    fig = plt.figure()
    #fig.suptitle('cap Radii', fontsize=12)
    #ax1 = fig.add_subplot(111)
    ax2 = fig.add_subplot(111)
    #ax1.set_title('M- network')
    #ax1.set_xlabel(r'radius/$\mu m$')
    #ax1.set_ylabel(r'$\rho$')
    ax2.set_title('Adapted- network')
    ax2.set_xlabel(r'radius/$\mu m$')
    ax2.set_ylabel(r'$\rho$')
    
    for (i,t) in enumerate(typelist.split()):
      print('Capillary radius for type: %s' % t)
      filteredFiles = filter( lambda fn: t in fn,filenames) 
      files = [h5files.open(fn, 'r+') for fn in filteredFiles]
  
      if(len(files)==0):#that means no file of dedicated type is left after filter
        continue        
        #h_nA, h_yA, x_edges_nA, x_edges_yA = generate_capillary_hist(dataman, files, destination_group, t[:-1])
      #data = generate_murray_data_plot3(dataman, files, destination_group, t[:-1])      
      else:#do the shit
        h_nA,h_yA,xedges_nA,xedges_yA = generate_capillary_hist_radius(dataman, files, destination_group, t[:-1] + '_cap')
        xedges_nA = np.asarray(xedges_nA)
        xedges_yA = np.asarray(xedges_yA)
        h1 =np.asarray(h_nA)
        h2 =np.asarray(h_yA)
        #ax1.plot((xedges_nA[1:]+xedges_nA[:-1])/2, h1,  label = type_to_label_dict[t], marker= vesselTypeMarkers[i], color=colors[i])
        ax2.plot((xedges_yA[1:]+xedges_yA[:-1])/2, h2,  label = type_to_label_dict[t], marker= vesselTypeMarkers[i], color=colors[i])
          
    fontP = FontProperties()
    fontP.set_size('x-small')
    ax2.legend(prop=fontP)
    ax2.semilogy()
    #ax1.semilogy()
    #ax1.text(-0.1, 1.15, 'A', transform=ax1.transAxes,fontsize=16, fontweight='bold', va='top', ha='right')
    ax2.text(-0.1, 1.15, 'B', transform=ax2.transAxes,fontsize=16, fontweight='bold', va='top', ha='right')
    plt.tight_layout()
    pdfpages.savefig(fig)

##### FLOW
@myutils.UsesDataManager
def get_capillary_flow(datamanager, destination_group, f):
  datanames = 'nA_capillary_flow yA_capillary_flow'.split()
  def process(grp):
    veins, arteries, capillaries = getVesselTypes(grp)
    edges, flows = krebsutils.read_vessels_from_hdf(grp,['flow'])
    return flows[capillaries]
    
  def read(gmeasure, groupname):
    gmeasure = gmeasure[groupname]
    return [gmeasure[name][()] for name in datanames ]
    
  def write(gmeasure, groupname):
    gmeasure=gmeasure.create_group(groupname)
    group_without_adaption = f['adaption/recomputed']
    group_with_adaption = f['adaption/vessels_after_adaption']
    for (group,name) in zip([group_without_adaption,group_with_adaption],datanames):
      flows = process(group)
      gmeasure.create_dataset(name, data=flows)
  ret = myutils.hdf_data_caching(read, write, destination_group, (f.filename), (1, ))
  return ret
@myutils.UsesDataManager
def generate_capillary_hist(dataman, inputfiles, destination_group, destination_name):
  def process(inputfiles):
    #groups_with_adaption = [f['adaption/vessels_after_adaption'] for f in files]
    #groups_without_adaption = [f['adaption/recomputed'] for f in files]
    
    all_capillary_flows_nA = []
    all_capillary_flows_yA = []
    for f in inputfiles:
      nA_flow, yA_flow = get_capillary_flow(dataman,destination_group,f)
      all_capillary_flows_nA = np.hstack((nA_flow,all_capillary_flows_nA))
      all_capillary_flows_yA = np.hstack((yA_flow,all_capillary_flows_yA))
    return (all_capillary_flows_nA,all_capillary_flows_yA)
  def write(gmeasure,groupname):
    gmeasure = gmeasure.create_group(groupname)
    x_edges_nA = np.logspace(1,5.5,100)
    x_edges_yA = np.logspace(1,5.5,100)
    nA_flows,yA_flows, = process(inputfiles)
    print(nA_flows.shape)
    h1, x_edges_nA = np.histogram(nA_flows,bins=x_edges_nA,density=True)
    h2, x_edges_yA = np.histogram(yA_flows,bins=x_edges_yA,density=True)
    gmeasure.create_dataset('h_nA', data=h1)
    gmeasure.create_dataset('h_yA', data=h2)
    gmeasure.create_dataset('x_edges_nA', data= x_edges_nA)
    gmeasure.create_dataset('x_edges_yA', data= x_edges_yA)
  def read(gmeasure, groupname):
    gmeasure = gmeasure[groupname]
    r1 = gmeasure['h_nA']
    r2 = gmeasure['h_yA']
    r3 = gmeasure['x_edges_nA']
    r4 = gmeasure['x_edges_yA']
    return (r1,r2,r3,r4)
  '''
  cache number needs to be higher than other murray stuff
  since this will open the cache file for the second time
  '''
  ret = myutils.hdf_data_caching(read, write, destination_group, (destination_name,), (2,))
  return ret

def plotCapillaryFlowHistogram(dataman, f_measure, filenames, options, pdfpages):
    filenames = adaption.get_files_with_successful_adaption(filenames)
    files = [h5files.open(fn, 'r') for fn in filenames]
    destination_group = f_measure.require_group('cap_flow')
    if(options.two):
      typelist = 'typeD- typeE- typeG- typeH-'
    else:
      typelist = 'typeA- typeB- typeC- typeF- typeI-'
    if(options.all_types):
      typelist = 'typeA- typeB- typeC- typeD- typeE- typeF- typeG- typeH- typeI-'
      
    if(options.single):
      typelist = 'typeF- '
      
    fig = plt.figure()
    #fig.suptitle('cap Flows', fontsize=12)
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)
    ax1.set_title('M- network')
    ax1.set_xlabel(r'flow/$\frac{\mu m^3}{s}$')
    ax1.set_ylabel(r'$\rho$')
    ax2.set_title('Adapted- network')
    ax2.set_xlabel(r'flow/$\frac{\mu m^3}{s}$')
    ax2.set_ylabel(r'$\rho$')
    
    for (i,t) in enumerate(typelist.split()):
      print('Capillary flow for type: %s' % t)
      filteredFiles = filter( lambda fn: t in fn,filenames) 
      files = [h5files.open(fn, 'r+') for fn in filteredFiles]
  
      if(len(files)==0):#that means no file of dedicated type is left after filter
        continue        
        #h_nA, h_yA, x_edges_nA, x_edges_yA = generate_capillary_hist(dataman, files, destination_group, t[:-1])
      #data = generate_murray_data_plot3(dataman, files, destination_group, t[:-1])      
      else:
        h_nA,h_yA,xedges_nA,xedges_yA = generate_capillary_hist(dataman, files, destination_group, t[:-1] + '_cap')
        xedges_nA = np.asarray(xedges_nA)
        xedges_yA = np.asarray(xedges_yA)
        h1 =np.asarray(h_nA)
        h2 =np.asarray(h_yA)
        ax1.plot((xedges_nA[1:]+xedges_nA[:-1])/2, h1,  label = type_to_label_dict[t], marker= vesselTypeMarkers[i], color=colors[i])
        ax2.plot((xedges_yA[1:]+xedges_yA[:-1])/2, h2,  label = type_to_label_dict[t], marker= vesselTypeMarkers[i], color=colors[i])
          
    fontP = FontProperties()
    fontP.set_size('x-small')
    ax1.legend(prop=fontP)
    if double_log:
      ax2.loglog()
      ax1.loglog()
    else:
      ax2.semilogx()
      ax1.semilogx()
    ax1.text(-0.1, 1.15, 'A', transform=ax1.transAxes,fontsize=16, fontweight='bold', va='top', ha='right')
    ax2.text(-0.1, 1.15, 'B', transform=ax2.transAxes,fontsize=16, fontweight='bold', va='top', ha='right')
    plt.tight_layout()
    pdfpages.savefig(fig)
    
@myutils.UsesDataManager
def get_radius(datamanager, destination_group, f):
  datanames = 'nA_radius yA_radius'.split()
  def process(grp):
    #veins, arteries, capillaries = getVesselTypes(grp)
    edges, radii = krebsutils.read_vessels_from_hdf(grp,['radius'])
    return radii
    
  def read(gmeasure, groupname):
    gmeasure = gmeasure[groupname]
    return [gmeasure[name][()] for name in datanames ]
    
  def write(gmeasure, groupname):
    gmeasure=gmeasure.create_group(groupname)
    group_without_adaption = f['adaption/recomputed']
    group_with_adaption = f['adaption/vessels_after_adaption']
    for (group,name) in zip([group_without_adaption,group_with_adaption],datanames):
      radii = process(group)
      gmeasure.create_dataset(name, data=radii)
  ret = myutils.hdf_data_caching(read, write, destination_group, (f.filename), (1, ))
  return ret
@myutils.UsesDataManager
def generate_hist_radius(dataman, inputfiles, destination_group, destination_name):
  def process(inputfiles):
    #groups_with_adaption = [f['adaption/vessels_after_adaption'] for f in files]
    #groups_without_adaption = [f['adaption/recomputed'] for f in files]
    
    all_radii_nA = []
    all_radii_yA = []
    for f in inputfiles:
      nA_radii, yA_radii = get_radius(dataman,destination_group,f)
      all_radii_nA = np.hstack((nA_radii,all_radii_nA))
      all_radii_yA = np.hstack((yA_radii,all_radii_yA))
    return (all_radii_nA,all_radii_yA)
  def write(gmeasure,groupname):
    gmeasure = gmeasure.create_group(groupname)
    x_edges_nA = np.linspace(2,8,100)
    x_edges_yA = np.linspace(2,8,100)
    nA_radii,yA_radii, = process(inputfiles)
    
    h1, x_edges_nA = np.histogram(nA_radii,bins=x_edges_nA, density=True)
    h2, x_edges_yA = np.histogram(yA_radii,bins=x_edges_yA, density=True)
    gmeasure.create_dataset('h_nA_radii', data=h1)
    gmeasure.create_dataset('h_yA_radii', data=h2)
    gmeasure.create_dataset('x_edges_nA_radii', data= x_edges_nA)
    gmeasure.create_dataset('x_edges_yA_radii', data= x_edges_yA)
  def read(gmeasure, groupname):
    gmeasure = gmeasure[groupname]
    r1 = gmeasure['h_nA_radii']
    r2 = gmeasure['h_yA_radii']
    r3 = gmeasure['x_edges_nA_radii']
    r4 = gmeasure['x_edges_yA_radii']
    return (r1,r2,r3,r4)
  '''
  cache number needs to be higher than other murray stuff
  since this will open the cache file for the second time
  '''
  ret = myutils.hdf_data_caching(read, write, destination_group, (destination_name,), (2,))
  return ret

def plotRadiusHistogram(dataman, f_measure, filenames, options, pdfpages):
    filenames = adaption.get_files_with_successful_adaption(filenames)
    files = [h5files.open(fn, 'r') for fn in filenames]
    destination_group = f_measure.require_group('allRadius')
    if(options.two):
      typelist = 'typeD- typeE- typeG- typeH-'
    else:
      typelist = 'typeA- typeB- typeC- typeF- typeI-'
    if(options.all_types):
      typelist = 'typeA- typeB- typeC- typeD- typeE- typeF- typeG- typeH- typeI-'
      
    if(options.single):
      typelist = 'typeF- '
      
    fig = plt.figure()
    #fig.suptitle('all Radii', fontsize=12)
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)
    ax1.set_title('M- network')
    ax1.set_xlabel(r'radius/$\mu m$')
    ax1.set_ylabel(r'$\rho$')
    ax2.set_title('Adapted- network')
    ax2.set_xlabel(r'radius/$\mu m$')
    ax2.set_ylabel(r'$\rho$')
    
    for (i,t) in enumerate(typelist.split()):
      print('Capillary radius for type: %s' % t)
      filteredFiles = filter( lambda fn: t in fn,filenames) 
      files = [h5files.open(fn, 'r+') for fn in filteredFiles]
  
      if(len(files)==0):#that means no file of dedicated type is left after filter
        continue        
        #h_nA, h_yA, x_edges_nA, x_edges_yA = generate_capillary_hist(dataman, files, destination_group, t[:-1])
      #data = generate_murray_data_plot3(dataman, files, destination_group, t[:-1])      
      else:#do the shit
        h_nA,h_yA,xedges_nA,xedges_yA = generate_hist_radius(dataman, files, destination_group, t[:-1] + '_allRadii')
        xedges_nA = np.asarray(xedges_nA)
        xedges_yA = np.asarray(xedges_yA)
        h1 =np.asarray(h_nA)
        h2 =np.asarray(h_yA)
        ax1.plot((xedges_nA[1:]+xedges_nA[:-1])/2, h1,  label = type_to_label_dict[t], marker= vesselTypeMarkers[i], color=colors[i])
        ax2.plot((xedges_yA[1:]+xedges_yA[:-1])/2, h2,  label = type_to_label_dict[t], marker= vesselTypeMarkers[i], color=colors[i])
          
    fontP = FontProperties()
    fontP.set_size('x-small')
    ax1.legend(prop=fontP)
    ax2.semilogy()
    ax1.semilogy()
    ax1.text(-0.1, 1.15, 'A', transform=ax1.transAxes,fontsize=16, fontweight='bold', va='top', ha='right')
    ax2.text(-0.1, 1.15, 'B', transform=ax2.transAxes,fontsize=16, fontweight='bold', va='top', ha='right')
    plt.tight_layout()
    pdfpages.savefig(fig)

##### FLOW
@myutils.UsesDataManager
def get_flow(datamanager, destination_group, f):
  datanames = 'nA_flow yA_flow'.split()
  def process(grp):
    #veins, arteries, capillaries = getVesselTypes(grp)
    edges, flows = krebsutils.read_vessels_from_hdf(grp,['flow'])
    #capillaries = np.asarray(capillaries)
    return flows
    
  def read(gmeasure, groupname):
    gmeasure = gmeasure[groupname]
    return [gmeasure[name][()] for name in datanames ]
    
  def write(gmeasure, groupname):
    gmeasure=gmeasure.create_group(groupname)
    group_without_adaption = f['adaption/recomputed']
    group_with_adaption = f['adaption/vessels_after_adaption']
    for (group,name) in zip([group_without_adaption,group_with_adaption],datanames):
      flows = process(group)
      gmeasure.create_dataset(name, data=flows)
  ret = myutils.hdf_data_caching(read, write, destination_group, (f.filename), (1, ))
  return ret
@myutils.UsesDataManager
def generate_flow_hist(dataman, inputfiles, destination_group, destination_name):
  def process(inputfiles):
    #groups_with_adaption = [f['adaption/vessels_after_adaption'] for f in files]
    #groups_without_adaption = [f['adaption/recomputed'] for f in files]
    
    all_flows_nA = []
    all_flows_yA = []
    for f in inputfiles:
      nA_flow, yA_flow = get_flow(dataman,destination_group,f)
      all_flows_nA = np.hstack((nA_flow,all_flows_nA))
      all_flows_yA = np.hstack((yA_flow,all_flows_yA))
    return (all_flows_nA,all_flows_yA)
  def write(gmeasure,groupname):
    gmeasure = gmeasure.create_group(groupname)
    x_edges_nA = np.logspace(0.5,6.5,100)
    x_edges_yA = np.logspace(0.5,6.5,100)
    nA_flows,yA_flows, = process(inputfiles)
    print(nA_flows.shape)
    h1, x_edges_nA = np.histogram(nA_flows,bins=x_edges_nA, density=True)
    h2, x_edges_yA = np.histogram(yA_flows,bins=x_edges_yA, density=True)
    gmeasure.create_dataset('h_nA', data=h1)
    gmeasure.create_dataset('h_yA', data=h2)
    gmeasure.create_dataset('x_edges_nA', data= x_edges_nA)
    gmeasure.create_dataset('x_edges_yA', data= x_edges_yA)
  def read(gmeasure, groupname):
    gmeasure = gmeasure[groupname]
    r1 = gmeasure['h_nA']
    r2 = gmeasure['h_yA']
    r3 = gmeasure['x_edges_nA']
    r4 = gmeasure['x_edges_yA']
    return (r1,r2,r3,r4)
  '''
  cache number needs to be higher than other murray stuff
  since this will open the cache file for the second time
  '''
  ret = myutils.hdf_data_caching(read, write, destination_group, (destination_name,), (2,))
  return ret

def plotFlowHistogram(dataman, f_measure, filenames, options, pdfpages):
    filenames = adaption.get_files_with_successful_adaption(filenames)
    files = [h5files.open(fn, 'r') for fn in filenames]
    destination_group = f_measure.require_group('allFlow')
    if(options.two):
      typelist = 'typeD- typeE- typeG- typeH-'
    else:
      typelist = 'typeA- typeB- typeC- typeF- typeI-'
    if(options.all_types):
      typelist = 'typeA- typeB- typeC- typeD- typeE- typeF- typeG- typeH- typeI-'
      
    if(options.single):
      typelist = 'typeF- '
      
    fig = plt.figure()
    #fig.suptitle('all Flows', fontsize=12)
    #fig.set_title('all Flows')
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)
    ax1.set_title('M- network')
    ax1.set_xlabel(r'flow/$\frac{\mu m^3}{s}$')
    ax1.set_ylabel(r'$\rho$')
    ax1.text(-0.1, 1.15, 'A', transform=ax1.transAxes,fontsize=16, fontweight='bold', va='top', ha='right')
    ax2.set_title('Adapted- network')
    ax2.set_xlabel(r'flow/$\frac{\mu m^3}{s}$')
    ax2.set_ylabel(r'$\rho$')
    ax2.text(-0.1, 1.15, 'B', transform=ax2.transAxes,fontsize=16, fontweight='bold', va='top', ha='right')
    
    for (i,t) in enumerate(typelist.split()):
      print('Capillary flow for type: %s' % t)
      filteredFiles = filter( lambda fn: t in fn,filenames) 
      files = [h5files.open(fn, 'r+') for fn in filteredFiles]
  
      if(len(files)==0):#that means no file of dedicated type is left after filter
        continue        
        #h_nA, h_yA, x_edges_nA, x_edges_yA = generate_capillary_hist(dataman, files, destination_group, t[:-1])
      #data = generate_murray_data_plot3(dataman, files, destination_group, t[:-1])      
      else:
        h_nA,h_yA,xedges_nA,xedges_yA = generate_flow_hist(dataman, files, destination_group, t[:-1] + '_allFlow')
        xedges_nA = np.asarray(xedges_nA)
        xedges_yA = np.asarray(xedges_yA)
        h1 =np.asarray(h_nA)
        h2 =np.asarray(h_yA)
        ax1.plot((xedges_nA[1:]+xedges_nA[:-1])/2, h1,  label = type_to_label_dict[t], marker= vesselTypeMarkers[i], color=colors[i])
        ax2.plot((xedges_yA[1:]+xedges_yA[:-1])/2, h2,  label = type_to_label_dict[t], marker= vesselTypeMarkers[i], color=colors[i])
          
    fontP = FontProperties()
    fontP.set_size('x-small')
    ax1.legend(prop=fontP)
    if double_log:
      ax2.loglog()
      ax1.loglog()
    else:
      ax2.semilogx()
      ax1.semilogx()
    
    plt.tight_layout()
    pdfpages.savefig(fig)

def DoIt(filenames, options):
  fn_measure = basename(commonprefix(filenames))
  fn_measure = myutils.strip_from_end(fn_measure, '.h5')
  fn_measure = myutils.strip_from_end(fn_measure, '-type')

  f_measure = h5files.open('adaption_common.h5', 'a', search = False)
  
  files = [h5files.open(fn, 'r') for fn in filenames]
  groups_without_adaption = [f['/adaption/recomputed'] for f in files]
  
  groups_with_adaption = [f['/adaption/vessels_after_adaption'] for f in files]


  with mpl_utils.PdfWriter(fn_measure+'caps.pdf') as pdfpages:
    import analyzeGeneral
    dataman = myutils.DataManager(20, [ analyzeGeneral.DataBasicVessel(), analyzeGeneral.DataVesselSamples(), analyzeGeneral.DataVesselGlobal()])
#    vesselgroups_without = groups_without_adaption
#    vesselgroups_with = groups_with_adaption
    
#    geometric_data_before = getGeometricData(groups_without_adaption)
#    perfusion_data_before = getTotalPerfusion(groups_without_adaption)*60
    
#    geometric_data_after = getGeometricData(groups_with_adaption)
#    perfusion_data_after = getTotalPerfusion(groups_with_adaption)*60
    if 1:
      plotFlowHistogram(dataman, f_measure, filenames, options, pdfpages)
      
    if 1:
      plotRadiusHistogram(dataman, f_measure, filenames, options, pdfpages)    
    
    if 1:
      plotCapillaryFlowHistogram(dataman, f_measure, filenames, options, pdfpages)
      
    if 1:
      plotCapillaryRadiusHistogram(dataman, f_measure, filenames, options, pdfpages)

def DoIt_single(filenames, options):
  fn_measure = basename(commonprefix(filenames))
  fn_measure = myutils.strip_from_end(fn_measure, '.h5')
  fn_measure = myutils.strip_from_end(fn_measure, '-type')

  f_measure = h5files.open('adaption_common.h5', 'a', search = False)
  
  files = [h5files.open(fn, 'r') for fn in filenames]
  groups_without_adaption = [f['/adaption/recomputed'] for f in files]
  
  groups_with_adaption = [f['/adaption/vessels_after_adaption'] for f in files]

  import analyzeGeneral
  dataman = myutils.DataManager(20, [ analyzeGeneral.DataBasicVessel(), analyzeGeneral.DataVesselSamples(), analyzeGeneral.DataVesselGlobal()])
  with mpl_utils.PdfWriter(fn_measure+'flow_hist.pdf') as pdfpages:
    if 1:
      plotFlowHistogram(dataman, f_measure, filenames, options, pdfpages)
  with mpl_utils.PdfWriter(fn_measure+'radius_hist.pdf') as pdfpages:
    if 1:
      plotRadiusHistogram(dataman, f_measure, filenames, options, pdfpages)    
  with mpl_utils.PdfWriter(fn_measure+'cap_flow_hist.pdf') as pdfpages:
    if 1:
      plotCapillaryFlowHistogram(dataman, f_measure, filenames, options, pdfpages)
  with mpl_utils.PdfWriter(fn_measure+'cap_radius_hist.pdf') as pdfpages:    
    if 1:
      plotCapillaryRadiusHistogram(dataman, f_measure, filenames, options, pdfpages)


if __name__ == "__main__":
  import optparse  #Note: Deprecated since version 2.7. Use argparse instead
  parser = optparse.OptionParser()
  parser.add_option("-O","--with-o2", dest="with_o2", help="look at detailed o2 data", default=False, action="store_true")
  parser.add_option("-T","--only_two_root", dest="two", help="flag to change the considered types", default=False, action="store_true")  
  parser.add_option("-a","--with_all_types", dest="all_types", help="take all types",default=True, action="store_true")  
  parser.add_option("-s","--singel_type", dest="single", help="", default=False, action="store_true")  
  parser.add_option("-p","--publication_mode", dest="is_for_publication", help="", default=False, action="store_true")  
  options, args = parser.parse_args()

  #filenames, pattern = args[:-1], args[-1]
  filenames = args[:]
  if options.is_for_publication:
    DoIt_single(filenames, options)
  else:
    DoIt(filenames, options)