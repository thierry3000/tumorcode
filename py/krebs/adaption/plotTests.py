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

import numpy as np
import krebsutils as ku
import math
import collections
import itertools
import myutils
import os

import krebs.adaption as adap
import krebsjobs.parameters.parameterSetsAdaption

from matplotlib.backends.backend_pdf import PdfPages
import h5py
import matplotlib.pyplot as plt

def plot_pressure_relation():
    print("later")
def plot_asymetric(fn,pp):
    adaption_grp = fn['Asymetric/vessels_after_adaption/vessels_after_adaption']
    index_of_artery = np.bitwise_and(np.asarray(adaption_grp['edges/flags']), ku.ARTERY)>0
    index_of_capillary = np.bitwise_and(np.asarray(adaption_grp['edges/flags']), ku.CAPILLARY)>0
    index_of_vein = np.bitwise_and(np.asarray(adaption_grp['edges/flags']), ku.VEIN)>0
    
    shearforce = adaption_grp['edges/shearforce']
    shearforce = np.multiply(shearforce,10000)
    diameter = np.multiply(adaption_grp['edges/radius'],2)
    node_a_index =adaption_grp['edges/node_a_index']
    node_b_index =adaption_grp['edges/node_b_index']
        
    pressure_at_nodes = adaption_grp['nodes/pressure']
    
    
    pressure_at_vessel=[]
        
    for NodeA, NodeB in itertools.izip(node_a_index,node_b_index):
        pressure_at_vessel.append((pressure_at_nodes[NodeA]+pressure_at_nodes[NodeB])/2 * 7.5)
    pressure_at_vessel = np.array(pressure_at_vessel)
    fig = plt.figure()
    plt.subplot(2,1,1)
    a = pressure_at_vessel[index_of_artery]
    b = shearforce[index_of_artery]
    plt.plot(a,b,'or')
    a = pressure_at_vessel[index_of_capillary]
    b = shearforce[index_of_capillary]
    plt.plot(a,b,'yD')
    a = pressure_at_vessel[index_of_vein]
    b = shearforce[index_of_vein]
    plt.plot(a,b,'bv')
#    plt.semilogy(pressure_at_vessel[index_of_artery],shearforce[index_of_artery],'ro')
#    plt.semilogy(pressure_at_vessel[index_of_capillary],shearforce[index_of_capillary],'yD')
#    plt.semilogy(pressure_at_vessel[index_of_vein],shearforce[index_of_vein],'bv')
    plt.grid()
    plt.yscale('log')
    plt.xlabel('pressure/ mmHg')
    plt.xlim([0,80])
    #plt.ylim([1,100])
    plt.ylabel('shearstress dyne/cm^2')
    
    plt.subplot(2,1,2)
    plt.plot(pressure_at_vessel[index_of_artery],diameter[index_of_artery],'ro')
    plt.plot(pressure_at_vessel[index_of_capillary],diameter[index_of_capillary],'yD')
    plt.plot(pressure_at_vessel[index_of_vein],diameter[index_of_vein],'bv')
    plt.grid()
    plt.legend(['ART','CAP','VEN'])
    plt.xlabel('pressure/ mmHg')
    plt.ylabel('diameter/ mum')
    plt.xlim([1,80])
    plt.ylim([5,30])
    
    pp.savefig()
    
    plt.savefig('All_signals_CONVERGENT_asymetric')
def plot_symetricB(fn,pp):
    adaption_grp = fn['symetricB/vessels_after_adaption/vessels_after_adaption']
    index_of_artery = np.bitwise_and(np.asarray(adaption_grp['edges/flags']), ku.ARTERY)>0
    index_of_capillary = np.bitwise_and(np.asarray(adaption_grp['edges/flags']), ku.CAPILLARY)>0
    index_of_vein = np.bitwise_and(np.asarray(adaption_grp['edges/flags']), ku.VEIN)>0
    
    shearforce = adaption_grp['edges/shearforce']
    shearforce = np.multiply(shearforce,10000)
    diameter = np.multiply(adaption_grp['edges/radius'],2)
    node_a_index =adaption_grp['edges/node_a_index']
    node_b_index =adaption_grp['edges/node_b_index']
        
    pressure_at_nodes = adaption_grp['nodes/pressure']
    
    
    pressure_at_vessel=[]
        
    for NodeA, NodeB in itertools.izip(node_a_index,node_b_index):
        pressure_at_vessel.append((pressure_at_nodes[NodeA]+pressure_at_nodes[NodeB])/2 * 7.5)
    pressure_at_vessel = np.array(pressure_at_vessel)
    fig = plt.figure()
    plt.subplot(2,1,1)
    a = pressure_at_vessel[index_of_artery]
    b = shearforce[index_of_artery]
    plt.plot(a,b,'or')
    a = pressure_at_vessel[index_of_capillary]
    b = shearforce[index_of_capillary]
    plt.plot(a,b,'yD')
    a = pressure_at_vessel[index_of_vein]
    b = shearforce[index_of_vein]
    plt.plot(a,b,'bv')
#    plt.semilogy(pressure_at_vessel[index_of_artery],shearforce[index_of_artery],'ro')
#    plt.semilogy(pressure_at_vessel[index_of_capillary],shearforce[index_of_capillary],'yD')
#    plt.semilogy(pressure_at_vessel[index_of_vein],shearforce[index_of_vein],'bv')
    plt.grid()
    plt.yscale('log')
    plt.xlabel('pressure/ mmHg')
    #plt.xlim([0,80])
    #plt.ylim([1,100])
    plt.ylabel('shearstress dyne/cm^2')
    
    plt.subplot(2,1,2)
    plt.plot(pressure_at_vessel[index_of_artery],diameter[index_of_artery],'ro')
    plt.plot(pressure_at_vessel[index_of_capillary],diameter[index_of_capillary],'yD')
    plt.plot(pressure_at_vessel[index_of_vein],diameter[index_of_vein],'bv')
    plt.grid()
#    plt.legend(['ART','CAP','VEN'],loc='upper center')
    plt.xlabel('pressure/ mmHg')
    plt.ylabel('diameter/ mum')
    #plt.xlim([1,80])
    #plt.ylim([5,30])
    
    #pp.savefig()
    
    plt.savefig('All_signals_CONVERGENT_symetricB')
def plot_symetricA(fn):
    adaption_grp = fn['symetricA/vessels_after_adaption/vessels_after_adaption']
    index_of_artery = np.bitwise_and(np.asarray(adaption_grp['edges/flags']), ku.ARTERY)>0
    index_of_capillary = np.bitwise_and(np.asarray(adaption_grp['edges/flags']), ku.CAPILLARY)>0
    index_of_vein = np.bitwise_and(np.asarray(adaption_grp['edges/flags']), ku.VEIN)>0
    
    shearforce = adaption_grp['edges/shearforce']
    shearforce = np.multiply(shearforce,10000)
    diameter = np.multiply(adaption_grp['edges/radius'],2)
    node_a_index =adaption_grp['edges/node_a_index']
    node_b_index =adaption_grp['edges/node_b_index']
        
    pressure_at_nodes = adaption_grp['nodes/pressure']
    
    
    pressure_at_vessel=[]
        
    for NodeA, NodeB in itertools.izip(node_a_index,node_b_index):
        pressure_at_vessel.append((pressure_at_nodes[NodeA]+pressure_at_nodes[NodeB])/2 * 7.5)
    pressure_at_vessel = np.array(pressure_at_vessel)
    fig = plt.figure()
    plt.subplot(2,1,1)
    a = pressure_at_vessel[index_of_artery]
    b = shearforce[index_of_artery]
    plt.plot(a,b,'or')
    a = pressure_at_vessel[index_of_capillary]
    b = shearforce[index_of_capillary]
    plt.plot(a,b,'yD')
    a = pressure_at_vessel[index_of_vein]
    b = shearforce[index_of_vein]
    plt.plot(a,b,'bv')
#    plt.semilogy(pressure_at_vessel[index_of_artery],shearforce[index_of_artery],'ro')
#    plt.semilogy(pressure_at_vessel[index_of_capillary],shearforce[index_of_capillary],'yD')
#    plt.semilogy(pressure_at_vessel[index_of_vein],shearforce[index_of_vein],'bv')
    plt.grid()
    plt.yscale('log')
    plt.xlabel('pressure/ mmHg')
    #plt.xlim([0,80])
    #plt.ylim([1,100])
    plt.ylabel('shearstress dyne/cm^2')
    
    plt.subplot(2,1,2)
    plt.plot(pressure_at_vessel[index_of_artery],diameter[index_of_artery],'ro')
    plt.plot(pressure_at_vessel[index_of_capillary],diameter[index_of_capillary],'yD')
    plt.plot(pressure_at_vessel[index_of_vein],diameter[index_of_vein],'bv')
    plt.grid()
    plt.legend(['ART','CAP','VEN'],loc='upper center')
    plt.xlabel('pressure/ mmHg')
    plt.ylabel('diameter/ mum')
    #plt.xlim([1,80])
    #plt.ylim([5,30])
    
    #pp.savefig()
    plt.savefig('All_signals_CONVERGENT_symetricA.png')

    print("hello")

def plot_symetricIrregular(fn,pp):
    adaption_grp = fn['symetricIrregular/vessels_after_adaption/vessels_after_adaption']
    index_of_artery = np.bitwise_and(np.asarray(adaption_grp['edges/flags']), ku.ARTERY)>0
    index_of_capillary = np.bitwise_and(np.asarray(adaption_grp['edges/flags']), ku.CAPILLARY)>0
    index_of_vein = np.bitwise_and(np.asarray(adaption_grp['edges/flags']), ku.VEIN)>0
    
    shearforce = adaption_grp['edges/shearforce']
    shearforce = np.multiply(shearforce,10000)
    diameter = np.multiply(adaption_grp['edges/radius'],2)
    node_a_index =adaption_grp['edges/node_a_index']
    node_b_index =adaption_grp['edges/node_b_index']
        
    pressure_at_nodes = adaption_grp['nodes/pressure']
    
    
    pressure_at_vessel=[]
        
    for NodeA, NodeB in itertools.izip(node_a_index,node_b_index):
        pressure_at_vessel.append((pressure_at_nodes[NodeA]+pressure_at_nodes[NodeB])/2 * 7.5)
    pressure_at_vessel = np.array(pressure_at_vessel)
    fig = plt.figure()
    plt.subplot(2,1,1)
    a = pressure_at_vessel[index_of_artery]
    b = shearforce[index_of_artery]
    plt.plot(a,b,'or')
    a = pressure_at_vessel[index_of_capillary]
    b = shearforce[index_of_capillary]
    plt.plot(a,b,'yD')
    a = pressure_at_vessel[index_of_vein]
    b = shearforce[index_of_vein]
    plt.plot(a,b,'bv')
#    plt.semilogy(pressure_at_vessel[index_of_artery],shearforce[index_of_artery],'ro')
#    plt.semilogy(pressure_at_vessel[index_of_capillary],shearforce[index_of_capillary],'yD')
#    plt.semilogy(pressure_at_vessel[index_of_vein],shearforce[index_of_vein],'bv')
    plt.grid()
    plt.yscale('log')
    plt.xlabel('pressure/ mmHg')
    plt.xlim([0,80])
    #plt.ylim([1,100])
    plt.ylabel('shearstress dyne/cm^2')
    
    plt.subplot(2,1,2)
    plt.plot(pressure_at_vessel[index_of_artery],diameter[index_of_artery],'ro')
    plt.plot(pressure_at_vessel[index_of_capillary],diameter[index_of_capillary],'yD')
    plt.plot(pressure_at_vessel[index_of_vein],diameter[index_of_vein],'bv')
    plt.grid()
    plt.legend(['ART','CAP','VEN'])
    plt.xlabel('pressure/ mmHg')
    plt.ylabel('diameter/ mum')
    plt.xlim([1,80])
    plt.ylim([5,30])
    
    pp.savefig()
    
    plt.savefig('All_signals_CONVERGENT_irregular')
    


if __name__ == '__main__':


  fn = 'test_configs.h5'
  f = h5py.File(fn)

  
#look at stuff
  common_filename = os.path.splitext(fn)[0]
  pp = PdfPages(common_filename + '_results_test_configs.pdf')
  
#  plot_symetricA(f)
#  plot_symetricB(f,pp)
#  plot_symetricIrregular(f,pp)
  plot_asymetric(f,pp)