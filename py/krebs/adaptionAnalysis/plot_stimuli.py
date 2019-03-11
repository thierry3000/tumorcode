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

if __name__ == '__main__':
  import os.path, sys
  sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../..'))

import matplotlib
import matplotlib.pyplot as plt
import identifycluster
if (identifycluster.getname()=='snowden' or identifycluster.getname()=='durga'):
  matplotlib.use('agg')
import h5py
import itertools
import numpy as np
import mpl_utils
import os
import krebsutils as ku


def plot_shearstress(vessel_grp,pp):
    index_of_artery = np.bitwise_and(np.asarray(vessel_grp['edges/flags']), ku.ARTERY)>0
    index_of_artery = index_of_artery[:,0]
    index_of_capillary = np.bitwise_and(np.asarray(vessel_grp['edges/flags']), ku.CAPILLARY)>0
    index_of_capillary = index_of_capillary[:,0]
    index_of_vein = np.bitwise_and(np.asarray(vessel_grp['edges/flags']), ku.VEIN)>0
    index_of_vein = index_of_vein[:,0]
  
    shearforce = np.asarray(vessel_grp['edges/shearforce'])
    shearforce = shearforce[:,0]
    shearforce = np.multiply(shearforce,10000)
    diameter = np.multiply(np.asarray(vessel_grp['edges/radius']),2)
    diameter = diameter[:,0]
    
    node_a_index = np.asarray(vessel_grp['edges/node_a_index'])
    node_a_index = node_a_index[:,0]
    node_b_index = np.asarray(vessel_grp['edges/node_b_index'])
    node_b_index = node_b_index[:,0]
        
    pressure_at_nodes = np.asarray(vessel_grp['nodes/pressure'])
    pressure_at_nodes = pressure_at_nodes[:,0]
    
    
    pressure_at_vessel=[]
        
    for NodeA, NodeB in itertools.izip(node_a_index,node_b_index):
        pressure_at_vessel.append((pressure_at_nodes[NodeA]+pressure_at_nodes[NodeB])/2 * 1/0.1333)
    
    fig = matplotlib.figure.Figure()
    ax1 = fig.add_subplot(111)
    ax1.tick_params(
    axis='both',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    labelsize = tickfontsize)
    
    #ax1.loglog(pressure_at_vessel,shearforce,'*')
    pressure_at_vessel=np.asarray(pressure_at_vessel)
    ax1.loglog(pressure_at_vessel[index_of_artery],shearforce[index_of_artery],'o', markerfacecolor='r', markeredgecolor='r')
    ax1.loglog(pressure_at_vessel[index_of_capillary],shearforce[index_of_capillary],'D',markerfacecolor='none', markeredgecolor='coral' )
    ax1.loglog(pressure_at_vessel[index_of_vein],shearforce[index_of_vein],'v', markerfacecolor='none', markeredgecolor='b')
    
    
    
    ax1.set_ylabel(r'shearstress dyne/cm$^2$', fontsize=18)
    ax1.yaxis.labelpad = -10
    if goodArguments.apj:
      ax1.set_xlim([12,100])
      ax1.set_ylim([1,1000])
      ax1.yaxis.set_major_locator(plt.NullLocator())
      ax1.yaxis.set_minor_locator(plt.NullLocator())
      ax1.set_yticks([1,10,100,1000])
      ax1.set_yticklabels(['1','10','100','1000'],fontsize=tickfontsize)
      ax1.xaxis.set_major_locator(plt.NullLocator())
      ax1.xaxis.set_minor_locator(plt.NullLocator())
      ax1.set_xticks([10,30,100])
      ax1.set_xticklabels(['10','30','100'],fontsize=tickfontsize)
    else:
      ax1.set_xlim([20,52])
      ax1.set_ylim([1,250])
      ax1.yaxis.set_major_locator(plt.NullLocator())
      ax1.yaxis.set_minor_locator(plt.NullLocator())
      ax1.set_yticks([1,50,250])
      ax1.set_yticklabels(['1','50','250'],fontsize=tickfontsize)
      ax1.xaxis.set_major_locator(plt.NullLocator())
      ax1.xaxis.set_minor_locator(plt.NullLocator())
      ax1.set_xticks([20,25,50])
      ax1.set_xticklabels(['20','25','50'],fontsize=tickfontsize)
    
    
    ax1.grid()
    
    pp.savefig(fig, 'hydorodynamic_charicteristics')
    
def plot_diameter(vessel_grp,pp):
    index_of_artery = np.bitwise_and(np.asarray(vessel_grp['edges/flags']), ku.ARTERY)>0
    index_of_artery = index_of_artery[:,0]
    index_of_capillary = np.bitwise_and(np.asarray(vessel_grp['edges/flags']), ku.CAPILLARY)>0
    index_of_capillary = index_of_capillary[:,0]
    index_of_vein = np.bitwise_and(np.asarray(vessel_grp['edges/flags']), ku.VEIN)>0
    index_of_vein = index_of_vein[:,0]
  
    shearforce = np.asarray(vessel_grp['edges/shearforce'])
    shearforce = shearforce[:,0]
    shearforce = np.multiply(shearforce,10000)
    diameter = np.multiply(np.asarray(vessel_grp['edges/radius']),2)
    diameter = diameter[:,0]
    
    node_a_index = np.asarray(vessel_grp['edges/node_a_index'])
    node_a_index = node_a_index[:,0]
    node_b_index = np.asarray(vessel_grp['edges/node_b_index'])
    node_b_index = node_b_index[:,0]
        
    pressure_at_nodes = np.asarray(vessel_grp['nodes/pressure'])
    pressure_at_nodes = pressure_at_nodes[:,0]
    
    
    pressure_at_vessel=[]
        
    for NodeA, NodeB in itertools.izip(node_a_index,node_b_index):
        pressure_at_vessel.append((pressure_at_nodes[NodeA]+pressure_at_nodes[NodeB])/2 * 1/0.1333)
    
    fig = matplotlib.figure.Figure()
    ax2 = fig.add_subplot(111)
    
    #ax1.loglog(pressure_at_vessel,shearforce,'*')
    pressure_at_vessel=np.asarray(pressure_at_vessel)
    
    
    #ax2.loglog(pressure_at_vessel,diameter,'r*')
    if goodArguments.apj:
      ax2.loglog(pressure_at_vessel[index_of_artery],diameter[index_of_artery],'o', markerfacecolor='r', markeredgecolor='r', label='ART')
      ax2.loglog(pressure_at_vessel[index_of_capillary],diameter[index_of_capillary],'D',markerfacecolor='none', markeredgecolor='coral', label='CAP' )
      ax2.loglog(pressure_at_vessel[index_of_vein],diameter[index_of_vein],'v', markerfacecolor='none', markeredgecolor='b', label='VEN')
    else:
      ax2.loglog(pressure_at_vessel[index_of_artery],diameter[index_of_artery],'o', markerfacecolor='r', markeredgecolor='r', label='ART')
      ax2.loglog(pressure_at_vessel[index_of_capillary],diameter[index_of_capillary],'D',markerfacecolor='none', markeredgecolor='coral', label='CAP' )
      ax2.loglog(pressure_at_vessel[index_of_vein],diameter[index_of_vein],'v', markerfacecolor='none', markeredgecolor='b', label='VEN')
      ax2.set_ylim([5,140])
      ax2.set_xlim([20,52])
    ax2.set_ylabel(r'diameter/ $\mu m$', fontsize=18)
    ax2.set_xlabel(r'pressure/ $mmHg$', fontsize=18)
    if goodArguments.apj:
      ax2.set_ylim([5,100])
      ax2.set_xlim([10,100])
      ax2.xaxis.set_major_locator(plt.NullLocator())
      ax2.xaxis.set_minor_locator(plt.NullLocator())
      ax2.set_xticks([10,30,100])
      ax2.set_xticklabels(['10','30','100'],fontsize=tickfontsize)
      ax2.yaxis.set_major_locator(plt.NullLocator())
      ax2.yaxis.set_minor_locator(plt.NullLocator())
      ax2.set_yticks([10,100])
      ax2.set_yticklabels(['10','100'],fontsize=tickfontsize)
    else:
      ax2.set_ylim([4.5,120])
      ax2.set_xlim([20,55])
      ax2.xaxis.set_major_locator(plt.NullLocator())
      ax2.xaxis.set_minor_locator(plt.NullLocator())
      ax2.set_xticks([20,25,50])
      ax2.set_xticklabels(['20','25','50'],fontsize=tickfontsize)
      ax2.yaxis.set_major_locator(plt.NullLocator())
      ax2.yaxis.set_minor_locator(plt.NullLocator())
      ax2.set_yticks([5,10,100])
      ax2.set_yticklabels(['5','10','100'],fontsize=tickfontsize)
    
    ax2.grid()
    #labels = ['ART', 'CAP', 'VEN']
    #dummies = [ax2.plot([], [], ls='o', c=c)[0] for c in colors]
    ax2.legend(loc='upper center')
    pp.savefig(fig, 'hydorodynamic_charicteristics')

def plot_hydrodynamic_stimuli(vessel_grp,pp):
    index_of_artery = np.bitwise_and(np.asarray(vessel_grp['edges/flags']), ku.ARTERY)>0
    index_of_artery = index_of_artery[:,0]
    index_of_capillary = np.bitwise_and(np.asarray(vessel_grp['edges/flags']), ku.CAPILLARY)>0
    index_of_capillary = index_of_capillary[:,0]
    index_of_vein = np.bitwise_and(np.asarray(vessel_grp['edges/flags']), ku.VEIN)>0
    index_of_vein = index_of_vein[:,0]
    
    shearforce = np.asarray(vessel_grp['edges/shearforce'])
    shearforce = shearforce[:,0]
    shearforce = np.multiply(shearforce,10000) #kpa to dyne
    node_a_index =np.asarray(vessel_grp['edges/node_a_index'])
    node_a_index = node_a_index[:,0]
    node_b_index =np.asarray(vessel_grp['edges/node_b_index'])
    node_b_index= node_b_index[:,0]
        
    pressure_at_nodes = np.asarray(vessel_grp['nodes/pressure'])
    pressure_at_nodes = pressure_at_nodes[:,0]
    
    
    pressure_at_vessel=[]
        
    for NodeA, NodeB in itertools.izip(node_a_index,node_b_index):
        pressure_at_vessel.append((pressure_at_nodes[NodeA]+pressure_at_nodes[NodeB])/2 * 1/0.1333)
    pressure_at_vessel = np.asarray(pressure_at_vessel)
    pressure_stimuli = 100 - 86 *np.exp(-5000*np.log10(np.log10(pressure_at_vessel))**5.4)
    fig = matplotlib.figure.Figure()
    ax = fig.add_subplot(111)
    ax.tick_params(
    axis='both',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    labelsize = tickfontsize)
    #ax.semilogx(pressure_at_vessel,np.log10(shearforce),'*')
    ax.semilogx(pressure_at_vessel[index_of_artery],np.log10(shearforce[index_of_artery]),'o', markerfacecolor='r', markeredgecolor='r', label='ART')
    ax.semilogx(pressure_at_vessel[index_of_capillary],np.log10(shearforce[index_of_capillary]),'D',markerfacecolor='none', markeredgecolor='coral', label='CAP' )
    ax.semilogx(pressure_at_vessel[index_of_vein],np.log10(shearforce[index_of_vein]),'v', markerfacecolor='none', markeredgecolor='b', label='VEN')
    
    #ax.semilogx(pressure_at_vessel,-np.log10(pressure_stimuli),'r*')
    ax.semilogx(pressure_at_vessel[index_of_artery],-np.log10(pressure_stimuli[index_of_artery]),'o', markerfacecolor='r', markeredgecolor='r')
    ax.semilogx(pressure_at_vessel[index_of_capillary],-np.log10(pressure_stimuli[index_of_capillary]),'D',markerfacecolor='none', markeredgecolor='coral' )
    ax.semilogx(pressure_at_vessel[index_of_vein],-np.log10(pressure_stimuli[index_of_vein]),'v', markerfacecolor='none', markeredgecolor='b')
    #ax.legend([r'$log_{10}($ shearforce $)$', r'$-log_{10}($ expected wallstress $)$'])
    ax.legend(loc='center right')
    ax.grid()
    ax.set_xlabel('pressure/ mmHg', fontsize=18)
    if goodArguments.apj:
      print("goodArguments.FileNames: %s" % goodArguments.FileNames[0].name)
      ax.set_ylim([-2.2,4])
      if not 'initial' in goodArguments.FileNames[0].name:
        print('first')
        ax.xaxis.set_major_locator(plt.NullLocator())
        ax.xaxis.set_minor_locator(plt.NullLocator())
        ax.set_xlim([12,120])
        ax.set_xticks([12,30,120])
        ax.set_xticklabels(['12','30','120'],fontsize=tickfontsize)
        
      else:
        print('else')
        ax.xaxis.set_major_locator(plt.NullLocator())
        ax.xaxis.set_minor_locator(plt.NullLocator())
        ax.set_xlim([10,100])
        ax.set_xticks([10,30,100])
        ax.set_xticklabels(['10','30','100'],fontsize=tickfontsize)
    else:
      ax.set_xlim([15,55])
      ax.set_ylim([-2.2,2.5])
      ax.xaxis.set_major_locator(plt.NullLocator())
      ax.xaxis.set_minor_locator(plt.NullLocator())
      ax.set_xticks([15,30,50])
      ax.set_xticklabels(['15','30','50'],fontsize=tickfontsize)
      #ax2.set_ylim([3,100])
    #ax.xaxis.set_label_coords(1.05, -20.025)
    #ax.set_xlim([10,100])
    #ax.set_ylim([-2.2,4])
    #ax.set_ylabel('stimuli')
    
    pp.savefig(fig, 'Hydrodynamic stimuli')
    #plt.show() 
    
def plot_conductive_stimuli(adaption_grp,pp):
    index_of_artery = np.bitwise_and(np.asarray(adaption_grp['edges/flags']), ku.ARTERY)>0
    index_of_artery = index_of_artery[:,0]
    index_of_capillary = np.bitwise_and(np.asarray(adaption_grp['edges/flags']), ku.CAPILLARY)>0
    index_of_capillary = index_of_capillary[:,0]
    index_of_vein = np.bitwise_and(np.asarray(adaption_grp['edges/flags']), ku.VEIN)>0
    index_of_vein = index_of_vein[:,0]
    metabolic = np.asarray(adaption_grp['edges/metabolicSignal'])
    metabolic = metabolic[:,0]
    conductive = np.asarray(adaption_grp['edges/conductivitySignal'])
    conductive = conductive[:,0]
    flow = np.asarray(adaption_grp['edges/flow'])
    flow =flow[:,0]
    flow = np.multiply(flow,60./1000000.)
    
    node_a_index = np.asarray(adaption_grp['edges/node_a_index'])
    node_a_index = node_a_index[:,0]
    node_b_index = np.asarray(adaption_grp['edges/node_b_index'])
    node_b_index = node_b_index[:,0]
    
    pressure_at_nodes = np.asarray(adaption_grp['nodes/pressure'])
    pressure_at_nodes = pressure_at_nodes[:,0]
    
    pressure_at_vessel=[]
    
    for NodeA, NodeB in itertools.izip(node_a_index,node_b_index):
        pressure_at_vessel.append((pressure_at_nodes[NodeA]+pressure_at_nodes[NodeB])/2 * 1/0.1333)
    
    fig2 = matplotlib.figure.Figure()
    ax = fig2.add_subplot(111)
    #ax.set_title('Stimuli')
    ax.semilogx(flow[index_of_artery],metabolic[index_of_artery],'o', markerfacecolor ='r', markeredgecolor='r')
    ax.semilogx(flow[index_of_capillary],metabolic[index_of_capillary], 'o', markerfacecolor ='coral', markeredgecolor='coral')
    ax.semilogx(flow[index_of_vein],metabolic[index_of_vein],'v', markerfacecolor ='b', markeredgecolor='b' )
    
    ax.semilogx(flow[index_of_artery],conductive[index_of_artery],'o', markerfacecolor ='none', markeredgecolor='r')
    ax.semilogx(flow[index_of_capillary],conductive[index_of_capillary],'o', markerfacecolor ='none', markeredgecolor='coral')
    ax.semilogx(flow[index_of_vein],conductive[index_of_vein],'v', markerfacecolor ='none', markeredgecolor='b')
    datasets = ax.get_lines()
    legend1=ax.legend([datasets[i] for i in [0,1,2]],['metabolic/ artery','metabolic/ capillary', 'metabolic/ vein'], loc=1)
    legend2=ax.legend([datasets[i] for i in [3,4,5]],['conductive/ artery','conductive/ capillary', 'conductive/ vein'], loc=2)
    #legend1=ax.legend(['metabolic/ artery','metabolic/ capillary', 'metabolic/ vein'])
    #legend2=ax.legend(['bld', 'blub'])
    ax.add_artist(legend1)
    ax.grid()
    ax.set_xlabel("flow/nl/min", fontsize=18)
    ax.set_ylabel("stimuli", fontsize=18)
    if goodArguments.apj:
      ax.set_xlim([0.01,1000])
      ax.set_ylim([0,4])
      ax.yaxis.set_major_locator(plt.NullLocator())
      ax.yaxis.set_minor_locator(plt.NullLocator())
      ax.set_yticks([0,1,2,3,4])
      ax.set_yticklabels(['0','1','2','3','4'],fontsize=tickfontsize)
      ax.xaxis.set_major_locator(plt.NullLocator())
      ax.xaxis.set_minor_locator(plt.NullLocator())
      ax.set_xticks([0.01,0.1,1,10,100,1000])
      ax.set_xticklabels(['0.01','0.1','1','10','100','1000'],fontsize=tickfontsize)
    else:
      ax.set_xlim([0.01,1500])
      ax.set_ylim([0,3.2])
      ax.yaxis.set_major_locator(plt.NullLocator())
      ax.yaxis.set_minor_locator(plt.NullLocator())
      ax.set_yticks([0,1,2,3])
      ax.set_yticklabels(['0','1','2','3'],fontsize=tickfontsize)
      ax.xaxis.set_major_locator(plt.NullLocator())
      ax.xaxis.set_minor_locator(plt.NullLocator())
      ax.set_xticks([0.01,0.1,1,10,100,1000])
      ax.set_xticklabels(['0.01','0.1','1','10','100','1000'],fontsize=tickfontsize)
    pp.savefig(fig2,'conductive stimuli')
    #plt.show()
def plot_movie(f,pp):
  #all_grp = f['Asymetric/vessels_after_adaption']
  all_grp = f['adaption']
    
  max_iterlength=len(all_grp.keys())-5
  i=1
  k=0
  while i<max_iterlength :
    #adaption_grp= f['/Asymetric/vessels_after_adaption/vessels_after_adaption_%i'%i]
    adaption_grp= f['/adaption/vessels_after_adaption_%i'%i]
    index_of_artery = np.bitwise_and(np.asarray(adaption_grp['edges/flags']), ku.ARTERY)>0
    index_of_capillary = np.bitwise_and(np.asarray(adaption_grp['edges/flags']), ku.CAPILLARY)>0
    index_of_vein = np.bitwise_and(np.asarray(adaption_grp['edges/flags']), ku.VEIN)>0
    
    shearforce = adaption_grp['edges/shearforce']
    shearforce = np.multiply(shearforce,10000)
    shearforce = np.log10(shearforce)
    diameter = np.multiply(adaption_grp['edges/radius'],2)
    node_a_index =adaption_grp['edges/node_a_index']
    node_b_index =adaption_grp['edges/node_b_index']
        
    pressure_at_nodes = adaption_grp['nodes/pressure']
    
    
    pressure_at_vessel=[]
        
    for NodeA, NodeB in itertools.izip(node_a_index,node_b_index):
        pressure_at_vessel.append((pressure_at_nodes[NodeA]+pressure_at_nodes[NodeB])/2 * 7.5)
    pressure_at_vessel = np.array(pressure_at_vessel)
    fig = plt.figure()
    plt.subplot(5,1,1)
    a = pressure_at_vessel[index_of_artery]
    sig1a = shearforce[index_of_artery]
    plt.plot(a,sig1a,'or')
    a = pressure_at_vessel[index_of_capillary]
    sig1c = shearforce[index_of_capillary]
    plt.plot(a,sig1c,'yD')
    a = pressure_at_vessel[index_of_vein]
    sig1v = shearforce[index_of_vein]
    plt.plot(a,sig1v,'bv')
#    plt.semilogy(pressure_at_vessel[index_of_artery],shearforce[index_of_artery],'ro')
#    plt.semilogy(pressure_at_vessel[index_of_capillary],shearforce[index_of_capillary],'yD')
#    plt.semilogy(pressure_at_vessel[index_of_vein],shearforce[index_of_vein],'bv')
    plt.grid()
    plt.xlabel('pressure/ mmHg')
    plt.xlim([10,80])
    plt.ylim([1.2,3])
    plt.ylabel('sig1')
    
    plt.subplot(5,1,2)
    signal2 = 100-86*np.exp(-5000*np.log10(np.log10(pressure_at_vessel))**5.4)
    signal2 = -np.log10(signal2)
    pa = pressure_at_vessel[index_of_artery]
    sig2a = signal2[index_of_artery]
    plt.plot(pa,sig2a,'or')
    pc = pressure_at_vessel[index_of_capillary]
    sig2c = signal2[index_of_capillary]
    plt.plot(pc,sig2c,'yD')
    pv = pressure_at_vessel[index_of_vein]
    sig2v = signal2[index_of_vein]
    plt.plot(pv,sig2v,'bv')
    #plt.subplot(4,1,3)
#    plt.plot(signal2[index_of_artery],diameter[index_of_artery],'ro')
#    plt.plot(signal2[index_of_capillary],diameter[index_of_capillary],'yD')
#    plt.plot(signal2[index_of_vein],diameter[index_of_vein],'bv')
    plt.grid()
    #plt.legend(['ART','CAP','VEN'])
    plt.xlabel('pressure/ mmHg')
    plt.xlim([10,80])
    plt.ylim([-2.,-1.1])
    plt.ylabel('sig2')
    #### signal 3
    plt.subplot(5,1,3)
    metabolic = adaption_grp['edges/metabolicSignal']
    conductive = adaption_grp['edges/conductivitySignal']
    flow = adaption_grp['edges/flow']
    flow = np.multiply(flow,60./1000000.)
#    plt.semilogx(flow[index_of_artery],metabolic[index_of_artery],'or')
#    plt.semilogx(flow[index_of_capillary],metabolic[index_of_capillary],'yD')
#    plt.semilogx(flow[index_of_vein],metabolic[index_of_vein],'bv')
    plt.plot(pa,metabolic[index_of_artery],'or')
    plt.plot(pc,metabolic[index_of_capillary],'yD')
    plt.plot(pv,metabolic[index_of_vein],'bv')
    
    plt.ylabel('sig3')
    plt.xlim([10,80])
    plt.ylim([0,1.2])
    plt.grid()
    
    #### signal 4
    plt.subplot(5,1,4)
#    plt.semilogx(flow[index_of_artery],conductive[index_of_artery],'or')
#    plt.semilogx(flow[index_of_capillary],conductive[index_of_capillary],'yD')
#    plt.semilogx(flow[index_of_vein],conductive[index_of_vein],'bv')
    plt.plot(pa,conductive[index_of_artery],'or')
    plt.plot(pc,conductive[index_of_capillary],'yD')
    plt.plot(pv,conductive[index_of_vein],'bv')
     
    plt.ylabel('sig4')
    plt.xlabel('pressure/mmHg')
    plt.xlim([10,80])
    plt.ylim([0,1.4])
    plt.grid()
    ### sum
    plt.subplot(5,1,5)
    plt.plot(pa,sig1a+sig2a+metabolic[index_of_artery]+conductive[index_of_artery],'or')
    plt.plot(pc,sig1c+sig2c+metabolic[index_of_capillary]+conductive[index_of_capillary],'yD')
    plt.plot(pv,sig1v+sig2v+metabolic[index_of_vein]+conductive[index_of_vein],'bv')
    
    plt.ylabel("sum")
    plt.xlim([10,80])
    plt.ylim([0,3.])
    plt.grid()
    plt.savefig("mov_%03i.png"%k)
    plt.close()
    #os.system("python2 /localdisk/thierry/tumorcode/py/krebs/povrayRenderVessels.py ../test_configs.h5 /Asymetric/vessels_after_adaption/vessels_after_adaption_%03i"%i)
    os.system("python2 /daten/tumorcode/py/krebs/povrayRenderVessels.py /daten/localdisk/adaption_project/vessels-q2d-8mm-P6-typeE-9x3L130-sample00_adption_p_typeE.h5 /adaption/vessels_after_adaption_%i"%i)
    i=i+100
    k=k+1
    #plt.show()
def plot_movie_typeE(f,pp):
  #all_grp = f['Asymetric/vessels_after_adaption']
  all_grp = f['adaption']
  x_min=20
  x_max=50
  max_iterlength=len(all_grp.keys())-5
  i=1
  k=0
  while i<max_iterlength :
    #adaption_grp= f['/Asymetric/vessels_after_adaption/vessels_after_adaption_%i'%i]
    adaption_grp= f['/adaption/vessels_after_adaption_%i'%i]
    index_of_artery = np.bitwise_and(np.asarray(adaption_grp['edges/flags']), ku.ARTERY)>0
    index_of_capillary = np.bitwise_and(np.asarray(adaption_grp['edges/flags']), ku.CAPILLARY)>0
    index_of_vein = np.bitwise_and(np.asarray(adaption_grp['edges/flags']), ku.VEIN)>0
    
    shearforce = adaption_grp['edges/shearforce']
    shearforce = np.multiply(shearforce,10000)
    shearforce = np.log10(shearforce)
    diameter = np.multiply(adaption_grp['edges/radius'],2)
    node_a_index =adaption_grp['edges/node_a_index']
    node_b_index =adaption_grp['edges/node_b_index']
        
    pressure_at_nodes = adaption_grp['nodes/pressure']
    
    
    pressure_at_vessel=[]
        
    for NodeA, NodeB in itertools.izip(node_a_index,node_b_index):
        pressure_at_vessel.append((pressure_at_nodes[NodeA]+pressure_at_nodes[NodeB])/2 * 7.5)
    pressure_at_vessel = np.array(pressure_at_vessel)
    fig = plt.figure()
    plt.subplot(5,1,1)
    a = pressure_at_vessel[index_of_artery]
    sig1a = shearforce[index_of_artery]
    plt.plot(a,sig1a,'or')
    a = pressure_at_vessel[index_of_capillary]
    sig1c = shearforce[index_of_capillary]
    plt.plot(a,sig1c,'yD')
    a = pressure_at_vessel[index_of_vein]
    sig1v = shearforce[index_of_vein]
    plt.plot(a,sig1v,'bv')
#    plt.semilogy(pressure_at_vessel[index_of_artery],shearforce[index_of_artery],'ro')
#    plt.semilogy(pressure_at_vessel[index_of_capillary],shearforce[index_of_capillary],'yD')
#    plt.semilogy(pressure_at_vessel[index_of_vein],shearforce[index_of_vein],'bv')
    plt.grid()
    plt.xlabel('pressure/ mmHg')
    plt.xlim([x_min,x_max])
    plt.ylim([-2,2.5])
    plt.ylabel('sig1')
    
    plt.subplot(5,1,2)
    signal2 = 100-86*np.exp(-5000*np.log10(np.log10(pressure_at_vessel))**5.4)
    signal2 = -np.log10(signal2)
    pa = pressure_at_vessel[index_of_artery]
    sig2a = signal2[index_of_artery]
    plt.plot(pa,sig2a,'or')
    pc = pressure_at_vessel[index_of_capillary]
    sig2c = signal2[index_of_capillary]
    plt.plot(pc,sig2c,'yD')
    pv = pressure_at_vessel[index_of_vein]
    sig2v = signal2[index_of_vein]
    plt.plot(pv,sig2v,'bv')
    #plt.subplot(4,1,3)
#    plt.plot(signal2[index_of_artery],diameter[index_of_artery],'ro')
#    plt.plot(signal2[index_of_capillary],diameter[index_of_capillary],'yD')
#    plt.plot(signal2[index_of_vein],diameter[index_of_vein],'bv')
    plt.grid()
    #plt.legend(['ART','CAP','VEN'])
    plt.xlabel('pressure/ mmHg')
    plt.xlim([x_min,x_max])
    plt.ylim([-1.9,-1.3])
    plt.ylabel('sig2')
    #### signal 3
    plt.subplot(5,1,3)
    metabolic = adaption_grp['edges/metabolicSignal']
    conductive = adaption_grp['edges/conductivitySignal']
    flow = adaption_grp['edges/flow']
    flow = np.multiply(flow,60./1000000.)
#    plt.semilogx(flow[index_of_artery],metabolic[index_of_artery],'or')
#    plt.semilogx(flow[index_of_capillary],metabolic[index_of_capillary],'yD')
#    plt.semilogx(flow[index_of_vein],metabolic[index_of_vein],'bv')
    plt.plot(pa,metabolic[index_of_artery],'or')
    plt.plot(pc,metabolic[index_of_capillary],'yD')
    plt.plot(pv,metabolic[index_of_vein],'bv')
    
    plt.ylabel('sig3')
    plt.xlim([x_min,x_max])
    plt.ylim([0,3.5])
    plt.grid()
    
    #### signal 4
    plt.subplot(5,1,4)
#    plt.semilogx(flow[index_of_artery],conductive[index_of_artery],'or')
#    plt.semilogx(flow[index_of_capillary],conductive[index_of_capillary],'yD')
#    plt.semilogx(flow[index_of_vein],conductive[index_of_vein],'bv')
    plt.plot(pa,conductive[index_of_artery],'or')
    plt.plot(pc,conductive[index_of_capillary],'yD')
    plt.plot(pv,conductive[index_of_vein],'bv')
     
    plt.ylabel('sig4')
    plt.xlabel('pressure/mmHg')
    plt.xlim([x_min,x_max])
    plt.ylim([0,3])
    plt.grid()
    ### sum
    plt.subplot(5,1,5)
    plt.plot(pa,sig1a+sig2a+metabolic[index_of_artery]+conductive[index_of_artery],'or')
    plt.plot(pc,sig1c+sig2c+metabolic[index_of_capillary]+conductive[index_of_capillary],'yD')
    plt.plot(pv,sig1v+sig2v+metabolic[index_of_vein]+conductive[index_of_vein],'bv')
    
    plt.ylabel("sum")
    plt.xlim([x_min,x_max])
    plt.ylim([-3,5.])
    plt.grid()
    plt.savefig("mov_%03i.png"%k)
    plt.close()
    #os.system("python2 /localdisk/thierry/tumorcode/py/krebs/povrayRenderVessels.py ../test_configs.h5 /Asymetric/vessels_after_adaption/vessels_after_adaption_%03i"%i)
    os.system("python2 /daten/tumorcode/py/krebs/povrayRenderVessels.py /daten/localdisk/adaption_project/vessels-q2d-8mm-P6-typeE-9x3L130-sample00_adption_p_typeE.h5 /adaption/vessels_after_adaption_%i"%i)
    i=i+100
    k=k+1
    #plt.show()
    
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Plot adaption stimuli')  
    parser.add_argument('FileNames', nargs='*', type=argparse.FileType('r'), default=sys.stdin, help='file containing a proper adaption')   
    parser.add_argument('--apj', help='sets the axis acording to Secombs apj paper', default=False,  action='store_true')
    goodArguments, otherArguments = parser.parse_known_args()    
    #create filename due to former standards
    filenames=[]
    for fn in goodArguments.FileNames:
      filenames.append(fn.name)   
    
    rc = matplotlib.rc
    rc('font', size = 8.)
    rc('axes', titlesize = 10., labelsize = 8.)
    tickfontsize=14
    for fn in filenames:
      with h5py.File(fn) as f:
        no_of_iterations = ''
        no_of_iterations = 'plot'
        vesselgrp = f['vessels_after_adaption']
        common_filename = os.path.splitext(os.path.basename(fn))[0]
      
        with mpl_utils.PdfWriter(no_of_iterations + '_' + common_filename + '_hydrodynamic_stimulies.pdf') as pp:
          hydrodynamic_fig = plot_hydrodynamic_stimuli(vesselgrp, pp)
        with mpl_utils.PdfWriter(no_of_iterations + '_' + common_filename + '_conductive_stimulies.pdf') as pp:
          plot_conductive_stimuli(vesselgrp,pp)
        with mpl_utils.PdfWriter(no_of_iterations + '_' + common_filename + '_diameter_hydrodynamics.pdf') as pp:
          plot_diameter(vesselgrp,pp)
        with mpl_utils.PdfWriter(no_of_iterations + '_' + common_filename + '_shearstress_hydrodynamics.pdf') as pp:
          plot_shearstress(vesselgrp,pp)
          
  
      #plot_movie(f,pp=None)
      #plot_movie_typeE(f,pp=None)
      
      

    
    
