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
"""
Created on Mon Sep 14 10:42:14 2015

@author: thierry
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 11:16:44 2015

@author: thierry
"""

import h5py
import numpy as np
import krebsutils as ku
import os

import myutils
from matplotlib.pylab import * 
import matplotlib.pyplot as plt

from analyzeGeneral   import DataBasicVessel
from analyzeBloodFlowResistance import ComputeVascularTreeBloodFlowResistances

def DoReadIn(filenames, pattern, fn_measure):
    #read in lots of stuff
    files = [h5files.open(fn, 'r') for fn in filenames]
    output_f = h5py.File(fn_measure)
    for afile in files:
        #edges
        index_of_circulated = np.bitwise_and(np.asarray(afile[pattern + '/edges/flags']), ku.CIRCULATED)>0
        index_of_arteries = np.bitwise_and(np.asarray(afile[pattern + '/edges/flags']), ku.ARTERY)>0
        index_of_veins = np.bitwise_and(np.asarray(afile[pattern + '/edges/flags']), ku.VEIN)>0
        index_of_capillaries = np.bitwise_and(np.asarray(afile[pattern + '/edges/flags']), ku.CAPILLARY)>0
        index_of_boundary = np.bitwise_and(np.asarray(afile[pattern + '/edges/flags']), ku.BOUNDARY)>0
        indeces_of_circulated_arteries = np.bitwise_and(index_of_circulated, index_of_arteries)
        indeces_of_circulated_veins = np.bitwise_and(index_of_circulated, index_of_veins)
        indeces_of_circulated_capillaries = np.bitwise_and(index_of_circulated, index_of_capillaries)
        
        radii_of_capillaries = np.asarray(afile[pattern + '/edges/radius'])[indeces_of_circulated_capillaries]
        radii_of_arteries = np.asarray(afile[pattern + '/edges/radius'])[indeces_of_circulated_arteries]
        radii_of_veins = np.asarray(afile[pattern + '/edges/radius'])[indeces_of_circulated_veins]
        radii_all = np.asarray(afile[pattern + '/edges/radius'])[index_of_circulated]
        flow_of_capillaries = np.asarray(afile[pattern + '/edges/flow'])[indeces_of_circulated_capillaries]
        flow_of_arteries = np.asarray(afile[pattern + '/edges/flow'])[indeces_of_circulated_arteries]
        flow_of_veins = np.asarray(afile[pattern + '/edges/flow'])[indeces_of_circulated_veins]
        flow_all = np.asarray(afile[pattern + '/edges/flow'])[index_of_circulated]
        #nodes
        index_of_circulated_nodes = np.bitwise_and(np.asarray(afile[pattern + '/nodes/nodeflags']), ku.CIRCULATED)>0
        indeces_of_roots = np.asarray(afile[pattern + '/nodes/roots'])
        artery_roots = []
        artery_roots_radii =[]
        artery_roots_pressure = []
        venous_roots = []
        venous_roots_radii = []
        venous_roots_pressure = []
        
        nodeflags = np.asarray(afile[pattern + '/nodes/nodeflags'])
        for aIndex in indeces_of_roots:
            if index_of_circulated_nodes[aIndex] and np.bitwise_and(nodeflags[aIndex], ku.ARTERY):
                artery_roots.append(aIndex)
                artery_roots_pressure.append(np.asarray(afile[pattern + '/nodes/pressure'])[aIndex])
            if index_of_circulated_nodes[aIndex] and np.bitwise_and(nodeflags[aIndex], ku.VEIN):
                venous_roots.append(aIndex)
                venous_roots_pressure.append(np.asarray(afile[pattern + '/nodes/pressure'])[aIndex])
        NodeAIndex = np.asarray(afile[pattern + '/edges/node_a_index'])[index_of_boundary]
        NodeBIndex = np.asarray(afile[pattern + '/edges/node_b_index'])[index_of_boundary]
        radii_of_boundary = np.asarray(afile[pattern + '/edges/radius'])[index_of_boundary]
        for (i,aNodeAIndex) in enumerate(NodeAIndex):
            if (aNodeAIndex in artery_roots) and not (aNodeAIndex in NodeBIndex):
                artery_roots_radii.append(radii_of_boundary[i])
        for (i,aNodeBIndex) in enumerate(NodeBIndex):
            if aNodeBIndex in artery_roots:
                artery_roots_radii.append(radii_of_boundary[i])
                
        for (i,aNodeAIndex) in enumerate(NodeAIndex):
            if aNodeAIndex in venous_roots and not (aNodeAIndex in NodeBIndex):
                venous_roots_radii.append(radii_of_boundary[i])
        for (i,aNodeBIndex) in enumerate(NodeBIndex):
            if aNodeBIndex in venous_roots:
                venous_roots_radii.append(radii_of_boundary[i])
              
        if 1:
            print("Circulated: %i " % np.sum(index_of_circulated))
            print("Circulated Arteries: %i " % np.sum(indeces_of_circulated_arteries))
            print("Circulated Veins: %i " % np.sum(indeces_of_circulated_veins))
            print("Circulated Capillaries: %i " % np.sum(indeces_of_circulated_capillaries))
        if 1:
            print("# Circulated Nodes: %i " % np.sum(index_of_circulated_nodes))
            print("# Arterious roots: %i " % len(artery_roots))
            print("# Venous roots: %i " % len(venous_roots))
            #print("Circulated Veins: %i " % np.sum(indeces_of_circulated_veins))
            #print("Circulated Capillaries: %i " % np.sum(indeces_of_circulated_capillaries))
        
        #write to file
        if 1:
            print("Working on file: %s"% afile.filename)
            output_grp = output_f.create_group(os.path.basename(afile.filename))
            #output_grp.attrs.create("#circulated")
            output_grp.attrs['#circulated'] = np.sum(index_of_circulated)
            N = np.sum(index_of_circulated)
            output_grp.attrs['#artery'] = np.sum(indeces_of_circulated_arteries)
            output_grp.attrs['#veins'] = np.sum(indeces_of_circulated_veins)
            output_grp.attrs['#capillary'] = np.sum(indeces_of_circulated_capillaries)
            l = afile[pattern + '/lattice'].attrs['SCALE']
            output_grp.attrs['length'] = l
            if(np.var(radii_of_capillaries))>0.0001:
                break
            r_0=np.mean(radii_of_capillaries)
            output_grp.attrs['c_radius'] = r_0
            output_grp.create_dataset('a_roots', data = artery_roots)
            output_grp.create_dataset('a_roots_radius', data = artery_roots_radii)
            output_grp.create_dataset('a_roots_pressure', data = artery_roots_pressure)
            output_grp.create_dataset('v_roots', data = venous_roots)
            output_grp.create_dataset('v_roots_radius', data = venous_roots)
            output_grp.create_dataset('v_roots_pressure', data = venous_roots_pressure)
            
            #since v2
            output_grp.create_dataset('radii_of_arteries', data = radii_of_arteries)
            output_grp.create_dataset('volumes_of_arteries', data = [ 2*np.pi*r**2*l for r in radii_of_arteries])#calculate tissue volume from radii
            output_grp.create_dataset('radii_of_veins', data = radii_of_veins)
            output_grp.create_dataset('volumes_of_veins', data = [ 2*np.pi*r**2*l for r in radii_of_veins])#calculate tissue volume from radii
            output_grp.create_dataset('radii_of_capillaries', data = radii_of_capillaries)
            output_grp.create_dataset('volumes_of_capillaries', data = [ 2*np.pi*r**2*l for r in radii_of_capillaries])#calculate tissue volume from radii
            output_grp.create_dataset('radii_all', data = radii_all)
            output_grp.create_dataset('volumes_all', data = [ 2*np.pi*r**2*l for r in radii_all])#calculate tissue volume from radii
            
            output_grp.create_dataset('flow_of_arteries', data = flow_of_arteries)
            output_grp.create_dataset('flow_of_veins', data = flow_of_veins)
            output_grp.create_dataset('flow_of_capillaries', data = flow_of_capillaries)
            output_grp.create_dataset('flow_all', data = flow_all)
            
            
            dataman = myutils.DataManager(100, [DataBasicVessel()])
            vessels = dataman.obtain_data('vessel_graph', afile['/' + pattern + '/'], ['flow', 'pressure', 'flags', 'radius','nodeflags'])
            conductivities, avgVenousPressure, avgArterialPressure, totalFlow = ComputeVascularTreeBloodFlowResistances(vessels)
            output_grp.attrs['avgVenousPressure'] = avgVenousPressure
            output_grp.attrs['avgArterialPressure'] = avgArterialPressure
            output_grp.attrs['totalFlow_from_BloodFlowResistance'] = totalFlow
            
def DoCalculate(fn_measure):
    f_measure = h5files.open(fn_measure,'r+')
    if 'means' in f_measure:
        del f_measure['/means']
        print('Delete!')
    #for graphics
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
    color_index = 0
    f_measure.create_group('means')
    
    if 0:
        fig1 = matplotlib.pylab.figure(1)
        gs = matplotlib.gridspec.GridSpec(2,2)
        ax1 = fig1.add_subplot(gs[0,:])
        ax1.grid()
        plt.title('normalized rBV')
        plt.xlim([-10,400])
        plt.ylabel('rBV_norm')
        plt.xlabel('mean_number_of_root_nodes')
        
        ax2 = fig1.add_subplot(gs[1,0])
        ax2.grid()
        plt.xlim([-10,400])
        plt.ylabel('rBV_norm (only arterial)')
        plt.xlabel('number_of_a_root_nodes')
        ax3 = fig1.add_subplot(gs[1,1])
        ax3.grid()
        plt.xlim([-10,400])
        plt.ylabel('rBV_norm (only venous)')
        plt.xlabel('number_of_v_root_nodes')
        
        fig2 = matplotlib.pylab.figure(2)
        axtwo1 = fig2.add_subplot(gs[0,:])
        axtwo1.grid()
        #plt.title('normalized rBV')
        plt.ylabel('rBF [ 1/min]')
        plt.xlabel('rBV_norm (only arterial)')
        axtwo2 = fig2.add_subplot(gs[1,0])
        axtwo2.grid()
        #plt.title('normalized rBV')
        plt.ylabel('rBF [ 1/min]')
        plt.xlabel('rBV_norm (only venous)')
        axtwo3 = fig2.add_subplot(gs[1,1])
        axtwo3.grid()
        #plt.title('normalized rBV')
        plt.ylabel('rBF [ 1/min]')
        plt.xlabel('rBV_norm (mean)')
        for t in 'typeA typeB typeC typeD typeE typeF typeG typeH typeI'.split():
        #for t in 'typeE typeG typeH'.split():
        #for t in 'typeA'.split():
            print(t)
            groups_of_that_type = []
            number_of_arteries = []
            number_of_veins = []
            number_of_circulated = []
            number_of_a_roots = []
            number_of_v_roots = []
            cnt_mean_roots = []
            TV_of_type = []
            TV_a_of_type = []
            TV_v_of_type = []
            rBV_of_type = []
            rBV_a_of_type = []
            rBV_v_of_type = []
            rBV_of_type_normalized = []
            rBV_a_of_type_normalized = []
            rBV_v_of_type_normalized = []
            TF_of_type = []
            rBF_of_type = []
            for key in f_measure.keys():
                if t in key:
                    groups_of_that_type.append(key)
                    number_of_arteries.append(f_measure[key].attrs['#artery'])
                    number_of_veins.append(f_measure[key].attrs['#veins'])
                    number_of_circulated.append(f_measure[key].attrs['#circulated'])
                    
                    p_in = np.asarray(f_measure[key+'/a_roots_pressure'])
                    p_out = np.asarray(f_measure[key+'/v_roots_pressure'])                
                    
                    number_of_a_roots.append(len(p_in))
                    number_of_v_roots.append(len(p_out))
                    mean_R = (len(p_in)+len(p_out))/2
                    cnt_mean_roots.append(mean_R)
                    
                    #according to equation D6 
                    alpha = 3.
                    
                    N = f_measure[key].attrs['#circulated']/mean_R/2
                    N_a = f_measure[key].attrs['#artery']/len(p_in)
                    N_v = f_measure[key].attrs['#veins']/len(p_out)
                    l = f_measure[key].attrs['length']
                    r_0 = f_measure[key].attrs['c_radius']
                    rho_vessel = f_measure[key].attrs['#circulated']/512e9
                    V_single_cap = np.pi*r_0**2*l   
                    TV = 2*np.pi*l*r_0**2*(((N+1)-np.power(N+1,2./alpha))/(2-np.power(2,2./alpha)))
                    f_measure[key].attrs['TV'] = TV
                    rBV = TV /(512e9)
                    f_measure[key].attrs['rBV'] = rBV
                    rBV_of_type.append(rBV)
                    rBV_normalized = rBV/(rho_vessel*V_single_cap)
                    rBV_of_type_normalized.append(rBV_normalized)
                    
                    TV_a = 2*np.pi*l*r_0**2*(((N_a+1)-np.power(N_a+1,2./alpha))/(2-np.power(2,2./alpha)))
                    TV_a_of_type.append(TV_a)
                    rBV_a = TV_a /(512e9)
                    f_measure[key].attrs['rBV_a'] = rBV_a
                    rBV_a_of_type.append(rBV_a)
                    rBV_a_normalized = rBV_a/(rho_vessel*V_single_cap)
                    rBV_a_of_type_normalized.append(rBV_a_normalized)
                    
                    TV_v = 2*np.pi*l*r_0**2*(((N_v+1)-np.power(N_v+1,2./alpha))/(2-np.power(2,2./alpha)))
                    TV_v_of_type.append(TV_v)
                    rBV_v = TV_v /(512e9)
                    f_measure[key].attrs['rBV_v'] = rBV_v
                    rBV_v_of_type.append(rBV_v)
                    rBV_v_normalized = rBV_v/(rho_vessel*V_single_cap)
                    rBV_v_of_type_normalized.append(rBV_v_normalized)
                    
                    #now flow
                    #using flow weighted average pressure
                    
                    if (len(p_in)!=len(p_out)):
                        print("Warning: not equal roots")
                    TF = 1/2.
                    TF = TF*(np.power(2,4/alpha)-2)/(2*np.power(2,4/alpha))
                    TF = TF*((N+1)*np.power(N+1,4/alpha))/(np.power(N+1,4/alpha)-(N+1))
    
                    avgVenousPressure = f_measure[key].attrs['avgVenousPressure']
                    avgArterialPressure = f_measure[key].attrs['avgArterialPressure']
                    
                    TF = TF*np.pi*r_0**4/(8*8*4e-6*l)*( avgArterialPressure-avgVenousPressure)*60. #from sec to min
                    #TF = 1/2.*(np.power(2,4/alpha)-2)/(2*np.power(2,4/alpha))*((N+1)*np.power(N+1,4/alpha))/(np.power(N+1,4/alpha)-(N+1))*np.pi*r_0**4/(8*4e-6*l)*(artery_roots_pressure[0]-venous_roots_pressure[0])
                    f_measure[key].attrs['TF'] = TF
                    f_measure[key].attrs['rBF'] = TF/(512e9)
                    TV_of_type.append(TV)
                    
                    TF_of_type.append(TF)
                    rBF_of_type.append(TF/(512e9))
            #groups_of_that_type = [value for key, value in f_measure.items() if t in key]
            f_measure['/means'].create_group(t)
            f_measure['/means/' + t].attrs['artery'] = np.mean(number_of_arteries)
            f_measure['/means/' + t].attrs['artery_std'] = np.std(number_of_arteries)
            f_measure['/means/' + t].attrs['veins'] = np.mean(number_of_veins)
            f_measure['/means/' + t].attrs['veins_std'] = np.std(number_of_veins)
            f_measure['/means/' + t].attrs['circulated'] = np.mean(number_of_circulated)
            f_measure['/means/' + t].attrs['circulated_std'] = np.std(number_of_circulated)
            f_measure['/means/' + t].attrs['TV'] = np.mean(TV_of_type)
            f_measure['/means/' + t].attrs['TV_std'] = np.std(TV_of_type)
            f_measure['/means/' + t].attrs['rBV'] = np.mean(rBV_of_type)
            f_measure['/means/' + t].attrs['rBV_std'] = np.std(rBV_of_type)
            f_measure['/means/' + t].attrs['TF'] = np.mean(TF_of_type)
            f_measure['/means/' + t].attrs['TF_std'] = np.std(TF_of_type)
            f_measure['/means/' + t].attrs['rBF'] = np.mean(rBF_of_type)
            f_measure['/means/' + t].attrs['rBF_std'] = np.std(rBF_of_type)
            f_measure['/means/' + t].attrs['number_of_a_roots'] = np.mean(number_of_a_roots)
            f_measure['/means/' + t].attrs['number_of_v_roots'] = np.mean(number_of_v_roots)
            f_measure['/means/' + t].attrs['cnt_mean_roots'] = np.mean(cnt_mean_roots)
            
            matplotlib.pylab.figure(1)
            ax1.scatter(cnt_mean_roots,rBV_of_type_normalized,marker = '<>*osd^+x'[color_index], c = colors[color_index], facecolor='none', edgecolor = colors[color_index], linewidth=1., s = 50)
            ax1.set_yscale('log')
            #ax1.set_xscale('log')        
            ax2.scatter(number_of_a_roots,rBV_a_of_type_normalized,marker = '<>*osd^+x'[color_index], c = colors[color_index], facecolor='none', edgecolor = colors[color_index], linewidth=1., s = 50)
            ax2.set_yscale('log')
            #ax2.set_xscale('log')         
            ax3.scatter(number_of_v_roots,rBV_v_of_type_normalized,marker = '<>*osd^+x'[color_index], c = colors[color_index], facecolor='none', edgecolor = colors[color_index], linewidth=1., s = 50)
            ax3.set_yscale('log')
            #ax3.set_xscale('log')         
            matplotlib.pylab.figure(2)    
            axtwo1.scatter(rBV_of_type_normalized,rBF_of_type,marker = '<>*osd^+x'[color_index], c = colors[color_index], facecolor='none', edgecolor = colors[color_index], linewidth=1., s = 50)        
            axtwo1.set_yscale('log')
            axtwo1.set_xscale('log')
            axtwo2.scatter(rBV_a_of_type_normalized,rBF_of_type,marker = '<>*osd^+x'[color_index], c = colors[color_index], facecolor='none', edgecolor = colors[color_index], linewidth=1., s = 50)        
            axtwo2.set_yscale('log')
            axtwo2.set_xscale('log')        
            axtwo3.scatter(rBV_v_of_type_normalized,rBF_of_type,marker = '<>*osd^+x'[color_index], c = colors[color_index], facecolor='none', edgecolor = colors[color_index], linewidth=1., s = 50)                
            axtwo3.set_yscale('log')
            axtwo3.set_xscale('log')        
            color_index = color_index + 1
        
        matplotlib.pylab.figure(1)
        fig1.tight_layout()
        #ax1.legend('RC1 RC2 RC3 RC4 RC5 RC6 RC7 RC8 RC9'.split(), loc="center right")
        ax1.legend('RC1 RC2 RC3 RC4 RC5 RC6 RC7 RC8 RC9'.split(), bbox_to_anchor=(0.75, 1), loc=2, borderaxespad=0.)
        plt.savefig("rBV_vs_roots.png")
        
        matplotlib.pylab.figure(2)
        fig2.tight_layout()
        axtwo1.legend('RC1 RC2 RC3 RC4 RC5 RC6 RC7 RC8 RC9'.split(), bbox_to_anchor=(0.9, 0.9), loc=2, borderaxespad=0.)
        plt.savefig("rBF_vs_rBV.png")
    if 1:
        #gamma fitting
        x_0=[]
        x_1=[]
        x_2=[]
        x_3=[]
        rBV_for_fit_gamma = []
        rBV_for_fit_gamma_normalized = []
        fig1 = matplotlib.pylab.figure(1)
        gs1 = matplotlib.gridspec.GridSpec(2,2)
        ax1 = fig1.add_subplot(gs1[0,:])
        ax1.grid()
        plt.title('rBV (normalized to single capillary volume)')
        plt.xlim([-10,400])
        plt.ylabel('rBV_norm')
        plt.xlabel('mean #roots (#a_root+#v_root)/2')
        
        ax2 = fig1.add_subplot(gs1[1,0])
        ax2.grid()
        plt.xlim([-10,400])
        plt.ylabel('rBV_norm (arterial)')
        plt.xlabel('number_of_a_root_nodes')
        ax3 = fig1.add_subplot(gs1[1,1])
        ax3.grid()
        plt.xlim([-10,400])
        plt.ylabel('rBV_norm (venous)')
        plt.xlabel('number_of_v_root_nodes')
        
        fig2 = matplotlib.pylab.figure(2)
        gs = matplotlib.gridspec.GridSpec(2,2)
        axtwo1 = fig2.add_subplot(gs[0,:])
        axtwo1.grid()
        plt.ylim([0.01,0.06])
        plt.xlim([0.0069,0.00758])
        #plt.title('normalized rBV')
        plt.ylabel('rBF [ 1/min]')
        plt.xlabel('rBV')
       # plt.yscale('log')
        #plt.xscale('log')
        axtwo2 = fig2.add_subplot(gs[1,0])
        axtwo2.grid()
        plt.ylim([0.01,0.06])
        plt.xlim([0.0022,0.0028])
        #plt.title('normalized rBV')
        plt.ylabel('rBF [ 1/min]')
        plt.xlabel('rBV (only artery)')
        axtwo3 = fig2.add_subplot(gs[1,1])
        axtwo3.grid()
        plt.ylim([0.01,0.06])
        plt.xlim([0.0045,0.005])
        #plt.title('normalized rBV')
        plt.ylabel('rBF [ 1/min]')
        plt.xlabel('rBV (only venous)')
        
        fig3 = matplotlib.pylab.figure(3)
#        axthree = fig3.add_plot()
#        axthree.grid()
        plt.ylabel('rBF [1/min]')
        plt.xlabel('mean #roots (#a_root+#v_root)/2')
        
        for t in 'typeA typeB typeC typeD typeE typeF typeG typeH typeI'.split():
        #for t in 'typeE typeG typeH'.split():
        #for t in 'typeA'.split():
            print(t)
            velocity_of_that_type = []
            velocity_of_that_type_std = []
            flow_of_that_type = []
            flow_of_that_type_std = []
            groups_of_that_type = []
            number_of_arteries = []
            number_of_veins = []
            number_of_circulated = []
            number_of_a_roots = []
            number_of_v_roots = []
            cnt_total_roots = []
            cnt_mean_roots = []
            TV_of_type = []
            TV_a_of_type = []
            TV_v_of_type = []
            TV_cap_of_type = []
            rBV_of_type = []
            rBV_a_of_type = []
            rBV_v_of_type = []
            rBV_cap_of_type = []
            rBV_of_type_normalized = []
            rBV_a_of_type_normalized = []
            rBV_v_of_type_normalized = []
            rBV_cap_of_type_normalized = []
            TF_of_type = []
            rBF_of_type = []
            for key in f_measure.keys():
                if t in key:
                    groups_of_that_type.append(key)
                    number_of_arteries.append(f_measure[key].attrs['#artery'])
                    number_of_veins.append(f_measure[key].attrs['#veins'])
                    number_of_circulated.append(f_measure[key].attrs['#circulated'])
                    
                    current_radii = np.asarray(f_measure[key+'/radii_all'])
                    current_radii_artery = np.asarray(f_measure[key+'/radii_of_arteries'])
                    current_radii_vein = np.asarray(f_measure[key+'/radii_of_veins'])
                    current_radii_capillary = np.asarray(f_measure[key+'/radii_of_capillaries'])
                    
                    current_volume = np.asarray(f_measure[key+'/volumes_all'])
                    current_volume_a = np.asarray(f_measure[key+'/volumes_of_arteries'])
                    current_volume_v = np.asarray(f_measure[key+'/volumes_of_veins'])
                    current_volume_cap = np.asarray(f_measure[key+'/volumes_of_capillaries'])
                    
                    
                    current_flow = np.asarray(f_measure[key+'/flow_all'])
                    current_flow_artery = np.asarray(f_measure[key+'/flow_of_arteries'])
                    current_flow_vein = np.asarray(f_measure[key+'/flow_of_veins'])
                    current_flow_capillary = np.asarray(f_measure[key+'/flow_of_capillaries'])
                    
                    current_velocities_all = current_flow/(np.pi*np.power(current_radii,2))
                    if 0:
                        current_velocities_all = current_flow/(np.pi*np.power(current_radii,2))
                        velocity_of_that_type.append(np.mean(current_velocities_all))
                        velocity_of_that_type_std.append(np.std(current_velocities_all))
                    if 0:
                        current_velocities_cap = current_flow_capillary/(np.pi*np.power(current_radii_capillary,2))
                        velocity_of_that_type.append(np.mean(current_velocities_cap))
                        velocity_of_that_type_std.append(np.std(current_velocities_cap))
                    if 0:
                        flow_of_that_type.append(np.mean(current_flow))
                        flow_of_that_type_std.append(np.std(current_flow))
                    if 1:
                        flow_of_that_type.append(np.mean(current_flow_capillary))
                        flow_of_that_type_std.append(np.std(current_flow_capillary))
                    p_in = np.asarray(f_measure[key+'/a_roots_pressure'])
                    p_out = np.asarray(f_measure[key+'/v_roots_pressure'])                
                    
                    number_of_a_roots.append(len(p_in))
                    number_of_v_roots.append(len(p_out))
                    mean_R = (len(p_in)+len(p_out))/2
                    total_no_roots = len(p_in)+len(p_out)
                    cnt_mean_roots.append(mean_R)
                    cnt_total_roots.append(total_no_roots)
                    
                    #according to equation D6 
                    alpha = 3.
                    
                    N = f_measure[key].attrs['#circulated']/mean_R/2
                    N_a = f_measure[key].attrs['#artery']/len(p_in)
                    N_v = f_measure[key].attrs['#veins']/len(p_out)
                    l = f_measure[key].attrs['length']
                    r_0 = f_measure[key].attrs['c_radius']
                    rho_vessel = f_measure[key].attrs['#circulated']/512e9
                    x_0.append(rho_vessel)
                    x_1.append(f_measure[key].attrs['#circulated'])
                    x_2.append(len(p_in)+len(p_out))
                    x_3.append(N)
                    V_single_cap = np.pi*r_0**2*l
                    
                    #current_TV = [ 2*np.pi*r**2*l for r in current_radii] #calculate tissue volume from radii
                    TV = np.sum(current_volume)
                    TV_of_type.append(TV)
                    f_measure[key].attrs['TV'] = TV
                    rBV = TV /(512e9)
                    rBV_for_fit_gamma.append(rBV)
                    f_measure[key].attrs['rBV'] = rBV
                    rBV_of_type.append(rBV)
                    rBV_normalized = rBV/(rho_vessel*V_single_cap)
                    rBV_of_type_normalized.append(rBV_normalized)
                    rBV_for_fit_gamma_normalized.append(rBV_normalized)
                    
                    #current_TV_a = [ 2*np.pi*r**2*l for r in current_radii_artery]
                    TV_a = np.sum(current_volume_a)
                    TV_a_of_type.append(TV_a)
                    rBV_a = TV_a /(512e9)
                    f_measure[key].attrs['rBV_a'] = rBV_a
                    rBV_a_of_type.append(rBV_a)
                    rBV_a_normalized = rBV_a/(rho_vessel*V_single_cap)
                    rBV_a_of_type_normalized.append(rBV_a_normalized)
                    
                    #current_TV_v = [ 2*np.pi*r**2*l for r in current_radii_vein]
                    TV_v = np.sum(current_volume_v)
                    TV_v_of_type.append(TV_v)
                    rBV_v = TV_v /(512e9)
                    f_measure[key].attrs['rBV_v'] = rBV_v
                    rBV_v_of_type.append(rBV_v)
                    rBV_v_normalized = rBV_v/(rho_vessel*V_single_cap)
                    rBV_v_of_type_normalized.append(rBV_v_normalized)
                    
                    #current_TV_cap = [ 2*np.pi*r**2*l for r in current_radii_capillary]
                    TV_cap = np.sum(current_volume_cap)
                    TV_cap_of_type.append(TV_cap)
                    rBV_cap = TV_cap /(512e9)
                    f_measure[key].attrs['rBV_cap'] = rBV_cap
                    rBV_cap_of_type.append(rBV_cap)
                    rBV_cap_normalized = rBV_v/(rho_vessel*V_single_cap)
                    rBV_cap_of_type_normalized.append(rBV_cap_normalized)
                    
                    #now flow
                    #using flow weighted average pressure
                    
#                    if (len(p_in)!=len(p_out)):
#                        print("Warning: not equal roots")
#                    TF = 1/2.
#                    TF = TF*(np.power(2,4/alpha)-2)/(2*np.power(2,4/alpha))
#                    TF = TF*((N+1)*np.power(N+1,4/alpha))/(np.power(N+1,4/alpha)-(N+1))
    
                    avgVenousPressure = f_measure[key].attrs['avgVenousPressure']
                    avgArterialPressure = f_measure[key].attrs['avgArterialPressure']
                    
#                    TF = TF*np.pi*r_0**4/(8*8*4e-6*l)*( avgArterialPressure-avgVenousPressure)*60. #from sec to min
                    #TF = 1/2.*(np.power(2,4/alpha)-2)/(2*np.power(2,4/alpha))*((N+1)*np.power(N+1,4/alpha))/(np.power(N+1,4/alpha)-(N+1))*np.pi*r_0**4/(8*4e-6*l)*(artery_roots_pressure[0]-venous_roots_pressure[0])
#                    current_flow = np.asarray(f_measure[key+'/flow_all'])
#                    current_flow_artery = np.asarray(f_measure[key+'/flow_of_arteries'])
#                    current_flow_vein = np.asarray(f_measure[key+'/flow_of_veins'])
#                    current_flow_capillary = np.asarray(f_measure[key+'/flow_of_capillaries'])                    
                    TF = np.sum(current_flow)
                    f_measure[key].attrs['TF'] = TF
                    f_measure[key].attrs['rBF'] = TF/(512e9)
                    
                    
                    TF_of_type.append(TF)
                    rBF_of_type.append(TF/(512e9))
            #groups_of_that_type = [value for key, value in f_measure.items() if t in key]
            f_measure['/means'].create_group(t)
            f_measure['/means/' + t].attrs['artery'] = np.mean(number_of_arteries)
            f_measure['/means/' + t].attrs['artery_std'] = np.std(number_of_arteries)
            f_measure['/means/' + t].attrs['veins'] = np.mean(number_of_veins)
            f_measure['/means/' + t].attrs['veins_std'] = np.std(number_of_veins)
            f_measure['/means/' + t].attrs['circulated'] = np.mean(number_of_circulated)
            f_measure['/means/' + t].attrs['circulated_std'] = np.std(number_of_circulated)
            f_measure['/means/' + t].attrs['TV'] = np.mean(TV_of_type)
            f_measure['/means/' + t].attrs['TV_std'] = np.std(TV_of_type)
            f_measure['/means/' + t].attrs['rBV'] = np.mean(rBV_of_type)
            f_measure['/means/' + t].attrs['rBV_std'] = np.std(rBV_of_type)
            f_measure['/means/' + t].attrs['TF'] = np.mean(TF_of_type)
            f_measure['/means/' + t].attrs['TF_std'] = np.std(TF_of_type)
            f_measure['/means/' + t].attrs['rBF'] = np.mean(rBF_of_type)
            f_measure['/means/' + t].attrs['rBF_std'] = np.std(rBF_of_type)
            f_measure['/means/' + t].attrs['number_of_a_roots'] = np.mean(number_of_a_roots)
            f_measure['/means/' + t].attrs['number_of_v_roots'] = np.mean(number_of_v_roots)
            f_measure['/means/' + t].attrs['cnt_mean_roots'] = np.mean(cnt_mean_roots)
            
            matplotlib.pylab.figure(1)
            ax1.scatter(cnt_mean_roots,rBV_of_type_normalized,marker = '<>*osd^+x'[color_index], c = colors[color_index], facecolor='none', edgecolor = colors[color_index], linewidth=1., s = 50)
            #ax1.set_yscale('log')
            #ax1.set_xscale('log')
            if 0:
                rBV_from_gamma= [a for a in cnt_mean_roots]
                ax1.plot(cnt_mean_roots,rBV_from_gamma)
            ax2.scatter(number_of_a_roots,rBV_a_of_type_normalized,marker = '<>*osd^+x'[color_index], c = colors[color_index], facecolor='none', edgecolor = colors[color_index], linewidth=1., s = 50)
            #ax2.set_yscale('log')
            #ax2.set_xscale('log')         
            ax3.scatter(number_of_v_roots,rBV_v_of_type_normalized,marker = '<>*osd^+x'[color_index], c = colors[color_index], facecolor='none', edgecolor = colors[color_index], linewidth=1., s = 50)
            #ax3.set_yscale('log')
            #ax3.set_xscale('log')         
            matplotlib.pylab.figure(2)    
            axtwo1.scatter(rBV_of_type,rBF_of_type,marker = '<>*osd^+x'[color_index], c = colors[color_index], facecolor='none', edgecolor = colors[color_index], linewidth=1., s = 50)        
                        
            #axtwo1.set_yscale('log')
            #axtwo1.set_xscale('log')
            axtwo2.scatter(rBV_a_of_type,rBF_of_type,marker = '<>*osd^+x'[color_index], c = colors[color_index], facecolor='none', edgecolor = colors[color_index], linewidth=1., s = 50)        
            #axtwo2.set_yscale('log')
            #axtwo2.set_xscale('log')
            axtwo3.scatter(rBV_v_of_type,rBF_of_type,marker = '<>*osd^+x'[color_index], c = colors[color_index], facecolor='none', edgecolor = colors[color_index], linewidth=1., s = 50)                
            #axtwo3.set_yscale('log')
            #axtwo3.set_xscale('log')
            
            matplotlib.pylab.figure(3)
            plt.scatter(cnt_mean_roots,rBF_of_type,marker = '<>*osd^+x'[color_index], c = colors[color_index], facecolor='none', edgecolor = colors[color_index], linewidth=1., s = 50)

            if 0:
                matplotlib.pylab.figure(4)
                plt.errorbar(range(len(velocity_of_that_type)),velocity_of_that_type, yerr = velocity_of_that_type_std, marker= '<>*osd^+x'[color_index], c = colors[color_index])        
                plt.xlabel("Index of sample")
                plt.ylabel(r"Blood flow velocity/ $\mu m/sec$")
            
            if 1:
                matplotlib.pylab.figure(4)
                plt.errorbar(range(len(flow_of_that_type)),flow_of_that_type, yerr = flow_of_that_type_std, marker= '<>*osd^+x'[color_index], c = colors[color_index])        
                plt.xlabel("Index of sample")
                plt.ylabel(r"Blood flow / $(\mu m)^s3/sec$")
                color_index = color_index + 1
        
        matplotlib.pylab.figure(1)
        fig1.tight_layout()
        matplotlib.pyplot.autoscale(enable=True, axis=u'both', tight=None)
        #ax1.legend('RC1 RC2 RC3 RC4 RC5 RC6 RC7 RC8 RC9'.split(), loc="center right")
        ax1.legend('RC1 RC2 RC3 RC4 RC5 RC6 RC7 RC8 RC9'.split(), bbox_to_anchor=(0.9, 1.1), loc=2, borderaxespad=0., fontsize= 'x-small')
        
        
        matplotlib.pylab.figure(2)
        fig2.tight_layout()
        #matplotlib.pyplot.autoscale(enable=True, axis=u'both', tight=None)
        axtwo1.legend('RC1 RC2 RC3 RC4 RC5 RC6 RC7 RC8 RC9'.split(), bbox_to_anchor=(0.9, 0.7), loc=2, borderaxespad=0., fontsize= 'x-small')
        plt.savefig("rBF_vs_rBV.png")
        
        matplotlib.pylab.figure(3)
        plt.legend('RC1 RC2 RC3 RC4 RC5 RC6 RC7 RC8 RC9'.split(), bbox_to_anchor=(0.9, 0.7), loc=2, borderaxespad=0., fontsize= 'x-small')        
        plt.grid()        
        plt.savefig("rBF_vs_root.png")
        
        #rieger review
        matplotlib.pylab.figure(4)
        
        plt.legend('RC1 RC2 RC3 RC4 RC5 RC6 RC7 RC8 RC9'.split(), bbox_to_anchor=(0.9, 0.7), loc=2, borderaxespad=0., fontsize= 'x-small')
        plt.grid()        
        plt.savefig("velocities.png")        
        
        #fitting
        from scipy.optimize import curve_fit

        #x_0 is vessel density
        #x_1 is M number circulated
        #x_2 is number of root nodes        
        #N>>1
        def my_fit_func(x,p_0):
            first = np.pi*np.power(2.5,2)*130.*x[0]/2
            second = 1-(np.power(x[1],p_0)/x[1])*(2*x[2]/(np.power(2*x[2],p_0)))
            third = 1.-(np.power(2,p_0)/2)
            
            return (first*second)/third
        #all N
        #x_3 is N
        def my_fit_func2(x,p_0):
            first = (1+1/x[0])
            second = 1-(np.power(x[0]+1,p_0)/(x[0]+1))
            third = 2-np.power(2,p_0)
            
            return (first*second)/third
        
        p0=0.8
        
        #for my_fit_func
        #x_4d = np.asarray([x_0,x_1,x_2])
        
        #for my_fit_func2
        x_4d = np.asarray([x_3])
        
        #fitParams, fitCovariance = curve_fit(my_fit_func,xdata=x_4d,ydata=rBV_for_fit_gamma, p0=p0)
        fitParams, fitCovariance = curve_fit(my_fit_func2,xdata=x_4d,ydata=rBV_for_fit_gamma_normalized, p0=p0)
        print( ' fit coefficients:\n', fitParams )
        print( ' fit Covariance:\n', fitCovariance )
        #use data from Rinneberg 28.8.2015 plot curve 
        R_Rinneberg = np.arange(1,400,20)
        N_Rinneberg = [260000/(2*thisR) for thisR in R_Rinneberg]
        rBV_norm_Rinneberg = [(1+1/NN)*((1-(np.power(NN+1,fitParams[0]))/(NN+1)))/(2-np.power(2,fitParams[0])) for NN in N_Rinneberg]
        rBV_norm_Rinneberg2 = [(1+1/NN)/(2*np.log(2))*np.log(NN+1) for NN in N_Rinneberg]
        #rBV_norm_Rinneberg = [(1 + 1/NN)*(1 - ((M/(2*float(R)) + 1)**(0.7)/((M/(2*float(R))) + 1)))/(2 - 2**(0.7)) for R in R_Rinneberg]
        matplotlib.pylab.figure(1)
        ax1.scatter(R_Rinneberg,rBV_norm_Rinneberg, c = 'b', facecolor='none', linewidth=1., s = 50)
        plt.savefig("rBV_vs_roots.png")
        #plt.show()
           
        
if __name__ == "__main__":
  import optparse
  parser = optparse.OptionParser()
  parser.add_option("-O","--with-o2", dest="with_o2", help="look at detailed o2 data", default=False, action="store_true")
  options, args = parser.parse_args()

  filenames, pattern = args[:-1], args[-1]
  fn_measure = "vs_compartment_v2.h5"
  #DoReadIn(filenames, pattern,  fn_measure)
  DoCalculate(fn_measure)
