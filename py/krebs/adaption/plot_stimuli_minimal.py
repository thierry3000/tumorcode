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

import matplotlib.pyplot as plt
import h5py
import itertools
import numpy as np
import myutils
import os
from matplotlib.backends.backend_pdf import PdfPages
import krebsutils as ku


def plot_hydorodynamic_charicteristics(vessel_grp,pp):
    shearforce = vessel_grp['edges/shearforce']
    shearforce = np.multiply(shearforce,10000)
    diameter = np.multiply(vessel_grp['edges/radius'],2)
    node_a_index =vessel_grp['edges/node_a_index']
    node_b_index =vessel_grp['edges/node_b_index']
        
    pressure_at_nodes = vessel_grp['nodes/pressure']
    
    
    pressure_at_vessel=[]
        
    for NodeA, NodeB in itertools.izip(node_a_index,node_b_index):
        pressure_at_vessel.append((pressure_at_nodes[NodeA]+pressure_at_nodes[NodeB])/2 * 1/0.1333)
    fig = plt.figure()
    plt.subplot(2,1,1)
    plt.semilogy(pressure_at_vessel,shearforce,'*')
    plt.grid()
    plt.xlabel('pressure/ mmHg')
    plt.ylabel('shearstress dyne/cm^2')
    plt.subplot(2,1,2)
    plt.semilogy(pressure_at_vessel,diameter,'r*')
    plt.grid()
    plt.xlabel('pressure/ mmHg')
    plt.ylabel('diameter/ mum')
    pp.savefig()
    plt.show()

def plot_hydrodynamic_stimuli(vessel_grp,pp):
    shearforce = vessel_grp['edges/shearforce']
    shearforce = np.multiply(shearforce,10000) #kpa to dyne
    node_a_index =vessel_grp['edges/node_a_index']
    node_b_index =vessel_grp['edges/node_b_index']
        
    pressure_at_nodes = vessel_grp['nodes/pressure']
    
    
    pressure_at_vessel=[]
        
    for NodeA, NodeB in itertools.izip(node_a_index,node_b_index):
        pressure_at_vessel.append((pressure_at_nodes[NodeA]+pressure_at_nodes[NodeB])/2 * 1/0.1333)
    pressure_stimuli = 100 - 86 *np.exp(-5000*np.log10(np.log10(pressure_at_vessel))**5.4)
    fig = plt.figure()
    plt.semilogx(pressure_at_vessel,np.log10(shearforce),'*')
    plt.semilogx(pressure_at_vessel,-np.log10(pressure_stimuli),'r*')
    plt.xlim([10,100])
    plt.ylim([-2.2,4])
    plt.legend(['Hydrodynamic stimuli', 'Pressure'])
    plt.grid()
    plt.xlabel('pressure/ mmHg')
    plt.ylabel('stimuli')
    pp.savefig()
    plt.show() 
    
def plot_conductive_stimuli(adaption_grp,pp):
    metabolic = adaption_grp['edges/metabolicSignal']
    conductive = adaption_grp['edges/conductivitySignal']
    flow = adaption_grp['edges/flow']
    flow = np.multiply(flow,60./1000000.)
    
    node_a_index =adaption_grp['edges/node_a_index']
    node_b_index =adaption_grp['edges/node_b_index']
    
    pressure_at_nodes = adaption_grp['nodes/pressure']
    
    pressure_at_vessel=[]
    
    for NodeA, NodeB in itertools.izip(node_a_index,node_b_index):
        pressure_at_vessel.append((pressure_at_nodes[NodeA]+pressure_at_nodes[NodeB])/2 * 1/0.1333)
    
    fig2 = plt.figure()
    plt.title('conductive stimuli')
    plt.semilogx(flow,metabolic,'*r')
    plt.semilogx(flow,conductive,'*b')
    plt.legend(['metabolic','conductive'])
    plt.grid()
    plt.xlabel("flow/nl/min")
    pp.savefig()
    plt.show()
    
def plot_wall_shear_stress(vessel_grp,pp):
    index_of_artery = np.bitwise_and(np.asarray(vessel_grp['edges/flags']), ku.ARTERY)>0
    index_of_capillary = np.bitwise_and(np.asarray(vessel_grp['edges/flags']), ku.CAPILLARY)>0
    index_of_vein = np.bitwise_and(np.asarray(vessel_grp['edges/flags']), ku.VEIN)>0
    shearforce = vessel_grp['edges/shearforce']
    shearforce = np.multiply(shearforce,10000) #kpa to dyne
    node_a_index =vessel_grp['edges/node_a_index']
    node_b_index =vessel_grp['edges/node_b_index']
        
    pressure_at_nodes = vessel_grp['nodes/pressure']
    
    
    pressure_at_vessel=[]
        
    for NodeA, NodeB in itertools.izip(node_a_index,node_b_index):
        pressure_at_vessel.append((pressure_at_nodes[NodeA]+pressure_at_nodes[NodeB])/2 * 1/0.1333)
    fig = plt.figure()
    pressure_at_vessel = np.asarray(pressure_at_vessel)
    plt.loglog(pressure_at_vessel[index_of_artery],shearforce[index_of_artery],'r*')
    plt.loglog(pressure_at_vessel[index_of_capillary],shearforce[index_of_capillary],'y*')
    plt.loglog(pressure_at_vessel[index_of_vein],shearforce[index_of_vein],'b*')
    plt.legend(['ART','CAP','VEN'])
    plt.grid()
    plt.xlabel('pressure/ mmHg')
    plt.ylabel('wall shear stress/ dyne/cm^2')
    plt.xlim([10,100])
    plt.ylim([1,1000])
    pp.savefig()
    plt.show()
    
def plot_diameters(vessel_grp,pp):
    index_of_artery = np.bitwise_and(np.asarray(vessel_grp['edges/flags']), ku.ARTERY)>0
    index_of_capillary = np.bitwise_and(np.asarray(vessel_grp['edges/flags']), ku.CAPILLARY)>0
    index_of_vein = np.bitwise_and(np.asarray(vessel_grp['edges/flags']), ku.VEIN)>0
    radii = vessel_grp['edges/radius']
    diameter = np.multiply(radii,2)
    node_a_index =vessel_grp['edges/node_a_index']
    node_b_index =vessel_grp['edges/node_b_index']
        
    pressure_at_nodes = vessel_grp['nodes/pressure']
    
    
    pressure_at_vessel=[]
        
    for NodeA, NodeB in itertools.izip(node_a_index,node_b_index):
        pressure_at_vessel.append((pressure_at_nodes[NodeA]+pressure_at_nodes[NodeB])/2 * 1/0.1333)
    pressure_at_vessel = np.asarray(pressure_at_vessel)
    fig = plt.figure()
    plt.loglog(pressure_at_vessel[index_of_artery],diameter[index_of_artery],'r*')
    plt.loglog(pressure_at_vessel[index_of_capillary],diameter[index_of_capillary],'y*')
    plt.loglog(pressure_at_vessel[index_of_vein],diameter[index_of_vein],'b*')
    plt.legend(['ART', 'CAP', 'VEN'])
    plt.xlim([10,100])
    plt.ylim([3,100])
    plt.grid()
    plt.xlabel('pressure/ mmHg')
    plt.ylabel('diameter/mum')
    pp.savefig()
    plt.show() 
    
if __name__ == '__main__':
#    filename1= 'symetric_world_adption_p_adaption_symetric_world.h5'
#    filename2= 'Asymetric_world_adption_p_adaption_asymetric_world.h5'
#    f_sym1 = h5py.File(filename1,'r')
#    f_asym = h5py.File(filename2,'r')
#    
#    vesselgrp_sym1 = f_sym1['/adaption/vessels_after_adaption_40']
#    vesselgrp_asym = f_asym['/adaption/vessels_after_adaption_238']
#    common_filename = 'minimal_plots'
#    
#    pp = PdfPages( common_filename + '_stimulies.pdf')
#    
#    plot_wall_shear_stress(vesselgrp_sym1, pp)    
#    plot_diameters(vesselgrp_sym1,pp)
#    
#    plot_wall_shear_stress(vesselgrp_asym, pp)    
#    plot_diameters(vesselgrp_asym,pp)
#    
#    #hydrodynamic_fig = plot_hydrodynamic_stimuli(vesselgrp, pp)
#    #plot_conductive_stimuli(vesselgrp,pp)
#    #plot_hydorodynamic_charicteristics(vesselgrp,pp)    
#    f_sym1.close()
#    f_asym.close()
#    pp.close()
    
    
    filename1= '/localdisk/thierry/output_adaption/secomb/mesentry_secomb_546_adption_p_adaption_default.h5'
    f = h5py.File(filename1,'r')
    
    
    vesselgrp = f['/adaption/vessels_after_adaption']
    vesselgrp = f['/vessels/recomputed_flow']
    
    common_filename = 'minimal_plots'
    
    pp = PdfPages( common_filename + '_stimulies.pdf')
    
    plot_wall_shear_stress(vesselgrp, pp)    
    plot_diameters(vesselgrp,pp)
    plot_hydrodynamic_stimuli(vesselgrp,pp)
    
    #hydrodynamic_fig = plot_hydrodynamic_stimuli(vesselgrp, pp)
    #plot_conductive_stimuli(vesselgrp,pp)
    #plot_hydorodynamic_charicteristics(vesselgrp,pp)    
    f.close()
    
    pp.close()
    
    