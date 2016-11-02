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
    plt.loglog(pressure_at_vessel,shearforce,'*')
    plt.grid()
    plt.xlabel('pressure/ mmHg')
    plt.xlim([10,100])
    plt.ylim([1,1000])
    plt.ylabel('shearstress dyne/cm^2')
    plt.subplot(2,1,2)
    plt.loglog(pressure_at_vessel,diameter,'r*')
    plt.grid()
    plt.xlabel('pressure/ mmHg')
    plt.ylabel('diameter/ mum')
    plt.xlim([10,100])
    plt.ylim([3,100])
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
    plt.legend(['Hydrodynamic stimuli', 'Pressure'])
    plt.grid()
    plt.xlabel('pressure/ mmHg')
    plt.xlim([10,100])
    plt.ylim([-2.2,4])
    plt.ylabel('stimuli')
    pp.savefig()
    plt.show() 
    
def plot_conductive_stimuli(adaption_grp,pp):
    index_of_artery = np.bitwise_and(np.asarray(adaption_grp['edges/flags']), ku.ARTERY)>0
    index_of_capillary = np.bitwise_and(np.asarray(adaption_grp['edges/flags']), ku.CAPILLARY)>0
    index_of_vein = np.bitwise_and(np.asarray(adaption_grp['edges/flags']), ku.VEIN)>0
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
    plt.semilogx(flow[index_of_artery],metabolic[index_of_artery],'*r')
    plt.semilogx(flow[index_of_capillary],metabolic[index_of_capillary],'*y')
    plt.semilogx(flow[index_of_vein],metabolic[index_of_vein],'*b')
    
    plt.semilogx(flow[index_of_artery],conductive[index_of_artery],'Dr')
    plt.semilogx(flow[index_of_capillary],conductive[index_of_capillary],'Dy')
    plt.semilogx(flow[index_of_vein],conductive[index_of_vein],'Db')
    plt.legend(['metabolic','conductive'])
    plt.grid()
    plt.xlabel("flow/nl/min")
    plt.xlim([0.01,1000])
    plt.ylim([0,4])
    pp.savefig()
    plt.show()
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
    filename= 'mesentry_secomb_546_adption_p_mesentry_subset_vary.h5'
    filename= '/daten/localdisk/adaption_project/vessels-q2d-8mm-P6-typeE-9x3L130-sample00_adption_p_typeE.h5'
    filename= 'mesentry_secomb_546_v3_adption_p_mesentry_play.h5'    
    f = h5py.File(filename)
    no_of_iterations = '1'
    #vesselgrp = f['/adaption/vessels_after_adaption_' + no_of_iterations]
    vesselgrp = f['/adaption/vessels_after_adaption']
    #vesselgrp = f['/Asymetric/vessels_after_adaption/vessels_after_adaption']
    common_filename = os.path.splitext(filename)[0]
    
    pp = PdfPages(no_of_iterations + '_' + common_filename + '_stimulies.pdf')
    hydrodynamic_fig = plot_hydrodynamic_stimuli(vesselgrp, pp)
    plot_conductive_stimuli(vesselgrp,pp)
    plot_hydorodynamic_charicteristics(vesselgrp,pp)
    #plot_movie(f,pp=None)
    #plot_movie_typeE(f,pp=None)
    
    f.close
    #pp.close()
    
    
#    common_filename = myutils.sanitize_posixpath(splitext(commonprefix(filename))[0])
#    mpl_utils.SinglePageWriter.
#    with mpl_utils.SinglePageWriter(common_filename+'.pdf') as pdfwriter:
#        pdfwriter(hydrodynamic_fig)