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

# -*- coding: utf-8 -*-
import krebsutils as ku
import h5files
import math
import numpy as np
import scipy.spatial.distance as sdist
import sys
import os

def get_network_from_file(vesselgrp):
    filename = os.path.dirname(sys.argv[0]) + '/apj_network_546_segments.txt'
    f=open(filename,'r')
    data = f.readlines()
    cropped_edges = data[7:553]

    label=[]
    mw_vessel_flag= []
    node_a_index=[]
    node_b_index=[]
    radii = []
    length = []
    hema = []
    velocity = []
    
    #to be calculated
    flow = []

    for line in cropped_edges:
        k=1
        for value in line.split():
            if k==1:
                label.append(int(value))
            if k==2:
                flag=0
                flag=np.bitwise_or(flag,ku.CIRCULATED)
                if(int(value)==1):
                    mw_vessel_flag.append(np.bitwise_or(flag, ku.ARTERY))
                elif(int(value)==2):
                    mw_vessel_flag.append(np.bitwise_or(flag, ku.CAPILLARY))
                elif(int(value)==3):
                    mw_vessel_flag.append(np.bitwise_or(flag, ku.VEIN))
                else:
                    print("error: unkown vessel type")             
            if k==3:
                node_a_index.append(int(value))
            if k==4:
                node_b_index.append(int(value))
            if k==5:
                radii.append(float(value)/2)
            if k==6:
                length.append(float(value))
            if k==7:
                hema.append(float(value))
            if k==8:
                velocity.append(float(value))
            k=k+1
    
    for (i,v) in enumerate(velocity):
        if(v!=0):
            flow.append((np.pi*radii[i]**2)*length[i]*velocity[i]*1000)
        
    #edge stuff
    N_edges= 546
    #N_edges=2666
    edgegrp = vesselgrp.create_group("edges")
    edgegrp.attrs.create('COUNT',N_edges)
    ds_nodeA = edgegrp.create_dataset('node_a_index', data=node_a_index)
    ds_nodeB = edgegrp.create_dataset('node_b_index', data= node_b_index)
    ds_radius = edgegrp.create_dataset('radius', data=radii)
    ds_hema = edgegrp.create_dataset('hematocrit', data=hema)
    ds_flow = edgegrp.create_dataset('flow', data=flow)
    ds_flags = edgegrp.create_dataset('flags', data= mw_vessel_flag)
    ds_vessel_label = edgegrp.create_dataset('vessel_label', data= label)
    f.close()
        
        
        
def get_nodes_from_file(vesselgrp):
    filename = os.path.dirname(sys.argv[0]) + '/apj_network_546_segments_nodes.txt'
    f=open(filename,'r')
    data = f.readlines()
    cropped_nodes = data[1140:2112]


    node_index2label = []
    node_label2index = dict()
    #positions_of_nodes = np.zeros([ 972, 3])
    positions_of_nodes = []
    max_x=0
    max_y=0
    
    for line in cropped_nodes:
        k=0
        coordinate = []
        #print line
        for value in line.split():
            #print value
            if k==0:
                node_index2label.append(int(value))
                node_label2index[int(value)]=len(node_index2label)-1
            if k==1:
                coordinate.append(float(value))
                if(float(value)>max_x):
                    max_x=float(value)
            if k==2:
                coordinate.append(float(value))
                if(float(value)>max_y):
                    max_y=float(value)
            if k==3:
                #coordinate.append(int(value))
                #set to z=0 plane
                coordinate.append(float(0))
            k=k+1
        positions_of_nodes.append(coordinate)

    cropped_roots = data[2114:2150]
    indeces_of_roots=[]
    bctyp_of_roots=[]
    value_of_bc=[]
    for line in cropped_roots:
        k=0
        for value in line.split():
            if k==0:
                indeces_of_roots.append(node_label2index[int(value)])
            if k==1:
                bctyp_of_roots.append(int(value))
            if k==2:
                value_of_bc.append(float(value))
            k=k+1
            
    for (i,bctyp) in enumerate(bctyp_of_roots):
        if(bctyp == 2):#if it is a flow condition
            #from nl/min to mul/s
            #change signs
            value_of_bc[i] = -value_of_bc[i]/60.*1000000 /2
        if(bctyp == 0):#if it is a pressure condition
    #        #from mmhg to kpa
            if(value_of_bc[i]>0): 
                #michael works with possitve pressures only
                value_of_bc[i] = value_of_bc[i]*0.1333
            else:
                sys.exit("bad misstake to use negative pressure!")

    #node stuff
    N_nodes=972
    #N_nodes = 5585
    nodegrp = f3['vessels'].create_group("nodes")
    nodegrp.attrs.create('COUNT', N_nodes)
    ds_world_pos = nodegrp.create_dataset('world_pos', data = np.array(positions_of_nodes))
    ds_value_of_bc = nodegrp.create_dataset('bc_value', data=value_of_bc)        
    ds_bctyp_of_roots = nodegrp.create_dataset('bc_type', data= bctyp_of_roots)
    ds_bc_node_index = nodegrp.create_dataset('bc_node_index', data = indeces_of_roots)
    ds_roots = nodegrp.create_dataset('roots', data = indeces_of_roots)
    f.close()
    return node_label2index

def correct_vessel_indeces(node_label2index, vesselgrp):
    va = vesselgrp['edges/node_a_index']
    vb = vesselgrp['edges/node_b_index']
    del vesselgrp['edges/node_a_index']
    del vesselgrp['edges/node_b_index']
    va_new = []
    vb_new = []
    
    for old_index in va:
        va_new.append(node_label2index[old_index])
    for old_index in vb:
        vb_new.append(node_label2index[old_index])
        
    vesselgrp.create_dataset('edges/node_a_index', data= va_new)
    vesselgrp.create_dataset('edges/node_b_index', data= vb_new)
    
    print("done... renumbering")


if __name__ == '__main__':
    #write the data to a h5 file
    fn = 'apj.h5'
    f3 = h5files.open(fn,'w')
    vesselgrp = f3.create_group('vessels')
    vesselgrp.attrs.create('CLASS','REALWORLD')
    
    get_network_from_file(vesselgrp)
    node_label2index = get_nodes_from_file(vesselgrp)
    
    correct_vessel_indeces(node_label2index, vesselgrp)
    
    
    import krebsjobs.parameters.parameterSetsAdaption
    adaptionParams = getattr(krebsjobs.parameters.parameterSetsAdaption, 'apj_1')
    ##CALCULATE!!!!
    pressure, flow, force, hema = ku.calc_vessel_hydrodynamics_world(f3['vessels'], False, False, None, adaptionParams['calcflow'])
    
    f3.close


##node stuff
#N_nodes=len(positions_of_nodes)
#nodegrp = f3['vessels'].create_group("nodes")
#nodegrp.attrs.create('COUNT', N_nodes)
#ds_roots = f3['vessels/nodes'].create_dataset('roots', data = indeces_of_roots)
#
##dummy pressure for compatiblility
#ds_pressure = f3['vessels/nodes'].create_dataset('roots_pressure', data = roots_pressure)
#ds_value_of_bc = f3['vessels/nodes'].create_dataset('value_of_boundary_condition', data=value_of_bc)        
#ds_bctyp_of_roots = f3['vessels/nodes'].create_dataset('bctyp_of_roots', data= bctyp_of_roots)
#ds_world_pos = f3['vessels/nodes'].create_dataset('world_pos', data = np.array(positions_of_nodes))
#
#import krebs.adaption.parameterSets
#adaptionParams = getattr(krebs.adaption.parameterSets, 'default_rats')
##CALCULATE!!!!
#pressure, flow, force, hema = ku.calc_vessel_hydrodynamics_world(f3['vessels'], False, False, None, adaptionParams['calcflow'])
#
#diameter = np.multiply(radii,2)
#
#import matplotlib.pyplot as plot
#import itertools
#pressure_at_vessel=[]
#
#for NodeA, NodeB in itertools.izip(node_a_index,node_b_index):
#    pressure_at_vessel.append((pressure[NodeA]+pressure[NodeB])/2 * 1/0.1333)
#
#plot.loglog(pressure_at_vessel,diameter,'*')
#plot.xlabel("pressure")
#plot.ylabel("diameter")
##plot.plot(pressure_at_vessel)
#plot.grid()
#plot.show()
#
#print "ende"