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

''' 
    This program uses measured data on rat mesentry network
    available at 
    http://ajpheart.physiology.org/highwire/filestream/95027/field_highwire_adjunct_files/1/network_546_segments.txt
    
    For convenience a copy is provide however the credits
    for measurement go to the authors secomb et al.
    
    1. we read the data from txt file an generate the file apj.h5
       which is in 'tumorcode' format.
    2. the 'ku.calc_vessel_hydrodynamics' command calculated
       all properties relevant for further proceding with 
       that data in tumorcode e.g. vizualization, 
'''
if __name__ == '__main__':
  import os.path, sys
  sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../..'))

import krebsutils as ku
import math
import numpy as np
import scipy.spatial.distance as sdist
import sys
import os

def correct_bc_type_from_secomb_to_MW(anArray):
  for i,v in enumerate(anArray):
    if v == 0:
      anArray[i] = 1
  return anArray

def get_network_from_file(vesselgrp, path):
    #filename = os.path.dirname(sys.argv[0]) + '/apj_network_546_segments.txt'
    filename = path + '/apj_network_546_segments.txt'    
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
            #flow.append((np.pi*radii[i]**2)*length[i]*velocity[i]*1000)
            #flow.append((np.pi*2*radii[i])*length[i]*velocity[i]*1000)
            flow.append((np.pi*radii[i]**2)*velocity[i]*1000)
    #edge stuff
    N_edges= 546
    ''' there are double entries, which maybe spoil the computation '''
#    matrix = np.asarray([node_a_index,node_b_index])
#    #unique_matrix = np.unique(matrix, axis=1, return_counts=True)
#    mt = matrix.transpose()
#    how_often_present=(mt[:, np.newaxis] == mt).all(axis=2).sum(axis=1)
#    repeated_edges=how_often_present>1;
#    single_edges = np.logical_not(repeated_edges);
#    N_edges=np.sum(single_edges)

#    list_to_change = ['node_a_index','node_b_index','radii','hema', 'flow', 'mw_vessel_flag','label']
#    for element in list_to_change:
#      exec("%s=np.asarray(%s)[single_edges]" % (element,element))
    edgegrp = vesselgrp.create_group("edges")
    edgegrp.attrs.create('COUNT',N_edges)
    ds_nodeA = edgegrp.create_dataset('node_a_index', data= node_a_index)
    ds_nodeB = edgegrp.create_dataset('node_b_index', data= node_b_index)
    ds_radius = edgegrp.create_dataset('radius', data=radii)
    ds_hema = edgegrp.create_dataset('hematocrit', data=hema)
    ds_flow = edgegrp.create_dataset('flow', data=flow)
    ds_flags = edgegrp.create_dataset('flags', data= mw_vessel_flag)
    ds_vessel_label = edgegrp.create_dataset('vessel_label', data= label)
    f.close()
        
        
        
def get_nodes_from_file(vesselgrp, path):
    #filename = os.path.dirname(sys.argv[0]) + '/apj_network_546_segments_nodes.txt'
    filename = path + '/apj_network_546_segments_nodes.txt'    
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
            value_of_bc[i] = -value_of_bc[i]/60.*1000000
        if(bctyp == 1):#if it is a pressure condition
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
    ds_bctyp_of_roots = nodegrp.create_dataset('bc_type', data= correct_bc_type_from_secomb_to_MW(bctyp_of_roots))
    ds_bc_node_index = nodegrp.create_dataset('bc_node_index', data = indeces_of_roots)
    ds_bc_conductivity_value = nodegrp.create_dataset('bc_conductivity_value', data=np.zeros_like(value_of_bc)) 
    ds_roots = nodegrp.create_dataset('roots', data = indeces_of_roots)
    ds_nodeflags = nodegrp.create_dataset('nodeflags', data = np.zeros(N_nodes,dtype='int32'))
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

''' I found that there are some edges twice, 
    delete them
    '''
def delete_double_edges(vesselgrp):
  #176 deleted for convergenz issues
  by_hand_identified_doubles = [279,283,284,294,312,325,328]
  va = np.asarray(vesselgrp['edges/node_a_index'])
  vb = np.asarray(vesselgrp['edges/node_b_index'])
  flags = np.asarray(vesselgrp['edges/flags'])
  radii = np.asarray(vesselgrp['edges/radius'])
  hema = np.asarray(vesselgrp['edges/hematocrit'])
  no_edges = len(va)
  good_indeces = np.ones(no_edges)
  for i in by_hand_identified_doubles:
    good_indeces[i] = 0
  good_indeces = good_indeces.astype(bool)
  del vesselgrp['edges/node_a_index']
  del vesselgrp['edges/node_b_index']
  del vesselgrp['edges/flags']
  del vesselgrp['edges/radius']
  del vesselgrp['edges/hematocrit']
#  matrix = np.asarray([va,vb])
#  unique_matrix = np.unique(matrix, axis=1)
#  va_new = unique_matrix[0,:]
#  vb_new = unique_matrix[1,:]
  va_new = va[good_indeces]
  vb_new = vb[good_indeces]
  vesselgrp.create_dataset('edges/node_a_index', data= va_new)
  vesselgrp.create_dataset('edges/node_b_index', data= vb_new)
  N_edges= len(va_new)
  edgegrp = vesselgrp['edges']
  edgegrp.attrs['COUNT']= N_edges
  flags_new = flags[good_indeces]
  vesselgrp.create_dataset('edges/flags', data= flags_new)
  radii_new = radii[good_indeces]
  vesselgrp.create_dataset('edges/radius', data= radii_new)
  hema_new = hema[good_indeces]
  vesselgrp.create_dataset('edges/hematocrit', data= hema_new)
  
  print("done... deleting double entries")

if __name__ == '__main__':
    #write the data to a h5 file
    fn = 'apj.h5'
    f3 = h5files.open(fn,'w')
    vesselgrp = f3.create_group('vessels')
    vesselgrp.attrs.create('CLASS','REALWORLD')
    
    '''***** INPUT *****'''
    '''path to apj_network_546_segments.txt
    this might change during your installation'''
    path_to_txt = '/home/usersHR/thierry/tumorcode/py/krebs/adaption'
    get_network_from_file(vesselgrp, path_to_txt)
    node_label2index = get_nodes_from_file(vesselgrp, path_to_txt)
    delete_double_edges(vesselgrp)
    correct_vessel_indeces(node_label2index, vesselgrp)
    
    '''***** tumorcode stuff *****'''
    import krebsjobs.parameters.parameterSetsAdaption
    adaptionParams = getattr(krebsjobs.parameters.parameterSetsAdaption, 'apj')
    ##CALCULATE!!!!
    #pressure, flow, force, hema, flags = ku.calc_vessel_hydrodynamics(f3['vessels'], False, False, None, adaptionParams['calcflow'],storeCalculationInHDF=True)
    dd = ku.calc_vessel_hydrodynamics(f3['vessels'], return_flags = True, bloodflowparams = adaptionParams['calcflow'],storeCalculationInHDF=True)
    f3.close

