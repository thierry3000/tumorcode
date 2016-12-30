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
Needs the file Network_New_version.dat
provided in the tumorcode repository.
measurement data provided by Secomb et. al.
'''
import krebsutils as ku
import h5files
import math
import numpy as np
import scipy.spatial.distance as sdist
import sys

def correct_bc_type_from_secomb_to_MW(anArray):
  for i,v in enumerate(anArray):
    if v == 0:
      anArray[i] = 1
  return anArray
#measurement file
f_more=open('Network_New_version.dat','r')

#input stream
data = f_more.readlines()
cropped_edges = data[8:453]
secomb_type=[]
node_a_index=[]
node_b_index=[]
radii = []
flow = []
hema = []
schrott =[]

uncirculated_vessels=[194,195,196,197,198,102]

for line in cropped_edges:
    k=0
    for value in line.split("\t"):
#        if (k==6 and float(value) ==0):
#            break
        if k==2:
            secomb_type.append(int(value))
        if k==3:
            node_a_index.append(int(value)-1)
        if k==4:
            node_b_index.append(int(value)-1)
        if k==5:
            radii.append(float(value)/2)
        if k==6:
            flow.append(float(value)*60./1000000)
            #flow.append(float(value)/60*1000000) 
        if k==7:
            hema.append(float(value))
        k=k+1


mw_vessel_flags=[]
for (i, aSecombType) in enumerate(secomb_type):
    #everything is suposed to be circulated
#    if( aSecombType == 10):#maybe it is better to ignore these ones?
#        node_a_index.pop(i+1)
#        node_b_index.pop(i+1)
#        radii.pop(i+1)
#        flow.pop(i+1)
#        hema.pop(i+1)
    flag_value = 0
    if( aSecombType != 10):
        flag_value = np.bitwise_or(flag_value, ku.CIRCULATED)
    #mw_vessel_flags.append(np.bitwise_or(ku.CIRCULATED,0))
    if(aSecombType == 5):
        #mw_vessel_flags.append(np.bitwise_or(ku.ARTERY,0))
        flag_value = np.bitwise_or(flag_value, ku.ARTERY)
    if(aSecombType == 4):
        #mw_vessel_flags.append(np.bitwise_or(ku.VEIN,0))
        flag_value = np.bitwise_or(flag_value, ku.VEIN)
#    if(aSecombType == 10):
#        mw_vessel_flags[i] = np.bitwise_xor(mw_vessel_flags[i],ku.CIRCULATED)
    mw_vessel_flags.append(flag_value)

''' defining capillaries 
according to: 
[1]A. R. Pries, T. W. Secomb, and P. Gaehtgens, “Structural adaptation and stability of microvascular networks: theory and simulations,” American Journal of Physiology - Heart and Circulatory Physiology, vol. 275, no. 2, pp. H349–H360, Aug. 1998.

there are 167 capillaries
due to the lack of knowledge we are guessing
'''
if 0:
    manually_adding=[59,60,101,102,103,
                     121, 122, 123, 380, 381, 382, 441, 442, 332, 331,
                     57, 58, 138, 139, 140, 391, 392, 393, 394, 395, 390, 389]
    for i in range(0,len(hema)):
        if( (hema[i] <= 0.1 and radii[i]<1) or (flow[i]<= 0.001563 and radii[i]<10) ): 
            #guess it is a capillary
            mw_vessel_flags[i] = np.bitwise_or(mw_vessel_flags[i], ku.CAPILLARY)
            #delte artery
            mw_vessel_flags[i] = np.bitwise_xor(mw_vessel_flags[i], ku.ARTERY)
    #    elif( i in manually_adding):
    #        mw_vessel_flags[i] = np.bitwise_xor(mw_vessel_flags[i], ku.CAPILLARY)
    #        mw_vessel_flags[i] = np.bitwise_xor(mw_vessel_flags[i], ku.ARTERY)
    

negative_flow_indexes_secomb = np.where(np.array(flow)<0)[0]
print("%s have negative flows!" % negative_flow_indexes_secomb)

#making negative flows positive
for ind in negative_flow_indexes_secomb:
#    mybuffer = node_a_index[ind]
#    node_a_index[ind] = node_b_index[ind]
#    node_b_index[ind] = mybuffer
    flow[ind] = math.fabs(flow[ind])
    
cropped_nodes = data[455:867]
positions_of_nodes = []
max_x=0
max_y=0

for line in cropped_nodes:
    k=0
    coordinate = []
    #print line
    for value in line.split("\t"):
        #print value
        if k==2:
            coordinate.append(float(value))
            if(float(value)>max_x):
                max_x=float(value)
        if k==3:
            coordinate.append(float(value))
            if(float(value)>max_y):
                max_y=float(value)
        if k==4:
            #coordinate.append(int(value))
            #set to z=0 plane
            coordinate.append(0)
        k=k+1
    positions_of_nodes.append(coordinate)

cropped_roots = data[869:874]
indeces_of_roots=[]
bctyp_of_roots=[]
value_of_bc=[]
for line in cropped_roots:
    k=0
    for value in line.split("\t"):
        if k==1:
            indeces_of_roots.append(int(value)-1)
        if k==2:
            bctyp_of_roots.append(int(value))
        if k==3:
            value_of_bc.append(float(value))
        k=k+1
f_more.close()

dummy_pressure=np.zeros(len(positions_of_nodes))
roots_pressure=np.zeros(len(positions_of_nodes))
for (i,bctyp) in enumerate(bctyp_of_roots):
    if(bctyp == 2):#if it is a flow condition
        #from nl/min to mul/s
        #value_of_bc[i] = math.fabs(value_of_bc[i])*60/1000000
        value_of_bc[i] = value_of_bc[i]*60./1000
    if(bctyp == 0):#if it is a pressure condition
#        #from mmhg to kpa
        if(value_of_bc[i]>0): 
            #michael works with possitve pressures only
            value_of_bc[i] = value_of_bc[i]*0.1333
            dummy_pressure[indeces_of_roots[i]] = value_of_bc[i]
            roots_pressure[indeces_of_roots[i]] = value_of_bc[i]
        else:
            sys.exit("bad misstake to use negative pressure!")

#calculate pressure at boundary nodes

#for vessel_index in indeces_of_boundary_vessels:
#    real_world_a = np.array(positions_of_nodes[node_a_index[vessel_index]])
#    real_world_b = np.array(positions_of_nodes[node_b_index[vessel_index]])
#    l = sdist.pdist([real_world_a.transpose(), real_world_b.transpose()])
#    print("distance: %f" % l)
#note 15 boundary artery, 23 boundary vein
#border_vessel_flags=[15,23,23,15,15,23,23,23,23]   
#press=[5, 2.5, 2, 3, 3.5, 1, 1.4,1.5,1.6]

#one negative flow in roots
#roots[3] = 102-1
print("secombs max extension in x and y: %i, %i " % (max_x,max_y))
print("found: %i node points @ secombs file "% len(positions_of_nodes))
print("found: %i edges @ secombs file "% len(node_a_index))
print("found: %i root points " % len(indeces_of_roots))
print("indeces_of_roots: %s " % indeces_of_roots)
print("bctyp: %s " % bctyp_of_roots)
print("with value mul/s: %s" % value_of_bc)

print("*********")
print("*** start creating vessel file")
print("*********")
fn = 'mesentry_secomb_test_more_data.h5'
f3 = h5files.open(fn,'w')
f3.create_group('vessels')
N_edges=len(radii)
N_nodes=len(positions_of_nodes)
edgegrp = f3['vessels'].create_group("edges")
edgegrp.attrs.create('COUNT',N_edges)

nodegrp = f3['vessels'].create_group("nodes")
nodegrp.attrs.create('COUNT', N_nodes)
ds_nodeA = edgegrp.create_dataset('node_a_index', data=node_a_index)
ds_nodeB = edgegrp.create_dataset('node_b_index', data= node_b_index)
ds_radius = edgegrp.create_dataset('radius', data=radii)
ds_hema = edgegrp.create_dataset('hematocrit', data=hema)
ds_flow = edgegrp.create_dataset('flow', data=flow)

dummy_vessel_flags = np.zeros(len(node_a_index))
#ds_vesselflags = edgegrp.create_dataset('flags', data = dummy_vessel_flags)
ds_vesselflags = edgegrp.create_dataset('flags', data = mw_vessel_flags)


ds_roots = f3['vessels/nodes'].create_dataset('roots', data = indeces_of_roots)
#dummy pressure
#ds_pressure = f3['vessels/nodes'].create_dataset('roots_pressure', data = roots_pressure)
ds_node_index = f3['vessels/nodes'].create_dataset('bc_node_index', data = indeces_of_roots)
ds_value_of_bc = f3['vessels/nodes'].create_dataset('bc_value', data=value_of_bc)        
ds_bctyp_of_roots = f3['vessels/nodes'].create_dataset('bc_type', data=correct_bc_type_from_secomb_to_MW(bctyp_of_roots))
ds_world_pos = f3['vessels/nodes'].create_dataset('world_pos', data = np.array(positions_of_nodes))
ds_cond = f3['vessels/nodes'].create_dataset('bc_conductivity_value', data=np.zeros_like(value_of_bc))



print("*********")
print("*** vessel file created ")
print("*********")
#f3.close()

import krebsjobs.parameters.parameterSetsAdaption
adaptionParams = getattr(krebsjobs.parameters.parameterSetsAdaption, 'secomb_mesentry')

f3['vessels'].attrs.create('CLASS','WORLD')
pressure, flow, force, hema = ku.calc_vessel_hydrodynamics(f3['vessels'], False, False, None, adaptionParams['calcflow'],storeCalculationInHDF=True)



#print(dd[1])
#fn = 'mesentry_secomb.h5'
#import krebs.adaption.parameterSets
#adaptionParams = getattr(krebs.adaption.parameterSets, 'default')
  #this creates a ordinary vessel group
#if 0:
#  f2 = h5files.open(fn,'w')
#  ku.vesselgen_generate_grid_no_flow(f2, sz, 2. , 'fcc', 3., 1., 5., adaptionParams['calcflow'] )
#else:
#  f2 = h5files.open(fn,'r+')  
#  print('we are loading...')
#edges, real_world_positions = ku.read_vessels_from_hdf(f2['vessels'],
#                                                       ('position'),return_graph=True)
#real_world_positions = ku.read_vessel_positions_from_hdf(f2['vessels'])



print("ende")