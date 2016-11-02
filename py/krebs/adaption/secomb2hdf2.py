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
import krebsutils as ku
import h5files
import math
import numpy as np
import scipy.spatial.distance as sdist
import sys

f=open('Network_secomb_mesentry.dat','r')

data = f.readlines()

cropped_edges = data[8:453]

node_a_index=[]
node_b_index=[]
radii = []
flow = []
hema = []
schrott =[]

for line in cropped_edges:
    k=0
    for value in line.split("\t"):
        if k==2:
            node_a_index.append(int(value)-1)
        if k==3:
            node_b_index.append(int(value)-1)
        if k==4:
            radii.append(float(value)/2)
        if k==5:
            #flow.append(float(value)*60)
            flow.append(float(value)/60*1000000)
        if k==6:
            hema.append(float(value))
        k=k+1

    
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
        if k==1:
            coordinate.append(int(value))
            if(int(value)>max_x):
                max_x=int(value)
        if k==2:
            coordinate.append(int(value))
            if(int(value)>max_y):
                max_y=int(value)
        if k==3:
            #coordinate.append(int(value))
            #set to z=0 plane
            coordinate.append(0)
        k=k+1
    positions_of_nodes.append(coordinate)

cropped_roots = data[869:879]
indeces_of_roots=[]
bctyp_of_roots=[]
value_of_bc=[]
for line in cropped_roots:
    k=0
    for value in line.split("\t"):
        if k==0:
            indeces_of_roots.append(int(value)-1)
        if k==1:
            if int(value)==2:#if secomb flow make mw flow
                bctyp_of_roots.append(int(value)-1)
            else:
                bctyp_of_roots.append(int(value))
        if k==2:
            value_of_bc.append(float(value))
        k=k+1
f.close

#create artificial flags
#everything is VEIN
flags=np.multiply(np.ones(len(node_a_index)),ku.VEIN)

for i in range(0,44):
    flags[i]=ku.ARTERY
for i in range(44,50):
    flags[i]=ku.ARTERY
    
flags[50] = ku.ARTERY
flags[51] = ku.ARTERY
flags[52] = ku.ARTERY
flags[53] = ku.ARTERY
for i in range(54,61):
    flags[i]=ku.ARTERY
for i in range(61,76):
    flags[i]=ku.ARTERY
for i in range(76,89):
    flags[i]=ku.ARTERY
for i in range(89,94):
    flags[i]=ku.ARTERY
for i in range(94,101):
    flags[i]=ku.ARTERY
for i in range(101,104):
    flags[i]=ku.ARTERY
    

for i in range(104,113):
    flags[i]=ku.CAPILLARY
for i in range(133,138):
    flags[i]=ku.CAPILLARY

for i in range(169,180):
    flags[i]=ku.CAPILLARY
for i in range(159,169):
    flags[i]=ku.CAPILLARY
for i in range(184,194):
    flags[i]=ku.CAPILLARY
for i in range(226,232):
    flags[i]=ku.CAPILLARY
for i in range(244,248):
    flags[i]=ku.CAPILLARY
flags[248] = ku.CAPILLARY
flags[249] = ku.CAPILLARY
flags[250] = ku.CAPILLARY

for i in range(104,261):
    flags[i]=ku.CAPILLARY




            
#read pressures manually 0,94,193,207,256,308,330,340,407
if 0:
    for i in range(len(bctyp_of_roots)):
        bctyp_of_roots[i]=1#all pressure
    value_of_bc[0] = 7.79685
    value_of_bc[1] = 4.19766
    value_of_bc[2] = 3.7617
    value_of_bc[3] = 3.60944
    value_of_bc[4] = 1.9995
    value_of_bc[5] = 2.11303
    value_of_bc[6] = 2.681
    value_of_bc[7] = 3.32947
    value_of_bc[8] = 2.51393
    


print("secombs max extension in x and y: %i, %i " % (max_x,max_y))
print("found: %i node points @ secombs file "% len(positions_of_nodes))
print("found: %i edges @ secombs file "% len(node_a_index))
print("found: %i root points " % len(indeces_of_roots))
print("indeces_of_roots: %s " % indeces_of_roots)
print("bctyp: %s " % bctyp_of_roots)
print("with value mul/s: %s" % value_of_bc)


#write the data to a h5 file
fn = '/localdisk/thierry/output_adaption/secomb/mesentry_secomb.h5'
f3 = h5files.open(fn,'w')
f3.create_group('vessels')

#edge stuff
N_edges=len(radii)
edgegrp = f3['vessels'].create_group("edges")
edgegrp.attrs.create('COUNT',N_edges)
ds_nodeA = edgegrp.create_dataset('node_a_index', data=node_a_index)
ds_nodeB = edgegrp.create_dataset('node_b_index', data= node_b_index)
ds_radius = edgegrp.create_dataset('radius', data=radii)
ds_hema = edgegrp.create_dataset('hematocrit', data=hema)
ds_flow = edgegrp.create_dataset('flow', data=flow)

#mw_vessel_flags = get_flags_from_other_file()

#ds_vesselflags = edgegrp.create_dataset('flags', data = get_flags_from_other_file())

#node stuff
N_nodes=len(positions_of_nodes)
nodegrp = f3['vessels'].create_group("nodes")
nodegrp.attrs.create('COUNT', N_nodes)
ds_roots = f3['vessels/nodes'].create_dataset('roots', data = indeces_of_roots)

ds_bc_value = f3['vessels/nodes'].create_dataset('bc_value', data=value_of_bc)        
ds_bc_type = f3['vessels/nodes'].create_dataset('bc_type', data= bctyp_of_roots)
ds_bc_node_index = f3['vessels/nodes'].create_dataset('bc_node_index', data= indeces_of_roots)
ds_world_pos = f3['vessels/nodes'].create_dataset('world_pos', data = np.array(positions_of_nodes))

ds_vesselflags = edgegrp.create_dataset('flags', data = flags)

import krebs.adaption.parameterSets
adaptionParams = getattr(krebs.adaption.parameterSets, 'adaption_default')
#CALCULATE!!!!
dd = ku.calc_vessel_hydrodynamics_world(f3['vessels'], False, False, None, adaptionParams['calcflow'])
#pressure, flow, force, hema = ku.calc_vessel_hydrodynamics_world(f3['vessels'], False, False, None, adaptionParams['calcflow'])

diameter = np.multiply(radii,2)

if 0:
    import matplotlib.pyplot as plot
    import itertools
    pressure_at_vessel=[]
    
    for NodeA, NodeB in itertools.izip(node_a_index,node_b_index):
        pressure_at_vessel.append((pressure[NodeA]+pressure[NodeB])/2 * 1/0.1333)
    
    plot.loglog(pressure_at_vessel,diameter,'*')
    plot.xlabel("pressure")
    plot.ylabel("diameter")
    #plot.plot(pressure_at_vessel)
    plot.grid()
    plot.show()

print "ende"