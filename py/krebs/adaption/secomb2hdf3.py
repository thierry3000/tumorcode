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
import matplotlib.pyplot as plot
import itertools
import subprocess
import os

script_location = os.path.dirname(os.path.realpath(__file__))
f=open('%s/mesent_28_10_net.dat' % script_location,'r')

data = f.readlines()

cropped_edges = data[8:1138]

node_a_index=[]
node_b_index=[]
radii = []
flow = []
hema = []
schrott =[]

my_index_a = 0
my_index_b = 0
for line in cropped_edges:
    k=0
    for value in line.split():
        if k==2:
            node_a_index.append(int(value))
        if k==3:
            node_b_index.append(int(value))
        if k==4:
            radii.append(float(value)/2)
        if k==5:
            #flow.append(float(value)*60)
            flow.append(float(value)/60*1e6)
        if k==6:
            hema.append(float(value))
        k=k+1

    
negative_flow_indexes_secomb = np.where(np.array(flow)<0)[0]
print("%s have negative flows!" % negative_flow_indexes_secomb)

#making negative flows positive
#there are no negative flows except in at boundary
if 0:
  for ind in negative_flow_indexes_secomb:
      flow[ind] = math.fabs(flow[ind])
    
cropped_nodes = data[1140:2112]
positions_of_nodes = []
max_x=0
max_y=0

secomb_to_mw_index = dict()
mw_counter=0
for line in cropped_nodes:
    k=0
    coordinate = []
    #print line
    for value in line.split():
        #print value
        if k==0:
            secomb_to_mw_index[int(value)] = mw_counter
            mw_counter = mw_counter+1
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
            coordinate.append(0)
        k=k+1
    positions_of_nodes.append(coordinate)

#also create inverse map
mw2secomb = {v: k for k, v in secomb_to_mw_index.items()}
#now we have the map so that we can rename the nodes
node_a_index_new =[]
node_b_index_new =[]

for aindex in node_a_index:
    node_a_index_new.append(secomb_to_mw_index[aindex])
for aindex in node_b_index:
    node_b_index_new.append(secomb_to_mw_index[aindex])

cropped_roots = data[2114:2150]
indeces_of_roots=[]
bctyp_of_roots=[]
value_of_bc=[]
for line in cropped_roots:
    k=0
    for value in line.split():
        if k==0:
            indeces_of_roots.append(secomb_to_mw_index[int(value)])
        if k==1:
            bctyp_of_roots.append(int(value))
        if k==2:
            value_of_bc.append(float(value))
        k=k+1
f.close

#delete some boundary nodes ;-)
#will_delete_this_index = []
print mw2secomb[180]
#will_delete_this_index.append(23)
##will_delete_this_index.append(np.where(np.array(indeces_of_roots)==179)[0][0])
#for aindex in will_delete_this_index:
#    indeces_of_roots = np.delete(indeces_of_roots, aindex)
#    bctyp_of_roots = np.delete(bctyp_of_roots,aindex)
#    value_of_bc = np.delete(value_of_bc, aindex)

#adapt boundray
indeces_of_flow = np.where(np.array(bctyp_of_roots)==2)[0]#still secomb flow
indeces_of_pressure = np.where(np.array(bctyp_of_roots)==1)[0]
for aindex in indeces_of_flow:
    value_of_bc[aindex] = value_of_bc[aindex]/60 * 1000000
    #value_of_bc[aindex] = -value_of_bc[aindex]
    #value_of_bc[aindex] = np.fabs(value_of_bc[aindex])
    #bctyp_of_roots[aindex] = 1
for aindex in indeces_of_pressure:
    value_of_bc[aindex] = value_of_bc[aindex]*0.1333
    
#create artificial flags
#everything is VEIN
flags=np.multiply(np.ones(len(node_a_index), dtype = np.int32),ku.VEIN)
#looks better
for i in range(0,224):
    flags[i]=ku.ARTERY
for i in range(297,714):
    flags[i]=ku.CAPILLARY
for i in range(856,860):
    flags[i]=ku.CAPILLARY

#looks not so meaningfull
#for indA, indB in itertools.izip(node_a_index,node_b_index):
#        if indA<800:
#            flags[secomb_to_mw_index[indA]] = ku.ARTERY
#        if indB<800:
#            flags[secomb_to_mw_index[indA]] = ku.ARTERY
#        if indA>2000 and indA<5000:
#            flags[secomb_to_mw_index[indA]] = ku.CAPILLARY
#        if indB>2000 and indB<5000:
#            flags[secomb_to_mw_index[indA]] = ku.CAPILLARY

#print("secomb_to_mw_index[820]: %i" % secomb_to_mw_index[820])
#print("secomb_to_mw_index[825]: %i" % secomb_to_mw_index[825])
#print("secomb_to_mw_index[830]: %i" % secomb_to_mw_index[830])

print("mw2secomb[175]: %i" % mw2secomb[175])
print("mw2secomb[181]: %i" % mw2secomb[181])
print("mw2secomb[185]: %i" % mw2secomb[185])
print("mw2secomb[187]: %i" % mw2secomb[187])
print("mw2secomb[188]: %i" % mw2secomb[188])




            

print("secombs max extension in x and y: %i, %i " % (max_x,max_y))
print("found: %i node points @ secombs file "% len(positions_of_nodes))
print("found: %i edges @ secombs file "% len(node_a_index))
print("found: %i root points " % len(indeces_of_roots))
print("indeces_of_roots: %s " % indeces_of_roots)
print("bctyp: %s " % bctyp_of_roots)
print("with value mul/s: %s" % value_of_bc)


#write the data to a h5 file
fn = 'mesentry_secomb_546_v3.h5'
f3 = h5files.open(fn,'w')
g = f3.create_group('vessels')
g.attrs['CLASS'] = 'REALWORLD'

#edge stuff
N_edges=len(radii)
edgegrp = f3['vessels'].create_group("edges")
edgegrp.attrs.create('COUNT',N_edges)
ds_nodeA = edgegrp.create_dataset('node_a_index', data=node_a_index_new)
ds_nodeB = edgegrp.create_dataset('node_b_index', data= node_b_index_new)
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
adaptionParams = getattr(krebs.adaption.parameterSets, 'mesentry_paper')
#CALCULATE!!!!
pressure, flow, force, hema, again_flags = ku.calc_vessel_hydrodynamics(f3['vessels'], True, True, None, adaptionParams['calcflow'])
#pressure, flow, force, hema = ku.calc_vessel_hydrodynamics_world(f3['vessels'], False, False, None, adaptionParams['calcflow'])
#use all the nice things from mw
graph = ku.read_vessels_from_hdf(f3['vessels'], ['position','flags','pressure'], return_graph = True, return_not_found = False)
read_flags = graph.edges['flags']
ds_nodeflags = f3['vessels/nodes'].create_dataset('nodeflags',
data = ku.edge_to_node_property(len(mw2secomb),graph.edgelist,read_flags, 'or')
)
ds_pressure = f3['vessels/nodes'].create_dataset('pressure', data=graph.nodes['pressure'] )
#hopefully these flags are nicer
del edgegrp['flags']
ds_vesselflags2 = edgegrp.create_dataset('flags', data = again_flags)


diameter = np.multiply(radii,2)

if 0:
    
    pressure_at_vessel=[]
    
    for NodeA, NodeB in itertools.izip(node_a_index,node_b_index):
        pressure_at_vessel.append((pressure[NodeA]+pressure[NodeB])/2 * 1/0.1333)
    
    plot.loglog(pressure_at_vessel,diameter,'*')
    plot.xlabel("pressure")
    plot.ylabel("diameter")
    #plot.plot(pressure_at_vessel)
    plot.grid()
    plot.show()

if 0:
    h5files.closeall
    #subprocess.call(['rm', 'mesentry_secomb_546_adption_p_adaption_default.h5'])
    #subprocess.call(['python2', '/localdisk/thierry/spitzenforschung/branches/adaption/py/krebsjobs/submit-adaption.py' , 'adaption_default', 'mesentry_secomb_546.h5', '/vessels' ])
    subprocess.Popen(['python2 /localdisk/thierry/tumorcode/py/krebsjobs/submit-adaption.py mesentry_subset_vary /localdisk/thierry/output_adaption/secomb/mesentry_secomb_546.h5 /vessels'  ],stdout = subprocess.PIPE,stderr = subprocess.PIPE,shell=True)    
    subprocess.Popen(['python2 /localdisk/thierry/tumorcode/py/krebs/povrayRenderVessels.py mesentry_secomb_546_adption_p_mesentry_subset_vary.h5 /vessels/recomputed_flow'  ],stdout = subprocess.PIPE,stderr = subprocess.PIPE,shell=True)        
    subprocess.Popen(['python2 /localdisk/thierry/tumorcode/py/krebs/hdfvessels2vtk.py mesentry_secomb_546_adption_p_mesentry_subset_vary.h5 noconduc /vessels/recomputed_flow'  ],stdout = subprocess.PIPE,stderr = subprocess.PIPE,shell=True)
#python2 /localdisk/thierry/spitzenforschung/branches/adaption/py/krebs/hdfvessels2vtk.py mesentry_secomb_546_adption_p_adaption_default.h5 'noconduc' /vessels/recomputed_flow    
    #os.system("rm mesentry_secomb_546_adption_p_adaption_default.h5")
    #os.system('python2 /localdisk/thierry/spitzenforschung/branches/adaption/py/krebsjobs/submit-adaption.py 'adaption_default mesentry_secomb_546.h5 /vessels'')
    #os.system("python2 /localdisk/thierry/spitzenforschung/branches/adaption/py/krebs/povrayRenderVessels.py 'mesentry_secomb_546_adption_p_adaption_default.h5' '/vessels/recomputed_flow'")

print "ende"