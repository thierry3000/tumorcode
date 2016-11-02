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

f=open('/localdisk/thierry/output_adaption/secomb/Network_secomb_mesentry.dat','r')

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
            bctyp_of_roots.append(int(value))
        if k==2:
            value_of_bc.append(float(value))
        k=k+1
f.close

#create artificial flags

flags=np.multiply(np.ones(len(node_a_index)),ku.VEIN)
for i in range(0,42):
    flags[i]=ku.ARTERY

flags[65] = ku.ARTERY
flags[66] = ku.ARTERY
flags[67] = ku.ARTERY
flags[68] = ku.ARTERY
flags[69] = ku.ARTERY
flags[70] = ku.ARTERY
flags[143] = ku.ARTERY
flags[144] = ku.ARTERY
flags[145] = ku.ARTERY
flags[146] = ku.ARTERY
flags[147] = ku.ARTERY

flags[89] = ku.ARTERY
flags[90] = ku.ARTERY
flags[91] = ku.ARTERY
flags[92] = ku.ARTERY
flags[93] = ku.ARTERY

flags[54] = ku.ARTERY
flags[55] = ku.ARTERY
flags[56] = ku.ARTERY
flags[57] = ku.ARTERY
flags[58] = ku.ARTERY
flags[138] = ku.ARTERY
flags[139] = ku.ARTERY
flags[140] = ku.ARTERY

flags[123] = ku.ARTERY
flags[122] = ku.ARTERY
flags[121] = ku.ARTERY

flags[52] = ku.ARTERY
flags[53] = ku.ARTERY

flags[133] = ku.CAPILLARY
flags[134] = ku.CAPILLARY
flags[135] = ku.CAPILLARY
flags[136] = ku.CAPILLARY
flags[137] = ku.CAPILLARY

flags[104] = ku.ARTERY
flags[105] = ku.ARTERY
flags[106] = ku.ARTERY
flags[107] = ku.ARTERY
flags[108] = ku.ARTERY

flags[112] = ku.ARTERY
flags[113] = ku.ARTERY
flags[114] = ku.ARTERY
flags[115] = ku.ARTERY
flags[116] = ku.ARTERY
flags[117] = ku.ARTERY

flags[76] = ku.ARTERY
flags[77] = ku.ARTERY
flags[78] = ku.ARTERY
flags[79] = ku.ARTERY
flags[80] = ku.ARTERY
flags[81] = ku.ARTERY
flags[82] = ku.ARTERY
flags[83] = ku.ARTERY
flags[84] = ku.ARTERY
flags[85] = ku.ARTERY

flags[169] = ku.ARTERY
flags[170] = ku.ARTERY
flags[169] = ku.ARTERY
flags[171] = ku.ARTERY
flags[172] = ku.ARTERY
flags[173] = ku.ARTERY
flags[174] = ku.ARTERY
flags[175] = ku.ARTERY

flags[159] = ku.ARTERY
flags[160] = ku.ARTERY
flags[161] = ku.ARTERY
flags[162] = ku.ARTERY
flags[163] = ku.ARTERY
flags[164] = ku.ARTERY

flags[94] = ku.ARTERY
flags[95] = ku.ARTERY
flags[96] = ku.ARTERY
flags[97] = ku.ARTERY
flags[98] = ku.ARTERY
flags[99] = ku.ARTERY
flags[100] = ku.ARTERY

flags[207] = ku.CAPILLARY
flags[255] = ku.CAPILLARY
flags[256] = ku.CAPILLARY
flags[257] = ku.CAPILLARY
flags[258] = ku.CAPILLARY
flags[367] = ku.CAPILLARY
flags[368] = ku.CAPILLARY

#flags[340] = ku.ARTERY
#flags[341] = ku.ARTERY
#flags[342] = ku.ARTERY
#flags[343] = ku.ARTERY
#flags[344] = ku.ARTERY
#flags[337] = ku.ARTERY
#flags[338] = ku.ARTERY

#flags[353] = ku.ARTERY
#flags[354] = ku.ARTERY
#flags[355] = ku.ARTERY
#flags[358] = ku.ARTERY
#flags[357] = ku.ARTERY
#flags[356] = ku.ARTERY


flags[136] = ku.CAPILLARY

#flags[339] = ku.CAPILLARY






flags[248] = ku.CAPILLARY
flags[249] = ku.CAPILLARY
flags[250] = ku.CAPILLARY


flags[148] = ku.CAPILLARY
flags[235] = ku.CAPILLARY
flags[236] = ku.CAPILLARY
flags[237] = ku.CAPILLARY
flags[71] = ku.CAPILLARY
flags[72] = ku.CAPILLARY
flags[73] = ku.CAPILLARY
flags[149] = ku.CAPILLARY
    
flags[188]=ku.CAPILLARY
flags[187]=ku.CAPILLARY
flags[186]=ku.CAPILLARY
flags[185]=ku.CAPILLARY
flags[184]=ku.CAPILLARY

flags[194]=ku.CAPILLARY
flags[195]=ku.CAPILLARY
flags[196]=ku.CAPILLARY
flags[197]=ku.CAPILLARY
flags[198]=ku.CAPILLARY

flags[193]=ku.CAPILLARY
flags[192]=ku.CAPILLARY
flags[191]=ku.CAPILLARY
flags[190]=ku.CAPILLARY
flags[189]=ku.CAPILLARY

flags[247]=ku.CAPILLARY
flags[246]=ku.CAPILLARY
flags[245]=ku.CAPILLARY
flags[244]=ku.CAPILLARY

flags[64]=ku.CAPILLARY
flags[63]=ku.CAPILLARY
flags[62]=ku.CAPILLARY
flags[61]=ku.CAPILLARY

flags[226]=ku.CAPILLARY
flags[227]=ku.CAPILLARY
flags[228]=ku.CAPILLARY
flags[229]=ku.CAPILLARY
flags[230]=ku.CAPILLARY
flags[231]=ku.CAPILLARY

flags[336]=ku.CAPILLARY
flags[335]=ku.CAPILLARY
flags[338]=ku.CAPILLARY
flags[337]=ku.CAPILLARY
flags[132]=ku.CAPILLARY
flags[131]=ku.CAPILLARY
flags[51]=ku.CAPILLARY
flags[50]=ku.CAPILLARY
flags[223]=ku.CAPILLARY
flags[224]=ku.CAPILLARY
flags[225]=ku.CAPILLARY
flags[130]=ku.CAPILLARY
flags[129]=ku.CAPILLARY
flags[222]=ku.CAPILLARY
flags[221]=ku.CAPILLARY
flags[49]=ku.CAPILLARY
flags[48]=ku.CAPILLARY
flags[47]=ku.CAPILLARY
flags[46]=ku.CAPILLARY
flags[45]=ku.CAPILLARY
flags[44]=ku.CAPILLARY
flags[42]=ku.CAPILLARY
flags[43]=ku.CAPILLARY
flags[211]=ku.CAPILLARY
flags[212]=ku.CAPILLARY
flags[213]=ku.CAPILLARY
flags[214]=ku.CAPILLARY
flags[215]=ku.CAPILLARY
flags[216]=ku.CAPILLARY
flags[217]=ku.CAPILLARY
flags[218]=ku.CAPILLARY
flags[219]=ku.CAPILLARY
flags[220]=ku.CAPILLARY
flags[127]=ku.CAPILLARY
flags[128]=ku.CAPILLARY
flags[124]=ku.CAPILLARY
flags[125]=ku.CAPILLARY
flags[126]=ku.CAPILLARY
flags[320]=ku.CAPILLARY
flags[321]=ku.CAPILLARY
flags[322]=ku.CAPILLARY

flags[141]=ku.CAPILLARY
flags[142]=ku.CAPILLARY
flags[391]=ku.CAPILLARY
flags[392]=ku.CAPILLARY
flags[393]=ku.CAPILLARY
flags[394]=ku.CAPILLARY
flags[395]=ku.CAPILLARY
flags[389]=ku.CAPILLARY
flags[390]=ku.CAPILLARY
flags[260]=ku.CAPILLARY
flags[259]=ku.CAPILLARY
flags[208]=ku.CAPILLARY
flags[209]=ku.CAPILLARY
flags[210]=ku.CAPILLARY
flags[103]=ku.CAPILLARY
flags[102]=ku.CAPILLARY
flags[101]=ku.CAPILLARY
flags[59]=ku.CAPILLARY
flags[60]=ku.CAPILLARY
flags[232]=ku.CAPILLARY
flags[233]=ku.CAPILLARY
flags[234]=ku.CAPILLARY
flags[380]=ku.CAPILLARY
flags[381]=ku.CAPILLARY
flags[382]=ku.CAPILLARY
flags[441]=ku.CAPILLARY
flags[442]=ku.CAPILLARY


flags[421]=ku.CAPILLARY
flags[420]=ku.CAPILLARY
flags[111]=ku.CAPILLARY
flags[110]=ku.CAPILLARY
flags[109]=ku.CAPILLARY
flags[120]=ku.CAPILLARY
flags[119]=ku.CAPILLARY
flags[118]=ku.CAPILLARY

flags[168]=ku.CAPILLARY
flags[167]=ku.CAPILLARY
flags[166]=ku.CAPILLARY
flags[165]=ku.CAPILLARY

flags[179]=ku.CAPILLARY
flags[178]=ku.CAPILLARY
flags[177]=ku.CAPILLARY
flags[176]=ku.CAPILLARY

flags[437]=ku.CAPILLARY
flags[436]=ku.CAPILLARY
flags[183]=ku.CAPILLARY
flags[182]=ku.CAPILLARY
flags[181]=ku.CAPILLARY
flags[180]=ku.CAPILLARY
flags[243]=ku.CAPILLARY
flags[242]=ku.CAPILLARY
flags[241]=ku.CAPILLARY
flags[240]=ku.CAPILLARY
flags[239]=ku.CAPILLARY
flags[238]=ku.CAPILLARY
flags[88]=ku.CAPILLARY
flags[87]=ku.CAPILLARY
flags[86]=ku.CAPILLARY


flags[203]=ku.CAPILLARY
flags[204]=ku.CAPILLARY
flags[205]=ku.CAPILLARY
flags[206]=ku.CAPILLARY

flags[251]=ku.CAPILLARY
flags[252]=ku.CAPILLARY
flags[253]=ku.CAPILLARY
flags[254]=ku.CAPILLARY

#flags[361]=ku.CAPILLARY
#flags[362]=ku.CAPILLARY
#flags[363]=ku.CAPILLARY
#flags[364]=ku.CAPILLARY
#flags[365]=ku.CAPILLARY
#flags[366]=ku.CAPILLARY
#flags[371]=ku.CAPILLARY
#flags[370]=ku.CAPILLARY
#flags[369]=ku.CAPILLARY
#flags[375]=ku.CAPILLARY
#flags[374]=ku.CAPILLARY
#flags[373]=ku.CAPILLARY
#flags[372]=ku.CAPILLARY

flags[199]=ku.CAPILLARY
flags[200]=ku.CAPILLARY
flags[201]=ku.CAPILLARY
flags[202]=ku.CAPILLARY

#flags[379]=ku.CAPILLARY
#flags[378]=ku.CAPILLARY
#flags[377]=ku.CAPILLARY
#flags[376]=ku.CAPILLARY



#needed for compatiblility with old code
roots_pressure=np.zeros(len(positions_of_nodes))
for (i,bctyp) in enumerate(bctyp_of_roots):
    if(bctyp == 2):#if it is a flow condition
        #from nl/min to mul/s
        #value_of_bc[i] = np.fabs(value_of_bc[i]*1000000/60)
        value_of_bc[i] = value_of_bc[i]*1000000/60
    if(bctyp == 1):#if it is a pressure condition
#        #from mmhg to kpa
        if(value_of_bc[i]>0): 
            #michael works with possitve pressures only
            value_of_bc[i] = value_of_bc[i]*0.1333
        else:
            sys.exit("bad misstake to use negative pressure!")
            
#mw_vessel_flags=[]
#for (i, aflow) in enumerate(flow):
#    #everything is suposed to be circulated
#    flag_value = 0
#    flag_value = np.bitwise_or(flag_value, ku.CIRCULATED)
#    mw_vessel_flags.append(flag_value)


print("secombs max extension in x and y: %i, %i " % (max_x,max_y))
print("found: %i node points @ secombs file "% len(positions_of_nodes))
print("found: %i edges @ secombs file "% len(node_a_index))
print("found: %i root points " % len(indeces_of_roots))
print("indeces_of_roots: %s " % indeces_of_roots)
print("bctyp: %s " % bctyp_of_roots)
print("with value mul/s: %s" % value_of_bc)


#write the data to a h5 file
fn = 'mesentry_secomb.h5'
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

import krebsjobs.parameters.parameterSetsAdaption
adaptionParams = getattr(krebsjobs.parameters.parameterSetsAdaption, 'adaption_default')
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