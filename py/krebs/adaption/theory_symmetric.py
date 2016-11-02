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

import numpy as np
import h5py
import itertools

def read_vessel_lengths(vessel_grp):
    world_pos= np.asarray(vessel_grp['nodes/world_pos'])
    node_a_index = np.asarray(vessel_grp['edges/node_a_index'])
    node_b_index = np.asarray(vessel_grp['edges/node_b_index'])
    world_vecs = []
    for (a,b) in itertools.izip(node_a_index,node_b_index):
        world_vecs.append(world_pos[a,:]-world_pos[b,:])
    
    l_ges = np.linalg.norm(world_pos[33,:]-world_pos[32,:])
    return np.linalg.norm(world_vecs, axis=1), l_ges


#
#p_1 = 8.347
#
##first level
#p_2 = p_1 -(((Q_ges/2)*(8*etha*560))/np.pi*23.*np.power(2,-1/.2))
#
#print(p_2)

if __name__ == '__main__':
    filename1= '/localdisk/thierry/output_adaption/own_networks/test_configs.h5'
    f = h5py.File(filename1,'r')
    vessel_grp = f['/symetric/vessels_after_adaption/vessels_after_adaption/']
    
    vessel_lengths, l_ges = read_vessel_lengths(vessel_grp)
    print(vessel_lengths)
    
    one_total_sequenz = [37, 8, 7, 5, 4, 2, 1, 0, 36]
    delta_p_for_sequenz = []
    total_circulated_length = 0
    for number in one_total_sequenz:
        total_circulated_length = total_circulated_length + vessel_lengths[number]
    for number in one_total_sequenz:
        delta_p_for_sequenz.append(8.*0.3*vessel_lengths[number]/total_circulated_length)
    print total_circulated_length
    print l_ges
    dp_ges = 8 * 0.3 #kPa
    r_0 = 23./2.
    etha= 8.*1.2e-6
    Q_ges = dp_ges*np.pi*r_0**4/(8.*etha*l_ges)

    #zero level (just du make sure)
    under_root= (Q_ges*8*etha*vessel_lengths[one_total_sequenz[0]])/(delta_p_for_sequenz[0]*np.pi)
    r0 = np.power(under_root,1/4.)
    r00 = vessel_lengths[one_total_sequenz[0]]/50./delta_p_for_sequenz[0]
    print(r0)
    print("From const. shear: %f" % np.round(r00,1))
    
    #first level
    under_root= (Q_ges/2.*8*etha*vessel_lengths[one_total_sequenz[1]])/(delta_p_for_sequenz[1]*np.pi)
    r1 = np.power(under_root,1/4.)
    r11 = vessel_lengths[one_total_sequenz[1]]/50./delta_p_for_sequenz[1]
    print(np.round(r1,1))
    print("From const. shear: %f" % np.round(r11,1))
    
    #second level
    under_root= (Q_ges/4.*8*etha*vessel_lengths[one_total_sequenz[2]])/(delta_p_for_sequenz[2]*np.pi)
    r2 = np.power(under_root,1/4.)
    r22 = vessel_lengths[one_total_sequenz[2]]/50./delta_p_for_sequenz[2]
    print(np.round(r2,1))
    print("From const. shear: %f" % np.round(r22,1))
    
    #third level
    under_root= (Q_ges/(2**3)*8*etha*vessel_lengths[one_total_sequenz[3]])/(delta_p_for_sequenz[3]*np.pi)
    r3 = np.power(under_root,1/4.)
    r33 = vessel_lengths[one_total_sequenz[3]]/50./delta_p_for_sequenz[3]
    print(np.round(r3,1))
    print("From const. shear: %f" % np.round(r33,1))
    
    #fourth level
    under_root= (Q_ges/(2**4)*8*etha*vessel_lengths[one_total_sequenz[4]])/(delta_p_for_sequenz[4]*np.pi)
    r4 = np.power(under_root,1/4.)
    r44 = vessel_lengths[one_total_sequenz[4]]/50./delta_p_for_sequenz[4]
    print(np.round(r4,1))
    print("From const. shear: %f" % np.round(r44,1))
    