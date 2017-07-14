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
This script is used to test and developt the structural aption in microvessels
/*
 * Adaption According to  Pries, Secomb, Gaehtgens 1998
 * "Structural adaption and stability of microvascular networks: theory and simulations"
 */
"""
if __name__ == '__main__':
  import os.path, sys
  sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'..'))
  
import numpy as np
import krebsutils as ku
import math
import collections
import itertools
import h5files
import myutils
import os
import sys

import h5py


def createMostBasic_world(fn):    
    
    points=[]
    points.append([-5,0,0])
    p1=[0,0,0]
    p2=[2.5,3,0]
    p3=[7.5,3,0]
    
    p4=[10,0,0]
    
    p5=[2.5,-3,0]
    p6=[7.5,-3,0]
    p7=[15,0,0]
    
    points.append(p1)
    points.append(p2)
    points.append(p3)
    points.append(p4)
    points.append(p5)
    points.append(p6)
    points.append(p7)
    
    
    node_a_index=[0,1,2,3,1,5,6,4]
    node_b_index=[1,2,3,4,5,6,4,7]
    flags=[]
    radii=[]

    flags.append(int(ku.ARTERY))#0
    radii.append(7)
    flags.append(int(ku.ARTERY))#1
    radii.append(5)
    flags.append(int(ku.CAPILLARY))#2
    radii.append(3)
    flags.append(int(ku.VEIN))#3
    radii.append(6)
    flags.append(int(ku.ARTERY))#4
    radii.append(5)
    flags.append(int(ku.CAPILLARY))#5
    radii.append(3)
    flags.append(int(ku.VEIN))#6
    radii.append(6)
    flags.append(int(ku.VEIN))#7
    radii.append(10)

    f2 = fn
    f2.create_group('mostBasic_world/vessels')
    
    #edge stuff
    N_edges=len(radii)
    edgegrp = f2['mostBasic_world/vessels'].create_group("edges")
    edgegrp.attrs.create('COUNT',N_edges)
    ds_nodeA = edgegrp.create_dataset('node_a_index', data=node_a_index)
    ds_nodeB = edgegrp.create_dataset('node_b_index', data= node_b_index)
    ds_radius = edgegrp.create_dataset('radius', data=radii)
    #ds_hema = edgegrp.create_dataset('hematocrit', data=hema)
    #ds_flow = edgegrp.create_dataset('flow', data=flow)
    
    ds_vesselflags = edgegrp.create_dataset('flags', data = flags)
    
    #node stuff
    N_nodes=len(points)
    nodegrp = f2['mostBasic_world/vessels'].create_group("nodes")
    nodegrp.attrs.create('COUNT', N_nodes)
    ds_roots = f2['mostBasic_world/vessels/nodes'].create_dataset('roots', data = [0,7])
    ds_bc_node_index = f2['mostBasic_world/vessels/nodes'].create_dataset('bc_node_index', data = [0,7])
    
       
    ds_value_of_bc = f2['mostBasic_world/vessels/nodes'].create_dataset('bc_value', data=[5, 3])        
    ds_bctyp_of_roots = f2['mostBasic_world/vessels/nodes'].create_dataset('bc_type', data= [1,1])
    ds_world_pos = f2['mostBasic_world/vessels/nodes'].create_dataset('world_pos', data = np.array(points))
    f2['mostBasic_world/vessels'].attrs.create('CLASS','REALWORLD')
    
    '''get some hydrodynamics parameters'''
    import krebsjobs.parameters.parameterSetsVesselGen
    set_name = 'default'
    vesselGenParams = getattr(krebsjobs.parameters.parameterSetsVesselGen, set_name)
    #CALCULATE!!!!
    myutils.hdf_write_dict_hierarchy(f2['/'], 'mostBasic_world/parameters', vesselGenParams['calcflow'])
    f2['mostBasic_world/parameters'].attrs.create('name', set_name)
    f2['/'].file.flush()
    #calc_vessel_hydrodynamics(vesselgroup, calc_hematocrit=False, return_flags=False, override_hematocrit = None, bloodflowparams = dict(),storeCalculationInHDF=False)
    dd = ku.calc_vessel_hydrodynamics(f2['mostBasic_world/vessels'], bloodflowparams = vesselGenParams['calcflow'], storeCalculationInHDF = True)
    
    print('Calcflow done')
    


if __name__ == '__main__':

#create stuff  
  fn = 'test_configs.h5'
  f = h5py.File(fn)
  if not 'mostBasic_world' in f:
    createMostBasic_world(f)
  else:
    input_var = raw_input('you want to redo mostBasic_wold: ')
    #print(str(input_var))
    #imput_var = sys.stdin.readline()
    if str(input_var) == 'y':
      del f['mostBasic_world']
      createMostBasic_world(f)
  
