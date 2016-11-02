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
  sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../..'))
  
import numpy as np
import krebsutils as ku
import math
import collections
import itertools
import h5files
import myutils
import os
import sys

import krebs.adaption as adap
import krebsjobs.parameters.parameterSetsAdaption

from matplotlib.backends.backend_pdf import PdfPages
import h5py
import matplotlib.pyplot as plt

def scale_points(points):
      points2 = []
      scale = 100
      for point in points:
        points2.append([point[0]*scale,point[1]*scale,point[2]*scale])
      return points2

def createAsymetric_world(fn):
    points=[]
    points.append([0,0,0])
    points.append([200,0,0])
    points.append([400,0,0])
    points.append([600,0,0])
    points.append([800,0,0])
    points.append([1000,0,0])
    points.append([1200,0,0])
    points.append([1400,0,0])
    points.append([1450,120,0])
    points.append([1500,240,0])
    points.append([1500,360,0])
    points.append([1450,480,0])
    points.append([1400,600,0])
    points.append([1350,480,0])
    points.append([1300,360,0])
    points.append([1300,240,0])
    points.append([1350,120,0])
    points.append([1200,600,0])
    points.append([1000,600,0])
    points.append([800,600,0])
    points.append([600,600,0])
    points.append([400,600,0])
    points.append([200,600,0])
    points.append([0,600,0])
    points.append([1200,360,0])
    points.append([1000,360,0])
    points.append([800,360,0])
    points.append([600,360,0])
    points.append([400,360,0])
    points.append([200,360,0])
    points.append([1200,240,0])
    points.append([1000,240,0])
    points.append([800,240,0])
    points.append([600,240,0])
    points.append([400,240,0])
    points.append([200,240,0])
    
    def scale_points(points):
      points2 = []
      scale = 1
      for point in points:
            points2.append([-point[0]*scale,-point[1]*scale,point[2]]*scale)
      return points2
    
    #points = reverse_xy(points)    
    points = scale_points(points)    
    
    node_a_index=[]
    node_b_index=[]
    flags=[]
    radii = []
    
    #1
    node_a_index.append(0)
    node_b_index.append(1)
    flags.append(int(ku.ARTERY))
    radii.append(7.5)
    
    #2
    node_a_index.append(1)
    node_b_index.append(2)
    flags.append(int(ku.ARTERY))
    radii.append(7.)
    #3
    node_a_index.append(2)
    node_b_index.append(3)
    flags.append(int(ku.ARTERY))
    radii.append(6.5)
    #4
    node_a_index.append(3)
    node_b_index.append(4)
    flags.append(int(ku.ARTERY))
    radii.append(6.0)
    #5
    node_a_index.append(4)
    node_b_index.append(5)
    flags.append(int(ku.ARTERY))
    radii.append(5.8)
    #6
    node_a_index.append(5)
    node_b_index.append(6)
    flags.append(int(ku.ARTERY))
    radii.append(5.6)
    #7
    node_a_index.append(6)
    node_b_index.append(7)
    flags.append(int(ku.ARTERY))
    radii.append(5.4)
    #8
    node_a_index.append(7)
    node_b_index.append(8)
    flags.append(int(ku.ARTERY))
    radii.append(5.2)
    #9
    node_a_index.append(8)
    node_b_index.append(9)
    flags.append(int(ku.ARTERY))
    radii.append(5.)
    #10
    node_a_index.append(9)
    node_b_index.append(10)
    flags.append(int(ku.CAPILLARY))
    radii.append(4.5)
    
    #11
    node_a_index.append(10)
    node_b_index.append(11)
    flags.append(int(ku.VEIN))
    radii.append(7.)
    #12
    node_a_index.append(11)
    node_b_index.append(12)
    flags.append(int(ku.VEIN))
    radii.append(7.)
    #13
    node_a_index.append(12)
    node_b_index.append(13)
    flags.append(int(ku.VEIN))
    radii.append(7.)
    #14
    node_a_index.append(13)
    node_b_index.append(14)
    flags.append(int(ku.VEIN))
    radii.append(7.)
    #15
    node_a_index.append(14)
    node_b_index.append(15)
    flags.append(int(ku.CAPILLARY))
    radii.append(4.5)
    #16
    node_a_index.append(15)
    node_b_index.append(16)
    flags.append(int(ku.ARTERY))
    radii.append(7.)
    #17
    node_a_index.append(16)
    node_b_index.append(7)
    flags.append(int(ku.ARTERY))
    radii.append(7.)
    #18
    node_a_index.append(12)
    node_b_index.append(17)
    flags.append(int(ku.VEIN))
    radii.append(7.)
    #19
    node_a_index.append(17)
    node_b_index.append(18)
    flags.append(int(ku.VEIN))
    radii.append(7.)
    #20
    node_a_index.append(18)
    node_b_index.append(19)
    flags.append(int(ku.VEIN))
    radii.append(7.)
    #21
    node_a_index.append(19)
    node_b_index.append(20)
    flags.append(int(ku.VEIN))
    radii.append(7.)
    #22
    node_a_index.append(20)
    node_b_index.append(21)
    flags.append(int(ku.VEIN))
    radii.append(7.)
    #23
    node_a_index.append(21)
    node_b_index.append(22)
    flags.append(int(ku.VEIN))
    radii.append(7.)
    #24
    node_a_index.append(22)
    node_b_index.append(23)
    flags.append(int(ku.VEIN))
    radii.append(7.)
    #25
    node_a_index.append(17)
    node_b_index.append(24)
    flags.append(int(ku.VEIN))
    radii.append(7.)
    #26
    node_a_index.append(18)
    node_b_index.append(25)
    flags.append(int(ku.VEIN))
    radii.append(7.)
    #27
    node_a_index.append(19)
    node_b_index.append(26)
    flags.append(int(ku.VEIN))
    radii.append(7.)
    #28
    node_a_index.append(20)
    node_b_index.append(27)
    flags.append(int(ku.VEIN))
    radii.append(7.)
    #29
    node_a_index.append(21)
    node_b_index.append(28)
    flags.append(int(ku.VEIN))
    radii.append(7.)
    #30
    node_a_index.append(22)
    node_b_index.append(29)
    flags.append(int(ku.VEIN))
    radii.append(7.)
    #31
    node_a_index.append(24)
    node_b_index.append(30)
    flags.append(int(ku.CAPILLARY))
    radii.append(4.5)
    #32
    node_a_index.append(25)
    node_b_index.append(31)
    flags.append(int(ku.CAPILLARY))
    radii.append(4.5)
    #33
    node_a_index.append(26)
    node_b_index.append(32)
    flags.append(int(ku.CAPILLARY))
    radii.append(4.5)
    #34
    node_a_index.append(27)
    node_b_index.append(33)
    flags.append(int(ku.CAPILLARY))
    radii.append(4.5)
    #35
    node_a_index.append(28)
    node_b_index.append(34)
    flags.append(int(ku.CAPILLARY))
    radii.append(4.5)
    #36
    node_a_index.append(29)
    node_b_index.append(35)
    flags.append(int(ku.CAPILLARY))
    radii.append(4.5)
    #37
    node_a_index.append(30)
    node_b_index.append(6)
    flags.append(int(ku.ARTERY))
    radii.append(7.)
    #38
    node_a_index.append(31)
    node_b_index.append(5)
    flags.append(int(ku.ARTERY))
    radii.append(7.)
    #39
    node_a_index.append(32)
    node_b_index.append(4)
    flags.append(int(ku.ARTERY))
    radii.append(7.)
    #40
    node_a_index.append(33)
    node_b_index.append(3)
    flags.append(int(ku.ARTERY))
    radii.append(7.)
    #41
    node_a_index.append(34)
    node_b_index.append(2)
    flags.append(int(ku.ARTERY))
    radii.append(7.)
    #42
    node_a_index.append(35)
    node_b_index.append(1)
    flags.append(int(ku.ARTERY))
    radii.append(7.)
    
    #add circulated flag
    flags2=[]
    for aflag in flags:
        flags2.append(np.bitwise_or(ku.CIRCULATED,aflag))
        
    #radii=30*np.ones(len(node_a_index))
    #write the data to a h5 file
    #f3 = h5files.open(fn,'w')
    f3 = fn
    f3.create_group('Asymetric/vessels')
    
    #edge stuff
    N_edges=len(radii)
    edgegrp = f3['Asymetric/vessels'].create_group("edges")
    edgegrp.attrs.create('COUNT',N_edges)
    ds_nodeA = edgegrp.create_dataset('node_a_index', data=node_a_index)
    ds_nodeB = edgegrp.create_dataset('node_b_index', data= node_b_index)
    ds_radius = edgegrp.create_dataset('radius', data=radii)
    #ds_hema = edgegrp.create_dataset('hematocrit', data=hema)
    #ds_flow = edgegrp.create_dataset('flow', data=flow)
    
    #mw_vessel_flags = get_flags_from_other_file()
    
    ds_vesselflags = edgegrp.create_dataset('flags', data = flags2)
    
    #node stuff
    N_nodes=len(points)
    nodegrp = f3['Asymetric/vessels'].create_group("nodes")
    nodegrp.attrs.create('COUNT', N_nodes)
    ds_bc_node_index = f3['Asymetric/vessels/nodes'].create_dataset('bc_node_index', data = [0,23])
    ds_roots = f3['Asymetric/vessels/nodes'].create_dataset('roots', data = [0,23])

    #dummy pressure for compatiblility
    #ds_pressure = f3['vessels/nodes'].create_dataset('roots_pressure', data = roots_pressure)
    ds_value_of_bc = f3['Asymetric/vessels/nodes'].create_dataset('bc_value', data=[70*0.133,2.5e6])
    """ Michael 1: press, 2, flow, 
        Secomb 0: press, 2, flow
    """        
    ds_bctyp_of_roots = f3['Asymetric/vessels/nodes'].create_dataset('bc_type', data= [1,2])
    ds_world_pos = f3['Asymetric/vessels/nodes'].create_dataset('world_pos', data = np.array(points))

    f3['Asymetric/vessels'].attrs.create('CLASS','REALWORLD')    
    
    import krebsjobs.parameters.parameterSetsAdaption
    set_name = 'adaption_asymetric_world'
    adaptionParams = getattr(parameterSets, set_name)
    #CALCULATE!!!!
    myutils.hdf_write_dict_hierarchy(f3['/'], 'Asymetric/parameters', adaptionParams)
    f3['Asymetric/parameters'].attrs.create('name', set_name)
    f3['/'].file.flush()
    dd = ku.calc_vessel_hydrodynamics(f3['Asymetric/vessels'],adaptionParams['calcflow']['includePhaseSeparationEffect'] , False, None, adaptionParams['calcflow'])
    print('Calcflow done')
    
    #r = adap.computeAdaption(f3['vessels'], None, adaptionParams, adaptionParams['calcflow'], f3['/'] )
    dd = adap.computeAdaption_(f3['Asymetric'].create_group('vessels_after_adaption'), f3['Asymetric/vessels'], None, adaptionParams)
    
    
    #os.system("python2 /localdisk/thierry/tumorcode/py/krebs/hdfvessels2vtk.py test_configs.h5 'conduc' '/Asymetric/vessels_after_adaption/vessels_after_adaption'")
    #os.system("python2 /localdisk/thierry/tumorcode/py/krebs/povrayRenderVessels.py test_configs.h5 /Asymetric/vessels_after_adaption/vessels_after_adaption")

def createSymetric_worldA(fn):    
    
    points=[]
    points.append([-10,0,0])
    p1=[-6,4,0]
    p2=[-3,6,0]
    p3=[-1,7,0]
    
    p4=[-1,5,0]
    
    p5=[-3,2,0]
    p6=[-1,3,0]
    p7=[-1,1,0]
    
    points.append(p1)
    points.append(p2)
    points.append(p3)
    points.append(p4)
    points.append(p5)
    points.append(p6)
    points.append(p7)
    
    #mirror y axis
    for i in range(len(points)):
        points.append([-points[i][0],points[i][1],points[i][2]])
    
    node_a_index=[0,1,2,2, 3,11,12,10,9,1,5, 6,14,5, 7,15,13, 4, 0,17,18,18,19,27,28,26,25,17,21,22,30,21,23,31,29,20,32,8]
    node_b_index=[1,2,3,4,11,10,10, 9,8,5,6,14,13,7,15,13,9 ,12,17,18,19,20,27,26,26,25, 8,21,22,30,29,23,31,29,25,28,0,33]
    flags=[]
    radii=[]
#    level0 = 12.
#    level1 = 10.1
#    level2 = 8.5
#    level3 = 7.1
    level0 = 11.95
    level1 = 9.55
    level2 = 7.55
    level3 = 6.0
    flags.append(int(ku.ARTERY))#0
    radii.append(level1)
    flags.append(int(ku.ARTERY))#1
    radii.append(level2)
    flags.append(int(ku.ARTERY))#2
    radii.append(level3)
    flags.append(int(ku.ARTERY))#3
    radii.append(level3)
    flags.append(int(ku.CAPILLARY))#4
    radii.append(level3)
    flags.append(int(ku.VEIN))#5
    radii.append(level3)
    flags.append(int(ku.VEIN))#6
    radii.append(level3)
    flags.append(int(ku.VEIN))#7
    radii.append(level2)
    flags.append(int(ku.VEIN))#8
    radii.append(level1)
    flags.append(int(ku.ARTERY))#9
    radii.append(level2)
    flags.append(int(ku.ARTERY))#10
    radii.append(level3)
    flags.append(int(ku.CAPILLARY))#11
    radii.append(level3)
    flags.append(int(ku.VEIN))#12
    radii.append(level3)
    flags.append(int(ku.ARTERY))#13
    radii.append(level3)
    flags.append(int(ku.CAPILLARY))#14
    radii.append(level3)
    flags.append(int(ku.VEIN))#15
    radii.append(level3)
    flags.append(int(ku.VEIN))#16
    radii.append(level2)
    flags.append(int(ku.CAPILLARY))#17
    radii.append(level3)
    
    #mirror x axis
    for i in range(len(points)):
        points.append([points[i][0],-points[i][1],points[i][2]])
               
        
    #double flags
    for i in range(len(flags)):
        flags.append(flags[i])
        radii.append(radii[i]) 
    
    #add entrence
    flags.append(int(ku.ARTERY))
    radii.append(level0)
    points.append([-13,0,0])
    flags.append(int(ku.VEIN))
    radii.append(level0)
    points.append([13,0,0])
    
    def reverse_xy(points):
        points2 = []
        scale = 100
        for point in points:
            points2.append([point[1]*scale,point[0]*scale,point[2]*scale])
        return points2
    
    #points = reverse_xy(points)    
    points = scale_points(points)
    #add circulated flag
    flags2=[]
    for aflag in flags:
        flags2.append(np.bitwise_or(ku.CIRCULATED,aflag))
       
    #this datalist is after successfull adaption with all signals   
    datalist=[13.301738,9.563003,7.1691427,7.108207,5.7596083,5.6054125,5.584573,6.330745,8.065095,9.596836,7.166868,
              5.7524166,5.6175594,7.166868,5.7524166,5.6175594,6.354387,5.704988,13.301738,9.563003,7.1691427,7.108207,
              5.7596083,5.6054125,5.584573,6.330745,8.065095,9.596836,7.166868,5.7524166,5.6175594,7.166868,5.7524166,
              5.6175594,6.354387,5.704988,17.057396,9.972254]
    #radii = np.array(datalist)        
    #radii=10*np.ones(len(node_a_index))
    #radii=6.5*np.ones(len(node_a_index))
    #write the data to a h5 file
    #f2 = h5files.open(fn,'w')
    f2 = fn
    f2.create_group('symetricA/vessels')
    
    #edge stuff
    N_edges=len(radii)
    edgegrp = f2['symetricA/vessels'].create_group("edges")
    edgegrp.attrs.create('COUNT',N_edges)
    ds_nodeA = edgegrp.create_dataset('node_a_index', data=node_a_index)
    ds_nodeB = edgegrp.create_dataset('node_b_index', data= node_b_index)
    ds_radius = edgegrp.create_dataset('radius', data=radii)
    #ds_hema = edgegrp.create_dataset('hematocrit', data=hema)
    #ds_flow = edgegrp.create_dataset('flow', data=flow)
    
    #mw_vessel_flags = get_flags_from_other_file()
    
    ds_vesselflags = edgegrp.create_dataset('flags', data = flags2)
    
    #node stuff
    N_nodes=len(points)
    nodegrp = f2['symetricA/vessels'].create_group("nodes")
    nodegrp.attrs.create('COUNT', N_nodes)
    ds_roots = f2['symetricA/vessels/nodes'].create_dataset('roots', data = [32,33])
    ds_bc_node_index = f2['symetricA/vessels/nodes'].create_dataset('bc_node_index', data = [32,33])
    
    #dummy pressure for compatiblility
    #ds_pressure = f3['vessels/nodes'].create_dataset('roots_pressure', data = roots_pressure)
    #this works    
    #ds_value_of_bc = f2['symetricA/vessels/nodes'].create_dataset('bc_value', data=[70*0.133,10*0.133])        
    ds_value_of_bc = f2['symetricA/vessels/nodes'].create_dataset('bc_value', data=[-8e5, 15*0.133])        
    ds_bctyp_of_roots = f2['symetricA/vessels/nodes'].create_dataset('bc_type', data= [2,1])
    ds_world_pos = f2['symetricA/vessels/nodes'].create_dataset('world_pos', data = np.array(points))
    f2['symetricA/vessels'].attrs.create('CLASS','REALWORLD')
    
    import krebsjobs.parameters.parameterSetsAdaption
    #set_name = 'adaption_symetric_world'
    set_name = 'adaption_symetricA_world_break'
    #set_name = 'adaption_symetric_world'
    adaptionParams = getattr(parameterSets, set_name)
    #CALCULATE!!!!
    myutils.hdf_write_dict_hierarchy(f2['/'], 'symetricA/parameters', adaptionParams)
    f2['symetricA/parameters'].attrs.create('name', set_name)
    f2['/'].file.flush()
    dd = ku.calc_vessel_hydrodynamics(f2['symetricA/vessels'],adaptionParams['calcflow']['includePhaseSeparationEffect'] , False, None, adaptionParams['calcflow'])
    
    print('Calcflow done')
    
    #r = adap.computeAdaption(f3['vessels'], None, adaptionParams, adaptionParams['calcflow'], f3['/'] )
    #dst_grp = f2['/symetric'].create_group('symetric/vessels_after_adaption')
    dd = adap.computeAdaption_(f2['/symetricA'].create_group('vessels_after_adaption'), f2['symetricA/vessels'], None, adaptionParams)

#    a = os.system("python2 /localdisk/thierry/tumorcode/py/krebs/hdfvessels2vtk.py test_configs.h5 'conduc' /symetricA/vessels_after_adaption/vessels_after_adaption")
#    os.system("python2 /localdisk/thierry/tumorcode/py/krebs/povrayRenderVessels.py test_configs.h5 /symetricA/vessels_after_adaption/vessels_after_adaption") 
 

def createSymetric_worldB(fn):    
    
    points=[]
    points.append([-10,0,0])
    p1=[-6,4,0]
    p2=[-3,6,0]
    p3=[-1,7,0]
    
    p4=[-1,5,0]
    
    p5=[-3,2,0]
    p6=[-1,3,0]
    p7=[-1,1,0]
    
    points.append(p1)
    points.append(p2)
    points.append(p3)
    points.append(p4)
    points.append(p5)
    points.append(p6)
    points.append(p7)
    
    #mirror y axis
    for i in range(len(points)):
        points.append([-points[i][0],points[i][1],points[i][2]])
    
    node_a_index=[0,1,2,2, 3,11,12,10,9,1,5, 6,14,5, 7,15,13, 4, 0,17,18,18,19,27,28,26,25,17,21,22,30,21,23,31,29,20,32,8]
    node_b_index=[1,2,3,4,11,10,10, 9,8,5,6,14,13,7,15,13,9 ,12,17,18,19,20,27,26,26,25, 8,21,22,30,29,23,31,29,25,28,0,33]
    flags=[]
    radii=[]
#    level0 = 12.
#    level1 = 10.1
#    level2 = 8.5
#    level3 = 7.1
    level0 = 11.95
    level1 = 9.55
    level2 = 7.55
    level3 = 6.0
    flags.append(int(ku.ARTERY))#0
    radii.append(level1)
    flags.append(int(ku.ARTERY))#1
    radii.append(level2)
    flags.append(int(ku.ARTERY))#2
    radii.append(level3)
    flags.append(int(ku.ARTERY))#3
    radii.append(level3)
    flags.append(int(ku.CAPILLARY))#4
    radii.append(level3)
    flags.append(int(ku.VEIN))#5
    radii.append(level3)
    flags.append(int(ku.VEIN))#6
    radii.append(level3)
    flags.append(int(ku.VEIN))#7
    radii.append(level2)
    flags.append(int(ku.VEIN))#8
    radii.append(level1)
    flags.append(int(ku.ARTERY))#9
    radii.append(level2)
    flags.append(int(ku.ARTERY))#10
    radii.append(level3)
    flags.append(int(ku.CAPILLARY))#11
    radii.append(level3)
    flags.append(int(ku.VEIN))#12
    radii.append(level3)
    flags.append(int(ku.ARTERY))#13
    radii.append(level3)
    flags.append(int(ku.CAPILLARY))#14
    radii.append(level3)
    flags.append(int(ku.VEIN))#15
    radii.append(level3)
    flags.append(int(ku.VEIN))#16
    radii.append(level2)
    flags.append(int(ku.CAPILLARY))#17
    radii.append(level3)
    
    #mirror x axis
    for i in range(len(points)):
        points.append([points[i][0],-points[i][1],points[i][2]])
               
        
    #double flags
    for i in range(len(flags)):
        flags.append(flags[i])
        radii.append(radii[i]) 
    
    #add entrence
    flags.append(int(ku.ARTERY))
    radii.append(level0)
    points.append([-13,0,0])
    flags.append(int(ku.VEIN))
    radii.append(level0)
    points.append([13,0,0])
    
    def reverse_xy(points):
        points2 = []
        scale = 100
        for point in points:
            points2.append([point[1]*scale,point[0]*scale,point[2]]*scale)
        return points2
    
    #points = reverse_xy(points)    
    points = scale_points(points)
    #add circulated flag
    flags2=[]
    for aflag in flags:
        flags2.append(np.bitwise_or(ku.CIRCULATED,aflag))
       
    #this datalist is after successfull adaption with all signals   
    datalist=[13.301738,9.563003,7.1691427,7.108207,5.7596083,5.6054125,5.584573,6.330745,8.065095,9.596836,7.166868,
              5.7524166,5.6175594,7.166868,5.7524166,5.6175594,6.354387,5.704988,13.301738,9.563003,7.1691427,7.108207,
              5.7596083,5.6054125,5.584573,6.330745,8.065095,9.596836,7.166868,5.7524166,5.6175594,7.166868,5.7524166,
              5.6175594,6.354387,5.704988,17.057396,9.972254]
    #radii = np.array(datalist)        
    #radii=10*np.ones(len(node_a_index))
    #radii=6.5*np.ones(len(node_a_index))
    #write the data to a h5 file
    #f2 = h5files.open(fn,'w')
    f2 = fn
    f2.create_group('symetricB/vessels')
    
    #edge stuff
    N_edges=len(radii)
    edgegrp = f2['symetricB/vessels'].create_group("edges")
    edgegrp.attrs.create('COUNT',N_edges)
    ds_nodeA = edgegrp.create_dataset('node_a_index', data=node_a_index)
    ds_nodeB = edgegrp.create_dataset('node_b_index', data= node_b_index)
    ds_radius = edgegrp.create_dataset('radius', data=radii)
    #ds_hema = edgegrp.create_dataset('hematocrit', data=hema)
    #ds_flow = edgegrp.create_dataset('flow', data=flow)
    
    #mw_vessel_flags = get_flags_from_other_file()
    
    ds_vesselflags = edgegrp.create_dataset('flags', data = flags2)
    
    #node stuff
    N_nodes=len(points)
    nodegrp = f2['symetricB/vessels'].create_group("nodes")
    nodegrp.attrs.create('COUNT', N_nodes)
    ds_roots = f2['symetricB/vessels/nodes'].create_dataset('roots', data = [32,33])
    ds_bc_node_index = f2['symetricB/vessels/nodes'].create_dataset('bc_node_index', data = [32,33])
    
    #dummy pressure for compatiblility
    #ds_pressure = f3['vessels/nodes'].create_dataset('roots_pressure', data = roots_pressure)
    ds_value_of_bc = f2['symetricB/vessels/nodes'].create_dataset('bc_value', data=[70*0.133,0.5e6])        
    ds_bctyp_of_roots = f2['symetricB/vessels/nodes'].create_dataset('bc_type', data= [1,2])
    ds_world_pos = f2['symetricB/vessels/nodes'].create_dataset('world_pos', data = np.array(points))
    f2['symetricB/vessels'].attrs.create('CLASS','REALWORLD')
    
    import krebsjobs.parameters.parameterSetsAdaption
    set_name = 'adaption_symetricB_world'
    adaptionParams = getattr(parameterSets, set_name)
    #CALCULATE!!!!
    myutils.hdf_write_dict_hierarchy(f2['/'], 'symetricB/parameters', adaptionParams)
    f2['symetricB/parameters'].attrs.create('name', set_name)
    f2['/'].file.flush()
    dd = ku.calc_vessel_hydrodynamics(f2['symetricB/vessels'],adaptionParams['calcflow']['includePhaseSeparationEffect'] , False, None, adaptionParams['calcflow'])
    
    print('Calcflow done')
    
    #r = adap.computeAdaption(f3['vessels'], None, adaptionParams, adaptionParams['calcflow'], f3['/'] )
    #dst_grp = f2['/symetric'].create_group('symetric/vessels_after_adaption')
    dd = adap.computeAdaption_(f2['/symetricB'].create_group('vessels_after_adaption'), f2['symetricB/vessels'], None, adaptionParams)
    

    a = os.system("python2 /localdisk/thierry/tumorcode/py/krebs/hdfvessels2vtk.py test_configs.h5 /symetricB/vessels_after_adaption/vessels_after_adaption")
    os.system("python2 /localdisk/thierry/tumorcode/py/krebs/povrayRenderVessels.py test_configs.h5 /symetricB/vessels_after_adaption/vessels_after_adaption") 

def createSymetricIrregular_world(fn):
    points=[]
    points.append([-10,0,0])
    p1=[-6,4,0]
    p2=[-3,6,0]
    p3=[-1,7,0]
    
    p4=[-1,5,0]
    
    p5=[-3,2,0]
    p6=[-1,3,0]
    p7=[-1,1,0]
    
    points.append(p1)
    points.append(p2)
    points.append(p3)
    points.append(p4)
    points.append(p5)
    points.append(p6)
    points.append(p7)
    
    #mirror y axis
    for i in range(len(points)):
        points.append([-points[i][0],points[i][1],points[i][2]])
    
    node_a_index=[0,1,2,2, 3,11,12,10,9,1,5, 6,14,5, 7,15,13, 4, 0,17,18,18,19,27,28,26,25,17,21,22,30,21,23,31,29,20,32,8]
    node_b_index=[1,2,3,4,11,10,10, 9,8,5,6,14,13,7,15,13,9 ,12,17,18,19,20,27,26,26,25, 8,21,22,30,29,23,31,29,25,28,0,33]
    flags=[]
    radii=[]
    flags.append(int(ku.ARTERY))#0
    radii.append(9.5)
    flags.append(int(ku.ARTERY))#1
    radii.append(8)
    flags.append(int(ku.ARTERY))#2
    radii.append(7.5)
    flags.append(int(ku.ARTERY))#3
    radii.append(7.)
    flags.append(int(ku.CAPILLARY))#4
    radii.append(6.5)
    flags.append(int(ku.VEIN))#5
    radii.append(7.)
    flags.append(int(ku.VEIN))#6
    radii.append(7.5)
    flags.append(int(ku.VEIN))#7
    radii.append(8)
    flags.append(int(ku.VEIN))#8
    radii.append(9.5)
    flags.append(int(ku.ARTERY))#9
    radii.append(8)
    flags.append(int(ku.ARTERY))#10
    radii.append(7.5)
    flags.append(int(ku.CAPILLARY))#11
    radii.append(6.5)
    flags.append(int(ku.VEIN))#12
    radii.append(7.5)
    flags.append(int(ku.ARTERY))#13
    radii.append(7.5)
    flags.append(int(ku.CAPILLARY))#14
    radii.append(6.5)
    flags.append(int(ku.VEIN))#15
    radii.append(7.5)
    flags.append(int(ku.VEIN))#16
    radii.append(8)
    flags.append(int(ku.CAPILLARY))#17
    radii.append(6.5)
    
    #mirror x axis
    for i in range(len(points)):
        points.append([points[i][0],-points[i][1],points[i][2]])
               
        
    #double flags
    for i in range(len(flags)):
        flags.append(flags[i])
        radii.append(radii[i]) 
    
    #add entrence
    flags.append(int(ku.ARTERY))
    radii.append(11.5)
    points.append([-13,0,0])
    flags.append(int(ku.VEIN))
    radii.append(11.5)
    points.append([13,0,0])
    
    def reverse_xy(points):
        points2 = []
        scale = 100
        for point in points:
            points2.append([point[1]*scale,point[0]*scale,point[2]]*scale)
        return points2
    def scale_points(points):
      points2 = []
      scale = 100
      for point in points:
            points2.append([point[0]*scale,point[1]*scale,point[2]]*scale)
      return points2
    
    #points = reverse_xy(points)    
    points = scale_points(points)

    
    #add circulated flag
    flags2=[]
    for aflag in flags:
        flags2.append(np.bitwise_or(ku.CIRCULATED,aflag))
        
    
    #radii=30*np.ones(len(node_a_index))
    #radii=10*np.ones(len(node_a_index))
    #now we will randomize the points a little bit
    def screw_points(points,variation):
        points2 = []
        for point in points:
            points2.append([point[0]+point[0]*np.random.normal(0,variation),
            point[1]+point[1]*np.random.normal(0,2*variation),
            point[2]+point[2]*np.random.normal(0,variation)])
        return points2
    
    if os.path.isfile('screwed_points.h5'):
      f_out = h5py.File('screwed_points.h5','r')
      points=np.asarray(f_out['screwed_points'])
    else:
      points = screw_points(points,0.05)
      f_out = h5py.File('screwed_points.h5','w')       
      f_out.create_dataset('screwed_points', data = points)
      f_out.close()
    #write the data to a h5 file
    #f2 = h5files.open(fn,'w')
    f2 = fn
    f2.create_group('symetricIrregular/vessels')
    
    #edge stuff
    N_edges=len(radii)
    edgegrp = f2['symetricIrregular/vessels'].create_group("edges")
    edgegrp.attrs.create('COUNT',N_edges)
    ds_nodeA = edgegrp.create_dataset('node_a_index', data=node_a_index)
    ds_nodeB = edgegrp.create_dataset('node_b_index', data= node_b_index)
    ds_radius = edgegrp.create_dataset('radius', data=radii)
    #ds_hema = edgegrp.create_dataset('hematocrit', data=hema)
    #ds_flow = edgegrp.create_dataset('flow', data=flow)
    
    #mw_vessel_flags = get_flags_from_other_file()
    
    ds_vesselflags = edgegrp.create_dataset('flags', data = flags2)
    
    #node stuff
    N_nodes=len(points)
    nodegrp = f2['symetricIrregular/vessels'].create_group("nodes")
    nodegrp.attrs.create('COUNT', N_nodes)
    ds_roots = f2['symetricIrregular/vessels/nodes'].create_dataset('roots', data = [32,33])
    ds_bc_node_index = f2['symetricIrregular/vessels/nodes'].create_dataset('bc_node_index', data = [32,33])
    
    #dummy pressure for compatiblility
    #ds_pressure = f3['vessels/nodes'].create_dataset('roots_pressure', data = roots_pressure)
    ds_value_of_bc = f2['symetricIrregular/vessels/nodes'].create_dataset('bc_value', data=[70*0.133,2.5e6])        
    ds_bctyp_of_roots = f2['symetricIrregular/vessels/nodes'].create_dataset('bc_type', data= [1,2])
    ds_world_pos = f2['symetricIrregular/vessels/nodes'].create_dataset('world_pos', data = np.array(points))
    f2['symetricIrregular/vessels'].attrs.create('CLASS','REALWORLD')
    
    import krebsjobs.parameters.parameterSetsAdaption
    set_name = 'adaption_symetricIrregular_world'
    adaptionParams = getattr(parameterSets, set_name)
    #CALCULATE!!!!
    myutils.hdf_write_dict_hierarchy(f2['/'], 'symetricIrregular/parameters', adaptionParams)
    f2['symetricIrregular/parameters'].attrs.create('name', set_name)
    f2['/'].file.flush()
    dd = ku.calc_vessel_hydrodynamics(f2['symetricIrregular/vessels'],adaptionParams['calcflow']['includePhaseSeparationEffect'] , False, None, adaptionParams['calcflow'])
    
    print('Calcflow done')
    
    #r = adap.computeAdaption(f3['vessels'], None, adaptionParams, adaptionParams['calcflow'], f3['/'] )
    #dst_grp = f2['/symetric'].create_group('symetric/vessels_after_adaption')
    dd = adap.computeAdaption_(f2['/symetricIrregular'].create_group('vessels_after_adaption'), f2['symetricIrregular/vessels'], None, adaptionParams)

    a = os.system("python2 /localdisk/thierry/tumorcode/py/krebs/hdfvessels2vtk.py test_configs.h5 /symetricIrregular/vessels_after_adaption/vessels_after_adaption")
    os.system("python2 /localdisk/thierry/tumorcode/py/krebs/povrayRenderVessels.py test_configs.h5 /symetricIrregular/vessels_after_adaption/vessels_after_adaption") 


if __name__ == '__main__':

#create stuff  
  fn = 'test_configs.h5'
  f = h5py.File(fn)
  if not 'symetricA' in f:
    createSymetric_worldA(f)
    print('Symetric_worldA done!')
  if not 'symetricB' in f:
    createSymetric_worldB(f)
    print('Symetric_worldB done')
  if not 'symetricIrregular' in f: 
    createSymetricIrregular_world(f)
    print('SymetricIrregular_world done')
  if not 'Asymetric' in f:
    createAsymetric_world(f)
    print('Asymetric done')
  
#look at stuff
  common_filename = os.path.splitext(fn)[0]
  pp = PdfPages(common_filename + '_results_test_configs.pdf')


  print("shit")
  #plot_pressure_relation()  
  
#  a = os.system("python2 /localdisk/thierry/tumorcode/py/krebs/hdfvessels2vtk.py test_configs.h5 'conduc' /symetric/vessels_after_adaption/vessels_after_adaption")
#  os.system("echo 'hello'")
#  os.system("python2 /localdisk/thierry/tumorcode/py/krebs/hdfvessels2vtk.py test_configs.h5 'conduc' '/Asymetric/vessels_after_adaption/vessels_after_adaption'")
#  os.system("python2 /localdisk/thierry/tumorcode/py/krebs/povrayRenderVessels.py test_configs.h5 /symetric/vessels_after_adaption/vessels_after_adaption") 
#  os.system("python2 /localdisk/thierry/tumorcode/py/krebs/povrayRenderVessels.py test_configs.h5 /Asymetric/vessels_after_adaption/vessels_after_adaption")
