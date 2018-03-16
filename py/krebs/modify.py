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
#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 13:16:06 2016

@author: thierry
"""

if __name__ == '__main__':
  import os.path, sys
  sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'..'))

import myutils
import krebsutils
import itertools
import sys
import numpy as np
import qsub



#enlarge_vessels(float(this_enlargement_factor),out_grp, in_grp, bloodflowparams)
#def worker_on_client(fn, pattern, enlargement_factor):
  

def create_Y_sample(file_to_save):
  radii = [2,2,np.power(2**3+2**3,1/3.)]
  radii = [radius * 1.2 for radius in radii]
  N_edges=len(radii)
  file_to_save.create_group('vessels')
  edgegrp = file_to_save['vessels'].create_group("edges")
  edgegrp.attrs.create('COUNT',N_edges)
  ds_nodeA = edgegrp.create_dataset('node_a_index', data=[0, 1, 3])
  ds_nodeB = edgegrp.create_dataset('node_b_index', data=[3, 3, 2])
  ds_radius = edgegrp.create_dataset('radius', data=radii)
  ds_flags = edgegrp.create_dataset('flags', data=[int(krebsutils.VEIN),int(krebsutils.VEIN),int(krebsutils.VEIN)])
  
  points=[]
  points.append([-100,+100,0])
  points.append([-100,-100,0])
  points.append([0,0,100])
  points.append([0,0,0])

  N_nodes=len(points)
  nodegrp = file_to_save['vessels'].create_group("nodes")
  nodegrp.attrs.create('COUNT', N_nodes)
  ds_roots = file_to_save['vessels/nodes'].create_dataset('roots', data = [0,1,2])
  ds_bc_node_index = file_to_save['vessels/nodes'].create_dataset('bc_node_index', data = [0,1,2])
  
  #dummy pressure for compatiblility
  #ds_pressure = f3['vessels/nodes'].create_dataset('roots_pressure', data = roots_pressure)
  #this works    
  #ds_value_of_bc = f2['symetricA/vessels/nodes'].create_dataset('bc_value', data=[70*0.133,10*0.133])        
  ds_value_of_bc = file_to_save['vessels/nodes'].create_dataset('bc_value', data=[2,2,1])
  ds_conductivity_of_bc = file_to_save['vessels/nodes'].create_dataset('bc_conductivity_value', data=[0,0,0])
  ds_bctyp_of_roots = file_to_save['vessels/nodes'].create_dataset('bc_type', data= [1,1,1])
  ds_world_pos = file_to_save['vessels/nodes'].create_dataset('world_pos', data = np.array(points))
  file_to_save['vessels'].attrs.create('CLASS','REALWORLD')
  
  #no bloodflowParams --> simple --> no hematocrit output
  (pressure, flow, shearforce, flags) = krebsutils.calc_vessel_hydrodynamics(file_to_save['vessels'],return_flags=True)
  # then we save the new data to complete the network copy
  edgegrp.create_dataset('flow'      , data = flow       , compression = 9)
  edgegrp.create_dataset('shearforce', data = shearforce , compression = 9)
  
  nodegrp.create_dataset('pressure'  , data = pressure   , compression = 9)

def enlarge_vessels(factor, gv_filename, bloodflowparams):
  '''gvdst = group where the data is placed in, does not create a 'vesse' folder in it but writes nodes, edges directly;
     gv   = source vessel group
  '''
  inFile = h5files.open(gv_filename, 'r')
  gv = inFile['vessels']
  fac_as_percent = int(np.ceil(float(factor)*100-100))
  dest_file = h5files.open(myutils.strip_from_end(gv_filename, '.h5')+'_growth_by_%02i.h5'%fac_as_percent, 'a')
  dest_file.attrs.create('enlargeFactor', data=fac_as_percent)  
  gvdst = dest_file.create_group('vessels')
  gvdst.attrs['CLASS'] = 'GRAPH'
  myutils.buildLink(gvdst, 'SOURCE', gv)
  # first we need to copy some of the vessel data
  gvedst = gvdst.create_group('edges')
  gvndst = gvdst.create_group('nodes')
  gvedst.attrs['COUNT'] = gv['edges'].attrs['COUNT']
  gvndst.attrs['COUNT'] = gv['nodes'].attrs['COUNT']
  gv.copy('lattice',gvdst)
  for name in ['lattice_pos', 
               'roots',
               'nodeflags','gf',
               'bc_conductivity_value',
               'bc_node_index',
               'bc_type',
               'bc_value']:
    gv['nodes'].copy(name, gvndst)
  for name in ['radius', 'node_a_index', 'node_b_index','flags']:
    if name=='radius':
      radii=gv['edges/'+name]
      radii = factor*np.asarray(radii)
      gvedst.create_dataset(name,data=radii)
    else:
      gv['edges'].copy(name, gvedst)
  # then we recompute blood flow because the alorithm has changed and we may or may not want hematocrit
  (pressure, flow, shearforce, hematocrit, flags) = krebsutils.calc_vessel_hydrodynamics(gvdst,return_flags=True, bloodflowparams=bloodflowparams)
  # then we save the new data to complete the network copy
  gvedst.create_dataset('flow'      , data = flow       , compression = 9)
  gvedst.create_dataset('shearforce', data = shearforce , compression = 9)
  gvedst.create_dataset('hematocrit', data = hematocrit , compression = 9)
  gvndst.create_dataset('pressure'  , data = pressure   , compression = 9)
  
def DoIt(inputFileNames, pattern, options):
  for this_enlargement_factor in options.enlarge_factor:
    inFiles = [h5files.open(fn, 'r') for fn in inputFileNames]
    
    
    inGroups = list(itertools.chain.from_iterable(myutils.walkh5(f, pattern, return_h5objects=True) for f in inFiles))
    if len(inGroups)<=0:
      print 'no matching groups in hdf file(s)'
      sys.exit(0)
    for in_file in inputFileNames:
      bloodflowparams = krebsutils.pyDictFromParamGroup(h5files.open(in_file,'r')['parameters/calcflow'])
      #enlarge_vessels(float(this_enlargement_factor),in_file, bloodflowparams)
      
      qsub.submit(qsub.func(enlarge_vessels, float(this_enlargement_factor), in_file,bloodflowparams),
                  name = 'job_modify_enlarge_'+str(this_enlargement_factor)+'_',
                  num_cpus = 6,
                  days = 5,  # about one day per thread, the idea being that number of threads is proportional to systems size and runtime is 
                  mem = '%iMB' % (1000),
                  change_cwd = True)
if __name__=='__main__':
  import argparse
  parser = argparse.ArgumentParser(description='')  
  parser.add_argument('-f',
                      action='store',
                      help = 'fator by which the vessels are enlarged', 
                      dest='enlarge_factor',
                      nargs='+', type=float                      
                      )
  parser.add_argument('args', nargs='*')
  argv = qsub.parse_args(sys.argv)
  parseResult = parser.parse_args(argv[1:])
  if 1:
    filenames, pattern = parseResult.args[0:-1], parseResult.args[-1]  
    DoIt(filenames, pattern, parseResult)
  if 0:#realy simple setting  
    f = h5files.open('test_set12.h5', 'a')
    create_Y_sample(f)