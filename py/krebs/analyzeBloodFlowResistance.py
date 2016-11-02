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
#!/usr/bin/env python
# -*- coding: utf-8 -*-
if __name__ == '__main__':
  import os.path, sys
  sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'..'))
  
import os,sys
import h5py
import h5files
import numpy as np
import extensions # for asarray with h5py support
import krebsutils
import posixpath
import math
import glob
import itertools
from pprint import pprint

import myutils

from analyzeGeneral   import DataBasicVessel

def ComputeVascularTreeBloodFlowResistances(vessels):
  edgelist, flags, pressure, flow, nodeflags = vessels.edgelist, vessels['flags'], vessels['pressure'], vessels['flow'], vessels['nodeflags']
  
  arterial = {}
  venous = {}
  totalFlow = 0.
  if sys.flags.debug:
      num_circ = 0
      num_artery = 0
      num_vein = 0
      num_boundary = 0
      
  for i, (a,b) in enumerate(edgelist):
    #print("i: %i, a: %i, b: %i" %(i,a,b))
    if sys.flags.debug:
        if( flags[i] & krebsutils.CIRCULATED ):
            num_circ = num_circ+1
        if( flags[i] & krebsutils.ARTERY ):
            num_artery = num_artery+1
        if( flags[i] & krebsutils.VEIN ):
            num_vein = num_vein+1
        if( nodeflags[a] & krebsutils.BOUNDARY or nodeflags[b] & krebsutils.BOUNDARY):
            num_boundary = num_boundary+1
    if not (flags[i] & krebsutils.CIRCULATED): continue
    if (nodeflags[a] & krebsutils.BOUNDARY): continue
    if (nodeflags[a] & krebsutils.BOUNDARY): continue   
    if (flags[i] & krebsutils.ARTERY):
      arterial[i] = (0.5*(pressure[a]+pressure[b]), flow[i])
      totalFlow   += flow[i]
      if sys.flags.debug:
          print("pressure[a]: %f, pressure[b]: %f, flow[i]: %f" % (pressure[a],pressure[b],flow[i]))
    if (flags[i] & krebsutils.CIRCULATED):
      venous[i] = (0.5*(pressure[a]+pressure[b]), flow[i])
  if sys.flags.debug:
      print("num_circ: %i" % num_circ)
      print("num_artery: %i" % num_artery)
      print("num_vein: %i" % num_vein)
      print("num_boundary: %i" % num_boundary)
      print(len(venous))
      print(len(arterial))
  avgVenousPressure = np.average([p for (p,q) in venous.values()], weights = [q for (p,q) in venous.values()])
  avgArterialPressure = np.average([p for (p,q) in arterial.values()], weights = [q for (p,q) in arterial.values()])
  conductivities = {} # Q = S*dp
  for k, (p, q) in arterial.iteritems():
    conductivities[k] = q/(p-avgVenousPressure), q, (p-avgVenousPressure)
  return conductivities, avgVenousPressure, avgArterialPressure, totalFlow


def AddConductivityBoundaryConditions(vessels, vesselgroup, avgConductivity):
  print '---- resistor BCs to %s -----' % str(vesselgroup)
  fudgeFactor = 100.
  rootNodes = vessels.roots
  pressures = vessels['pressure']
  numNodes = len(pressures)
  nodeFlags = krebsutils.edge_to_node_property(numNodes, vessels.edgelist, vessels['flags'], 'or')
  nodeRadi  = krebsutils.edge_to_node_property(numNodes, vessels.edgelist, vessels['radius'], 'avg')
  minp = np.amin(pressures)
  # storage area for h5 datasets
  dataArrays = []
  # use radius^4 weighted average conductivity
  weights = dict()
  for rootNode in rootNodes:
    weights[rootNode] = nodeRadi[rootNode]**4
    assert nodeRadi[rootNode]>10. or pressures[rootNode]>minp*1.1 or (nodeFlags[rootNode]&krebsutils.VEIN)
  weightSum = np.sum(weights.values())  
  for rootNode in rootNodes:
    weight = weights[rootNode]/weightSum
    #weight = 1
    conductivity = fudgeFactor * avgConductivity * weight
    dataArrays.append((rootNode, krebsutils.FLOWBC_RESIST, pressures[rootNode], conductivity))
    print 'node %i -> S=%f, w=%f' % (rootNode, conductivity, weight)
  # write to h5
  dataArrays = zip(*dataArrays)
  names = ['bc_node_index', 'bc_type', 'bc_value', 'bc_conductivity_value']
  gnodes = vesselgroup['nodes']
  for name, values in zip(names, dataArrays):
    if name in gnodes:
      gnodes[name][...] = values
    else:
      gnodes.create_dataset(name, data = values)



if __name__ == "__main__":
  import optparse  #Note: Deprecated since version 2.7. Use argparse instead
  parser = optparse.OptionParser()
  parser.add_option("--add-resistor-bc", dest="add_resistor_bc", help="add resistor flow boundary condition, inplace!", default=False, action="store_true")
  parser.add_option("--force-flow-recompute", dest="force_flow_recompute", help="recompute blood flow", default=False, action="store_true")
  parser.add_option("--write-flow", dest="write_flow", help="writes the recomputed flow; note: hematocrit is completely ignored!", default=False, action="store_true")
  options, args = parser.parse_args()
  if options.write_flow:
    options.force_flow_recompute = True
  
  dataman = myutils.DataManager(100, [DataBasicVessel()])

  filenames, pattern = args[:-1], args[-1]
  files = [h5files.open(fn, 'a' if options.add_resistor_bc else 'r+') for fn in filenames]
  groups = list(itertools.chain.from_iterable(myutils.walkh5(f, pattern, return_h5objects=True) for f in files))

  for vesselgroup in groups:
    if options.force_flow_recompute and not options.add_resistor_bc:
      (pressure, flow, shearforce) = krebsutils.calc_vessel_hydrodynamics(vesselgroup)
      vessels = dataman.obtain_data('vessel_graph', vesselgroup, ['flags', 'radius'])
      vessels.edges['flow'] = flow
      vessels.nodes['pressure'] = pressure
    else:
      vessels = dataman.obtain_data('vessel_graph', vesselgroup, ['flow', 'pressure', 'flags', 'radius'])
      
    def DoBC():
      conductivities, avgVenousPressure, avgArterialPressure, totalFlow = ComputeVascularTreeBloodFlowResistances(vessels)
      avgConductivity = (totalFlow/(avgArterialPressure-avgVenousPressure))
      print 'avgVenousPressure', avgVenousPressure
      print 'avgArterialPressure', avgArterialPressure
      print 'totalFlow', totalFlow
      print 'avgConductivity', avgConductivity
      return avgConductivity
    avgConductivity = DoBC()
    
    if options.add_resistor_bc:
      AddConductivityBoundaryConditions(vessels, vesselgroup, avgConductivity)
      vesselgroup.file.flush()
#    else:
#      xy = []
#      for i, (s, q, dp) in conductivities.iteritems():
#        r = vessels['radius'][i]
#        #print 'S=%f, q=%f, dp=%f, r=%f' % (s, q, dp, r)
#        xy.append((r, s))
#      xy = sorted(xy, key=lambda (a,b): a)
#      xy = np.asarray(xy).transpose()
#      pyplot.plot(xy[0], xy[1])
      #pyplot.show()
    
    if options.force_flow_recompute and options.add_resistor_bc:
      (pressure, flow, shearforce) = krebsutils.calc_vessel_hydrodynamics(vesselgroup)
      vessels.edges['flow'] = flow
      vessels.nodes['pressure'] = pressure
      DoBC()      
      if options.write_flow:
        vesselgroup['nodes/pressure'][...] = pressure
        vesselgroup['edges/flow'][...] = flow
        try:
          vesselgroup['edges/shearforce'][...] = shearforce
        except:
          pass
        vesselgroup.file.flush()
        print 'written new blood flow data'
          