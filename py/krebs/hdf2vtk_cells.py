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
if __name__ == '__main__':
  import os.path, sys
  sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'..'))
  sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'.'))
  
import os, sys
from os.path import join, basename, dirname, splitext
import krebsutils
import extractVtkFields
import pprint
import numpy as np
import h5py
import identifycluster
if identifycluster.getname()=='snowden':
  print('Install vtk on snowden!!!')
  #return 0
else:
  from vtk import *
#from vtk import *
from vtkcommon import *
import extensions
import itertools
import myutils
import fnmatch


def ConvertMyHdfVesselsToVTKPolydata(graph, newflag_for_backward_compatibility, options):
  import posixpath
  edges = graph.edgelist
  pos = graph['position']
  pressure = graph['pressure']
  if newflag_for_backward_compatibility:
    nodeflags = graph['nodeflags']
  hema = graph['hematocrit']
  flow = graph['flow']
  rad = graph['radius']
  flags = graph['flags']
  shearforce = graph['shearforce']
  if 'metabolicSignal' in graph.keys():
    meta_sig = graph['metabolicSignal']
  if 'conductivitySignal' in graph.keys():
    cond_sig = graph['conductivitySignal']
  if 'S_tot' in graph.keys():
    S_tot = graph['S_tot']
  #print('WARNING!!!!! SCALING RADIUS')
  #rad = rad * 3
  
  print(pos)  
  
  num_verts = len(pos)
  num_edges = len(edges)
  print( "v %i e %i" % (num_verts, num_edges))

  pts = vtkPoints()
  pts.SetNumberOfPoints(len(pos))
  for i,(x,y,z) in enumerate(pos):
    pts.SetPoint(i, (float(x), float(y), float(z)))

  cells = asVtkCellArray(edges)

  polydata = vtkPolyData()
  polydata.SetPoints(pts)
  polydata.SetLines(cells)

  if 1:
    #radius
    polydata.GetCellData().AddArray(asVtkArray(rad, "radius", vtkFloatArray))
    #nodes
    node_rad = np.zeros(num_verts)
    node_n   = np.zeros(num_verts)
    for r,(a,b) in itertools.izip(rad,edges): # itertools, or it will blow up in yr face cuz ultra slow python lists and shit will be created
      node_rad[a] += r
      node_rad[b] += r
      node_n[a] += 1.
      node_n[b] += 1.
    node_rad /= np.maximum(1, node_n)
    polydata.GetPointData().AddArray(asVtkArray(node_rad, "point_radius", vtkFloatArray))
    
  if 1:
      #pressure
      polydata.GetPointData().AddArray(asVtkArray(pressure, "pressure_at_node", vtkFloatArray))
      #isnodeboundary = np.asarray(np.bitwise_and(nodeflags,krebsutils.BOUNDARY) > 0, dtype=np.int32)
      #polydata.GetPointData().AddArray(asVtkArray(isnodeboundary, "isNodeBoundary", vtkFloatArray))
  if 1:
    #flags
    iscirc = np.asarray(np.bitwise_and(flags,krebsutils.CIRCULATED) > 0, dtype=np.int32)
    polydata.GetCellData().AddArray(asVtkArray(iscirc, "isCirculated", vtkIntArray))
    
  if 1:
    #flags
    isartery = np.asarray(np.bitwise_and(flags,krebsutils.ARTERY) > 0, dtype=np.int32)
    polydata.GetCellData().AddArray(asVtkArray(isartery, "isArtery", vtkIntArray))
    
  if 1:
    #flags
    isvein = np.asarray(np.bitwise_and(flags,krebsutils.VEIN) > 0, dtype=np.int32)
    polydata.GetCellData().AddArray(asVtkArray(isvein, "isVein", vtkIntArray))

  if 1:
    #flags
    iscap = np.asarray(np.bitwise_and(flags,krebsutils.CAPILLARY) > 0, dtype=np.int32)
    polydata.GetCellData().AddArray(asVtkArray(iscap, "isCapillary", vtkIntArray))
    
  if 1 and newflag_for_backward_compatibility:
    #flags
    isnodeboundary = np.asarray(np.bitwise_and(nodeflags,krebsutils.BOUNDARY) > 0, dtype=np.int32)    
    polydata.GetPointData().AddArray(asVtkArray(isnodeboundary, "isBoundary", vtkFloatArray))
    polydata.GetCellData().AddArray(asVtkArray(graph['edge_boundary'], "isBoundaryVessel", vtkIntArray))
    #polydata.GetPointData().AddArray(asVtkArray(isnodeboundary, "isNodeBoundary", vtkFloatArray))
    #isboundary = np.asarray(np.bitwise_and(flags,krebsutils.BOUNDARY) > 0, dtype=np.int32)
    #polydata.GetCellData().AddArray(asVtkArray(isboundary, "isBoundary", vtkIntArray))
    
  if 1:
    #indices of edges
    print(len(edges))
    polydata.GetCellData().AddArray(asVtkArray(np.arange(0,len(edges)), "IndexOfEdges2", vtkIntArray))
    polydata.GetCellData().AddArray(asVtkArray(edges, "IndexOfEdges", vtkIntArray))
    polydata.GetCellData().AddArray(asVtkArray(hema, "hematocrit", vtkFloatArray))
    #flow to nl
    if options.out_nl:
      flow = np.multiply(flow,60./1000000.)
    polydata.GetCellData().AddArray(asVtkArray(flow, "flow", vtkFloatArray))
    #shearforce from kpa to dyne/cm
    shearforce = np.multiply(shearforce, 10000)
    polydata.GetCellData().AddArray(asVtkArray(shearforce, "shearforce", vtkFloatArray))

  if 'cond_sig' in locals():
    if np.sum(np.isnan(cond_sig))==0 and np.sum(np.isinf(cond_sig))==0:
      #normalize to 0,1
      #cond_sig = np.multiply(cond_sig,1/cond_sig.max())
      polydata.GetCellData().AddArray(asVtkArray(cond_sig, "conductive", vtkFloatArray))
    else:
      print("Warning: bad cond_sig value")
      print( np.where(np.isnan(cond_sig) ==1)[0])
      bad_cond_sig = np.isnan(cond_sig)
      cond_sig[bad_cond_sig] = 0
      polydata.GetCellData().AddArray(asVtkArray(cond_sig, "conductive", vtkFloatArray))
  if 'meta_sig' in locals():
    if np.sum(np.isnan(meta_sig))==0 and np.sum(np.isinf(meta_sig))==0:
      polydata.GetCellData().AddArray(asVtkArray(meta_sig, "metabolic", vtkFloatArray))
    else:
      print("warning: bad meta_sig value")
      print( np.where(np.isnan(meta_sig) ==1)[0])
      print( np.where(np.isinf(meta_sig) ==1)[0])
      bad_meta_sig_nan = np.isnan(meta_sig)
      bad_meta_sig_inf = np.isinf(meta_sig)
      meta_sig[bad_meta_sig_nan] = 0
      meta_sig[bad_meta_sig_inf] = 0
      polydata.GetCellData().AddArray(asVtkArray(meta_sig, "metabolic", vtkFloatArray))
  if 'S_tot' in locals():
    if np.sum(np.isnan(S_tot))==0 and np.sum(np.isinf(S_tot))==0:
      polydata.GetCellData().AddArray(asVtkArray(S_tot, "S_tot", vtkFloatArray))
    else:
      print("warning: bad S_tot at ")
      print( np.where(np.isnan(S_tot) ==1)[0])
      print( np.where(np.isinf(S_tot) ==1)[0])
      bad_S_tot_nan = np.isnan(S_tot)
      bad_S_tot_inf = np.isinf(S_tot)
      S_tot[bad_S_tot_nan]=0
      S_tot[bad_S_tot_inf]=0
      print("len(S_tot)")
      print(len(S_tot))
      polydata.GetCellData().AddArray(asVtkArray(S_tot, "S_tot", vtkFloatArray))
  if 'po2_node' in graph.keys():
    po2_node = graph['po2_node']    
    polydata.GetPointData().AddArray(asVtkArray(po2_node, "point_vessel_po2", vtkFloatArray))
  if 'po2_vessel' in graph.keys():
    po2_vessels = graph['po2_vessel']
    polydata.GetCellData().AddArray(asVtkArray(po2_vessels, "po2_vessel", vtkFloatArray))
  if 0:
    #vesselgroup = h5file[posixpath.join(grpname, 'vessels')]
    vesselgroup = h5file['/vessels']

    #everything else
    ignore = ['lattice_pos','roots','node_a_index','node_b_index','radius']
    g = vesselgroup['edges'].values() + vesselgroup['nodes'].values()
    for item in g:
      short_name = item.name.split('/')[-1]
      if short_name in ignore: continue
      n = item.shape
      if len(n) <> 1:
        print 'ignoring item %s because if unsupported shape %s' % (item.name, n)
        continue
      n = n[0]
      if n == num_verts:
        print 'vert data %s [%s,%s]' % (item.name, item.dtype, str(item.shape))
        polydata.GetPointData().AddArray(asVtkArray(np.asarray(item), short_name))
      elif n == num_edges:
        print 'edge data %s [%s,%s]' % (item.name, item.dtype, str(item.shape))
        polydata.GetCellData().AddArray(asVtkArray(np.asarray(item), short_name))
      else:
        print 'ignoring item %s because array size does not match' % (item.name)
  return polydata



# now done by filter graph
#def removeUncirculatedVessels(dataset):
#  flagarr = fromVtkArray(dataset.GetCellData().GetArray("flags"))
#  ids = [ i for i, f in enumerate(flagarr) if (f & krebs.CIRCULATED) ]
#  return copyPolyDataCells(dataset, ids)

def writeFields_(graph, options):
  fn = str(graph.from_fn)
  f = h5py.File(fn, 'r')
  fn, _ = myutils.splitH5PathsFromFilename(fn)
  #search_groups.append('po2/adaption/vessels_after_adaption')
  #search_groups = ['po2/adaption/vessels_after_adaption/po2field']
  search_groups = []
  search_groups.append(str(options.writeFields))
  e = extractVtkFields.Extractor(f, search_groups, recursive = True, lattice_path="field_ld") 
  print 'found field datasets:'
  pprint.pprint(e.getDatasetPaths())
  filename_append = 'fields_%s' % str(os.path.split(search_groups[0])[-1])
  e.write(options.outfn % filename_append)
  del e

def writeVessels_(graph, options):  
  polydata = ConvertMyHdfVesselsToVTKPolydata(graph, False, goodArguments);
  writer = vtkPolyDataWriter()
  print("use vtkVersion: %s" % vtkVersion.GetVTKVersion())
  if(int(vtkVersion.GetVTKVersion()[0])>5):
    writer.SetInputData(polydata)
  else:
    writer.SetInput(polydata)
  
  writer.SetFileName(options.outfn % 'vessels')  
  writer.Write()

def writeCells_(graph, options):
  fn = str(graph.from_fn)
  f = h5py.File(fn, 'r')
  #fn, _ = myutils.splitH5PathsFromFilename(fn)
  #search_groups.append('po2/adaption/vessels_after_adaption')
  #search_groups = ['po2/adaption/vessels_after_adaption/po2field']
  #search_groups = []
  #search_groups.append(str(options.writeFields))
  #search_groups.append("out0005")
  pos = np.asarray(f[str(options.grp_pattern)+'/cells/cell_center_pos'])
  
  pts = vtkPoints()
  pts.SetNumberOfPoints(len(pos))
  for i,(x,y,z) in enumerate(pos):
    pts.SetPoint(i, (float(x), float(y), float(z)))
    
  polydata = vtkPolyData()
  polydata.SetPoints(pts)
  
  cellGroup = f[str(options.grp_pattern)+'/cells']
  for aKey in cellGroup.keys():
    npReadOut = np.asarray(cellGroup[aKey])
    polydata.GetPointData().AddArray(asVtkArray(npReadOut, aKey, vtkFloatArray))
  
  writer = vtkPolyDataWriter()
  print("use vtkVersion: %s" % vtkVersion.GetVTKVersion())
  if(int(vtkVersion.GetVTKVersion()[0])>5):
    writer.SetInputData(polydata)
  else:
    writer.SetInput(polydata)
  
  writer.SetFileName(options.outfn % 'cells')  
  writer.Write()

def hdftumor2vtk(graph, options ):
  if True:# could be like write cells
    writeCells_(graph, options)
  if options.writeFields:
    writeFields_(graph, options)
  if options.writeVessels:
    writeVessels_(graph, options)


if __name__ == '__main__':

  import argparse
  parser = argparse.ArgumentParser(description='Export hdf data to ParaView.')  
  parser.add_argument('vesselFileNames', nargs='*', type=argparse.FileType('r'), default=sys.stdin)
  parser.add_argument('grp_pattern')
  parser.add_argument("-d","--data", dest="datalist", help="which data (pressure, flow, shearforce, hematocrit) as comma separated list", default='pressure', action="store")
  parser.add_argument("-f","--filter-uncirculated", dest="filteruncirculated", help="filter uncirculated vessels", default=False, action="store_true")
  parser.add_argument("--filter-radius-high-pass", dest="filterradiushighpass", type=float, default=-1.0)
  parser.add_argument("--no-overlay", dest="overlay", default = True, action="store_false")
  parser.add_argument("--dpi", dest="dpi", default=None, action="store")
  parser.add_argument("--format", dest="format", default=None, action="store")
  parser.add_argument("-n", dest="out_nl", default = False, action="store") 
  parser.add_argument("--outFilename", dest="outfn", default= None, type=str)
  '''this option come from the tumor side'''
  parser.add_argument("--writeVessels", help="when doing the tumor, export vesesls", default=True, action="store_true")  
  parser.add_argument("--writeFields", help="when doing the tumor, export Fields", default=None)  
  goodArguments, otherArguments = parser.parse_known_args()
  #create filename due to former standards
  filenames=[]
  for fn in goodArguments.vesselFileNames:
    filenames.append(fn.name)
    
  datalist = map(lambda s: s, map(str.strip, goodArguments.datalist.split(',')))
  pattern   = goodArguments.grp_pattern
  
  '''check input'''
  try:
    dirs = set()
    for fn in goodArguments.vesselFileNames:
      if not os.path.isfile(fn.name):
        raise AssertionError('The file %s is not present!'%fn)
      with h5py.File(fn.name, 'r') as f:
        d = myutils.walkh5(f, goodArguments.grp_pattern)
        if not len(d)>0:
          raise AssertionError('pattern "%s" not found in "%s"!' % (grp_pattern, fn))
        else:
          dirs = set.union(dirs,d)
  except Exception, e:
    print e.message
    sys.exit(-1)

  for fn in filenames:
    fn, _ = myutils.splitH5PathsFromFilename(fn)
    f = h5py.File(fn, 'r')
    dirs = myutils.walkh5(f['/'], pattern)
    if goodArguments.outfn:
      goodArguments.outfn = goodArguments.outfn + '_%s.vtk'
    else:
      goodArguments.outfn = outfn = "%s-%%s.vtk" % (os.path.splitext(os.path.basename(fn))[0]+('_%s'%pattern))
    print("you chose: %s as outfilename" % goodArguments.outfn)
    for d in dirs:
      if 'vessels' in d and 'po2' not in d:
        vesselgroup = f[join('/',d)]['.']
        new = False
        if new:
          graph = krebsutils.read_vessels_from_hdf(vesselgroup, ['position', 'radius', 'hematocrit', 'pressure', 'flow', 'flags','shearforce','nodeflags','edge_boundary'] + datalist, return_graph=True)
        else:
          graph = krebsutils.read_vessels_from_hdf(vesselgroup, ['position', 'radius', 'hematocrit', 'pressure', 'flow', 'flags','shearforce'] + datalist, return_graph=True)
        if goodArguments.filteruncirculated:
          graph = graph.get_filtered(edge_indices = myutils.bbitwise_and(graph['flags'], krebsutils.CIRCULATED))      
        writeVessels_(graph, goodArguments)

      elif 'out' in d:
        vesselgroup = f[join('/',d+'/vessels')]['.']
        if '/parameters' in f:
          useConstO2 = f['/parameters'].attrs['useConstO2']
          # comes as string
          useConstO2 = bool(useConstO2)
        else:
          useConstO2 = False
          
        if useConstO2:
          ''' the po2_node will not be present in this case '''
          graph = krebsutils.read_vessels_from_hdf(vesselgroup, ['position', 'radius', 'hematocrit', 'pressure', 'flow', 'flags','shearforce'] + datalist, return_graph=True)
        else:
          graph = krebsutils.read_vessels_from_hdf(vesselgroup, ['po2_node','position', 'radius', 'hematocrit', 'pressure', 'flow', 'flags','shearforce'] + datalist, return_graph=True)
        if goodArguments.filteruncirculated:
          graph = graph.get_filtered(edge_indices = myutils.bbitwise_and(graph['flags'], krebsutils.CIRCULATED))
        hdftumor2vtk(graph, goodArguments)
      elif 'po2' in d:
        #vesselgroup = f[join('/',d+'/vessels')]['.']
        if 'adaption' in d:
          vesselgroup = f['/recomputed_flow/adaption/vessels_after_adaption']
          new = False
          if new:
            graph = krebsutils.read_vessels_from_hdf(vesselgroup, ['position', 'radius', 'hematocrit', 'pressure', 'flow', 'flags','shearforce','nodeflags','edge_boundary'] + datalist, return_graph=True)
          else:
            graph = krebsutils.read_vessels_from_hdf(vesselgroup, ['position', 'radius', 'hematocrit', 'pressure', 'flow', 'flags','shearforce'] + datalist, return_graph=True)
          if goodArguments.filteruncirculated:
            graph = graph.get_filtered(edge_indices = myutils.bbitwise_and(graph['flags'], krebsutils.CIRCULATED))
          #amazing, out of the box this works for the o2 simulation as well.
          hdftumor2vtk(graph, goodArguments)
          print('bla2')
        else:
          vesselgroup = f['/recomputed_flow']
          new = False
          if new:
            graph = krebsutils.read_vessels_from_hdf(vesselgroup, ['position', 'radius', 'hematocrit', 'pressure', 'flow', 'flags','shearforce','nodeflags','edge_boundary'] + datalist, return_graph=True)
          else:
            graph = krebsutils.read_vessels_from_hdf(vesselgroup, ['position', 'radius', 'hematocrit', 'pressure', 'flow', 'flags','shearforce'] + datalist, return_graph=True)
          if goodArguments.filteruncirculated:
            graph = graph.get_filtered(edge_indices = myutils.bbitwise_and(graph['flags'], krebsutils.CIRCULATED))
          #amazing, out of the box this works for the o2 simulation as well.
          hdftumor2vtk(graph, goodArguments)
