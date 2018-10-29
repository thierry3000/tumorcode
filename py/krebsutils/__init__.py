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
import os, sys
import numpy as np
import h5py
import extensions # for efficient asarray with h5py
import posixpath
import sys

sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../../lib'))
'''overcomes some serious mpi issues!!!
https://github.com/baidu-research/tensorflow-allreduce/issues/4
http://users.open-mpi.narkive.com/BxD0j82y/ompi-users-problem-with-using-mpi-in-a-python-extension

NOTE: this issues should not be present in a properly configured system
'''
import platform
theUnameList = platform.uname()
isUbuntu=False
for entry in theUnameList:
    if 'Ubuntu' in entry:
        isUbuntu=True
if isUbuntu:        
    import ctypes
    ctypes.CDLL("libmpi.so", mode=ctypes.RTLD_GLOBAL)

# leaks a bit of memory each time it is imported!
from scipy.ndimage.interpolation import geometric_transform

if sys.flags.debug:
  libkrebs = __import__('libkrebs_d', globals(), locals())
else:
  libkrebs = __import__('libkrebs_', globals(), locals())

''' this call import often used cpp elements

SetupFieldLattice: 
  arguments: wbbox, dim, spacing, safety_spacing
    wbbox:  the world box for the lattice
    dim:    
    spacing:lattice spacing
    
'''
#HACK2018
#imports_ = [ f.strip() for f in
#    '\
#    LatticeData, \
#    read_lattice_data_from_hdf_by_filename, \
#    write_lattice_data_to_hdf_by_filename, \
#    export_network_for_povray, \
#    ClipShape, \
#    povray_clip_object_str, \
#    make_position_field, \
#    calcBulkTissueSourceTerm, \
#    make_vessel_volume_fraction_field, \
#    calc_vessel_boxcounts, \
#    set_num_threads, \
#    run_vesselgen, \
#    vesselgen_generate_grid, \
#    vesselgen_generate_single, \
#    vesselgen_generate_symmetric, \
#    GetHealthyVesselWallThickness, \
#    CalcRelativeViscosity, \
#    CalcFahraeusEffect, \
#    CalcIntervascularInterpolationField_, \
#    SetupFieldLattice, \
#    PressureRadiusRelation, \
#    SumIsoSurfaceIntersectionWithVessels_, \
#    get_Murray2, \
#    get_Murray_scale, \
#    CalcViscosities, \
#    CalcConductivities'.split(',')
#]

if libkrebs.is_vbl_used():
    print("VBL is used!")
    imports_ = [ f.strip() for f in
        '\
        LatticeData, \
        read_lattice_data_from_hdf_by_filename, \
        write_lattice_data_to_hdf_by_filename, \
        export_network_for_povray, \
        export_VBL_Cells_for_povray, \
        ClipShape, \
        povray_clip_object_str, \
        make_position_field, \
        calcBulkTissueSourceTerm, \
        make_vessel_volume_fraction_field, \
        calc_vessel_boxcounts, \
        run_vesselgen, \
        vesselgen_generate_grid, \
        vesselgen_generate_single, \
        vesselgen_generate_symmetric, \
        GetHealthyVesselWallThickness, \
        CalcRelativeViscosity, \
        CalcFahraeusEffect, \
        CalcIntervascularInterpolationField_, \
        SetupFieldLattice, \
        PressureRadiusRelation, \
        SumIsoSurfaceIntersectionWithVessels_, \
        get_Murray2, \
        get_Murray_scale, \
        CalcViscosities, \
        CalcConductivities'.split(',')
    ]
else:
    print("VBL is unused!")
    imports_ = [ f.strip() for f in
        '\
        LatticeData, \
        read_lattice_data_from_hdf_by_filename, \
        write_lattice_data_to_hdf_by_filename, \
        export_network_for_povray, \
        ClipShape, \
        povray_clip_object_str, \
        make_position_field, \
        calcBulkTissueSourceTerm, \
        make_vessel_volume_fraction_field, \
        calc_vessel_boxcounts, \
        run_vesselgen, \
        vesselgen_generate_grid, \
        vesselgen_generate_single, \
        vesselgen_generate_symmetric, \
        GetHealthyVesselWallThickness, \
        CalcRelativeViscosity, \
        CalcFahraeusEffect, \
        CalcIntervascularInterpolationField_, \
        SetupFieldLattice, \
        PressureRadiusRelation, \
        SumIsoSurfaceIntersectionWithVessels_, \
        get_Murray2, \
        get_Murray_scale, \
        CalcViscosities, \
        CalcConductivities'.split(',')
    ]

# CalcRelativeViscosityByTable, \
# load functions from libkrebs_ into local namespace
locals().update( (f,getattr(libkrebs, f)) for f in imports_)

    #DetailedO2Parameters, \
    #AllocateDetailedO2ParametersFromDict, \

imports_ = [ f.strip() for f in 'testCalcOxy, testCalcOxy2, testCalcOxy3, testCalcOxy4, testCalcOxy5'.split(',') ]
imports_ = [ f for f in imports_ if hasattr(libkrebs, f) ]
locals().update( (f,getattr(libkrebs, f)) for f in imports_)

''' arguments are: 
      vesselgroup
      return_flags          --- decide if  vessel flags are returned
      bloodflowparams       --- need for calcflow
      simple                --- if simple, no hematocrit is returned
      storeCalculationInHDF --- add a "recomputed" subfolder to the vesselfile
    return values:
      a python list object with elements:
    (pressure, flow, force, **hematocrit**, **flags**)
      ** values are optional
'''
calc_vessel_hydrodynamics_Ccode = libkrebs.calc_vessel_hydrodynamics
#read_vessel_positions_from_hdf_ = libkrebs.read_vessel_positions_from_hdf
read_vessel_positions_from_hdf_by_filename = libkrebs.read_vessel_positions_from_hdf_by_filename
#read_vessel_positions_from_hdf_edges_ = libkrebs.read_vessel_positions_from_hdf_edges
#read_vessel_positions_from_hdf_world_ = libkrebs.read_vessel_positions_from_hdf_world
flood_fill_ = libkrebs.flood_fill
distancemap_ = libkrebs.distancemap
curvature_ = libkrebs.calcCurvature
fill_with_smooth_delta_ = libkrebs.fill_with_smooth_delta

'''if flag av is false:
  returns asked property inside and outside tumor in 2 lists.
   if flag av is True:
  returns asked property inside tumor arteries and vein,
  and outside for arteries and veins.
  4 list:
    inside:
      arteries veins
    outside
      arteries veins
'''
#calculate_within_fake_tumor_lattice_based = libkrebs.calculate_within_fake_tumor_lattice_based

#globals
typelist = 'typeA typeB typeC typeD typeE typeF typeG typeH typeI'.split()

#----------------------------------------------------------------------------------#
#  utility routines
#----------------------------------------------------------------------------------#
#def get_Murray(vesselgrp):
#  #return get_Murray2_p(vesselgrp, alpha)
#  return get_Murray2(vesselgrp)

def get_full_tumor_executable_path(name):
  from os.path import join, abspath, normpath, dirname
  '''Get the full path using only the file name. It is guesstimated based on conventions.
     See get_tumor_executable_dir
     
     Determine the path where tumor sim executables reside. This depends on whether
     we run from the build dir or from an install dir and also on whether we have a
     debug build or not. The build dir follows the convention that release build dir
     is named "buildopt" and the debug build dir is named "build" and both reside in the
     tumorcode base directory'''
     
  path = normpath(os.environ['TC_INSTALL']) # probable base directory of the project
  path1 = join(path, 'bin') # install subdirectory
  path2 = join(path, 'lib') # put here when install PYTHON_LIBS_OUTPUT_IN_SOURCE_DIR is true
  for p in [path1, path2]:
    exe = abspath(join(p, name))
    if os.path.isfile(exe):
      print 'found executable in %s' % exe
      return exe
  raise RuntimeError('tried to find tumor simulation executable %s but failed' % name)


#----------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------#

def vector_lengths(a):
  """
    take n x m array and compute the 2-norm of each row
  """
  a = a.transpose()
  r = np.zeros_like(a[0])
  for x in a:
    r += x*x
  r = np.sqrt(r)
  return r


asarray = np.asarray


def errormsg_(cond, msg, *args):
  if not cond:
    raise RuntimeError(msg % args)


def get_krebssubroutine_by_type_(name, dtype):
  """
    Give a base name and numpy dtype to get the corresponding function object of the c implementation.
    The c implementation follows a naming convenction so it can be found.
  """
  postfix = {
    np.dtype(np.float32) : '_float',
    np.dtype(np.float64) : '_double',
    np.dtype(np.int32)   : '_int',
    np.dtype(np.uint32)   : '_uint',
    np.dtype(np.int8)    : '_char',
    np.dtype(np.uint8)   : '_uchar',
  }
  return getattr(libkrebs, name+postfix[dtype])



#----------------------------------------------------------------------------------#
#  lattice data
#----------------------------------------------------------------------------------#

def LatticeDataGetWorldBox(ld):
  return ld.GetWorldBoxRavel()
    
def LatticeDataGetBox(ld):
  return ld.GetBoxRavel()

def LatticeDataGetWorldSize(ld):
  """
    world size as 3 component float array
  """
  bbox = LatticeDataGetWorldBox(ld).reshape(3,2)
  return bbox[:,1]-bbox[:,0]


def LatticeDataRescale(ld, factor):
  """
    scale around the origin by factor
  """
  ld.SetScale(ld.scale*factor)
  x0, _, y0, _, z0, _ = ld.GetWorldBoxRavel()
  ld.SetOriginPosition(factor * np.asarray((x0,y0,z0)))

def LatticeDataGetShape(ld):
  """
    lateral size of the box
  """
  b = ld.GetBoxRavel().reshape(3,2)
  return b[:,1]-b[:,0]+1

def LatticeDataEqual(ld, other):
  return all(ld.GetBox() == other.GetBox()) and (ld.GetScale() == other.GetScale()) # TODO: cell centering

def LatticeDataGetCentered(ld):
  #org[i] = -0.5*size[i]+(origin[i]-wbb.min[i]);
  wbox = ld.worldBox.reshape(3,2)
  originOffset = ld.GetOriginPosition()
  compute = lambda ((l, r), o): -0.5*(r-l) + o - l
  originOffset = map(compute, zip(wbox, originOffset))
  ld = ld.Copy()
  ld.SetOriginPosition(originOffset)
  return ld
  
# probably not a good idea since world coordinates are not compared
#def LatticeDataHash(ld):
#  return hash((ld.scale, tuple(ld.GetBox()))) # TODO: cell centering
      
# in order to get get the bounding box in the format
# [(xmin, xmax)
#  (ymin, ymax)
#  (zmin, zmax)]
# do ld.box.reshape(3,2)

LatticeData.GetWorldSize = LatticeDataGetWorldSize
LatticeData.GetWorldBox  = LatticeDataGetWorldBox
LatticeData.GetBox       = LatticeDataGetBox
LatticeData.Rescale = LatticeDataRescale
LatticeData.shape = property(LatticeDataGetShape)
LatticeData.scale = property(LatticeData.GetScale)
LatticeData.box   = property(LatticeData.GetBox)
LatticeData.worldBox = property(LatticeData.GetWorldBox)
LatticeData.worldSize = property(LatticeData.GetWorldSize)
LatticeData.Copy  = lambda self : LatticeData(self)
LatticeData.GetCentered = LatticeDataGetCentered
LatticeData.__eq__ = LatticeDataEqual
#LatticeData.__hash__ = LatticeDataHash

LatticeDataQuad3d = lambda *args : LatticeData('QUAD3D',*args)
LatticeDataFCC    = lambda *args : LatticeData('FCC',*args)


#----------------------------------------------------------------------------------#
#  functions to read vessel graph data
#----------------------------------------------------------------------------------#

def filter_graph_byedge2( edges, edge_data, node_data, indices, return_indices=False):
    if isinstance(indices, tuple) and len(indices)==1:
      indices = indices[0]
    assert len(indices.shape)==1
    tmpedges = edges[indices,...]
    n = np.amax(edges)+1
    u = np.zeros((n,),dtype=np.int)
    np.put(u,tmpedges.ravel(),1)
    map, = np.where(u)
    invmap = -1*np.ones((n,),dtype=edges.dtype)
    np.put(invmap,map,np.arange(map.shape[0],dtype=edges.dtype))
    newedges = np.asarray(np.take(invmap,tmpedges), dtype=edges.dtype)
    new_edge_data = tuple( q[indices,...] for q in edge_data )
    new_node_data = tuple( q[map,...] for q in node_data )
    ''' difficult with position which is dict with 3 entries'''
    #new_node_data = node_data
    
    ret = (newedges, new_edge_data, new_node_data)
    if return_indices:
        ret = ret+(indices, map, invmap) # map  is array to obtain old node index from new index by: oldIndex=map[newIndex]
    return ret


def read_graph_(grp, *prop_names):
    """
    read graph structure:
        edges\
            edge properties: array
            ...
        nodes\
            node_a_index: 1d int array
            node_b_index: 1d int array
            node properties: array
            ...
    return a tuple with first the edge->node indices as nx2 int array, then properties in the order of input names
    """
    ge = grp['edges']
    gn = grp['nodes']
    node_a = asarray(ge['node_a_index'], dtype=np.int32)
    node_b = asarray(ge['node_b_index'], dtype=np.int32)
    nodes = np.column_stack((node_a,node_b))
    g = Graph(edgelist = nodes)
    not_found = set()
    for prop_name in prop_names:
        if prop_name in gn:
            if prop_name == 'pressure':
              g.nodes[prop_name] = 7.5*asarray(gn[prop_name]) # to mmHg
            else:
              g.nodes[prop_name] = asarray(gn[prop_name])
        elif prop_name in ge:
            g.edges[prop_name] = asarray(ge[prop_name])
        else:
            not_found.add(prop_name)
    return g, not_found


#def filter_graph_byedge( edges, edge_data, node_data, filter, names=None, return_indices=False ):
#    if names:
#        tmp = list(edge_data.__iter__()) + list(node_data.__iter__())
#        tmp.reverse()
#        kw = dict()
#        for name in names:
#            if name:
#                kw[name] = tmp.pop()
#            else:
#                tmp.pop()
#        indices, = np.where( filter(edges,kw) )
#        del tmp,kw
#    else:
#        indices, = np.where( filter(edges,*edge_data) )
#    return filter_graph_byedge2(edges, edge_data, node_data, indices, return_indices)


class Graph(object):
  def __init__(self, edgelist = []):
    self.nodes = {}
    self.edges = {}
    self.edgelist = edgelist
    self.roots = None
    self.from_fn = None
  def __getitem__(self, name):
    if name == 'edgelist':
      item = self.edgelist
    try:
      item = self.nodes[name]
    except KeyError:
      try:
        item = self.edges[name]
      except KeyError:
        raise KeyError('cannot find %s in neither in nodes or edges' % name)
    return item
  def __contains__(self, name):
    return name == 'edgelist' or name in self.nodes or name in self.edges
  def keys(self):
    return self.nodes.keys()+self.edges.keys()
  def values(self):
    return self.nodes.values()+self.edges.values()
  def items(self):
    return self.nodes.items()+self.edges.items()
  def get_filtered(self, edge_indices=None):
    #assert not self.roots, 'implement filtering of roots'
    if self.roots is not None:
      print 'WARNING: test Graph.get_filtered for tree roots (untested)!'
    try:
      tmp = np.asarray(edge_indices)
    except:
      pass
    else:
      if tmp.dtype == np.bool and tmp.shape[0] == self.edgelist.shape[0]:
        edge_indices = np.nonzero(tmp)[0]
    el = self.edgelist
    nnames, nprop = zip(*self.nodes.items()) if self.nodes else ([], [])
    enames, eprop = zip(*self.edges.items()) if self.edges else ([], [])
    new_el, new_eprop, new_nprop, _, toOldNode, toNewNode = filter_graph_byedge2(el, eprop, nprop, edge_indices, True)
    if self.roots is not None:
      self.roots = np.take(toNewNode, self.roots)
    g = Graph(edgelist = new_el)
    g.from_fn = self.from_fn
    g.nodes = dict(zip(nnames,new_nprop))
    g.edges = dict(zip(enames,new_eprop))
    return g
  def __str__(self):
    return '<Graph N: %i %s, E: %i %s>' % (self.edgelist.max(), self.nodes.keys(), len(self.edgelist), self.edges.keys())


def find_lattice_group_(vesselgroup):
  """
    since i am not consistent where to store the lattice data for the
    vessel network i need this helper to find it.

    returns - lattice data group
  """
  if 'lattice' in vesselgroup:
    ldgroup = vesselgroup['lattice']
  else:
    ldgroup = vesselgroup.file['/lattice']
  return ldgroup


class CannotComputeException(Exception):
  pass


def vessels_require_(vesselgroup, g, name):
  """
    internal routine to compute certain often used quantities and insert them in the graph structure
  """
  assert isinstance(name, str)

  if name in ['pressure', 'flow', 'shearforce','nodeflags']:
    if all(('pressure' in g.nodes, 'flow' in g.edges, 'shearforce' in g.edges, 'nodeflags' in g.nodes)): return
  if name in g.edges or name in g.nodes: return

  e = g.edgelist
  if 'position' == name:
    
    if "CLASS" in vesselgroup.attrs:
      # c++ site now manages this
      # I try to remove the this dependency, since h5py not always uses the same
      # library as c++
      #pos = read_vessel_positions_from_hdf_(vesselgroup).transpose()
      fn=str(vesselgroup.file.filename)
      path = str(vesselgroup.name)
      pos_x, pos_y, pos_z = read_vessel_positions_from_hdf_by_filename(fn, path)
      pos = np.asarray([pos_x, pos_y, pos_z])
      print("before:")
      print("max x: %f" % np.max(pos_x))
      print("min x: %f" % np.min(pos_x))
      print("max y: %f" % np.max(pos_y))
      print("min y: %f" % np.min(pos_y))
      print("max z: %f" % np.max(pos_z))
      print("min z: %f" % np.min(pos_z))
      print(pos.shape)
      print(pos)
      pos = pos.transpose()
      print("after:")
      print(pos.shape)
      print(pos)
      
      print("max x: %f" % np.max(pos[:,0]))
      print("min x: %f" % np.min(pos[:,0]))
      print("max y: %f" % np.max(pos[:,1]))
      print("min y: %f" % np.min(pos[:,1]))
      print("max z: %f" % np.max(pos[:,2]))
      print("min z: %f" % np.min(pos[:,2]))
      
      g.nodes['position'] = pos  
    else:
      print("WARNING")
      print("We assume this is a non lattice based structure!")
      #pos = read_vessel_positions_from_hdf_(vesselgroup).transpose()
      g.nodes['position'] = np.asarray(vesselgroup['nodes/world_pos'])
  elif 'edge_position' == name:
    
    if "CLASS" in vesselgroup.attrs:
      # c++ site now manages this
      pos = read_vessel_positions_from_hdf_edges_(vesselgroup).transpose()
      g.edges['edge_position'] = pos  
    else:
      print("unknown structure")

  elif name in ['pressure', 'flow', 'shearforce']:
      if "CLASS" in vesselgroup.attrs:
        dd = calc_vessel_hydrodynamics(vesselgroup)
        g.nodes['pressure'] = dd[0]
        g.edges['flow'] = dd[1]
        g.edges['shearforce'] = dd[2]
#          if vesselgroup.attrs['CLASS'] == 'LATTICE':
#            dd = calc_vessel_hydrodynamics(vesselgroup)
#            g.nodes['pressure'] = dd[0]
#            g.edges['flow'] = dd[1]
#            g.edges['shearforce'] = dd[2]
#          if vesselgroup.attrs['CLASS'] == 'REALWORLD':
#            dd = calc_vessel_hydrodynamics_world(vesselgroup)
#            g.nodes['pressure'] = dd[0]
#            g.edges['flow'] = dd[1]
#            g.edges['shearforce'] = dd[2]
      else:
          print("unknown structure")

  elif 'pressure_gradient' == name:
    vessels_require_(vesselgroup, g, 'pressure')
    vessels_require_(vesselgroup, g, 'length')
    dp = g.nodes['pressure']
    dp = np.abs(dp[e[:,1] - dp[e[:,1]]]) / g.edges['length']
    g.edges['pressure_gradient'] = dp

  elif 'hematocrit' == name:
    g.edges['hematocrit'] = 0.45 * np.ones_like(g.edges['radius'])

  elif 'conductivity' == name:
    vessels_require_(vesselgroup, g, 'length')
    vessels_require_(vesselgroup, g, 'hematocrit')
    r = g.edges['radius']
    l = g.edges['length']
    h = g.edges['hematocrit']
    c = calc_vessel_conductivities(r, l, h)
    g.edges['conductivity'] = c

  elif 'edge_pressure' == name:
    vessels_require_(vesselgroup, g, 'pressure')
    p = g.nodes['pressure']
    p = 0.5*(p[e[:,0]] + p[e[:,1]])
    g.edges['pressure'] = p
  elif 'edge_boundary' == name:
    isnodeboundary = np.asarray(np.bitwise_and(g.nodes['nodeflags'],BOUNDARY) > 0, dtype=np.int32)    
    g.edges['edge_boundary'] = np.bitwise_or(isnodeboundary[e[:,0]],isnodeboundary[e[:,1]])
  elif 'length' == name:
    vessels_require_(vesselgroup, g, 'position')
    pos = g.nodes['position']
    dp = pos[e[:,0],:] - pos[e[:,1],:]
    l = vector_lengths(dp)
    g.edges['length'] = l
  elif 'flags' == name:
    vessels_require_(vesselgroup, g, 'flags')
  elif 'po2_node' == name:
    #find detailed o2 group
    if 'recomputed' in vesselgroup.name:
      ''' usualy the case for deatailed o2 calculation '''
      ''' vesselgroup.parent.parent is po2 group'''
      po2_vessels = np.asarray(vesselgroup.parent.parent['po2/vessels/po2vessels'])
    else:
      ''' if we have the mts tumor things are slightly different'''
      po2_vessels = np.asarray(vesselgroup.parent['po2/po2vessels'])
    num_verts = len(g['position'])
    node_po2 = np.zeros(num_verts)
    ''' counts the vertices, they might apear several time
        so we take the average
    '''
    node_n = np.zeros(num_verts)
    for ((a,b),(po2_a,po2_b)) in zip(e, po2_vessels):
      node_po2[a] += po2_a
      node_po2[b] += po2_b
      node_n[a] += 1.
      node_n[b] += 1.
    node_po2 /= np.maximum(1, node_n)
    g.nodes['po2_node'] = node_po2
  
  elif 'po2_vessel' == name:
    #find detailed o2 group  -- works
    if 'recomputed' in vesselgroup.name:
      ''' usualy the case for deatailed o2 calculation '''
      ''' vesselgroup.parent.parent is po2 group'''
      po2_vessels = np.asarray(vesselgroup.parent.parent['po2/vessels/po2vessels'])
    else:
      po2_vessels = np.asarray(vesselgroup.parent['detailedPo2/po2_vessels'])
    po2AtVessel = np.average(po2_vessels,1)
    g.edges['po2_vessel'] = po2AtVessel 
  else:
    raise CannotComputeException('do not know how to compute %s' % name)



def read_vessels_from_hdf(f, prop_names, return_graph=False, return_not_found=False, compute_not_found=True):
  """
    This is the main routine to load a vessel graph and associate data
    prop_names  -    is a list of data names
    f           -   is a hdf5 group where the vessel data is found
    return : a Graph object
  """
  assert isinstance(f, (h5py.File, h5py.Group)), "read_vessels_from_hdf needs a h5py.File or h5py.Group instance"
  vesselgroup = f#f[posixpath.join(grpname,'vessels')]
  g, not_found = read_graph_(vesselgroup, *prop_names)

  roots = np.asarray(vesselgroup['nodes/roots'])
  g.roots = roots  

  if compute_not_found:
    for x in not_found.copy():
#      try:
        vessels_require_(vesselgroup, g, x)
#      except CannotComputeException:
#        pass
#      else:
        not_found.discard(x)

  if return_graph:
    g.from_fn = f.file.filename
    r = g
  else:
    r = [g.edgelist] + [ (g[name] if name in g else None) for name in prop_names ]
  if return_not_found:
    r = r, not_found
  return r


def read_vesselgraph(vesselgroup, datanames, return_not_found = False):
    return read_vessels_from_hdf(vesselgroup, datanames, return_graph = True, return_not_found = return_not_found)

#this assures backward compatibility
def calc_vessel_hydrodynamics(vesselgroup, calc_hematocrit=False, return_flags=False, override_hematocrit = None, bloodflowparams = dict(),storeCalculationInHDF=False):
    if 'CLASS' in vesselgroup.attrs:
      return calc_vessel_hydrodynamics_(vesselgroup, calc_hematocrit, return_flags, override_hematocrit, bloodflowparams, storeCalculationInHDF)
    else:
        os.error('Unknown vessel structure!')

def calc_vessel_hydrodynamics_(vesselgroup, calc_hematocrit, return_flags, override_hematocrit, bloodflowparams, storeCalculationInHDF):
  """
    calls the c++ routines to load a vessel network into a VesselList3d object,
    compute flowrates etc with calcflow, and return (pressure, flow, shearforce, hematocrit, flags)
    as numpy arrays.
  """
  # this is very messy now due to the new bloodflowparameters
  simple = True # keep hematocrit that is already set in teh vesselgraph
  if calc_hematocrit:
    simple = False
  ''' do not use the override_hematorit flag anymore'''
#  if override_hematocrit is not None:
#    bloodflowparams['inletHematocrit'] = float(override_hematocrit)
#  if calc_hematocrit:
#    bloodflowparams['includePhaseSeparationEffect'] = True
#  if ('inletHematocrit' in bloodflowparams) or ('includePhaseSeparationEffect' in bloodflowparams):
#    simple = False # assign new hematocrit
  if not bool(bloodflowparams):
    print('Warning: bloodflowparams are empty!!!')
    print('Using c++ default falue')
  #usage:
  #const py::object &vess_grp_obj ,bool return_flags, const BloodFlowParameters &bfparams, bool simple
  #return calc_vessel_hydrodynamics_Ccode(vesselgroup, return_flags, bloodflowparams, simple, storeCalculationInHDF)
  fn=str(vesselgroup.file.filename)
  vessel_path=str(vesselgroup.name)
  return calc_vessel_hydrodynamics_Ccode(fn, vessel_path, return_flags, bloodflowparams, simple, storeCalculationInHDF)

def calc_vessel_conductivities(rad, length, hema, bloodflowparams = dict()):
  rad    = np.ascontiguousarray(rad, dtype = np.float64)
  length = np.ascontiguousarray(length, dtype = np.float64)
  hema   = np.ascontiguousarray(hema, dtype = np.float64)
  visc = CalcViscosities(rad, hema, bloodflowparams)
  return CalcConductivities(rad, length, visc)



def read_vessel_positions_from_hdf(vesselgroup):
  """
    returns a n x 3 numpy array with node positions
  """
  if vesselgroup.attrs['CLASS']=='GRAPH':
      return read_vessel_positions_from_hdf_(vesselgroup)
  if vesselgroup.attrs['CLASS']=='REALWORLD':
      return read_vessel_positions_from_hdf_world_(vesselgroup)


def edge_to_node_property(num_nodes, edges, prop, combinefunc):
  """
    convert data from edges to nodes
  """
  cf = {
    'max' : 1,
    'min' : 2,
    'and' : 3,
    'or'  : 4,
    'avg' : 6
  }[combinefunc]
  func = get_krebssubroutine_by_type_('edge_to_node_property', prop.dtype)
  component_shape = prop.shape[1:]
  try:
    prop = prop.reshape((len(prop),-1))
  except ValueError:
    print('No elements of considered property is present')
  res = func(num_nodes, edges, prop, cf)
  res = res.reshape((num_nodes,)+component_shape)
  return res


BOUNDARY = 1
CIRCULATED = 2
CONNECTED = 4
ARTERY = 8
VEIN = 16
CAPILLARY = 32
WITHIN_TUMOR = 64

# flow boundary conditions
# numbers must match c++ side in calcflow.h!!!!!
FLOWBC_PIN = 1
FLOWBC_CURRENT = 2
FLOWBC_RESIST = 3

#----------------------------------------------------------------------------------#
#  here come functions to sample vessels and manipulate gridded data
#----------------------------------------------------------------------------------#
'''
arguments:
  pypos:      coordinate in 3D space, best the nodepoints
  pyedges:    begining and endpoint
  sample_len: distance to put sample points on the vessel
  
return:
  returns an array with the weight 
'''
sample_edges_weights = libkrebs.sample_edges_weights

class VesselSamplingFlags(object):
  DATA_PER_NODE = 1
  DATA_CONST = 2
  DATA_LINEAR = 4

"""
"""
def sample_edges(pos, edges, data, sample_len, mode_flags):
  """ creates sample points on edges
arguments:
  pos - world position of nodes
  edges - n x 2 node indices
  data - data associated with edges or nodes, depending on mode
  sample_len - average distance between samples. The samples are taken in regular intervals if this value is smaller than the length of vessels. Otherwise the sampling becomes random.
  mode - see Mode enum
returns:
  
  """
  edges = np.asarray(edges, dtype=np.int32)
  pos = np.asarray(pos, dtype=np.float32)
  #print 'sample_edges', pos.shape, pos.dtype, edges.shape, edges.dtype, data.shape, data.dtype
  func = get_krebssubroutine_by_type_('sample_edges', data.dtype)
  errormsg_(pos.shape[1] == 3, "pos array shape[1] must be 3")
  errormsg_(edges.shape[1] == 2, "edges array shape[1] must be 2")
  errormsg_(not ((mode_flags&VesselSamplingFlags.DATA_PER_NODE) and len(data) != len(pos)), "pos and data arrays must have the same first dimension size")
  errormsg_(not ((not mode_flags&VesselSamplingFlags.DATA_PER_NODE) and len(data) != len(edges)), "edge and data arrays must have the same first dimension size")
  if (mode_flags&VesselSamplingFlags.DATA_LINEAR) and not (mode_flags&VesselSamplingFlags.DATA_PER_NODE):
    errormsg_(len(data.shape)==3 and data.shape[2]==2, "data array len(shape) must be 3 and shape[2]==2")
  else:
    errormsg_(len(data.shape)<=2, "data array len(shape) must be <= 2")
  res = func(pos, edges, data, sample_len, mode_flags)
  if len(data.shape)>1:
    res = res.reshape((res.shape[0], data.shape[1]))
  else:
    res = res.reshape((res.shape[0],))
  return res


def sample_field(pos, field, ld, linear_interpolation=True, extrapolation_value = None):
  assert type(pos) is np.ndarray and len(pos.shape) == 2 and pos.shape[1]==3
  #assert pos.dtype == field.dtype
  func = get_krebssubroutine_by_type_('sample_field', field.dtype)
  field = np.atleast_3d(field)
  assert all(ld.shape[i] == field.shape[i] for i in xrange(3))
  use_extrapolation_value = extrapolation_value is not None
  if not use_extrapolation_value:
    extrapolation_value = 0
  res = func(pos, field, ld, linear_interpolation, use_extrapolation_value, extrapolation_value)
  return res



def flood_fill(field, pos, dtype=None):
  """
    fills regions of zeros starting from pos.
    returns a new array where filled regions are 1 others are 0.
    return type (dtype) is optional. It is np.uint8 if not given.
  """
  org_field = field
  field = np.asarray(np.atleast_3d(field), dtype=np.uint8)
  pos = np.asarray(pos, dtype=np.int32)
  pos.resize((3,))
  res = flood_fill_(field, pos)
  res = res.reshape(org_field.shape)
  if not dtype is None:
    res = np.asarray(res, dtype=dtype)
  return res



def make_radial_field(ld):
  f = make_position_field(ld)
  fr = np.sqrt(np.sum((np.square(f[...,0]), np.square(f[...,1]), np.square(f[...,2])), axis=0))
  del f
  return fr;



def distancemap(field):
  """
    compute the signed distance from the interface of a region.
    the region is defined by positive values and zeros elsewhere.
    the result is a float field, < 0 within the region and > 0 elsewhere.
  """
  org_field = field
  field = np.asarray(np.atleast_3d(field)>0, dtype=np.uint8)
  res = distancemap_(field)
  res = np.reshape(res, org_field.shape)
  return res



def field_gradient(field, spacing=1.):
  res = []
  org_field = field
  func = get_krebssubroutine_by_type_('diff_field', field.dtype)
  for ax in xrange(len(field.shape)):
    d = func(field, ax, 1./spacing)
#    print d.shape, d.min(), d.max()
    d = np.reshape(d, org_field.shape)
    res.append(d)
  return res


def radial_correlation(field1, field2, distance, bins_per_unit=2, return_raw_data = False, subtract_avg = True, mask = None):
  typesize = { np.dtype(np.float32): 1,
               np.dtype(np.float64): 2}
  if typesize[field1.dtype] > typesize[field2.dtype]:
    field1 = np.asarray(field1, dtype = field2.dtype)
  if typesize[field2.dtype] > typesize[field1.dtype]:
    field2 = np.asarray(field2, dtype = field1.dtype)
  if mask is not None:
    mask = np.asarray(mask, dtype = np.ubyte)
  func = get_krebssubroutine_by_type_('radial_correlation', field1.dtype)
  dim = max(len(field1.shape), len(field2.shape))
  #distance = int(bins[-1]+1)
  distance = (distance,)*dim + (0,)*(3-dim)
  f1 = np.atleast_3d(field1)
  f2 = np.atleast_3d(field2)
  assert len(f1.shape)==3 and len(f2.shape)==3
  r, n, c, s = func(f1, f2, distance, bins_per_unit, subtract_avg, mask)
  if return_raw_data == True:
    return r, n, c, s
  else:
    n = np.maximum(1, n)
    x    = r
    y    = c / n
    yerr = np.sqrt((s/n - np.power(y, 2.))/n)
    return x, y, yerr


def curvature(ld, distancemap, return_force = False, has_ghost_boundary = False):
  '''
    given a distance map, it computes the curvature and surface tension force.
    If has_ghost_boundary is true, then distancemap has an additonal boundarylayer (for b.c.s), else boundary values are extrapolated
  '''
  return curvature_(ld, np.asarray(distancemap, dtype=np.float32), has_ghost_boundary, return_force)


def smooth_delta(distancemap, width):
  a = np.array(distancemap, dtype=np.float32, copy=True)
  a = np.ravel(np.atleast_1d(a))
  fill_with_smooth_delta_(a, width)
  a = a.reshape(distancemap.shape)
  return a


def resample_field(data_src, bbox_src, shape_dst, bbox_dst, order = 2, mode = 'constant', cval=0.):
  #see http://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.interpolation.geometric_transform.htm
  bbox0, bbox1 = np.asarray(bbox_dst), np.asarray(bbox_src)
  o_ = (bbox0[0]-bbox1[0])/(bbox1[1]-bbox1[0])
  s_ = (bbox0[1]-bbox0[0])/(bbox1[1]-bbox1[0])
  s_ /= np.asarray(shape_dst)-1
  o_ *= np.asarray(data_src.shape)-1
  s_ *= np.asarray(data_src.shape)-1
#  def trafo(c_):
#    q = tuple(c*s+o for (c,s,o) in zip(c_, s_, o_))
#    return q
  lerp = libkrebs.PyLerp(o_, s_) # speedup of c++ to python implementation: ca x30!!
  #import time
  #t_ = time.time()
  ret = geometric_transform(data_src, lerp.apply, shape_dst, order=order, mode=mode, cval=cval)
  #print 'resample time: %f' % (time.time()-t_)
  return ret


def CalcIntervascularInterpolationField(edges, radius, pos, values, latticeData, source_factor = 100.):
  '''
   solve: laplace U + lambda (U0 - U) = 0
   where U = local volume average of values
   and lambda is a large factor >> 1.
   edges, radius, pos and values are properties of a graph. pos must be an Nx3 array,
   and values must by an Ex2 array, where N = number of nodes and E = number of edges.
   values is indexed by edges numbers and gives a value at the start and at the end of
   a vessel in the second dimension.
  '''
  edges  = np.asarray(edges, dtype=np.int32)
  pos    = np.asarray(pos  , dtype=np.float32)
  radius = np.asarray(radius, dtype=np.float32)
  values = np.asarray(values, dtype=np.float32)
  assert isinstance(latticeData, LatticeData)
  return CalcIntervascularInterpolationField_(pos, edges, radius, values, latticeData, source_factor, 5)


def SumIsoSurfaceIntersectionWithVessels(level, edges, pressure, flags, nodalLevel, dataValue):
  '''return tuple with (sum of dataValue for inflow, sum of dataValue for outflow)'''
  result = SumIsoSurfaceIntersectionWithVessels_(level,
                                                 np.asarray(edges, np.int32),
                                                 np.asarray(pressure, np.float32),
                                                 np.asarray(flags, np.int32),
                                                 np.asarray(nodalLevel, np.float32),
                                                 np.asarray(dataValue, np.float64))
  return result

def test():
  a = np.arange(12, dtype=np.float).reshape((4,3))
  b = resample_field(a, ((0.,0.),(2.,2.)), ((0.,0.),(1.,1.)), (6,3))
  print a
  print b

def GetWorldBox(vesselgroup):
    if( vesselgroup.attrs['CLASS'] == 'GRAPH' ):
      #ld = read_lattice_data_from_hdf(vesselgroup['lattice'])
      fn=str(vesselgroup.file.filename)
      path=str(vesselgroup.name)
      #ld = read_lattice_data_from_hdf_by_filename(fn, path)
      worldbox = ld.worldBox
      #worldbox = read_lattice_data_from_hdf(vesselgroup['lattice']).GetWorldBox()
    if( vesselgroup.attrs['CLASS'] == 'REALWORLD'):
      pos=vesselgroup['nodes/world_pos']
      x_min = np.min(pos[:,0])
      x_max = np.max(pos[:,0])
      y_min = np.min(pos[:,1])
      y_max = np.max(pos[:,1])
      z_min = np.min(pos[:,2])
      z_max = np.max(pos[:,2])
      worldbox=np.asarray([x_min, x_max, y_min, y_max, z_min, z_max])
    return worldbox
def pyDictFromParamGroup(grp):
  aDict = dict()
#  for key in grp.attrs.keys():
#    aDict[str(key)] = grp.attrs.get(key)
#  return aDict
  if 'calcflow' in grp.name:
    aDict['includePhaseSeparation'] = bool(grp.attrs.get('includePhaseSeparation'))
    aDict['inletHematocrit'] = float(grp.attrs.get('inletHematocrit'))
    aDict['rheology'] = str(grp.attrs.get('rheology'))
    aDict['viscosityPlasma'] = float(grp.attrs.get('viscosityPlasma'))
  return aDict

def test2():
  rad = np.linspace(4., 20., 5)  
  hema = np.linspace(0.45, 0.5, 5)
  length = np.linspace(100., 100., 5)
  print calc_vessel_conductivities(rad, length, hema)
  print CalcViscosities(rad, hema, dict())



#test2()
