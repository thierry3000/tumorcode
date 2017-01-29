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
from os.path import basename, dirname, join, splitext
import sys
import os
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'..'))
import h5py
import h5files

import numpy as np
import extensions # for hdf5 support in np.asarray
import krebsutils
import myutils
import posixpath
import math
import collections

def RemoveArteriovenousFlagsFromCapillaries(flags):
  capil_mask = myutils.bbitwise_and(flags, krebsutils.CAPILLARY)
  flags[capil_mask] = np.bitwise_and(flags[capil_mask], np.asarray(~(krebsutils.ARTERY | krebsutils.VEIN), dtype=flags.dtype))
  return flags
'''returns in \mu m^3 '''
def totalLdVolume(vesselgrp):
  if vesselgrp.attrs['CLASS'] == 'GRAPH':
    ld = krebsutils.read_lattice_data_from_hdf(vesselgrp['lattice'])
    #if we have 2d config one entry of sizeOfDim is zero!  
    sizeOfDim = ld.GetWorldSize()
    sizeOfDim = np.asarray(sizeOfDim)
    vol = np.prod(sizeOfDim[sizeOfDim.nonzero()])
    if vol>0:
      return vol
    else:
      print("to be implemented")
  if vesselgrp.attrs['CLASS'] == 'REALWORLD':
    pos=vesselgrp['nodes/world_pos']
    x_min = np.min(pos[:,0])
    x_max = np.max(pos[:,0])
    delta_x = np.fabs(x_max-x_min)
    #failsave for 2d
    if delta_x == 0:
        delta_x =1
    y_min = np.min(pos[:,1])
    y_max = np.max(pos[:,1])
    delta_y = np.fabs(y_max-y_min)
    if delta_y == 0:
        delta_y =1
    z_min = np.min(pos[:,2])
    z_max = np.max(pos[:,2])
    delta_z = np.fabs(z_max-z_min)
    if delta_z == 0:
        delta_z =1
    vol = delta_x*delta_y*delta_z
    return vol

'''developed with adaption an move here for general usage'''
def getGeometricData(vesselgroups):
  '''returns mvd in 1/mm^2 and blood vessel volume fractions'''
  from analyzeBloodVolumeSimple import cylinderCollectionVolumeDensity, cylinderCollectionLineDensity
  #from analyzeBloodVolumeSimple import computeMeanCapillaryDistance
  data = []
  for g in vesselgroups:
    #meanCapillaryDistance = dataman.obtain_data('basic_vessel_global', 'avg_cap_dist', g, cachelocation(g))
    #meanCapillaryDistance = computeMeanCapillaryDistance(datamanager,dest_group,g)
    rv, rv_a, rv_v, rv_c = cylinderCollectionVolumeDensity(g)
    mvd, mvd_a, mvd_v, mvd_c = cylinderCollectionLineDensity(g)
    radii = np.asarray(g['edges/radius'])
    mean_r = np.mean(radii)
    data.append((rv, rv_a, rv_v, rv_c, mvd, mvd_a, mvd_v, mvd_c, mean_r))
  data = np.asarray(data).transpose() # first dimension gives quantity
  return data

def getTotalPerfusion(vesselgroups):
  '''simply sums up the blood flow from arterial boundary nodes and divides by the system volume.
     Returns list of perfusion values in units of blood volume per volume and sec'''
  And = myutils.bbitwise_and
  data = []
  import analyzeBloodVolumeSimple as anaBloodV
  if 'edges' in vesselgroups.keys():
    g = vesselgroups # only a single group there
    graph = read_vessels_data(g, ['flow', 'flags'])
    #ld    = krebsutils.read_lattice_data_from_hdf(g['lattice'])
    roots = set(g['nodes/roots'][...])
    mask  = np.asarray(map(lambda (a,b): a in roots or b in roots, graph.edgelist), dtype = np.bool)
    vol   = anaBloodV.totalLdVolume(g)
    flow  = graph.edges['flow']
    flags = graph.edges['flags']
    mask  = mask & And(flags, krebsutils.ARTERY)
    flow  = flow[mask]
    perfusion = (np.sum(flow)/vol)
    data.append(perfusion)
  else:
    for g in vesselgroups:
      graph = read_vessels_data(g, ['flow', 'flags'])
      #ld    = krebsutils.read_lattice_data_from_hdf(g['lattice'])
      roots = set(g['nodes/roots'][...])
      mask  = np.asarray(map(lambda (a,b): a in roots or b in roots, graph.edgelist), dtype = np.bool)
      vol   = anaBloodV.totalLdVolume(g)
      flow  = graph.edges['flow']
      flags = graph.edges['flags']
      mask  = mask & And(flags, krebsutils.ARTERY)
      flow  = flow[mask]
      perfusion = (np.sum(flow)/vol)
      data.append(perfusion)
  return np.asarray(data)

def read_vessels_data(vesselgroup, datanames):
  return krebsutils.read_vessels_from_hdf(vesselgroup, datanames, return_graph = True, return_not_found = False)
  
@myutils.UsesDataManager
def generate_rBV_of_group(datamanager, destination_group, f):
  datanames = 'rbv'.split()
  # structure in HDF file:
  #            gmeasure/groupname/
  #                               rbv  <- dataset
  #                               a    <- dataset
  #                               v    <- dataset
  
  def read(gmeasure, groupname):
    gmeasure = gmeasure[groupname]
    return [gmeasure[name][()] for name in datanames ]
    
  def write(gmeasure, groupname):
    if 'adaption' in f.keys():
      vessel_grp = f['adaption/vessels_after_adaption']
    else:
      vessel_grp = f['vessels']
    #group_with_adaption = f['adaption/vessels_after_adaption']
    gmeasure = gmeasure.create_group(groupname)
    for name in datanames:
      geometric_data = getGeometricData([vessel_grp])
      rbv, a, v, c = geometric_data[:4]
      gmeasure.create_dataset(name, data = rbv)
#    
#    for name, data in zip(datanames, [rbv, a, v, c]):
#      gmeasure.create_dataset(name, data = data)
  
  ret = myutils.hdf_data_caching(read, write, destination_group, (f.filename), (1, ))  
  # so, gmeasure/groupname is "/measurements/adaption".(None, 1,) are version number specifications, 
  # one for each path component. None means no version number is checked. If version number is larger than
  # stored number, then data is recomputed instead of loaded.
  return ret
@myutils.UsesDataManager
def generate_adaption_data_average_rBV(datamanager, inputfiles, destination_group, destination_name):
  def process(inputfiles):
    tmp = []
    for f in inputfiles:
      rBV = generate_rBV_of_group(datamanager, destination_group, f)
      tmp.append(rBV)
    avg_rBV = np.average(tmp)
    std_rBV = np.std(tmp)    
#    avg_w = np.average(tmp[0], axis = 0)
#    std_w = np.std(tmp[0], axis = 0)
    return avg_rBV, std_rBV
  
  def write(gmeasure, groupname):
    gmeasure = gmeasure.create_group(groupname)    
  
    #groups_without_adaption = [f['adaption/recomputed'] for f in inputfiles]
    #groups_with_adaption = [f['adaption/vessels_after_adaption'] for f in inputfiles]
    avg_rBV, std_rBV = process(inputfiles)
    #for name in ['with_adaption', 'without_adaption']:
    #  avg_with, std_with, avg_without, std_without = process(groups)
    gmeasure.create_dataset('rBV_avg', data = avg_rBV)
    gmeasure.create_dataset('rBV_std', data = std_rBV)
    #gmeasure.create_dataset('without_adaption_rBV_avg', data = avg_without)
    #gmeasure.create_dataset('without_adaption_rBV_std', data = std_without)
      # will give datasets:
      #     with_adaption_avg, with_adaption_std, without_adaption_avg, and without_adaption_std.
      # Each is tuple containing (rbv, a, v, c) returned by generate_adaption_data_of_group(...)
  
  def read(gmeasure, groupname):
    gmeasure = gmeasure[groupname]
    r1 = gmeasure['rBV_avg']
    r2 = gmeasure['rBV_std']
    #r3 = gmeasure['without_adaption_rBV_avg']
    #r4 = gmeasure['without_adaption_rBV_std']
    return (r1, r2)#, r3, r4)
  
  ret = myutils.hdf_data_caching(read, write, destination_group, (destination_name,), (1,))
  return ret

class BinsSpecRange(object):
  def __init__(self, start, end, step):
    self.p    = start, end, step
  def __hash__(self):
    return hash(self.p)
  def __eq__(self, other):
    return self.p == other.p
  def arange(self):
    return np.arange(self.p[0], self.p[1], self.p[2])
  def __str__(self):
    return 'BinSpecRange'+str(self.p)


class BinsSpecArray(object):
  def __init__(self, arr):
    self.arr = tuple(arr)
  def __hash__(self):
    return hash(self.arr)
  def __eq__(self, other):
    return self.arr == other.arr
  def arange(self):
    return np.asarray(self.arr)
  def __str__(self):
    return 'BinSpecArray'+str(self.arr)


#def calcMVDAsLineDensity(vesselgroup, distancemap, ld, bins):
#  import plotVessels
#  sample_length = 10.
#  far = 1.e10
#  cellvol = ld.scale**3
#  vessels = krebsutils.read_vesselgraph(vesselgroup, ['position','flags'])
#  vessels = vessels.get_filtered(edge_indices = np.asarray(np.bitwise_and(vessels.edges['flags'], krebsutils.CIRCULATED), dtype=np.bool))
#  #s = plotVessels.generate_samples(vessels, 'position', 'nodes', sample_length)
#  s = krebsutils.sample_edges(graph.nodes['position'], graph.edgelist, graph.nodes['position'], sample_length, DATA_LINEAR | DATA_PER_NODE)
#  distancesamples = krebsutils.sample_field(s, distancemap, ld, linear_interpolation=True, extrapolation_value = far)
#  a = myutils.MeanValueArray.fromHistogram1d(bins, distancesamples, np.ones_like(distancesamples))
#  b = myutils.MeanValueArray.fromHistogram1d(bins, distancemap.ravel(), np.ones_like(distancemap.ravel()))
#  a.cnt = b.cnt.copy()
#  a.sum *= sample_length/cellvol
#  a.sqr *= a.sum**2
#  return a


def calc_distmap(ds, ld, level):
  distmap = np.asarray(ds) > level
  distmap = krebsutils.flood_fill(distmap, (0,0,0))
  distmap = np.logical_not(distmap)
  distmap = krebsutils.distancemap(distmap)*ld.scale
  return distmap


def generate_samples(graph, data, association, scale):
  '''A convenence function to call the c++ sampling routine for some data
  stored in my graph structure. Additionaly it converts nodal to edge
  properties because the sampling needs data associated with edges, but it
  can still interpolate between data for the endpoints for each vessel.
  Returns an array of the sampled data.
  data = the name of the data, or the data itself,
  association = 'edges' or 'nodes' or 'avg_edges'
  scale = length which each sample represents approximately. Larger value = less samples.
  '''
  DATA_LINEAR = krebsutils.VesselSamplingFlags.DATA_LINEAR
  DATA_CONST = krebsutils.VesselSamplingFlags.DATA_CONST
  DATA_PER_NODE = krebsutils.VesselSamplingFlags.DATA_PER_NODE
  if isinstance(data, str): data = getattr(graph, association)[data]
  if association == 'edges':
    #data = krebsutils.edge_to_node_property(int(np.amax(graph.edgelist)+1), graph.edgelist, graph.edges[name], 'avg')
    return krebsutils.sample_edges(graph.nodes['position'], graph.edgelist, data, scale, DATA_CONST)
  else:
    return krebsutils.sample_edges(graph.nodes['position'], graph.edgelist, data, scale, DATA_LINEAR | DATA_PER_NODE)


def tumor_path_(f, path):
  group = f[path]
  if 'tumor' in group and isinstance(group['tumor'], h5py.Group): # for simulations with blood vessel network
    return posixpath.join(path, 'tumor')
  else:
    return path

def tumor_group(f, path):
  return f[tumor_path_(f,path)]

def try_find_tumor_group_from_vesselgroup(vesselgroup):
  '''Look for a group named 'tumor'. This is needed to find the tumor group from copies of the tumorsim data.'''
  parent = vesselgroup.parent
  try:
    g = parent['tumor']
    return g
  except KeyError:
    pass
  try:
    g = h5files.openLink(vesselgroup, 'SOURCE')  # for when we have have copy of the original network
  except KeyError:
    pass
  else:
    return try_find_tumor_group_from_vesselgroup(g) # go look in the source which is hopefully the data from the tumor sim
  return None



def obtain_distmap_(dataman, tumorgroup, distance_distribution_name, ld = None):
  if ld is None or distance_distribution_name != 'radial':
    f = tumorgroup.file
    tumor_ld = dataman.obtain_data('ld', f)
  if ld is None:
    ld = tumor_ld        
  if distance_distribution_name == 'radial':    
    return dataman.obtain_data('distance_from_center_distribution', ld), ld
  elif  distance_distribution_name == 'levelset':
    if tumorgroup.attrs['TYPE'] == 'faketumor':
      radius = tumorgroup.attrs['TUMOR_RADIUS']
      field = krebsutils.make_radial_field(ld)
      field -= radius
      return field, ld
    elif tumorgroup.attrs['TYPE'].startswith('BulkTissue'):
      field = dataman.obtain_data('fieldvariable', f, 'dist_tumor', tumorgroup.name)
      if not (ld == tumor_ld):
        field = krebsutils.resample_field(field, tumor_ld.worldBox, ld.shape, ld.worldBox, mode = 'nearest')
        assert all(field.shape == ld.shape)
      return field, ld
    else:
      assert False            

#def resample_field(data_src, bbox_src, shape_dst, bbox_dst, order = 2, mode = 'constant', cval=0.):
class DataDistanceFromCenter(object):
    keywords = [
      'distance_from_center_distribution', 'distancemap_samples'
    ]

    def obtain_data(self, dataman, dataname, *args):
      if dataname == 'distance_from_center_distribution':
        ld, = args
        return krebsutils.make_radial_field(ld)        

      elif dataname == 'distancemap_samples':
        try:
          vesselgroup, tumorgroup, sample_length, distance_distribution_name, ld = args
        except ValueError:
          vesselgroup, tumorgroup, sample_length, distance_distribution_name = args
          ld = None
        pos_smpl   = dataman.obtain_data('basic_vessel_samples', 'position', vesselgroup, sample_length)
        #tumor_ld   = dataman.obtain_data('ld', tumorgroup.file)z
        distmap, ld = obtain_distmap_(dataman, tumorgroup, distance_distribution_name, ld)
        dist_smpl   = krebsutils.sample_field(pos_smpl, distmap, ld, linear_interpolation=True) #, extrapolation_value = 1.e10)
        bmin, bmax  = ld.worldBox.reshape(3,2).transpose()
        mask        = (pos_smpl>bmin) & (pos_smpl<bmax)
        mask        = (mask[:,0] & mask[:,1]) & mask[:,2]
        return dist_smpl, distmap, mask, ld


#radialDistance = myutils.JustHoldsAString('radial')
#levelsetDistance = myutils.JustHoldsAString('levelset')
#
#
#@myutils.UsesDataManager
#def GetRadialDistance(dataman, ld):
#    return krebsutils.make_radial_field(ld), ld
#
#
#def GetLevelsetDistance(dataman, tumorgroup, ld):
#    tumor_ld = dataman.obtain_data('ld', tumorgroup.file)
#    if tumorgroup.attrs['TYPE'] == 'faketumor':
#      radius = tumorgroup.attrs['TUMOR_RADIUS']
#      field = GetRadialDistance(dataman, ld)
#      field -= radius
#      return field, ld
#    elif tumorgroup.attrs['TYPE'].startswith('BulkTissue'):
#      field = dataman.obtain_data('fieldvariable', tumorgroup.file, 'dist_tumor', tumorgroup.name)
#      if not (ld == tumor_ld):
#        field = krebsutils.resample_field(field, tumor_ld.worldBox, ld.shape, ld.worldBox, mode = 'nearest')
#        assert all(field.shape == ld.shape)
#      return field, ld
#    else:
#      assert False  
# 
#
#@myutils.UsesDataManager
#def GetDistanceMap(dataman, tumorgroup, distanceType, ld = None):
#    f = tumorgroup.file
#    tumor_ld = dataman.obtain_data('ld', f)
#    if ld is None:
#      ld = tumor_ld    
#    if distanceType == radialDistance:
#      return GetRadialDistance(dataman, ld)
#    elif  distanceType == levelsetDistance:
#      return GetLevelsetDistance(dataman, tumorgroup, ld)
#
#
#@myutils.UsesDataManager
#def GetDistanceSamples(dataman, vesselgroup, tumorgroup, sample_length, distanceType):
#    pos_smpl    = dataman.obtain_data('basic_vessel_samples', 'position', vesselgroup, sample_length)
#    distmap, ld = GetDistanceMap(dataman, tumorgroup, distanceType)
#    dist_smpl   = krebsutils.sample_field(pos_smpl, distmap, ld, linear_interpolation=True)
#    bmin, bmax  = ld.worldBox.reshape(3,2).transpose()
#    mask        = (pos_smpl>bmin) & (pos_smpl<bmax)
#    mask        = (mask[:,0] & mask[:,1]) & mask[:,2] # position must be within bounds in each dimension
#    return dist_smpl, distmap, ld, mask  


def CalcPhiVessels(dataman, vesselgroup, ld, scaling, samples_per_cell = 5):
  '''samples per cell mean the lattice grid cell,
     the total number of samples for a vessel is determined by the ratio
     of its volume to the volume of a grid cell times the samples_per_cell'''
  # vessel volume fraction
  #vessels = krebsutils.read_vesselgraph(ds, ['radius','position', 'flags'])
  vessels = dataman.obtain_data('vessel_graph', vesselgroup, ['radius', 'position', 'flags'])
  mask = np.nonzero(vessels.edges['flags'] & krebsutils.CIRCULATED)[0]
  vessels = vessels.get_filtered(edge_indices = mask)
  vessels.nodes['position'] *= scaling
  vessels.edges['radius'] *= scaling
  vessel_fraction = krebsutils.make_vessel_volume_fraction_field(
    vessels.nodes['position'],
    vessels.edgelist,
    vessels.edges['radius'],
    ld,
    samples_per_cell)
  return vessel_fraction


class DataBasicVessel(object):
    keywords = [
      'vessel_graph', 'vessel_graph_property', 'vessel_system_length'
    ]

    def get_property(self, dataman, vesselgroup, association, property_name):
      if association == 'auto':
        if   property_name in vesselgroup['edges']: association = 'edges'
        elif property_name in vesselgroup['nodes']: association = 'nodes'
      if property_name == 'position':
        return krebsutils.read_vessel_positions_from_hdf_(vesselgroup).transpose(), 'nodes'
      elif property_name == 'length':
        graph = dataman.obtain_data('vessel_graph', vesselgroup, ['position'])
        pos = graph.nodes['position']
        dp = pos[graph.edgelist[:,0],:] - pos[graph.edgelist[:,1],:]
        l = krebsutils.vector_lengths(dp)
        return l, 'edges'
      elif property_name == 'velocity':
        r,_    = dataman.obtain_data('vessel_graph_property', vesselgroup, 'edges', 'radius')
        flow,_ = dataman.obtain_data('vessel_graph_property', vesselgroup, 'edges', 'flow')
        return (flow/(math.pi*r*r)), 'edges'
      #elif property_name == 'flags':
      #  flags,_    = dataman.obtain_data('vessel_graph_property', vesselgroup, 'edges', 'flags')
      else:
        return np.asarray(vesselgroup[association][property_name]), association

    def obtain_data(self, dataman, dataname, *args):
      if dataname == 'vessel_graph':
        vesselgroup, properties = args
        graph = krebsutils.read_vessels_from_hdf(vesselgroup, [], return_graph=True)
        for prop in properties:
          data, association = self.get_property(dataman, vesselgroup, 'auto', prop)
          getattr(graph, association)[prop] = data
        return graph

      elif dataname == 'vessel_graph_property':
        res, a = self.get_property(dataman, *args)
        return res, a

      elif dataname == 'vessel_system_length':
        group, = args
        def read(gmeasure, groupname):
          return np.asscalar(gmeasure[groupname][...])
        def write(gmeasure, groupname):
          l = np.sum(dataman.obtain_data('vessel_graph_property', group, 'edges', 'length')[0])
          gmeasure.create_dataset(groupname, data = l)
        return myutils.hdf_data_caching(read, write, group, ('vessel_system_length',), (1,))


def combineSamples(obtain_samples, arg_list, every):
  results = []
  for args in arg_list:
    smpl = obtain_samples(args)
    results.append(smpl[::every])
  return np.concatenate(results)


class DataVesselSamples(object):
    keywords = [
      'basic_vessel_samples', 'basic_vessel_samples_avg'
    ]

    def obtain_data(self, dataman, dataname, *args):
      if dataname == 'basic_vessel_samples':
        property_name, vesselgroup, sample_length = args
        #graph = dataman.obtain_data('vessel_graph', vesselgroup, ['position','flags'])
        graph = dataman.obtain_data('vessel_graph', vesselgroup, ['position'])
        if property_name == 'weight':
          smpl = krebsutils.sample_edges_weights(graph.nodes['position'], graph.edgelist, sample_length)
        else:
          data, association = dataman.obtain_data('vessel_graph_property', vesselgroup, 'auto', property_name)
          smpl = generate_samples(graph, data, association, sample_length)
        return smpl

      if dataname == 'basic_vessel_samples_avg':
        property_name, vesselgroups, sample_length, every = args
        def obtain_samples(group):
          return dataman.obtain_data('basic_vessel_samples', property_name, group, sample_length)
        return combineSamples(obtain_samples, vesselgroups, every)




radialAvgPerVolume = object()
radialAvgPerVessels = object()


def GenerateRadialDistributions(dataman, vesselgroup, tumorgroup, sample_length, bins_spec, distance_distribution_name, ld, sample_iterator):
  weight_smpl = dataman.obtain_data('basic_vessel_samples', 'weight', vesselgroup, sample_length)
  flags       = dataman.obtain_data('basic_vessel_samples', 'flags', vesselgroup, sample_length)
  # get teh radial distance function (either distance from tumor border or distance from center)
  dist_smpl, distmap, mask, tumor_ld   = dataman.obtain_data('distancemap_samples', vesselgroup, tumorgroup, sample_length, distance_distribution_name, ld)
  # note: tumor_ld might be another unrelated lattice

   #filter uncirculated
  mask = mask & myutils.bbitwise_and(flags, krebsutils.CIRCULATED)
  dist_smpl   = dist_smpl[mask]
  weight_smpl = weight_smpl[mask]

  res = []
  for (smpl, avg_mode) in sample_iterator:
    bins = bins_spec.arange()
    if avg_mode is radialAvgPerVolume:
      a = myutils.MeanValueArray.fromHistogram1d(bins, dist_smpl, smpl[mask], weight_smpl) # integral over length within bins
      b = myutils.MeanValueArray.fromHistogram1d(bins, distmap.ravel(), np.ones_like(distmap.ravel())) # how much tissue volume in the bins
      a.cnt = b.cnt.copy()
      a.sum *= 1./(tumor_ld.scale**3)
      a.sqr *= a.sum**2  # length integral per tissue volume
      #a *= 1.e6
    elif avg_mode is radialAvgPerVessels:
      a = myutils.MeanValueArray.fromHistogram1d(bins, dist_smpl, smpl[mask], weight_smpl)
    else:
      assert (avg_mode is radialAvgPerVolume or avg_mode is radialAvgPerVessels)
    res.append(a)
  return res


def HdfCacheRadialDistribution((read, write), property_name, bins_spec, distance_distribution_name, cachelocation, version):
  path1 = 'radial_vs_'+distance_distribution_name+'_bins'+str(hash(bins_spec))
  version_args = myutils.checksum(30., 13)
  return myutils.hdf_data_caching(read, write, cachelocation[0], (path1,cachelocation[1],property_name), (version_args,0,version))


#def GetMaskUncirculated(dataman, vesselgroup):
#  flags       = dataman.obtain_data('basic_vessel_samples', 'flags', vesselgroup, sample_length)
#  return myutils.bbitwise_and(flags, krebsutils.CIRCULATED)

def GetMaskArea(dataman, vesselgroup, tumorgroup, distance_distribution_name, ld, distmin, distmax, sample_length = 30.):
  dist_smpl, distmap, mask, _   = dataman.obtain_data('distancemap_samples', vesselgroup, tumorgroup, sample_length, distance_distribution_name, ld)
  mask2 = (dist_smpl>distmin) & (dist_smpl<distmax)
  return mask & mask2


def ComputeSampleHistogram(dataman, vesselgroup, binspec, samples, mask = None, sample_length = 30.):
  weight      = dataman.obtain_data('basic_vessel_samples', 'weight', vesselgroup, sample_length)
  bins = binspec.arange()
  if mask is not None:
    samples = samples[mask]
    weight  = weight[mask]
  return myutils.fromHistogram1d(bins, smpl, 1, weight)

def FixUnits_(name, data):
  if name == 'mvd': return data*1.e6
  return data

def MakeVersionId(inputData):
  if isinstance(inputData, (h5py.Group, h5py.Dataset)):
    try:
      v = inputData.attrs['VERSION']
      return v
    except:
      return 1
  return 1


class DataVesselRadial(object):
    keywords = [
      'basic_vessel_radial'
    ]

    def obtain_data(self, dataman, dataname, *args):
      if dataname == 'basic_vessel_radial':
        property_name, vesselgroup, tumorgroup, sample_length, bins_spec, distance_distribution_name, ld, cachelocation = args

        def write(gmeasure, groupname):
          if property_name == 'mvd':
            smpl = dataman.obtain_data('basic_vessel_samples', 'flags', vesselgroup, sample_length)
            smpl  = np.ones_like(smpl) # only need number of samples
            avg_mode = radialAvgPerVolume
          elif property_name == 'phi_vessels':
            smpl = dataman.obtain_data('basic_vessel_samples', 'radius', vesselgroup, sample_length)
            smpl = math.pi*np.power(smpl, 2.) # pi r^2 * length, where length comes in through GenerateRadialDistributions
            avg_mode = radialAvgPerVolume
          elif property_name in ['S_over_V', 'S_rho']: # alternative names for now ...
            smpl = dataman.obtain_data('basic_vessel_samples', 'radius', vesselgroup, sample_length)
            smpl = (math.pi*2.)*smpl # 2 pi r = surface area
            avg_mode = radialAvgPerVolume
          elif property_name == 'radial_distance':
            smpl, _, _ , _  = dataman.obtain_data('distancemap_samples', vesselgroup, tumorgroup, sample_length, distance_distribution_name, ld)
            avg_mode = radialAvgPerVessels
          else:
            smpl = dataman.obtain_data('basic_vessel_samples', property_name, vesselgroup, sample_length)
            avg_mode = radialAvgPerVessels
          a, = GenerateRadialDistributions(dataman,  vesselgroup, tumorgroup, sample_length, bins_spec, distance_distribution_name, ld, [(smpl, avg_mode)])
          ds = a.write(gmeasure, groupname)
          ds.attrs['BINS_SPEC'] = str(bins_spec)

        def read(gmeasure, groupname):
          assert groupname == property_name
          return FixUnits_(groupname, myutils.MeanValueArray.read(gmeasure, groupname))

        version = myutils.checksum(MakeVersionId(vesselgroup), MakeVersionId(tumorgroup))
        return HdfCacheRadialDistribution((read, write), property_name, bins_spec, distance_distribution_name, cachelocation, version)


class DataVesselGlobal(object):
  keywords = [
    'basic_vessel_global',
    'geometric_data'
  ]

  def obtain_data(self, dataman, dataname, *args):
    if dataname == 'geometric_data':
      property_name, vesselgroup, cachelocation = args
      
      property_list='phi_a phi_v phi_c mvd_a mvd_v mvd_c mean_r'
      def write_geometric(gmeasure, groupname):
        '''returns mvd in 1/mm^2 and blood vessel volume fractions'''
        from analyzeBloodVolumeSimple import cylinderCollectionVolumeDensity, cylinderCollectionLineDensity
        phi, phi_a, phi_v, phi_c = cylinderCollectionVolumeDensity(vesselgroup)
        mvd, mvd_a, mvd_v, mvd_c = cylinderCollectionLineDensity(vesselgroup)
        radii = np.asarray(vesselgroup['edges/radius'])
        mean_r = np.mean(radii)
        gmeasure = gmeasure.create_group(groupname)
        for name in property_list.split():
          gmeasure.create_dataset(name, data = locals()[name])
        #data.append((phi_a, rv_v, rv_c, mvd_a, mvd_v, mvd_c, mean_r))
        #data = np.asarray(data).transpose() # first dimension gives quantity
      def read_geometric(gmeasure, groupname):
        gmeasure = gmeasure[groupname]
        return [gmeasure[name][()] for name in property_list.split() ]
      return myutils.hdf_data_caching(read_geometric, write_geometric, cachelocation[0], ('global', cachelocation[1], 'vessel', property_name), (1,1,MakeVersionId(vesselgroup),1))    
    
    if dataname == 'basic_vessel_global':
      property_name, vesselgroup, cachelocation = args

      def write(gmeasure, groupname):
        ld = krebsutils.read_lattice_data_from_hdf(vesselgroup['lattice'])
        volume = np.prod(ld.GetWorldSize())
        if property_name == 'mvd_linedensity':
          graph = dataman.obtain_data('vessel_graph', vesselgroup, ['length'])
          l = np.asarray(graph['length'], dtype=np.float64)
          data  = (np.sum(l) / volume, 0.)
          data  = [d*1e6 for d in data] #from 1/mum^2 to 1/mm^2
        elif property_name == 'avg_cap_dist':
          vessels = krebsutils.read_vesselgraph(vesselgroup, ['flags', 'length','radius'])
          flags   = RemoveArteriovenousFlagsFromCapillaries(vessels['flags'])
          mask1 = myutils.bbitwise_and(flags, krebsutils.CIRCULATED)
          #totalvol = totalLdVolume(vesselgroup)          
          def compute(flagMask):
            if flagMask:
              mask = mask1 & myutils.bbitwise_and(flags, flagMask)
            else:
              mask = mask1
            length = np.asarray(vessels['length'][mask], dtype = np.float64)
            total = np.sum(length)
            return total/volume
          '''luedemann et. al. assume capillaries are vessels
          smaller than this:
            see manuscript draft
          '''
          mvd_cap = compute((vessels['radius']<=4.0).all())
          data = (np.sqrt(1/mvd_cap), 0.)
        elif property_name == 'phi_vessels':
          graph = dataman.obtain_data('vessel_graph', vesselgroup, ['length','radius'])
          l = np.asarray(graph['length'], dtype=np.float64)
          r = np.asarray(graph['radius'], dtype=np.float64)
          data = l*np.square(r)*math.pi
          data = (np.sum(data) / volume, 0.)
        elif property_name == 'total_perfusion':
          data = getTotalPerfusion(vesselgroup)*60 #to minutes
        else:
          data   = dataman.obtain_data('basic_vessel_samples', property_name, vesselgroup, 30.)
          weight = dataman.obtain_data('basic_vessel_samples', 'weight', vesselgroup, 30.)
          data = myutils.WeightedAverageStd(data, weights=weight)
        gmeasure.create_dataset(groupname, data = data)
      def read(gmeasure, groupname):
        #what did he here?-> well for radial average it makes sense to return the avg and the std
        # changed to only return the average
        return FixUnits_(groupname, gmeasure[groupname][0])
      return myutils.hdf_data_caching(read, write, cachelocation[0], ('global', cachelocation[1], 'vessel', property_name), (1,1,MakeVersionId(vesselgroup),1))


def ApproximateTumorRadius(dataman, tumorgroup):
  if tumorgroup.attrs.get('TYPE','')=='faketumor':
    return tumorgroup.attrs['TUMOR_RADIUS']
  else:
    tumor_ld    = dataman.obtain_data('ld', tumorgroup.file)
    radial      = dataman.obtain_data('distance_from_center_distribution', tumor_ld)
    distancemap = dataman.obtain_data('fieldvariable', tumorgroup.file, 'dist_tumor', tumorgroup.name)
    radius      = np.average(radial[np.nonzero(np.logical_and(distancemap>-2*tumor_ld.scale, distancemap<2*tumor_ld.scale))])
    return radius


class DataTumorTissueSingle(object):
  keywords = [
    'fieldvariable', 'fieldvariable_radial','ld','time', 'approximate_tumor_radius',
  ]

  def obtain_data(self, dataman, dataname, *args):
    if dataname == 'ld':
      f, = args
      ld = krebsutils.read_lattice_data_from_hdf(f['field_ld'])
      return ld

    ####
    if dataname == 'time':
      if len (args) == 1:
        group, = args
      else:
        group = args[0][args[1]] # f, group = args
      return group.attrs['time'] # in hours

    #####
    if dataname == 'fieldvariable':
      f, fieldname, group = args[:3]
      g = f[group] # group name to actual group
      tumor_path = tumor_path_(f, group)

      ld = dataman.obtain_data('ld', f)

      def get(name):
        return self.obtain_data(dataman, 'fieldvariable', f, name, group)

      if fieldname == 'phi_cells':
        data = f[tumor_path]['conc']
      elif fieldname == 'dist_tumor_':
        data = f[tumor_path]['ls']
      elif fieldname == 'theta_tumor':
        gtumor = f[tumor_path]
        if gtumor.attrs['TYPE'] == 'faketumor':
          data = gtumor['tc_density']
        else:
          data = gtumor['ptc']
      elif fieldname == 'phi_necro':
        data = f[tumor_path]['necro']
      elif fieldname in ('oxy'):
        try:
          data = g['oxy']
        except KeyError:
          data = g['fieldOxy']
      elif fieldname in ('gf'):
        data = g['fieldGf']
      elif fieldname in ('press', 'sources', 'vel'):
        data = f[tumor_path][fieldname]
      elif fieldname == 'phi_viabletumor':
        theta_tumor = get('theta_tumor')
        phi_cells = get('phi_cells')
        phi_necro = get('phi_necro')
        data = theta_tumor * (phi_cells - phi_necro)
      elif fieldname == 'phi_tumor':
        theta_tumor = get('theta_tumor')
        phi_cells = get('phi_cells')
        data = theta_tumor * phi_cells
      elif fieldname == 'dist_necro':
        phi_necro = get('phi_necro')
        phi_cells = get('phi_cells')
        tumor_contour_level = 0.5*np.average(phi_cells)
        data = calc_distmap(phi_necro, ld, tumor_contour_level)
      elif fieldname == 'dist_viabletumor':
        def read(gmeasure, dsname):
          return np.asarray(gmeasure[dsname])
        def write(gmeasure, dsname):
          dist_tumor = get('dist_tumor')
          dist_necro = get('dist_necro')
          data = np.maximum(dist_tumor, -dist_necro)
          gmeasure.create_dataset(dsname, data = data, compression = 9)
        fm = myutils.MeasurementFile(f, h5files)
        data = myutils.hdf_data_caching(read, write, fm, (tumor_path, 'dist_viabletumor'), (0,1))
      elif fieldname == 'phi_vessels':
        def read(gmeasure, name):
          return np.asarray(gmeasure[name])
        def write(gmeasure, name):
          phi_vessels = CalcPhiVessels(dataman, f[group]['vessels'], ld, scaling = 1.)
          gmeasure.create_dataset(name, data = phi_vessels, compression = 9)
        fm = myutils.MeasurementFile(f, h5files)
        data = myutils.hdf_data_caching(read, write, fm, (group, 'phi_vessels'), (0,1))
      elif fieldname == 'dist_tumor':
        def read(gmeasure, name):
          return np.asarray(gmeasure[name])
        def write(gmeasure, name):
          ls = get('dist_tumor_')
          ls = -ls
          dist = calc_distmap(np.asarray(ls < 0., dtype=np.float32), ld, 0.5)
          gmeasure.create_dataset(name, data = dist, compression = 9)
        fm = myutils.MeasurementFile(f, h5files)
        data = myutils.hdf_data_caching(read, write, fm, (tumor_path, 'dist_tumor_full'), (0,1))
      else:
        raise RuntimeError('unkown field %s' % fieldname)

      if len(args)>3 and args[3] == 'imslice':
        import plotBulkTissue
        return plotBulkTissue.imslice(data)
      else:
        return np.asarray(data)

    if dataname == 'fieldvariable_radial':
      property_name, tumorgroup, bins_spec, distance_distribution_name, cachelocation = args

      def write(gmeasure, groupname):
        distmap, ld = obtain_distmap_(dataman, tumorgroup, distance_distribution_name)
        data    = dataman.obtain_data('fieldvariable', tumorgroup.file, property_name, tumorgroup.name)
        bins    = bins_spec.arange()
        a = myutils.MeanValueArray.fromHistogram1d(bins, np.ravel(distmap), np.ravel(data))
        ds = a.write(gmeasure, groupname)
        ds.attrs['BINS_SPEC'] = str(bins_spec)

      def read(gmeasure, groupname):
        assert groupname == property_name
        return myutils.MeanValueArray.read(gmeasure, groupname)

      return HdfCacheRadialDistribution((read, write), property_name, bins_spec, distance_distribution_name, cachelocation, 1)

    if dataname == 'approximate_tumor_radius':
      tumorgroup, = args
      def write(gmeasure, groupname):
        rad = ApproximateTumorRadius(dataman, tumorgroup)
        ds = gmeasure.create_dataset(groupname, data = rad)
      def read(gmeasure, groupname):
        return np.asscalar(gmeasure[groupname][()])
      return myutils.hdf_data_caching(read, write, tumorgroup, ('approximate_tumor_radius',), (1,))
