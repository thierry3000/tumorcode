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
import os, sys
from os.path import join, basename, dirname, splitext
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../..'))

from os.path import basename
import os,sys
import h5py

import numpy as np
import extensions # for hdf5 support in np.asarray
import myutils
import posixpath
import math
from copy import deepcopy

from krebsjobs.parameters import parameterSetsO2
import krebsutils # import of this must come in front of import of detailedo2 libs because some required initialization stuff on the c++ side (see mainboost.cpp)
if sys.flags.debug:
  print 'detailedo2 DEBUG MODE'
  detailedo2current = __import__('libdetailedo2_d', globals(), locals())
  #detailedo2legacy = __import__('libdetailedo2legacy_d', globals(), locals())
else:
  detailedo2current = __import__('libdetailedo2_', globals(), locals())
  #detailedo2legacy = __import__('libdetailedo2legacy_', globals(), locals())


def pickDetailedO2Library(parameters):
  '''This is needed for compatibility. When old files are no longer relevant this could be deleted and the legacy version discarded'''
  # usesNewModel = 'conductivity_coeff1' in parameters
  return detailedo2current #if usesNewModel else detailedo2legacy

def HandleNdimensions(*argsToRavel):
  def wrapper(func):
    def wrapped(*args):
      newargs = list(args)
      for i in argsToRavel:
        a = np.ascontiguousarray(args[i], dtype = np.float32)
        shape = a.shape
        newargs[i] = a.ravel()
        del a
      result = func(*newargs)
      result = np.reshape(result, shape)
      return result
    return wrapped
  return wrapper

@HandleNdimensions(0)
def PO2ToSaturation(po2, parameters):
  '''call c++ also convert arrays to the expected format'''
  print('shape of po2: %s' % po2.shape)
  sat = detailedo2current.computeSaturation_(po2, parameters)
  print('calculated saturation:')
  print(sat)
  return sat
@HandleNdimensions(0, 1)
def ConcentrationToPO2(conc, hema, parameters):
  return pickDetailedO2Library(parameters).computePO2FromConc_(conc, hema, parameters)
  
@HandleNdimensions(0, 1)
def PO2ToConcentration(po2, hema, parameters):
  return pickDetailedO2Library(parameters).computeConcentration_(po2, hema, parameters)

@HandleNdimensions(0)
def computeMassTransferCoefficient(radius, parameters):
  return pickDetailedO2Library(parameters).computeMassTransferCoefficient_(radius, parameters)

def copyVesselnetworkAndComputeFlow(gvdst, gv, bloodflowparams):
  '''gdst = group where the data is placed in, does not create a 'vesse' folder in it but writes nodes, edges directly;
     gv   = source vessel group
  '''
  if gv.attrs['CLASS'] == 'GRAPH':
    gvdst.attrs['CLASS'] = 'GRAPH'
    myutils.buildLink(gvdst, 'SOURCE', gv)
    # first we need to copy some of the vessel data
    gvedst = gvdst.create_group('edges')
    gvndst = gvdst.create_group('nodes')
    gvedst.attrs['COUNT'] = gv['edges'].attrs['COUNT']
    gvndst.attrs['COUNT'] = gv['nodes'].attrs['COUNT']
    gv.copy('lattice',gvdst)
    for name in ['lattice_pos', 'roots','bc_conductivity_value','bc_node_index','bc_type','bc_value']:
      gv['nodes'].copy(name, gvndst)
    for name in ['radius', 'node_a_index', 'node_b_index']:
      gv['edges'].copy(name, gvedst)
  if gv.attrs['CLASS'] == 'REALWORLD':
    gvdst.attrs['CLASS'] = 'REALWORLD'
    myutils.buildLink(gvdst, 'SOURCE', gv)
    # first we need to copy some of the vessel data
    gvedst = gvdst.create_group('edges')
    gvndst = gvdst.create_group('nodes')
    gvedst.attrs['COUNT'] = gv['edges'].attrs['COUNT']
    gvndst.attrs['COUNT'] = gv['nodes'].attrs['COUNT']
    #gv.copy('lattice',gvdst)
    for name in ['world_pos', 'roots','bc_conductivity_value','bc_node_index','bc_type','bc_value']:
      gv['nodes'].copy(name, gvndst)
    for name in ['radius', 'node_a_index', 'node_b_index']:
      gv['edges'].copy(name, gvedst)
  # then we recompute blood flow because the alorithm has changed and we may or may not want hematocrit
  pressure, flow, shearforce, hematocrit, flags = krebsutils.calc_vessel_hydrodynamics_(gv, calc_hematocrit=True, return_flags = True, override_hematocrit= True, bloodflowparams = bloodflowparams,storeCalculationInHDF=False)
  # then we save the new data to complete the network copy
  gvedst.create_dataset('flow'      , data = flow       , compression = 9)
  gvedst.create_dataset('shearforce', data = shearforce , compression = 9)
  gvedst.create_dataset('hematocrit', data = hematocrit , compression = 9)
  gvedst.create_dataset('flags'     , data = flags      , compression = 9)
  gvndst.create_dataset('pressure'  , data = pressure   , compression = 9)


#def computePO2_(gdst, vesselgroup, tumorgroup, parameters):
#  myutils.buildLink(gdst, 'SOURCE_VESSELS', vesselgroup)
#  myutils.buildLink(gdst, 'SOURCE_TISSUE', tumorgroup)
#  gdst.attrs['PARAMETER_CHECKSUM'] = myutils.checksum(parameters)
#  gdst.attrs['UUID']               = myutils.uuidstr()
#  #writing parameters as acsii string. but i think it is better to store them as hdf datasets directly
#  #ds = gdst.create_dataset('PARAMETERS_AS_JSON', data = np.string_(json.dumps(parameters, indent=2)))
#  myutils.hdf_write_dict_hierarchy(gdst, 'parameters', parameters)
#  gdst.file.flush()
#  # now we can compute PO2. The c++ routine can a.t.m. only read from hdf so it has to use our new file
#  #pickDetailedO2Library(parameters).computePO2(vesselgroup, tumorgroup, parameters, parameters.get('calcflow'), gdst)
#  fn=str(vesselgroup.file.filename)
#  vessel_path=str(vesselgroup.name)
#  if( not tumorgroup):
#    tumor_path="not_found_tumor"
#  else:
#    tumor_path=str(tumorgroup.name)  
#  pickDetailedO2Library(parameters).computePO2(fn, vessel_path, tumor_path, parameters, parameters.get('calcflow'), str("bla"))
#  #r = pickDetailedO2Library(parameters).computePO2(vesselgroup, tumorgroup, parameters, parameters.get('calcflow'), gdst)
#  #return r


def readParameters(po2group):
  print('detailed o2 group found at %s' % po2group)
  if('simType' in po2group.parent.attrs.keys()):
    if(po2group.parent.attrs['simType'] == 'MTS'):
      print('assume o2 parameters at: /parameters/o2_sim')
      p = dict()
      for aKey in po2group.parent.parent.parent['parameters/o2_sim'].attrs.keys():
        temp = po2group.parent.parent.parent['parameters/o2_sim'].attrs[aKey]
        if (type(temp) == type(str())):
          p[aKey] = temp
        else:
          p[aKey] = np.asscalar(temp)
      #o2_paramset_name = po2group.parent.parent.parent['parameters/o2_sim'].attrs['detailedO2name']
      #p = getattr(parameterSetsO2, o2_paramset_name)
      #p=myutils.hdf_read_dict_hierarchy_attr(po2group.parent.parent['parameters/o2_params'])
  else:
    p = dict()
    for aKey in po2group.parent.parent['parameters/o2'].attrs.keys():
      temp = po2group.parent.parent['parameters/o2'].attrs[aKey]
      if (type(temp) == type(str())):
        p[aKey] = temp
      else:
        p[aKey] = np.asscalar(temp)
  return p
  

def sampleVessels(po2group, vesselgroup, tumorgroup, sample_length):
  po2vessels  = np.asarray(po2group['po2vessels'])
  if po2vessels.shape[0] <> 2: po2vessels = np.transpose(po2vessels)
  po2field  = np.asarray(po2group['po2field'])
  ld = krebsutils.read_lattice_data_from_hdf(po2group['field_ld'])
  parameters = readParameters(po2group)
  # call thee c++ stuff
  samples, fluxes = pickDetailedO2Library(parameters).sampleVessels(vesselgroup, tumorgroup, parameters, po2vessels, po2field, ld, sample_length)
  # convert units to mlO2/s, and also return dicts
  fluxes = dict((k, f*60.e-12) for k,f in zip(('Jin_root', 'Jout_root', 'Jout_tv', 'tv_cons'), fluxes))
  samples = dict(
      po2    = samples[0],
      extpo2 = samples[1],
      jtv    = samples[2]*60.e-4,  # 1.e-4 to convert from mu m^3 O2 / um^2 / s to ml O2 / cm^2 / min
      dS_dx  = samples[3]*1.e3, # in 1/mmm
    )
  return samples, fluxes


def computeO2Uptake(po2group, tumorgroup):
  po2field  = np.asarray(po2group['po2field'])
  po2ld = krebsutils.read_lattice_data_from_hdf(po2group['field_ld'])
  parameters = readParameters(po2group)
  uptake = pickDetailedO2Library(parameters).computeO2Uptake(po2field, po2ld, tumorgroup, parameters)
  uptake *= 60. # 1/s -> 1/min
  return uptake;


def getSourceRefs_(po2group):
  srcVessels = myutils.H5FileReference(po2group.attrs['SOURCE_VESSELS_FILE'],po2group.attrs['SOURCE_VESSELS_PATH'])
  srcTissue  = myutils.H5FileReference(po2group.attrs['SOURCE_TISSUE_FILE'],po2group.attrs['SOURCE_TISSUE_PATH'])
  return srcVessels, srcTissue

def OpenVesselAndTumorGroups(po2group):
  print('attempt to read from %s' % po2group)
  if('simType' in po2group.parent.attrs.keys()):
    if(po2group.parent.attrs['simType']=='MTS'):
      print('found MTS simulation')
      gvessels = po2group.parent.parent['vessels']
      gtumor = None
  else:
    if 'po2' in po2group.name: #pure o2 sim
      gvessels = po2group.parent.parent['recomputed_flow/vessels']
      gtumor = None
    else: #fakeTumor??
      vessel_input_file = h5py.File(po2group.attrs['SOURCE_VESSELS_FILE'],'r')
      gvessels = vessel_input_file[po2group.attrs['SOURCE_VESSELS_PATH']]
      gtumor = None
#    refVessels, refTumor = getSourceRefs_(po2group)
#    if refVessels.fn:
#      gvessels = h5files.open(refVessels.fn, 'r+', relatedObjectForSearch = po2group)[refVessels.path]
#    else:
#      gvessels = po2group.file[refVessels.path] # look in the same file
#      # expect it to return someting. cannot have no vessel network
#    if refTumor.fn in os.listdir('.'):
#      gtumor = h5files.open(refTumor.fn, 'r+', relatedObjectForSearch = po2group)[refTumor.path] if refTumor.fn else None # maybe there is no tumor
#    else:
#      print('Warning: did not find tumor in same folder!')
#      if refTumor.fn in os.listdir('../'):
#        gtumor = h5files.open('../'+refTumor.fn, 'r+', relatedObjectForSearch = po2group)[refTumor.path] if refTumor.fn else None # maybe there is no tumor
#      else:
#        print('Warning: also not in top level folder')
#        if refTumor.fn in os.listdir('../../'):
#          gtumor = h5files.open('../../'+refTumor.fn, 'r+', relatedObjectForSearch = po2group)[refTumor.path] if refTumor.fn else None # maybe there is no tumor
#        else:
#          print('also not in top top level')
#          gtumor = None
  return gvessels, gtumor


##############################################################################
### somewhat higher level function. Use on a time snapshot group from the tumor sim.
##############################################################################

# save what initial vessel network (type) was used
# WARNING: untested, expect some trouble here
def CopyInputFileInfo_(fdst, fsrc):
  names = ['VESSELFILE_ENSEMBLE_INDEX','VESSELFILE_MESSAGE']
  if names[0] in fsrc.attrs: # if is tumor file
    for name in names:
      fdst.attrs[name] = fsrc.attrs[name]
  elif 'ENSEMBLE_INDEX' in fsrc and 'MESSAGE' in fsrc: # if is initial vessel file
    d = [('ENSEMBLE_INDEX', names[0]), 
         ('MESSAGE', names[1])]
    if 'parameters' in fsrc:
      for n1, n2 in d:
        fdst.attrs[n2] = fsrc['parameters'].attrs[n1]


def computePO2(parameters):
  print("Computing o2 for file: %s" % parameters['input_file_name'])
  print("at group: %s" % parameters['input_group_path'])
  parameters['vessel_group_path'] = "recomputed_flow"
  parameters['output_group_path'] = "po2/" + parameters['input_group_path']
  output_buffer_name = basename(parameters['output_file_name']).rsplit('.h5')[0]
  parameters['output_file_name'] = "%s-%s.h5" % (output_buffer_name, parameters['input_group_path'])  
  print("storing in file: %s at %s" % (parameters['output_file_name'], parameters['output_group_path']))
  tumorgroup = None
  parameters['tumor_file_name'] = 'none'
  parameters['tumor_group_path'] = 'none'
  
  #f_out = h5files.open(parameters['output_file_name'], 'a', search = False)
  
  
      
  caching = False;
  if caching:
     
    #====  this is for recomputing flow =====#
    def read1(gmeasure, name):
      gmeasure = gmeasure[name]
      return myutils.H5FileReference(gmeasure.file.filename, gmeasure.name)
  
    def write1(gmeasure, name):
      gdst = gmeasure.create_group(name)
      f = h5py.File(parameters['input_file_name'])
      #f = h5files.open(parameters['input_file_name'], 'r')
      input_vessel_group = f[parameters['input_group_path']]
      copyVesselnetworkAndComputeFlow(gdst, input_vessel_group, parameters.get("calcflow"))
      f.close()
      
    new_flow_data_ref = myutils.hdf_data_caching(read1, write1, f_out, ('recomputed_flow'), (0, 1))
  
  #  #====  this is for po2 =====#
    def read2(gmeasure, name):
      gmeasure = gmeasure[name]
      return myutils.H5FileReference(gmeasure.file.filename, gmeasure.name)
  
    def write2(gmeasure, name):
      f_out.create_group("po2")
      detailedo2current.computePO2(parameters, parameters.get('calcflow')) 
      
  #
  #  #==== execute reading or computing and writing =====#
    with h5py.File(parameters['output_file_name'], 'a') as f_out:
      #f_out_bla = h5files.open("blub.h5", 'a', search = False)
      o2data_ref = myutils.hdf_data_caching(read2, write2, f_out, ('po2'), (0,1))
    #  #=== return filename and path to po2 data ====#
    #  #return o2data_re
  else:
    print("no caching!")
#    with h5py.File(parameters['input_file_name'], 'r') as f:
#      input_vessel_group = f[parameters['input_group_path']]
#      with h5py.File(parameters['output_file_name'], 'a') as f_out:
#        gdst = f_out.create_group('/recomputed_flow/' + parameters['input_group_path'])
#        copyVesselnetworkAndComputeFlow(gdst, input_vessel_group, parameters.get("calcflow"))
#        f_out.flush()
    #at this point all h5Files should be closed on python side
    detailedo2current.computePO2(parameters, parameters.get('calcflow')) 
  return parameters['output_file_name']

##############################################################################


def DiffRadiusToConsCoeff(rd, kd):
  return kd/(rd*rd)

chb_of_rbcs = 340./64458. # in mol/l, entspricht 34 g/dl


def doit(fn, pattern, (parameters, parameters_name)):
  #fnpath = dirname(fn)
  fnbase = basename(fn).rsplit('.h5')[0]
  #outfn_no_ext = join(fnpath, fnbase+'_detailedpo2')

  output_links = []

  with h5py.File(fn, 'r') as f_in_doit:
    dirs = myutils.walkh5(f_in_doit['.'], pattern)
  
  for group_path in dirs:
    #cachelocation = (outfn_no_ext+'.h5', group_path+'_'+parameters_name)
    #output_file_name = 'o2_' + fnbase+'_'+parameters_name+'.h5'
    #vessel_path = group_path
    parameters['input_file_name'] = fn
    parameters['input_group_path'] = group_path
    parameters['output_file_name']='o2_' + fnbase+'_'+parameters_name+'.h5'
    #cachelocation = ('o2_' + fnbase+'_'+parameters_name+'.h5', group_path)
    computePO2(parameters)
    print('computed po2 stored in: %s' % parameters['output_file_name'])
    #output_links.append(ref)
  #return output_links
  
if __name__ == '__main__':
  o2params = getattr(parameterSetsO2, 'default_o2')
  o2params['output_file_name'] = 'o2_vessels-default-typeI-15x19L130-sample00_default_o2-vessels.h5'
  o2params['vessel_group_path'] = 'recomputed_flow'
#  o2params['output_file_name'] = 'blub42.h5'
  detailedo2current.computePO2(o2params, o2params.get('calcflow')) 

