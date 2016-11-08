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
from os.path import basename
import os,sys
import h5py
import h5files
import numpy as np
import extensions # for hdf5 support in np.asarray
import myutils
import posixpath
import math
from copy import deepcopy

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
  return pickDetailedO2Library(parameters).computeSaturation_(po2, parameters)

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
  # then we recompute blood flow because the alorithm has changed and we may or may not want hematocrit
  pressure, flow, shearforce, hematocrit, flags = krebsutils.calc_vessel_hydrodynamics(gv, return_flags = True, bloodflowparams = bloodflowparams)
  # then we save the new data to complete the network copy
  gvedst.create_dataset('flow'      , data = flow       , compression = 9)
  gvedst.create_dataset('shearforce', data = shearforce , compression = 9)
  gvedst.create_dataset('hematocrit', data = hematocrit , compression = 9)
  gvedst.create_dataset('flags'     , data = flags      , compression = 9)
  gvndst.create_dataset('pressure'  , data = pressure   , compression = 9)


def computePO2_(gdst, vesselgroup, tumorgroup, parameters):
  myutils.buildLink(gdst, 'SOURCE_VESSELS', vesselgroup)
  myutils.buildLink(gdst, 'SOURCE_TISSUE', tumorgroup)
  gdst.attrs['PARAMETER_CHECKSUM'] = myutils.checksum(parameters)
  gdst.attrs['UUID']               = myutils.uuidstr()
  #writing parameters as acsii string. but i think it is better to store them as hdf datasets directly
  #ds = gdst.create_dataset('PARAMETERS_AS_JSON', data = np.string_(json.dumps(parameters, indent=2)))
  myutils.hdf_write_dict_hierarchy(gdst, 'parameters', parameters)
  gdst.file.flush()
  # now we can compute PO2. The c++ routine can a.t.m. only read from hdf so it has to use our new file
  pickDetailedO2Library(parameters).computePO2(vesselgroup, tumorgroup, parameters, parameters.get('calcflow'), gdst)
  #r = pickDetailedO2Library(parameters).computePO2(vesselgroup, tumorgroup, parameters, parameters.get('calcflow'), gdst)
  #return r


def readParameters(po2group):
  p = myutils.hdf_read_dict_hierarchy(po2group['parameters'])
  for oldName, newName in [("kD_tissue", 'D_tissue'),
                           ("alpha_t", "solubility_tissue"),
                           ("alpha_p", "solubility_plasma")]:
    if not newName in p:
      p[newName] = p[oldName]
      del p[oldName]
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
  refVessels, refTumor = getSourceRefs_(po2group)
  if refVessels.fn:
    gvessels = h5files.open(refVessels.fn, 'r+', relatedObjectForSearch = po2group)[refVessels.path]
  else:
    gvessels = po2group.file[refVessels.path] # look in the same file
    # expect it to return someting. cannot have no vessel network
  if refTumor.fn in os.listdir('.'):
    gtumor = h5files.open(refTumor.fn, 'r+', relatedObjectForSearch = po2group)[refTumor.path] if refTumor.fn else None # maybe there is no tumor
  else:
    print('Warning: did not find tumor in same folder!')
    if refTumor.fn in os.listdir('../'):
      gtumor = h5files.open('../'+refTumor.fn, 'r+', relatedObjectForSearch = po2group)[refTumor.path] if refTumor.fn else None # maybe there is no tumor
    else:
      print('Warning: also not in top level folder')
      if refTumor.fn in os.listdir('../../'):
        gtumor = h5files.open('../../'+refTumor.fn, 'r+', relatedObjectForSearch = po2group)[refTumor.path] if refTumor.fn else None # maybe there is no tumor
      else:
        print('also not in top top level')
        gtumor = None
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
  else: # if is initial vessel file
    d = [('ENSEMBLE_INDEX', names[0]), 
         ('MESSAGE', names[1])]
    for n1, n2 in d:
      fdst.attrs[n2] = fsrc['parameters'].attrs[n1]


def computePO2(f, group_path, parameters, cachelocation):
  if not isinstance(f, h5py.File):
    f = h5files.open(f, 'r+', search = False)  # we open with r+ because the cachelocation might be this file so we need to be able to write to it
    return computePO2(f, group_path, parameters, cachelocation) # recurse

  #==== is this from a tumor sim, or just a vessel network? get hdf group objects =====#
  group = f[group_path]
  tumorgroup = None
  vesselgroup1 = group
  if 'vessels' in group:
    vesselgroup1 = group['vessels']
    if 'tumor' in group:
      tumorgroup = group['tumor']

  #====  this is for recomputing flow =====#
  def read1(gmeasure, name):
    gmeasure = gmeasure[name]
    return myutils.H5FileReference(gmeasure.file.filename, gmeasure.name)

  def write1(gmeasure, name):
    gdst = gmeasure.create_group(name)
    copyVesselnetworkAndComputeFlow(gdst, vesselgroup1, parameters.get("calcflow"))

  #==== execute reading or computing and writing =====#
  fm = h5files.open(cachelocation[0], 'a', search = False) # this is the file where o2 stuff is stored
  CopyInputFileInfo_(fm, f)
  new_flow_data_ref = myutils.hdf_data_caching(read1, write1, fm, ('recomputed_flow', cachelocation[1]), (0, 1))
  fv = h5files.open(new_flow_data_ref.fn,'r+', search = False)
  vesselgroup2 = fv[new_flow_data_ref.path]

  #====  this is for po2 =====#
  def read2(gmeasure, name):
    gmeasure = gmeasure[name]
    return myutils.H5FileReference(gmeasure.file.filename, gmeasure.name)

  def write2(gmeasure, name):
    gdst = gmeasure.create_group(name)
    myutils.buildLink(gdst, 'SOURCE', group)
    computePO2_(gdst, vesselgroup2, tumorgroup, parameters)

  #==== execute reading or computing and writing =====#
  o2data_ref = myutils.hdf_data_caching(read2, write2, fm, ('po2',cachelocation[1]), (1,1))
  #=== return filename and path to po2 data ====#
  return o2data_ref


##############################################################################


def DiffRadiusToConsCoeff(rd, kd):
  return kd/(rd*rd)

chb_of_rbcs = 340./64458. # in mol/l, entspricht 34 g/dl


def doit(fn, pattern, (parameters, parameters_name)):
  #fnpath = dirname(fn)
  fnbase = basename(fn).rsplit('.h5')[0]
  #outfn_no_ext = join(fnpath, fnbase+'_detailedpo2')

  output_links = []

  f = h5files.open(fn, 'r')
  dirs = myutils.walkh5(f['.'], pattern)
  for group_path in dirs:
    #cachelocation = (outfn_no_ext+'.h5', group_path+'_'+parameters_name)
    cachelocation = ('o2_' + fnbase+'_'+parameters_name+'.h5', group_path)
    ref = computePO2(f, group_path, parameters, cachelocation)
    print 'computed po2 stored in:', ref
    output_links.append(ref)
  return output_links

