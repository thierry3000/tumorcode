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
import krebsutils as _ku # import of this must come in front of import of detailedo2 libs because some required initialization stuff on the c++ side (see mainboost.cpp)
from os.path import basename, dirname, join, splitext
import os, os.path, sys
import h5py
import h5files
import uuid
import numpy as np
#import extensions # for hdf5 support in np.asarray
import myutils
#import posixpath
#import math
import warnings
#from copy import deepcopy
#import krebs
from krebs.analyzeGeneral   import DataBasicVessel
from krebs.analyzeBloodFlowResistance import ComputeVascularTreeBloodFlowResistances


sys.path.append(os.path.join(os.path.dirname(__file__),'../../../lib'))

if sys.flags.debug:
  adaption_cpp = __import__('libadaption_d', globals(), locals())
else:
  adaption_cpp = __import__('libadaption_', globals(), locals())

def worker_on_client(fn, grp_pattern, adaptionParams, num_threads=1):
  print('Adaption on %s / %s / param: %s' % (fn, grp_pattern, adaptionParams['name']))
  #h5files.search_paths = [dirname(fn)] # so the plotting and measurement scripts can find the original tumor files using the stored basename alone
  _ku.set_num_threads(num_threads)
  
  #params['name'] = parameter_set_name
  adaptionParams['adaption'].update(
      vesselFileName = fn,
      vesselGroupName = grp_pattern,
      )

  doit( adaptionParams)
  
  #h5files.closeall() # just to be sure

#def buildLink(gdst, linkname, gsrc):
#  fn = myutils.sanitize_posixpath(basename(gsrc.file.filename)) if gsrc else ''
#  gdst.attrs[linkname+'_PATH'] = gsrc.name if gsrc else ''
#  gdst.attrs[linkname+'_FILE'] = fn
#  if gsrc:
#    gdst[linkname] = h5py.SoftLink(gsrc.name) if (gsrc.file == gdst.file) else h5py.ExternalLink(fn, gsrc.name)
def doit_optimize(vesselFileName,adaptParams,BfParams):
  returns = adaption_cpp.doAdaptionOptimization(vesselFileName,adaptParams,BfParams)
  print("should be optimized vaules:")  
  print(returns)
  print('from file: %s' %vesselFileName)
  f_results = h5files.open("optimize_results.h5", 'a')
  a_uuid = str(uuid.uuid4())
  g = f_results.create_group(vesselFileName + '_' + a_uuid)
  g.create_dataset('x_opt', data=returns)
  myutils.hdf_write_dict_hierarchy(g, 'adaptParams', adaptParams)
  myutils.hdf_write_dict_hierarchy(g, 'BfParams', BfParams)  
  g.file.flush()  
  #f_results.close()
  return returns
def copyVesselnetworkAndComputeFlow(gvdst, gv, bloodflowparams):
  '''gdst = group where the data is placed in, does not create a 'vesse' folder in it but writes nodes, edges directly;
     gv   = source vessel group
  '''
  myutils.buildLink(gvdst, 'CLASS', gv)
  myutils.buildLink(gvdst, 'SOURCE', gv)
  # first we need to copy some of the vessel data
  gvedst = gvdst.create_group('edges')
  gvndst = gvdst.create_group('nodes')
  gvedst.attrs['COUNT'] = gv['edges'].attrs['COUNT']
  gvndst.attrs['COUNT'] = gv['nodes'].attrs['COUNT']
  if( gv.attrs.__contains__('CLASS')):
      #special things if it is lattice based
      if (gv.attrs['CLASS']=='GRAPH'):
          gv.copy('lattice',gvdst)
          gvedst.attrs['CLASS']='GRAPH'
          for name in ['lattice_pos', 'roots']:
              gv['nodes'].copy(name, gvndst)
      if (gv.attrs['CLASS']=='REALWORLD'):
          gvedst.attrs['CLASS']='REALWORLD'
          for name in ['world_pos','roots', 'bc_node_index', 'bc_type', 'bc_value']:
              gv['nodes'].copy(name, gvndst)
          
  
  for name in ['radius', 'node_a_index', 'node_b_index']:
    gv['edges'].copy(name, gvedst)
  
  # then we recompute blood flow because the alorithm has changed and we may or may not want hematocrit
  if bloodflowparams['includePhaseSeparationEffect'] == True:
      pressure, flow, shearforce, hematocrit, flags = _ku.calc_vessel_hydrodynamics(gv, return_flags = True, bloodflowparams = bloodflowparams)
      gvedst.create_dataset('flow'      , data = flow       , compression = 9)
      gvedst.create_dataset('shearforce', data = shearforce , compression = 9)
      gvedst.create_dataset('hematocrit', data = hematocrit , compression = 9)
      gvedst.create_dataset('flags'     , data = flags      , compression = 9)
      gvndst.create_dataset('pressure'  , data = pressure   , compression = 9)
  else:
      print("Computing without phase seperation effect")
      pressure, flow, shearforce, hematocrit, flags = _ku.calc_vessel_hydrodynamics(gv, return_flags = True, bloodflowparams = bloodflowparams)
      gvedst.create_dataset('flow'      , data = flow       , compression = 9)
      gvedst.create_dataset('shearforce', data = shearforce , compression = 9)
      gvedst.create_dataset('hematocrit', data = hematocrit , compression = 9)
      gvedst.create_dataset('flags'     , data = flags      , compression = 9)
      gvndst.create_dataset('pressure'  , data = pressure   , compression = 9)

def computeAdaption_(gdst, vesselgroup, parameters):
  #myutils.buildLink(gdst, 'SOURCE_VESSELS', vesselgroup)
  #gdst.attrs['PARAMETER_CHECKSUM'] = myutils.checksum(parameters)
  #gdst.attrs['UUID']               = myutils.uuidstr()
  #writing parameters as acsii string. but i think it is better to store them as hdf datasets directly
  #ds = gdst.create_dataset('PARAMETERS_AS_JSON', data = np.string_(json.dumps(parameters, indent=2)))
  #myutils.hdf_write_dict_hierarchy(gdst, 'parameters', parameters)
  #gdst.file.flush()
  # now we can compute Adaption. The c++ routine can a.t.m. only read from hdf so it has to use our new file
  
  r = adaption_cpp.computeAdaption(vesselgroup, parameters['adaption'], parameters['calcflow'], gdst)
  myutils.hdf_write_dict_hierarchy(gdst, 'parameters', parameters)
  gdst.file.flush()
  
#  r = adaption.computeAdaption(vesselgroup, tumorgroup, parameters, None,gdst)

#def computeAdaption(f, group_path, parameters, cachelocation): 
#  if not isinstance(f, h5py.File):
#    f = h5files.open(f, 'r+', search = False)  # we open with r+ because the cachelocation might be this file so we need to be able to write to it
#    return computeAdaption(f, group_path, parameters, cachelocation) # recurse
#
#  #==== is this from a tumor sim, or just a vessel network? get hdf group objects =====#
#  group = f[group_path]
#  tumorgroup = None
#  vesselgroup1 = group
#  if 'vessels' in group:
#    vesselgroup1 = group['vessels']
#    if 'tumor' in group:
#      tumorgroup = group['tumor']
#
#  #====  this is for recomputing flow =====#
#  def read1(gmeasure, name):
#    gmeasure = gmeasure[name]
#    return myutils.H5FileReference(gmeasure.file.filename, gmeasure.name)
#
#  def write1(gmeasure, name):
#    gdst = gmeasure.create_group(name)
#    copyVesselnetworkAndComputeFlow(gdst, vesselgroup1, parameters.get("calcflow"))
#
#  #==== execute reading or computing and writing =====#
#  fm = h5files.open(cachelocation[0], 'a', search = False)
#  if 'addAverageConductivity' in parameters['adaption']:
#    avgCond = parameters['adaption'].pop('addAverageConductivity')
#  else:
#    avgCond = True
#  if avgCond:
#    parameters = AddAverageConductivity(vesselgroup1, parameters)
#  myutils.hdf_write_dict_hierarchy(fm['/'], 'parameters', parameters) 
#
#  #====  this is for Adaption =====#
#  def read2(gmeasure, name):
#    gmeasure = gmeasure[name]
#    return myutils.H5FileReference(gmeasure.file.filename, gmeasure.name)
#
#  def write2(gmeasure, name):
#    gdst = gmeasure.create_group(name)
#    myutils.buildLink(gdst, 'SOURCE', group)
#    gdst.attrs['CLASS'] = group.attrs['CLASS']
#    computeAdaption_(gdst, vesselgroup1, parameters)
#
#  #==== execute reading or computing and writing =====#
#  adaption_data_ref = myutils.hdf_data_caching(read2, write2, fm, 'adaption', (1,1))
#  #=== return filename and path to po2 data ====#
#  return adaption_data_ref


  
  #params['name'] = parameter_set_name
  
def doit_optimize_deap(individual):
  
  _ku.set_num_threads(1)
  if sys.flags.debug:
    print("individual in doit_optimize_deap")
    print(individual)
    print(individual.adaptionParameters['adaption'])
  
    print('starting doit in python')
  ''' update parameters '''
  individual.adaptionParameters['adaption'].update(
      k_c = individual[0],
      k_m = individual[1],
      k_s = individual[2],
      )
  returnState, mean = adaption_cpp.computeAdaption(individual.adaptionParameters['adaption'],individual.adaptionParameters['calcflow'], False)
  if 0: #hardcore debugging  
    print('mean: %f' % mean)
    print('param: %f' % individual.adaptionParameters['optimization']['desired_cap_flow'])
    print((mean-individual.adaptionParameters['optimization']['desired_cap_flow'])**2,)
  if returnState == 0:
    if sys.flags.debug:
      print("adaption succesful with mean: %f" % mean)
#  if not returnState == 0:
#    warnings.warn("adation broken", RuntimeWarning)
  
  if sys.flags.debug:
    print('mean: %f' % mean)
#  return (mean-individual.adaptionParameters['optimization']['desired_cap_flow'])**2,
  ''' I know it is wired, but now mean is variance 
      if this works out I have to rename it!!!!'''
  return mean,


def doit(parameters):
  _ku.set_num_threads(1)
  if sys.flags.debug:
    print(parameters)
#  fn = parameters['adaption']['vesselFileName']
#  pattern = parameters['adaption']['vesselGroupName']
  
  print('starting doit in python ... paramset :%s' % parameters['name'])

#  f_ = h5files.open(fn, 'r')
#  dirs = myutils.walkh5(f_['.'], pattern)
#  f_.close()
#  print(dirs)


#  for group_path in dirs:
    #cachelocation = (outfn_no_ext+'.h5', group_path+'_'+parameters_name)
    #cachelocation = (fnbase+'_adption_p_'+ parameters['name'] +'.h5', group_path)
    #ref = adaption_cpp.computeAdaption(f, group_path, parameters['adaption'],parameters['calcflow'], cachelocation)
  returnState, mean = adaption_cpp.computeAdaption(parameters['adaption'],parameters['calcflow'], True)
  if returnState == 0:
    print("adaption succesful with mean: %f" % mean)
  if not returnState == 0:
    warnings.warn("adation broken", RuntimeWarning)
#    print 'computed Adaption stored in:', ref
#    output_links.append(ref)
  return returnState, mean
  
def test():
    adaption_cpp.testAdaption()
    
def get_files_with_successful_adaption(filenames):
  good_files=[]
  for fn in filenames:
    with h5py.File(fn, 'r') as f:
      e = 'adaption/vessels_after_adaption' in f
      print e, 'in %s' % f.filename
      if e:
        good_files.append(fn)
  return good_files

def AddAverageConductivity(vesselgroup, parameters):
  dataman = myutils.DataManager(100, [DataBasicVessel()])
  vessels = dataman.obtain_data('vessel_graph', vesselgroup, ['flow', 'pressure', 'flags', 'radius','nodeflags'])
  def DoBC():
      conductivities, avgVenousPressure, avgArterialPressure, totalFlow = ComputeVascularTreeBloodFlowResistances(vessels)
      avgConductivity = (totalFlow/(avgArterialPressure-avgVenousPressure))
      print 'avgVenousPressure', avgVenousPressure
      print 'avgArterialPressure', avgArterialPressure
      print 'totalFlow', totalFlow
      print 'avgConductivity', avgConductivity
      return avgConductivity
  avgConductivity = DoBC()
  parameters['adaption']['avgRootNodeConductivity'] = avgConductivity
  return parameters
  
def getVesselTypes(vessel_groups):
  veins = []
  arteries = []
  capillaries = []
  if ('edges' in vessel_groups.keys()):# in fact it is only one group!
    g = vessel_groups    
    flags = np.array(_ku.read_vessels_from_hdf(g,['flags'])[1])
    circulated = np.bitwise_and(flags,_ku.CIRCULATED)
    circulated_indeces = np.nonzero(circulated)
    veins = np.bitwise_and(flags,_ku.VEIN)
    veins = veins[circulated_indeces]
    arteries = np.bitwise_and(flags,_ku.ARTERY)
    arteries = arteries[circulated_indeces]
    capillaries = np.bitwise_and(flags,_ku.CAPILLARY)
    capillaries = capillaries[circulated_indeces]
    veins = np.nonzero(veins)
    arteries =np.nonzero(arteries)
    capillaries = np.nonzero(capillaries)
  else:
    for g in vessel_groups:
      flags = _ku.read_vessels_from_hdf(g,['flags'])
      goods = np.sum(np.bitwise_and(flags[1],_ku.VEIN))
      veins.append(goods)
      goods = np.sum(np.bitwise_and(flags[1],_ku.ARTERY))
      arteries.append(goods)
      goods = np.sum(np.bitwise_and(flags[1],_ku.CAPILLARY))
      capillaries.append(goods)
    
  return veins, arteries, capillaries

## like worker_on_client but fixed filename
def do_simple_adaption( vfile_name, grp_name, paramset):
  #factory = getattr(parameterSetsAdaption, "pagmo_test")
  #vfile_name = "/localdisk/thierry/vessel_trees_better/my_chosen/PSO_data_vessels-large_2d-typeE-17x1L600-sample05_adption_p_human_guess.h5"
  #grp_name = "adaption/vessels_after_adaption"
  paramset['adaption'].update(
      vesselFileName = vfile_name,
      vesselGroupName = grp_name,
      )
  if sys.flags.debug:
    print(paramset)
  doit(paramset)