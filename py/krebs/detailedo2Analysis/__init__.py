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
#import os,sys
#from os.path import basename, dirname, join, splitext, commonprefix
import h5py
import h5files
import numpy as np
import extensions # for hdf5 support in np.asarray
import krebsutils
import myutils
#import posixpath
import math
import collections
#import pprint

#import matplotlib
#import matplotlib.cm
#import matplotlib.pyplot as pyplot
#import mpl_utils

#from krebs.analyzeGeneral import DataVesselGlobal, DataTumorTissueSingle, DataDistanceFromCenter, DataBasicVessel, DataVesselSamples, DataVesselRadial, BinsSpecRange, BinsSpecArray, obtain_distmap_, generate_samples, combineSamples, HdfCacheRadialDistribution, CalcPhiVessels, calc_distmap
#from krebs.analyzeBloodFlow import DataTumorBloodFlow
from krebs import detailedo2
from krebs import analyzeGeneral

#from myutils import f2l, f2s

def MakeSampleLocation(po2group):
  cachelocation = (po2group.file, po2group.name.lstrip('/')) #join('samples',po2group.name.lstrip('/').lstrip('po2').lstrip('/')))
  return cachelocation


def WriteSamplesToDisk(po2group):
  dataman = myutils.DataManager(5, [ DataDetailedPO2()])
  cachelocation = MakeSampleLocation(po2group)
  dataman.obtain_data('detailedPO2_total_fluxes', None, po2group, 30., None, cachelocation)


def getuuid_(h5grp):
  try:
    return h5grp.attrs['UUID']
  except KeyError:
    return ''

#readParameters(po2group)
class DataDetailedPO2(object):
  '''cache management'''
  keywords = [
    'detailedPO2Parameters', 'detailedPO2', 'detailedPO2_consumption', 'detailedPO2_samples', 'detailedPO2_total_fluxes', 'detailedPO2_radial', 'detailedPO2_global'
  ]

  def obtain_data(self, dataman, dataname, *args):
    if dataname == 'detailedPO2Parameters':
      po2group, = args
      return detailedo2.readParameters(po2group)
    
    if dataname == 'detailedPO2':
      po2group, = args
      a  = np.asarray(po2group['po2vessels'])
      if a.shape[0] <> 2: a = np.transpose(a)
      po2field  = po2group['po2field']
      ld = krebsutils.read_lattice_data_from_hdf(po2group['field_ld'])
      parameters = dataman.obtain_data('detailedPO2Parameters', po2group)
      return a, ld, po2field, parameters

    if dataname == 'detailedPO2_consumption':
      po2group, tumorgroup = args
      return detailedo2.computeO2Uptake(po2group, tumorgroup)

    if dataname == 'detailedPO2_samples' or dataname == 'detailedPO2_total_fluxes':
      prop, po2group, sample_length, every, cachelocation = args

      def read(gmeasure, name):
        gmeasure = gmeasure[name]
        if dataname == 'detailedPO2_samples':
          if prop == 'gtv':
            parameters = dataman.obtain_data('detailedPO2Parameters', po2group)
            #po2 = dataman.obtain_data('detailedPO2_samples','po2', *args[1:])
            #extpo2 = dataman.obtain_data('detailedPO2_samples','extpo2', *args[1:])
            #return (po2-extpo2)/(parameters['grid_lattice_const']*parameters.get('transvascular_ring_size',0.5))
            gvessels, gtumor = detailedo2.OpenVesselAndTumorGroups(po2group)
            jtv = dataman.obtain_data('detailedPO2_samples','jtv', *args[1:])
            #radius = dataman.obtain_data('basic_vessel_samples', 'radius', gvessels, sample_length)
            #radius = radius * (60.e-4*math.pi*2.* parameters['kD_tissue']*parameters['alpha_t'])
            radius = (60.e-4*parameters.get('D_tissue', parameters.get('kD_tissue', None))*parameters.get('solubility_tissue',parameters.get('alpha_t', None)))
            gtv = jtv/radius
            return gtv
          elif prop == 'sat':
            po2 = dataman.obtain_data('detailedPO2_samples','po2', *args[1:])
            parameters = dataman.obtain_data('detailedPO2Parameters', po2group)
            return detailedo2.PO2ToSaturation(po2, parameters)
          else:
            ds = gmeasure['smpl_'+prop]
            ds = ds[...] if every is None else ds[::every]
            return ds
        else:
          keys = filter(lambda k: k.startswith('flux_'), gmeasure.keys())
          fluxes = dict(map(lambda k: (k[5:], gmeasure[k][()]), keys))
          fluxes['e1'] = abs(100.*(fluxes['Jin_root']-fluxes['Jout_root']-fluxes['Jout_tv'])/fluxes['Jin_root'])
          fluxes['e2'] = abs(100.*(fluxes['Jin_root']-fluxes['Jout_root']-fluxes['Jout_cons']-fluxes.get('tv_cons',0.))/fluxes['Jin_root'])
          fluxes['e3'] = abs(100.*(fluxes['Jout_tv']-fluxes['Jout_cons']-fluxes.get('tv_cons',0.))/fluxes['Jout_tv'])
          return fluxes

      def write(gmeasure, name):
        # do the sampling
        gvessels, gtumor = detailedo2.OpenVesselAndTumorGroups(po2group)
        smpl, fluxes = detailedo2.sampleVessels(po2group, gvessels, gtumor, sample_length)
        #cons = dataman.detailedPO2_consumption(po2group, gtumor)
        cons = dataman.obtain_data('detailedPO2_consumption', po2group,gtumor)
        # get the raw data, just for the consumption flux
        po2vessels, ld, po2field, parameters = dataman.obtain_data('detailedPO2', po2group)
        avgCons = np.asarray(myutils.largeDatasetAverage(cons), dtype=np.float64)
        fluxes['Jout_cons'] = avgCons*np.product(cons.shape)*(ld.scale*1.e-4)**3
        del po2vessels, po2field, ld

        gmeasure = gmeasure.create_group(name)
        for k, v in smpl.iteritems():
          gmeasure.create_dataset('smpl_'+k, data = v, compression = 9, dtype = np.float32)
        for k, v in fluxes.iteritems():
          gmeasure.create_dataset('flux_'+k, data = v) # scalar dataset

      version_id = myutils.checksum(sample_length, 3, getuuid_(po2group))
      ret =  myutils.hdf_data_caching(read, write, cachelocation[0], (cachelocation[1],'samples_and_fluxes'), (None,version_id))
      return ret

    if dataname == 'detailedPO2_global':
      prop, po2group, sample_length, cachelocation = args
      samplelocation = MakeSampleLocation(po2group)

      def write(gmeasure, measurename):
        assert prop == measurename
        gvessels, gtumor = detailedo2.OpenVesselAndTumorGroups(po2group)
        if prop in ['po2','sat','gtv', 'jtv']:
          w = dataman.obtain_data('basic_vessel_samples', 'weight', gvessels, sample_length)
          d = dataman.obtain_data('detailedPO2_samples', prop, po2group, sample_length, 1, samplelocation)
          gmeasure.create_dataset(prop, data = myutils.WeightedAverageStd(d, weights=w))
        elif prop in ['sat_vein', 'sat_capi', 'sat_art']:
          w = dataman.obtain_data('basic_vessel_samples', 'weight', gvessels, sample_length)
          d = dataman.obtain_data('detailedPO2_samples', 'sat', po2group, sample_length, 1, samplelocation)
          f = dataman.obtain_data('basic_vessel_samples', 'flags', gvessels, sample_length)
          mask = ~myutils.bbitwise_and(f, krebsutils.WITHIN_TUMOR) & myutils.bbitwise_and(f, krebsutils.CIRCULATED)
          m = { 'sat_vein' : krebsutils.VEIN, 'sat_capi' : krebsutils.CAPILLARY, 'sat_art' : krebsutils.ARTERY }
          mask &= myutils.bbitwise_and(f, m[prop])
          d, w = d[mask], w[mask]
          gmeasure.create_dataset(prop, data = myutils.WeightedAverageStd(d, weights=w))
        elif prop in ['e1', 'e2', 'e3', 'Jin_root', 'Jout_root', 'Jout_tv', 'tv_cons', 'Jout_cons']:
          d = dataman.obtain_data('detailedPO2_total_fluxes', prop, po2group, sample_length, 1, samplelocation)
          gmeasure.create_dataset(prop, data = [d[prop], 0])
        elif prop == 'po2_tissue':
          _, po2ld, po2field, parameters = dataman.obtain_data('detailedPO2', po2group)
          d = myutils.largeDatasetAverageAndStd(po2field)
          gmeasure.create_dataset(prop, data = d)
        elif prop == 'mro2':
          uptakefield = detailedo2.computeO2Uptake(po2group, gtumor)
          d = myutils.largeDatasetAverageAndStd(uptakefield)
          gmeasure.create_dataset(prop, data = d)
        elif prop in ('sat_via_hb_ratio', 'vfhb_oxy', 'vfhb_deoxy', 'vfhb'):
          weight = dataman.obtain_data('basic_vessel_samples', 'weight', gvessels, sample_length)
          flags  = dataman.obtain_data('basic_vessel_samples', 'flags', gvessels, sample_length)
          mask = myutils.bbitwise_and(flags, krebsutils.CIRCULATED)
          hema = dataman.obtain_data('basic_vessel_samples', 'hematocrit', gvessels, sample_length)[mask]
          sat  = dataman.obtain_data('detailedPO2_samples' , 'sat', po2group, sample_length, None, samplelocation)[mask]
          rad  = dataman.obtain_data('basic_vessel_samples', 'radius', gvessels, sample_length)[mask]
          weight = weight[mask]
          hbvolume = weight*rad*rad*math.pi*hema
          ld = krebsutils.read_lattice_data_from_hdf(po2group['field_ld'])
          volume = np.product(ld.GetWorldSize())
          if prop == 'sat_via_hb_ratio':
            result = np.sum(hbvolume*sat) / np.sum(hbvolume)
          elif prop == 'vfhb_oxy':
            result = np.sum(hbvolume*sat)/volume
          elif prop == 'vfhb_deoxy':
            result = np.sum(hbvolume*(1.-sat))/volume
          elif prop == 'vfhb':
            result = np.sum(hbvolume)/volume
          gmeasure.create_dataset(prop, data = [result, 0.])
        elif prop in ('chb_oxy', 'chb_deoxy', 'chb'):
          m = { 'chb_oxy' : 'vfhb_oxy', 'chb_deoxy' : 'vfhb_deoxy', 'chb':'vfhb'}
          result = dataman.obtain_data('detailedPO2_global', m[prop], po2group, sample_length, cachelocation)
          result = result*detailedo2.chb_of_rbcs;
          gmeasure.create_dataset(prop, data = [result, 0.])
        elif prop == 'mro2_by_j':
          fluxes = dataman.obtain_data('detailedPO2_total_fluxes', prop, po2group, sample_length, 1)
          ld = krebsutils.read_lattice_data_from_hdf(po2group['field_ld'])
          worldbb = ld.worldBox
          result = fluxes['Jout_tv']/np.prod(worldbb[1]-worldbb[0])*1.e12
          gmeasure.create_dataset(prop, data = [result, 0.])
        elif prop == 'oef':
          fluxes = dataman.obtain_data('detailedPO2_total_fluxes', prop, po2group, sample_length, 1, samplelocation)
          result = (fluxes['Jin_root']-fluxes['Jout_root'])/fluxes['Jin_root']
          gmeasure.create_dataset(prop, data = [result, 0.])
        
        else:
          assert False
      def read(gmeasure, measurename):
        d = np.asarray(gmeasure[measurename])
        if measurename in ('chb_oxy', 'chb', 'chb_deoxy'):
          d *= 1.e6
        return d[0] # its a tuple (avg, std), we want avg now.
      version_num = collections.defaultdict(lambda : 3)
      version_id = myutils.checksum(sample_length, version_num[prop], getuuid_(po2group))
      #version_id = myutils.checksum(sample_length, (2 if prop in  else 1))
      return myutils.hdf_data_caching(read, write, cachelocation[0], ('global', cachelocation[1], prop), (1,1,version_id))

    if dataname == 'detailedPO2_radial':
      po2group, sample_length, bins_spec, distance_distribution_name, cachelocation = args
      samplelocation =  MakeSampleLocation(po2group)
      # we assume that there is a tumor. without this measurement makes little sense

      def read(gmeasure, name):
        d = dict(gmeasure[name].items())
        d = dict((k, myutils.MeanValueArray.read(v)) for k,v in d.items())
        hbo = d['vfhb_oxy']
        hbd = d['vfhb_deoxy']
        hb  = myutils.MeanValueArray(hbo.cnt, hbo.sum+hbd.sum, hbo.sqr+hbd.sqr)
        d['vfhb'] = hb
        sat = hbo.avg / hb.avg
        d['sat_via_hb_ratio'] = myutils.MeanValueArray(np.ones_like(hb.cnt), sat, sat*sat)
        d['chb_oxy'] = d['vfhb_oxy']*detailedo2.chb_of_rbcs*1.e6
        d['chb_deoxy'] = d['vfhb_deoxy']*detailedo2.chb_of_rbcs*1.e6
        d['chb'] = d['vfhb']*detailedo2.chb_of_rbcs*1.e6
        return d

      def write(gmeasure, name):
        gvessels, gtumor = detailedo2.OpenVesselAndTumorGroups(po2group)
        weight_smpl = dataman.obtain_data('basic_vessel_samples', 'weight', gvessels, sample_length)
        flags       = dataman.obtain_data('basic_vessel_samples', 'flags', gvessels, sample_length)
        # get teh radial distance function (either distance from tumor border or distance from center)
        dist_smpl, distmap, mask, tumor_ld   = dataman.obtain_data('distancemap_samples', gvessels, gtumor, sample_length, distance_distribution_name, None)
        # tumor_ld might actually be a unrelated lattice

        #filter uncirculated
        mask = mask & myutils.bbitwise_and(flags, krebsutils.CIRCULATED)
        dist_smpl = dist_smpl[mask]
        weight_smpl = weight_smpl[mask]

        bins = bins_spec.arange()
        gmeasure = gmeasure.create_group(name)

        for name in ['po2','extpo2','jtv','sat','gtv','dS_dx']:
          smpl = dataman.obtain_data('detailedPO2_samples', name, po2group, sample_length, None, samplelocation)
          myutils.MeanValueArray.fromHistogram1d(bins, dist_smpl, smpl[mask], w = weight_smpl).write(gmeasure, name)
        del smpl

        _, po2ld, po2field, parameters = dataman.obtain_data('detailedPO2', po2group)
        po2field = krebsutils.resample_field(np.asarray(po2field), po2ld.worldBox, tumor_ld.shape, tumor_ld.worldBox, order=1, mode='nearest')
        myutils.MeanValueArray.fromHistogram1d(bins, distmap.ravel(), po2field.ravel()).write(gmeasure, 'po2_tissue')
        del po2field

        uptakefield = detailedo2.computeO2Uptake(po2group, gtumor)
        uptakefield = krebsutils.resample_field(uptakefield, po2ld.worldBox, tumor_ld.shape, tumor_ld.worldBox, order=1, mode='nearest')
        myutils.MeanValueArray.fromHistogram1d(bins, distmap.ravel(), uptakefield.ravel()).write(gmeasure, 'mro2')
        del uptakefield

        hema = dataman.obtain_data('basic_vessel_samples', 'hematocrit', gvessels, sample_length)[mask]
        sat  = dataman.obtain_data('detailedPO2_samples' , 'sat', po2group, sample_length, None, samplelocation)[mask]
        rad  = dataman.obtain_data('basic_vessel_samples', 'radius', gvessels, sample_length)[mask]
        hbvolume = weight_smpl*rad*rad*math.pi*hema
        vol_per_bin = myutils.MeanValueArray.fromHistogram1d(bins, distmap.ravel(), np.ones_like(distmap.ravel())).cnt*(tumor_ld.scale**3)
        tmp = myutils.MeanValueArray.fromHistogram1d(bins, dist_smpl, hbvolume*sat)
        tmp.cnt = vol_per_bin.copy()
        tmp.write(gmeasure, 'vfhb_oxy')
        tmp = myutils.MeanValueArray.fromHistogram1d(bins, dist_smpl, hbvolume*(1.-sat))
        tmp.cnt = vol_per_bin.copy()
        tmp.write(gmeasure, 'vfhb_deoxy')
        del tmp, hbvolume, vol_per_bin

      version = getuuid_(po2group)
      ret = analyzeGeneral.HdfCacheRadialDistribution((read, write), 'po2', bins_spec, distance_distribution_name, cachelocation, version)
      return ret
    assert False




class DataDetailedO2Peff(object):
  """
   P = mass transfer coefficient /oxygen solubility in blood
   P = gamma/alpha=Sauerstoff-Diffusionskoeffizient* Nusselt number/Radius P=D_p*nu/r
   P_eff = P*beta/(1+beta)
   beta = sol. O2 (Plasma) / bound O2 (Hemoglobin)
   
   apparently beta/(1+beta) = c(O2, Plasma) / c(O2,Total), hier beta_factor genannt
  """
  keywords = [
    'detailedPO2_peff_samples', 'detailedPO2_peff_radial', 'detailedPO2_peffSrho_radial'
  ]

  def obtain_data(self, dataman, dataname, *args):
    if dataname == 'detailedPO2_peff_samples':
      po2group, sample_length, every, cachelocation = args
      gvessels, gtumor = detailedo2.OpenVesselAndTumorGroups(po2group)
      
      def read(gmeasure, name):
        ds = gmeasure[name]
        return ds[...] if every is None else ds[::every]
      
      def write(gmeasure, name):
        samplelocation = MakeSampleLocation(po2group)
        parameters = dataman.obtain_data('detailedPO2Parameters', po2group)
        rad  = dataman.obtain_data('basic_vessel_samples', 'radius', gvessels, sample_length)
        hema = dataman.obtain_data('basic_vessel_samples', 'hematocrit', gvessels, sample_length)
        po2  = dataman.obtain_data('detailedPO2_samples', 'po2', po2group, sample_length, 1, samplelocation)
        mtc  = detailedo2.computeMassTransferCoefficient(rad, parameters)
        blood_solubility = parameters.get('solubility_plasma', parameters.get('alpha_p', None))
        c_o2_total = detailedo2.PO2ToConcentration(po2, hema, parameters)
        c_o2_plasma = detailedo2.PO2ToConcentration(po2, np.zeros_like(hema), parameters) # to get the o2 conc. in plasma i just set the amount of RBCs to zero
        c_o2_plasma *= (1.0 - hema) # and reduce the concentration according to the volume fraction of plasma
        beta_factor = c_o2_plasma / c_o2_total
        peff = (1.0 / blood_solubility) * mtc * beta_factor
        
#        print 'po2', np.average(po2)
#        print 'mtc', np.average(mtc)
#        print 'hema', np.average(hema)
#        print 'blood_solubility', blood_solubility
#        print 'c_o2_total', np.average(c_o2_total)
#        print 'c_o2_plasma', np.average(c_o2_plasma)
#        print 'peff', np.average(peff)
#        print 'beta_factor', np.average(beta_factor)
        
        gmeasure.create_dataset(name, data = peff, compression = 9)
    
      version_id = myutils.checksum(sample_length, 1, getuuid_(po2group))
      ret =  myutils.hdf_data_caching(read, write, cachelocation[0], (cachelocation[1],'samples_and_fluxes', 'Peff'), (None,None,version_id))
      return ret
    
    elif dataname == 'detailedPO2_peff_radial':
      po2group, sample_length, bins_spec, distance_distribution_name, cachelocation = args
      # we assume that there is a tumor. without this measurement makes little sense

      def read(gmeasure, name):
        return myutils.MeanValueArray.read(gmeasure[name])

      def write(gmeasure, name):
        samplelocation = MakeSampleLocation(po2group)
        gvessels, gtumor = detailedo2.OpenVesselAndTumorGroups(po2group)
        smpl = dataman.obtain_data('detailedPO2_peff_samples', po2group, sample_length, None, samplelocation)
        
        data, = analyzeGeneral.GenerateRadialDistributions(dataman, gvessels, gtumor, sample_length, bins_spec, distance_distribution_name, None,
                                                           [(smpl, analyzeGeneral.radialAvgPerVessels)])
        data.write(gmeasure, name)

      version = myutils.checksum(2, getuuid_(po2group))
      ret = analyzeGeneral.HdfCacheRadialDistribution((read, write), 'Peff', bins_spec, distance_distribution_name, cachelocation, version)
      return ret
    
    # IDEA: store config stuff like sample length and distance distribution as member of the data handler for less call args
    # or store it in a common config instance.
    # Then allow to obtain a temporary handle on data which knows its GUID for a quick accurate check if data has changed.
    # The handle probably would have to obtain the GUID from disk but this is it.
    elif dataname == 'detailedPO2_peffSrho_radial':
      po2group, sample_length, bins_spec, distance_distribution_name, cachelocation = args

      def read(gmeasure, name):
        return myutils.MeanValueArray.read(gmeasure[name])

      def write(gmeasure, name):
        samplelocation = MakeSampleLocation(po2group)
        gvessels, gtumor = detailedo2.OpenVesselAndTumorGroups(po2group)        
        peff = dataman.obtain_data('detailedPO2_peff_samples', po2group, sample_length, None, samplelocation)
        rad  = dataman.obtain_data('basic_vessel_samples', 'radius', gvessels, sample_length)        
        peff = peff * math.pi*2.*rad        
        
        data, = analyzeGeneral.GenerateRadialDistributions(dataman, gvessels, gtumor, sample_length, bins_spec, distance_distribution_name, None,
                                                           [(peff, analyzeGeneral.radialAvgPerVolume)])
        data.write(gmeasure, name)

      version = myutils.checksum(2, getuuid_(po2group)) # this should actually depend on the samples
      ret = analyzeGeneral.HdfCacheRadialDistribution((read, write), 'PeffSrho', bins_spec, distance_distribution_name, cachelocation, version)
      return ret

    assert False


def ObtainOxygenExtractionFraction(dataman, po2group, cachelocation):
    def read(gmeasure, groupname):
      return dict((k, float(v[...])) for (k,v) in gmeasure[groupname].iteritems())

    def write(gmeasure, groupname):
      gvessels, gtumor = detailedo2.OpenVesselAndTumorGroups(po2group)
      parameters = dataman.obtain_data('detailedPO2Parameters', po2group)
      vessels = dataman.obtain_data('vessel_graph', gvessels, ['flow', 'pressure', 'position', 'flags', 'hematocrit'])
      pos = vessels['position']
      press = vessels['pressure']
      flow = vessels['flow']
      flags = vessels['flags']
      hema = vessels['hematocrit']
      po2vessels  = np.asarray(po2group['po2vessels'])
      hema = np.column_stack((hema, hema))
      conc = detailedo2.PO2ToConcentration(po2vessels, hema, parameters)
      #sat  = detailedo2.PO2ToSaturation(po2vessels, parameters)
      conc_diff = np.abs(conc[:,1]-conc[:,0])
      conc      = np.average(conc, axis=1)
      flow_diff = 0.5*flow*conc_diff
      flow      = flow*conc
      roots = set(gvessels['nodes/roots'][...])
      del hema
      del conc
      del conc_diff
      del po2vessels

      if gtumor:
        ldtumor = dataman.obtain_data('ld', gtumor.file)
        dist = dataman.obtain_data('fieldvariable', gtumor.file, 'theta_tumor', gtumor.name)
        dist = krebsutils.sample_field(pos, dist, ldtumor, linear_interpolation=True)
        dist = dist - 0.5 # is greater 0 inside the tumor?!

      total_flow_in, total_flow_out = 0., 0.
      flow_in, flow_out = 0., 0.
      flow_in_err, flow_out_err, total_flow_in_err, total_flow_out_err = 0., 0., 0., 0.
      for i, (a,b) in enumerate(vessels.edgelist):
        if not (flags[i] & krebsutils.CIRCULATED): continue
        if gtumor and dist[a]<0 and dist[b]>0: # b is in the tumor
          if press[a]<press[b]:
            flow_out += flow[i]
            flow_out_err += flow_diff[i]
          else:
            flow_in += flow[i]
            #print 'tumor inflow of sat', sat[i]
            flow_in_err += flow_diff[i]
        elif gtumor and dist[a]>0 and dist[b]<0: # a is in the tumor
          if press[a]>press[b]:
            flow_out += flow[i]
            flow_out_err += flow_diff[i]
          else:
            flow_in += flow[i]
            #print 'tumor inflow of sat', sat[i]
            flow_in_err += flow_diff[i]
        if (a in roots or b in roots):
            if flags[i] & krebsutils.ARTERY:
              #print 'total inflow of sat', sat[i]
              total_flow_in += flow[i]
              total_flow_in_err += flow_diff[i]
            elif flags[i] & krebsutils.VEIN:
              #print 'total inflow of sat', sat[i]
              total_flow_out += flow[i]
              total_flow_out_err += flow_diff[i]

      Err = lambda a,b,da,db: abs(1.0/a - (a-b)/a/a)*da + abs(1.0/a)*db

      res = dict(
        tumor_o2_in  = flow_in, 
        tumor_o2_out = flow_out, 
        total_o2_out = total_flow_out, 
        total_o2_in  = total_flow_in,
        oef_total = (total_flow_in-total_flow_out)/total_flow_in,
        oef_tumor = (flow_in-flow_out)/flow_in,
        flow_in_err = flow_in_err, 
        flow_out_err = flow_out_err,
        total_flow_in_err = total_flow_in_err,
        total_flow_out_err = total_flow_out_err,
        oef_total_err = Err(total_flow_in, total_flow_out, total_flow_in_err, total_flow_out_err),
        oef_err       = Err(flow_in, flow_out, flow_in_err, flow_out_err)
      )
      #print 'computed oef: '
      #pprint.pprint(res)
      g = gmeasure.create_group(groupname)
      for k,v in res.iteritems():
        g.create_dataset(k, data = v)
    ret = myutils.hdf_data_caching(read, write, cachelocation[0], ('global', cachelocation[1], 'oxygen_extraction'), (None, None, 8,))
    return ret

O2DataHandlers = [DataDetailedPO2, DataDetailedO2Peff]