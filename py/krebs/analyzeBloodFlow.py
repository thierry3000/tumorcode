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
if __name__ == '__main__':
  import os.path, sys
  sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'..'))
  
import os,sys
from os.path import basename, splitext
import h5py
import numpy as np
import extensions # for asarray with h5py support
import krebsutils
import posixpath
import math
import glob
import collections
from pprint import pprint
import itertools

import matplotlib
import identifycluster
if (identifycluster.getname()=='snowden' or identifycluster.getname()=='durga'):
  matplotlib.use('agg')
import matplotlib.pyplot as pyplot
import mpl_utils

import krebs.quantities
import myutils
import krebs.analyzeGeneral


def ComputeIsosurfaceBloodFlow(dataman, vesselgroup, nodalLevel, level):
    vessels = dataman.obtain_data('vessel_graph', vesselgroup, ['flow', 'pressure', 'flags'])
    result = krebsutils.SumIsoSurfaceIntersectionWithVessels(level, vessels.edgelist, vessels['pressure'], vessels['flags'], nodalLevel, vessels['flow'])    
    return dict(flow_in = result[0], flow_out = result[1])


def ComputeIsosurfaceRegionalBloodFlow(dataman, vesselgroup, nodalLevel, ld, distancemap, level):
    data = ComputeIsosurfaceBloodFlow(dataman, vesselgroup, nodalLevel, level)    
    volume = np.count_nonzero(distancemap.ravel() < level)
    volume *= ld.scale**3
    #volume = math.pow(level, 3.) *4./3.*math.pi
    return data['flow_in'] / volume
    
#### some debug output
#        loc_in = np.asarray(loc_in)
#        loc_out = np.asarray(loc_out)
#        fig, ax = pyplot.subplots(1,1)
#        ax.plot(loc_in[:,0], loc_in[:,1], lw = 0, ms=3, marker = 'x', color = 'r')
#        ax.plot(loc_out[:,0], loc_out[:,1], lw = 0, ms=3, marker = 'x', color = 'b')
#        fig.text(0.1,0.1,'total_flow: %f |  %f' % (total_flow_in, total_flow_out))
#        fig.text(0.2,0.1,'tumor_flow: %f |  %f' % (flow_in, flow_out))
#        pyplot.show()


def ComputeIsosurfaceAvgRadius(dataman, vesselgroup, nodalLevel, ld, distancemap, level):
    vessels = dataman.obtain_data('vessel_graph', vesselgroup, ['flags','radius', 'pressure'])
    radiusSumIn, radiusSumOut = krebsutils.SumIsoSurfaceIntersectionWithVessels(level, vessels.edgelist, vessels['pressure'], vessels['flags'], nodalLevel, vessels['radius'])
    countIn, countOut = krebsutils.SumIsoSurfaceIntersectionWithVessels(level, vessels.edgelist, vessels['pressure'], vessels['flags'], nodalLevel, np.ones(vessels.edgelist.shape[0], dtype = np.float64))    
    return radiusSumIn/countIn


#def ComputeIsosurfaceAvgFlow(dataman, vesselgroup, nodalLevel, ld, distancemap, level):
#    vessels = dataman.obtain_data('vessel_graph', vesselgroup, ['flags','flow'])
#    data = vessels['flow']
#    flags = vessels['flags']
#
#    dataList = []
#    for i, (a,b) in enumerate(vessels.edgelist):
#      if not (flags[i] & krebsutils.CIRCULATED): continue
#      if nodalLevel[a]<level and nodalLevel[b]>level: # b is in the tumor
#        dataList.append(data[i])
#      elif nodalLevel[a]>level and nodalLevel[b]<level: # a is in the tumor
#        dataList.append(data[i])
#    res = np.average(dataList)
#    return res


def ComputeSphereVesselDensity(dataman, vesselgroup, nodalLevel, ld, distancemap, level):
    vessels = dataman.obtain_data('vessel_graph', vesselgroup, ['flags','radius', 'pressure'])
    countIn, countOut = krebsutils.SumIsoSurfaceIntersectionWithVessels(level, vessels.edgelist, vessels['pressure'], vessels['flags'], nodalLevel, np.ones(vessels.edgelist.shape[0], dtype = np.float64))
    surfaceArea = 4.*math.pi*level**2
    return (countIn + countOut) / surfaceArea


def ComputeIsosurfaceRadialCurve(dataman, vesselgroup, tumorgroup, bins_spec, distance_distribution_name, ld, cachelocation, WorkerFunction):
    print 'compute', WorkerFunction.__name__, vesselgroup.name
    vessels = dataman.obtain_data('vessel_graph', vesselgroup, ['position', 'flags'])
    pos = vessels['position']
    distancemap, ld = krebs.analyzeGeneral.obtain_distmap_(dataman, tumorgroup, distance_distribution_name, ld)
    dist = krebsutils.sample_field(pos, distancemap, ld, linear_interpolation=True)
    bounds = np.amin(distancemap), np.amax(distancemap)
    bins = bins_spec.arange()
    results = []
    for level in bins:
      if level < bounds[0] or level > bounds[1]:
        data = np.nan
      else:
        data = WorkerFunction(dataman, vesselgroup, dist, ld, distancemap, level)
        print 'level %f, data %f' % (level, data)
      results.append(data)
    return results, bins


def GetRootVesselData(vessels, data):
    flags = vessels['flags']
    roots = set(vessels.roots)
    arterialData = []
    venousData = []
    for i, (a,b) in enumerate(vessels.edgelist):
      if not (flags[i] & krebsutils.CIRCULATED): continue
      if (a in roots or b in roots):
          if flags[i] & krebsutils.ARTERY:
            arterialData.append(data[i])
          elif flags[i] & krebsutils.VEIN:
            venousData.append(data[i])
    return arterialData, venousData


@myutils.UsesDataManager
def ComputeRootFlowInfo(dataman, vesselgroup, cachelocation):
    def write(gmeasure, groupname):
      vessels = dataman.obtain_data('vessel_graph', vesselgroup, ['flags', 'flow', 'radius'])
      edgeData = zip(vessels['flow'], vessels['radius'])      
      arterialData, venousData = GetRootVesselData(vessels, edgeData)
      arterialFlow, arterialRadius = zip(*arterialData)
      venousFlow, venousRadius = zip(*venousData)
      d = dict(
        totalFlow         = np.sum(arterialFlow),
        arterialCount     = len(arterialData),
        venousCount       = len(venousData),
        avgArterialRadius = np.average(arterialRadius),
        avgVenousRadius   = np.average(venousRadius),
      )
      g = gmeasure.create_group(groupname)
      for k, v in d.iteritems():
        g.attrs[k] = v
    
    def read(gmeasure, groupname):
      return dict(gmeasure[groupname].attrs.items())
    
    return myutils.hdf_data_caching(read, write, cachelocation[0], ('global', cachelocation[1], 'rootNodeFlowData'), (None, None, 1))


@myutils.UsesDataManager
def ComputeTumorDistanceMapAndVolumeForPerfusion_(dataman, tumorGroup):
  distancemap, ld = krebs.analyzeGeneral.obtain_distmap_(dataman, tumorGroup, 'levelset')
  volume = np.count_nonzero(distancemap.ravel() < 0)
  volume *= ld.scale**3
  return distancemap, ld, volume


@myutils.UsesDataManager
def ComputeIsoTumorSpherePerfusion(dataman, vesselgroup, tumorGroup, cachelocation):
  def write(gmeasure, groupname):
    distancemap, ld, volume = ComputeTumorDistanceMapAndVolumeForPerfusion_(dataman, tumorGroup)
    vessels = dataman.obtain_data('vessel_graph', vesselgroup, ['position', 'flags'])
    distSamples = krebsutils.sample_field(vessels['position'], distancemap, ld, linear_interpolation=True)
    result = ComputeIsosurfaceBloodFlow(dataman, vesselgroup, distSamples, 0)
    result = result['flow_in']        
    result /= volume
    ds = gmeasure.create_dataset(groupname, data = result)
    ds.attrs['tumorGroup'] = str(tumorGroup)
    ds.attrs['volume']     = volume
    ds.attrs['unit']       = '1 / s'
  
  def read(gmeasure, groupname):
    return gmeasure[groupname]

  return myutils.hdf_data_caching(read, write, cachelocation[0], ('global', cachelocation[1], 'isoTumorSphereRBF'), (None, None, 1))


@myutils.UsesDataManager
def ComputeSystemPerfusion(dataman, vesselgroup, cachelocation):
  def write(gmeasure, groupname):
    vessels = dataman.obtain_data('vessel_graph', vesselgroup, ['flags', 'flow'])
    arterialFlow, _ = GetRootVesselData(vessels, vessels['flow'])
    arterialFlow = np.sum(arterialFlow)
    ldvessels = krebsutils.read_lattice_data_from_hdf(vesselgroup['lattice'])
    totalVolume = np.cumprod(ldvessels.GetWorldSize())[2]
    perfusion = arterialFlow / totalVolume
    ds = gmeasure.create_dataset(groupname, data = perfusion)
    ds.attrs['unit'] = '1 / s'

  def read(gmeasure, groupname):
    return gmeasure[groupname]  
  
  return myutils.hdf_data_caching(read, write, cachelocation[0], ('global', cachelocation[1], 'systemRBF'), (None, None, 1))


@myutils.UsesDataManager
def ComputeIsoTumorSpherePerfusionScaleFactor(dataman, vesselGroup, tumorVesselGroup, tumorGroup, cachelocation, cachelocationTumor):
  def write(gmeasure, groupname):
    systemPerfusion   = ComputeSystemPerfusion(dataman, vesselGroup, cachelocation)
    isoTumorPerfusion = ComputeIsoTumorSpherePerfusion(dataman, vesselGroup, tumorGroup, cachelocation)
    rescaledPerfusion = systemPerfusion[...] / isoTumorPerfusion[...]
    ds = gmeasure.create_dataset(groupname, data = rescaledPerfusion)
    ds.attrs['unit'] = '1'
 
  def read(gmeasure, groupname):
    return gmeasure[groupname]
    
  return myutils.hdf_data_caching(read, write, cachelocation[0], ('global', cachelocation[1], 'perfusion_scaling_factor'), (None, None, 1))


@myutils.UsesDataManager
def ComputeIsoTumorSphereRescaledPerfusion(dataman, vesselGroup, tumorVesselGroup, tumorGroup, cachelocation, cachelocationTumor):
  def write(gmeasure, groupname):
    tumorPerfusion    = ComputeIsoTumorSpherePerfusion(dataman, tumorVesselGroup, tumorGroup, cachelocationTumor)
    scalingFactor     = ComputeIsoTumorSpherePerfusionScaleFactor(dataman, vesselGroup, tumorVesselGroup, tumorGroup, cachelocation, cachelocationTumor)
    rescaledPerfusion = scalingFactor[...] * tumorPerfusion[...]
    ds = gmeasure.create_dataset(groupname, data = rescaledPerfusion)
    ds.attrs['unit'] = tumorPerfusion.attrs['unit']
  
  def read(gmeasure, groupname):
    return gmeasure[groupname]
    
  return myutils.hdf_data_caching(read, write, cachelocation[0], ('global', cachelocation[1], 'scaled_RBF'), (None, None, 1))


## ---------------- ------- -----------------------
## obtain tumor flow flow in in/out flow per tissue volume
## ---------------- ------- -----------------------
class DataTumorBloodFlow(object):
  keywords = [
    'blood_flow', 'blood_flow_rbf', 'blood_flow_resistances', 'cum_rbf_radial', 'avg_surf_vessel_rad_radial', 'sphere_vessel_density'
  ]

  @staticmethod
  def FixUnit_((name, data)):
    if 'rBF' in name: data = data*60. # 1/s -> 1/min
    if 'flow' in name: data = data*60.*1.e-12 # um/s -> ml/min
    if 'volume' in name: data = data*1.e-9
    return (name, data)

  @staticmethod
  def ComputeTotalBloodFlow_(vessels):
      arterialFlow, venousFlow = GetRootVesselData(vessels, vessels['flow'])
      res = dict(total_flow_out = np.sum(venousFlow), total_flow_in = np.sum(arterialFlow))
      return res

  def obtain_data(self, dataman, dataname, *args):
    if dataname == 'blood_flow':
      vesselgroup, tumorgroup, cachelocation = args
      has_tumor = tumorgroup is not None
      
      def read(gmeasure, groupname):
        return dict((k, float(v[...])) for (k,v) in gmeasure[groupname].iteritems())

      def write(gmeasure, groupname):
        #vessels = krebsutils.read_vesselgraph(vesselgroup, ['flow', 'pressure', 'position', 'flags'])
        vessels = dataman.obtain_data('vessel_graph', vesselgroup, ['position', 'flags', 'flow'])
        pos = vessels['position']
        
        if has_tumor:
          ldtumor = dataman.obtain_data('ld', tumorgroup.file)
          dist = dataman.obtain_data('fieldvariable', tumorgroup.file, 'theta_tumor', tumorgroup.name)
          dist = krebsutils.sample_field(pos, dist, ldtumor, linear_interpolation=True)
          res = ComputeIsosurfaceBloodFlow(dataman, vesselgroup, dist, 0.5)
          res.update(
            DataTumorBloodFlow.ComputeTotalBloodFlow_(vessels)
          )
        else:
          res = DataTumorBloodFlow.ComputeTotalBloodFlow_(vessels)

        g = gmeasure.create_group(groupname)
        for k,v in res.iteritems():
          g.create_dataset(k, data = v)

      #fm = myutils.MeasurementFile(f, h5files)
      ret = myutils.hdf_data_caching(read, write, cachelocation[0], ('global', cachelocation[1], 'tissue', 'blood_flow'), (1, 1, 1, 3))
      return ret

    if dataname == 'blood_flow_rbf':
      vesselgroup, tumorgroup, cachelocation = args

      has_tumor = tumorgroup is not None
      data      = dataman.obtain_data('blood_flow', vesselgroup, tumorgroup, cachelocation).copy()
      ldvessels = krebsutils.read_lattice_data_from_hdf(vesselgroup['lattice'])
      total_flow = data['total_flow_in']
      total_volume = np.cumprod(ldvessels.GetWorldSize())[2]
      total_flow_p_volume = total_flow/total_volume
      data['rBF_total'] = total_flow_p_volume
      data['total_volume'] = total_volume
      if has_tumor:
        ldtumor = dataman.obtain_data('ld', tumorgroup.file)
        theta_tumor = dataman.obtain_data('fieldvariable', tumorgroup.file, 'theta_tumor', tumorgroup.name)
        tumor_volume = np.sum(theta_tumor)*(ldtumor.scale**3)
        tumor_flow = data['flow_in']
        tumor_flow_p_volume = tumor_flow/tumor_volume
        data['rBF_tumor'] = tumor_flow_p_volume
        data['tumor_volume'] = tumor_volume
        #print 'estimated tumor volume:', tumor_volume
        #print 'tumor flow:', tumor_flow
        #print 'rBF:', tumor_flow_p_volume*60.
      data = dict(map(DataTumorBloodFlow.FixUnit_, data.items()))
      return data
    
    if dataname in ('cum_rbf_radial', 'avg_surf_vessel_rad_radial', 'sphere_vessel_density'):
      vesselgroup, tumorgroup, bins_spec, distance_distribution_name, ld, cachelocation = args
      WorkerFunction = {
        'cum_rbf_radial' : ComputeIsosurfaceRegionalBloodFlow,
        'avg_surf_vessel_rad_radial':ComputeIsosurfaceAvgRadius,
        'sphere_vessel_density': ComputeSphereVesselDensity,
      }[dataname]
      version = {
        'cum_rbf_radial': 4,
        'avg_surf_vessel_rad_radial': 3,
        'sphere_vessel_density': 1,
      }[dataname]
      
      def read(gmeasure, groupname):
        gmeasure = gmeasure[groupname]
        return gmeasure['values'], gmeasure['bins']
      
      def write(gmeasure, groupname):
        values, bins = ComputeIsosurfaceRadialCurve(dataman, vesselgroup, tumorgroup, bins_spec, distance_distribution_name, ld, cachelocation, WorkerFunction)
        gmeasure = gmeasure.create_group(groupname)
        gmeasure.create_dataset('values', data = values)
        gmeasure.create_dataset('bins', data = bins)
      return krebs.analyzeGeneral.HdfCacheRadialDistribution((read,write), dataname, bins_spec, distance_distribution_name, cachelocation, version)


## ---------------- ------- -----------------------
## util for average blood flow
## ---------------- ------- -----------------------
def obtain_averaged_blood_flow(dataman, files, group, cachelocation):
  l = collections.defaultdict(list)
  for f in files:
    vesselgroup   = f[group]
    tumorgroup    = krebs.analyzeGeneral.try_find_tumor_group_from_vesselgroup(vesselgroup)
    rdata = dataman.obtain_data('blood_flow_rbf', vesselgroup, tumorgroup, cachelocation(f[group])).copy()
    data = dataman.obtain_data('blood_flow', vesselgroup, tumorgroup, cachelocation(f[group])).copy()
    for k, v in data.iteritems():
      l[k].append(v)
    for k, v in rdata.iteritems():
      l[k].append(v)
  for k, v in l.iteritems():
    l[k] = np.average(v), np.std(v)
  return l


def obtain_radial_curves(dataman, groups, cachelocation):
  '''tumor groups required. TODO better finding of groups'''
  curves = []
  binspec = krebs.analyzeGeneral.BinsSpecRange(500., 5000., 100.)
  for group in groups:  
    vesselgroup = group
    tumorgroup  = tumorgroup = krebs.analyzeGeneral.try_find_tumor_group_from_vesselgroup(vesselgroup)
    data, bins = dataman.obtain_data('cum_rbf_radial', vesselgroup, tumorgroup, binspec, 'radial', None, cachelocation(group))
    fuck, _    = dataman.obtain_data('avg_surf_vessel_rad_radial', vesselgroup, tumorgroup, binspec, 'radial', None, cachelocation(group))
    density, _ = dataman.obtain_data('sphere_vessel_density', vesselgroup, tumorgroup, binspec, 'radial', None, cachelocation(group))
    curves.append(
      (np.asarray(data), np.asarray(bins), np.asarray(fuck), np.asarray(density))
    )
  return curves


def GetTimeLabel(vesselgroup):
  try:
    tumorgroup = krebs.analyzeGeneral.try_find_tumor_group_from_vesselgroup(vesselgroup)
    return 't = %i h' % int(tumorgroup.parent.attrs['time'])
  except:
    return vesselgroup.name

# note: use insertVesselConfigInO2File.py script (in scripts repository) to copy vesselfile message over. It is the type name in it.
def GetVesselTypeLabel(group):
  import re
  try:
    msg = group.file.attrs['VESSELFILE_MESSAGE']
  except KeyError:
    return group.name
  m = re.search('type(\w)', msg)
  if not m:
    return 'unkown'
  m = m.group(1)
  # specific naming scheme for paper
  from detailedo2Analysis.plotsForPaper import RewriteVesselLabel
  return RewriteVesselLabel(m)


## ---------------- ------- -----------------------
## ---------------- ------- -----------------------
if __name__ == '__main__':
  filenames = sys.argv[1:-1]
  pattern   = sys.argv[-1]
    
  files = [ h5py.File(fn, 'r') for fn in filenames ]
  
  dataman = myutils.DataManager(100, [krebs.analyzeGeneral.DataTumorTissueSingle(), 
                                      krebs.analyzeGeneral.DataBasicVessel(), 
                                      krebs.analyzeGeneral.DataDistanceFromCenter(),
                                      DataTumorBloodFlow()])


  groupnames = myutils.walkh5(files[0], pattern, return_h5objects=False)
  allgroups = list(itertools.chain.from_iterable((f[g] for g in groupnames) for f in files))

  # determination of storage file, copied from analyzeVesselsBulkTumor
  prefix, suffix = myutils.splitcommonpresuffix(map(lambda s: basename(s), filenames))
  outputbasename, _ = splitext(prefix+suffix)
  fn_measure = 'common-radial-cache.h5'
  f_measure = h5files.open(fn_measure, 'a')
  def cachelocation(g):
    path = posixpath.join('FileCS_'+myutils.checksum(basename(g.file.filename)), g.name.strip(posixpath.sep))
    return (f_measure, path)

  with mpl_utils.PdfWriter('bloodFlow-%s.pdf' % outputbasename) as pdfpages:
    def bloodflow(group):
      bv = obtain_averaged_blood_flow(dataman, files, group, cachelocation)
      print 'of ', group, 'bf = ', bv
      def printit(id, label, unit):
        s = r'$%s$ = $%s \pm %s$ [$%s$]' % (
          label, myutils.f2s(bv[id][0],latex=True), myutils.f2s(bv[id][1], latex=True), unit)
        return s
      s = '\n'.join([
        GetTimeLabel(files[0][group]),
        #printit('flow_in', r'\tilde{BF}_{tumor}', '\mu m^3\,Blood\,/\,min'),
        printit('total_flow_in', r'BF_{total}', '\mu m^3\,Blood\,/\,min'),
        printit('total_flow_out', r'BF_{total}', '\mu m^3\,Blood\,/\,min'),
        #printit('rBF_tumor', r'\tilde{rBF}_{tumor}', 'ml\,Blood\,/\,(ml\,Tissue\,min)'),
        printit('rBF_total', r'rBF_{total}', 'ml\,Blood\,/\,(ml\,Tissue\,min)'),
        printit('total_volume', r'vol', 'ml')
      ])
      return s

    sall = [ bloodflow(g) for g in groupnames ]
    fig = pyplot.figure(figsize = (8,6))
    fig.text(0.1, 0.5, '\n'.join(sall))
    pdfpages.savefig(fig)
    
    colors = [
       '#e50000', #red
       '#9a0eea', #violet
       '#0343df', #blue
       '#f97306', #orange
       '#677a04', #olive green
       '#ceb301', #mustard
       '#04d8b2', #aquamarine
       '#06470c', #forest greeen
       '#840000', #dark red
       '#607c8e', #blue grey
    ]  
    
    if 0:
      fig, axes = pyplot.subplots(3,1, figsize = mpl_utils.a4size*np.array([0.5, 0.6]))
      groupsByTime = collections.defaultdict(list)
      for g in allgroups:
        groupsByTime[g.name].append(g)
      dataByTime = []
      for groupname, groups in groupsByTime.iteritems():
        curves = obtain_radial_curves(dataman, groups, cachelocation)
        curves = np.average(curves, axis = 0)
        data, bins, fuck, density = curves
        dataByTime.append((groupname, data))
        data = 60.*data # 1/s -> 1/min
        density *= 1.e6
        axes[0].plot(bins, data, label = GetTimeLabel(groups[0]))
        axes[1].plot(bins, fuck)
        axes[2].plot(bins, density)
      dataByTime = sorted(dataByTime, key = lambda (t, v): t)
      if len(dataByTime)>1:
        ax0twin = axes[0].twinx()
        ax0twin.plot(bins, dataByTime[1][1] / dataByTime[0][1], color = 'k', label = 'ratio')
        ax0twin.set(ylabel = r'$\tilde{rBF}(tumor)/\tilde{rBF}(normal)$')
        ax0twin.legend(loc = mpl_utils.loc.upper_left)
      axes[0].set(ylabel = r'$\tilde{rBF}$ [$ml / ml / min$]')
      axes[1].set(ylabel = r'$Vascular Radius\quad<r_v>$')
      axes[2].set(ylabel = r'$Density\quad n / S$ [$mm^{-2}$]', xlabel = '$r$ [$\mu m$]')
      axes[0].legend()
      pdfpages.savefig(fig)
      #pyplot.show()
    
    if 0:
      groupsByTime = collections.defaultdict(list)
      for g in allgroups:
        groupsByTime[g.name].append(g)
      for groupname, groups in groupsByTime.iteritems():
        fig, ax = pyplot.subplots(1,1, figsize = mpl_utils.a4size*np.array([0.5, 0.2]))
        curves = obtain_radial_curves(dataman, groups, cachelocation)
        for data, bins, _, _ in curves:
          data = 60.*data # 1/s -> 1/min
          ax.plot(bins, data, color = (0.5,0.5,0.5,0.3))
        data, bins, _, _ = np.average(curves, axis = 0)
        data = 60.*data
        ax.plot(bins, data, color = 'k', label = 'average')
        ax.legend()
        ax.set(title = GetTimeLabel(groups[0]), ylabel = r'$\tilde{rBF}$ [$ml / ml / min$]', xlabel = '$r$ [$\mu m$]')
        pdfpages.savefig(fig)

    if 1:
      groupsByTimeAndType = collections.defaultdict(lambda: collections.defaultdict(list))
      for g in allgroups:
        groupsByTimeAndType[g.name][GetVesselTypeLabel(g)].append(g)
        
    if 0:
      for timeName, byType in sorted(groupsByTimeAndType.items(), key = lambda (a,b): a):
        fig, ax = pyplot.subplots(1,1, figsize = mpl_utils.a4size*np.array([0.5, 0.2]))
        for typeName, groups in sorted(byType.items(), key = lambda (a,b): a):
          curves = obtain_radial_curves(dataman, groups, cachelocation)
          curves = np.average(curves, axis=0)
          data, bins, fuck, density = curves
          data = 60.*data # 1/s -> 1/min
          ax.plot(bins, data, label = typeName)
        ax.legend()
        ax.set(title = GetTimeLabel(groups[0]), ylabel = r'$\tilde{rBF}$ [$ml / ml / min$]', xlabel = '$r$ [$\mu m$]')
        pdfpages.savefig(fig)
    
    if 1:
      for timeName, byType in sorted(groupsByTimeAndType.items(), key = lambda (a,b): a):
        fig, axes = pyplot.subplots(3, 1, figsize = mpl_utils.a4size*np.array([0.5, 0.6]))
        for i, (typeName, groups) in enumerate(sorted(byType.items(), key = lambda (a,b): a)):
          flowInfos = [ ComputeRootFlowInfo(dataman, group, cachelocation(group)) for group in groups ]
          flowInfos = myutils.zipListOfDicts(flowInfos) # -> dicts of numpy arrays
          print timeName, typeName#,'\n', flowInfos
          axes[0].plot(flowInfos['arterialCount'], flowInfos['totalFlow'], colors[i], lw = 0., marker = '>', markeredgecolor = colors[i], label = typeName)
          axes[1].plot(flowInfos['avgArterialRadius'], flowInfos['totalFlow'], colors[i], lw = 0., marker = '>', markeredgecolor = colors[i], label = typeName)
          axes[2].plot(flowInfos['arterialCount'], flowInfos['avgArterialRadius'], colors[i], lw = 0., marker = '>', markeredgecolor = colors[i], label = typeName)
        axes[0].set(xlabel = 'arterial node count', ylabel = r'$q$ total $[\mu m^3/s]$', xscale = 'log')
        axes[2].legend(bbox_to_anchor=(0.85, 1.1), loc=2, borderaxespad=0.)
        axes[1].set(xlabel = 'arterial $<r>$', ylabel = r'$q$ total $[\mu m^3/s]$')
        axes[2].set(xlabel = 'arterial node count', ylabel = 'arterial $<r>$', xscale = 'log')
        fig.suptitle(GetTimeLabel(groups[0]))
        pyplot.tight_layout()
        pdfpages.savefig(fig)
    
    if 1:
      tommHg = krebs.quantities.kPa.asNumber(krebs.quantities.mmHg)
      x = np.linspace(-200., 200., 400)
      y = map(lambda x: tommHg*krebsutils.PressureRadiusRelation(x, False), x)
      fig, ax = pyplot.subplots(1, 1, figsize = mpl_utils.a4size*np.array([0.5, 0.3]))
      ax.plot(x, y)
      ax.set(xlabel = 'vessel radius $r$ [$\mu m$]', ylabel = 'pressure boundary condition $p^{(BC)}$ [$mmHg$]')
      ax.yaxis.grid(True, which='major', color = '#f0f0f0', ls = '-')
      ax.xaxis.grid(True, which='major', color = '#f0f0f0', ls = '-')
      pyplot.tight_layout()
      pdfpages.savefig(fig)
      