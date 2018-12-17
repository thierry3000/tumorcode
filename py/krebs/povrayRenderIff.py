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
import os,sys
from os.path import basename, dirname, join, splitext, commonprefix
if __name__=='__main__': sys.path.append(join(dirname(__file__),'..'))
import h5py
import os,sys
import numpy as np
import extensions # for hdf5 support in np.asarray
import krebsutils
import myutils
import posixpath
from copy import deepcopy
import math
#from mpl_utils import *

#from krebs.povrayRenderVessels import  addVesselTree
from krebs.povrayEasy import *
from krebs.detailedo2 import PO2ToSaturation, OpenVesselAndTumorGroups, chb_of_rbcs
from krebs.detailedo2Analysis import DataDetailedPO2
from krebs.analyzeGeneral import DataBasicVessel

import matplotlib
import matplotlib.cm
#import matplotlib.pyplot as pyplot #not good on cluster without graphics!!!


colors = [(0, 0.1,0,0.4), (0.7, 0.8,0,0.05), (1., 0.9, 0.7, 0.7) ]
#cm_po2 = matplotlib.colors.LinearSegmentedColormap('', {
#  'red'   : [(x, r, r) for (x,r,g,b) in colors],
#  'green' : [(x, g, g) for (x,r,g,b) in colors],
#  'blue'  : [(x, b, b) for (x,r,g,b) in colors],
#}, N = 256, gamma = 1.)
cm_po2 = matplotlib.cm.jet


def InsertGraphColors(vesselgraph, ifffield, data_name, automatic_scale=False, max_conc=1.0):
  edges = vesselgraph.edgelist
  num_nodes = len(vesselgraph.nodes['position'])

#  if data_name in vesselgraph.edges:
#    edgedata = data = vesselgraph.edges[data_name]
#    nodedata = krebsutils.edge_to_node_property(num_nodes, edges, data, 'avg')
#  else:
#    nodedata = data = vesselgraph.nodes[data_name]
#    edgedata = np.average((data[edges[:,0]], data[edges[:,1]]), axis=0)
  #note everything is uncirculated anyway
  #edgedata = data = vesselgraph.edges['radius']
  #nodedata = krebsutils.edge_to_node_property(num_nodes, edges, data, 'avg')

  if data_name == 'auc':
    #p1 = np.amax(data)
    #if po2field is not None:
    #  p1 = max(p1, np.amax(po2field))
    #p1 = math.ceil(p1/10.0)*10.0  # round to powers of something
    
    #p1 = 100.0
    #value_range = (np.amin(ifffield), np.amax(ifffield))
    color_norm = matplotlib.colors.LogNorm().autoscale(ifffield)
    cm = matplotlib.cm.ScalarMappable(cmap = 'Reds', norm=color_norm)
    #value_range = (0, 37000)
    value_range = (np.amin(ifffield), np.amax(ifffield))
    if value_range[0]==0 and value_range[1] == 0:
      value_range = (-1,1)
    #cm = matplotlib.cm.ScalarMappable(cmap = 'Reds')
    cm.set_clim(*value_range)
    #cm = matplotlib.colors.LinearSegmentedColormap.from_list('my_cmap',[(0,0,0,0),'red'],256)
    #cm = plt.get_cmap('Reds')
  elif data_name == 'iff_pressure':
    cm = matplotlib.cm.ScalarMappable(cmap = matplotlib.cm.spectral)
    value_range = (np.amin(ifffield), np.amax(ifffield))
    cm.set_clim(*value_range)
  elif data_name == 'drug_conc':
    cm = matplotlib.cm.ScalarMappable(cmap = matplotlib.cm.Spectral, norm=matplotlib.colors.LogNorm())
    if automatic_scale:
      value_range = (np.amin(ifffield), np.amax(ifffield))
    else:
      value_range= (1e-3,0.5e0)
    # new global variable
    if not max_conc == 1.0:
      value_range = (1e-3, max_conc)
    else:
      value_range = (np.amin(ifffield), np.amax(ifffield))
    print(value_range)
    if value_range[0]==0 and value_range[1] == 0:
      value_range = (1,2)
    #value_range = (0, 0.42)
    #color_norm = matplotlib.colors.LogNorm(*value_range)
    #color_norm = matplotlib.colors.LogNorm().autoscale(ifffield)
    #cm = matplotlib.cm.ScalarMappable(cmap = 'Reds', norm=color_norm)
    #cm = matplotlib.cm.ScalarMappable(cmap = 'Reds', norm=color_norm)
    #if np.amin(ifffield) <0:
    #  value_range = (0, np.amax(ifffield))
    #else:
    #value_range = (np.amin(ifffield), np.amax(ifffield))
    #value_range = (1e-7, 0.5)
    
    #value_range = (0, 0.42)
    
    #value_range = (0.01, 0.3)
    cm.set_clim(*value_range)
    #cm.autoscale(ifffield)
    #cm.LogNorm()
  else:
    cm = matplotlib.cm.ScalarMappable(cmap = matplotlib.cm.gnuplot)
    p1 = math.ceil(np.amax(data))
    value_range = (0., p1)
    
  

  #colors = lambda arr: np.power(cm.to_rgba(arr)[:,:3], 2.4)
  edgedata = vesselgraph.edges['radius']
  nodedata = krebsutils.edge_to_node_property(num_nodes, edges, edgedata, 'avg')
  
  if 1:
    print('custom cm:')
    mycm = matplotlib.colors.Colormap('Reds')
    
    theEdgeColors = np.ones((len(edgedata),3))
    theNodeColors = np.ones((len(nodedata),3))
    blue=matplotlib.colors.colorConverter.to_rgba('blue')
    normalizedEdgeData = matplotlib.colors.Normalize(vmin=np.amin(edgedata),vmax=np.amax(edgedata))(edgedata)
    normalizedNodeData = matplotlib.colors.Normalize(vmin=np.amin(edgedata),vmax=np.amax(edgedata))(nodedata)
    #normalizedEdgeData = matplotlib.colors.LogNorm(vmin=np.amin(edgedata),vmax=np.amax(edgedata))(edgedata)
    #normalizedNodeData = matplotlib.colors.LogNorm(vmin=np.amin(edgedata),vmax=np.amax(edgedata))(nodedata)
    #normalizedEdgeData = matplotlib.colors.LogNorm()(edgedata)
    #normalizedNodeData = matplotlib.colors.LogNorm()(nodedata)
    for (no, intensity) in enumerate(normalizedEdgeData):
      theEdgeColors[no,:] = np.multiply(blue[0:3],intensity)
    for (no, intensity) in enumerate(normalizedNodeData):
      theNodeColors[no,:] = np.multiply(blue[0:3],intensity)
    vesselgraph.edges['colors'] = theEdgeColors
    vesselgraph.nodes['colors'] = theNodeColors
  return cm


def renderSliceWithDistribution((vessel_ld, vessel_graph, data_name), (volume_ld, volumedata), imagefn, label, options, max_conc=None):
  wbbox = volume_ld.worldBox
  options.wbbox=wbbox
  trafo = calc_centering_normalization_trafo(wbbox)
  volume_ld = trafo.transform_ld( volume_ld)
  vessel_ld = trafo.transform_ld( vessel_ld)
  cm = InsertGraphColors(vessel_graph, volumedata, data_name, automatic_scale = options.auto_colorscale, max_conc=max_conc)  

  def DoTheRendering(fn, options):
    with EasyPovRayRender(options) as epv:
      epv.setBackground(options.background)
      central_position = wbbox[4]+planeZCoord/100.*(wbbox[5]-wbbox[4])
      central_position = central_position*trafo.w
      
      cam_fov = 60.
      cam_distance_factor = ComputeCameraDistanceFactor(cam_fov, options.res, wbbox)
      epv.setCamera((0,0,cam_distance_factor*1.0), lookat = (0,0,0), fov = cam_fov, up = 'y')
  
      epv.addLight(10.*Vec3(1,0.5,2), 1.2)
      options.vessel_clip=('zslice', -150*trafo.w, +150*trafo.w)      
  
      pvcm = matplotlibColormapToPovray('DATACOLORMAP', cm)
      epv.declareColorMap(pvcm)

      if not options.not_render_volume:
        epvvol = epv.declareVolumeData(volumedata, volume_ld.GetWorldBox())
        epv.addVolumeDataSlice(epvvol, (0,0,planeZCoord), (0, 0, 1), pvcm)        
        
      if not options.not_render_vessels:
        epv.addVesselTree2(epv, vessel_graph, trafo = trafo, options=options)   
      CallPovrayAndOptionallyMakeMPLPlot(epv, cm, label, options)
  
  planeZCoord = 0
  DoTheRendering('slice_at_%0.1f_percent'%planeZCoord+imagefn, options)
    

def renderSlice((vessel_ld, vessel_graph, data_name), (volume_ld, volumedata), imagefn, label, kwargs):
  kwargs = deepcopy(kwargs)

  wbbox = vessel_ld.worldBox
  trafo = calc_centering_normalization_trafo(wbbox)

  cm = InsertGraphColors(vessel_graph, volumedata, data_name)
  
  with EasyPovRayRender(**kwargs) as epv:
    epv.setBackground(kwargs.pop('background',1.0))
    cam_fov = 60.
    cam_distance_factor = ComputeCameraDistanceFactor(cam_fov, kwargs['res'], wbbox)
    epv.setCamera((0,0,cam_distance_factor*1.05), lookat = (0,0,0), fov = cam_fov, up = 'y')
    epv.addLight(10.*Vec3(1,0.5,2), 1.2)
    if (wbbox[1]-wbbox[0]) < (wbbox[5]-wbbox[4])*2.:
      kwargs.update(vessel_clip =('zslice', -300*trafo.w, +300*trafo.w))

    addVesselTree2(epv, vessel_graph, trafo = trafo, **kwargs)  
    
    CallPovrayAndOptionallyMakeMPLPlot(epv, imagefn, cm, label, **kwargs)


def renderVasculatureWTumor((vessel_ld, vessel_graph, data_name), gtumor, imagefn, label, kwargs):
  kwargs = deepcopy(kwargs)

  wbbox = vessel_ld.worldBox
  trafo = calc_centering_normalization_trafo(wbbox)

  cm = InsertGraphColors(vessel_graph, None, data_name)
  
  with EasyPovRayRender(**kwargs) as epv:
    epv.setBackground(kwargs.pop('background',1.0))
    cam_fov = 60.
    cam_distance_factor = ComputeCameraDistanceFactor(cam_fov, kwargs['res'], wbbox)
    epv.setCamera((0,0,cam_distance_factor*1.05), lookat = (0,0,0), fov = cam_fov, up = 'y')
    epv.addLight(10.*Vec3(1,0.5,2), 3.)
    if (wbbox[1]-wbbox[0]) < (wbbox[5]-wbbox[4])*2.:
      kwargs.update(vessel_clip =('zslice', -500*trafo.w, +500*trafo.w))
      kwargs.update(tumor_clip =('zslice', -250*trafo.w, 250*trafo.w))

    epv.addVesselTree2(epv, vessel_graph, trafo = trafo, **kwargs)
    
    # totally ad-hoc BS
    # only for simple sphere model for now
    tumRadius = gtumor.attrs['TUMOR_RADIUS']
    shellHalfThickness = 10 # microns
    style = '''
      texture {
        pigment {
          color rgb<1,0.5, 0.1>
        }
        finish {
          specular 0.1
          ambient 1.
        }
      }'''
    t1 = pv.Sphere((0,0,0), (tumRadius-shellHalfThickness)*trafo.w,
        style
    )
    t2 = pv.Sphere((0,0,0), (tumRadius+shellHalfThickness)*trafo.w,
        style
    )
    o = pv.Difference(t2, t1)
    clipStyle = "pigment { color rgb<%0.2f,%0.2f,%0.2f> }" % defaultCrossSectionColor
    clip = clipFactory(kwargs.get('tumor_clip', None))
    s = krebsutils.povray_clip_object_str(clip, clipStyle)
    if s:
      o = pv.Intersection(o, s)
    o.write(epv.pvfile)
    
    CallPovrayAndOptionallyMakeMPLPlot(epv, imagefn, cm, label, **kwargs)




def renderScene(drug_grp, imagefn, options): 
  imagefn, ext = splitext(imagefn)
  ext = '.' + options.format
  
  options.imageFileName = imagefn+'_iff_drug'+ext
  
  max_conc = getMaxConcentration(drug_grp.file)    
  dataman = myutils.DataManager(2, [DataBasicVessel()])
  
  timepoint = drug_grp.attrs['time'] #comes in seconds
  timepoint = timepoint/3600.
  #gvessels = drug_grp.parent['iff/vessels']
  gvessels = drug_grp.parent['vessels']
  iff_pressure_field = drug_grp.parent['iff/iff_pressure']
  drug_conc_field = drug_grp['conc']
  cell_drug_conc_field = drug_grp['conc_cell']
  ex_drug_conc_field = drug_grp['conc_ex']
  
  #ex_drug_auc_field = drug_grp.parent['measurements/drug_local_integral']['auc_ex']
  #in_drug_auc_field = drug_grp.parent['measurements/drug_local_integral']['auc_in']
  
  #iff_ld  = krebsutils.read_lattice_data_from_hdf_by_filename(drug_grp.parent['field_ld'])
  iff_ld  = krebsutils.read_lattice_data_from_hdf_by_filename(str(drug_grp.file.filename),str(drug_grp.parent.name)+'field_ld')
  
  #ld = iffgroup['lattice']

  #po2vessels, po2field_ld, po2field, parameters = dataman('detailedPO2', po2group)
  #po2vessels = np.average(po2vessels, axis=0)
  #print 'po2vessels:', po2vessels.min(), po2vessels.max()
  print 'ifpfield:', np.amin(iff_pressure_field), np.amax(iff_pressure_field)
  print 'drug_conc_field:', np.amin(drug_conc_field), np.amax(drug_conc_field)

  #vessel_ld = krebsutils.read_lattice_data_from_hdf(gvessels['lattice'])
  vessel_ld = krebsutils.read_lattice_data_from_hdf_by_filename(str(gvessels.file.filename),str(gvessels.name)+'/lattice')
  vessel_graph = dataman('vessel_graph', gvessels, ['position', 'flags', 'radius'])  
    
  #vessel_graph.edges['po2vessels'] = po2vessels
  #vessel_graph.edges['saturation'] = PO2ToSaturation(po2vessels, parameters)
  #vessel_graph.edges['hboconc'] = vessel_graph.edges['saturation']*vessel_graph.edges['hematocrit']*chb_of_rbcs*1.0e3
  vessel_graph = vessel_graph.get_filtered(edge_indices = myutils.bbitwise_and(vessel_graph['flags'], krebsutils.CIRCULATED))
  if options.filterradiuslowpass>0.0:
    print("lowpass filter activated:")
    vessel_graph = vessel_graph.get_filtered(edge_indices = vessel_graph['radius']< kwargs['filterradiuslowpass'])
  imagefn, ext = splitext(imagefn)
  ext = '.' + options.format
  
  if 1:
    if timepoint==0:
      renderSliceWithDistribution((vessel_ld, vessel_graph, 'iff_pressure'), (iff_ld, iff_pressure_field), (imagefn+'_iff_pressure_t%0.1fh'%timepoint )+ext, 'IF pressure t=%.1f h'%timepoint, options, max_conc=max_conc)
    renderSliceWithDistribution((vessel_ld, vessel_graph, 'drug_conc'), (iff_ld, drug_conc_field), (imagefn+'_iff_drug_t%0.1fh'%timepoint )+ext, '%s t=%.1f h'%(drug_grp.file.attrs.get('MESSAGE'),timepoint), options, max_conc=max_conc)
    #renderSliceWithDistribution((vessel_ld, vessel_graph, 'drug_conc'), (iff_ld, cell_drug_conc_field), (imagefn+'_iff_drug_incell_t%0.1fh'%timepoint)+ext, 'Tr. intr. t=%.1f h'%timepoint, kwargs)
    #renderSliceWithDistribution((vessel_ld, vessel_graph, 'drug_conc'), (iff_ld, ex_drug_conc_field), (imagefn+'_iff_drug_excell_t%0.1fh'%timepoint)+ext, 'Tr. extr. t=%.1f h'%timepoint, kwargs)
    if options.plot_auc:
      renderSliceWithDistribution((vessel_ld, vessel_graph, 'auc'), (iff_ld, ex_drug_auc_field), imagefn+'_iff_ex_drug_auc'+ext, 'Tr. extr. t=%.1f h'%timepoint, options, max_conc=max_conc)
    #renderSliceWithDistribution((vessel_ld, vessel_graph, 'auc'), (iff_ld, in_drug_auc_field), imagefn+'_iff_in_drug_auc'+ext, 'Tr. intr. t=%.1f h'%timepoint, kwargs)
  if 0:
    renderSlice((vessel_ld, vessel_graph, 'radius'), (iff_ld, iff_pressure_field), imagefn+'_iff_pressure'+ext, '', options)
    renderSlice((vessel_ld, vessel_graph, 'radius'), (iff_ld, drug_conc_field), imagefn+'_iff_drug'+ext, '', options)
    renderSlice((vessel_ld, vessel_graph, 'radius'), (iff_ld, cell_drug_conc_field), imagefn+'_iff_drug_incell'+ext, '', options)
    renderSlice((vessel_ld, vessel_graph, 'radius'), (iff_ld, ex_drug_conc_field), imagefn+'_iff_drug_excell'+ext, '', options)
  
  #renderSlice((vessel_ld, vessel_graph, 'saturation'), (None, None), imagefn+'_saturation'+ext, '', kwargs)
  #renderSlice((vessel_ld, vessel_graph, 'hboconc'), (None, None), imagefn+'_hboconc'+ext, 'HbO [mmol/l blood]', kwargs)
  #renderVasculatureWTumor((vessel_ld, vessel_graph, 'po2vessels'), gtumor, imagefn+'_po2vt'+ext, '', kwargs)


def getMaxConcentration(f):
  if 'max_conc' not in f.attrs.keys():
    max_conc=0.0
    paths = myutils.walkh5(f['.'], 'out*')
    for path in paths:
      if 'conc' in f[path]:
        current_max = np.amax(f[path+'/conc'])
        if current_max > max_conc:
          max_conc=current_max
    print('found max conc over time: %f'%max_conc)
    #f.attrs.create('max_conc',max_conc)
  return max_conc
    
  

