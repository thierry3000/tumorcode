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
import h5files
import os,sys
import numpy as np
import extensions # for hdf5 support in np.asarray
import krebsutils
import myutils
import posixpath
from copy import deepcopy
import math

from krebs.povrayRenderVessels import  addVesselTree
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


def InsertGraphColors(vesselgraph, po2field, data_name):
  edges = vesselgraph.edgelist
  num_nodes = len(vesselgraph.nodes['position'])

  if data_name in vesselgraph.edges:
    edgedata = data = vesselgraph.edges[data_name]
    nodedata = krebsutils.edge_to_node_property(num_nodes, edges, data, 'avg')
  else:
    nodedata = data = vesselgraph.nodes[data_name]
    edgedata = np.average((data[edges[:,0]], data[edges[:,1]]), axis=0)

  if data_name == 'po2vessels':
    try:
      p1 = np.amax(data)
    except ValueError:
      print("p1 not found")
      pass
    if po2field is not None:
      p1 = max(p1, np.amax(po2field))
    p1 = math.ceil(p1/10.0)*10.0  # round to powers of something
    #p1 = 100.0
    value_range = (0., p1)
    cm = matplotlib.cm.ScalarMappable(cmap = cm_po2)
  elif data_name == 'saturation':
    cm = matplotlib.cm.ScalarMappable(cmap = matplotlib.cm.spectral)
    value_range = (0,1.)
  elif data_name == 'hboconc':
    cm = matplotlib.cm.ScalarMappable(cmap = matplotlib.cm.gnuplot)
    p1 = math.ceil(np.amax(data))
    value_range = (0., p1)
  cm.set_clim(*value_range)

  colors = lambda arr: np.power(cm.to_rgba(arr)[:,:3], 2.4)
  if data_name in vesselgraph.edges:
    edgecolors = colors(data)
    nodecolors = colors(nodedata)
  else:
    edgecolors = colors(edgedata)
    nodecolors = colors(data)

  flags = vesselgraph.edges['flags']
  nflags = krebsutils.edge_to_node_property(num_nodes, edges, flags, 'or')
  is_not_set = lambda flags_,flag: np.bitwise_not(np.asarray(np.bitwise_and(flags_, flag), np.bool))
  gray = np.asarray((0.3,0.3,0.3))
  uncirculated = is_not_set(flags,krebsutils.CIRCULATED)
  nuncirculated = is_not_set(nflags,krebsutils.CIRCULATED)
  edgecolors[uncirculated] = gray
  nodecolors[nuncirculated] = gray

  print 'colormap range ', cm.get_clim()

  vesselgraph.edges['colors'] = edgecolors
  vesselgraph.nodes['colors'] = nodecolors
  return cm


#def CallPovrayAndOptionallyMakeMPLPlot(epv, imagefn, cm, label, **kwargs):
#    if kwargs.get('overlay', True):
#      RenderImageWithOverlay(epv, imagefn, cm, label, **kwargs)      
#    else:
#      epv.render(imagefn)


def renderSliceWithDistribution((vessel_ld, vessel_graph, data_name), (volume_ld, volumedata), imagefn, label, options):
  #kwargs = deepcopy(kwargs_in)
  wbbox = volume_ld.worldBox
#  kwargs['wbbox']=wbbox
  options.wbbox = wbbox
  trafo = calc_centering_normalization_trafo(wbbox)
  volume_ld = transform_ld(trafo, volume_ld)
  vessel_ld = transform_ld(trafo, vessel_ld)
  cm = InsertGraphColors(vessel_graph, volumedata, data_name)  

  def DoTheRendering(fn, options):
    with EasyPovRayRender(options) as epv:
      epv.setBackground(options.background)
  
      cam_fov = 60.
      cam_distance_factor = ComputeCameraDistanceFactor(cam_fov, options.res, wbbox)
      epv.setCamera((0,0,cam_distance_factor*1.05), lookat = (0,0,0), fov = cam_fov, up = 'y')
  
      epv.addLight(10.*Vec3(1,0.5,2), 1.2)
      options.vessel_clip=('zslice', -150*trafo.w, +150*trafo.w)
      #kwargs.update(vessel_clip =('zslice', -150*trafo.w, +150*trafo.w))
  
      pvcm = matplotlibColormapToPovray('DATACOLORMAP', cm)
      epv.declareColorMap(pvcm)
#      print vessel_ld
#      print volume_ld
#      print vessel_ld.worldBox
#      print volume_ld.worldBox
      if options.render_volume:
        epvvol = epv.declareVolumeData(volumedata, volume_ld.GetWorldBox())
        epv.addVolumeDataSlice(epvvol, (0,0,planeZCoord), (0, 0, 1.), pvcm)
      if options.render_vessels:
        addVesselTree(epv, vessel_graph, trafo = trafo, options=options)
#      if kwargs.get('render_volume_', True):
#        epvvol = epv.declareVolumeData(volumedata, volume_ld.GetWorldBox())
#        epv.addVolumeDataSlice(epvvol, (0,0,planeZCoord), (0, 0, 1.), pvcm)
#      if kwargs.get('render_vessels_', True):
#        addVesselTree(epv, vessel_graph, trafo = trafo, **kwargs)
        
      CallPovrayAndOptionallyMakeMPLPlot(epv, fn, cm, label, options)
  
#  if options.projection_plot:
  if False:
    planeZCoord = vessel_ld.worldBox[4]+0.1*vessel_ld.scale
    basefn, ext = splitext(imagefn)
    tf1 = mkstemp.File(suffix='.png', prefix='mwpov_', text=False, keep=True)
    tf2 = mkstemp.File(suffix='.png', prefix='mwpov_', text=False, keep=True)
    fn1 = tf1.filename
    fn2 = tf2.filename
    kwargs['overlay'] = False
    kwargs['render_volume_'] = False
    DoTheRendering(fn1, kwargs)
    kwargs['render_volume_'] = True
    kwargs['render_vessels_'] = False
    DoTheRendering(fn2, kwargs)
    import subprocess
    if kwargs_in.get('overlay', True):
      tf3 = mkstemp.File(suffix='.png', prefix='mwpov_', text=False, keep=True)
      subprocess.call(['convert', '-page', '+0+0', fn2, '-page', '+0+0', fn1, '-flatten', tf3.filename])
      plotsettings = dict(myutils.iterate_items(kwargs, ['dpi','fontcolor'], skip=True))
      OverwriteImageWithColorbar(tf3.filename, cm, label, output_filename = imagefn, **plotsettings)
    else:
      subprocess.call(['convert', '-page', '+0+0', fn2, '-page', '+0+0', fn1, '-flatten', imagefn])
    return
  else:
    planeZCoord = 0.
    DoTheRendering(imagefn, options)


def renderSlice((vessel_ld, vessel_graph, data_name), (volume_ld, volumedata), imagefn, label, options):
  #kwargs = deepcopy(kwargs)

  wbbox = vessel_ld.worldBox
  options.wbbox = wbbox
  #kwargs['wbbox']=wbbox
  trafo = calc_centering_normalization_trafo(wbbox)

  cm = InsertGraphColors(vessel_graph, volumedata, data_name)
  
  with EasyPovRayRender(options) as epv:
    epv.setBackground(options.background)
    cam_fov = 60.
    cam_distance_factor = ComputeCameraDistanceFactor(cam_fov, options.res, wbbox)
    epv.setCamera((0,0,cam_distance_factor*1.05), lookat = (0,0,0), fov = cam_fov, up = 'y')
    epv.addLight(10.*Vec3(1,0.5,2), 1.2)
    if (wbbox[1]-wbbox[0]) < (wbbox[5]-wbbox[4])*2.:
      options.vessel_clip=('zslice', -300*trafo.w, +300*trafo.w)
      #kwargs.update(vessel_clip =('zslice', -300*trafo.w, +300*trafo.w))

    addVesselTree(epv, vessel_graph, trafo = trafo, options=options)  
    
    CallPovrayAndOptionallyMakeMPLPlot(epv, imagefn, cm, label, options)


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

    addVesselTree(epv, vessel_graph, trafo = trafo, **kwargs)
    
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




def renderScene(po2group, imagefn, options):
  #kwargs = myutils.updated(default_parameters, kwargs)

  dataman = myutils.DataManager(2, [DataDetailedPO2(), DataBasicVessel()])

  gvessels, gtumor = OpenVesselAndTumorGroups(po2group)

  po2vessels, po2field_ld, po2field, parameters = dataman('detailedPO2', po2group)
  po2vessels = np.average(po2vessels, axis=0)
  print 'po2vessels:', po2vessels.min(), po2vessels.max()
  print 'po2field:', np.amin(po2field), np.amax(po2field)

  vessel_ld = krebsutils.read_lattice_data_from_hdf(gvessels['lattice'])
  vessel_graph = dataman('vessel_graph', gvessels, ['position', 'flags', 'radius', 'hematocrit'])  
    
  vessel_graph.edges['po2vessels'] = po2vessels
  vessel_graph.edges['saturation'] = PO2ToSaturation(po2vessels, parameters)
  vessel_graph.edges['hboconc'] = vessel_graph.edges['saturation']*vessel_graph.edges['hematocrit']*chb_of_rbcs*1.0e3
  vessel_graph = vessel_graph.get_filtered(edge_indices = myutils.bbitwise_and(vessel_graph['flags'], krebsutils.CIRCULATED))
  if options.filterradiuslowpass>0:
    print("lowpass filter activated:")
    vessel_graph = vessel_graph.get_filtered(edge_indices = vessel_graph['radius']< options.filterradiuslowpass)
    filenamepostfix = '_rlp'
#  if 'filterradiuslowpass' in kwargs.keys():
#    if kwargs['filterradiuslowpass'] >0:
#      print("lowpass filter activated:")
#      vessel_graph = vessel_graph.get_filtered(edge_indices = vessel_graph['radius']< kwargs['filterradiuslowpass'])
#      filenamepostfix = '_rlp'

  imagefn, ext = splitext(imagefn)
  ext = '.' + options.format
  renderSliceWithDistribution((vessel_ld, vessel_graph, 'po2vessels'), (po2field_ld, po2field), imagefn+'_po2vessels'+ext, '', options)
  renderSlice((vessel_ld, vessel_graph, 'saturation'), (None, None), imagefn+'_saturation'+ext, '', options)
  renderSlice((vessel_ld, vessel_graph, 'hboconc'), (None, None), imagefn+'_hboconc'+ext, 'HbO [mmol/l blood]', options)
  #renderVasculatureWTumor((vessel_ld, vessel_graph, 'po2vessels'), gtumor, imagefn+'_po2vt'+ext, '', kwargs)


#def doit(fn, pattern, parameters = dict()):
#  f = h5files.open(fn, 'r+')
#  paths = myutils.walkh5(f['.'], pattern)
#  for path in paths:
#    po2group = f[path]
#    imagefn = '-'.join([splitext(basename(f.filename))[0], po2group.attrs.get('SOURCE_PATH','').strip(posixpath.sep).replace(posixpath.sep,'-')])+'.'+parameters.pop('format', 'png')
#    renderScene(po2group, imagefn, parameters)


#if __name__ == '__main__':
#  import optparse #Note: Deprecated since version 2.7: The optparse module is deprecated and will not be developed further; development will continue with the argparse module.
#  parser = optparse.OptionParser()
#  parser.add_option("--dpi", dest="dpi", default=None, action="store")
#  parser.add_option("--format", dest="format", default=None, action="store")
#  parser.add_option("--no-overlay", dest="overlay", default=True, action="store_false")
#  options, args = parser.parse_args()  
#  
#  krebsutils.set_num_threads(4)
#
#  import povrayRenderSettings
#  settings = deepcopy(povrayRenderSettings.image)
#  settings.update(povrayRenderSettings.tumor)
#
#  if options.dpi: 
#    settings['dpi'] = float(options.dpi)
#  if options.format:
#    settings['format'] = options.format
#  settings['overlay'] = options.overlay    
#  
#  settings['projection_plot'] = False
#    
#  filenames, pattern = args[:-1], args[-1]
#  for fn in filenames:
#    doit(fn, pattern, settings)
