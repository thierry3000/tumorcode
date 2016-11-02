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

import os,sys
from os.path import join, basename, dirname, splitext
if __name__=='__main__': sys.path.append(join(dirname(__file__),'..'))
import time
#import plottools
import krebsutils
import h5py
import numpy as np
import matplotlib.cm
import myutils
import math

import povray as pv
from krebs.povrayEasy import *

#see http://assorted-experience.blogspot.com/2007/07/custom-colormaps.html
cm_redblue = matplotlib.colors.LinearSegmentedColormap('redblue', {
    'red' : ((0,0,0), (1,1,1)),
    'green' : ((0,0,0) ,(1,0,0)),
    'blue' : ((0,1,1), (1,0,0))
  })

filter_circulated = lambda e,kw: (kw['flags']&krebsutils.CIRCULATED)

def filter_small_rad(rmin):
    return lambda e,kw: np.bitwise_and(kw['rad']>rmin,kw['flags']&krebsutils.CIRCULATED)


def filter_tumorvessels1(e,kw):
    flags = kw['flags']
    rad = kw['rad']
    return np.logical_and(
        flags&krebs.CIRCULATED,
        np.logical_or(
            flags&krebs.WITHIN_TUMOR,
            rad>=10.0
        )
    )



def make_pressure_color_arrays(vesselgraph):
    edges = vesselgraph.edgelist
    flags = vesselgraph.edges['flags']
    data = vesselgraph.nodes['pressure']
    num_nodes = len(vesselgraph.nodes['position'])
    nflags = krebsutils.edge_to_node_property(num_nodes, edges, flags, 'or')

    is_set = lambda flags_,flag: np.asarray(np.bitwise_and(flags_, flag), np.bool)

    gray = np.asarray((0.1,0.1,0.1))
    lgray = np.asarray((0.5, 0.5, 0.5))

    circulated = is_set(flags,krebsutils.CIRCULATED)
    ncirculated = is_set(nflags,krebsutils.CIRCULATED)
    capillary = is_set(flags, krebsutils.CAPILLARY)
    ncapillary = is_set(nflags, krebsutils.CAPILLARY)
    p0 = np.amin(data[ncirculated])
    p1 = np.amax(data[ncirculated])
    cm = matplotlib.cm.ScalarMappable(cmap=cm_redblue)
    cm.set_clim(p0, p1)

    edgedata = np.average((data[edges[:,0]], data[edges[:,1]]), axis=0)
    edgecolors = cm.to_rgba(edgedata)[:,:3]
    edgecolors[~circulated] = gray
    #edgecolors[capillary] = lgray

    nodecolors = cm.to_rgba(data)[:,:3]
    nodecolors[~ncirculated] = gray
    #nodecolors[ncapillary] = lgray

    vesselgraph.edges['colors'] = edgecolors
    vesselgraph.nodes['colors'] = nodecolors



colors = [(0, 0.2,0,0.8), (0.45, 1,0,0), (0.6, 1, 1,0), (1., 1, 1, 1) ]
cm_hematocrit = matplotlib.colors.LinearSegmentedColormap('', {
  'red'   : [(x, r, r) for (x,r,g,b) in colors],
  'green' : [(x, g, g) for (x,r,g,b) in colors],
  'blue'  : [(x, b, b) for (x,r,g,b) in colors],
}, N = 256, gamma = 1.)


def make_any_color_arrays(vesselgraph, data_name):
  edges = vesselgraph.edgelist
  num_nodes = len(vesselgraph.nodes['position'])
  flags = vesselgraph.edges['flags']
  nflags = krebsutils.edge_to_node_property(num_nodes, edges, flags, 'or')

  mask = myutils.bbitwise_and(flags,krebsutils.CIRCULATED)
  nmask = myutils.bbitwise_and(nflags,krebsutils.CIRCULATED)

  if data_name in vesselgraph.edges:
    edgedata = vesselgraph.edges[data_name]
    nodedata = krebsutils.edge_to_node_property(num_nodes, edges, edgedata, 'avg')
  else:
    nodedata = vesselgraph.nodes[data_name]
    edgedata = np.average((nodedata[edges[:,0]], nodedata[edges[:,1]]), axis=0)

  gray = np.asarray((0.1,0.1,0.1))
  edgecolors = np.repeat(gray.reshape(1,-1), len(edgedata), axis=0)
  nodecolors = np.repeat(gray.reshape(1,-1), len(nodedata), axis=0)
  #colors = lambda arr: cm.to_rgba(arr)[:,:3]
  colors = lambda arr: np.power(cm.to_rgba(arr)[:,:3], 2.4)

  if data_name == 'hematocrit':
    cm = matplotlib.cm.ScalarMappable(cmap = cm_hematocrit)
    cm.set_clim(0, 1)
    unmapped_range = (0.,1.)
    edgecolors[mask] = colors(edgedata[mask])
    nodecolors[nmask] = colors(nodedata[nmask])
  elif data_name == 'pressure':
    p0 = np.amin(nodedata)
    p1 = np.amax(nodedata)
    unmapped_range = (p0, p1)
    cm = matplotlib.cm.ScalarMappable(cmap=cm_redblue)
    cm.set_clim(p0, p1)
    edgecolors[mask] = colors(edgedata[mask])
    nodecolors[nmask] = colors(nodedata[nmask])
  elif data_name == 'shearforce':
    mask = mask & (edgedata>0)
    nmask = nmask & (nodedata>0)
    edgedata = edgedata[mask]
    nodedata = nodedata[nmask]
    unmapped_range = edgedata.min(), edgedata.max()
    edgedata = np.log10(edgedata)
    nodedata = np.log10(nodedata)
    p0 = -4#np.amin(edgedata)
    p1 = -1#np.amax(edgedata)
    cm = matplotlib.cm.ScalarMappable(cmap=matplotlib.cm.spectral)
    cm.set_clim(p0, p1)
    edgecolors[mask] = colors(edgedata)
    nodecolors[nmask] = colors(nodedata)
  elif data_name == 'flow':
    mask = mask & (edgedata>0)
    nmask = nmask & (nodedata>0)
    edgedata = edgedata[mask]
    nodedata = nodedata[nmask]
    unmapped_range = edgedata.min(), edgedata.max()
    edgedata = np.log10(edgedata)
    nodedata = np.log10(nodedata)
    p0 = np.floor(np.amin(edgedata))
    p1 = np.ceil(np.amax(edgedata))
    cm = matplotlib.cm.ScalarMappable(cmap=matplotlib.cm.jet)
    cm.set_clim(p0, p1)
    edgecolors[mask] = colors(edgedata)
    nodecolors[nmask] = colors(nodedata)
  elif data_name == 'flags':
    unmapped_range = (0.,1.)
    mask = mask & (edgedata>0)
    nmask = nmask & (nodedata>0)
    edgedata = edgedata[mask]
    nodedata = nodedata[nmask]
    edgedata = np.bitwise_and(edgedata,krebsutils.ARTERY)
    cm = matplotlib.cm.ScalarMappable(cmap=matplotlib.cm.jet)
    edgecolors[mask] = colors(edgedata)
    nodecolors[nmask] = colors(nodedata)
  vesselgraph.edges['colors'] = edgecolors
  vesselgraph.nodes['colors'] = nodecolors
  return cm, unmapped_range



def addVesselTree(epv, vesselgraph, trafo, **kwargs):
    edges = vesselgraph.edgelist
    pos = vesselgraph['position']
    rad = vesselgraph.edges['radius']
    #print 'positions = ', np.amin(pos, axis=0), np.amax(pos, axis=0)
    pos = transform_position(trafo, pos)
    rad = transform_scalar(trafo, rad)
    #print 'positions after trafo = ', np.amin(pos, axis=0), np.amax(pos, axis=0)

    clip_vessels = clipFactory(kwargs.pop('vessel_clip', None))

    edgecolors = vesselgraph.edges['colors']
    nodecolors = vesselgraph.nodes['colors']

    class Styler(object):
      @staticmethod
      def edge_style(i, a, b):
        c = edgecolors[i]
        return "MkStyle(<%0.2f,%0.2f,%0.2f>)" % tuple(c)
      @staticmethod
      def node_style(i):
        c = nodecolors[i]
        return "MkStyle(<%0.2f,%0.2f,%0.2f>)" % tuple(c)

    if kwargs.pop('colored_slice',False):
      ClipStyler = Styler
    else:
      class ClipStyler(object):
        styleString = "pigment { color rgb<%0.2f,%0.2f,%0.2f> }" % (0.9, 0.9, 0.9)
        @staticmethod
        def edge_style(i, a, b):
          return ClipStyler.styleString
        @staticmethod
        def node_style(i):
          return ClipStyler.styleString
#	  ambient 0.05 * c
    epv.pvfile.write("""
    #macro MkStyle(c)
    texture {
	pigment {
	  color rgb c
	}
	finish {
	  specular 0.3
      ambient 1.
	}
    }
    #end
    """)

    #print pos.shape, rad.shape, edges.shape, clip_vessels
    epv.addVesselTree(np.asarray(edges, dtype=np.int32),
                      np.asarray(pos, dtype=np.float32),
                      np.asarray(rad, dtype=np.float32),
                      Styler, clip_vessels, ClipStyler)
    #del Styler; del ClipStyler
    #print 'write time: ', time.time(



def renderScene(vesselgroup, imagefn, **kwargs):
    vess_ldgroup = vesselgroup['lattice']
    graph = krebsutils.read_vessels_from_hdf(vesselgroup, ['position', 'flags', 'radius', 'pressure'], return_graph=True)
    colorfactory = kwargs.pop('colorfactory', make_pressure_color_arrays)
    colorfactory(graph)

    wbbox = krebsutils.read_lattice_data_from_hdf(vess_ldgroup).worldBox
    trafo = calc_centering_normalization_trafo(wbbox)
    zsize = (wbbox[5]-wbbox[4])

    with EasyPovRayRender(**kwargs) as epv:
      epv.setBackground(kwargs.pop('background',0.0))
      cam = kwargs.pop('cam','topdown') # this stuff is copy pasted everywhere, should be refactored in a extra function
      if cam in ('topdown', 'topdown_slice'):
        cam_fov = 60.
        cam_distance_factor = ComputeCameraDistanceFactor(cam_fov, kwargs['res'], wbbox)
        epv.addLight(10*Vec3(1.7,1.2,2), 1., area=(4, 4, 3, 3), jitter=True)
        if cam == 'topdown_slice':
          kwargs.update(vessel_clip =('zslice', -201*trafo.w, 201*trafo.w),
                        tumor_clip =('zslice', -100*trafo.w, 100*trafo.w))
          epv.setCamera((0,0,cam_distance_factor*0.5*(200.*trafo.w+2.)), (0,0,0), cam_fov, up = 'y')
        else:
          epv.setCamera((0,0,cam_distance_factor*0.5*(zsize*trafo.w+2.)), (0,0,0), cam_fov, up = 'y')
      else:
        basepos = np.asarray((0.6,0.7,0.7))*(1./1.4)
        epv.setCamera(basepos, (0,0,0), 90, up = (0,0,1))
        num_samples_large_light = 10
        num_samples_small_light = 3
        epv.addLight(10*Vec3(0.7,1.,0.9), 0.8, area=(1., 1., num_samples_small_light, num_samples_small_light), jitter=True)
        epv.addLight(10*Vec3(0.5,0.5,0.5), 0.6, area=(5., 5., num_samples_large_light, num_samples_large_light), jitter=True)
        kwargs.update(vessel_clip = ('pie', 0.),
                      tumor_clip = ('pie', -20*trafo.w))


      addVesselTree(epv, graph, trafo = trafo, **kwargs)
      epv.render(imagefn)


def CreateScene2(epv, wbbox, graph, imagefn, **kwargs):
  #wbbox = ld.GetWorldBox()
  trafo = calc_centering_normalization_trafo(wbbox)
  zsize = (wbbox[5]-wbbox[4])

  epv.setBackground(kwargs.pop('background',1.0))
  cam = kwargs.pop('cam','topdown')
  if cam in ('topdown', 'topdown_slice'):
    cam_fov = 60.
    cam_distance_factor = ComputeCameraDistanceFactor(cam_fov, kwargs['res'], wbbox)
    #epv.addLight(10*Vec3(1.7,1.2,2), 1., area=(4, 4, 3, 3), jitter=True)
    epv.addLight(10.*Vec3(1,0.5,2), 1.2)
    if cam == 'topdown_slice':
      kwargs.update(vessel_clip =('zslice', -301*trafo.w, 301*trafo.w),
                    tumor_clip =('zslice', -100*trafo.w, 100*trafo.w))
      epv.setCamera((0,0,cam_distance_factor*1.05), lookat = (0,0,0), fov = cam_fov, up = 'y')
    else:
      epv.setCamera((0,0,cam_distance_factor*0.5*(zsize*trafo.w+2.)), (0,0,0), cam_fov, up = 'y')
  else:
    basepos = np.asarray((0.6,0.7,0.7))*(1./1.4)
    epv.setCamera(basepos, (0,0,0), 90, up = (0,0,1))
    num_samples_large_light = 10
    num_samples_small_light = 3
    epv.addLight(10*Vec3(0.7,1.,0.9), 0.8, area=(1., 1., num_samples_small_light, num_samples_small_light), jitter=True)
    epv.addLight(10*Vec3(0.5,0.5,0.5), 0.6, area=(5., 5., num_samples_large_light, num_samples_large_light), jitter=True)
    kwargs.update(vessel_clip = ('pie', 0.),
                  tumor_clip = ('pie', -20*trafo.w))
  addVesselTree(epv, graph, trafo = trafo, **kwargs)




if (__name__ == '__main__'):
    import optparse
    parser = optparse.OptionParser()
    parser.add_option("-d","--data", dest="datalist", help="which data (pressure, flow, shearforce, hematocrit) as comma separated list", default='pressure', action="store")
    parser.add_option("-f","--filter-uncirculated", dest="filteruncirculated", help="filter uncirculated vessels", default=False, action="store_true")
    parser.add_option("--filter-radius-high-pass", dest="filterradiushighpass", action="store", type="float", default = -1)
    parser.add_option("--no-overlay", dest="overlay", default = True, action="store_false")
    options, args = parser.parse_args()

    filenames = args[:-1]
    pattern   = args[-1]
    datalist = map(lambda s: s, map(str.strip, options.datalist.split(',')))

    labels = {
      'flow' : 'Flow Rate /nl',
      'shearforce': '$log_{10}$ Shear Force',
      'hematocrit' : 'Hematocrit',
      'pressure' : 'Blood Pressure',
      'flags' : ' Flags',
    }

    for fn in filenames:
      f = h5py.File(fn,'r')
      dirs = myutils.walkh5(f['.'], pattern)
      for d in dirs:
        filenamepostfix = ''
        vesselgroup = f[d]
        if( vesselgroup.attrs['CLASS'] == 'GRAPH' ):
          worldbox = krebsutils.read_lattice_data_from_hdf(vesselgroup['lattice']).GetWorldBox()
        if( vesselgroup.attrs['CLASS'] == 'REALWORLD'):
          pos=vesselgroup['nodes/world_pos']
          x_min = np.min(pos[:,0])
          x_max = np.max(pos[:,0])
          y_min = np.min(pos[:,1])
          y_max = np.max(pos[:,1])
          z_min = np.min(pos[:,2])
          z_max = np.max(pos[:,2])
          worldbox=np.asarray([x_min, x_max, y_min, y_max, z_min, z_max])
          
        graph = krebsutils.read_vessels_from_hdf(vesselgroup, ['position', 'flags', 'radius'] + datalist, return_graph=True)
        if options.filteruncirculated:
          graph = graph.get_filtered(edge_indices = myutils.bbitwise_and(graph['flags'], krebsutils.CIRCULATED))
        if options.filterradiushighpass>0:
          graph = graph.get_filtered(edge_indices = graph['radius']> options.filterradiushighpass)
          filenamepostfix = '_rhp'
        for data_name in datalist:
          cm, (datamin, datamax) = make_any_color_arrays(graph, data_name)
          imagefn = splitext(basename(fn))[0]+myutils.sanitize_posixpath(vesselgroup.name).replace('/','-')+'_'+data_name+filenamepostfix+'.png'
          kwargs = dict(
            #res=(2048,2048), aa=4, dpi = 320,
            res=(1024,1024), aa=2, dpi = 180,
            #clip_pie=0,
            #clip_zslice=(0.4, 0.6),
            #cam='corner',
            #cam='topdown_slice',
            cam = 'topdown',
            colored_slice = False,
            out_alpha = True,
            num_threads = 2,
            background = 0.8,
            #ambient_color = (0.5, 0.5, 0.5),
            ambient_color = (0.3, 0.3, 0.3),
          )
          with EasyPovRayRender(**kwargs) as epv:
            CreateScene2(epv, worldbox, graph, imagefn, **kwargs)
            if options.overlay:
              RenderImageWithOverlay(epv, imagefn, cm, labels[data_name], **kwargs)
            else:
              epv.render(imagefn)
    #        z = krebsutils.read_lattice_data_from_hdf(vesselgroup['lattice']).GetWorldSize()[2]
    #        renderScene(vesselgroup, outname,
    #          res=(1200,1200), aa=2,
    #          #clip_pie=0,
    #          #clip_zslice=(0.4, 0.6),
    #          #cam='corner',
    #          cam='topdown',
    #          colored_slice=False,
    #          out_alpha=False,
    #          num_threads=2,
    #          debug=False
    #          )
