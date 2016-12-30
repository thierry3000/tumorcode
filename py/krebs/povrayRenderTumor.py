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
if __name__ == '__main__':
  import os.path, sys
  sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'..'))


import os,sys
from os.path import basename, splitext

from povrayRenderVessels import *
import myutils
import math
import copy


def addBulkTissueTumor(epv, tumorgroup, trafo, options):
    ld = krebsutils.read_lattice_data_from_hdf(tumorgroup.file[tumorgroup['conc'].attrs['LATTICE_PATH']])
    ld = transform_ld(trafo, ld)

    ds_necro    = np.asarray(tumorgroup['necro'])
    data        = np.clip(np.asarray(tumorgroup['conc']), 0., 1.) # - np.asanyarray(tumorgroup['necro']), 0., 1.)

    ds_levelset = -np.minimum(np.asarray(tumorgroup['ls']), 0.4 - ds_necro)
    ds_levelset = krebsutils.distancemap(ds_levelset) * ld.scale

#    import matplotlib
#    matplotlib.use('Qt4Agg')
#    import matplotlib.pyplot as pyplot
#    pyplot.imshow(ds_levelset[:,:,8])
#    pyplot.contour(ds_levelset[:,:,8],[0.])
#    pyplot.show()
    if 'tumor_clip' in options:
      clip = clipFactory(options.tumor_clip)
    else:
      clip = clipFactory('None')

    voldata_ls    = epv.declareVolumeData(ds_levelset, ld.GetWorldBox())
    voldata_cells = epv.declareVolumeData(data, ld.GetWorldBox())

    value_bounds = voldata_cells.value_bounds
    style = """
      texture {
        pigment {
          function { %f + %f*%s(x,y,z) }
          color_map {
            [0.0  color <0.3,0,0>]
            [0.5  color <1,0.8, 0.3>]
            [0.8  color <1,1,0.1>]
          }
        }
        finish { 
          specular 0.3
        }
      }""" % (value_bounds[0], (value_bounds[1]-value_bounds[0]), voldata_cells.name)
    #style = " texture { pigment { color rgb<1,0.8,0.3> }  finish { specular 0.3 }}"
    epv.addIsosurface(voldata_ls, 0., lambda : style, clip, style)



def renderScene(vesselgroup, tumorgroup, imagefn, options):
    if vesselgroup is not None:
      wbbox = krebsutils.read_lattice_data_from_hdf(vesselgroup['lattice']).worldBox
    else:
      wbbox = krebsutils.read_lattice_data_from_hdf(tumorgroup.file['field_ld']).worldBox
    trafo = calc_centering_normalization_trafo(wbbox)
    zsize = (wbbox[5]-wbbox[4])
    
    vessel_ld = krebsutils.read_lattice_data_from_hdf(vesselgroup['lattice'])  
    options.wbbox = vessel_ld.GetWorldBox()     

    with EasyPovRayRender(options) as epv:
      epv.setBackground(options.background)
      cam = options.cam
      if cam in ('topdown', 'topdown_slice'):
        cam_fov = 60.
        cam_distance_factor = options.cam_distance_multiplier * ComputeCameraDistanceFactor(cam_fov, options.res, wbbox)
        epv.addLight(10*Vec3(1.7,1.2,2), 1., area=(4, 4, 3, 3), jitter=True)
        if cam == 'topdown_slice':
          imagefn += '_slice'
          options.vessel_clip =('zslice', -201*trafo.w, 201*trafo.w)
          options.tumor_clip=('zslice', -100*trafo.w, 100*trafo.w)
          epv.setCamera((0,0,cam_distance_factor*0.5*(200.*trafo.w+2.)), (0,0,0), cam_fov, up = 'y')
        else:
          imagefn += '_top'
          epv.setCamera((0,0,cam_distance_factor*0.5*(zsize*trafo.w+2.)), (0,0,0), cam_fov, up = 'y')
      else:
        imagefn += '_pie'
        cam_fov = 60.
        basepos = options.cam_distance_multiplier * np.asarray((0.6,0.7,0.55))*(1./1.4)*math.tan(math.pi*0.25)/math.tan(math.pi/180./2.*cam_fov)
        epv.setCamera(basepos, (0,0,0), cam_fov, up = (0,0,1))
        num_samples_large_light = 4
        num_samples_small_light = 2
        epv.addLight(10*Vec3(0.7,1.,0.9), 0.8, area=(1., 1., num_samples_small_light, num_samples_small_light), jitter=True)
        epv.addLight(10*Vec3(0.5,0.5,0.5), 0.6, area=(5., 5., num_samples_large_light, num_samples_large_light), jitter=True)
        options.vessel_clip =('pie', 0.)
        options.tumor_clip = ('pie', -50*trafo.w)

      if vesselgroup is not None:
        graph = krebsutils.read_vessels_from_hdf(vesselgroup, ['position', 'flags', 'radius', 'pressure'], return_graph=True)
        if options.filteruncirculated:
          graph = graph.get_filtered(edge_indices = myutils.bbitwise_and(graph['flags'], krebsutils.CIRCULATED))        
        if 'colorfactory' in options:
          colorfactory = options.colorfactory
        else:
          colorfactory = make_pressure_color_arrays
        colorfactory(graph)
        #addVesselTree(epv, graph, trafo, vesselgroup = vesselgroup, options)
        addVesselTree(epv, graph, trafo, options)
        
      if tumorgroup is not None and 'conc' in tumorgroup:
        addBulkTissueTumor(epv, tumorgroup, trafo, options)
      
      if (tumorgroup is not None and tumorgroup.attrs['TYPE'] == 'faketumor'):
        print('nix')
      
      if options.noOverlay:
        epv.render(imagefn+'.png')
      else:
        RenderImageWithOverlay(epv, imagefn+'.png', None, 'tumor', options)
      



#if __name__ == '__main__':
#    print("Note: works only with vessels!")
#    print("Usage: filenames group")
#    filenames = sys.argv[1:-1]
#    grp_pattern = sys.argv[-1]
#    import povrayRenderSettings
#    settings = copy.deepcopy(povrayRenderSettings.image)
#    settings.update(povrayRenderSettings.tumor)
#    for fn in filenames:
#      with h5py.File(fn, 'r') as f:
#        dirs = myutils.walkh5(f['.'], grp_pattern)
#        outname = splitext(basename(fn))[0]
#        for d in dirs:
#          print 'rendering %s %s' % (f.filename, d)
#          renderScene(f[d]['vessels'], f[d]['tumor'], '%s-%s' % (outname,d), **settings)
