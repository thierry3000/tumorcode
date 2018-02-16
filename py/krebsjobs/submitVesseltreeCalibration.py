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

"""
  The purpose of this module is to experimentally (in-silico) determine a relationship
  between the lattice constant and the histological MVD and the line density S_D.
  In addition, it allows for conveniently running simulation with different
  parameters for a given lattice constant in order to determine configurations
  with desired rBV and rBF. Corresponding code for analysis can be found here, too.
"""
if __name__ == '__main__':
  import os.path, sys
  sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'..'))

import os, sys

import qsub
import copy
import myutils
import krebsutils
import krebs.analyzeBloodFlow
import h5py
import numpy as np

import krebsjobs.parameters.parameterSetsVesselGen


def ObtainDataOfVesselFile(f):
  dataman = myutils.DataManager(20, [krebs.analyzeGeneral.DataTumorTissueSingle(), 
                                      krebs.analyzeGeneral.DataVesselRadial(), 
                                      krebs.analyzeGeneral.DataDistanceFromCenter(),
                                      krebs.analyzeBloodFlow.DataTumorBloodFlow(),
                                      krebs.analyzeGeneral.DataBasicVessel(),
                                      krebs.analyzeGeneral.DataVesselSamples(),
                                      ])
  vesselgroup = f['vessels'] 
  ld = krebsutils.read_lattice_data_from_hdf(f['field_ld'])
  #print 'field_box = ', ld.worldBox
  #bins_spec   = krebs.analyzeGeneral.BinsSpecRange(100., 1000., 100.)
  bins_spec   = krebs.analyzeGeneral.BinsSpecRange(100., 1000., 100.)
  mvd, mvd_bins = dataman.obtain_data('sphere_vessel_density',  vesselgroup, None, bins_spec, 'radial', ld, (f, 'data'))
  rbf, rbf_bins = dataman.obtain_data('cum_rbf_radial', vesselgroup, None, bins_spec, 'radial', ld, (f, 'data'))
  rbv = dataman.obtain_data('basic_vessel_radial', 'phi_vessels', vesselgroup, None, 50., bins_spec, 'radial', ld, (f, 'data'))
  # rbv returns myutils.MeanValueArray
  scale = f['vessels/lattice'].attrs['SCALE']
  message = f['parameters'].attrs['MESSAGE']
  return dict(mvd = np.asarray(mvd), rbf = np.asarray(rbf), rbv = rbv.avg, bins = bins_spec.arange(), scale = scale, message = message)


def CenterTheLattice(f, h5_path):
    """
      The histological MVD is determined for intersections of
      the vascular network with a number of concentric spheres.
      Therefore it is convenient to first move the system
      so that it is centered at the coordinate origin.
    """
    ld = krebsutils.read_lattice_data_from_hdf(f[h5_path])   
    del f[h5_path]
    ld = ld.GetCentered()
    fn=str(f.file.filename)
    krebsutils.write_lattice_data_to_hdf_by_filename(fn, h5_path, ld)
    

def RunsOnClient(configstring_file, workdir, vesselfilename):
  os.chdir(workdir)
  if not os.path.exists(vesselfilename):
    print 'GENERATING %s' % vesselfilename
    with open(configstring_file, 'r') as f:
      configstring = f.read()
    krebsutils.run_vesselgen(configstring)
  with h5py.File(vesselfilename, 'r+') as f:
    # centering is needed because quantities are analyzed in dependence on the distance from the system origin!!
    CenterTheLattice(f, 'field_ld')
    CenterTheLattice(f, 'vessels/lattice')
    f.flush()
    ObtainDataOfVesselFile(f)


if not qsub.is_client:
  from krebsjobs.submitVesselgeneration import run_config_samples, configuration_type_factories, write_config
  from krebs.vesselgenerator import fix_worldsize, VD   
  
  
  def SubmitEnsemble(type, filename, options_):
     '''submits an ensemble of vessel networks for generation'''
     def Factory(i):
       '''prepares the network generation parameters for submission of sample i'''
       options = copy.deepcopy(options_)
       num_points = options.pop('num_points')
       scale      = options.pop('scale')
       hiter      = options.pop('hiter')
       vd = VD(shape=(num_points, num_points, num_points), scale=scale, latticetype='FCC', num_hierarchical_iterations=hiter, ensemble_index=i, **options)
       vd = fix_worldsize(vd)
       vd.name        = options.get('name', filename)  #this goes into the message attribute of the output file, default to filename
       vd.outfilename = filename
       vd.num_threads = 8
       vd = configuration_type_factories[type](vd, i)
       return vd
     #either submit or generate a config file (for debugging)
     if 0: # config file
       write_config(Factory)
     else:
       index_range = range(3) # two samples per config
       run_config_samples(Factory, index_range, RunsOnClient)      

  def SubmitOtherVariationJobsSwine():
    cube_width=2000
    hiter = 2
    scale = 50
    nums_points = int(cube_width / (2**hiter * scale) + 1)
    baseparams = dict(
      tip_radius_arterial = 2.0,
      tip_radius_capi = 2.0,
      tip_radius_vein = 3.72, #MW scale relation
      murray_alpha_vein = 3.,
      murray_alpha_artery = 3.,
      max_sprout_radius_artery = 8.,
      max_sprout_radius_vein = 8.,
      calcflow = dict(
        viscosityPlasma = 1.2e-6,
        rheology = 'RheologySecomb2005', # WARNING: changed to new model of secomb
        inletHematocrit = 0.45,
        includePhaseSeparationEffect = False,
      ),
      scale = scale, #to obtain 300 MVD from scale jobs
      name = 'base',
      hiter = hiter,
      num_points = nums_points,
    )
    variant1 = copy.deepcopy(baseparams)
    variant1['max_sprout_radius_artery'] = 4.
    variant1['max_sprout_radius_vein'] = 4.
    variant1['name'] = 'lower_sprout_radius_limit'
    
    variant2 = copy.deepcopy(baseparams)
    variant2['murray_alpha_vein'] = 2.7
    variant2['murray_alpha_artery'] = 2.7
    variant2['name'] = 'lower_murray_alpha'
    
    variant3 = copy.deepcopy(baseparams)
    variant3['generate_more_capillaries'] = True
    variant3['name'] = 'more_capillaries'
    
    variant4 = copy.deepcopy(baseparams)
    variant4['murray_alpha_vein'] = 3.5
    variant4['murray_alpha_artery'] = 3.0
    variant4['name'] = 'higher_murray_alpha'
    
    for t in 'typeI typeA typeB typeC typeD typeE typeF typeG typeH'.split():
      for params in [baseparams, variant1, variant2, variant3, variant4]:
        filename = 'vess_calibr_%s_%s' % (params['name'], t)
        SubmitEnsemble(t, filename, params)
  
  def SubmitOtherVariationJobs():
    baseparams = dict(
      tip_radius_arterial = 2.5,
      tip_radius_capi = 2.5,
      tip_radius_vein = 3.4,
      murray_alpha_vein = 3.,
      murray_alpha_artery = 3.,
      max_sprout_radius_artery = 8.,
      max_sprout_radius_vein = 8.,
      calcflow = dict(
        viscosityPlasma = 1.2e-6,
        rheology = 'RheologySecomb2005', # WARNING: changed to new model of secomb
        inletHematocrit = 0.45,
        includePhaseSeparationEffect = False,
      ),
      scale = 90.,
      name = 'base',
      hiter = 1,
      num_points = 13,
    )
    variant1 = copy.deepcopy(baseparams)
    variant1['max_sprout_radius_artery'] = 4.
    variant1['max_sprout_radius_vein'] = 4.
    variant1['name'] = 'lower_sprout_radius_limit'
    
    variant2 = copy.deepcopy(baseparams)
    variant2['murray_alpha_vein'] = 2.7
    variant2['murray_alpha_artery'] = 2.7
    variant2['name'] = 'lower_murray_alpha'
    
    variant3 = copy.deepcopy(baseparams)
    variant3['generate_more_capillaries'] = True
    variant3['name'] = 'more_capillaries'
    
    variant4 = copy.deepcopy(variant1)
    variant4['hiter'] = 2
    variant4['name'] = 'lower_sprlimit_large'
    
    for t in 'typeI typeA typeB typeC typeD typeE typeF typeG typeH'.split():
      for params in [variant4]: #[baseparams, variant1, variant2]:
        filename = 'vess_calibr_%s_%s' % (params['name'], t)
        SubmitEnsemble(t, filename, params)
    
    
  
  def SubmitScaleVariationJobs():
    params = dict(
      tip_radius_arterial = 2.5,
      tip_radius_capi = 2.5,
      tip_radius_vein = 3.8,
      # estimating from relative volume of blood in veins 65%, 
      # one obtains a radius scale factor of r(vein) = 1.36 * r(arterial)
      murray_alpha_vein = 3.,
      murray_alpha_artery = 3.,
      max_sprout_radius_artery = 8.,
      max_sprout_radius_vein = 8.,
      calcflow = dict(
        viscosityPlasma = 1.2e-6,
        #rheology = 'RheologySecomb2005', # WARNING: changed to new model of secomb
        rheology = 'RheologyForRats', # WARNING: changed to new model of secomb
        inletHematocrit = 0.45,
        includePhaseSeparationEffect = False,
      ),
    )
    # @scale = 120, we have an histological MVD of ca. 60-70  = MVD1
    # MVD scales with scale (h) like MVD1 / MVD2 = h2^2 / h1^2
    # thus to obtain an MVD2 that is ca. 3 * MVD1, we scale h by sqrt(1/3) = 0.57
    if 0:
      hiters      = [ 1,    1,   1,    1,    1,   1,   1 ]
      scales      = [60., 70., 80., 90., 100., 110., 120.]
      cube_width  = 2200.
    if 0:
      hiters      = [ 1,    1,   1,    1,]
      scales      = [60., 80., 100., 120.]
      cube_width  = 1100.
    '''try to match luedemann et al'''
    if 0:
      params = getattr(krebsjobs.parameters.parameterSetsVesselGen, 'paramset24')
      params['calcflow'].update(
        rheology = 'RheologySecomb2005' #irrelevant since includePhaseSeparationEffect = False
      )
      hiters      = [ 2,    2,   2,    2,   2,]
      scales      = [50., 75., 100., 125.,150,]
      cube_width  = 1100.
    if 1:
      params = getattr(krebsjobs.parameters.parameterSetsVesselGen, 'paramset24')
      params['calcflow'].update(
        rheology = 'RheologySecomb2005' #irrelevant since includePhaseSeparationEffect = False
      )
      hiters      = [ 2,    2,   2,    2,   2,]
      scales      = [40., 45., 50.,  55.,  60,]
      cube_width  = 2000.
    # estimated size:
    # 2^hiter * scale * (num_points - 1)  =  S
    # num_points = S / 2^hiter / scale + 1
    # num_points >! 5
    # -> (10-1) * 120 * 1**1 = 1920 um!!
    ''' try to match
    1.Kyle, A. H., Baker, J. H. E., Gandolfo, M.-J., Reinsberg, S. A. & Minchinton, A. I. Tissue Penetration and Activity of Camptothecins in Solid Tumor Xenografts. Mol Cancer Ther 13, 2727–2737 (2014).
    supplementary Fig. S2
    6mm^2 tumor, try 1cm^2 tissue
    MVD \aprox 70
    Result: seems like scale 110 is doing a decent job!
    '''
    if 0:
#      hiters      = [   2,  2,   2,    2,]
#      scales      = [140.,80.,100., 120.]
      hiters      = [   2,  2,   2,    2,]
      scales      = [110.,90.,100., 105.]
      cube_width  = 1000.
    nums_points = map(lambda (hiter, scale): int(cube_width / (2**hiter * scale) + 1), zip(hiters, scales))
    #### main loop that does the submitting of various configurations #####
    for t in 'typeI typeA typeB typeC typeD typeE typeF typeG typeH'.split():
    #for t in ['typeI']:
      for (hiter, scale, num_points) in zip(hiters, scales, nums_points):  
        filename = 'vess_calibr_s%03i_ra%03i_rc%03i_rv%03i_%s' % (scale, params['tip_radius_arterial'], params['tip_radius_capi'], params['tip_radius_vein'], t)
        params['hiter'] = hiter
        params['scale'] = scale
        params['num_points'] = num_points
        SubmitEnsemble(t, filename, params)


  def ComputeSingleNumberAvgStd(datasets):
      '''
        input dict: name -> list of radial curves

        List of radial curves is summarized into a single value.
        Care must be taken in the selection of a suitable range to
        avoid boundary effects. The thresholds x_min and x_max
        should be within the minimal extends of the system bounding obx.
        
        retuns dict of tuples
      '''
      x_min = 200.
      x_max = 500.
      bins = datasets['bins']
      x     = bins
      x_rbv = 0.5*(x[:,1:]+x[:,:-1])
      mask1 = np.bitwise_and(x < x_max, x > x_min)
      mask2 = np.bitwise_and(x_rbv < x_max, x_rbv > x_min)
      masks = {
        'rbv' : mask2,
        'mvd' : mask1,
        'rbf' : mask1,
        'scale' : None,
      }
      result = dict()
      for k, v in datasets.items():
        if k == 'bins': continue
        mask = masks[k]
        if mask is not None:
#            print 'k = ', k
#            print mask.shape, mask
#            print v.shape, v
          tmp = v[mask].ravel()
        else:
          tmp = v
        avg = np.average(tmp)
        std = np.std(tmp)/np.sqrt(len(tmp))
        result[k] = (avg, std)
      return result


  def AnalyzeScaleVariation(filenames):
    import matplotlib.pyplot as pyplot
    import mpl_utils
    import krebs.quantities as Q
    import collections

    data = []
    for fn in filenames:
      with h5py.File(fn, 'r+') as f:
        sample = ObtainDataOfVesselFile(f)
        del sample['message']
        data.append(sample)

    byScale = collections.defaultdict(list)
    for d in data:
      scale = d['scale']
      byScale[scale].append(d)
    for k, v in byScale.items():
      byScale[k] = myutils.zipListOfDicts(v) 
          
    curves = collections.defaultdict(list)
    for k, v in byScale.items():
      res = ComputeSingleNumberAvgStd(v)
      for name, (std, avg) in res.items():
        curves[name].append((std, avg))
    order = np.argsort(np.asarray(curves['scale'])[:,0])
    for k, v in curves.items():
      curves[k] = np.asarray(v).transpose()[:, order]

    scales = {
      'mvd': (Q.um**-2).asNumber(Q.mm**-2),
      'rbv': 100.,
      'rbf': 60.,
    }

    with mpl_utils.PageWriter('vessel-calibration-analysis', fileformats=['pdf']) as pdfwriter:
      fig, axes = pyplot.subplots(3,1, figsize = mpl_utils.a4size*np.asarray([0.4, 0.5]))
      for scale, data in byScale.items():
        bins = np.average(data['bins'], axis=0) # sanity check, all bins arrays have equal size, so just try to average the bin boundaries, even if it makes no real sense
        x    = bins
        x_rbv    = 0.5*(x[1:]+x[:-1]) # bin center for rBV        
        # plot things
        ax = axes[0]
        ya = data['mvd']
        ax.errorbar(x, scales['mvd']*np.average(ya, axis=0), yerr = scales['mvd']*np.std(ya, axis=0), label = ('h = %0.f' % scale))
        legend = ax.legend(loc=4, fontsize='xx-small')        
        ax = axes[1]
        scale = scales['rbv']
        ya = data['rbv']
        ax.errorbar(x_rbv, scale*np.average(ya, axis=0), yerr = scale*np.std(ya, axis=0))        
        ax = axes[2]
        ya = data['rbf']
        scale = scales['rbf']
        ax.errorbar(x, scale*np.average(ya, axis=0), yerr = scale*np.std(ya, axis=0))
        axes[0].set(ylabel = 'mvd [$mm^{-1}$]')
        axes[1].set(ylabel = 'rbv [$\%$]')
        axes[2].set(ylabel = 'rbf [$min^{-1}$]', xlabel = '$|x| [\mu m]$')
      pyplot.tight_layout()
      pyplot.legend()
      pdfwriter.savefig(fig)
      
      fig, axes = pyplot.subplots(3,1, figsize = mpl_utils.a4size*np.asarray([0.4, 0.5]))

      for ax in axes:
        ax.grid(linestyle=':', linewidth=0.5, color='#aaaaaa')
      
      x = curves['scale'][0,:]
      ax = axes[0]
      y, yerr = scales['mvd']*curves['mvd']
      ax.errorbar(x, y, yerr = yerr)
      ax = axes[1]
      y, yerr = scales['rbv']*curves['rbv']
      ax.errorbar(x, y, yerr = yerr)
      ax = axes[2]
      y, yerr = scales['rbf']*curves['rbf']
      ax.errorbar(x, y, yerr = yerr)      
      
      axes[0].set(ylabel = 'mvd [$mm^{-1}$]')
      axes[1].set(ylabel = 'rbv [$\%$]')
      axes[2].set(ylabel = 'rbf [$min^{-1}$]', xlabel = '$h [\mu m]$')      
      
      pyplot.tight_layout()
      pdfwriter.savefig(fig)

  def AnalyzeOtherVariations(filenames):
    import krebs.quantities as Q
    import collections
    scales = {
      'mvd': (Q.um**-2).asNumber(Q.mm**-2),
      'rbv': 100.,
      'rbf': 60.,
    }

    dataByMsg = collections.defaultdict(list)
    for fn in filenames:
      with h5py.File(fn, 'r+') as f:
        sample = ObtainDataOfVesselFile(f)
        msg = sample.pop('message')
        dataByMsg[msg].append(sample)

    for msg, data in dataByMsg.items():
      data = myutils.zipListOfDicts(data)
      data = ComputeSingleNumberAvgStd(data)
      fmt = lambda k: (scales[k]*data[k][0], scales[k]*data[k][1])
      print('MSG = %s' % msg)
      print('  MVD = %f +/- %f' % fmt('mvd'))
      print('  rBV = %f +/- %f' % fmt('rbv'))
      print('  rBF = %f +/- %f' % fmt('rbf'))
      # rBF should be about 0.06 ml/ml/min !!!


if (not qsub.is_client) and __name__ == '__main__':
  import argparse   # for example see sync.py
  parser = argparse.ArgumentParser(description = 'Run a number of vessel tree generation jobs and analyze the resulting MVD, rBV and rBF')
  parser.add_argument('-a', '--analyze', help = 'analyze data, supply filenames', action='store_true')
  parser.add_argument('-n', '--no_scale_variation', help = 'No scale variation', default=False, action='store_true',)  
  parser.add_argument('FILES',  nargs='*')
  goodArguments, otherArguments = parser.parse_known_args()
  qsub.parse_args(otherArguments)
  
  if goodArguments.no_scale_variation:
    #submitFunction = SubmitOtherVariationJobs
    submitFunction = SubmitOtherVariationJobsSwine
    analyzeFunction = AnalyzeOtherVariations
  else:
    submitFunction = SubmitScaleVariationJobs
    analyzeFunction = AnalyzeScaleVariation

  if goodArguments.analyze:
    filenames = goodArguments.FILES
    analyzeFunction(filenames)
  else:
    submitFunction()
