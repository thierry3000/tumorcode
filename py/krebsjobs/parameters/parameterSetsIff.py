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

import numpy as np
from copy import deepcopy
import myutils
cluster_threads = myutils.cluster_threads
trastuzumab_factor = np.sqrt(540/145531.5)#see wikipedia
topotecan_factor = np.sqrt(540./421.)
irinotecan_factor = np.sqrt(540./587.)

#from drug_variant42-inj   
iff_default = dict(
  num_threads = 6,
  tumor = dict(
    #out_intervall = 100,  #fake tumor only??
    #tumor_speed = 2.,
    #tumor_radius = 250.,
    ),
  vessels = dict(
    probCollapse = 0.5,
    forceCollapse = 0.1 * 1.e-3, #3.e-3 *  4.3/8.8 * 0.5, # scaled from original ifp paper config to match lower shear stress in P10 vessel configs
    distSproutMin = 3,
    radMax = 25.,
    #radius under max pressure is = reference_r * vesselCompressionFactor
    vesselCompressionFactor = 1.,
    gfVessProl = 0.0005,
    timeProlEcSprout  = 2,  # 2h per added segment
    timeEcRadiusInflate = 300, # after 360 hours, vessel radius should be ca 27micron
    timeProlEcSproutLifetime = 50,  # the time until sprouts regress
    sproutDelay = 24, # the time during which a tumor vessel can still generate sprouts (hours)
    bRelativeShearforceCollapse = False,
    radInit   = 3.5,
    #images in (holash et al.) indicate approximately 0.04,
    #based on the time it takes to destroy a 25micron *radius* vessel
    #maturation of 25 micron vessel ca 16, 5 micron vessel ca 4.5
    dematuration_rate = 0.05,
    #initial wall thickness should be around 4 micron, so this is maturation_crit, so capillaries should immediately collapse
    #note that also vessels within normal tissue can collapse if their maturation value is < maturation_crit
    maturation_crit = 0.,
    #bArterialSprouting = True, # old code
    sproutMaxVesselWallThicknessArterial = 50., # replacement for arterial sprouting = true
    sproutMaxVesselWallThickness = 50.,
    bShearStressControlledDilatationAndRadiusDependentCollapse = False,
    
    bRadiusSmoothing = False,
    radiusSmoothingDiffusionCoeff = 1.,
  ),

  calcflow = dict(
    viscosityPlasma = 1.2e-6,
    rheology = 'RheologyForHuman',
    inletHematocrit = 0.37,
    includePhaseSeparationEffect = False,
  ),
  iff = dict(
#      iff_cond_tissue = 4 / .2,  ### FUUUUU jain 2007 used 187, 751 to 7510 for normal skin!
#      iff_cond_tumor  = 4 / .2,
#      iff_cond_necro  = 40 / .2,
      iff_cond_tissue = 20.,
      iff_cond_tumor  = 20.,
      iff_cond_necro  = 200.,
      ecm_friction    = 0.05,
      cell_friction   = 0.005,
      iff_capillary_permeability_tissue = 0.01,
      iff_capillary_permeability_tumor = 1.,
      lymph_surf_per_vol = 3. * (2. * 5. * 3.1415)*(1. / 100.**2), #0.0094
      lymph_press = -0.5,
      lymph_permeability = 0.02,
      capillary_oncotic_pressure = 2.7,
      interstitial_oncotic_pressure = 1.33,
      osmotic_reflection_coeff = 0.9,
      osmotic_reflection_coeff_tumor = 0.,
    ),
  ift = dict(
    inject_mode = 'DF_INJECT_MODE_EXP',
    uptake_mode = 'DF_UPTAKE_LINEAR',#this is default
    kdiff = 16.,
    inject_t = 1. * 3600., #10 * 3600,
    inject_max = 1.,
    stepper = 'vsimexbdf2',
    stepper_compartments = 'vsimexbdf2',
    capillary_permeability_normal = 0.00017,
    capillary_permeability_tumor = 0.017,
    comprates_k12 = 0.1,
    comprates_k21 = 0.001
  ),
  ift_measure = dict(),
  # out_times has to be a list for json to work
  # timepoints are in seconds!!!
  out_times = list(3600.*np.concatenate((np.arange(0., .5, 0.1),))),
  message = "iff_default",
  out_intervall = 10. * 60.,
  fn_out = "set me",
  fn_tumor = "set me",
  #tend = -1,
  h5_path_vessel = "set me",
  h5_path_lattice = "set me",
  h5_path_tumor = "set me",
#  tend = 300
)
ifp_paper_config = deepcopy(iff_default)
iff_defaultconfig2 = deepcopy(ifp_paper_config)
myutils.UpdateHierarchical(iff_defaultconfig2, dict(
  ift = dict(
    inject_mode = 'DF_INJECT_MODE_EXP',
    uptake_mode = 'DF_UPTAKE_LINEAR',#this is default
    kdiff = 16.,
    inject_t = 1. * 3600., #10 * 3600,
    inject_max = 1.,
    stepper = 'vsimexbdf2',
    stepper_compartments = 'vsimexbdf2',
    capillary_permeability_normal = 0.00017,
    capillary_permeability_tumor = 0.017,
    comprates_k12 = 0.1,
    comprates_k21 = 0.001
  ),
  ift_measure = dict(),
  # out_times has to be a list for json to work
  out_times = list(3600.*np.concatenate((np.arange(0., 1., 0.25), np.arange(1., 6., 1.), np.arange(6., 36., 3.), np.arange(36., 97., 6.)))),
  message = "iff_defaultconfig2",
  num_threads = 16,
  )
)
iff_defaultconfig2_JUMP = deepcopy(iff_defaultconfig2)
myutils.UpdateHierarchical(iff_defaultconfig2_JUMP, dict(
  ift = dict(
    inject_mode = 'DF_INJECT_MODE_JUMP',
    inject_t = 24. * 3600., #10 * 3600,
    )
  )
)
iff_MW_guess = deepcopy(iff_defaultconfig2)
myutils.UpdateHierarchical(iff_MW_guess, dict(
  ift = dict(
    inject_mode = 'DF_INJECT_MODE_EXP',
    capillary_permeability_normal=1.7e-6,
    capillary_permeability_tumor=1.7e-4,
    comprates_k12 = 0.0,
    comprates_k21 = 0.0,
    inject_t = 1. * 3600., #10 * 3600,
    kdiff = 0.16,
    )
  )
)
iff_MW_guess_linear = deepcopy(iff_MW_guess)
myutils.UpdateHierarchical(iff_MW_guess_linear, dict(
  # out_times has to be a list for json to work
  #out_times = list(3600.*np.concatenate((np.arange(0., 1., 0.25), np.arange(1., 6., 1.), np.arange(6., 36., 3.), np.arange(36., 97., 6.)))),
  ift = dict(
    inject_mode = 'DF_INJECT_MODE_EXP',
    capillary_permeability_normal=1.7e-6,
    capillary_permeability_tumor=1.7e-4,
    comprates_k12 = 0.0,
    comprates_k21 = 0.0,
    inject_t = 1. * 3600., #10 * 3600,
    kdiff = 0.16,
    ),  
  out_times = list(3600.*(np.arange(0., 96., 0.5))),
  message = "iff_MW_guess_linear",
  num_threads = 16,
  )
)
iff_MW_guess_linear_JUMP = deepcopy(iff_MW_guess)
myutils.UpdateHierarchical(iff_MW_guess_linear_JUMP, dict(
  # out_times has to be a list for json to work
  #out_times = list(3600.*np.concatenate((np.arange(0., 1., 0.25), np.arange(1., 6., 1.), np.arange(6., 36., 3.), np.arange(36., 97., 6.)))),
  ift = dict(
    inject_mode = 'DF_INJECT_MODE_JUMP',
    inject_t = 28. * 24.* 3600., #10 * 3600,
    ),  
  out_times = list(3600.*(np.arange(0., 96., 0.5))),
  message = "iff_MW_guess_linear_JUMP",
  num_threads = 16,
  )
)
iff_MW_guess_trastuzumab = deepcopy(iff_MW_guess_linear)
myutils.UpdateHierarchical(iff_MW_guess_trastuzumab, dict(
  # out_times has to be a list for json to work
  #out_times = list(3600.*np.concatenate((np.arange(0., 1., 0.25), np.arange(1., 6., 1.), np.arange(6., 36., 3.), np.arange(36., 97., 6.)))),
  ift = dict(
    inject_mode = 'DF_INJECT_MODE_EXP',
    capillary_permeability_normal=1.7e-6 * trastuzumab_factor,
    capillary_permeability_tumor=1.7e-4 * trastuzumab_factor,
    comprates_k12 = 0.0,
    comprates_k21 = 0.0,
    inject_t = 28 * 24. * 3600., #10 * 3600,
    kdiff = 0.16 * trastuzumab_factor,
    ),  
  out_times = list(3600.*(np.arange(0., 96., 0.5))),
  message = "iff_MW_guess_trastuzumab",
  num_threads = 16,
  )
)

iff_topotecan_1 = deepcopy(iff_default)
myutils.UpdateHierarchical(iff_topotecan_1, dict(
  # out_times has to be a list for json to work
  #out_times = list(3600.*np.concatenate((np.arange(0., 1., 0.25), np.arange(1., 6., 1.), np.arange(6., 36., 3.), np.arange(36., 97., 6.)))),
  ift = dict(
    inject_mode = 'DF_INJECT_MODE_EXP',
    capillary_permeability_normal=1.7e-6 * topotecan_factor,
    capillary_permeability_tumor=1.7e-4 * topotecan_factor,
    comprates_k12 = 0.0,
    comprates_k21 = 0.0,
    inject_t = 1200, #about 0.3 hrs
    kdiff = 0.16 * topotecan_factor,
    ),  
  out_times = list(3600.*(np.arange(0., 24., 0.5))),
  message = "iff_topotecan_1",
  num_threads = 16,
  )
)
iff_topotecan_2 = deepcopy(iff_topotecan_1)
iff_topotecan_2.update(message = "iff_topotecan_2",)
iff_topotecan_2['ift'].update(
  comprates_k12 = 0.2,
  comprates_k21 = 0.002,
)
iff_topotecan_3 = deepcopy(iff_topotecan_1)
iff_topotecan_3.update(message = "iff_topotecan_3",)
iff_topotecan_3['ift'].update(
  comprates_k12 = 0.05,
  comprates_k21 = 0.0,
)
iff_topotecan_4 = deepcopy(iff_topotecan_1)
iff_topotecan_4.update(
  message = "iff_topotecan_4",
  out_times = list(3600.*np.concatenate((np.arange(0.,2.0,0.2),np.arange(2., 24., 0.5)))),
)
iff_topotecan_4['ift'].update(
  comprates_k12 = 0.1,
  comprates_k21 = 0.0,
  kdiff = 0.16 * topotecan_factor *0.5,
)
iff_topotecan_5 = deepcopy(iff_topotecan_1)
iff_topotecan_5.update(
  message = "iff_topotecan_5",
  out_times = list(3600.*np.concatenate((np.arange(0.,2.0,0.1),np.arange(2., 24., 0.5)))),
)
iff_topotecan_5['ift'].update(
  comprates_k12 = 0.01,
  comprates_k21 = 0.0,
  kdiff = 0.16 * topotecan_factor,
)
iff_topotecan_6 = deepcopy(iff_topotecan_5)
iff_topotecan_6.update(
  message = "iff_topotecan_6",
)
iff_topotecan_6['ift'].update(
  comprates_k12 = 0.01 * topotecan_factor,
  comprates_k21 = 0.005 * topotecan_factor,
  kdiff = 0.16 * topotecan_factor,
)
iff_topotecan_7 = deepcopy(iff_topotecan_5)
iff_topotecan_7.update(
  message = "iff_topotecan_7",
)
iff_topotecan_7['ift'].update(
  comprates_k12 = 0.01 * topotecan_factor,
  comprates_k21 = 0.001 * topotecan_factor,
  kdiff = 0.16 * topotecan_factor,
)
iff_topotecan_8 = deepcopy(iff_topotecan_7)
iff_topotecan_8.update(
  message = "iff_topotecan_8",
  out_times = list(3600.*np.concatenate((np.arange(0.,2.0,0.1),np.arange(2., 73., 0.5)))),
)
iff_irinotecan_1 = deepcopy(iff_default)
myutils.UpdateHierarchical(iff_irinotecan_1, dict(
  # out_times has to be a list for json to work
  #out_times = list(3600.*np.concatenate((np.arange(0., 1., 0.25), np.arange(1., 6., 1.), np.arange(6., 36., 3.), np.arange(36., 97., 6.)))),
  ift = dict(
    inject_mode = 'DF_INJECT_MODE_EXP',
    capillary_permeability_normal=1.7e-6 * irinotecan_factor,
    capillary_permeability_tumor=1.7e-4 * irinotecan_factor,
    comprates_k12 = 0.01 * irinotecan_factor,
    comprates_k21 = 0.001 * irinotecan_factor,
    inject_t = 3240, #about 0.9 hrs
    kdiff = 0.16 * irinotecan_factor,
    ),  
  out_times = list(3600.*np.concatenate((np.arange(0.,2.0,0.1),np.arange(2., 24., 0.5)))),
  message = "iff_irinotecan_1",
  num_threads = 16,
  )
)
iff_irinotecan_2 = deepcopy(iff_irinotecan_1)
iff_irinotecan_2.update(
  message = "iff_irinotecan_2",
  out_times = list(3600.*np.concatenate((np.arange(0.,2.0,0.1),np.arange(2., 73., 0.5)))),
)
iff_small = deepcopy(iff_defaultconfig2)
myutils.UpdateHierarchical(iff_small, dict(
  ift = dict(
    kdiff = 0.16,
  ),
  # out_times has to be a list for json to work
  out_times = list(3600.*np.concatenate((np.arange(0., 1., 0.25), np.arange(1., 6., 1.), np.arange(6., 36., 3.), np.arange(36., 97., 6.)))),
  message = "iff_small",
  num_threads = 12,
  )
)
iff_small_highdiff = deepcopy(iff_small)
myutils.UpdateHierarchical(iff_small_highdiff, dict(
  ift = dict(
    kdiff = 50.,
  ),
  # out_times has to be a list for json to work
  out_times = list(3600.*np.concatenate((np.arange(0., 1., 0.25), np.arange(1., 6., 1.), np.arange(6., 36., 3.), np.arange(36., 97., 6.)))),
  message = "iff_small",
  num_threads = 12,
  )
)
iff_highdiff_JUMP = deepcopy(iff_defaultconfig2)
myutils.UpdateHierarchical(iff_highdiff_JUMP, dict(
  ift = dict(
    inject_mode = 'DF_INJECT_MODE_JUMP',
    kdiff = 50.,
    inject_t = 24. * 3600., #10 * 3600,
    )
  )
)
#iff_trastuzumab = deepcopy(iff_defaultconfig2)
iff_trastuzumab = deepcopy(iff_defaultconfig2)
myutils.UpdateHierarchical(iff_trastuzumab, dict(
  ift = dict(
    inject_mode = 'DF_INJECT_MODE_EXP',
    kdiff = 16.*trastuzumab_factor,
    inject_t = 1. * 3600., #10 * 3600,
    inject_max = 1.,
    stepper = 'vsimexbdf2',
    stepper_compartments = 'vsimexbdf2',
    capillary_permeability_normal = 0.00017*trastuzumab_factor,
    capillary_permeability_tumor = 0.017*trastuzumab_factor,
    comprates_k12 = 0.1*trastuzumab_factor,
    comprates_k21 = 0.001*trastuzumab_factor,
  ),
  ift_measure = dict(),
  # out_times has to be a list for json to work
  out_times = list(3600.*np.concatenate((np.arange(0., 1., 0.25), np.arange(1., 6., 1.), np.arange(6., 36., 3.), np.arange(36., 97., 6.)))),
  message = "drug test 1",
  num_threads = 16,
  )
)
iff_trastuzumab_JUMP = deepcopy(iff_defaultconfig2)
myutils.UpdateHierarchical(iff_trastuzumab_JUMP, dict(
  ift = dict(
    inject_mode = 'DF_INJECT_MODE_JUMP',
    kdiff = 16.*trastuzumab_factor,
    inject_t = 2. * 3600., #10 * 3600,
    inject_max = 1.,
    stepper = 'vsimexbdf2',
    stepper_compartments = 'vsimexbdf2',
    capillary_permeability_normal = 0.00017*trastuzumab_factor,
    capillary_permeability_tumor = 0.017*trastuzumab_factor,
    comprates_k12 = 0.1*trastuzumab_factor,
    comprates_k21 = 0.001*trastuzumab_factor,
  ),
  ift_measure = dict(),
  # out_times has to be a list for json to work
  out_times = list(3600.*np.concatenate((np.arange(0., 1., 0.25), np.arange(1., 6., 1.), np.arange(6., 36., 3.)))),
  message = "trastuzumab_JUMP",
  num_threads = 12,
  )
)
iff_trastuzumab2 = deepcopy(iff_trastuzumab)
myutils.UpdateHierarchical(iff_trastuzumab2, dict(
  ift = dict(
    inject_t = 1. * 3600., #10 * 3600,
  ),
  num_threads = 12,
  message = "iff_trastuzumab2",
  )
)
'''
Clinical Pharmacology of Trastuzumab
Dominique Levêque 1,* , Luc Gigou 1 and Jean Pierre Bergerat 2
see 28 days
'''
iff_trastuzumab3 = deepcopy(iff_trastuzumab)
myutils.UpdateHierarchical(iff_trastuzumab3, dict(
  ift = dict(
    inject_t = 28. * 24. * 3600., #10 * 3600,
  ),
  num_threads = 12,
  message = "iff_trastuzumab3",
  )
)
iff_trastuzumab3_JUMP = deepcopy(iff_trastuzumab3)
myutils.UpdateHierarchical(iff_trastuzumab3_JUMP, dict(
  ift = dict(
    inject_mode = 'DF_INJECT_MODE_JUMP',
  ),
  num_threads = 12,
  message = "iff_trastuzumab3_JUMP",
  )
)

iff_trastuzumab2_JUMP = deepcopy(iff_trastuzumab2)
myutils.UpdateHierarchical(iff_trastuzumab2_JUMP, dict(
  ift = dict(
    inject_mode = 'DF_INJECT_MODE_JUMP',
    inject_t = 1. * 3600., #10 * 3600,
  ),
  ift_measure = dict(),
  # out_times has to be a list for json to work
  out_times = list(3600.*np.concatenate((np.arange(0., 1., 0.25), np.arange(1., 6., 1.), np.arange(6., 36., 3.)))),
  message = "trastuzumab2_JUMP",
  num_threads = 12,
  )
)
iff_trastuzumab4 = deepcopy(iff_defaultconfig2)
myutils.UpdateHierarchical(iff_trastuzumab4, dict(
  ift = dict(
    inject_mode = 'DF_INJECT_MODE_EXP',
    kdiff = 16.*trastuzumab_factor,
    inject_t = 28. * 24. * 3600., #10 * 3600,
    inject_max = 1.,
    stepper = 'vsimexbdf2',
    stepper_compartments = 'vsimexbdf2',
    capillary_permeability_normal = 0.00017, #*trastuzumab_factor,
    capillary_permeability_tumor = 0.017, #*trastuzumab_factor,
    comprates_k12 = 0.1, #*trastuzumab_factor,
    comprates_k21 = 0.001, #*trastuzumab_factor,
  ),
  ift_measure = dict(),
  # out_times has to be a list for json to work
  out_times = list(3600.*np.concatenate((np.arange(0., 1., 0.25), np.arange(1., 6., 1.), np.arange(6., 36., 3.), np.arange(36., 97., 6.)))),
  message = "iff_trastuzumab4",
  num_threads = 12,
  )
)  
iff_trastuzumab5 = deepcopy(iff_trastuzumab4)
myutils.UpdateHierarchical(iff_trastuzumab5, dict(
  ift = dict(
    capillary_permeability_normal = 0.00017*trastuzumab_factor,
    capillary_permeability_tumor = 0.017*trastuzumab_factor,
    comprates_k12 = 0.1, #*trastuzumab_factor,
    comprates_k21 = 0.001, #*trastuzumab_factor,
  ),
  message = "iff_trastuzumab5",
  num_threads = 16,
  )
)
iff_trastuzumab6 = deepcopy(iff_trastuzumab4)
myutils.UpdateHierarchical(iff_trastuzumab6, dict(
  ift = dict(
    capillary_permeability_normal = 0.00017,#*trastuzumab_factor,
    capillary_permeability_tumor = 0.017,#*trastuzumab_factor,
    comprates_k12 = 0.1*trastuzumab_factor,
    comprates_k21 = 0.001*trastuzumab_factor,
  ),
  message = "iff_trastuzumab6",
  num_threads = 16,
  )
)
iff_trastuzumab7 = deepcopy(iff_trastuzumab4)
myutils.UpdateHierarchical(iff_trastuzumab7, dict(
  ift = dict(
    capillary_permeability_normal = 0.00017*trastuzumab_factor,
    capillary_permeability_tumor = 0.017*trastuzumab_factor,
    comprates_k12 = 0.1*trastuzumab_factor,
    comprates_k21 = 0.001*trastuzumab_factor,
  ),
  message = "iff_trastuzumab7",
  num_threads = 16,
  )
)
irinotecan_factor = np.sqrt(540/586.68)#see wikipedia
iff_irinotecan1 = deepcopy(iff_MW_guess_trastuzumab)
myutils.UpdateHierarchical(iff_irinotecan1, dict(
  # out_times has to be a list for json to work
  #out_times = list(3600.*np.concatenate((np.arange(0., 1., 0.25), np.arange(1., 6., 1.), np.arange(6., 36., 3.), np.arange(36., 97., 6.)))),
  ift = dict(
    inject_mode = 'DF_INJECT_MODE_EXP',
    capillary_permeability_normal=1.7e-6 * irinotecan_factor,
    capillary_permeability_tumor=1.7e-4 * irinotecan_factor,
    comprates_k12 = 0.0,
    comprates_k21 = 0.0,
    #in 'DF_INJECT_MODE_EXP' this means halftime
    # https://www.ncbi.nlm.nih.gov/pubmed/9342501, 5-27 hours--> mean 16
    # day hour seconds
    inject_t = 1 * 16. * 3600., 
    kdiff = 0.16 * irinotecan_factor,
    ),  
  out_times = list(3600.*(np.arange(0., 96., 0.5))),
  message = "iff_irinotecan1",
  num_threads = cluster_threads,
  )
)
iff_irinotecan1_quick = deepcopy(iff_irinotecan1)
myutils.UpdateHierarchical(iff_irinotecan1_quick, dict( 
  out_times = list(3600.*(np.arange(0., 3., 0.5))),
  message = "iff_irinotecan1_quick",
  )
)
if __name__=='__main__':
  print(myutils.namestr(irinotecan_factor,globals()))