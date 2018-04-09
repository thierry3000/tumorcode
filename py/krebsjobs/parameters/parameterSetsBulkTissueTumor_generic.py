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
''' Some details about the parameters
We use dict everywhere here.

vessels:
  if *bShearStressControlledDilatationAndRadiusDependentCollapse*
  than:
    vessels with radius smaller than *radMin* collapse with properbility *probCollapse*
  if NOT *bShearStressControlledDilatationAndRadiusDependentCollapse*
  than:
    parameters:forceCollapse, isShearForceEffectLinear and bRelativeShearforceCollapse 
    come to play. We interpolate linarly between 0 and *probCollapse*.
    We assign a StabilityAmountDueToShearStress in the following manner.
    
    if *isShearForceEffectLinear*
    than:
      if *bRelativeforceCollapse*
        current shearforce/ initial shearforce/ *forceCollapse*
      else
        current shearforce/ *forceCollapse*
    if NOT *isShearForceEffectLinear*
    than:
      if *bRelativeforceCollapse*
        ??
      else:
        if shearfoce< *forceCollapse --> we collapse the vessel with properbiltiy *probCollapse*
'''

circumferentialgrowth = dict(
  num_threads = 2,
  out_intervall = 100,
  tumor_speed = 2.,
  tumor_radius = 250.,
  vessels = dict(
    forceCollapse = 0.05e-3, # in kPa
    probCollapse = 0.1,
    distSproutMin = 8, # in number of lattice sites
    radMax = 20., # micrometer
    #radius under max pressure is = reference_r * vesselCompressionFactor
    vesselCompressionFactor = 1.,
    gfVessProl = 0.0005, # unitless concentration
    timeProlEcSprout  = 2,  # 2h per added segment
    timeEcRadiusInflate = 144, # after 360 hours, vessel radius should be ca 27micron
    timeProlEcSproutLifetime = 50,  # the time until sprouts regress
    sproutDelay = 12, # the time during which a tumor vessel can still generate sprouts (hours)
    bRelativeShearforceCollapse = False,
    radInit = 2.6,
    #images in (holash et al.) indicate approximately 0.04,
    #based on the time it takes to destroy a 25micron *radius* vessel
    #maturation of 25 micron vessel ca 16, 5 micron vessel ca 4.5
    dematuration_rate = 0.05, # micrometer per hour
    #initial wall thickness should be around 4 micron, so this is maturation_crit, so capillaries should immediately collapse
    #note that also vessels within normal tissue can collapse if their maturation value is < maturation_crit
    maturation_crit = 0., # micrometer
    #bArterialSprouting = True, # old code
    sproutMaxVesselWallThicknessArterial = 50., # replacement for arterial sprouting = true
    sproutMaxVesselWallThickness = 50.,
    bShearStressControlledDilatationAndRadiusDependentCollapse = False,
    isShearForceEffectLinear = False,
    bRadiusSmoothing = False,
    radiusSmoothingDiffusionCoeff = 1.,
    onlyReduceMaturationIfUnderperfused = False,
  ),
  calcflow = dict(
    viscosityPlasma = 1.2e-6,
    rheology = 'RheologyForHuman',
    inletHematocrit = 0.37,
    includePhaseSeparationEffect = False,
  ),
)


circumferentialgrowth2 = deepcopy(circumferentialgrowth)
circumferentialgrowth2['vessels'].update(
  radMax = 10.
)

circumferentialgrowth3 = deepcopy(circumferentialgrowth2)
circumferentialgrowth3['vessels'].update(
  probCollapse = 0.3, # 0.3 was kinda good with phi_v leveling off at ca 0.025
  forceCollapse = 0.075e-3,
)

circumferentialgrowth4 = deepcopy(circumferentialgrowth2)
circumferentialgrowth4['vessels'].update(
  probCollapse = 0.5, # 0.3 was kinda good with phi_v leveling off at ca 0.025
  forceCollapse = 0.1e-3,
)

circumferentialgrowth5 = deepcopy(circumferentialgrowth2)
circumferentialgrowth5['vessels'].update(
  probCollapse = 1.,
  forceCollapse = 0.1e-3,
  onlyReduceMaturationIfUnderperfused = True
)

circumferentialgrowth6 = deepcopy(circumferentialgrowth2)
circumferentialgrowth6['vessels'].update(
  probCollapse = 1.,
  forceCollapse = 0.5e-3,
  onlyReduceMaturationIfUnderperfused = True
)


circumferentialgrowth7 = deepcopy(circumferentialgrowth2)
circumferentialgrowth7['vessels'].update(
  probCollapse = 1.,
  forceCollapse = 0.25e-3,
  onlyReduceMaturationIfUnderperfused = True
)
trastuzumab_growth1 = deepcopy(circumferentialgrowth7)
trastuzumab_growth1.update(
  num_threads = 12,
  out_intervall = 30,
  #tumor_speed = 2., fake only
  #tumor_radius = 50., fake only
)
trastuzumab_growth1['vessels'].update(
  probCollapse = 0.3,
  forceCollapse = 0.25e-3,
  onlyReduceMaturationIfUnderperfused = True
)



flowctrldisintegrationconfig = dict(
  # This is a configuration where collapses are deterministic.
  # In this case the maturation variable decreases only when a vessel has
  # low shear stress. Collapse probability is set to 1.
  # Other changes are made as well. forceCollapse, sproutDelay and timeProlEcEnlarge.
  # In particular timeProlEcEnlarge was increased to match the experimental data again.
  num_threads = 2,
  out_intervall = 100,
  tumor_speed = 2.,
  tumor_radius = 200.,
  vessels = dict(
    gfVessProl = 0.0005,
    timeProlEcSprout = 2,
    timeProlEcEnlarge = 180,
    timeTillDeathEc = 50,
    sproutDelay = 24,
    probCollapse = 1,
    bRelativeShearforceCollapse = False,
    forceCollapse = 0.0005,
    forceCollapseRel = 0.2,
    radMax = 25,
    radInit = 4,
    distSproutMin = 8,
    dematuration_rate = 0.05,
    maturation_crit = 0.,
    bSproutModelSimple = False,
    vesselCompressionFactor = 1.,
    sproutMaxVesselWallThickness = 50,
    sproutMaxVesselWallThicknessArterial = 50., # replacement for arterial sprouting = true
    onlyReduceMaturationIfUnderperfused = True,
    isShearForceEffectLinear = False
  ),
)


newradiusdeflation = dict(
  num_threads = 2,
  out_intervall = 100,
  tumor_speed = 2.,
  tumor_radius = 200.,
  tissuePressureDistribution = 'shell',
  tissuePressureWidth = 500.,
  tissuePressureCenterFraction = 0.,
  vessels = dict(
    gfVessProl = 0.0005,
    timeProlEcSprout = 2,
    timeProlEcSproutLifetime = 50,
    timeEcRadiusInflate = 12,
    timeEcRadiusDeflate = 12,
    sproutDelay = 12,
    probCollapse = 1,
    bRelativeShearforceCollapse = False,
    forceCollapse = 0.1,
    forceEnlarge = 1.5,
    radMax = 12,
    distSproutMin = 8,
    dematuration_rate = 0.05,
    maturation_crit = 0.,
    bSproutModelSimple = False,
    vesselCompressionFactor = 1.,
    sproutMaxVesselWallThickness = 16,
    sproutMaxVesselWallThicknessArterial = 16., # replacement for arterial sprouting = true
    onlyReduceMaturationIfUnderperfused = True,
    isShearForceEffectLinear = True,
    bRadiusSmoothing = False,
    radiusSmoothingDiffusionCoeff = 1.,
    bShearStressControlledDilatationAndRadiusDependentCollapse = True,
    radMin = 2.4, # must bee lower than capillary radius (adjust if needed for initial networks)
    radInit = 2.6, # was 4
  ),
  calcflow = dict(
    viscosityPlasma = 1.2e-6,
    rheology = 'RheologyForHuman',
    inletHematocrit = 0.37,
    includePhaseSeparationEffect = False,
  )
)

newradiusdeflationInfinite = deepcopy(newradiusdeflation)
newradiusdeflationInfinite['vessels']['radMax'] = 200.

radiusdeflation2 = deepcopy(newradiusdeflation)
radiusdeflation2['vessels'].update(
  sproutMaxVesselWallThickness = 30.,
  sproutMaxVesselWallThicknessArterial = 5.,
  distSproutMin = 8,
  forceCollapse = 0.1,
  forceEnlarge = 3.0,
  radMax = 20,
  timeEcRadiusInflate = 24,
  timeEcRadiusDeflate = 24,
)

# this was used for breastv2 and breastv3 runs, with tumor sims labeled raddefl15
raddefl15 = deepcopy(radiusdeflation2)
raddefl15['vessels'].update(
  distSproutMin = 6,
  timeEcRadiusDeflate = 72, # not very sensitive against that
  forceEnlarge = 3.0, # lowering this *decreases the MVD
  radMax = 15,
)


test = deepcopy(circumferentialgrowth7)
test.update(
  tissuePressureDistribution = 'sphere',
  tissuePressureWidth = 100.,
  tumor_speed = 3.,
  num_threads = 2,
  tend = 1200.,
)
test['vessels'].update(
  vesselCompressionFactor = 0.7,
  radMax = 14.
)



def GenerateConfigurationsCompression(N):
  parameter_sets = []
  a = parameter_sets.append
  compressionFactors = np.random.permutation(myutils.random_stratified(0.2, 1.0, N))
  for cf in compressionFactors:
    c = deepcopy(circumferentialgrowth)
    c['vessels']['vesselCompressionFactor'] = cf
    c['vessels']['radMax'] = 25.
    c['vessels']['probCollapse'] = 0.05
    c['tissuePressureDistribution'] = 'sphere'
    c['tissuePressureWidth'] = 400.
    a(c)
  return parameter_sets  


def GenerateConfigurationsCompression2(N):
  parameter_sets = []
  a = parameter_sets.append
  compressionFactors = np.random.permutation(myutils.random_stratified(0.5, 1.0, N))
  for cf in compressionFactors:
    c = deepcopy(circumferentialgrowth7)
    c['vessels']['vesselCompressionFactor'] = cf
    c['vessels']['radMax'] = 14.
    c['tissuePressureDistribution'] = 'sphere'
    c['tissuePressureWidth'] = 400.
    a(c)
  return parameter_sets  
GenerateConfigurationsCompression2.name = 'cmpr2'


def GenerateConfigurationsRMAX(N):
  parameter_sets = []
  a = parameter_sets.append
  compressionFactors = np.random.permutation(myutils.random_stratified(0.5, 1.0, N))
  for cf in compressionFactors:
    c = deepcopy(circumferentialgrowth7)
    c['vessels']['radMax'] = 14. * cf
    a(c)
  return parameter_sets
GenerateConfigurationsRMAX.name = 'rmax'

""" PREZIOSI PARAMETER 

TODO
Attention needed to comment the O2 Parameter. 
The C++ Code is not net bulktissue and o2 compatible!

26.04.2016 --> now it is a little bit

"""

defaultconfig_bulktissue = dict(
  # this is the configuration which is used in submit-bulktissue-w-vessels
  lattice_size  = "set me to match the vessel domain",
  fn_out = "set_me",
  tend = 300.,
  out_intervall = 60,
  num_threads = cluster_threads,
  lattice_scale  = 30,

  #o2_rel_tumor_source_density = 0,
  #use_o2_source_decay = True,
  #o2_source_decay_time = 8,
  #adaption = appReal14['adaption'],
  prez_o2 = dict(
    o2_level_normal = 0.8,
    o2_range_tumor = 50,
    o2_range_necro = 100,
    o2_range_normal = 100,
    ),
  vessels = dict(
    gfVessProl = 0.0005,
    timeProlEcSprout = 2,
    #timeProlEcEnlarge = 180,
    #timeTillDeathEc = 50,
    sproutDelay = 24,
    probCollapse = 1,
    bRelativeShearforceCollapse = False,
    forceCollapse = 0.0005,
    #forceCollapseRel = 0.2,
    radMax = 25,
    radInit = 4,
    distSproutMin = 8,
    dematuration_rate = 0.05,
    maturation_crit = 0.,
    bSproutModelSimple = False,
    vesselCompressionFactor = 1.,
    sproutMaxVesselWallThickness = 50,
    sproutMaxVesselWallThicknessArterial = 50., # replacement for arterial sprouting = true
    onlyReduceMaturationIfUnderperfused = True,
    isShearForceEffectLinear = False
  ),
  
  calcflow = dict(
    viscosityPlasma = 1.2e-6,
    rheology = 'RheologyForHuman',
    inletHematocrit = 0.37,
    includePhaseSeparationEffect = False,
  ),

  tumor = dict(
    tumor_diameter = 500,
    init_shape = 'sphere',
    time_prol = 24,
    time_death = 240,
    time_necrosis = 48,
    time_death_tumor = 48000, # apoptosis

    ncells_norm = 0.4,
    ncells_tumor = 0.6,
    ncells_sigma = 0.1,
    ncells_ecm = 0.2,
    cell_mobility = 7.2e-1,
    cell_mobility_tumor = 7.2e-1, # micron^2 / (h kPa)
    cell_compression_modulus = 0.5e4, # kPa
    cell_mobility_ecmstar = -1,
    write_face_velocities = False,
    write_levelset_function = True,
    wrong_model = False,
    ecm_noise_density = -1,
    ecm_noise_std = -1,

    o2_necro_threshold = 0.03,
    o2_prol_threshold = 0.3,
    source_model_version = 3,
    use_necrotic_regions = True,

    surface_tension = -1.,
  ),
  simple_o2 = dict(
     o2_range_necro       = 100,
     o2_range_normal      = 100,
     o2_range_tumor       = 100,
     o2_cons_coeff_necro  = 0.0001,
     o2_cons_coeff_normal = 0.0001,
     o2_cons_coeff_tumor  = 0.0001,
     o2_rel_tumor_source_density = 0.2,
     o2_level_normal = 100,
     o2_source_decay_time = 8,
     hematocrit_init = 0.45,
     reference_intercapillary_distance = 80,
     capillary_wall_permeability = 0.013,
     test_obstacle = 0,
     use_o2_source_decay = 0,#bool
    ),
  reference_intercapillary_distance = 80,
  hematocrit_init = 0.45,
)
''' neat parameters for desktop pc'''
default = deepcopy(defaultconfig_bulktissue)
default.update(
    num_threads = 6,
    lattice_scale  = 50,
    out_intervall = 60, #1 min
    tend = 600., #5 min
)
default['tumor'].update(
    tumor_diameter = 250,
    )
merge = deepcopy(default)
merge['calcflow'].update(
    viscosityPlasma = 1.2e-6,
    rheology = 'RheologySecomb2005',
    inletHematocrit = 0.37,
    includePhaseSeparationEffect = False,
    )
bulktissue_small_2d_snowden = dict(
  # this is the configuration which is used in submit-bulktissue-w-vessels
  lattice_size  = "set me to match the vessel domain",
  fn_out = "set_me",
  tend = 2000.,
  out_intervall = 30,
  num_threads = cluster_threads,
  lattice_scale  = 30,

  #o2_rel_tumor_source_density = 0,
  #use_o2_source_decay = True,
  #o2_source_decay_time = 8,
  #adaption = appReal14['adaption'],
  prez_o2 = dict(
    o2_level_normal = 0.8,
    o2_range_tumor = 50,
    o2_range_necro = 100,
    o2_range_normal = 100,
    ),
  vessels = dict(
    gfVessProl = 0.0005,
    timeProlEcSprout = 2,
    #timeProlEcEnlarge = 180,
    #timeTillDeathEc = 50,
    sproutDelay = 24,
    probCollapse = 1,
    bRelativeShearforceCollapse = False,
    forceCollapse = 0.0005,
    #forceCollapseRel = 0.2,
    radMax = 25,
    radInit = 4,
    distSproutMin = 8,
    dematuration_rate = 0.05,
    maturation_crit = 0.,
    bSproutModelSimple = False,
    vesselCompressionFactor = 1.,
    sproutMaxVesselWallThickness = 50,
    sproutMaxVesselWallThicknessArterial = 50., # replacement for arterial sprouting = true
    onlyReduceMaturationIfUnderperfused = True,
    isShearForceEffectLinear = False
  ),
  
  calcflow = dict(
    viscosityPlasma = 1.2e-6,
    rheology = 'RheologySecomb2005',
    inletHematocrit = 0.37,
    includePhaseSeparationEffect = True,
  ),

  tumor = dict(
    tumor_diameter = 10,
    init_shape = 'sphere',
    time_prol = 24,
    time_death = 240,
    time_necrosis = 48,
    time_death_tumor = 48000, # apoptosis

    ncells_norm = 0.4,
    ncells_tumor = 0.6,
    ncells_sigma = 0.1,
    ncells_ecm = 0.2,
    cell_mobility = 7.2e-1,
    cell_mobility_tumor = 7.2e-1, # micron^2 / (h kPa)
    cell_compression_modulus = 0.5e4, # kPa
    cell_mobility_ecmstar = -1,
    write_face_velocities = False,
    write_levelset_function = True,
    wrong_model = False,
    ecm_noise_density = -1,
    ecm_noise_std = -1,

    o2_necro_threshold = 0.03,
    o2_prol_threshold = 0.3,
    source_model_version = 3,
    use_necrotic_regions = True,

    surface_tension = -1.,
  ),
)
oxy_try_2 = deepcopy(bulktissue_small_2d_snowden)
oxy_try_2.update(
    out_intervall = 100,
    tend = 200000.,
    num_threads = 16,
    )

bulktissue_tutorial = deepcopy(defaultconfig_bulktissue)
bulktissue_tutorial.update(
  tend = 101.,
  out_intervall = 50,
)
colorectal=deepcopy(defaultconfig_bulktissue)
colorectal.update(
  lattice_scale  = 50,
  num_threads = 8,
)
colorectal['tumor'].update(
  tumor_diameter = 500,
)
colorectal2 = deepcopy(colorectal)
colorectal2.update(
tend = 58000.,
)

colorectal2['calcflow'].update(
  rheology = 'RheologyForRats',
  includePhaseSeparationEffect = True,
)
prez2 = dict(
  # this is the configuration which is used in submit-bulktissue-w-vessels
  lattice_size  = "set me to match the vessel domain",
  fn_out = "set_me",
  tend = 1200.,
  out_intervall = 100,
  num_threads = cluster_threads,
  lattice_scale  = 30,

  #o2_rel_tumor_source_density = 0,
  #use_o2_source_decay = True,
  #o2_source_decay_time = 8,
  #adaption = appReal14['adaption'],
  prez_o2 = dict(
    o2_level_normal = 0.8,
    o2_range_tumor = 50,
    o2_range_necro = 100,
    o2_range_normal = 100,
    ),
  vessels = dict(
    gfVessProl = 0.0005,
    timeProlEcSprout = 2,
    #timeProlEcEnlarge = 180,
    #timeTillDeathEc = 50,
    sproutDelay = 24,
    probCollapse = 1,
    bRelativeShearforceCollapse = True,
    forceCollapse = 0.0005,
    radMax = 25,
    radInit = 2.5,
    distSproutMin = 8,
    dematuration_rate = 0.05,
    maturation_crit = 0.,
    bSproutModelSimple = False,
    vesselCompressionFactor = 1.,
    sproutMaxVesselWallThickness = 30,
    sproutMaxVesselWallThicknessArterial = 5., # replacement for arterial sprouting = true
    onlyReduceMaturationIfUnderperfused = True,
    isShearForceEffectLinear = True,
    bShearStressControlledDilatationAndRadiusDependentCollapse = True,
  ),
  
  calcflow = dict(
    viscosityPlasma = 1.2e-6,
    rheology = 'RheologyForHuman',
    inletHematocrit = 0.45,
    includePhaseSeparationEffect = False,
  ),

  tumor = dict(
    tumor_diameter = 500,
    init_shape = 'sphere',
    time_prol = 24,
    time_death = 240,
    time_necrosis = 48,
    time_death_tumor = 48000, # apoptosis

    ncells_norm = 0.4,
    ncells_tumor = 0.6,
    ncells_sigma = 0.1,
    ncells_ecm = 0.2,
    cell_mobility = 7.2e-1,
    cell_mobility_tumor = 7.2e-1, # micron^2 / (h kPa)
    cell_compression_modulus = 0.5e4, # kPa
    cell_mobility_ecmstar = -1,
    write_face_velocities = False,
    write_levelset_function = True,
    wrong_model = False,
    ecm_noise_density = -1,
    ecm_noise_std = -1,

    o2_necro_threshold = 0.03,
    o2_prol_threshold = 0.3,
    source_model_version = 3,
    use_necrotic_regions = True,

    surface_tension = -1.,
  ),
)
prez3 = deepcopy(prez2)
myutils.UpdateHierarchical(prez3,
dict( 
vessels = dict(
  forceEnlarge = 0.01,#0.01 is c++ default value
  forceCollapse = 0.2,
  probCollapse = 0.5,
  radInit = 3.5,
  #note: vessels with maturation above that threshold cannot collapse
  maturation_crit = 4., # micrometer
  radMin = 2.5,
  sproutMaxVesselWallThicknessArterial = 5., # replacement for arterial sprouting = true
  sproutMaxVesselWallThickness = 16.,
),
tumor = dict(
  tumor_diameter = 100,
)
))
prez4 = deepcopy(prez3)
myutils.UpdateHierarchical(prez4,
dict(
vessels = dict(
  radInit = 2.6,
  probCollapse = 1.0,
  #forceCollapse = 0.005,
  bShearStressControlledDilatationAndRadiusDependentCollapse = True,
  radMin = 2.5,
  isShearForceEffectLinear = True,
  bRelativeShearforceCollapse = True,
  forceCollapse = 1.1,
)
))
prez5 = deepcopy(prez3)
myutils.UpdateHierarchical(prez5,
dict(
vessels = dict(
  radInit = 2.6,
  maturation_crit = 2., # micrometer
  probCollapse = .5,
  forceCollapse = 0.0005,
  bShearStressControlledDilatationAndRadiusDependentCollapse = False,
  radMin = 2.5,
  isShearForceEffectLinear = False,
  bRelativeShearforceCollapse = False,
)
))
prez6 = deepcopy(prez3)
myutils.UpdateHierarchical(prez6,
dict(
tend = 24000.,
vessels = dict(
  radInit = 2.6,
  maturation_crit = 2., # micrometer
  probCollapse = .5,
  forceCollapse = 0.0005,
  bShearStressControlledDilatationAndRadiusDependentCollapse = False,
  radMin = 2.5,
  isShearForceEffectLinear = False,
  bRelativeShearforceCollapse = False,
)
))
prez_new=deepcopy(prez3)
prez_new.update(
    tend = 200,
    )
prez_new['tumor'].update(
    tumor_diameter = 100,
    time_prol = 1,
    )



# for tumor sim without vessels
realisticconfig_bulktissue = deepcopy(defaultconfig_bulktissue)
myutils.UpdateHierarchical(realisticconfig_bulktissue, dict(
  tumor = dict(
    cell_mobility = 7.2e-1,
    cell_mobility_tumor = 7.2e-1, # micron^2 / (h kPa)
    cell_compression_modulus = 0.5e4, # kPa
  )
))
trastuzumab_growth2 = deepcopy(defaultconfig_bulktissue)
#del trastuzumab_growth2['tumor_speed']
#del trastuzumab_growth2['tumor_radius']
trastuzumab_growth2.update(
  num_threads = 16,
  out_intervall = 30,
  lattice_scale  = 30,
)
trastuzumab_growth2['vessels'].update(
  probCollapse = 0.2,
  forceCollapse = 0.2e-3,
  onlyReduceMaturationIfUnderperfused = True
)
trastuzumab_growth2['tumor'] = dict(
  tumor_diameter = 100,
  init_shape = 'sphere',
  time_prol = 24,
  time_death = 240,
  time_necrosis = 48,
  time_death_tumor = 48000, # apoptosis

  ncells_norm = 0.4,
  ncells_tumor = 0.6,
  ncells_sigma = 0.1,
  ncells_ecm = 0.2,
  cell_mobility = 5000,
  cell_mobility_tumor = 5000,
  cell_mobility_ecmstar = -1,
  write_face_velocities = False,
  write_levelset_function = True,
  wrong_model = False,
  ecm_noise_density = -1,
  ecm_noise_std = -1,

  o2_necro_threshold = 0.03,
  o2_prol_threshold = 0.3,
  source_model_version = 3,
  use_necrotic_regions = True,
  surface_tension = -1.,
)
trastuzumab_growth3 = deepcopy(trastuzumab_growth2)
trastuzumab_growth3.update(
  lattice_scale  = 80,
)
raddefl_bulktissue = deepcopy(defaultconfig_bulktissue)
raddefl_bulktissue['vessels'] = deepcopy(raddefl15['vessels'])
myutils.UpdateHierarchical(raddefl_bulktissue, dict(
  vessel_volume_exclusion = False,
  reference_intercapillary_distance = 100,
  fake_height_2d = 100.,
  gf_production_threshold = 0.3,  # setting this too low was an error
  tumor = dict(
    cell_mobility = 5.e-1,
    cell_mobility_tumor = 5.e-1, # micron^2 / (h kPa)
    cell_compression_modulus = 0.5e4, # kPa
  )
))


adaption_master = dict(
  num_threads = 6,
  out_intervall = 10,
  tumor_speed = 2.,
  tumor_radius = 25.,
  apply_adaption_intervall = 10,
  vessels = dict(
    forceCollapse = 0.25e-3, # in kPa
    probCollapse = 1.0,
    distSproutMin = 8, # in number of lattice sites
    radMax = 20., # micrometer
    #radius under max pressure is = reference_r * vesselCompressionFactor
    vesselCompressionFactor = 1.,
    gfVessProl = 0.0005, # unitless concentration
    timeProlEcSprout  = 2,  # 2h per added segment
    timeEcRadiusInflate = 144, # after 360 hours, vessel radius should be ca 27micron
    timeProlEcSproutLifetime = 50,  # the time until sprouts regress
    sproutDelay = 12, # the time during which a tumor vessel can still generate sprouts (hours)
    bRelativeShearforceCollapse = False,
    radInit = 2.6,
    #images in (holash et al.) indicate approximately 0.04,
    #based on the time it takes to destroy a 25micron *radius* vessel
    #maturation of 25 micron vessel ca 16, 5 micron vessel ca 4.5
    dematuration_rate = 0.05, # micrometer per hour
    #initial wall thickness should be around 4 micron, so this is maturation_crit, so capillaries should immediately collapse
    #note that also vessels within normal tissue can collapse if their maturation value is < maturation_crit
    maturation_crit = 0., # micrometer
    #bArterialSprouting = True, # old code
    sproutMaxVesselWallThicknessArterial = 50., # replacement for arterial sprouting = true
    sproutMaxVesselWallThickness = 50.,
    bShearStressControlledDilatationAndRadiusDependentCollapse = False,
    isShearForceEffectLinear = False,
    bRadiusSmoothing = False,
    radiusSmoothingDiffusionCoeff = 1.,
    onlyReduceMaturationIfUnderperfused = True,
    badaption_on_off = True,
  ),
  calcflow = dict(
    viscosityPlasma = 1.2e-6,
    rheology = 'RheologySecomb2005',
    inletHematocrit = 0.37,
    includePhaseSeparationEffect = False,
  ),
  adaption = dict(
    k_m = 0.900009,
    k_c = 2.616314,
    k_s = 1.826391,
    cond_length = 2452.249030,
    Q_refdot = 40,
    S_0 = 20,
    #if this is 0 we iterate until other convergenze criteria is reached
    max_nun_iterations = 100,
    qdev = .0,
    #if starting_radii is 0. we use the values given in
    #the input file
    starting_radii = 0.,
    delta_t = 0.0,
    tum_manitulate_s1 = False,
    tum_manitulate_s2 = False,
    tum_manitulate_s3 = False,
    tum_manitulate_s4 = False,
    tum_manitulate_s5 = False,
  ),
)
adaption_master_phase = deepcopy(adaption_master)
adaption_master_phase['calcflow'].update(
  includePhaseSeparationEffect = True,
)

little_tumor_typeE_vary1 = dict(
  num_threads = 6,
  out_intervall = 20,
  tumor_speed = 2.,
  tumor_radius = 50.,
  apply_adaption_intervall = 10,
  vessels = dict(
    forceCollapse = 0.25e-3, # in kPa
    probCollapse = 1.0,
    distSproutMin = 8, # in number of lattice sites
    radMax = 20., # micrometer
    #radius under max pressure is = reference_r * vesselCompressionFactor
    vesselCompressionFactor = 1.,
    gfVessProl = 0.0005, # unitless concentration
    timeProlEcSprout  = 2,  # 2h per added segment
    timeEcRadiusInflate = 144, # after 360 hours, vessel radius should be ca 27micron
    timeProlEcSproutLifetime = 50,  # the time until sprouts regress
    sproutDelay = 12, # the time during which a tumor vessel can still generate sprouts (hours)
    bRelativeShearforceCollapse = False,
    radInit = 2.6,
    #images in (holash et al.) indicate approximately 0.04,
    #based on the time it takes to destroy a 25micron *radius* vessel
    #maturation of 25 micron vessel ca 16, 5 micron vessel ca 4.5
    dematuration_rate = 0.05, # micrometer per hour
    #initial wall thickness should be around 4 micron, so this is maturation_crit, so capillaries should immediately collapse
    #note that also vessels within normal tissue can collapse if their maturation value is < maturation_crit
    maturation_crit = 0., # micrometer
    #bArterialSprouting = True, # old code
    sproutMaxVesselWallThicknessArterial = 50., # replacement for arterial sprouting = true
    sproutMaxVesselWallThickness = 50.,
    bShearStressControlledDilatationAndRadiusDependentCollapse = False,
    isShearForceEffectLinear = False,
    bRadiusSmoothing = False,
    radiusSmoothingDiffusionCoeff = 1.,
    onlyReduceMaturationIfUnderperfused = True,
    badaption_on_off = True,
  ),
  calcflow = dict(
    viscosityPlasma = 1.2e-6,
    rheology = 'RheologySecomb2005',
    inletHematocrit = 0.37,
    includePhaseSeparationEffect = False,
  ),
  adaption = dict(
    k_m = 0.900009,
    k_c = 2.616314,
    k_s = 1.826391,
    cond_length = 2452.249030,
    Q_refdot = 40,
    S_0 = 20,
    #if this is 0 we iterate until other convergenze criteria is reached
    max_nun_iterations = 100,
    qdev = .0,
    #if starting_radii is 0. we use the values given in
    #the input file
    starting_radii = 0.,
    delta_t = 0.0,
    tum_manitulate_s1 = False,
    tum_manitulate_s2 = False,
    tum_manitulate_s3 = False,
    tum_manitulate_s4 = False,
    tum_manitulate_s5 = False,
  ),
)