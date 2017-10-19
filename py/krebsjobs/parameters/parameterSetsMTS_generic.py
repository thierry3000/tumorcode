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

from copy import deepcopy
import myutils
cluster_threads = myutils.cluster_threads
import parameterSetsO2

default = dict(
  num_threads = 6,
  out_intervall = 100,
  dt = 1,
  message = '',
  rGf = 200.,
  gf_production_threshold = 0.1,
  stopping_radius_fraction = 0.6,
  tumor_speed = 1.42,
  tumor_radius = 50.,
  tissuePressureDistribution = 'shell',
  tissuePressureWidth = 500.,
  tissuePressureCenterFraction = 0.,
  lattice_scale = 5.,
  lattice_size = 25,
  vessels = dict(
    gfVessProl = 0.0005,
    timeProlEcSprout = 2,
    timeProlEcSproutLifetime = 50,
    timeEcRadiusInflate = 12,
    timeEcRadiusDeflate = 12,
    sproutDelay = 12,
    probCollapse = 0.5,
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
    includePhaseSeparationEffect = 0,
  )
)
'''works as test, creates 4-5 cells '''
milotti_mts_test = deepcopy(default)
milotti_mts_test.update(
    out_intervall = 1,
    tend=200.,
    tumor_speed = 0.1,
    tumor_radius = 5.,
    lattice_scale = 50.,
    lattice_size = 20,
    tissuePressureDistribution = 'sphere',
    vessels = dict(
        timeProlEcSprout = 2,
        timeProlEcSproutLifetime = 50,
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
        )
    )
milotti_detailed = deepcopy(milotti_mts_test)
milotti_detailed['detailedo2'] = parameterSetsO2.milotti_o2
milotti_detailed.update(
    rGf = 100.,
    tend= 200, # usually starts to get lame over 200
    )

milotti_simple = deepcopy(milotti_detailed)
milotti_simple['detailedo2'].update(
    parameterSetsO2.milotti_o2_simple
    )
milotti_feedback = deepcopy(milotti_simple)
milotti_feedback.update(
    rGf = 100.,
    tend= 860000, # about 10 days
    )
milotti_test = deepcopy(milotti_detailed)
milotti_test['detailedo2'] = parameterSetsO2.milotti_o2_test
#milotti_detailed['detailedo2'] = parameterSetsO2.default_o2
#factory = getattr(parameterSets, tumorParameterName)
#milotti_detailed['detailedo2'].update(
#    max_iter = 500,
#    )
plos_One_2016 = dict(
    tumor_speed = 2.,
    tumor_radius = 250.,
    rGF = 200.,
    stopping_radius_fraction = 0.6,
    tend = 1000,
    tissuePressureCenterFraction = 0,
    tissuePressureDistribution = 'sphere',
    tissuePressureWith = 500,
    num_threads = 6,
    vessels = dict(
        bRadiusSmoothing = False,
        bRelativeShearforceCollapse = False,
        bShearStressControlledDilatationAndRadiusDependentCollapse = False,
        bSproutModelSimple = False,
        dematuration_rate = 0.05,
        distSproutMin = 8,
        forceCollapse = 0.00025,
        forceEnlarge = 0.01,
        gfVessProl = 0.0005,
        isShearForceEffectLinear = False,
        maturation_crit = 0,
        onlyReduceMaturationIfUnderperfused = True,
        pressMax = 13,
        pressMin = 0,
        probCollapse = 1,
        radInit = 2.6,
        radMax = 8.14502,
        radMin = 2.9,
        radiusSmoothingDiffusionCoeff = 1,
        seed = 12348,
        sproutDelay = 12,
        sproutMaxVesselWallThickness = 50,
        sproutMaxVesselWallThicknessArterial = 50,
        timeEcRadiusDeflate = 60,
        timeEcRadiusInflate = 144,
        timeProlEcSprout = 2,
        timeProlEcSproutLifetime = 50,
        vesselCompressionFactor = 1,
        ),
    calcflow = dict(
        includePhaseSeparationEffect=0,
        inletHematocrit = 0.37,
        rheology = 'RheologyForHuman',
        viscosityPlasma = 1.2e-6,
        )
    )
milotti_mts_1 = deepcopy(plos_One_2016)
milotti_mts_1.update(
    tumor_speed = 0.2,
    out_intervall = 100,
    tend=1400.,
    tissuePressureDistribution = 'sphere',
    )
forCatSim = deepcopy(default)
forCatSim.update(
  num_threads = cluster_threads,
  tumor_speed = 1.042,
  tumor_radius = 50.,
)
forCatSim['rheology'] = 'RheologySecomb2005'

#this creates a nice tumor of size 5.3 mm
#but in 1 month time --> experiments last 3 months
gero_3month_to_5mmb = deepcopy(default)
gero_3month_to_5mmb.update(
    num_threads = cluster_threads,
    out_intervall = 43200, # 3600*12 half a day
    tumor_speed = 0.001,
    tumor_radius = 20.,
    tend=7000000.,
    dt = 1000.,
    )
gero_3month_to_5mmc = deepcopy(default)
gero_3month_to_5mmc.update(
    num_threads = cluster_threads,
    out_intervall = 43200, # 3600*12 half a day
    tumor_speed = 0.003,
    tumor_radius = 20.,
    tend=7000000.,
    dt = 3000.,
    )
gero_3month_to_5mmd = deepcopy(default)
gero_3month_to_5mmd.update(
    num_threads = cluster_threads,
    out_intervall = 43200, # 3600*12 half a day
    tumor_speed = 0.000333333333,
    tumor_radius = 20.,
    tend=7000000.,
    dt = 3600.,
    )
if __name__ == '__main__':
  print(milotti_detailed)