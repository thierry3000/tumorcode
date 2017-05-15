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
import myutils
cluster_threads = myutils.cluster_threads

from copy import deepcopy

default = dict(
  num_threads = 6,
  out_intervall = 100,
  tumor_speed = 1.42,
  tumor_radius = 50.,
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
    includePhaseSeparationEffect = False,
  )
)
newradiusdeflation = deepcopy(default)
newradiusdeflation['num_threads'] = cluster_threads

gero_3month_to_5mm = deepcopy(newradiusdeflation)
gero_3month_to_5mm.update(
    tumor_speed = 0.001,
    tumor_radius = 20.,
    )
