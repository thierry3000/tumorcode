#!/usr/bin/env python2
# -*- coding: utf-8 -*-
'''
This file is part of tumorcode project.
(http://www.uni-saarland.de/fak7/rieger/homepage/research/tumor/tumor.html)

Copyright (C) 2018  Thierry Fredrich

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
import parameterSetsFakeTumor
import parameterSetsVesselGen

default = parameterSetsFakeTumor.default
default.update(
  out_intervall = 100,
  dt = 1,
  message = '',
  rGf = 200.,
  rO2Consumtion = 100.,
  gf_production_threshold = 0.1,
  stopping_radius_fraction = 0.6,
  lattice_scale = 5.,
  lattice_size = 25,
  useConstO2 = True,
)
default['detailedo2'] = parameterSetsO2.milotti_o2

'''works as test, creates 4-5 cells '''
milotti_mts = deepcopy(default)
milotti_mts['detailedo2'] = parameterSetsO2.milotti_o2
milotti_mts['detailedo2'].update(
    useCellBasedUptake=True
    )
milotti_mts['vessels'] = parameterSetsFakeTumor.milotti_mts_3['vessels']
milotti_mts.update(
    out_intervall = 1,
    tumor_speed = 0.2,
    tumor_radius = 5.,
    lattice_scale = 50.,
    lattice_size = 20,
    tissuePressureDistribution = 'sphere',
    rGf = 100.,
    rO2Consumtion = 10.,
    tend=1400.,
    useConstO2 = False,
    )
milotti_mts = deepcopy(milotti_mts)
milotti_mts['detailedo2'].update(
    solubility_tissue = 2.8e-5,
    extra_tissue_source_const = 0.0,
    )
milotti_mts_2 = deepcopy(milotti_mts)
milotti_mts_2.update(
    useConstO2=True,
    )
#note: one needs more than 5 iterations for the bisect to converge
'''
error in function boost::math::tools::bisect<d>: No change of sign in 
boost::math::tools::bisect, either there is no root to find, 
or there are multiple roots in the interval (f(min) = 7.2726904991920653)
'''
milotti_mts_3 = deepcopy(milotti_mts)
milotti_mts_3['detailedo2'].update(
    max_iter = 1,
    input_group_path = "vessels",
    )
milotti_mts_3.update(
    tend = 100,
    useConstO2 = True,
    max_iteration_per_rerun = 3,
    )

milotti_mts_4 = deepcopy(milotti_mts)
milotti_mts_4['detailedo2'].update(
    rO2Consumtion = 5.,
    )
milotti_mts_5 = deepcopy(milotti_mts)
milotti_mts_5.update(
    rO2Consumtion = 30.,
    )
milotti_mts_6 = deepcopy(milotti_mts)
milotti_mts_6.update(
    tumor_radius = 10.,
    tumor_speed = 1.,
    out_intervall = 1,
    tend=1000.,
    )
milotti_mts_6['vessels'].update(
    radInit = 2.6,
    radMin = 3.5,
    probCollapse = 0.1,
    forceCollapse = 0.00025,
    bShearStressControlledDilatationAndRadiusDependentCollapse = False,
)
milotti_mts_7 = deepcopy(milotti_mts_6)
milotti_mts_7['useTumorcodeVessels'] = True
milotti_mts_7['detailedo2'].update(
    max_iter = 50,
    loglevel=1,
    useCellBasedUptake=True,
    michaelis_menten_uptake = True,
    )
milotti_mts_7b = deepcopy(milotti_mts_7)
milotti_mts_7b.update(
	tend=9999.,
)
milotti_mts_7c = deepcopy(milotti_mts_7b)
milotti_mts_7c.update(
	rGf = 1.,
	gf_production_threshold = 0.01,
)
milotti_mts_8 = deepcopy(milotti_mts_7)
#milotti_mts_8['useTumorcodeVessels'] = False
milotti_mts_8.update(
    useConstO2 = True,
    )
milotti_mts_8b= deepcopy(milotti_mts_8)
milotti_mts_8b.update(
    tend=9999.,
    )
milotti_mts_8c= deepcopy(milotti_mts_8b)
milotti_mts_8c.update(
    rGf = 1.,
    gf_production_threshold = 0.01,
    )
''' changed g_env[k] to g_env[i]
    old geometry as sent from milotti by email
'''
'''
no const o2, skip convergence
'''
milotti_mts_9 = deepcopy(milotti_mts_7)
milotti_mts_9['detailedo2'].update(
    max_iter = 200,
    convergence_tolerance = 1e-2,
    )
milotti_mts_10 = deepcopy(milotti_mts_7)
milotti_mts_10.update(
    rO2Consumtion = 5.,
    )
milotti_mts_11 = deepcopy(milotti_mts_7)
milotti_mts_11.update(
    rO2Consumtion = 30.,
    )
milotti_mts_12 = deepcopy(milotti_mts_7)
milotti_mts_12.update(
    rO2Consumtion = 50.,
    )
milotti_mts_13 = deepcopy(milotti_mts_7c)
milotti_mts_13.update(
	rGf = 1.,
	gf_production_threshold = 1e8,
)
milotti_mts_14 = deepcopy(milotti_mts_13)
milotti_mts_14.update(
	useConstO2 = True,
)
milotti_mts_14['vessels'].update(
	gfVessProl = 1e8,
)
''' NOTE: lattice_size is not used for MTS
'''
vbl = dict(
    tissuePressureDistribution = 'sphere',
    gf_production_threshold =  0.01,
    lattice_size = 20,
    lattice_scale = 50.0,
    tumor_speed = 1.0,
    useTumorcodeVessels = True,
    rGf = 1.0,
    dt = 1,
    tend = 9999.0,
    max_iteration_per_rerun = 9999,
    message ='',
    tumor_radius = 10.0,
    rO2Consumtion = 10.0,
    tissuePressureWidth = 500.0,
    tissuePressureCenterFraction = 0.0,
    stopping_radius_fraction = 0.6,
    useConstO2 = False,
    out_intervall = 1,
    detailedo2 = dict(
        solubility_plasma = 3.1e-05,
        solubility_tissue = 2.8e-05,
        po2_mmcons_k_norm = 4.0,
        po2_mmcons_k_tum = 2.0,
        po2_mmcons_k_necro = 2.0,
        po2_mmcons_m0_norm = 6.1e-05,
        po2_mmcons_m0_tum = 0.000244,
        po2_mmcons_m0_necro = 0.0,
        po2init_r0 = 55.0,
        po2init_dr = 1.0,
        po2init_cutoff = 100.0,
        massTransferCoefficientModelNumber = 1,
        rd_necro = 150.0,
        rd_norm = 150.0,
        rd_tum = 50.0,
        D_plasma = 2750.0,
        extra_tissue_source_linear = 0.0,
        extra_tissue_source_const = 0.0,
        tissue_po2_boundary_condition = 'neumann',
        tissue_boundary_value = 0.0,
        sat_curve_p50 = 27.0,
        sat_curve_exponent = 2.7,
        haemoglobin_binding_capacity = 0.5,
        detailedO2name = 'vbl_o2',
        axial_integration_step_factor = 0.25,
        michaelis_menten_uptake = True,
        approximateInsignificantTransvascularFlux = True,
        conductivity_coeff1 = 8,
        conductivity_coeff2 = 4.7,
        conductivity_coeff3 = 0, 
        max_iter = 50,
        transvascular_ring_size = 0.5,
        c0 = 0.5,
        convergence_tolerance = 0.001,
        num_threads = 4,
        loglevel = 1,
        debug_zero_o2field = False,
        useCellBasedUptake = True,
        calcflow = dict(
            viscosityPlasma = 1.2e-06,
            rheology = 'RheologySecomb2005',
            includePhaseSeparationEffect = 1,
            inletHematocrit =  0.45
            ),
        debug_fn = 'none.h5',
        grid_lattice_const = 30.0,
        safety_layer_size = 60.0,
        input_group_path = 'vessels',
        ),
    vessels = dict(
        bRadiusSmoothing = False,
        bRelativeShearforceCollapse = False,
        bShearStressControlledDilatationAndRadiusDependentCollapse = False,
        bSproutModelSimple = False,
        dematuration_rate = 0.05,
        distSproutMin = 8,
        forceCollapse = 0.00025,
        forceEnlarge = 0.01,
        gfVessProl = 0.01,
        isShearForceEffectLinear = False,
        maturation_crit = 0,
        onlyReduceMaturationIfUnderperfused = True,
        pressMax = 13,
        pressMin = 0,
        probCollapse =  0.1,
        radInit = 2.6,
        radMax = 8.14502,
        radMin = 3.5,
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
#    calcflow = dict(
#        viscosityPlasma = 1.2e-06,
#        rheology = 'RheologyForHuman',
#        includePhaseSeparationEffect = False,
#        inletHematocrit = 0.37
#        ),
)
vbl['calcflow'] = vbl['detailedo2']['calcflow']

vbl_disable_vessel_remodelling = deepcopy(vbl) #former milotti_mts_14
vbl_disable_vessel_remodelling['vessels']['gfVessProl'] = 1e8

vbl_dvr_const_o2 = deepcopy(vbl_disable_vessel_remodelling)
vbl_dvr_const_o2['useConstO2'] = True

vbl_const_o2 = deepcopy(vbl)
vbl_const_o2['useConstO2'] = True

if __name__ == '__main__':
  #print(milotti_mts_13)
  aDict = milotti_mts_7c['vessels']
  print(vbl)
