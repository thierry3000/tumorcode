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
  gf_production_threshold = 0.1,
  stopping_radius_fraction = 0.6,
  lattice_scale = 5.,
  lattice_size = 25,
)
default['detailedo2'] = parameterSetsO2.milotti_o2

'''works as test, creates 4-5 cells '''
milotti_mts = deepcopy(default)
milotti_mts['detailedo2'] = parameterSetsO2.milotti_o2
milotti_mts.update(
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
milotti_detailed = deepcopy(milotti_mts)
milotti_detailed.update(
    rGf = 100.,
    tend= 200, # usually starts to get lame over 200
    )


milotti_mts_1 = deepcopy(default)
milotti_mts_1.update(
    rGf = 100.,
    tumor_speed = 0.2,
    out_intervall = 1,
    tend=1400.,
    tissuePressureDistribution = 'sphere',
    )
milotti_mts_1['detailedo2'] = parameterSetsO2.milotti_o2
milotti_mts_2 = deepcopy(milotti_mts_1)
milotti_mts_2['calcflow'] = parameterSetsVesselGen.no_phase
if __name__ == '__main__':
  print(milotti_detailed)