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
# -*- coding: utf-8 -*-
import numpy as np
from copy import deepcopy
import math

default = dict(
  num_threads = 8,
  adaption = dict(
    k_c = 1.5,
    k_m = 1.5,
    k_s = 1.,
    Q_refdot = 40,
    S_0 = 20,
    cond_length = 2500.,
    #if this is 0 we iterate until qdev is reached
    max_nun_iterations = 10,
    qdev = 0.01,
    #if starting_radii is 0. we use the values given in
    #the input file
    starting_radii = 0.,
    # 0.0 means the time step is dynamically chosen!
    delta_t = 0.0, 
    ),
  calcflow = dict(
    viscosityPlasma = 1.2e-6, #commented means using default for rats
    rheology = 'RheologyForRats',
    inletHematocrit = 0.40,
    includePhaseSeparationEffect =  False,
  ),
)
default_phase=deepcopy(default)
default_phase['calcflow'].update(
  includePhaseSeparationEffect =  True,
  rheology = 'RheologySecomb2005',
)
adaption_default = dict(
  adaption = dict(
    k_c = 2.11,
    k_m = 2.5,
    k_s = 1.7,
    Q_refdot = 20,
    S_0 = 10,
    cond_length = 2000.,
    #if this is 0 we iterate until qdev is reached
    max_nun_iterations = 10,
    qdev = 0.01,
    #if starting_radii is 0. we use the values given in
    #the input file
    starting_radii = 0.,
    delta_t = 0.001,
    ),
  calcflow = dict(
    viscosityPlasma = 1.2e-6, #commented means using default for rats
    rheology = 'RheologyForRats',
    inletHematocrit = 0.40,
    includePhaseSeparationEffect =  False,
  ),
)

adaption_typeE = dict(
  adaption = dict(
    k_c = 2.11,
    k_m = 2.5,
    k_s = 1.7,
    Q_refdot = 40,
    S_0 = 10,
    cond_length = 2000.,
    #if this is 0 we iterate until qdev is reached
    max_nun_iterations = 3000,
    qdev = 0.001,
    #if starting_radii is 0. we use the values given in
    #the input file
    starting_radii = 0.,
    delta_t = 0.01,
    ),
  calcflow = dict(
    viscosityPlasma = 1.2e-6, #commented means using default for rats
    rheology = 'RheologyForRats',
    inletHematocrit = 0.40,
    includePhaseSeparationEffect =  False,
  ),
)

size_v_typeE = dict(
  adaption = dict(
    k_c = 2.11,
    k_m = 2.5,
    k_s = 1.7,
    Q_refdot = 40,
    S_0 = 0.5, #guess log(no edges/no roots)
    cond_length = 2000.,
    #if this is 0 we iterate until qdev is reached
    max_nun_iterations = 3,
    qdev = 0.001,
    #if starting_radii is 0. we use the values given in
    #the input file
    starting_radii = 0.,
    delta_t = 0.01,
    ),
  calcflow = dict(
    viscosityPlasma = 1.2e-6, #commented means using default for rats
    rheology = 'RheologyForRats',
    inletHematocrit = 0.40,
    includePhaseSeparationEffect =  False,
  ),
)
size_v_typeEb = dict(
  adaption = dict(
    k_c = 1.5,
    k_m = 1.2,
    k_s = 1.7,
    Q_refdot = 40,
    S_0 = 0.5, #guess log(no edges/no roots)
    cond_length = 200.,
    #if this is 0 we iterate until qdev is reached
    max_nun_iterations = 2000,
    qdev = 0.02,
    #if starting_radii is 0. we use the values given in
    #the input file
    starting_radii = 10.,
    delta_t = 0.01,
    ),
  calcflow = dict(
    viscosityPlasma = 1.2e-6, #commented means using default for rats
    rheology = 'RheologySecomb2005',
    inletHematocrit = 0.40,
    includePhaseSeparationEffect =  True,
  ),
)
  
adaption_PLOS_2009 = dict(
  adaption = dict(
    k_c = 2.113,
    k_m = 0.6444,
    k_s = 1.69,
    Q_refdot = 1.98,
    S_0 = 20,
    cond_length = 2000.,
    #if this is 0 we iterate until qdev is reached
    max_nun_iterations = 500,
    qdev = 0.1,
    #if starting_radii is 0. we use the values given in
    #the input file
    starting_radii = 0.,
    delta_t = 0.01,
    ),
  calcflow = dict(
    viscosityPlasma = 1.2e-6, #commented means using default for rats
    rheology = 'RheologyForRats',
    inletHematocrit = 0.40,
    includePhaseSeparationEffect =  False,
  ),
)
  
adaption_typeB = dict(
  adaption = dict(
    k_c = 2.113,
    k_m = 0.6444,
    k_s = 1.69,
    Q_refdot = .98,
    S_0 = 20,
    cond_length = 2000.,
    #if this is 0 we iterate until qdev is reached
    max_nun_iterations = 500,
    qdev = 0.1,
    #if starting_radii is 0. we use the values given in
    #the input file
    starting_radii = 0.,
    delta_t = 0.01,
    ),
  calcflow = dict(
    viscosityPlasma = 6.e-6, #commented means using default for rats
    rheology = 'RheologyForRats',
    inletHematocrit = 0.40,
    includePhaseSeparationEffect =  False,
  ),
)
adaption_typeC = dict(
  adaption = dict(
    k_c = 2.113,
    k_m = 0.6444,
    k_s = .929,
    Q_refdot = .098,
    S_0 = 20,
    cond_length = 1500.,
    #if this is 0 we iterate until qdev is reached
    max_nun_iterations = 500,
    qdev = 0.1,
    #if starting_radii is 0. we use the values given in
    #the input file
    starting_radii = 0.,
    delta_t = 0.01,
    ),
  calcflow = dict(
    viscosityPlasma = 6.e-6, #commented means using default for rats
    rheology = 'RheologySecomb2005',
    inletHematocrit = 0.40,
    includePhaseSeparationEffect =  True,
  ),
)
  
apj_1 = dict(
  adaption = dict(
    k_c = 2.45,
    k_m = 0.7,
    k_s = 1.72,
    Q_refdot = 0.198,
    S_0 = 27.9,
    cond_length = 173000.,
    max_nun_iterations = 2,
    qdev = 0.00001,
    starting_radii = 3.,
    delta_t = 0.01,
    ),
  calcflow = dict(
    viscosityPlasma = 1.2e-6, #commented means using default for rats
    rheology = 'RheologyForRats',
    inletHematocrit = 0.40,
    includePhaseSeparationEffect =  True,
  ),
)
apj_2 = dict(
  adaption = dict(
    k_c = 0.7,
    k_m = 0.83,
    k_s = 1.79,
    Q_refdot = 40,
    S_0 = 20,
    cond_length = 1500.,
    max_nun_iterations = 0,
    qdev = 0.00001,
    starting_radii = 0.,
    delta_t = 0.01,
    ),
  calcflow = dict(
    viscosityPlasma = 1.2e-6, #commented means using default for rats
    rheology = 'RheologyForRats',
    inletHematocrit = 0.40,
    includePhaseSeparationEffect =  True,
  ),
)

adaption_symetric_world = dict(
  adaption = dict(
    k_c = 1.61,
    k_m = 0.5,
    k_s = 1.16,
    Q_refdot = 40.5,
    S_0 = 5,
    cond_length = 1700.,
    max_nun_iterations = 1000,
    qdev = 0.005,
    starting_radii = 0.,
    delta_t = .001,
    ),
  calcflow = dict(
    viscosityPlasma = 1.2e-6, #commented means using default for rats
    rheology = 'RheologyForRats',
    inletHematocrit = 0.40,
    includePhaseSeparationEffect =  True,
  ),
)

adaption_symetricA_world = dict(
  adaption = dict(
    k_c = 0.0,
    k_m = 0.0,
    k_s = 0.0,
    Q_refdot = 40.5,
    S_0 = 5,
    cond_length = 1700.,
    max_nun_iterations = 4000,
    qdev = 0.001,
    starting_radii = 10.,
    delta_t = .001,
    ),
  calcflow = dict(
    viscosityPlasma = 1.2e-6, #commented means using default for rats
    rheology = 'RheologyForRats',
    inletHematocrit = 0.40,
    includePhaseSeparationEffect =  True,
  ),
)
adaption_symetricA_world_break = dict(
  adaption = dict(
    k_c = 1.61,
    k_m = 1.5,
    k_s = 0.1,
    Q_refdot = 20.5,
    S_0 = 5,
    cond_length = 700.,
    max_nun_iterations = 100,
    qdev = 0.01,
    starting_radii = 0.,
    delta_t = .01,
    ),
  calcflow = dict(
    viscosityPlasma = 1.2e-6, #commented means using default for rats
    rheology = 'RheologyForRats',
    inletHematocrit = 0.40,
    includePhaseSeparationEffect =  True,
  ),
)
  
adaption_symetricB_world = dict(
  adaption = dict(
    k_c = 1.61,
    k_m = 2.5,
    k_s = 1.16,
    Q_refdot = 10.5,
    S_0 = 5,
    cond_length = 1700.,
    max_nun_iterations = 1000,
    qdev = 0.01,
    starting_radii = 10.,
    delta_t = .101,
    ),
  calcflow = dict(
    viscosityPlasma = 1.2e-6, #commented means using default for rats
    rheology = 'RheologyForRats',
    inletHematocrit = 0.40,
    includePhaseSeparationEffect =  True,
  ),
)
  
adaption_symetricIrregular_world = dict(
  adaption = dict(
    k_c = 0.0,
    k_m = 0.8,
    k_s = .57,
    Q_refdot = 40.5,
    S_0 = 5,
    cond_length = 1700.,
    max_nun_iterations = 3000,
    qdev = 0.01,
    starting_radii = 0.,
    delta_t = .01,
    ),
  calcflow = dict(
    viscosityPlasma = 1.2e-6, #commented means using default for rats
    rheology = 'RheologySecomb2005',
    inletHematocrit = 0.40,
    includePhaseSeparationEffect =  True,
  ),
)
adaption_symetricIrregular_world_bak = dict(
  adaption = dict(
    k_c = 0.0,
    k_m = 0.8,
    k_s = .57,
    Q_refdot = 40.5,
    S_0 = 5,
    cond_length = 1700.,
    max_nun_iterations = 1000,
    qdev = 0.01,
    starting_radii = 0.,
    delta_t = .01,
    ),
  calcflow = dict(
    viscosityPlasma = 1.2e-6, #commented means using default for rats
    rheology = 'RheologySecomb2005',
    inletHematocrit = 0.40,
    includePhaseSeparationEffect =  True,
  ),
)


adaption_asymetric_world = dict(
  adaption = dict(
    k_c = 2.8,
    k_m = 0.8,
    k_s = 1.75,
    Q_refdot = 40,
    S_0 = 20,
    cond_length = 1400.,
    max_nun_iterations = 00,
    qdev = 0.001,
    starting_radii = 0.,
    delta_t = .001,
    ),
  calcflow = dict(
    viscosityPlasma = 1.2e-6, #commented means using default for rats
    rheology = 'RheologyForRats',
    inletHematocrit = 0.40,
    includePhaseSeparationEffect =  True,
  ),
)
adaption_asymetric_world_bak = dict(
  adaption = dict(
    k_c = 1.61,
    k_m = 1.91,
    k_s = 1.6,
    Q_refdot = 40,
    S_0 = 20,
    cond_length = 1200.,
    max_nun_iterations = 1000,
    qdev = 0.0005,
    starting_radii = 0.,
    delta_t = .01,
    ),
  calcflow = dict(
    viscosityPlasma = 1.2e-6, #commented means using default for rats
    rheology = 'RheologyForRats',
    inletHematocrit = 0.40,
    includePhaseSeparationEffect =  True,
  ),
)
  
mesentry_subset = dict(
  adaption = dict(
    k_c = 1.0834,
    k_m = 0.73,
    k_s = 1.399449,
    Q_refdot = 20,
    S_0 = 10,
    cond_length = 1500.,
    max_nun_iterations = 0,
    qdev = 10.1,
    starting_radii = 0.,
    delta_t = .001,
    ),
  calcflow = dict(
    viscosityPlasma = 1.2e-6, #commented means using default for rats
    rheology = 'RheologyForRats',
    inletHematocrit = 0.40,
    includePhaseSeparationEffect =  True,
  ),
)

mesentry_subset_not_so_bad = dict(
  adaption = dict(
    k_c = 1.15,
    k_m = 1.40,
    k_s = 1.1,
    Q_refdot = 30,
    S_0 = 15,
    cond_length = 1500.,
    max_nun_iterations = 0,
    qdev = 1.1,
    starting_radii = 0.,
    delta_t = .01,
    ),
  calcflow = dict(
    viscosityPlasma = 1.2e-6, #commented means using default for rats
    rheology = 'RheologyForRats',
    inletHematocrit = 0.40,
    includePhaseSeparationEffect =  True,
  ),
)

mesentry_subset_vary = dict(
  adaption = dict(
    k_c = 2.9,
    k_m = 0.83,
    k_s = 1.7,
    Q_refdot = 40,
    S_0 = 20,
    cond_length = 1800.,
    max_nun_iterations = 4000,
    qdev = 1.1,
    starting_radii = 0.,
    delta_t = .001,
    ),
  calcflow = dict(
    viscosityPlasma = 1.2e-6, #commented means using default for rats
    rheology = 'RheologyForRats',
    inletHematocrit = 0.40,
    includePhaseSeparationEffect =  True,
  ),
)
  
typeA = dict(
  adaption = dict(
    k_c = 1.5,
    k_m = 1.15,
    k_s = 1.7,
    Q_refdot = 45,
    S_0 = 25,
    cond_length = 9000.,
    #if this is 0 we iterate until qdev is reached
    max_nun_iterations = 1000,
    qdev = 1.1,
    #if starting_radii is 0. we use the values given in
    #the input file
    starting_radii = 0.,
    delta_t = 0.0001,
    ),
  calcflow = dict(
    viscosityPlasma = 1.2e-6, #commented means using default for rats
    rheology = 'RheologyForRats',
    inletHematocrit = 0.40,
    includePhaseSeparationEffect =  False,
  ),
)
  
typeE = dict(
  adaption = dict(
    k_c = 2.91,
    k_m = .8,
    k_s = 1.9,
    Q_refdot = 45,
    S_0 = 20,
    cond_length = 2000.,
    #if this is 0 we iterate until qdev is reached
    max_nun_iterations = 4000,
    qdev = .01,
    #if starting_radii is 0. we use the values given in
    #the input file
    starting_radii = 0.,
    delta_t = 0.001,
    ),
  calcflow = dict(
    viscosityPlasma = 1.2e-6, #commented means using default for rats
    rheology = 'RheologySecomb2005',
    inletHematocrit = 0.40,
    includePhaseSeparationEffect =  True,
  ),
)
typeE_bak = dict(
  adaption = dict(
    k_c = 1.35,
    k_m = 2.,
    k_s = 1.85,
    Q_refdot = 30,
    S_0 = 20,
    cond_length = 8000.,
    #if this is 0 we iterate until qdev is reached
    max_nun_iterations = 500,
    qdev = .1,
    #if starting_radii is 0. we use the values given in
    #the input file
    starting_radii = 0.,
    delta_t = 0.05,
    ),
  calcflow = dict(
    viscosityPlasma = 1.2e-6, #commented means using default for rats
    rheology = 'RheologyForRats',
    inletHematocrit = 0.40,
    includePhaseSeparationEffect =  True,
  ),
)
typeE5 = dict(
  adaption = dict(
    k_c = 1.85,
    k_m = 0.90,
    k_s = 1.3,
    Q_refdot = 35,
    S_0 = 15,
    cond_length = 7000.,
    #if this is 0 we iterate until qdev is reached
    max_nun_iterations = 1000,
    qdev = 1.1,
    #if starting_radii is 0. we use the values given in
    #the input file
    starting_radii = 20.,
    delta_t = 0.1,
    ),
  calcflow = dict(
    viscosityPlasma = 1.2e-6, #commented means using default for rats
    rheology = 'RheologySecomb2005',
    inletHematocrit = 0.40,
    includePhaseSeparationEffect =  True,
  ),
)

typeE_break = dict(
  num_threads = 6,
  adaption = dict(
    k_c = 2.85,
    k_m = 1.5,
    k_s = 1.85,
    Q_refdot = 22.,
    S_0 = 20,
    cond_length = 1500.,
    #if this is 0 we iterate until qdev is reached
    max_nun_iterations = 10000,
    qdev = 1.0,
    #if starting_radii is 0. we use the values given in
    #the input file
    starting_radii = 0.,
    delta_t = 0.01,
    ),
  calcflow = dict(
    viscosityPlasma = 1.2e-6, #commented means using default for rats
    rheology = 'RheologyForRats',
    inletHematocrit = 0.40,
    includePhaseSeparationEffect =  True,
  ),
)
typeE_break_longer_cond = dict(
  num_threads = 6,
  adaption = dict(
    k_c = 2.85,
    k_m = 1.5,
    k_s = 1.85,
    Q_refdot = 22.,
    S_0 = 20,
    cond_length = 3000.,
    #if this is 0 we iterate until qdev is reached
    max_nun_iterations = 10000,
    qdev = 1.0,
    #if starting_radii is 0. we use the values given in
    #the input file
    starting_radii = 0.,
    delta_t = 0.01,
    ),
  calcflow = dict(
    viscosityPlasma = 1.2e-6, #commented means using default for rats
    rheology = 'RheologyForRats',
    inletHematocrit = 0.40,
    includePhaseSeparationEffect =  True,
  ),
)
typeE_break_longer_cond_less_km = dict(
  num_threads = 6,
  adaption = dict(
    k_c = 2.85,
    k_m = 0.7,
    k_s = 1.85,
    Q_refdot = 42.,
    S_0 = 20,
    cond_length = 3000.,
    #if this is 0 we iterate until qdev is reached
    max_nun_iterations = 3000,
    qdev = 1.0,
    #if starting_radii is 0. we use the values given in
    #the input file
    starting_radii = 0.,
    delta_t = 0.1,
    ),
  calcflow = dict(
    viscosityPlasma = 1.2e-6, #commented means using default for rats
    rheology = 'RheologySecomb2005',
    inletHematocrit = 0.40,
    includePhaseSeparationEffect =  True,
  ),
)
typeE_break_longer_cond_human_rheo = dict(
  num_threads = 8,
  adaption = dict(
    k_c = 2.85,
    k_m = 1.5,
    k_s = 1.85,
    Q_refdot = 22.,
    S_0 = 20,
    cond_length = 3000.,
    #if this is 0 we iterate until qdev is reached
    max_nun_iterations = 3000,
    qdev = 1.0,
    #if starting_radii is 0. we use the values given in
    #the input file
    starting_radii = 0.,
    delta_t = 0.01,
    ),
  calcflow = dict(
    viscosityPlasma = 1.2e-6, #commented means using default for rats
    rheology = 'RheologySecomb2005',
    inletHematocrit = 0.40,
    includePhaseSeparationEffect =  True,
  ),
)
typeE_break_longer_cond_less_Q_ref = dict(
  num_threads = 8,
  adaption = dict(
    k_c = 2.85,
    k_m = 1.5,
    k_s = 1.85,
    Q_refdot = 10.,
    S_0 = 20,
    cond_length = 3000.,
    #if this is 0 we iterate until qdev is reached
    max_nun_iterations = 2500,
    qdev = 1.0,
    #if starting_radii is 0. we use the values given in
    #the input file
    starting_radii = 0.,
    delta_t = 0.01,
    ),
  calcflow = dict(
    viscosityPlasma = 1.2e-6, #commented means using default for rats
    rheology = 'RheologyForRats',
    inletHematocrit = 0.40,
    includePhaseSeparationEffect =  True,
  ),
)
typeE_break_longer_cond_higher_Q_ref = dict(
  num_threads = 8,
  adaption = dict(
    k_c = 2.85,
    k_m = 1.5,
    k_s = 1.85,
    Q_refdot = 30.,
    S_0 = 20,
    cond_length = 3000.,
    #if this is 0 we iterate until qdev is reached
    max_nun_iterations = 2500,
    qdev = 1.0,
    #if starting_radii is 0. we use the values given in
    #the input file
    starting_radii = 0.,
    delta_t = 0.01,
    ),
  calcflow = dict(
    viscosityPlasma = 1.2e-6, #commented means using default for rats
    rheology = 'RheologyForRats',
    inletHematocrit = 0.40,
    includePhaseSeparationEffect =  True,
  ),
)


dep_test = dict(
  num_threads = 6,
  adaption = dict(
    k_c = 2.85,
    k_m = 1.5,
    k_s = 2.85,
    Q_refdot = 22.,
    S_0 = 20,
    cond_length = 1500.,
    #if this is 0 we iterate until qdev is reached
    max_nun_iterations = 10000,
    qdev = 1.0,
    #if starting_radii is 0. we use the values given in
    #the input file
    starting_radii = 0.,
    delta_t = 0.01,
    ),
  calcflow = dict(
    viscosityPlasma = 1.2e-6, #commented means using default for rats
    rheology = 'RheologySecomb2005',
    inletHematocrit = 0.40,
    includePhaseSeparationEffect =  True,
  ),
)

mesentry_paper = dict(
  adaption = dict(
    k_c = 2.74,
    k_m = 0.83,
    k_s = 1.59,
    Q_refdot = 60,
    S_0 = 20,
    cond_length = 1500,
    max_nun_iterations = 0,
    qdev = 0.0001,
    starting_radii = 0.,
    delta_t = .01,
    ),
  calcflow = dict(
    viscosityPlasma = 1.2e-6, #commented means using default for rats
    rheology = 'RheologyForRats',
    inletHematocrit = 0.40,
    includePhaseSeparationEffect =  True,
  ),
)
mesentry_paper_original = dict(
  adaption = dict(
    k_c = 2.74,
    k_m = 0.83,
    k_s = 1.79,
    Q_refdot = 40,
    S_0 = 20,
    cond_length = 1500,
    max_nun_iterations = 1000,
    qdev = 0.001,
    starting_radii = 0.,
    delta_t = .01,
    ),
  calcflow = dict(
    viscosityPlasma = 1.2e-6, #commented means using default for rats
    rheology = 'RheologySecomb2005',
    inletHematocrit = 0.40,
    includePhaseSeparationEffect =  True,
  ),
)
mesentry_play = dict(
  adaption = dict(
    k_c = 2.74,
    k_m = 0.83,
    k_s = 1.79,
    Q_refdot = 40,
    S_0 = 20,
    cond_length = 1500,
    max_nun_iterations = 2000,
    qdev = 0.1,
    starting_radii = 0.,
    delta_t = .01,
    ),
  calcflow = dict(
    viscosityPlasma = 1.2e-6, #commented means using default for rats
    rheology = 'RheologyForRats',
    inletHematocrit = 0.40,
    includePhaseSeparationEffect =  True,
  ),
)
RK4 = dict(
  num_threads = 4,
  adaption = dict(
    k_c = 2.85,
    k_m = 1.5,
    k_s = 1.85,
    Q_refdot = 22.,
    S_0 = 20,
    cond_length = 1500.,
    #if this is 0 we iterate until qdev is reached
    max_nun_iterations = 2000,
    qdev = .00001,
    #if starting_radii is 0. we use the values given in
    #the input file
    starting_radii = 0.,
    delta_t = 0.1,
    ),
  calcflow = dict(
    viscosityPlasma = 1.2e-6, #commented means using default for rats
    rheology = 'RheologySecomb2005',
    inletHematocrit = 0.40,
    includePhaseSeparationEffect =  True,
  ),
)
RK4_b = dict(
  num_threads = 4,
  adaption = dict(
    k_c = 2.85,
    k_m = 1.5,
    k_s = 1.85,
    Q_refdot = 22.,
    S_0 = 20,
    cond_length = 1500.,
    #if this is 0 we iterate until qdev is reached
    max_nun_iterations = 2000,
    qdev = .00001,
    #if starting_radii is 0. we use the values given in
    #the input file
    starting_radii = 0.,
    delta_t = 0.01,
    ),
  calcflow = dict(
    viscosityPlasma = 1.2e-6, #commented means using default for rats
    rheology = 'RheologySecomb2005',
    inletHematocrit = 0.40,
    includePhaseSeparationEffect =  True,
  ),
)
QDEV_watch = dict(
  num_threads = 4,
  adaption = dict(
    k_c = 2.85,
    k_m = 1.5,
    k_s = 1.85,
    Q_refdot = 22.,
    S_0 = 20,
    cond_length = 1500.,
    #if this is 0 we iterate until qdev is reached
    max_nun_iterations = 3000,
    qdev = .00001,
    #if starting_radii is 0. we use the values given in
    #the input file
    starting_radii = 0.,
    delta_t = 0.01,
    ),
  calcflow = dict(
    viscosityPlasma = 1.2e-6, #commented means using default for rats
    rheology = 'RheologySecomb2005',
    inletHematocrit = 0.40,
    includePhaseSeparationEffect =  False,
  ),
)
nQDEV = dict(
  num_threads = 6,
  adaption = dict(
    k_c = 2.85,
    k_m = 1.5,
    k_s = 1.85,
    Q_refdot = 22.,
    S_0 = 20,
    cond_length = 1500.,
    #if this is 0 we iterate until qdev is reached
    max_nun_iterations = 3000,
    qdev = .00001,
    #if starting_radii is 0. we use the values given in
    #the input file
    starting_radii = 0.,
    delta_t = 0.01,
    ),
  calcflow = dict(
    viscosityPlasma = 1.2e-6, #commented means using default for rats
    rheology = 'RheologySecomb2005',
    inletHematocrit = 0.40,
    includePhaseSeparationEffect =  True,
  ),
)
appReal1 = dict(
  num_threads = 6,
  adaption = dict(
    k_c = 2.85,
    k_m = 1.5,
    k_s = 2.85,
    Q_refdot = 42.,
    S_0 = 20,
    cond_length = 3500.,
    #if this is 0 we iterate until qdev is reached
    max_nun_iterations = 3000,
    qdev = .00001,
    #if starting_radii is 0. we use the values given in
    #the input file
    starting_radii = 0.,
    delta_t = 0.1,
    ),
  calcflow = dict(
    viscosityPlasma = 1.2e-6, #commented means using default for rats
    rheology = 'RheologySecomb2005',
    inletHematocrit = 0.40,
    includePhaseSeparationEffect =  False,
  ),
)
appReal2 = dict(
  num_threads = 6,
  adaption = dict(
    k_c = 2.85,
    k_m = 1.5,
    k_s = 1.85,
    Q_refdot = 22.,
    S_0 = 10,
    cond_length = 3500.,
    #if this is 0 we iterate until qdev is reached
    max_nun_iterations = 3000,
    qdev = .00001,
    #if starting_radii is 0. we use the values given in
    #the input file
    starting_radii = 0.,
    delta_t = 0.1,
    ),
  calcflow = dict(
    viscosityPlasma = 1.2e-6, #commented means using default for rats
    rheology = 'RheologySecomb2005',
    inletHematocrit = 0.40,
    includePhaseSeparationEffect =  False,
  ),
)
appReal3 = dict(
  num_threads = 6,
  adaption = dict(
    k_c = 1.85,
    k_m = 1.5,
    k_s = 1.85,
    Q_refdot = 22.,
    S_0 = 20,
    cond_length = 1500.,
    #if this is 0 we iterate until qdev is reached
    max_nun_iterations = 3500,
    qdev = .00001,
    #if starting_radii is 0. we use the values given in
    #the input file
    starting_radii = 0.,
    delta_t = 0.1,
    ),
  calcflow = dict(
    viscosityPlasma = 1.2e-6, #commented means using default for rats
    rheology = 'RheologySecomb2005',
    inletHematocrit = 0.40,
    includePhaseSeparationEffect =  False,
  ),
)
appReal4 = dict(
  num_threads = 6,
  adaption = dict(
    k_c = 2.85,
    k_m = 1.5,
    k_s = 1.85,
    Q_refdot = 42.,
    S_0 = 20,
    cond_length = 3500.,
    #if this is 0 we iterate until qdev is reached
    max_nun_iterations = 3000,
    qdev = .00001,
    #if starting_radii is 0. we use the values given in
    #the input file
    starting_radii = 0.,
    delta_t = 0.01,
    ),
  calcflow = dict(
    viscosityPlasma = 1.2e-6, #commented means using default for rats
    rheology = 'RheologySecomb2005',
    inletHematocrit = 0.40,
    includePhaseSeparationEffect =  False,
  ),
)
appReal5 = dict(
  num_threads = 6,
  adaption = dict(
    k_c = 2.85,
    k_m = 1.5,
    k_s = 1.95,
    Q_refdot = 40.,
    S_0 = 20,
    cond_length = 3500.,
    #if this is 0 we iterate until qdev is reached
    max_nun_iterations = 3000,
    qdev = .00001,
    #if starting_radii is 0. we use the values given in
    #the input file
    starting_radii = 0.,
    delta_t = 0.01,
    ),
  calcflow = dict(
    viscosityPlasma = 1.2e-6, #commented means using default for rats
    rheology = 'RheologySecomb2005',
    inletHematocrit = 0.40,
    includePhaseSeparationEffect =  False,
  ),
)
appReal6 = dict(
  num_threads = 6,
  adaption = dict(
    k_c = 2.85,
    k_m = 1.5,
    k_s = 1.95,
    Q_refdot = 10.,
    S_0 = 20,
    cond_length = 3500.,
    #if this is 0 we iterate until qdev is reached
    max_nun_iterations = 3000,
    qdev = .001,
    #if starting_radii is 0. we use the values given in
    #the input file
    starting_radii = 0.,
    delta_t = 0.01,
    ),
  calcflow = dict(
    viscosityPlasma = 1.2e-6, #commented means using default for rats
    rheology = 'RheologySecomb2005',
    inletHematocrit = 0.40,
    includePhaseSeparationEffect =  False,
  ),
)
secomb_mesentry=deepcopy(default)