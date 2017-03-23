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

boundary_Condition_handling_map = dict(
  KEEP = 0,
  VEIN_AS_FLOW_ARTERY_PRESSURE = 1,
  LARGE_2D = 2,
  LARGE_2D_2 = 3,
  LARGE_2D_3 = 4,
  LARGE_2D_like_paper = 5,
  VALUE = 42,
)

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

pyGmo = dict(
  num_threads = 1,
  adaption = dict(
    k_c = 1.5,
    k_m = 1.5,
    k_s = 1.,
    Q_refdot = 40,
    S_0 = 20,
    cond_length = 1500.,
    #if this is 0 we iterate until qdev is reached
    max_nun_iterations = 1000,
    qdev = 0.01,
    #if starting_radii is 0. we use the values given in
    #the input file
    starting_radii = 10.,
    # 0.0 means the time step is dynamically chosen!
    delta_t = 0.1, 
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
adaption_fix_r_10 = deepcopy(default_phase)
adaption_fix_r_10['adaption'].update(
  starting_radii = 10.,
)

adaption_fix_r_5 = deepcopy(default_phase)
adaption_fix_r_5['adaption'].update(
  starting_radii = 5.,
  delta_t = 0.01,
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

'''taken from 
1.Pries, A. R., Secomb, T. W. & Gaehtgens, P. 
Structural adaptation and stability of microvascular networks: theory and simulations. American Journal of Physiology - Heart and Circulatory Physiology 275, H349–H360 (1998).
'''
mesentry = dict(
  num_threads = 4,
  adaption = dict(
    k_c = 2.74,
    k_m = 0.83,
    k_s = 1.79,
    Q_refdot = 40.,
    S_0 = 20,
    cond_length = 1500.,
    #if this is 0 we iterate until qdev is reached
    max_nun_iterations = 1000,
    #qdev = 0, means local conditions!
    qdev = .0,
    #if starting_radii is 0. we use the values given in
    #the input file
    starting_radii = 0.,
    delta_t = 0.05,
    addAverageConductivity = False,
    radMin_for_kill = 1.5,
    boundary_Condition_handling = boundary_Condition_handling_map['KEEP'],
    ),
  calcflow = dict(
    viscosityPlasma = 1.2e-6, #commented means using default for rats
    rheology = 'RheologyForRats',
    inletHematocrit = 0.40,
    includePhaseSeparationEffect =  False,
  ),
)
mesentry_new_h = deepcopy(mesentry)
mesentry_new_h['calcflow'].update(
  includePhaseSeparationEffect =  True,
)
large_2d = dict(
  num_threads = 4,
  adaption = dict(
    k_c = 2.74,
    k_m = 0.83,
    k_s = 1.79,
    Q_refdot = 40.,
    S_0 = 20,
    cond_length = 1500.,
    #if this is 0 we iterate until qdev is reached
    max_nun_iterations = 1000,
    #qdev = 0, means local conditions!
    qdev = .0,
    #if starting_radii is 0. we use the values given in
    #the input file
    starting_radii = 10.,
    delta_t = 0.05,
    addAverageConductivity = False,
    radMin_for_kill = 1.5,
    boundary_Condition_handling = boundary_Condition_handling_map['LARGE_2D'],
    ),
  calcflow = dict(
    viscosityPlasma = 1.2e-6, #commented means using default for rats
    rheology = 'RheologyForRats',
    inletHematocrit = 0.40,
    includePhaseSeparationEffect =  True,
  ),
)
large_2d_2 = deepcopy(large_2d)
large_2d_2['adaption'].update(
  boundary_Condition_handling = boundary_Condition_handling_map['LARGE_2D_2'],
)
human_guess = dict(
  num_threads = 4,
  adaption = dict(
    k_c = 2.74,
    k_m = 0.83,
    k_s = 1.79,
    Q_refdot = 40.,
    S_0 = 20,
    cond_length = 1500.,
    #if this is 0 we iterate until qdev is reached
    max_nun_iterations = 1000,
    #qdev = 0, means local conditions!
    qdev = .0,
    #if starting_radii is 0. we use the values given in
    #the input file
    starting_radii = 20.,
    delta_t = 0.05,
    addAverageConductivity = False,
    radMin_for_kill = 3.5,
    boundary_Condition_handling = boundary_Condition_handling_map['LARGE_2D_2'],
    ),
  calcflow = dict(
    viscosityPlasma = 1.2e-6, #commented means using default for rats
    rheology = 'RheologySecomb2005',
    inletHematocrit = 0.40,
    includePhaseSeparationEffect =  True,
  ),
)
human_guess_2 = deepcopy(human_guess)
human_guess_2['adaption'].update(
  boundary_Condition_handling = boundary_Condition_handling_map['LARGE_2D_3'],
)
value_1 = deepcopy(human_guess)
value_1['adaption'].update(
  boundary_Condition_handling = boundary_Condition_handling_map['VALUE'],
  a_pressure = 2.,
  a_flow = 20000.,
  pop = 50,
  individuals = 300,
)
value_2 = deepcopy(human_guess)
value_2['adaption'].update(
  boundary_Condition_handling = boundary_Condition_handling_map['VALUE'],
  a_pressure = 2.,
  a_flow = 50000.,
  pop = 5,
  individuals = 50,
  opt_iter = 10,
  max_nun_iterations = 150,
  delta_t = 0.1,
)
pagmo_test = deepcopy(value_1)
pagmo_test['adaption'].update(
  boundary_Condition_handling = boundary_Condition_handling_map['VALUE'],
  a_pressure = 2.,
  a_flow = 20000.,
  pop = 3,
  individuals = 9,
  opt_iter = 10,
  max_nun_iterations = 10,
)
pagmo_test['calcflow'].update(
  includePhaseSeparationEffect =  False,
)

def _value():
  #import random
  pressures = [2.,2.5,3,3.5,4.]
  flows = [1e5,1e6,1e7,1e8]
  import itertools
  combined = list(itertools.product(pressures,flows))
  def mk1(cnt):
    p = deepcopy(value_1)
    p['name'] = 'value_%i' % cnt
    pressure = combined[cnt][0]
    flow = combined[cnt][1]
    p['adaption'].update(
      a_pressure = pressure,
      a_flow = flow,
    )
    return p
  return list(mk1(i) for i in xrange(len(combined)))

value_list = _value()

LARGE_2D_like_paper = dict(
  num_threads = 4,
  adaption = dict(
    k_c = 2.74,
    k_m = 0.83,
    k_s = 1.79,
    Q_refdot = 40.,
    S_0 = 20,
    cond_length = 1500.,
    #if this is 0 we iterate until qdev is reached
    max_nun_iterations = 1000,
    #qdev = 0, means local conditions!
    qdev = .0,
    #if starting_radii is 0. we use the values given in
    #the input file
    starting_radii = 10.,
    delta_t = 0.05,
    addAverageConductivity = False,
    radMin_for_kill = 1.5,
    boundary_Condition_handling = boundary_Condition_handling_map['LARGE_2D_like_paper'],
    ),
  calcflow = dict(
    viscosityPlasma = 1.2e-6, #commented means using default for rats
    rheology = 'RheologyForRats',
    inletHematocrit = 0.40,
    includePhaseSeparationEffect =  True,
  ),
)
if __name__ == '__main__':
  abc = _value()
  print abc