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
  BOTH_PRESSURE = 43,
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
    max_nun_iterations = 200,
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
    includePhaseSeparationEffect = 0,
  ),
)

default_phase=deepcopy(default)
default_phase['calcflow'].update(
  includePhaseSeparationEffect = 1,
  rheology = 'RheologySecomb2005',
)
apj = dict(
  num_threads = 1,
  name = 'apj',
  #this is converging, note k_s changed from 1.79 to 1.99
  adaption = dict(
    k_c = 2.74,
    k_m = 0.83,
    k_s = 1.99,
    Q_refdot = 40,
    S_0 = 20,
    cond_length = 1500.,
    #if this is 0 we iterate until qdev is reached
    max_nun_iterations = 250,
    qdev = .001,
    #if starting_radii is 0. we use the values given in
    #the input file
    starting_radii = 0.,
    # 0.0 means the time step is dynamically chosen!
    delta_t = .1, 
    boundary_Condition_handling = boundary_Condition_handling_map['KEEP'],
    ),
  calcflow = dict(
    viscosityPlasma = 1.2e-6, #commented means using default for rats
    rheology = 'RheologyForRats',
    inletHematocrit = 0.42,
    includePhaseSeparationEffect = 1,
  ),
)
apj_initial = deepcopy(apj)
apj_initial['adaption'].update(
    qdev = 10000., # other adaption will break since it is not convergent in this case
    max_mun_iterations = 2, # we need at least 2 iterations to get some data at all
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
    includePhaseSeparationEffect = 0,
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
    includePhaseSeparationEffect = 0,
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
    includePhaseSeparationEffect = 1,
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
    includePhaseSeparationEffect = 1,
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
    includePhaseSeparationEffect = 1,
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
    includePhaseSeparationEffect = 1,
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
    includePhaseSeparationEffect = 1,
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
    includePhaseSeparationEffect = 1,
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
    includePhaseSeparationEffect = 1,
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
    includePhaseSeparationEffect = 1,
  ),
)

mesentry_paper = dict(
  name = 'mesentry_paper',
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
    includePhaseSeparationEffect = 1,
  ),
)

'''taken from 
1.Pries, A. R., Secomb, T. W. & Gaehtgens, P. 
Structural adaptation and stability of microvascular networks: theory and simulations. American Journal of Physiology - Heart and Circulatory Physiology 275, H349–H360 (1998).
'''
mesentry = dict(
  num_threads = 4,
  name = 'mesentry',
  adaption = dict(
    k_c = 2.74,
    k_m = 0.83,
    k_s = 1.79,
    Q_refdot = 40.,
    S_0 = 20,
    cond_length = 1500.,
    #if this is 0 we iterate until qdev is reached
    max_nun_iterations = 300,
    #qdev = 0, means local conditions!
    qdev = .001,
    #if starting_radii is 0. we use the values given in
    #the input file
    starting_radii = 0.,
    delta_t = 0.1,
    addAverageConductivity = False,
    radMin_for_kill = 1.5,
    boundary_Condition_handling = boundary_Condition_handling_map['KEEP'],
    ),
  calcflow = dict(
    viscosityPlasma = 1.2e-6, #commented means using default for rats
    rheology = 'RheologyForRats',
    inletHematocrit = 0.40,
    includePhaseSeparationEffect = 0,
  ),
)
mesentry_new_h = deepcopy(mesentry)
mesentry_new_h['name'] = 'mesentry_new_h'
mesentry_new_h['calcflow'].update(
  includePhaseSeparationEffect = 1,
)
mesentry_new_h_initial = deepcopy(mesentry_new_h)
mesentry_new_h_initial['name'] = 'mesentry_new_h_initial'
mesentry_new_h_initial['adaption'].update(
    delta_t = 0.01,
    qdev = 179769313486231570814527423731704356798070567525844996598917476803157260780028538760589558632766878171540458953514382464234321326889464182768467546703537516986049910576551282076245490090389328944075868508455133942304583236903222948165808559332123348274797826204144723168738177180919299881250404026184124858368.000000,
    max_nun_iterations = 2,
    )
deap_test = dict(
  num_threads = 1,
  adaption = dict(
    k_c = 2.74,
    k_m = 0.83,
    k_s = 1.79,
    Q_refdot = 40.,
    S_0 = 20,
    cond_length = 1500.,
    #if this is 0 we iterate until qdev is reached
    max_nun_iterations = 500,
    #qdev = 0, means local conditions!
    qdev = 0.10,
    #if starting_radii is 0. we use the values given in
    #the input file
    starting_radii = 0.,
    delta_t = 0.05,
    addAverageConductivity = False,
    radMin_for_kill = 3.5,
    boundary_Condition_handling = boundary_Condition_handling_map['KEEP'],
    a_pressure = 2.,
    a_flow = 20000.,
    pop = 3,
    individuals = 9,
    opt_iter = 10,
    ),
  calcflow = dict(
    viscosityPlasma = 1.2e-6, #commented means using default for rats
    rheology = 'RheologySecomb2005',
    inletHematocrit = 0.40,
    includePhaseSeparationEffect = 0,
  ),
  optimization = dict(
      desired_cap_flow = 20000,      
      )
)
deap_test_1_value= deepcopy(deap_test)
deap_test_1_value['adaption'].update(
    a_pressure = 1.8,
    a_flow = -800000,
    boundary_Condition_handling = boundary_Condition_handling_map['VALUE'],
    )
deap_test_2 = dict(
  num_threads = 1,
  adaption = dict(
    k_c = 2.74,
    k_m = 0.83,
    k_s = 1.79,
    Q_refdot = 40.,
    S_0 = 20,
    cond_length = 1500.,
    #if this is 0 we iterate until qdev is reached
    max_nun_iterations = 500,
    #qdev = 0, means local conditions!
    qdev = 0.10,
    #if starting_radii is 0. we use the values given in
    #the input file
    starting_radii = 0.,
    delta_t = 0.05,
    addAverageConductivity = False,
    radMin_for_kill = 3.5,
    boundary_Condition_handling = boundary_Condition_handling_map['KEEP'],
    a_pressure = 2.,
    a_flow = 20000.,
    pop = 3,
    individuals = 9,
    opt_iter = 10,
    ),
  calcflow = dict(
    viscosityPlasma = 1.2e-6, #commented means using default for rats
    rheology = 'RheologySecomb2005',
    inletHematocrit = 0.40,
    includePhaseSeparationEffect = 1,
  ),
  optimization = dict(
      desired_cap_flow = 20000,      
      )
)
deap_test_3 = deepcopy(deap_test_2)
deap_test_3['adaption'].update(
    starting_radii=10.,
    )

deap_test_4 = deepcopy(deap_test_3)
deap_test_4['adaption'].update(
    starting_radii = 0.,
    boundary_Condition_handling = boundary_Condition_handling_map['VALUE'],
    a_pressure = 1.8,
    a_flow = -800000,
    )


# velocity in capillary about 1mm/s = 1000mu/s
# typical diameter 5mu e.g r=2.5mu
# 3.1415*2.5*2.5*1000 about 2e4 mu^3/s
def _value1():
  #import random
  pressures = [1.8,2.3, 3.]
  flows = [-1e5, -1e6, -1e7, ]
  import itertools
  combined = list(itertools.product(pressures,flows))
  def mk1(cnt):
    p = deepcopy(deap_test_4)
    p['name'] = 'value_%i' % cnt
    pressure = combined[cnt][0]
    flow = combined[cnt][1]
    p['adaption'].update(
      a_pressure = pressure,
      a_flow = flow,
    )
    return p
  return list(mk1(i) for i in xrange(len(combined)))

value_list1 = _value1()

def _value2():
  #import random
  pressures = [1.5, 2, 3,]
  flows = [-1e6, -1e8, -1e9, ]
  import itertools
  combined = list(itertools.product(pressures,flows))
  def mk1(cnt):
    p = deepcopy(deap_test_4)
    p['name'] = 'value_list2_%i' % cnt
    pressure = combined[cnt][0]
    flow = combined[cnt][1]
    p['adaption'].update(
      a_pressure = pressure,
      a_flow = flow,
      max_nun_iterations = 150,
    )
    return p
  return list(mk1(i) for i in xrange(len(combined)))
value_list2 = _value2()

def _value3():
  #import random
  pressures = [2.3, 2.8, 3.3,]
  flows = [-4.8e5, -4.8e6, -4.8e7,]
  import itertools
  combined = list(itertools.product(pressures,flows))
  def mk1(cnt):
    p = deepcopy(deap_test_4)
    p['name'] = 'value_list3_%i' % cnt
    pressure = combined[cnt][0]
    flow = combined[cnt][1]
    p['adaption'].update(
      a_pressure = pressure,
      a_flow = flow,
      max_nun_iterations = 150,
    )
    return p
  return list(mk1(i) for i in xrange(len(combined)))
value_list3 = _value3()

def _value4():
  #import random
  pressures = [2.3, 2.8, 3.3,]
  flows = [-4.8e5, -4.8e6, -4.8e7,]
  import itertools
  combined = list(itertools.product(pressures,flows))
  def mk1(cnt):
    p = deepcopy(deap_test_4)
    p['name'] = 'value_list4_%i' % cnt
    pressure = combined[cnt][0]
    flow = combined[cnt][1]
    p['adaption'].update(
      a_pressure = pressure,
      a_flow = flow,
      max_nun_iterations = 150,
    )
    return p
  return list(mk1(i) for i in xrange(len(combined)))
value_list4 = _value4()

def _value5():
  #import random
  pressures = [2.3, 2.8, 3.3,]
  flows = [-4.8e5, -4.8e6, -4.8e7,]
  import itertools
  combined = list(itertools.product(pressures,flows))
  def mk1(cnt):
    p = deepcopy(deap_test_4)
    p['name'] = 'value_list5_%i' % cnt
    pressure = combined[cnt][0]
    flow = combined[cnt][1]
    p['adaption'].update(
      a_pressure = pressure,
      a_flow = flow,
      max_nun_iterations = 150,
      Q_refdot = 10.,
    )
    return p
  return list(mk1(i) for i in xrange(len(combined)))
value_list5 = _value5()

def _value6():
  #import random
  pressures = [2.3, 2.8, 3.3,]
  flows = [-4.8e5, -4.8e6, -4.8e7,]
  import itertools
  combined = list(itertools.product(pressures,flows))
  def mk1(cnt):
    p = deepcopy(deap_test_4)
    p['name'] = 'value_list6_%i' % cnt
    pressure = combined[cnt][0]
    flow = combined[cnt][1]
    p['adaption'].update(
      a_pressure = pressure,
      a_flow = flow,
      max_nun_iterations = 150,
      #Q_refdot = 10.,
      cond_length = 2500.,
    )
    return p
  return list(mk1(i) for i in xrange(len(combined)))
value_list6 = _value6()

def _value7():
  #import random
  pressures = [2.3, 2.8, 3.3,]
  flows = [-4.8e5, -4.8e6, -4.8e7,]
  import itertools
  combined = list(itertools.product(pressures,flows))
  def mk1(cnt):
    p = deepcopy(deap_test_4)
    p['name'] = 'value_list7_%i' % cnt
    pressure = combined[cnt][0]
    flow = combined[cnt][1]
    p['adaption'].update(
      a_pressure = pressure,
      a_flow = flow,
      max_nun_iterations = 150,
      Q_refdot = 10.,
      cond_length = 2500.,
    )
    return p
  return list(mk1(i) for i in xrange(len(combined)))
value_list7 = _value7()

def _value8():
  #import random
  pressures = [3., 4, 5,]
  flows = [-4.8e6, -4.8e7, -4.8e8,]
  import itertools
  combined = list(itertools.product(pressures,flows))
  def mk1(cnt):
    p = deepcopy(deap_test_4)
    p['name'] = 'value_list8_%i' % cnt
    pressure = combined[cnt][0]
    flow = combined[cnt][1]
    p['adaption'].update(
      a_pressure = pressure,
      a_flow = flow,
      max_nun_iterations = 150,
      Q_refdot = 10.,
      cond_length = 2500.,
    )
    return p
  return list(mk1(i) for i in xrange(len(combined)))
value_list8 = _value8()

def _value9():
  #import random
  pressures = [3., 4, 5,]
  flows = [-4.8e6, -4.8e7, -4.8e8,]
  import itertools
  combined = list(itertools.product(pressures,flows))
  def mk1(cnt):
    p = deepcopy(deap_test_4)
    p['name'] = 'value_list9_%i' % cnt
    pressure = combined[cnt][0]
    flow = combined[cnt][1]
    p['adaption'].update(
      a_pressure = pressure,
      a_flow = flow,
      max_nun_iterations = 150,
      #Q_refdot = 10.,
      #cond_length = 2500.,
    )
    return p
  return list(mk1(i) for i in xrange(len(combined)))
value_list9 = _value9()

def _value10():
  #import random
  pressures = [6., 7, 8,]
  flows = [-4.8e8, -4.8e9, -4.8e10,]
  import itertools
  combined = list(itertools.product(pressures,flows))
  def mk1(cnt):
    p = deepcopy(deap_test_4)
    p['name'] = 'value_list10_%i' % cnt
    pressure = combined[cnt][0]
    flow = combined[cnt][1]
    p['adaption'].update(
      a_pressure = pressure,
      a_flow = flow,
      max_nun_iterations = 150,
      Q_refdot = 10.,
      cond_length = 2500.,
    )
    return p
  return list(mk1(i) for i in xrange(len(combined)))
value_list10 = _value10()

def _value11():
  #import random
  pressures = np.arange(3.5,4.5,0.1)
  flows = np.arange(-9.0,-1, 0.5) * 1e7
  import itertools
  combined = list(itertools.product(pressures,flows))
  def mk1(cnt):
    p = deepcopy(deap_test_4)
    p['name'] = 'value_list11_%i' % cnt
    pressure = combined[cnt][0]
    flow = combined[cnt][1]
    p['adaption'].update(
      a_pressure = pressure,
      a_flow = flow,
      max_nun_iterations = 150,
      #Q_refdot = 10.,
      #cond_length = 2500.,
    )
    return p
  return list(mk1(i) for i in xrange(len(combined)))
value_list11 = _value11()
def _value12():
  #import random
  pressures = np.arange(3.5,4.5,0.1)
  flows = np.arange(-9.0,-1, 0.5) * 1e7
  import itertools
  combined = list(itertools.product(pressures,flows))
  def mk1(cnt):
    p = deepcopy(deap_test_4)
    p['name'] = 'value_list12_%i' % cnt
    pressure = combined[cnt][0]
    flow = combined[cnt][1]
    p['adaption'].update(
      a_pressure = pressure,
      a_flow = flow,
      max_nun_iterations = 150,
      starting_radii = 15.,
      #Q_refdot = 10.,
      #cond_length = 2500.,
    )
    return p
  return list(mk1(i) for i in xrange(len(combined)))
value_list12 = _value12()
def _value13():
  #import random
  pressures = np.arange(1.8,4.5,0.1)
  flows = np.arange(-9.0,-1, 0.5) * 1e7
  flows2 = np.arange(-9.0,-1, 0.5) * 1e6
  flows = np.append(flows,flows2)
  import itertools
  combined = list(itertools.product(pressures,flows))
  def mk1(cnt):
    p = deepcopy(deap_test_4)
    p['name'] = 'value_list13_%i' % cnt
    pressure = combined[cnt][0]
    flow = combined[cnt][1]
    p['adaption'].update(
      a_pressure = pressure,
      a_flow = flow,
      max_nun_iterations = 150,
      #Q_refdot = 10.,
      #cond_length = 2500.,
    )
    return p
  return list(mk1(i) for i in xrange(len(combined)))
value_list13 = _value13()
def _value14():
  #import random
  pressures = np.arange(.5,9.5,0.5)
  flows = np.arange(-10.,-0.5, 0.5) * 1e7
  import itertools
  combined = list(itertools.product(pressures,flows))
  def mk1(cnt):
    p = deepcopy(deap_test_4)
    p['name'] = 'value_list14_%i' % cnt
    pressure = combined[cnt][0]
    flow = combined[cnt][1]
    p['adaption'].update(
      a_pressure = pressure,
      a_flow = flow,
      max_nun_iterations = 150,
      #Q_refdot = 10.,
      #cond_length = 2500.,
    )
    return p
  return list(mk1(i) for i in xrange(len(combined)))
value_list14 = _value14()
''' the a_pressure variable is the arterial pressure
    the a_flow variable is venous pressure
'''
default_pso = deepcopy(value_list11[142])
#this is from the optimization
default_pso['adaption'].update(
    k_c=3.27,
    k_m=1.31,
    k_s=1.85,
    )
def _value_pressure_1():
  #import random
  pressures1 = [6., 7., 8.,] # arteries
  pressures2 = [3., 4., 5.]  # veins
  import itertools
  combined = list(itertools.product(pressures1,pressures2))
  def mk1(cnt):
    p = deepcopy(deap_test_4)
    p['boundary_Condition_handling'] = boundary_Condition_handling_map['BOTH_PRESSURE'],
    p['name'] = 'value_list_pressure_1_%i' % cnt
    pressure = combined[cnt][0]
    flow = combined[cnt][1]
    p['adaption'].update(
      a_pressure = pressure,
      a_flow = flow,
      max_nun_iterations = 150,
      Q_refdot = 10.,
      cond_length = 2500.,
    )
    return p
  return list(mk1(i) for i in xrange(len(combined)))
value_list_pressure_1 = _value_pressure_1()
def _value_pressure_2():
  #import random
  pressures1 = [3., 3.5, 4.,] # arteries
  pressures2 = [1.5, 2., 2.5]  # veins
  import itertools
  combined = list(itertools.product(pressures1,pressures2))
  def mk1(cnt):
    p = deepcopy(deap_test_4)
    p['boundary_Condition_handling'] = boundary_Condition_handling_map['BOTH_PRESSURE'],
    p['name'] = 'value_list_pressure_2_%i' % cnt
    pressure = combined[cnt][0]
    flow = combined[cnt][1]
    p['adaption'].update(
      a_pressure = pressure,
      a_flow = flow,
      max_nun_iterations = 150,
      Q_refdot = 40.,
      cond_length = 1700.,
    )
    return p
  return list(mk1(i) for i in xrange(len(combined)))
value_list_pressure_2 = _value_pressure_2()

deap_id5 = deepcopy(deap_test_4)
deap_id5.update(
   boundary_Condition_handling = boundary_Condition_handling_map['BOTH_PRESSURE'],
   name = 'deap_id5',
   )
deap_id5['adaption'].update(
    a_pressure = 2.81,
    a_flow = 6.87,
    max_nun_iterations = 150,
    Q_refdot = 40.,
    cond_length = 1700.,
    )
   
final_I = dict(
  num_threads = 1,
  name = 'final_I',
  adaption = dict(
    k_c = 2.74,
    k_m = 0.83,
    k_s = 1.79,
    Q_refdot = 40.,
    S_0 = 20,
    cond_length = 1500.,
    #if this is 0 we iterate until qdev is reached
    max_nun_iterations = 500,
    #qdev = 0, means local conditions!
    qdev = 0.10,
    #if starting_radii is 0. we use the values given in
    #the input file
    starting_radii = 0.,
    delta_t = 0.05,
    addAverageConductivity = False,
    radMin_for_kill = 3.5,
    boundary_Condition_handling = boundary_Condition_handling_map['KEEP'],
    a_pressure = 2.8,
    a_flow = 20000000.,
    pop = 3,
    individuals = 9,
    opt_iter = 10,
    ),
  calcflow = dict(
    viscosityPlasma = 1.2e-6, #commented means using default for rats
    rheology = 'RheologySecomb2005',
    inletHematocrit = 0.40,
    includePhaseSeparationEffect = 0,
  ),
  optimization = dict(
      desired_cap_flow = 20000,      
      )
)
final_II = deepcopy(final_I)
final_II['adaption']['k_m'] = 0.63
final_II['adaption']['k_c'] = 3.19
final_II['adaption']['k_s'] = 1.57
final_II['name'] = 'final_II'

final_III = deepcopy(final_II)
final_III['adaption']['k_m'] = 0.97
final_III['adaption']['k_c'] = 2.85
final_III['adaption']['k_s'] = 1.79
final_III['name'] = 'final_III'
final_III['adaption']['max_nun_iterations'] = 1000

#this is clearly in the successfull adapted range
final_IV = deepcopy(final_II)
final_IV['adaption']['k_m'] = 1.76
final_IV['adaption']['k_c'] = 2.8
final_IV['adaption']['k_s'] = 2.33
final_IV['adaption']['max_nun_iterations'] = 1000
final_IV['name'] = 'final_IV'
final_IV['calcflow'].update(
    includePhaseSeparationEffect=True,
    )
final_V = deepcopy(final_I)
final_V['calcflow'].update(
    includePhaseSeparationEffect=True,
    )
final_apj = deepcopy(final_IV)
final_apj['adaption']['k_m'] = 0.83
final_apj['adaption']['k_c'] = 2.74
final_apj['adaption']['k_s'] = 1.79
if __name__ == '__main__':
#  index = 142
#  print_c = deepcopy(value_list14)
#  print('starting radii: %f ' % print_c[index]['adaption']['starting_radii'])
#  print('a_pressure: %f ' % print_c[index]['adaption']['a_pressure'])
#  print('a_flow: %0.1g ' % print_c[index]['adaption']['a_flow'])
#  print('a_flow: %f ' % print_c[index]['adaption']['a_flow'])
#  print('handling: %s ' % print_c[index]['adaption']['boundary_Condition_handling'])
#  print(len(print_c))
#  print(final_IV)
  
  print(value_list_pressure_2)
