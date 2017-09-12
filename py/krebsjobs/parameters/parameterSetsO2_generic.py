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
import myutils
cluster_threads = myutils.cluster_threads
#c1_, c2_, c3_ = (0.0052*1.1, 16, 0.0078/math.exp(-4./16)*1.6) # hand adjusted parameters based on secomb
''' see S1 from 
    https://doi.org/10.1371/journal.pone.0161267.s001
'''
c1_, c2_, c3_ = (8, 4.7, 0) 

default_o2 = dict(
  num_threads = cluster_threads,
  po2init_r0 = 55., #  mmHg;  po2_init = min(po2init_cutoff, po2_init_r0 + v->r * po2_init_dr)
  po2init_dr = 1., #  mmHg / um, these values are from yaseen 2011 for rat brain
  po2init_cutoff = 100., # mmHg; maximal attained value
  solubility_plasma = 3.1e-5,
  c0 = 0.5,
  sat_curve_exponent = 2.7,
  sat_curve_p50 = 27.,
  D_plasma = 2750.,
  solubility_tissue = 2.8e-5,
  rd_norm = 150.,
  rd_tum  = 50.,
  rd_necro = 150.,
  max_iter = 50,
  convergence_tolerance = 0.03,
  axial_integration_step_factor = .1,
  debug_zero_o2field = False,
  grid_lattice_const = 40.,
  calcflow = dict(
    viscosityPlasma = 1.2e-6,
    rheology = 'RheologySecomb2005',
    inletHematocrit = 0.45,
    includePhaseSeparationEffect = 1,
  ),
  michaelis_menten_uptake = 1,
  mmcons_m0_norm = 6.2e-5, # = 3.7e-3 ml O2 / ml / min; # fairly low like estimates by rinneberg and beany
  mmcons_m0_tum = 6.2e-5 * 4,
  massTransferCoefficientModelNumber = 1,
  conductivity_coeff1 = c1_,
  conductivity_coeff2 = c2_,
  conductivity_coeff3 = c3_,
  detailedO2name = 'default_o2',
)
small_size = deepcopy(default_o2)
small_size.update(
    axial_integration_step_factor = 0.05,
    grid_lattice_const = 10,
    )
small_size2 = deepcopy(default_o2)
small_size2.update(
    axial_integration_step_factor = 0.1,
    grid_lattice_const = 30,
    po2init_r0 = 25., #  mmHg;  po2_init = min(po2init_cutoff, po2_init_r0 + v->r * po2_init_dr)
    po2init_dr = .5, #  mmHg / um, these values are from yaseen 2011 for rat brain
    po2init_cutoff = 75., # mmHg; maximal attained value
    convergence_tolerance = 1e-3,
)
no_mm = deepcopy(default_o2)
no_mm['michaelis_menten_uptake'] = 0
no_mm['D_plasma'] = 2000.

milotti_o2 = deepcopy(default_o2)
milotti_o2.update(
    detailedO2name = 'milotti_o2',
    axial_integration_step_factor = 0.25,
    grid_lattice_const = 50.,
    max_iter = 150,
    convergence_tolerance = 1e-3,
    loglevel=2,
    )
milotti_o2_simple = deepcopy(milotti_o2)
milotti_o2_simple['calcflow'].update(
    includePhaseSeparationEffect = 0,
    )
exp = deepcopy(default_o2)
exp.update(
    detailedO2name = 'exp',
    axial_integration_step_factor = 0.25,
    grid_lattice_const = 100.,
    max_iter = 150,
    po2init_r0 = 55., #  mmHg;  po2_init = min(po2init_cutoff, po2_init_r0 + v->r * po2_init_dr)
    po2init_dr = 1., #  mmHg / um, these values are from yaseen 2011 for rat brain
    po2init_cutoff = 200., # mmHg; maximal attained value
    )

lowo2 = dict(
  po2init_r0 = 55., #  mmHg;  po2_init = min(po2init_cutoff, po2_init_r0 + v->r * po2_init_dr)
  po2init_dr = 1., #  mmHg / um, these values are from yaseen 2011 for rat brain
  po2init_cutoff = 100., # mmHg; maximal attained value
  solubility_plasma = 3.1e-5,
  c0 = 0.5,
  S_n = 2.7,
  S_p50 = 27.,
  D_plasma = 2000.,
  solubility_tissue = 2.8e-5,
  #dM = 0.05,
  rd_norm = 200.,
  rd_tum  = 100.,
  rd_necro = 500.,
  max_iter = 50,
  axial_integration_step_factor = 0.25,
  debug_zero_o2field = False,
  grid_lattice_const = 50,
  calcflow = dict(
    includePhaseSeparationEffect = True,
  )
)

medo2 = dict(
  po2init_r0 = 55., #  mmHg;  po2_init = min(po2init_cutoff, po2_init_r0 + v->r * po2_init_dr)
  po2init_dr = 1., #  mmHg / um, these values are from yaseen 2011 for rat brain
  po2init_cutoff = 100., # mmHg; maximal attained value
  solubility_plasma = 3.1e-5,
  c0 = 0.5,
  S_n = 2.7,
  S_p50 = 27.,
  D_plasma = 2000.,
  solubility_tissue = 2.8e-5,
  #dM = 0.05,
  rd_norm = 200.,
  rd_tum  = 50.,
  rd_necro = 500.,
  max_iter = 50,
  axial_integration_step_factor = 0.25,
  debug_zero_o2field = False,
  grid_lattice_const = 50.,
  calcflow = dict(
    includePhaseSeparationEffect = True,
  ),
  massTransferCoefficientModelNumber = 1,
  conductivity_coeff1 = c1_,
  conductivity_coeff2 = c2_,
  conductivity_coeff3 = c3_,
  detailedO2name = 'medo2',
)


michaelismenten_consumption = dict(
  num_threads = 6,
  po2init_r0 = 55., #  mmHg;  po2_init = min(po2init_cutoff, po2_init_r0 + v->r * po2_init_dr)
  po2init_dr = 1., #  mmHg / um, these values are from yaseen 2011 for rat brain
  po2init_cutoff = 100., # mmHg; maximal attained value; radius at cutoff is 45 um
  solubility_plasma = 3.1e-5,
  c0 = 0.5,
  sat_curve_exponent = 2.7,
  sat_curve_p50 = 27.,
  D_plasma = 2750.,
  solubility_tissue = 2.8e-5, # ml O2 / (ml tissue mmHg)
  #dM = 0.05,
  rd_norm = 150.,
  rd_tum  = 50.,
  rd_necro = 150.,
  max_iter = 50,
  axial_integration_step_factor = 0.25,
  debug_zero_o2field = False,
  grid_lattice_const = 50.,
  michaelis_menten_uptake = True,
  mmcons_k_norm = 4.,
  mmcons_k_tum = 2.,
  mmcons_k_necro = 2.,
  mmcons_m0_norm = 6.1e-5, # value for breast (ml O2 / ml Tissue / s)
  mmcons_m0_tum = 6.1e-5 * 4., # 4 times normal uptake
  mmcons_m0_necro = 0.,
  calcflow = dict(
    includePhaseSeparationEffect = True,
  ),
  massTransferCoefficientModelNumber = 1,
  conductivity_coeff1 = c1_,
  conductivity_coeff2 = c2_,
  conductivity_coeff3 = c3_,
  detailedO2name = 'michaelismenten_consumption',
)
pso = deepcopy(michaelismenten_consumption)
pso.update(
  grid_lattice_const = 10.,
  #grid_lattice_size = (30,30,30),
)
breast = deepcopy(michaelismenten_consumption)
breast.update(
  mmcons_m0_norm = 13.0e-5,
  mmcons_m0_tum  = 13.0e-5 * 4.,
  calcflow = dict(
    inletHematocrit = 0.45,
    includePhaseSeparationEffect = True,
  ),
  detailedO2name = 'breast',
)

breastlow = deepcopy(michaelismenten_consumption)
breastlow.update(
  mmcons_m0_norm = 7.3e-5,
  mmcons_m0_tum  = 7.3e-5 * 4.,
  calcflow = dict(
    inletHematocrit = 0.45,
    includePhaseSeparationEffect = True,
  ),
)

breastfixed = deepcopy(michaelismenten_consumption)
breastfixed.update(
  mmcons_m0_norm = 7.3e-5, # fairly low like estimates by rinneberg and beany
  mmcons_m0_tum  = 7.3e-5 * 4.,
  calcflow = dict(
    viscosityPlasma = 1.2e-6,
    rheology = 'RheologyForHuman',
    inletHematocrit = 0.45,
    includePhaseSeparationEffect = True,
  ),
  # new, now with hopefully better transvascular conductivity model
  detailedO2name = 'breastfixed',
)


breastv2 = deepcopy(michaelismenten_consumption)
breastv2.update(
  mmcons_m0_norm = 13.2e-5, # 7.9e-3 ml O2 / ml / min
  mmcons_m0_tum = 13.2e-5 * 4,
  calcflow = dict(
    viscosityPlasma = 1.2e-6,
    rheology = 'RheologyForHuman',
    inletHematocrit = 0.45,
    includePhaseSeparationEffect = True,
  ),
  # new, now with hopefully better transvascular conductivity model
  conductivity_coeff1 = c1_,
  conductivity_coeff2 = c2_,
  conductivity_coeff3 = c3_,
  grid_lattice_const = 40.,
  axial_integration_step_factor = 0.1,
  convergence_tolerance = 0.01,
  max_iter = 200, # experience values
  detailedO2name = 'breastv2',
)


breastv3 = deepcopy(breastv2)
breastv3.update(            #3.7e-3/60 = 6.2e-5 min to second!
  mmcons_m0_norm = 6.2e-5, # = 3.7e-3 ml O2 / ml / min; # fairly low like estimates by rinneberg and beany
  mmcons_m0_tum = 6.2e-5 * 4,
  detailedO2name = 'breastv3',
  convergence_tolerance = 0.03,
  max_iter = 400, # some need more
  num_threads = cluster_threads,
)

breast_pso = deepcopy(breastv3)
breast_pso['calcflow'].update(
  rheology = 'RheologySecomb2005',
)

thigh = dict(
  num_threads = cluster_threads,
  po2init_r0 = 55., #  mmHg;  po2_init = min(po2init_cutoff, po2_init_r0 + v->r * po2_init_dr)
  po2init_dr = 1., #  mmHg / um, these values are from yaseen 2011 for rat brain
  po2init_cutoff = 100., # mmHg; maximal attained value; radius at cutoff is 45 um
  solubility_plasma = 3.1e-5,
  c0 = 0.5,
  S_n = 2.7,
  S_p50 = 27.,
  D_plasma = 2000.,
  solubility_tissue = 2.8e-5, # ml O2 / (ml tissue mmHg)
  #dM = 0.05,
  rd_norm = 150.,
  rd_tum  = 50.,
  rd_necro = 150.,
  max_iter = 50,
  axial_integration_step_factor = 0.25,
  debug_zero_o2field = False,
  grid_lattice_const = 50.,
  michaelis_menten_uptake = True,
  mmcons_k_norm = 4.,
  mmcons_k_tum = 2.,
  mmcons_k_necro = 2.,
  #mmcons_m0_norm = 6.1e-5, # value for breast (ml O2 / ml Tissue / s)
  mmcons_m0_norm = 6.1e-4, # as assumption for a muscle I increase that
  mmcons_m0_tum = 6.1e-4 * 4., # 4 times normal uptake
  mmcons_m0_necro = 0.,
  calcflow = dict(
    includePhaseSeparationEffect = True,
  ),
  conductivity_coeff1 = c1_,
  conductivity_coeff2 = c2_,
  conductivity_coeff3 = c3_,
)
thigh2 = deepcopy(thigh)
thigh2.update(
  S_n = 2.7,
  S_p50 = 17.,
)
thigh3 = deepcopy(thigh2)
thigh3.update(
  S_n = 2.7,
  S_p50 = 17.,
  mmcons_k_norm = 0.9,
  mmcons_m0_norm = 2 * 6.1e-5,
)
thigh4 = deepcopy(thigh2)
thigh4.update(
  S_n = 2.7,
  S_p50 = 27.,
  mmcons_k_norm = 0.9,
  mmcons_m0_norm = 3 * 6.1e-5,
)
thigh5 = deepcopy(thigh2)
thigh5.update(
  S_n = 2.7,
  S_p50 = 37.,
  mmcons_k_norm = 1.0,
  mmcons_m0_norm = 3 * 6.1e-5,
)
thigh6 = deepcopy(thigh2)
thigh6.update(
  S_n = 2.7,
  S_p50 = 5.3,
  mmcons_k_norm = 1.0,
  mmcons_m0_norm = 3 * 6.1e-5,
)
thigh7 = deepcopy(thigh2)
thigh7.update(
  S_n = 2.7,
  S_p50 = 5.3,
  mmcons_k_norm = 1.0,
  mmcons_m0_norm = 4 * 6.1e-5,
)
numtest = deepcopy(breastv3)
numtest.update(
  max_iter = 2,
  detailedO2name = 'numtest',
)


secombComparison = deepcopy(breastv3)
secombComparison.update(
  massTransferCoefficientModelNumber = 1,
  conductivity_coeff1 = 5.10370371,  # fitted parameters to Nusselt numbers used by Secomb et al.
  conductivity_coeff2 = 4.78847198,
  conductivity_coeff3 = 0.02647925,
  mmcons_m0_norm = 6.2e-6, # EEE to pow of ----6
  detailedO2name = 'secombComparison',
  safety_layer_size = 500.,
)

secombComparison2 = deepcopy(breastv3)
secombComparison2.update(
  massTransferCoefficientModelNumber = 1,
  conductivity_coeff1 = 6.07978304,
  conductivity_coeff2 = 3.9202903,
  conductivity_coeff3 = 0.01977186,
  detailedO2name = 'secombComparison2',
  mmcons_m0_norm = 13.2e-5,# with these parameters the tissue PO2 should decay over a lengthscale of 82 um!
  safety_layer_size = 400.,
)

#  Distribute consumption rate lognormally.
#  In maple the distribution can be seen like so:
#
#  with(Statistics);
#  X := RandomVariable(LogNormal(mu, sigma));
#  f := PDF(X, u);
#  plot(eval(f, {mu = log(1), sigma = .4}), u = 0 .. 6);
#  Variance(X);
#
#  The variance is exp(2 mu + sigma^2)*(exp(sigma^2) - 1).
#  Sigma is approximately the standard deviation of the distribution.
#  Here std(the distribution w. sigma=0.3) = 32.1%
def breastv4(count):
  import random
  def mk1():
    p = deepcopy(breastv3)
    value = random.lognormvariate(mu = math.log(6.2e-5 * 4), sigma = 0.3)
    p.update(
      mmcons_m0_tum = value,
      detailedO2name = 'breastv4'
    )
    return p
  return list(mk1() for i in xrange(count))


breastv3const = deepcopy(breastv3)
breastv3const.update(
  detailedO2name = 'breastv3const',
  approximateInsignificantTransvascularFlux = True,
  po2init_r0 = 39.,
  po2init_dr = 0,
  po2init_cutoff = 39.,
)
breastv3const['calcflow']['includePhaseSeparationEffect'] = False


breastv3hconst = deepcopy(breastv3)
breastv3hconst['calcflow']['includePhaseSeparationEffect'] = False
breastv3hconst['name'] = 'breastv3hconst'

# use RheologySecomb2005 model. It will lower the viscosity for
# thin vessels (r < 10 micron)
fixedv3 = deepcopy(breastv3)  # not yet used, coefficients need adjustment to new model
fixedv3.update(
  massTransferCoefficientModelNumber = 1,
  conductivity_coeff1 = 8.0,
  conductivity_coeff2 = 4.7,
  conductivity_coeff3 = 0.0,
  calcflow = dict(
    viscosityPlasma = 1.2e-6,
    rheology = 'RheologySecomb2005',
    inletHematocrit = 0.45,
    includePhaseSeparationEffect = True,
  ),
  detailedO2name = 'fixedv3',
)

fixed_mtc3 = deepcopy(breastv3) 
fixed_mtc3.update(
  massTransferCoefficientModelNumber = 1,
  conductivity_coeff1 = 8.0,
  conductivity_coeff2 = 4.7,
  conductivity_coeff3 = 0.0,
  detailedO2name = 'fixed_mtc3',
)

mesentry = dict(
  detailedO2name = 'mesentry',
  po2init_r0 = 55., #  mmHg;  po2_init = min(po2init_cutoff, po2_init_r0 + v->r * po2_init_dr)
  po2init_dr = 1., #  mmHg / um, these values are from yaseen 2011 for rat brain
  po2init_cutoff = 100., # mmHg; maximal attained value; radius at cutoff is 45 um
  solubility_plasma = 3.1e-5,
  c0 = 0.5,
  S_n = 3.0,
  S_p50 = 38.,
  D_plasma = 3000.,
  solubility_tissue = 2.0e-5, # ml O2 / (ml tissue mmHg)
  #dM = 0.05,
  rd_norm = 150.,
  rd_tum  = 50.,
  rd_necro = 150.,
  max_iter = 70,
  axial_integration_step_factor = 0.01,
  debug_zero_o2field = False,
  grid_lattice_const = 20.,
  michaelis_menten_uptake = True,
  mmcons_k_norm = 1.,
  mmcons_k_tum = 2.,
  mmcons_k_necro = 2.,
  mmcons_m0_norm = 0.0004, # value for breast (ml O2 / ml Tissue / s)
  mmcons_m0_tum = 6.1e-5 * 4., # 4 times normal uptake
  mmcons_m0_necro = 0.,
  calcflow = dict(
    includePhaseSeparationEffect = True,
    rheology = 'RheologyForRats',
  ),
  massTransferCoefficientModelNumber = 1,
  conductivity_coeff1 = 7.2,
  conductivity_coeff2 = 4.0,
  conductivity_coeff3 = c3_,
)



v3bigass = deepcopy(breastv3)
v3bigass.update(
  grid_lattice_const = 20,
  detailedO2name = 'v3bigass',
  num_threads = 6,
)


breastv6 = deepcopy(breastv3)
breastv6.update(
  mmcons_m0_norm = 6.2e-5, # = 3.7e-3 ml O2 / ml / min; # fairly low like estimates by rinneberg and beany
  mmcons_m0_tum = 41.e-5, # ml O2 / ml / s   = 25 ul / ml / min
  detailedO2name = 'breastv6',
)

breastv7 = deepcopy(breastv3)
breastv7.update(
  mmcons_m0_norm = 6.2e-5, # = 3.7e-3 ml O2 / ml / min; # fairly low like estimates by rinneberg and beany
  mmcons_m0_tum = 82.e-5, # ml O2 / ml / s   = 50 ul / ml / min
  detailedO2name = 'breastv7',
)
#Rinneberg:  Müssten die BASE Rechnungen nochmals mit einem deutlich höheren  
# M _0 als 14.9 microliter O2/ml/min durchgeführt werden
# Ich hatte unabhängig von den Simulationsrechnungen aus dem Kompartimentmodell 
# eine Sauerstoff-Verbrauchsrate von 19 - 26 microliter O2/ml/min für die Brust-Tumoren abgeschätzt,
# Ich habe vor wenigen Tagen eine Arbeit gefunden, die die 
# Sauerstoff-Verbrauchsrate von wachsenden (also lebenden), 
# multizellulären  Tumor-Spheroiden (Grimes et al "A method for 
# estimating the oxygen consumption rate in multicellular tumour 
# spheroids" J R Soc Interface (2013), 11: 20131124,  siehe Anhang) zu 
# 43.7 +/- 8.4 microliter O2/ml/min angibt (human colorectal carcinoma 
# (DLD1) cells). Der PET-Wert für Brusttumore von 6.6 microliter 
# O2/ml/min (s. Tab. 4) ist nach Lammertsma (Br. J. Radiol, 65: 697-700 
# (1992), Fig. 1A-1C, alte Referenz [23]) grob unterschätzt, vermutlich 
# um einen Faktor von ca.  2-3, entsprechend Verbrauchsraten von 
# 13 - 20 microliter O2/ml/min. Es ist wahrscheinlich, dass die 
# Sauerstoff-Verbrauchsrate im Tumor zu klein angesetzt wurde.

# strange, i used 4x normal consumption which is about 24e-5 / s = 14.88 ul / ml / min
colorectal = deepcopy(fixedv3)
colorectal.update(
    mmcons_m0_norm = 7.3e-4,  #43.7e-3/60
    mmcons_mo_tum = 7.3e-4
    )
swine1 = dict(
  po2init_r0 = 55., #  mmHg;  po2_init = min(po2init_cutoff, po2_init_r0 + v->r * po2_init_dr)
  po2init_dr = 1., #  mmHg / um, these values are from yaseen 2011 for rat brain
  po2init_cutoff = 100., # mmHg; maximal attained value; radius at cutoff is 45 um
  solubility_plasma = 3.1e-5,
  c0 = 0.5,
  S_n = 2.7,
  S_p50 = 27.,
  D_tissue = 2000.,
  solubility_tissue = 2.8e-5, # ml O2 / (ml tissue mmHg)
  #dM = 0.05,
  rd_norm = 150.,
  rd_tum  = 50.,
  rd_necro = 150.,
  max_iter = 50,
  axial_integration_step_factor = 0.25,
  debug_zero_o2field = False,
  grid_lattice_const = 50.,
  michaelis_menten_uptake = True,
  tumor_is_isv = True,
  mmcons_k_norm = 4.,
  mmcons_k_tum = 2.,
  mmcons_k_necro = 2.,
  mmcons_m0_norm = 6.1e-5, # value for breast (ml O2 / ml Tissue / s)
  mmcons_m0_tum = 6.1e-5 * 4., # 4 times normal uptake
  mmcons_m0_necro = 0.,
  calcflow = dict(
    includePhaseSeparationEffect = True,
  ),
)
if __name__ == '__main__':
  abc = breastv4(3)
  print breastv4(3)
