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
import myutils
import math
from copy import copy, deepcopy
from collections import namedtuple
import krebsutils
from mystruct import Struct


def PressureGradientFromFlowRate(q, r, paramsbf):
  eta = paramsbf['viscosityPlasma']
  eta *= krebsutils.CalcRelativeViscosity(r, paramsbf['inletHematocrit'], paramsbf['rheology'])
  return q*8*eta/(math.pi*(r**4))


_c1, _c2, _c3 = (0.0052*1.1, 16, 0.0078/math.exp(-4./16)*1.6)  # hand picked coefficients based on secomb
_paramsbf = dict(
  viscosityPlasma = 1.2e-6,
  rheology = 'RheologyForHuman',
  inletHematocrit = 0.3,
  includePhaseSeparationEffect = False,
)
_paramspo2 = dict(
  po2init_r0 = 100.,
  po2init_dr = 0,
  po2init_cutoff = 100.,
  c0 = 0.5,
  S_n = 2.7,
  S_p50 = 27.,
  D_tissue = 2410.,
  D_plasma = 2750.,# 2.75 x 10^-5 cm^2 / s
  solubility_tissue = 3.4e-5, # ml O2 / (ml tissue mmHg)
  solubility_plasma = 2.9e-5,
  debug_zero_o2field = False,
  michaelis_menten_uptake = True,
  mmcons_k_norm = 4.,
  mmcons_k_tum = 2.,
  mmcons_k_necro = 2.,
  mmcons_m0_norm = 6.2e-5, #13.2e-5,
  mmcons_m0_tum = 6.2e-5, #13.2e-5,
  mmcons_m0_necro = 0.,
  #transvascular_ring_size = 0.1,
  max_iter = 200,
  axial_integration_step_factor = 0.1,
  safety_layer_size = -0.1,
  grid_lattice_const = 40.,
  conductivity_coeff1 = _c1,
  conductivity_coeff2 = _c2,
  conductivity_coeff3 = _c3,
  loglevel = 1,
)

_paramstube = dict(
  size = (401, 41, 41),
  scale = 10.,  # gives 4 x 0.4 x 0.4 mm^3 length
  pgrad = PressureGradientFromFlowRate(3.3e6, 50., _paramsbf),
  r = 50.,
  ld_type = 'quad',
  direction_mode = 0,
)

basecase = namedtuple('Params', 'paramsbf, paramspo2, paramsTube')(_paramsbf, _paramspo2, _paramstube)

basecase.paramspo2['massTransferCoefficientModelNumber'] = 1
basecase.paramspo2['conductivity_coeff1'] = 7.2
basecase.paramspo2['conductivity_coeff2'] = 4.0
basecase.paramspo2['conductivity_coeff3'] = 0.0

literature_cases = Struct()

# this is supposed to reproduce the cases from nair 1989
# where o2 diffuses from a perfused pipe through a block
# of rubber.
nair_release = deepcopy(basecase)
nair_release.paramsTube['size'] = (401, 19, 19)
nair_release.paramsTube['r'] = 27*0.5   # FUUUU ... paper gives diameter, must convert to radius
nair_release.paramsTube['pgrad'] = PressureGradientFromFlowRate(3.3e6, 27.*0.5, _paramsbf)  # flow rate of 12 ul /hr
nair_release.paramspo2['tissue_po2_boundary_condition'] = 'dirichlet_yz'
nair_release.paramspo2['po2init_r0'] = 160
nair_release.paramspo2['po2init_cutoff'] = 160
nair_release.paramspo2['mmcons_m0_tum'] = 0
nair_release.paramspo2['mmcons_m0_norm'] = 0
nair_release.paramspo2['solubility_tissue'] = 1
nair_release.paramspo2['D_tissue'] = 0.94
literature_cases['nair_release'] = nair_release

nair_uptake = deepcopy(nair_release)
nair_uptake.paramspo2['po2init_r0'] = 0.
nair_uptake.paramspo2['po2init_cutoff'] = 0.
nair_uptake.paramspo2['tissue_boundary_value']= 160.
literature_cases['nair_uptake'] = nair_uptake

# this should reproduce the base cases in moschandreou 2011
# where Po2 and Sat curves are shown for different radii.
# The flow velocity thereby is held constant at ca 5mm/s.
# see maple script moschandreou_flow for the computation of
# the flow rate from the parameters they presented.
moschandreou_base_case = deepcopy(basecase)
moschandreou_base_case.paramspo2.update(
  extra_tissue_source_linear = -0.625,
  extra_tissue_source_const = 30,
  mmcons_m0_norm = 13.2e-5,
  mmcons_m0_tum = 13.2e-5,
)

_moschandreou_parameter_variation = [ # radius, flow
  (6., 6.6e5),
  (10., 1.8e6),
  (13.5, 3.3e6),
  (27., 1.33e7),
  (40., 2.93e7),
  (50., 4.56e7),
  (100., 5000*100*100*math.pi**2),
]

for i, (r, q) in enumerate(_moschandreou_parameter_variation):
  moschandreou_case = deepcopy(moschandreou_base_case)
  moschandreou_case.paramsTube['size'] = (401, 85, 85)
  moschandreou_case.paramsTube['r'] = r
  moschandreou_case.paramsTube['pgrad'] = PressureGradientFromFlowRate(q, r, _paramsbf)
  literature_cases['moschandreou_case%02i' % i] = moschandreou_case
  globals()['moschandreou_case%02i' % i] = moschandreou_case
  

moschandreou_extra_long = deepcopy(moschandreou_base_case) # extra long
moschandreou_extra_long.paramsTube['size'] = (5001, 85, 85)
moschandreou_extra_long.paramsTube['r'] = 13.5
moschandreou_extra_long.paramsTube['pgrad'] = PressureGradientFromFlowRate(3.3e6, 13.5, _paramsbf)

moschandreou_diag_base = deepcopy(moschandreou_base_case)
moschandreou_diag_base.paramsTube['size'] = (401, 401, 85)
moschandreou_diag_base.paramsTube['ld_type'] = 'fcc'
moschandreou_diag_base.paramsTube['pgrad'] = PressureGradientFromFlowRate(3.3e6, 13.5, _paramsbf)
moschandreou_diag_base.paramsTube['r'] = 13.5

moschandreou_diag = deepcopy(moschandreou_diag_base)
moschandreou_diag.paramsTube['direction_mode'] = 1

#plot_single_capillary(dataman, f['moschandreou_diag'], useInsets = False)


#    ---- thierrys compare case ---
#    this parameter are set according to the program of secomb Version 3
#    http://www.physiology.arizona.edu/people/secomb/greens_c3
#    Secomb uses in inflow of 2nl/min = 1/3*10^5 (mu m)^3/s
#    Note that this vessel is perpendicular to the y-z plane
#    while Secomb's vessel is perpendicular to the x-y plane
#    
#    6.0e-10  tissue D*alpha in cm^3 O2/cm/s/mmHg
#       equaling 0.06 um^3 / um / s / mmHg
thierry_case = deepcopy(nair_release)
#thierry_case.paramsbf['rheology'] = 'RheologySecomb2005'
thierry_case.paramsbf['inletHematocrit'] = 0.4
thierry_case.paramsTube['size']=(12,10,10)
thierry_case.paramsTube['scale'] = 10.
thierry_case.paramsTube['pgrad'] = PressureGradientFromFlowRate(0.033e5, 5., thierry_case.paramsbf)
thierry_case.paramspo2['grid_lattice_const'] = 2.
thierry_case.paramsTube['r'] = 5.
thierry_case.paramspo2['mmcons_m0_norm'] = 0.002
thierry_case.paramspo2['mmcons_k_norm'] = 1
thierry_case.paramspo2['po2init_r0'] = 100
thierry_case.paramspo2['po2init_dr'] = 0
thierry_case.paramspo2['po2init_cutoff'] = 100
thierry_case.paramspo2['S_p50'] = 38
thierry_case.paramspo2['S_n'] = 3.0
thierry_case.paramspo2['c0'] = 0.5
thierry_case.paramspo2['solubility_plasma'] = 3.1e-5
thierry_case.paramspo2['solubility_tissue'] = 3.0e-5
thierry_case.paramspo2['D_tissue'] = 2000
thierry_case.paramspo2['tissue_po2_boundary_condition'] = 'neumann'
thierry_case.paramspo2['conductivity_coeff1'] = 0.01
thierry_case.paramspo2['conductivity_coeff2'] = 16.0 # length scale
thierry_case.paramspo2['conductivity_coeff3'] = 0


test_straight = deepcopy(basecase)
test_straight.paramsTube['ld_type'] = 'fcc'
test_straight.paramsTube['size'] = (401, 401, 41)

test_diag = deepcopy(test_straight)
test_diag.paramsTube['direction_mode'] = 1


def MakeGridSizeSeries(basecase, basename, isDiag):
  def sz(l, h): return (int(math.ceil(l/h))|1)*h/10.+1.0
    
  # system size: 2000 x 400 x 400 micron
  # tube system: 200 x 20 x 20 baseline with a o2 lattice that has 40 as divisor
  ldsize = [
    (200, 41), # 10
    (100, 21), # 20
    (67, 15), # 30
    (50, 11), # 40
    (40, 9), # 50
    (34, 7) # 60
  ] #math.ceil(2000 / 50.0)*50 / 10 + 1
  for i, scale in enumerate([ 10, 20, 30, 40, 50, 60]):
    series = deepcopy(basecase)
    lx, lz = 101, 11
    glx, glz = ldsize[i]
    glz = glz*2+1
    series.paramspo2['grid_lattice_const'] = scale
    series.paramspo2['grid_lattice_size']  = (glx, glx if isDiag else glz, glz)
    series.paramsTube['size'] = (lx, lx if isDiag else lz, lz)
    
    name = '%s-num%02i' % (basename, i)
    globals()[name] = series

hseriesbase1_ = deepcopy(basecase)
hseriesbase1_.paramspo2['max_iter'] = 50
hseriesbase1_.paramspo2['convergence_tolerance'] = 1.e-6
#hseriesbase1_.paramspo2['approximateInsignificantTransvascularFlux'] = True
hseriesbase1_.paramspo2['michaelis_menten_uptake'] = False
hseriesbase1_.paramspo2.update(
  rd_norm = 100.,
  rd_tum  = 100.,
  rd_necro = 100.,
  axial_integration_step_factor = 0.1
)
#hseriesbase1_.paramspo2['debug_zero_o2field'] = True
MakeGridSizeSeries(hseriesbase1_, 'hseries1', isDiag = False)

hseriesbase1diag_ = deepcopy(hseriesbase1_)
hseriesbase1diag_.paramsTube['ld_type'] = 'fcc'
hseriesbase1diag_.paramsTube['direction_mode'] = 1
MakeGridSizeSeries(hseriesbase1diag_, 'hseries1diag', isDiag = True)
