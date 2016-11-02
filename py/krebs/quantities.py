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
import unum
import unum.units

if not hasattr(unum.units, 'unitless'):
  unum.units.unitless = unum.Unum.unit('unitless', 1)  
  
  unum.units.mlO2  = unum.Unum.unit('mlO2', 0)
  unum.units.ulO2 = unum.Unum.unit('ulO2', 1.e-3 * unum.units.mlO2)
  unum.units.nlO2 = unum.Unum.unit('nlO2', 1.e-6 * unum.units.mlO2)
  
  unum.units.kPa  = unum.Unum.unit('kPa', 1000. * unum.units.Pa)
  unum.units.mmHg = unum.Unum.unit('mmHg', 133.322 * unum.units.Pa)
  
  unum.units.l    = unum.units.L
  unum.units.ml   = unum.Unum.unit('ml', 1.e-3 * unum.units.L)

  unum.units.mlBlood = unum.Unum.unit('mlBlood', 0)
  unum.units.ulBlood = unum.Unum.unit('ulBlood', 1.e-3 * unum.units.mlBlood)
  
  unum.units.mlBlood_per_ml = unum.Unum.unit('mlBlood_per_ml', unum.units.mlBlood / unum.units.ml)
  unum.units.ulBlood_per_ml = unum.Unum.unit('ulBlood_per_ml', unum.units.ulBlood / unum.units.ml)

  unum.units.mlO2_per_ml = unum.Unum.unit('mlO2_per_ml', unum.units.mlO2 / unum.units.ml)
  unum.units.ulO2_per_ml = unum.Unum.unit('ulO2_per_ml', unum.units.ulO2 / unum.units.ml)

  unum.units.percent   = unum.Unum.unit('percent', 0.01)
from unum.units import *

#note:
#  1 ml = 1 cm^3


def FormatUnumLatex(u):
  '''output a latex string for a given Unum number or unit.'''
  u = u._unit
  tbl = dict(
    um = r'\mu m',
    mlO2_per_ml = r'ml\,O_2\,ml^{-1}',
    ulO2_per_ml = r'\mu l\,O_2\,ml^{-1}',
    nlO2_per_ml = r'nl\,O_2\,ml^{-1}',
    mlO2        = r'ml\,O_2',
    ulO2        = r'\mu l\,O_2',
    nlO2        = r'nl\,O_2',
    ml_per_ml = r'ml\,ml^{-1}',
    ul_per_ml = r'\mu l\,ml^{-1}',
    nl_per_ml = r'nl\,ml^{-1}',
    unitless  = '',
    percent   = r'\%',
  )
  def fmtExponent(e):
    if e == 1:
      return ''
    return r'^{%i}' % e
  l = [] # collect string pieces
  items = sorted(u.items(), key = lambda (name, e): -e)
  for name, e in items:
    l.append(
      '%s%s' % (tbl.get(name, name),fmtExponent(e)))
  return '\\,'.join(l)


class Prettyfier(object):
  '''pretty print math symbols and values'''
  known_stuff = [
    ('pressure',        'Blood pressure',             'p',               'kPa'),
    ('mvd',             'Microvascular density',      'MVD', r'mm^{-2}'),
    ('mvd_linedensity', 'Linedensity or MVD',      'MVD_{line}', r'mm^{-2}'),
    ('mvd_exp',         'Spherical sampled MVD',   'MVD_{exp}', r'mm^{-2}'),
    ('mvd_a',           'Arterial MVD', 'MVD_a', r'mm^{-2}'),
    ('mvd_v',           'Venous MVD', 'MVD_v', r'mm^{-2}'),
    ('mvd_c',           'Capillary MVD', 'MVD_c', r'mm^{-2}'),
    ('rBV',             'Regional blood volume', r'rBV', mlBlood / ml),
    ('phi_vessels',     'Regional blood volume', r'rBV', unitless),
    ('phi_a',           'Artery Vessel volume fraction', r'rBV_a', unitless),
    ('phi_v',           'Vein Vessel volume fraction', r'rBV_v', unitless),
    ('phi_c',           'Capillary volume fraction', r'rBV_c', unitless),
    ('rBV_a',           'Arterial $rBV$', r'rBV_a', ''),
    ('rBV_v',           'Venous $rBV$', r'rBV_v', ''),
    ('rBV_c',           'Capillary $rBV$', r'rBV_c', ''),
    ('meanCapillaryDistance', 'Mean Distance between capillaries (from density)', r'meanCapillaryDistance','\mu m'),
    ('mean_r', 'Mean radius', r'mean_r','\mu m'),
    ('tubeoflux',       'Blood O$_2$ flux', r'\iota', 'ml\,O_2\,min^{-1}'),
    ('venous_rBV_fraction', 'Venous $rBV$ fraction', 'vrBV', ''),
    ('flow',            'Blood flow rate', 'q', r'\mu m^3 / s'),
    ('radius',          'Vessel radius', 'r', r'\mu m'),
    ('velocity',        'Blood flow velocity', 'v', r'\mu m / s'),
    ('velocity_mm',     'Blood flow velocity', 'v', r'mm / s'),
    ('shearforce',      'Vessel wall shear stress', 'f', r'Pa'),
    ('maturation',      'Maturation', 'w', ''),
    ('po2',            r'\mpoxy', 'P', r'mmHg'),
    ('hema',            'Hematocrit', 'H', ''),
    ('hematocrit',      'Hematocrit', 'H', ''),
    ('oxyconc',         'Blood oxygen concentration', 'c', 'ml\,O_2 / ml'),
    ('oxyconct',        'Tissue oxygen concentration', 'c_t', 'ml\,O_2 / ml'),
    ('sat',             'Hemoglobin saturation', 'S', r''),
    ('lattice_const',   'Lattice constant', 'h', '\mu m'),
    ('lateral_size',    'Lateral size of the simulation cube', 'L', 'mm'),
    #('mmcons_k'   , r'M-M Half \mpoxy', 'P_{MP50}'   , 'mmHg'),
    #('mmcons_M0'   , r'M-M Max. Rate', 'M_0'   , 'ml\,O2\,ml^{-1}\,s^{-1}'),
    ('sat_via_hb_ratio', 'Tissue blood O$_2$ saturation', 'Y', ''),
    ('po2_tissue',     r'Tissue \mpoxy', 'P_t', r'mmHg'),
    ('extpo2',         r'Near Vessel Tissue \mpoxy', 'P_{t,v}', r'mmHg'),
    ('mro2',            'Metabolic rate of O$_2$ consumption', r'MRO_2', mlO2_per_ml/min), #r'ml\,O_2\,ml^{-1}\,min^{-1}'),
    ('mro2_by_j',       'O$_2$ outflow per volume', 'J_{vout}/V_{total}', r'ml O_2\,ml^{-1}\,min^{-1}'),
    ('local_mro',       'Local rate of O$_2$ consumption', 'M', r'ml\,O_2\,ml^{-1}\,min^{-1}'),
    ('jtv',             'Transvascular flux', r'j_{tv}', mlO2 / cm**2 / min),#r'ml O_2 cm^{-2} min^{-1}'),
    ('gtv',             'Transvascular \mpoxy gradient', r'(P - P_t)/w)', r'mmHg / \mu m'),
    ('vfhb_oxy',        'Oxy RBC vol. frac.',r'\phi_{Ho}', ''),
    ('vfhb_deoxy',      'Deoxy RBC vol. frac.',r'\phi_{Hd}', ''),
    ('vfhb',            'RBC vol. frac.', r'\phi_{H}', ''),
    ('e1',              'Error (tv)', '(J_{vin}-J_{vout}-J_{tv})/J_{vin}', r'\%'),
    ('e2',              'Error (cons)', '(J_{vin}-J_{vout}-J_{cons})/J_{vin}', r'\%'),
    ('e3',              'Error (tv-cons)', '(J_{tv}-J_{cons})/J_{tv})', r'\%'),
    ('Jin_root',        'O$_2$ influx', 'J_{vin}', r'ml O_2 / min'),
    ('Jout_root',       'O$_2$ outflux' , 'J_{vout}', r'ml O_2 / min'),
    ('Jout_tv',         'Total transvascular flux', 'J_{tv}', r'ml O_2 / min'),
    ('Jout_cons',       'Total consumption', 'J_{cons}', r'ml O_2 / min'),
    ('rJin',            'Regional O2 inflow', 'rJ_{in}', ulO2_per_ml/min), 
    ('scaled_rJin',     'Scaled Regional O2 inflow', 'rJ_{in,scaled}', ulO2_per_ml/min), 
    ('tv_cons','',      'M_{tv}', r'ml O_2 / min'),
    ('approximate_tumor_radius', 'Tumor radius', r'R_{tumor}', r'mm'),
    ('dS_dx',           'dS/dx', r'dS/dx','1/mm'),
    ('rBF',             'Regional blood flow','rBF',r'ml\,g^{-1}\,min^{-1}'),
    ('total_perfusion', 'Total perfusion (in-out)', 'rBF_{tot}',r'ml\,g^{-1}\,min^{-1}'),
    ('scaled_rBF',     'Rescaled rBF', 'rBF_{scaled}', r'ml\,g^{-1}\,min^{-1}'),
    ('total_flow_in',   'Blood inflow', 'BF', r'ml/min^{-1}'),
    ('chb',             'Hemoglobin concentration', 'c_{Hb}', '\mu mol/l'),
    ('chb_oxy',         'Oxyhemoglobin concentration', 'c_{HbO}', '\mu mol/l'),
    ('chb_deoxy',       'Deoxyhemoglobin concentration', 'c_{HbD}', '\mu mol/l'),
    ('oef',             'O$_2$ extraction fraction', 'OEF', ''),
    ('scaled_oef',     'rescaled O$_2$ extraction fraction', 'OEF_{scaled}', ''),
    ('sat_vein',        'Saturation (veins)', 'S_v', ''),
    ('sat_capi',        'Saturation (capillaries)', 'S_c', ''),
    ('sat_art',         'Saturation (arteries)', 'S_a', ''),
    ('sat_estimated_by_acv', 'Estimated saturation', r'Y_{(from\, S_i \cdot rBV_i)}', ''),
    ('mtt', 'Mean Transit Time', 'MTT', 's'),
    ('S_rho', 'Surf. area p. vol.', r'S_D', um**-1),
    ('S_rho_over_rBV', 'Surf. area. p. Vessel vol.', r'S_D / rBV', um**-1),
    ('Sin', 'Inflow Saturation', r'\langle S_{in}\rangle_q', 1),
    ('Peff', 'Peff', r'P_{eff}', um/s),
    ('PeffSrho', 'Peff * S * rho', 'P_{eff}\,S\,rho', 1/s),
    ('PeffSrhoMTT', 'Peff * S * rho*MTT', 'P_{eff}\,S_D\,MTT', 1),
    ('kExpFun', 'Sin * (1-exp(-k))/k', 'S_{in} (1-exp(-k))/k', 1),
    ('Y_plus_oef', 'Y + OEF' , r'Y + OEF', 1),
    # vessel generation parameters
    ('vgp_lattice_size', 'Initial lattice size', r'N', ''),
    ('vgp_lattice_spacing', 'Lattice spacing', r'L', r'\mu m'),
    ('vgp_num_hierarchical_iterations', 'Num. refinement steps', r'nr', ''),
    ('vgp_max_sprout_radius_artery', 'Max capillary connection radius (Artery)', r'r_{ac}', r'\mu m'),
    ('vgp_max_sprout_radius_vein', 'Max capillary connection radius (Vein)', r'r_{vc}', r'\mu m'),
    ('vgp_radius_vein', 'Vein tip radius', 'r_{v}', r'\mu m'),
    ('vgp_radius_artery', 'Arterial tip radius', r'r_{a}', r'\mu m'),
    ('vgp_radius_capi', 'Capillary radius', r'r_{v}', r'\mu m'),
    ('vgp_murray_alpha', "Murray's alpha", r'\alpha', ''),
    ('vgp_murray_alpha_artery', 'murray_alpha_artery', '', ''),
    ('vgp_murray_alpha_vein', 'murray_alpha_vein', '', ''),
    ('vgp_generate_more_capillaries', 'generate_more_capillaries', '', ''),
    ('vgp_capillariesUntilLevel', 'capillariesUntilLevel', '', ''),
    # oxygen
    ('o2p_S_n', 'Hill exponent', r'n', ''),
    ('o2p_S_p50', r'Hill \mpoxy at blood saturation $S=0.5$', r'P_{S50}', 'mmHg'),
    ('o2p_alpha_p', 'O$_2$ solubility (plasma)', r'\alpha', r'ml\,O2\,ml^{-1}\,mmHg^{-1}'),
    ('o2p_alpha_t', 'O$_2$ solubility (tissue)', r'\alpha_t', r'ml\,O2\,ml^{-1}\,mmHg^{-1}'),
    ('o2p_c0', 'O$_2$ carrying capacity of RBCs', r'c_0', 'ml\,O_2\,/ml'),
    ('o2p_kD_tissue', 'O$_2$ diffusion constant (tissue)', r'D_t', r'\mu m^2/s'),
    ('o2p_kD_plasma', 'O$_2$ diffusion constant (plasma)', r'D_p', r'\mu m^2/s'),
    ('o2p_mmcons_k_necro' , r'M-M \mpoxy at 50\% of $M_0$ (necrotic)'      , 'P_{M50} (necrotic)', 'mmHg'),
    ('o2p_mmcons_k_norm'  , r'M-M \mpoxy at 50\% of $M_0$ (normal)'        , 'P_{M50} (normal)'  , 'mmHg'),
    ('o2p_mmcons_k_tum'   , r'M-M \mpoxy at 50\% of $M_0$ (tumor)'         , 'P_{M50} (tumor)'   , 'mmHg'),
    ('o2p_mmcons_k'       , 'M-M \mpoxy at 50\% of $M_0$'         , 'P_{M50}'   , 'mmHg'),
    ('o2p_mmcons_m0_necro', 'M-M maximal rate of O$_2$ consumption (necrotic)', 'M_{0} (necrotic)', 'ml\,O2\,ml^{-1}\,s^{-1}'),
    ('o2p_mmcons_m0_norm' , 'M-M maximal rate of O$_2$ consumption (normal)'  , 'M_{0} (normal)'  , 'ml\,O2\,ml^{-1}\,s^{-1}'),
    ('o2p_mmcons_m0_tum'  , 'M-M maximal rate of O$_2$ consumption (tumor)'   , 'M_{0} (tumor)'   , 'ml\,O2\,ml^{-1}\,s^{-1}'),
    ('o2p_mmcons_m0_p_min'            , 'M-M maximal rate of O$_2$ consumption', 'M_{0}', 'ml\,O2\,ml^{-1}\,min^{-1}'),
    ('o2p_mmcons_m0_necro_p_min', 'M-M maximal rate of O$_2$ consumption (necrotic)', 'M_{0} (necrotic)', 'ml\,O2\,ml^{-1}\,min^{-1}'),
    ('o2p_mmcons_m0_norm_p_min' , 'M-M maximal rate of O$_2$ consumption (normal)'  , 'M_{0} (normal)'  , 'ml\,O2\,ml^{-1}\,min^{-1}'),
    ('o2p_mmcons_m0_tum_p_min'  , 'M-M maximal rate of O$_2$ consumption (tumor)'   , 'M_{0} (tumor)'   , 'ml\,O2\,ml^{-1}\,min^{-1}'),
    ('o2p_po2init_cutoff', r'Root \mpoxy maximum', r'P^{(BC)}_c', r'mmHg'),  #P^{(BC)}_0 + \Delta P^{(BC)} r, P^{(BC)}_c
    ('o2p_po2init_dr'    , r'Root \mpoxy slope', r'\Delta P^{(BC)}', r'mmHg / \mu m'),
    ('o2p_po2init_r0'    , r'Root \mpoxy y-intercept', r'P^{(BC)}_0', 'mmHg'),
    ('o2p_rd_necro'     , r'Diffusion radius (necrotic)', 'rD (necro)', r'\mu m'),
    ('o2p_rd_norm'     , r'Diffusion radius (normal)', 'rD (normal)', r'\mu m'),
    ('o2p_rd_tum'     , r'Diffusion radius (tumor)', 'rD (tumor)', r'\mu m'),
    ('o2p_axial_integration_step_factor', r'Axial integration step factor', r'\Delta x / h', r''),
    ('o2p_grid_lattice_const', r'tissue grid constant', r'h', r'\mu m'),
    #('o2p_transvascular_ring_size', r'Transvascular Layer Breadth', r'w/h', r''),
    ('o2p_blood_press_bc', 'Blood pressure boundary values', r'p^{(BC)}(r)', 'kPa'),
    ('o2p_hema_bc', 'Hematocrit at inlets', r'H^{(BC)}', ''),
    ('o2p_po2_bc', r'\mpoxy at inlets', r'P^{(BC)})(r)', 'mmHg'),
    ('o2p_vessel_int_step', r'Blood \mpoxy interation step length', 'h_{v}', '\mu m'),
    ('o2p_tv_coeff', 'O$_2$ mass transfer coefficient', '\gamma', r'\mu m^3\,O2\,\mu m^{-2}\,s^{-1}\,mmHg^{-1}'),
    #('oxymtc',          'O$_2$ Mass Transfer Coefficient', r'\gamma', r'\mu m^3\,O2\,\mu m^{-2} s^{-1} mmHg^{-1}'),
    ('o2p_vesselCompressionFactor', 'Vessel compression factor', r'\xi_{cpr}', ''),
    ('o2p_rmax', 'R Max', r'r_{max}', um),
    ('avg_cap_dist', 'average distance between capillaries', r'\vec{x}_{cap_i}-\vec{x}_{cap_j}','\mu m'),
  ]
  known_stuff = dict((x[0], x[1:]) for x in known_stuff)
#  multipliers = collections.defaultdict(lambda : float(1.))
#  multipliers.update({
#    'mro2' : 60.,
#    'mro2_by_j' : 60.,
#    'jtv' : 60.,
#    'Jin_root' : 60.,
#    'Jout_root' : 60.,
#    'Jout_tv' : 60.,
#    'Jout_cons' : 60.,
#  })

  autoconversion_to_unit = dict(
    mro2 = ulO2_per_ml / min,
    S_rho = 1/mm,
  )

  @staticmethod
  def MathMode(f):
    @staticmethod
    def wrapper(*args, **kwargs):
      u = f(*args, **kwargs)
      return (r'$%s$' % u) if u else ''
    return wrapper

  @staticmethod
  def format_value_list(stuff, formater):
    known_stuff = Prettyfier.known_stuff
    for k, value in stuff:
      name, sym, unit = known_stuff[k]
      formater(name, sym, unit, value)

  @staticmethod
  def get_label(key):
    '''value label in clear text'''
    return Prettyfier.known_stuff[key][0]

  @staticmethod
  def get_sym(key):
    '''value symbol'''
    return Prettyfier.known_stuff[key][1]

  @staticmethod
  def get_unit(key):
    _, unitstr = Prettyfier.get_value_unit(key, 1)
    return unitstr

  @staticmethod
  def get_value_unit(key, value, mathMode = False, brackets = False):
    u =  Prettyfier.known_stuff[key][2]
    if type(u) == unum.Unum: # support physical unit representation by Unum
      displayUnit = Prettyfier.autoconversion_to_unit.get(key, None)
      if displayUnit is not None:
        #print 'conv %s -> %s' % (str(u), str(displayUnit)),
        multiplier = u.asNumber(displayUnit)
        unitstr    = FormatUnumLatex(displayUnit)
        #print '= %s' % (multiplier)
      else:
        multiplier = 1.
        unitstr    = FormatUnumLatex(u)
    else:
      if u is 1 or u is None:
        multiplier, unitstr = 1, ''
      else:
        multiplier, unitstr = 1, u
    
    if mathMode:
      unitstr = Prettyfier.mm(unitstr)
    if brackets:
      unitstr = Prettyfier.br(unitstr)
    value = multiplier * value
    return value, unitstr

  @staticmethod
  def get_bunit(key):
    '''value unit with brackets (no math mode req.)'''
    u = Prettyfier.get_unit(key)
    return (r'[$%s$]' % u) if u else ''

  @staticmethod
  def mm(s):
    return (r'$%s$' % s) if s else ''

  @staticmethod
  def br(s):
    return (r'[%s]' % s) if s else ''

#  @staticmethod
#  def get_multiplier(key):
#    return Prettyfier.multipliers[key]

Prettyfier.get_munit = Prettyfier.MathMode(Prettyfier.get_unit)
Prettyfier.get_msym  = Prettyfier.MathMode(Prettyfier.get_sym)


if __name__ == '__main__':
  # do some testing
  unum.AUTO_NORM = False
  #thunit = unum.units.um**3 / unum.units.s / unum.units.um**3
  #print thunit
  #thesecondunit = thunit.asUnit(unum.units.mlO2_per_ml / unum.units.min)
  #print thesecondunit
  #print FormatUnumLatex(thesecondunit)
  moreUnits = unum.units.mlO2
  print FormatUnumLatex(moreUnits)
  moreUnits = unum.units.mlO2 / unum.units.ml
  print moreUnits
  print moreUnits.asUnit(unum.units.mlO2_per_ml)
  u4 = unum.units.unitless.asUnit(unum.units.percent)
  print u4
  print (1*u4).asNumber()
  print unum.units.unitless.asNumber(unum.units.percent)
  #u1 = 1./unum.units.m 
  #u2 = 1./unum.units.km 
  u1 = unum.units.mlO2_per_ml / unum.units.min
  u2 = unum.units.ulO2_per_ml / unum.units.min
  u3 = unum.units.mlO2_per_ml / unum.units.s
  u4 = unum.units.mlO2 / unum.units.ml / unum.units.min
  print 'u1=',u1, FormatUnumLatex(u1)
  print 'u2=',u2, FormatUnumLatex(u2)
  print 'u3=',u3, FormatUnumLatex(u3)
  print 'u4=',u4, FormatUnumLatex(u4)
  print 'u1 as u2 = ',u1.asUnit(u2)
  print 'u2 as u1 = ',u2.asUnit(u1)
  print 'u3 as u1 = ',u3.asUnit(u1)
  print 'u4 as u1 = ',u4.asUnit(u1)

  u = mlO2_per_ml/ min
  _, mro2unit =  Prettyfier.get_value_unit('mro2', 1,  mathMode = True, brackets = True)
  s2 = r'[$ml\,O_2\,ml^{-1}\,min^{-1}$]'
  
  print (unum.units.um**-1).asNumber(unum.units.cm**2 / unum.units.ml)  
  print FormatUnumLatex(u)
  print mro2unit  
  print s2
  print r'$\mu m\,O2$'

  import matplotlib.pyplot as pyplot
  fig, ax = pyplot.subplots(1,1)
  ax.text(0.5, 0.5, mro2unit)
  pyplot.show()