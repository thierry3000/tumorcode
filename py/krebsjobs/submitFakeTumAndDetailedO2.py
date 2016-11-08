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

if __name__ == '__main__':
  import os.path, sys
  sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'..'))

import os, sys
import qsub
import glob
import h5py
import h5files
import numpy as np
import identifycluster
import krebsutils
from copy import deepcopy
from dicttoinfo import dicttoinfo, Vec
import myutils

import importlib
import imp

import krebs.detailedo2 as detailedo2
import krebsjobs.parameters.parameterSetsO2
import krebsjobs.submitFakeTum
import krebsjobs.parameters.tumorParameters
import krebs.povrayRenderOxygenDetailed
import krebs.analyzeTissueOxygen

import krebs.tumors


def worker_on_client(vessel_fn, tumor_parameters, o2_parameter_set_name, num_threads):
  krebsutils.set_num_threads(num_threads)
  tumor_fn = tumor_parameters['fn_out']
  tend = tumor_parameters['tend']
  pattern1 = 'out0000'
  pattern2 = 'out%04i' % int(round(tend/tumor_parameters['out_intervall']))
  pattern = '|'.join([pattern1,pattern2])
  print tumor_fn, pattern
  #os.system("%s -s '%s'" %  (krebsjobs.submitTumOnlyVessels.exe,dicttoinfo(tumor_parameters)))
  os.system("%s -s '%s'" %  (krebs.tumors.run_faketum,dicttoinfo(tumor_parameters)))
  o2_refs = detailedo2.doit(tumor_fn, pattern, (getattr(krebs.detailedo2Analysis.parameterSetsO2, o2_parameter_set_name), o2_parameter_set_name))
  for ref in o2_refs:
    po2group = h5files.open(ref.fn)[ref.path]
    krebs.analyzeTissueOxygenDetailed.WriteSamplesToDisk(po2group)
  krebs.povrayRenderOxygenDetailed.doit(o2_refs[1].fn, o2_refs[1].path)
  h5files.closeall()


if __name__ == '__main__' and not qsub.is_client:
  num_threads = 2
  argv = qsub.parse_args(sys.argv)

  tumorParameterName = argv[1]
  o2ParameterName = argv[2]
  filenames = argv[3:]

  try:
    if not tumorParameterName in dir(krebsjobs.parameters.tumorParameters):
      raise AssertionError('Unknown tumor parameter set %s!' % tumorParameterName)
    if not o2ParameterName in dir(krebsjobs.parameters.parameterSetsO2):
      raise AssertionError('Unknown o2 parameter set %s!' % o2ParameterName)
    for fn in filenames:
        if not os.path.isfile(fn):
            raise AssertionError('The file %s is not present!'%fn)
  except Exception, e:
    print e.message
    sys.exit(-1)

  #tumor_parameter_sets = [('defaultconfig', tum_only_vessels.defaultconfig)]
  #o2_parameter_set_name = 'breastv2'

  tumorParameterSet = getattr(krebsjobs.parameters.tumorParameters, tumorParameterName)
  o2ParameterSet = getattr(krebsjobs.parameters.parameterSetsO2, o2ParameterName)
  for vessel_fn in filenames:
    #for name, tumor_parameters in tumor_parameter_sets:
    job_name, tumorParameterSet = krebsjobs.submitFakeTum.PrepareConfigurationForSubmission(vessel_fn, tumorParameterName, 'fakeTumDetailed', tumorParameterSet)
    qsub.submit(qsub.func(worker_on_client, vessel_fn, tumorParameterSet, o2ParameterName, num_threads),
                name = ('job_%s_%s' % (tumorParameterName, o2ParameterName)),
                num_cpus = num_threads,
                days = 4.,
                mem = '2500MB',
                change_cwd = True)
