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
import os
from copy import deepcopy
from os.path import join,splitext, basename
import h5files

def MakeVesselFilenamePart(fn):
  with h5files.open(fn, mode='a') as f:
    if 'parameters' in f:
      if 'MESSAGE' in f['parameters'].attrs:
        msg = f['parameters'].attrs['MESSAGE']
        ensemble_index = f['parameters'].attrs['ENSEMBLE_INDEX']
        if msg.startswith('vessels-'): msg=msg[len('vessels-'):]
    if 'msg' not in locals():
      msg = "hm"
      ensemble_index = 1
      f['parameters'].attrs['MESSAGE'] = msg 
      f['parameters'].attrs['ENSEMBLE_INDEX'] = ensemble_index
    
  name = '%s-sample%02i' % (msg, ensemble_index)
  return name

def PrepareConfigurationForSubmission(vessel_fn, name, prepend, config_):
  if vessel_fn is not None: 
    dstdir = os.getcwd()  
    c = deepcopy(config_)
    vessel_fn_part = MakeVesselFilenamePart(vessel_fn)
    out_fn = join(dstdir, '%s-%s-%s.h5' % (prepend,vessel_fn_part, name))
    print 'generating tumor run with'
    print '  vessel:', vessel_fn
    print '  output:', out_fn
    print ' paramset:', name
    c['fn_out'] = out_fn
    c['fn_vessel'] = vessel_fn
    c['paramset_name'] = name
    name = splitext(basename(out_fn))[0]
  else:
    dstdir = os.getcwd()  
    c = deepcopy(config_)
    #vessel_fn_part = MakeVesselFilenamePart(vessel_fn)
    out_fn = join(dstdir, '%s-%s.h5' % (prepend,name))
    print 'generating tumor run no vessels'
    print '  output:', out_fn
    print ' paramset:', name
    c['fn_out'] = out_fn
    c['paramset_name'] = name
    name = splitext(basename(out_fn))[0]
  return name, c