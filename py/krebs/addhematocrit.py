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
"""
Created on Fri Jan 25 14:56:49 2013

@author: -
"""
if __name__ == '__main__':
  import os.path, sys
  sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'..'))

import h5py
import os,sys
import numpy as np
import extensions # for hdf5 support in np.asarray
import krebsutils
import myutils
import posixpath



class DataHematocritSolver(object):
  keywords = [
    'flow_w_hematocrit',
  ]

  def obtain_data(self, dataman, dataname, *args):
    f, args = args[0], args[1:]
    #obtain_data = lambda *args: dataman.obtain_data(args[0], f, args[1:])

    if dataname == 'flow_w_hematocrit':
      vesselgroup = f[args[0]]

      def read(gmeasure, name):
        grp = gmeasure[name]
        d = dict(
          hema = grp['edges/h_hema'],
          flow = grp['edges/h_flow'],
          force = grp['edges/h_force'],
          press = grp['nodes/h_press'],
        )
        for k, v in d.iteritems():
          d[k] = np.asarray(v)
        return d

      def write(gmeasure, name):
        press, flow, force, hema = krebsutils.calc_vessel_hydrodynamics(vesselgroup, True)
        grp = gmeasure.create_group(name)
        egrp = grp.create_group('edges')
        ngrp = grp.create_group('nodes')
        ngrp.create_dataset('h_press', data = press, compression = 9)
        egrp.create_dataset('h_flow', data = flow, compression = 9)
        egrp.create_dataset('h_hema', data = hema, compression = 9)
        egrp.create_dataset('h_force', data = force, compression = 9)

      fm = myutils.MeasurementFile(f, h5files ,prefix='hematocrit_')
      ret = myutils.hdf_data_caching(read, write, fm, args[0].split(posixpath.pathsep), (1,))
      return ret


if __name__ == '__main__':
  fn = sys.argv[1]
  pattern = sys.argv[2]
  with myutils.H5Files() as h5files:
    dataman = myutils.DataManager(2, [DataHematocritSolver(h5files)])
    f = h5files.open(fn, 'r')
    data = dataman('flow_w_hematocrit', f, pattern)
