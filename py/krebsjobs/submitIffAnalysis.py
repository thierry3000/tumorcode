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
from os.path import basename
import qsub


def runs_on_client(filename):
  import krebs.plotIff
  import myutils
  import h5py
  dataman = myutils.DataManager(100, krebs.plotIff.GetDataManDataInstances())  
  with h5py.File(filename,'r+') as iff_file:
    krebs.plotIff.ComputeIfpVsIffCorrelationDataLocal(dataman, iff_file)


if not qsub.is_client:
  filenames = sys.argv[1:]
  for i,filename in enumerate(filenames):
    qsub.submit(qsub.func(runs_on_client, filename),
                name = 'job_'+basename(filename),
                num_cpus = 1,
                days = 0.1,
                change_cwd = True)