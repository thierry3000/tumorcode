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
#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os,sys
import h5py
import numpy as np
import vtk
import vtkcommon
import extensions # for asarray with h5py support
import extractVtkFields

def print_array(g):
  print g.shape
  for i in np.ndindex(g.shape):
    print '%s = %03i' % (i, g[i])

f = h5py.File('test.h5','r')
#g = f['test_field']

print '-----------------'
print 'the hdf field as array:'
print_array(np.asarray(f['test_field']))
print '-----------------'  
print 'the dataset made from the array'
(dataset,) = extractVtkFields.Extractor(f, ['test_field']).asVtkDataSets()
print dataset
for i, t in enumerate(
  vtkcommon.vtkIterTuples(dataset.GetCellData().GetArray(0))
  ):
    a = np.zeros((6,), dtype=np.float32)
    dataset.GetCellBounds(i, a)
    a = 0.5 * (a[::2] + a[1::2])
    print '%03i' % t, a
print '-----------------'
print 'the dataset converted back to numpy'
(img,) = vtkcommon.vtkStructuredGridToNumpy(dataset, ['test_field'])
print_array(img)