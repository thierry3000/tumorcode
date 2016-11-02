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
from os.path import basename, dirname, join, splitext, commonprefix
if __name__=='__main__': sys.path.append(join(dirname(__file__),'../py'))
import h5py

if 0:
  f1 = h5py.File('t1.h5','w')
  f2 = h5py.File('t2.h5','w')

  f1.create_dataset('fubar', data = 1)
  f1.create_dataset('baz', data = 2)

  f2['l1'] = h5py.ExternalLink('t1.h5','fubar')
  f2['l2'] = h5py.ExternalLink('t1.h5','baz')

  f1.close()
  f2.close()

  #f1 = h5py.File('t1.h5','r+')
  f2 = h5py.File('t2.h5','r')
  ds1 = f2['l1']
  ds2 = f2['l2']
  print f2,
  print ds1
  print ds1.file
  print ds2
  print ds2.file
  f1 = h5py.File('t1.h5','r+')
  print f1

import h5files
with h5files.open('t1.h5', 'a') as f1:
  print f1
  with h5files.open('t1.h5', 'a') as f2:
    print f2
  print f1
print f1