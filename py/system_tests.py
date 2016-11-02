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

def test_multi():
  ''''
  got:
    Exception KeyError: KeyError(140317151476224,) in 'h5py._objects.ObjectID.__dealloc__' ignored
Exception KeyError: KeyError(140317151476224,) in 'h5py._objects.ObjectID.__dealloc__' ignored
Process PoolWorker-1:
Exception KeyError: Process PoolWorker-2:
KeyError(140317151476224,) in 'Traceback (most recent call last):
h5py._objects.ObjectID.__deall  File "/usr/lib64/python2.7/multiprocessing/process.py", line 258, in _bootstrap
ocTraceback (most recent call last):
__' ignored
  File "/usr/lib64/python2.7/multiprocessing/process.py", line 258, in _bootstrap

in pso.py  
  ''''
  import multiprocessing
  mp_pool = multiprocessing.Pool(2)
  fx = np.array(mp_pool.map(obj, x))
  
if __name__=='__main__':
  
  test_multi()