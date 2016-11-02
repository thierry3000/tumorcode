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
"""
Created on Mon Aug 22 12:57:28 2016

@author: thierry
"""


from os.path import basename, dirname, join, splitext
import os,sys
import h5py
import h5files
import numpy as np
import extensions # for hdf5 support in np.asarray
import myutils
import posixpath
import math
from copy import deepcopy
from krebs.analyzeGeneral   import DataBasicVessel
from krebs.analyzeBloodFlowResistance import ComputeVascularTreeBloodFlowResistances

import krebsutils # import of this must come in front of import of detailedo2 libs because some required initialization stuff on the c++ side (see mainboost.cpp)
if sys.flags.debug:
  iff_cpp = __import__('libiff_d', globals(), locals())
else:
  iff_cpp = __import__('libiff_', globals(), locals())
