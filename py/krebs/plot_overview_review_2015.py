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
import matplotlib.pyplot as plt
import h5py

vessel_file = h5py.File('/localdisk1/thierry/output_count_for_rinneberg/vessels-3d-8mm-P10/vessels-3d-8mm-P10-typeA-49x61L130-sample00.h5','r')

length_prob = np.asarray(vessel_file['data/lengths_prob'])

print('loaded file')

expectation_value=0
for i in range(len(length_prob[0,:,0])):
    expectation_value= expectation_value + length_prob[0,i,0]*length_prob[1,i,0]

print('Expectation value: %f' % expectation_value)

plt.semilogy(length_prob[0,:,0],length_prob[1,:,0],'x')
plt.show()