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

import random
import operator
from scoop import futures
bound= 100000
data = [random.randint(-1000, 1000) for r in range(bound)]


if __name__ == '__main__':
    # Python's standard serial function
    serialSum = sum(map(abs, data))

    # SCOOP's parallel function
    parallelSum = futures.mapReduce(abs, operator.add, data)

    assert serialSum == parallelSum