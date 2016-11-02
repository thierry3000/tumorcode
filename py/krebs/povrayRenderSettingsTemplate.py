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
Created on Thu Aug 11 16:09:05 2016

@author: mwelter
"""

image = dict(
      #res=(3600,3600), # brauche ca. 3600 fuer 150mm bei 600 dpi
      #res = (1800, 1800), # 75 mm 600 dpi
      #aa=4,
      res=(1024, 1024), aa=1,
      background = 1.,
      out_alpha=False,
      num_threads=3,
)

tumor = dict(
      #cam='topdown_slice',
      #cam='corner',
      #clip_zslice=(0.4, 0.6),
      cam_distance_multiplier = 1.0,
      cam='topdown_slice',
      debug=False,
      colored_slice = True,
)

vessels = dict(
      #clip_pie=0,
      #clip_zslice=(0.4, 0.6),
      #cam='corner',
      cam='topdown_slice',
      #cam = 'topdown',
      colored_slice = True,
      background = 1.,
      #ambient_color = (0.5, 0.5, 0.5),
      ambient_color = (0.3, 0.3, 0.3),
)