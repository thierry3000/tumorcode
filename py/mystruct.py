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
#!/usr/bin/env python2
import copy

class Struct(dict):
    __getattr__ = dict.__getitem__
    def __setattr__(self,k,v): dict.__setitem__(self,k,v)
    def __deepcopy__(self, memo):
      return Struct(copy.deepcopy(dict(self)))
    def updated(self, other):
      dict.update(self, other)
      return self
    def copy(self):
      return Struct(dict.copy(self))
    def __getstate__(self):
      return dict(self)
    def __setstate__(self, d):
      dict.update(self, d)