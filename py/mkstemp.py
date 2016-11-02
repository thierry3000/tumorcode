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
#! /usr/bin/env python
import tempfile,os
from tempfile import *

class File:
    """
      A temporary file which is deleted when the File object is destroyed.
      Diffierence to the ordinary mkstemp call:
	The file descriptor is immediately closed, but the file remains.
      Therefore subprocessed can open the file again to do IO,
    """
    def __init__(self,*args,**kw):
        self.keep = kw.pop('keep',False)
        self.fd,self.filename = tempfile.mkstemp(*args,**kw)
        os.close(self.fd)
    def remove(self):
      if self.filename and not self.keep:
        os.unlink(self.filename)
        self.filename = self.fd = None
    def __del__(self):
      self.remove()